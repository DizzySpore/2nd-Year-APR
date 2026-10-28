#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# separate_strain_bins.py
#
# Strain separation by contig reallocation for bins representing the same
# organism binned differently by different binners.
#
# Context:
# When two bins share the same lineage and similar quality scores but contain
# largely non-overlapping contig sets, they likely represent two co-occurring
# strains that different binners separated differently. Each binner drew bin
# boundaries independently based on coverage and tetranucleotide frequencies,
# resulting in two bins that represent the same genomic space but with
# different contig assignments.
#
# Unlike dereplication (which removes contigs appearing in multiple bins
# simultaneously), this script reassigns contigs from a pooled candidate set
# to the bin they most likely belong to. Every contig in both input bins is
# evaluated and assigned to one of the two output bins -- no contigs are lost.
#
# Decision hierarchy (same weight structure as smart_dereplicate_contigs.py):
#
#   1. Marker excess detection (GCN2+) -- hard penalty -3.0
#      If placing this contig in a bin would create a third copy of a marker
#      that bin already has twice, that is a reliable contamination signal.
#
#   2. Coverage similarity -- primary signal when no excess detected (weight 1.5)
#      secondary confirmation when excess detected (weight 0.3)
#      Contig mean coverage reflects relative strain abundance. A contig at
#      95x in a sample where strain A has 100x median and strain B has 30x
#      almost certainly belongs to strain A. Coverage non-uniformity at
#      strain-specific SNP/gene regions is captured implicitly -- reads from
#      the other strain map poorly to unique regions, pulling down the mean.
#
#   3. Marker gap-filling (GCN0) -- soft bonus 0.75
#      Placing a contig in a bin that is missing its marker is helpful but
#      not diagnostic of strain identity. Many shared housekeeping markers
#      will be missing from both bins and the coverage signal is more reliable.
#
#   4. Composite CheckM1 score -- final fallback
#      completeness + 5*(100-contamination)
#
# Outputs:
#   - {output_folder}/{bin_a}.fa  : contigs assigned to strain A
#   - {output_folder}/{bin_b}.fa  : contigs assigned to strain B
#   - {output_folder}/separation_report.txt : per-contig assignment details
#
# After separation, run CheckM on both output bins. Successful separation
# should show:
#   - Maintained or improved completeness in both bins
#   - Reduced contamination in both bins
#   - Reduced strain heterogeneity in both bins
# If one bin collapses below quality thresholds, the two were not separable
# by coverage/marker content and deduplication should remove the worse one.
#
# Usage:
#   separate_strain_bins.py bin_a_file bin_b_file checkm_dir coverage_file
#                           output_folder [min_completeness] [max_contamination]
#
# Arguments:
#   bin_a_file        : path to first bin .fa file
#   bin_b_file        : path to second bin .fa file
#   checkm_dir        : CheckM lineage_wf output dir for binsM
#                       must contain storage/bin_stats_ext.tsv and bins/{bin}/
#   coverage_file     : jgi_summarize_bam_contig_depths output
#                       columns: contigName contigLen totalAvgDepth ...
#   output_folder     : directory to write separated bins
#   min_completeness  : quality threshold for reporting (default 70)
#   max_contamination : quality threshold for reporting (default 10)

import sys
import os
import re
import ast
import math
import shutil
import statistics

# ---- Arguments --------------------------------------------------------------

bin_a_file        = sys.argv[1]
bin_b_file        = sys.argv[2]
checkm_dir        = sys.argv[3]
coverage_file     = sys.argv[4]
output_folder     = sys.argv[5]
min_comp          = float(sys.argv[6]) if len(sys.argv) > 6 else 70.0
max_cont          = float(sys.argv[7]) if len(sys.argv) > 7 else 10.0

bin_a_name = os.path.basename(bin_a_file)[:-3]   # strip .fa
bin_b_name = os.path.basename(bin_b_file)[:-3]

print("Strain separation: {} vs {}".format(bin_a_name, bin_b_name))
print("CheckM dir  : {}".format(checkm_dir))
print("Coverage    : {}".format(coverage_file))
print("Output      : {}".format(output_folder))

# ---- Weights (same as smart_dereplicate_contigs.py) -------------------------

WEIGHT_EXCESS               = -3.0
WEIGHT_MISSING              =  0.75
WEIGHT_SINGLE               =  0.25
COVERAGE_WEIGHT_PRIMARY     =  1.5
COVERAGE_WEIGHT_SECONDARY   =  0.3
LOG2_ZERO                   =  2.0

# ---- Composite score --------------------------------------------------------

def composite_score(completeness, contamination):
    return completeness + 5.0 * (100.0 - contamination)

# ---- Load marker copy numbers from CheckM bin_stats_ext.tsv -----------------

print("\nLoading marker GCN data from CheckM output...")

bin_stats_ext_file = os.path.join(checkm_dir, 'storage', 'bin_stats_ext.tsv')

# bin_marker_gcn[bin_name][marker_family] = copy_count
bin_marker_gcn = {}

if os.path.isfile(bin_stats_ext_file):
    with open(bin_stats_ext_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            tab_pos  = line.index('\t')
            bin_name = line[:tab_pos]
            dict_str = line[tab_pos + 1:]
            try:
                stats_dict = ast.literal_eval(dict_str)
            except Exception as e:
                print("WARNING: could not parse stats for {}: {}".format(
                    bin_name, e))
                continue
            bin_marker_gcn[bin_name] = {}
            for gcn_key in ['GCN0', 'GCN1', 'GCN2', 'GCN3', 'GCN4', 'GCN5+']:
                copy_count = int(gcn_key.replace('GCN', '').replace('+', ''))
                for marker in stats_dict.get(gcn_key, []):
                    bin_marker_gcn[bin_name][marker] = copy_count

    print("Loaded GCN data for {} bins".format(len(bin_marker_gcn)))
    if bin_a_name not in bin_marker_gcn:
        print("WARNING: {} not found in bin_stats_ext.tsv".format(bin_a_name))
    if bin_b_name not in bin_marker_gcn:
        print("WARNING: {} not found in bin_stats_ext.tsv".format(bin_b_name))
else:
    print("WARNING: bin_stats_ext.tsv not found -- "
          "marker decisions unavailable, using coverage only")

# ---- Build contig -> marker mapping from CheckM per-bin gene files ----------
#
# CheckM bins/{bin_name}/genes.gff  -- GFF3; col 0 = contig name
# CheckM bins/{bin_name}/hmmer.analyze.txt -- HMMER3 tabular;
#   col 0 = gene/protein (target), col 2 = HMM profile/marker family (query)

print("Building contig -> marker mapping...")

contig_markers = {}   # contig_name -> set of marker family names

bins_dir = os.path.join(checkm_dir, 'bins')

if os.path.isdir(bins_dir):
    for bin_name in [bin_a_name, bin_b_name]:
        bin_checkm_path = os.path.join(bins_dir, bin_name)
        gff_file        = os.path.join(bin_checkm_path, 'genes.gff')
        hmmer_file      = os.path.join(bin_checkm_path, 'hmmer.analyze.txt')

        if not os.path.isfile(gff_file) or not os.path.isfile(hmmer_file):
            print("WARNING: gene files not found for {} -- "
                  "marker mapping unavailable for this bin".format(bin_name))
            continue

        # GFF3 -> gene_id to contig_name
        gene_to_contig = {}
        with open(gff_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                contig_name = parts[0]
                attrs       = parts[8]
                id_match    = re.search(r'ID=([^;]+)', attrs)
                if id_match:
                    gene_to_contig[id_match.group(1)] = contig_name
                pos_key = '{}_{}_{}'.format(contig_name, parts[3], parts[4])
                gene_to_contig[pos_key] = contig_name

        # HMMER tabular -> gene_id to markers
        gene_to_markers = {}
        with open(hmmer_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) < 3:
                    continue
                gene_id    = parts[0]
                marker_fam = parts[2]
                if gene_id not in gene_to_markers:
                    gene_to_markers[gene_id] = set()
                gene_to_markers[gene_id].add(marker_fam)

        # Combine gene->contig and gene->markers to get contig->markers
        for gene_id, markers in gene_to_markers.items():
            contig = gene_to_contig.get(gene_id, None)
            if contig is None:
                stripped = re.sub(r'_\d+$', '', gene_id)
                if stripped in gene_to_contig.values():
                    contig = stripped
            if contig is None:
                continue
            if contig not in contig_markers:
                contig_markers[contig] = set()
            contig_markers[contig].update(markers)

    print("Built marker mapping for {} contigs".format(len(contig_markers)))
else:
    print("WARNING: CheckM bins directory not found -- "
          "marker decisions unavailable")

# ---- Load contig coverage ---------------------------------------------------

print("Loading contig coverage from {}...".format(coverage_file))

contig_coverage = {}   # contig_name -> float depth

if os.path.isfile(coverage_file):
    with open(coverage_file) as f:
        for line in f:
            if line.startswith('contigName') or line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            try:
                contig_coverage[parts[0]] = float(parts[2])
            except (ValueError, IndexError):
                continue
    print("Loaded coverage for {} contigs".format(len(contig_coverage)))
else:
    print("WARNING: coverage file not found -- "
          "coverage signal unavailable")

# ---- Load both bin FASTA files ----------------------------------------------

def load_fasta(filepath):
    """Load a FASTA file into {contig_name: sequence} dict."""
    contigs = {}
    current_name = None
    current_seq  = []
    with open(filepath) as f:
        for line in f:
            if line.startswith('>'):
                if current_name:
                    contigs[current_name] = ''.join(current_seq)
                current_name = line[1:].strip().split()[0]
                current_seq  = []
            else:
                current_seq.append(line.strip())
        if current_name:
            contigs[current_name] = ''.join(current_seq)
    return contigs

print("\nLoading bin FASTA files...")

contigs_a = load_fasta(bin_a_file)
contigs_b = load_fasta(bin_b_file)

print("  {} : {} contigs  total length {:.2f} Mbp".format(
    bin_a_name, len(contigs_a),
    sum(len(s) for s in contigs_a.values()) / 1e6))
print("  {} : {} contigs  total length {:.2f} Mbp".format(
    bin_b_name, len(contigs_b),
    sum(len(s) for s in contigs_b.values()) / 1e6))

# Pool all contigs from both bins
all_contigs = {}
all_contigs.update(contigs_a)
all_contigs.update(contigs_b)

# Track original bin membership for reporting
original_bin = {}
for c in contigs_a:
    original_bin[c] = bin_a_name
for c in contigs_b:
    # If a contig appears in both (unlikely but possible), bin_b overwrites
    # -- these true shared contigs will be handled by coverage/marker signal
    original_bin[c] = bin_b_name if c not in contigs_a else 'shared'

print("  Total pooled contigs: {}".format(len(all_contigs)))

shared_count = sum(1 for c in contigs_a if c in contigs_b)
if shared_count > 0:
    print("  {} contigs present in both bins "
          "(will be assigned by coverage/marker)".format(shared_count))

# ---- Compute per-bin median coverage from non-overlapping contigs -----------
#
# Use contigs unique to each bin as the "characteristic" coverage reference.
# These unambiguously belong to each bin and reflect its true strain abundance.

unique_to_a = set(contigs_a.keys()) - set(contigs_b.keys())
unique_to_b = set(contigs_b.keys()) - set(contigs_a.keys())

def compute_median_coverage(contig_set, coverage_dict, fallback_set=None):
    depths = [coverage_dict[c] for c in contig_set if c in coverage_dict]
    if len(depths) < 3 and fallback_set is not None:
        depths = [coverage_dict[c] for c in fallback_set if c in coverage_dict]
    if len(depths) >= 3:
        return statistics.median(depths)
    elif depths:
        return statistics.mean(depths)
    return None

median_a = compute_median_coverage(unique_to_a, contig_coverage, 
                                    fallback_set=contigs_a.keys())
median_b = compute_median_coverage(unique_to_b, contig_coverage,
                                    fallback_set=contigs_b.keys())

print("\nCharacteristic coverage from unique contigs:")
print("  {} : {:.2f}x  ({} unique contigs used)".format(
    bin_a_name,
    median_a if median_a is not None else 0.0,
    len([c for c in unique_to_a if c in contig_coverage])))
print("  {} : {:.2f}x  ({} unique contigs used)".format(
    bin_b_name,
    median_b if median_b is not None else 0.0,
    len([c for c in unique_to_b if c in contig_coverage])))

bin_median_coverage = {}
if median_a is not None:
    bin_median_coverage[bin_a_name] = median_a
if median_b is not None:
    bin_median_coverage[bin_b_name] = median_b

# ---- Affinity functions (identical to smart_dereplicate_contigs.py) ---------

def has_excess_signal(markers, bins):
    """True if any competing bin has GCN2+ for any marker on this contig."""
    if not markers or not bin_marker_gcn:
        return False
    for bin_name in bins:
        gcn = bin_marker_gcn.get(bin_name, {})
        if any(gcn.get(m, 0) >= 2 for m in markers):
            return True
    return False


def marker_affinity(bin_name, markers):
    if not markers or bin_name not in bin_marker_gcn:
        return 0.0
    gcn   = bin_marker_gcn[bin_name]
    score = 0.0
    for marker in markers:
        copy_count = gcn.get(marker, 0)
        if copy_count == 0:
            score += WEIGHT_MISSING
        elif copy_count == 1:
            score += WEIGHT_SINGLE
        else:
            score += WEIGHT_EXCESS
    return score


def coverage_affinity(bin_name, contig_name, weight):
    bin_median  = bin_median_coverage.get(bin_name, None)
    contig_cov  = contig_coverage.get(contig_name, None)
    if bin_median is None or contig_cov is None:
        return 0.0
    if bin_median <= 0 or contig_cov <= 0:
        return 0.0
    log2_ratio = abs(math.log2(contig_cov / bin_median))
    similarity = max(0.0, 1.0 - log2_ratio / LOG2_ZERO)
    return weight * similarity

# ---- Assign all pooled contigs to bin_a or bin_b ---------------------------

print("\nAssigning {} pooled contigs...".format(len(all_contigs)))

competing_bins = [bin_a_name, bin_b_name]

# Assignment tracking
assigned_to_a   = []
assigned_to_b   = []
report_lines    = []

decision_counts = {
    'excess'   : 0,
    'coverage' : 0,
    'gap_fill' : 0,
    'score'    : 0,
    'unchanged': 0
}

# Load composite scores from bin_stats_ext for fallback
bin_composite = {}
if os.path.isfile(bin_stats_ext_file):
    with open(bin_stats_ext_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            tab_pos  = line.index('\t')
            bin_name = line[:tab_pos]
            dict_str = line[tab_pos + 1:]
            try:
                d    = ast.literal_eval(dict_str)
                comp = d.get('Completeness', 0)
                cont = d.get('Contamination', 0)
                bin_composite[bin_name] = composite_score(comp, cont)
            except Exception:
                pass

for contig_name, sequence in all_contigs.items():
    markers       = contig_markers.get(contig_name, set())
    excess_avail  = has_excess_signal(markers, competing_bins)
    cov_weight    = COVERAGE_WEIGHT_SECONDARY if excess_avail \
                    else COVERAGE_WEIGHT_PRIMARY

    aff_a = (marker_affinity(bin_a_name, markers) +
             coverage_affinity(bin_a_name, contig_name, cov_weight))
    aff_b = (marker_affinity(bin_b_name, markers) +
             coverage_affinity(bin_b_name, contig_name, cov_weight))

    # Determine decision driver
    m_aff_a = marker_affinity(bin_a_name, markers)
    m_aff_b = marker_affinity(bin_b_name, markers)
    c_aff_a = coverage_affinity(bin_a_name, contig_name, cov_weight)
    c_aff_b = coverage_affinity(bin_b_name, contig_name, cov_weight)

    m_diff = abs(m_aff_a - m_aff_b)
    c_diff = abs(c_aff_a - c_aff_b)

    if aff_a > aff_b:
        winner = bin_a_name
    elif aff_b > aff_a:
        winner = bin_b_name
    else:
        # Tie -- use composite score as fallback
        score_a = bin_composite.get(bin_a_name, 0)
        score_b = bin_composite.get(bin_b_name, 0)
        winner  = bin_a_name if score_a >= score_b else bin_b_name

    # Classify decision for reporting
    if aff_a != aff_b:
        if excess_avail and m_diff > 0:
            decision = 'excess'
            decision_counts['excess'] += 1
        elif c_diff > 0:
            decision = 'coverage'
            decision_counts['coverage'] += 1
        elif m_diff > 0:
            decision = 'gap_fill'
            decision_counts['gap_fill'] += 1
        else:
            decision = 'score'
            decision_counts['score'] += 1
    else:
        decision = 'score'
        decision_counts['score'] += 1

    # Track if contig moved relative to original bin
    original = original_bin.get(contig_name, 'unknown')
    moved    = (original != winner and original != 'shared')
    if not moved:
        decision_counts['unchanged'] += 1

    contig_cov  = contig_coverage.get(contig_name, None)
    cov_str     = "{:.1f}x".format(contig_cov) if contig_cov else "n/a"
    contig_len  = len(sequence)

    report_lines.append("{}\t{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}\t{}\t{}".format(
        contig_name,
        contig_len,
        cov_str,
        original,
        winner,
        aff_a,
        aff_b,
        decision,
        'MOVED' if moved else 'kept'
    ))

    if winner == bin_a_name:
        assigned_to_a.append((contig_name, sequence))
    else:
        assigned_to_b.append((contig_name, sequence))

total = len(all_contigs)
moved_total = total - decision_counts['unchanged']

print("\nAssignment complete:")
print("  Assigned to {} : {} contigs  ({:.2f} Mbp)".format(
    bin_a_name, len(assigned_to_a),
    sum(len(s) for _, s in assigned_to_a) / 1e6))
print("  Assigned to {} : {} contigs  ({:.2f} Mbp)".format(
    bin_b_name, len(assigned_to_b),
    sum(len(s) for _, s in assigned_to_b) / 1e6))
print("  Contigs moved from original bin: {} ({:.1f}%)".format(
    moved_total, 100.0 * moved_total / max(total, 1)))
print("\nDecision breakdown:")
print("  Marker excess detection  : {:5d}  ({:.1f}%)".format(
    decision_counts['excess'],
    100 * decision_counts['excess']   / max(total, 1)))
print("  Coverage (primary mode)  : {:5d}  ({:.1f}%)".format(
    decision_counts['coverage'],
    100 * decision_counts['coverage'] / max(total, 1)))
print("  Marker gap-fill          : {:5d}  ({:.1f}%)".format(
    decision_counts['gap_fill'],
    100 * decision_counts['gap_fill'] / max(total, 1)))
print("  Score fallback           : {:5d}  ({:.1f}%)".format(
    decision_counts['score'],
    100 * decision_counts['score']    / max(total, 1)))

# ---- Write output -----------------------------------------------------------

if os.path.isdir(output_folder):
    shutil.rmtree(output_folder)
os.makedirs(output_folder)

def write_fasta(filepath, contigs):
    with open(filepath, 'w') as f:
        for contig_name, sequence in contigs:
            f.write('>{}\n'.format(contig_name))
            for i in range(0, len(sequence), 60):
                f.write(sequence[i:i+60] + '\n')

out_a = os.path.join(output_folder, bin_a_name + '.fa')
out_b = os.path.join(output_folder, bin_b_name + '.fa')

write_fasta(out_a, assigned_to_a)
write_fasta(out_b, assigned_to_b)

# Write per-contig report
report_path = os.path.join(output_folder, 'separation_report.txt')
with open(report_path, 'w') as f:
    f.write("# Strain separation report: {} vs {}\n".format(
        bin_a_name, bin_b_name))
    f.write("# Median coverage: {} = {:.2f}x  {} = {:.2f}x\n".format(
        bin_a_name, median_a if median_a else 0,
        bin_b_name, median_b if median_b else 0))
    f.write("# Contigs moved: {}/{} ({:.1f}%)\n".format(
        moved_total, total, 100.0 * moved_total / max(total, 1)))
    f.write("#\n")
    f.write("contig\tlength\tcoverage\toriginal_bin\tassigned_bin\t"
            "affinity_a\taffinity_b\tdecision\tstatus\n")
    for line in sorted(report_lines):
        f.write(line + '\n')

print("\nOutput written:")
print("  {}".format(out_a))
print("  {}".format(out_b))
print("  {}".format(report_path))
print("\nNext steps:")
print("  Run CheckM on output bins to evaluate separation quality.")
print("  Successful separation: maintained completeness + reduced")
print("  contamination and strain heterogeneity in both bins.")
print("  If one bin collapses below quality thresholds, the strains")
print("  could not be reliably separated -- remove the worse bin.")