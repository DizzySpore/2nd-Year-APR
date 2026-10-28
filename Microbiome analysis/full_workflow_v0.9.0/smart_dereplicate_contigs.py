#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# smart_dereplicate_contigs.py
#
# Marker-aware, coverage-primary replacement for metaWRAP's
# dereplicate_contigs_in_bins.py
#
# For each contig shared between two or more bins, decides which bin it
# truly belongs to based on:
#   (1) Marker excess detection  [hard penalty, always strong]
#   (2) Contig coverage          [primary strain identity signal]
#   (3) Marker gap-filling       [helpful but not strain-diagnostic]
#   (4) Composite CheckM1 score  [final fallback]
#
# ---- Rationale for weight hierarchy -----------------------------------------
#
# The original score-based approach assigned shared contigs to the highest-
# scoring bin regardless of what those contigs contributed. When two bins
# represent different strains of the same species, both bins are nearly
# identical in quality score, and both bins are missing the same marker genes
# (same GCN0 profile). A naive gap-filling approach would then place the same
# marker-carrying contig in whichever bin "needs" it most but since both
# bins have the same gaps, this produces a coin-flip result where patchworks
# of each strain's markers end up in the wrong bin.
#
# Coverage directly measures the relative abundance of the organism each contig
# originated from. A contig with 95x coverage in a sample where strain A has
# 100x median and strain B has 30x median almost certainly belongs to strain A,
# regardless of whether it carries a marker gene. Coverage is therefore
# elevated to the primary signal when marker evidence is only gap-filling.
#
# Marker excess detection (GCN2+) remains a strong negative signal because it
# is reliable: placing a contig in a bin that already has two copies of its
# marker creates unambiguous contamination. Coverage is demoted to secondary
# confirmation in this case since the marker evidence is already clear.
#
# ---- Weight structure -------------------------------------------------------
#
# WEIGHT_EXCESS   = -3.0   hard penalty for creating GCN2+ (reliable signal)
# WEIGHT_MISSING  =  0.75  soft bonus for filling GCN0 gap (helpful, not diagnostic)
# WEIGHT_SINGLE   =  0.25  minimal bonus for GCN1->GCN2 (minor contamination risk)
#
# COVERAGE_WEIGHT_PRIMARY   = 1.5  used when markers only show gap-filling preference
#                                   (no bin has GCN2+ for any marker on this contig)
#                                   strong enough to override gap-filling ties but
#                                   below the excess penalty magnitude (3.0)
#
# COVERAGE_WEIGHT_SECONDARY = 0.3  used when at least one bin has GCN2+ for a marker
#                                   on this contig excess detection already gives
#                                   a clear signal, coverage just confirms
#
# Coverage affinity (per bin) = weight * max(0, 1 - |log2(contig_cov/bin_median)| / 2)
# Maximum coverage contribution in primary mode: 1.5 (at perfect coverage match)
# Maximum coverage contribution in secondary mode: 0.3
#
# ---- Decision hierarchy for shared contigs ----------------------------------
#
#   total_affinity = marker_affinity(excess+gap) + coverage_affinity(adaptive weight)
#   Winner = bin with highest total_affinity
#   Tiebreak: composite CheckM1 score = completeness + 5*(100-contamination)
#
# Contigs with no marker genes use coverage as the primary signal directly.
# If no coverage data available, falls back to composite score (original behaviour).
#
# ---- Expected outcomes ------------------------------------------------------
#
# For two bins of the same species (patchwork strains scenario):
# - Both bins have GCN0 for the same markers -> coverage becomes primary
# - Contigs with coverage close to strain A's median go to strain A
# - Contigs with coverage close to strain B's median go to strain B
# - Each strain ends up with a coherent marker set rather than a patchwork
# - Strain heterogeneity should decrease; completeness should be maintained
#
# ---- Usage ------------------------------------------------------------------
#
#   smart_dereplicate_contigs.py stats_file bin_folder checkm_dir output_folder
#                                [coverage_file]
#
# Arguments:
#   stats_file    : summarize_checkm.py output for binsM
#                   columns: bin completeness contamination GC lineage N50 size
#   bin_folder    : directory containing binsM .fa files
#   checkm_dir    : CheckM lineage_wf output directory for binsM
#                   must contain storage/bin_stats_ext.tsv and bins/{bin}/
#   output_folder : output directory for dereplicated bins (binsO)
#   coverage_file : optional jgi_summarize_bam_contig_depths output
#                   columns: contigName contigLen totalAvgDepth ...

import sys
import os
import re
import ast
import math
import shutil
import statistics

# ---- Arguments --------------------------------------------------------------

stats_file    = sys.argv[1]
bin_folder    = sys.argv[2]
checkm_dir    = sys.argv[3]
output_folder = sys.argv[4]
coverage_file = sys.argv[5] if len(sys.argv) > 5 else None

# ---- Weights ----------------------------------------------------------------

# Marker weights
WEIGHT_EXCESS  = -3.0   # creating GCN2+: reliable contamination signal, hard penalty
WEIGHT_MISSING =  0.75  # filling GCN0:   helpful but not strain-diagnostic
WEIGHT_SINGLE  =  0.25  # GCN1->GCN2:     minor contamination risk, low weight

# Coverage weights adaptive based on quality of marker signal
COVERAGE_WEIGHT_PRIMARY   = 1.5  # used when markers only show gap-filling
                                  # primary strain identity signal
                                  # must be > max gap-filling benefit (0.75)
                                  # and < excess penalty magnitude (3.0)
COVERAGE_WEIGHT_SECONDARY = 0.3  # used when excess detection already gives clear signal
                                  # coverage just confirms the decision

# log2 scale for coverage: contribution falls to 0 at LOG2_ZERO doublings difference
LOG2_ZERO = 2.0   # 4x coverage difference = no coverage contribution

# ---- Helper: composite score ------------------------------------------------

def composite_score(completeness, contamination):
    return completeness + 5.0 * (100.0 - contamination)

# ---- Load bin stats ---------------------------------------------------------

print("Loading bin stats from {}".format(stats_file))

bin_scores = {}   # bin_name -> composite score
bin_stats  = {}   # bin_name -> (completeness, contamination)

with open(stats_file) as f:
    for line in f:
        if 'completeness' in line.lower():
            continue
        parts = line.strip().split('\t')
        if len(parts) < 3:
            continue
        name = parts[0]
        try:
            comp = float(parts[1])
            cont = float(parts[2])
            bin_scores[name] = composite_score(comp, cont)
            bin_stats[name]  = (comp, cont)
        except (ValueError, IndexError):
            continue

print("Loaded stats for {} bins".format(len(bin_scores)))

# ---- Load marker copy numbers from CheckM bin_stats_ext.tsv -----------------
#
# Format: bin_name \t python-repr dict with keys GCN0 GCN1 GCN2 GCN3 GCN4 GCN5+
# Each value is a list of marker family names at that copy number.
# GCN5+ treated as copy count 5 for scoring.

print("Loading marker copy numbers from CheckM output...")

bin_stats_ext_file = os.path.join(checkm_dir, 'storage', 'bin_stats_ext.tsv')

# bin_marker_gcn[bin_name][marker_family] = copy_count (int)
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
    print("Loaded marker GCN data for {} bins".format(len(bin_marker_gcn)))
else:
    print("WARNING: bin_stats_ext.tsv not found"
          "marker-based decisions unavailable, using coverage+score fallback")

# ---- Build contig -> marker mapping from CheckM per-bin files ---------------
#
# CheckM per-bin files under {checkm_dir}/bins/{bin_name}/:
#   genes.gff          GFF3 gene predictions; col 0 = contig name
#   hmmer.analyze.txt  HMMER3 tabular (--tblout);
#                      col 0 = gene/protein (target)
#                      col 2 = HMM profile / marker family (query)
#
# Prodigal gene IDs in GFF attributes: ID={contig_name}_{gene_number}
# Stripping the trailing _{digits} recovers the contig name as fallback.

print("Building contig -> marker mapping from CheckM gene files...")

contig_markers = {}   # contig_name -> set of marker family names

bins_dir = os.path.join(checkm_dir, 'bins')

if os.path.isdir(bins_dir):
    for bin_name in os.listdir(bins_dir):
        bin_checkm_path = os.path.join(bins_dir, bin_name)
        gff_file        = os.path.join(bin_checkm_path, 'genes.gff')
        hmmer_file      = os.path.join(bin_checkm_path, 'hmmer.analyze.txt')

        if not os.path.isfile(gff_file) or not os.path.isfile(hmmer_file):
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
                id_match = re.search(r'ID=([^;]+)', attrs)
                if id_match:
                    gene_id = id_match.group(1)
                    gene_to_contig[gene_id] = contig_name
                # Also index by position as fallback
                pos_key = '{}_{}_{}'.format(contig_name, parts[3], parts[4])
                gene_to_contig[pos_key] = contig_name

        # HMMER tabular -> gene_id to marker families
        gene_to_markers = {}
        with open(hmmer_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) < 3:
                    continue
                gene_id    = parts[0]   # target = protein/gene
                marker_fam = parts[2]   # query  = HMM profile (marker family)
                if gene_id not in gene_to_markers:
                    gene_to_markers[gene_id] = set()
                gene_to_markers[gene_id].add(marker_fam)

        # Combine to get contig -> markers
        for gene_id, markers in gene_to_markers.items():
            contig = gene_to_contig.get(gene_id, None)
            if contig is None:
                # Fallback: strip trailing _{digits}
                stripped = re.sub(r'_\d+$', '', gene_id)
                if stripped in gene_to_contig.values():
                    contig = stripped
            if contig is None:
                continue
            if contig not in contig_markers:
                contig_markers[contig] = set()
            contig_markers[contig].update(markers)

    n_marker_contigs = sum(1 for m in contig_markers.values() if m)
    print("Built marker mapping: {} contigs carry at least one marker gene".format(
        n_marker_contigs))
else:
    print("WARNING: CheckM bins directory not found"
          "marker decisions unavailable")

# ---- Load contig coverage ---------------------------------------------------
#
# From jgi_summarize_bam_contig_depths output:
#   contigName  contigLen  totalAvgDepth  sample.bam  sample.bam-var  ...
#
# totalAvgDepth (col 2) is used the mean depth across the whole contig.
#
# Coverage non-uniformity at strain-specific SNP/gene regions is captured
# implicitly: reads from the other strain map poorly to unique regions,
# pulling down the contig mean. This is intentional a shared contig whose
# mean coverage is closer to bin A's median is soft evidence it belongs to
# bin A's strain, and becomes the primary signal when marker evidence is
# only gap-filling.

contig_coverage = {}   # contig_name -> float depth

if coverage_file and os.path.isfile(coverage_file):
    print("Loading contig coverage from {}".format(coverage_file))
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
elif coverage_file:
    print("WARNING: coverage file {} not found"
          "coverage signal unavailable".format(coverage_file))
else:
    print("No coverage file provided using marker+score decisions only")

# ---- Load all bin FASTA files -----------------------------------------------

print("Loading bin FASTA files from {}...".format(bin_folder))

bins = {}   # bin_name -> {contig_name: sequence}

for bin_file in sorted(f for f in os.listdir(bin_folder) if f.endswith('.fa')):
    bin_name = bin_file[:-3]
    bin_path = os.path.join(bin_folder, bin_file)
    bins[bin_name] = {}
    current_name = None
    current_seq  = []
    with open(bin_path) as f:
        for line in f:
            if line.startswith('>'):
                if current_name:
                    bins[bin_name][current_name] = ''.join(current_seq)
                current_name = line[1:].strip().split()[0]
                current_seq  = []
            else:
                current_seq.append(line.strip())
        if current_name:
            bins[bin_name][current_name] = ''.join(current_seq)
    print("  {} : {} contigs".format(bin_name, len(bins[bin_name])))

# ---- Find shared contigs ----------------------------------------------------

contig_owners = {}
for bin_name, contigs in bins.items():
    for contig_name in contigs:
        if contig_name not in contig_owners:
            contig_owners[contig_name] = []
        contig_owners[contig_name].append(bin_name)

shared_contigs = {
    c: owners for c, owners in contig_owners.items() if len(owners) > 1
}

print("Found {} shared contigs".format(len(shared_contigs)))

def compute_bin_median_coverage(bin_name, bins, shared_contigs, contig_coverage):
    """
    Compute median coverage for a bin using non-shared contigs preferentially.
    Falls back to all contigs if no unique contigs are available, since
    after consolidation some bins may share all their contigs with others.
    """
    non_shared = [
        contig_coverage[c]
        for c in bins[bin_name]
        if c not in shared_contigs and c in contig_coverage
    ]
    if len(non_shared) >= 3:
        return statistics.median(non_shared), 'unique_contigs'

    # Fallback: use all contigs including shared ones
    all_depths = [
        contig_coverage[c]
        for c in bins[bin_name]
        if c in contig_coverage
    ]
    if len(all_depths) >= 3:
        return statistics.median(all_depths), 'all_contigs_fallback'
    elif all_depths:
        return statistics.mean(all_depths), 'all_contigs_fallback'
    return None, 'no_coverage_data'


# ---- Compute per-bin median coverage from non-shared contigs ----------------
#
# Non-shared contigs are the "characteristic" contigs of each bin they
# unambiguously belong to that bin and reflect its true abundance.
# Shared contigs are excluded from the median calculation because they are
# ambiguous by definition.

bin_median_coverage = {}   # bin_name -> median depth

if contig_coverage:
    for bin_name in bins:
        median, source = compute_bin_median_coverage(
            bin_name, bins, shared_contigs, contig_coverage)
        if median is not None:
            bin_median_coverage[bin_name] = median
            print("  {} : {:.2f}x  [{}]".format(bin_name, median, source))
        else:
            print("  {} : no coverage data available".format(bin_name))

    print("Computed median coverage for {}/{} bins:".format(
        len(bin_median_coverage), len(bins)))
    for b, med in sorted(bin_median_coverage.items()):
        print("  {} : {:.2f}x".format(b, med))

# ---- Affinity functions -----------------------------------------------------

def has_excess_signal(owners, markers):
    """
    Returns True if ANY competing bin would have GCN2+ for any marker
    on this contig meaning marker excess detection is available and
    reliable. When this is True, use secondary (low) coverage weight.
    """
    if not markers or not bin_marker_gcn:
        return False
    for bin_name in owners:
        gcn = bin_marker_gcn.get(bin_name, {})
        if any(gcn.get(m, 0) >= 2 for m in markers):
            return True
    return False


def marker_affinity(bin_name, markers):
    """
    Marker-based affinity score for placing this contig in bin_name.
    Excess penalty is strong and reliable.
    Gap-filling bonus is soft helpful but not strain-diagnostic.
    """
    if not markers or bin_name not in bin_marker_gcn:
        return 0.0
    gcn   = bin_marker_gcn[bin_name]
    score = 0.0
    for marker in markers:
        copy_count = gcn.get(marker, 0)
        if copy_count == 0:
            score += WEIGHT_MISSING    # gap-fill: helpful but weak
        elif copy_count == 1:
            score += WEIGHT_SINGLE     # second copy: minimal contamination risk
        else:
            score += WEIGHT_EXCESS     # excess: reliable, hard penalty
    return score


def coverage_affinity(bin_name, contig_name, weight):
    """
    Coverage-based affinity: similarity of contig coverage to bin's
    characteristic median coverage using log2 scale.

    Contribution decays linearly from `weight` at perfect match to 0
    at LOG2_ZERO doublings difference (default 4x).

    Uses contig mean depth the non-uniformity caused by strain-specific
    regions pulling down the mean is the intended signal, not a noise source.
    """
    bin_median = bin_median_coverage.get(bin_name, None)
    contig_cov = contig_coverage.get(contig_name, None)

    if bin_median is None or contig_cov is None:
        return 0.0
    if bin_median <= 0 or contig_cov <= 0:
        return 0.0

    log2_ratio = abs(math.log2(contig_cov / bin_median))
    similarity = max(0.0, 1.0 - log2_ratio / LOG2_ZERO)
    return weight * similarity


def total_affinity(bin_name, contig_name, markers, use_excess_signal):
    """
    Combined affinity score for placing this contig in bin_name.
    Coverage weight is adaptive: primary (1.5) when marker signal is
    only gap-filling, secondary (0.3) when excess detection is available.
    """
    m_aff  = marker_affinity(bin_name, markers)
    cov_wt = COVERAGE_WEIGHT_SECONDARY if use_excess_signal \
             else COVERAGE_WEIGHT_PRIMARY
    c_aff  = coverage_affinity(bin_name, contig_name, cov_wt)
    return m_aff + c_aff

# ---- Assign shared contigs --------------------------------------------------

print("\nAssigning shared contigs...")
print("  Weight structure:")
print("    Marker excess (GCN2+)  : {:.2f}  [hard penalty, reliable]".format(
    WEIGHT_EXCESS))
print("    Marker gap-fill (GCN0) : {:.2f}  [soft bonus, not strain-diagnostic]".format(
    WEIGHT_MISSING))
print("    Coverage (primary)     : {:.2f}  [when only gap-filling markers present]".format(
    COVERAGE_WEIGHT_PRIMARY))
print("    Coverage (secondary)   : {:.2f}  [when excess detection available]".format(
    COVERAGE_WEIGHT_SECONDARY))

contigs_to_remove = {bin_name: set() for bin_name in bins}

excess_led     = 0   # excess detection was the deciding factor
coverage_led   = 0   # coverage was the deciding factor (no excess signal)
gap_led        = 0   # gap-filling broke the tie (coverage unavailable/tied)
score_led      = 0   # composite score was the final fallback

for contig_name, owners in shared_contigs.items():
    markers          = contig_markers.get(contig_name, set())
    excess_available = has_excess_signal(owners, markers)

    affinities = {
        b: total_affinity(b, contig_name, markers, excess_available)
        for b in owners
    }

    max_aff    = max(affinities.values())
    candidates = [b for b, a in affinities.items() if a == max_aff]

    if len(candidates) == 1:
        winner = candidates[0]

        # Determine what drove the decision for reporting
        marker_affs   = {b: marker_affinity(b, markers) for b in owners}
        coverage_affs = {
            b: coverage_affinity(
                b, contig_name,
                COVERAGE_WEIGHT_SECONDARY if excess_available
                else COVERAGE_WEIGHT_PRIMARY
            )
            for b in owners
        }
        m_diff = max(marker_affs.values()) - min(marker_affs.values())
        c_diff = max(coverage_affs.values()) - min(coverage_affs.values())

        if excess_available and m_diff > 0:
            excess_led += 1
        elif c_diff > 0:
            coverage_led += 1
        elif m_diff > 0:
            gap_led += 1
        else:
            score_led += 1
    else:
        # Affinity tie use composite CheckM1 score as final fallback
        winner = max(candidates, key=lambda b: bin_scores.get(b, 0.0))
        score_led += 1

    for bin_name in owners:
        if bin_name != winner:
            contigs_to_remove[bin_name].add(contig_name)

total = excess_led + coverage_led + gap_led + score_led
print("\nDecision breakdown ({} shared contigs):".format(total))
print("  Marker excess detection  : {:5d}  ({:.1f}%)".format(
    excess_led,   100 * excess_led   / max(total, 1)))
print("  Coverage (primary mode)  : {:5d}  ({:.1f}%)".format(
    coverage_led, 100 * coverage_led / max(total, 1)))
print("  Marker gap-fill          : {:5d}  ({:.1f}%)".format(
    gap_led,      100 * gap_led      / max(total, 1)))
print("  Score fallback           : {:5d}  ({:.1f}%)".format(
    score_led,    100 * score_led    / max(total, 1)))

# Guard: if dereplication would completely empty any bin,
# these bins are complete subsets -- revert all assignments
# and let lineage-based selection handle them in Step 6D
bins_that_would_empty = [
    b for b in bins
    if len(bins[b]) > 0 and
    all(c in contigs_to_remove[b] for c in bins[b])
]

if bins_that_would_empty:
    print("\nWARNING: Dereplication would completely empty: {}".format(
        ', '.join(bins_that_would_empty)))
    print("These bins are complete subsets of other bins.")
    print("Reverting all assignments -- Step 6D will select the best version.")
    contigs_to_remove = {b: set() for b in bins}
    excess_led = coverage_led = gap_led = score_led = 0

# ---- Write output -----------------------------------------------------------

print("\nWriting dereplicated bins to {}...".format(output_folder))

if os.path.isdir(output_folder):
    shutil.rmtree(output_folder)
os.makedirs(output_folder)

for bin_name in sorted(bins.keys()):
    remove_set = contigs_to_remove[bin_name]
    out_path   = os.path.join(output_folder, bin_name + '.fa')
    kept = removed = 0

    with open(out_path, 'w') as out_f:
        for contig_name, sequence in bins[bin_name].items():
            if contig_name in remove_set:
                removed += 1
            else:
                out_f.write('>{}\n'.format(contig_name))
                for i in range(0, len(sequence), 60):
                    out_f.write(sequence[i:i+60] + '\n')
                kept += 1

    print("  {} : kept {} contigs, removed {} shared".format(
        bin_name, kept, removed))

print("\nDereplication complete.")
print("Output written to {}".format(output_folder))
print("\nExpected validation signals after re-running CheckM:")
print("  - Similar or improved completeness (contigs assigned to correct strain)")
print("  - Reduced contamination (excess markers resolved)")
print("  - Reduced strain heterogeneity (patchwork bins resolved by coverage)")
print("  - High coverage-led % in decision breakdown = coverage signal was useful")
print("  - Low score-fallback % = decisions were evidence-based, not arbitrary")