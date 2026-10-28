#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# deduplicate_bins.py
#
# Removes duplicate bins from a single bin set by comparing all pairs for
# contig overlap, keeping only the best-scoring representative where overlap
# exceeds the minimum threshold.
#
# This is needed after iterative consolidation because consolidate_two_sets_of_bins.py
# only compares bins BETWEEN two sets, never within a single set. Duplicate
# genomes from different binners can therefore both survive into binsM if they
# entered from different hybrid combinations during separate merge iterations.
# Running after dereplicate_contigs_in_bins.py ensures shared contigs have
# already been stripped so overlap calculations are based on unique content.
#
# Decision logic for each overlapping pair:
#   1. If completeness differs by < 2% OR contamination differs by < 2%:
#      use CheckM2 as primary decider (metrics too close for CheckM1 alone)
#   2. Otherwise use CheckM1 composite score directly
#   3. If CheckM2 also within 2% on both metrics (or unavailable): use N50
#   4. If N50 identical: use genome size (larger preferred)
#
# Scoring: SCORE = completeness + 5*(100-contamination)
#
# Usage:
# deduplicate_bins.py bin_folder stats_file output_folder
#                     min_completeness max_contamination
#                     [min_overlap_pct] [checkm2_report]
#
# Arguments:
#   bin_folder        : directory containing .fa bin files
#   stats_file        : summarize_checkm.py output (bin comp cont GC lineage N50 size)
#   output_folder     : directory to write deduplicated bins
#   min_completeness  : minimum completeness to consider a bin (e.g. 70)
#   max_contamination : maximum contamination to consider a bin (e.g. 10)
#   min_overlap_pct   : minimum % contig overlap to call two bins duplicates (default 50)
#   checkm2_report    : optional CheckM2 quality_report.tsv for tiebreaking

import sys
import os
import shutil

# ---- Parse arguments --------------------------------------------------------

bin_folder     = sys.argv[1]
stats_file     = sys.argv[2]
output_folder  = sys.argv[3]
min_comp       = float(sys.argv[4])
max_cont       = float(sys.argv[5])
min_overlap    = float(sys.argv[6]) if len(sys.argv) > 6 else 50.0
checkm2_report = sys.argv[7] if len(sys.argv) > 7 else None

# CheckM2 is used as tiebreaker when either completeness or contamination
# differs by less than this threshold between two bins
METRIC_TOLERANCE = 2.0

# ---- Scoring function -------------------------------------------------------

def composite_score(completeness, contamination):
    return completeness + 5.0 * (100.0 - contamination)

# ---- Load CheckM2 data (optional) -------------------------------------------

# Stored as bin_file -> (completeness, contamination, composite_score)
checkm2_data = {}

if checkm2_report and os.path.isfile(checkm2_report):
    with open(checkm2_report) as f:
        for line in f:
            if line.startswith('Name'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                try:
                    name = parts[0]
                    comp = float(parts[1])
                    cont = float(parts[2])
                    checkm2_data[name + '.fa'] = (
                        comp,
                        cont,
                        composite_score(comp, cont)
                    )
                except (ValueError, IndexError):
                    pass

print("Loaded CheckM2 scores for {} bins".format(len(checkm2_data)))

# ---- Load stats file --------------------------------------------------------

print("Loading stats from {}".format(stats_file))

bin_stats   = {}   # bin_file -> (completeness, contamination, n50, size)
bin_summary = {}   # bin_file -> original stats line
header_line = ""

with open(stats_file) as f:
    for line in f:
        if 'completeness' in line.lower():
            header_line = line
            continue
        parts = line.strip().split('\t')
        if len(parts) < 7:
            continue
        bin_name = parts[0]
        bin_file = bin_name + '.fa'
        try:
            comp = float(parts[1])
            cont = float(parts[2])
            n50  = int(parts[5])
            size = int(parts[6])
        except (ValueError, IndexError):
            continue
        bin_stats[bin_file]  = (comp, cont, n50, size)
        bin_summary[bin_file] = line

# Filter to bins meeting quality thresholds
good_bins = {
    b for b, (comp, cont, n50, size) in bin_stats.items()
    if comp >= min_comp and cont <= max_cont
}

print("Total bins in stats: {}  Meeting thresholds: {}".format(
    len(bin_stats), len(good_bins)))

# ---- Load contig sequences for each good bin --------------------------------

print("Loading contig sequences for good bins...")

bin_contigs = {}   # bin_file -> {contig_name: length}

for bin_file in good_bins:
    bin_path = os.path.join(bin_folder, bin_file)
    if not os.path.isfile(bin_path):
        print("WARNING: {} not found, skipping".format(bin_path))
        continue
    bin_contigs[bin_file] = {}
    contig_name = ""
    contig_len  = 0
    with open(bin_path) as f:
        for line in f:
            if line.startswith('>'):
                if contig_name:
                    bin_contigs[bin_file][contig_name] = contig_len
                    contig_len = 0
                contig_name = line[1:].strip().split()[0]
            else:
                contig_len += len(line.strip())
        if contig_name:
            bin_contigs[bin_file][contig_name] = contig_len

# ---- Calculate all pairwise overlaps ----------------------------------------

print("Calculating pairwise contig overlaps for {} bins...".format(
    len(bin_contigs)))

bin_list = sorted(bin_contigs.keys())
overlaps = {}   # (bin_a, bin_b) -> max overlap %

for i, bin_a in enumerate(bin_list):
    for bin_b in bin_list[i+1:]:
        contigs_a = bin_contigs[bin_a]
        contigs_b = bin_contigs[bin_b]

        match_a    = sum(contigs_b[c] for c in contigs_a if c in contigs_b)
        mismatch_a = sum(contigs_a[c] for c in contigs_a if c not in contigs_b)
        match_b    = sum(contigs_a[c] for c in contigs_b if c in contigs_a)
        mismatch_b = sum(contigs_b[c] for c in contigs_b if c not in contigs_a)

        denom_a = match_a + mismatch_a
        denom_b = match_b + mismatch_b
        ratio_a = 100.0 * match_a / denom_a if denom_a > 0 else 0.0
        ratio_b = 100.0 * match_b / denom_b if denom_b > 0 else 0.0

        max_overlap = max(ratio_a, ratio_b)
        if max_overlap >= min_overlap:
            overlaps[(bin_a, bin_b)] = max_overlap

print("Found {} pairs with >= {}% overlap".format(len(overlaps), min_overlap))

# ---- Decision function ------------------------------------------------------

def pick_winner(bin_a, bin_b):
    """
    Return (winner, decision_description) for an overlapping pair.

    Uses CheckM2 as primary decider whenever either completeness or
    contamination differs by less than METRIC_TOLERANCE between the two bins,
    as CheckM1 alone is not reliable enough to distinguish them at that scale.
    """
    comp_a, cont_a, n50_a, size_a = bin_stats[bin_a]
    comp_b, cont_b, n50_b, size_b = bin_stats[bin_b]

    comp_diff = abs(comp_a - comp_b)
    cont_diff = abs(cont_a - cont_b)

    # If either metric is close, defer to CheckM2
    if comp_diff < METRIC_TOLERANCE or cont_diff < METRIC_TOLERANCE:
        c2_a = checkm2_data.get(bin_a, None)
        c2_b = checkm2_data.get(bin_b, None)

        if c2_a is not None and c2_b is not None:
            c2_comp_a, c2_cont_a, c2_score_a = c2_a
            c2_comp_b, c2_cont_b, c2_score_b = c2_b

            c2_comp_diff = abs(c2_comp_a - c2_comp_b)
            c2_cont_diff = abs(c2_cont_a - c2_cont_b)

            # CheckM2 metrics clearly different — use CheckM2 score
            if c2_comp_diff >= METRIC_TOLERANCE or c2_cont_diff >= METRIC_TOLERANCE:
                winner = bin_a if c2_score_a >= c2_score_b else bin_b
                return winner, "CheckM2 scores ({:.1f} vs {:.1f})".format(
                    c2_score_a, c2_score_b)

            # CheckM2 also too close — fall through to N50
            fallback_reason = "CheckM1 and CheckM2 within {}%, ".format(
                METRIC_TOLERANCE)
        else:
            fallback_reason = "CheckM2 unavailable, "
    else:
        # Metrics clearly different — use CheckM1 composite score directly
        score_a = composite_score(comp_a, cont_a)
        score_b = composite_score(comp_b, cont_b)
        winner = bin_a if score_a >= score_b else bin_b
        return winner, "CheckM1 scores ({:.2f} vs {:.2f})".format(
            score_a, score_b)
        fallback_reason = ""

    # N50 tiebreaker
    if n50_a != n50_b:
        winner = bin_a if n50_a > n50_b else bin_b
        return winner, "{}N50 ({} vs {})".format(fallback_reason, n50_a, n50_b)

    # Final tiebreaker: genome size
    winner = bin_a if size_a >= size_b else bin_b
    return winner, "{}genome size ({} vs {})".format(
        fallback_reason, size_a, size_b)

# ---- Greedy deduplication ---------------------------------------------------
# Process pairs sorted by overlap descending so clearest duplicates resolved
# first. Once a bin is marked for removal it cannot win a later comparison.

removed = set()

for (bin_a, bin_b), overlap_pct in sorted(
        overlaps.items(), key=lambda x: x[1], reverse=True):

    if bin_a in removed or bin_b in removed:
        continue

    comp_a, cont_a, n50_a, _ = bin_stats[bin_a]
    comp_b, cont_b, n50_b, _ = bin_stats[bin_b]

    winner, decision = pick_winner(bin_a, bin_b)
    loser = bin_b if winner == bin_a else bin_a

    print("  Overlap {:.1f}%: {} (comp={:.2f}% cont={:.2f}%) vs "
          "{} (comp={:.2f}% cont={:.2f}%) -> keeping {} [{}]".format(
        overlap_pct,
        bin_a, comp_a, cont_a,
        bin_b, comp_b, cont_b,
        winner, decision
    ))
    removed.add(loser)

kept = [b for b in bin_list if b not in removed]

print("\nKept {} bins, removed {} duplicates".format(len(kept), len(removed)))
if removed:
    print("Removed:")
    for b in sorted(removed):
        comp, cont, n50, size = bin_stats[b]
        print("  {} (comp={:.2f}% cont={:.2f}% N50={})".format(
            b, comp, cont, n50))

# Write identified duplicate pairs for strain separation
pairs_path = output_folder + ".strain_pairs.txt"
with open(pairs_path, 'w') as f:
    f.write("# Overlapping bin pairs identified as potential strain duplicates\n")
    f.write("# Format: bin_a\tbin_b\toverlap_pct\tcheckm1_score_a\tcheckm1_score_b\n")
    for (bin_a, bin_b), overlap_pct in sorted(
            overlaps.items(), key=lambda x: x[1], reverse=True):
        score_a = bin_scores.get(bin_a[:-3], 0)  # strip .fa
        score_b = bin_scores.get(bin_b[:-3], 0)
        f.write("{}\t{}\t{:.2f}\t{:.2f}\t{:.2f}\n".format(
            bin_a[:-3], bin_b[:-3], overlap_pct, score_a, score_b))
print("Strain pairs written to {}".format(pairs_path))

# ---- Write output -----------------------------------------------------------

if os.path.isdir(output_folder):
    shutil.rmtree(output_folder)
os.makedirs(output_folder)

bin_ct = 1
new_summary = header_line if header_line else \
    "bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\tbinner\n"

for bin_file in sorted(kept):
    src = os.path.join(bin_folder, bin_file)
    if not os.path.isfile(src):
        print("WARNING: {} missing from bin folder, skipping".format(src))
        continue
    dst = os.path.join(output_folder, "bin.{}.fa".format(bin_ct))
    shutil.copy2(src, dst)

    orig_line = bin_summary.get(bin_file, "")
    parts = orig_line.strip().split('\t')
    parts[0] = "bin.{}".format(bin_ct)
    new_summary += "\t".join(parts) + "\n"
    bin_ct += 1

with open(output_folder + ".stats", 'w') as f:
    f.write(new_summary)

print("Output written to {} ({} bins)".format(output_folder, bin_ct - 1))
print("Stats written to {}.stats".format(output_folder))