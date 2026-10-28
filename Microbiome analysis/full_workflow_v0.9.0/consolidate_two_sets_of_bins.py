# -*- coding: utf-8 -*-
#!/usr/bin/env python3

# Modified from metaWRAP's consolidate_two_sets_of_bins.py
# Updated for the full metagenomics workflow:
# - Updated to Python 3
# - Integer division made explicit for Python 3 compatibility
# - Print statements updated to functions
# - Scoring formula made consistent with choose_best_bin.py:
#   SCORE = completeness + 5*(100-contamination)
#   (original used completeness - contamination*5 which is functionally
#    equivalent but less readable)
# - Added CheckM2 score integration as secondary tiebreaker
#
# Usage:
# ./script bin_folder_1 bin_folder_2 stats_file_1 stats_file_2 output_folder
#          min_completion max_contamination [checkm2_file_1] [checkm2_file_2]
#
# CheckM2 files are optional — if provided they are used as a tiebreaker
# when CheckM1 scores are identical between two bin versions

import sys
import os

c = float(sys.argv[6])
x = float(sys.argv[7])

# Optional CheckM2 quality report files for tiebreaking
checkm2_file_1 = sys.argv[8] if len(sys.argv) > 8 else None
checkm2_file_2 = sys.argv[9] if len(sys.argv) > 9 else None

################################################################################
# Load CheckM2 scores if provided
################################################################################

checkm2_scores_1 = {}
checkm2_scores_2 = {}

def load_checkm2(filepath, scores_dict):
    if filepath is None or not os.path.isfile(filepath):
        return
    with open(filepath) as f:
        for line in f:
            if line.startswith('Name'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                bin_name = parts[0] + '.fa'
                try:
                    completeness = float(parts[1])
                    contamination = float(parts[2]) if len(parts) > 2 else 0.0
                    scores_dict[bin_name] = completeness + 5 * (100 - contamination)
                except (ValueError, IndexError):
                    pass

load_checkm2(checkm2_file_1, checkm2_scores_1)
load_checkm2(checkm2_file_2, checkm2_scores_2)

################################################################################
# Load good bins from stats files
################################################################################

print("Loading list of good bins (comp>{}%, cont<{}%)".format(c, x))

good_bins_1 = {}
good_bins_2 = {}

for line in open(sys.argv[3]):
    if "completeness" in line:
        continue
    cut = line.strip().split('\t')
    if float(cut[1]) > c and float(cut[2]) < x:
        good_bins_1[cut[0] + '.fa'] = None

for line in open(sys.argv[4]):
    if "completeness" in line:
        continue
    cut = line.strip().split('\t')
    if float(cut[1]) > c and float(cut[2]) < x:
        good_bins_2[cut[0] + '.fa'] = None

################################################################################
# Load contig information for each good bin
################################################################################

bins_1 = {}
bins_2 = {}

print("Loading contig info for good bins...")

for bin_file in good_bins_1:
    bins_1[bin_file] = {}
    contig_len = 0
    contig_name = ""
    for line in open(sys.argv[1] + '/' + bin_file):
        if not line.startswith('>'):
            contig_len += len(line.strip())
        else:
            if contig_name != "":
                bins_1[bin_file][contig_name] = contig_len
                contig_len = 0
            contig_name = line[1:].strip()
    if contig_name != "":
        bins_1[bin_file][contig_name] = contig_len

for bin_file in good_bins_2:
    bins_2[bin_file] = {}
    contig_len = 0
    contig_name = ""
    for line in open(sys.argv[2] + '/' + bin_file):
        if not line.startswith('>'):
            contig_len += len(line.strip())
        else:
            if contig_name != "":
                bins_2[bin_file][contig_name] = contig_len
                contig_len = 0
            contig_name = line[1:].strip()
    if contig_name != "":
        bins_2[bin_file][contig_name] = contig_len

################################################################################
# Compare all bin pairs by contig overlap
################################################################################

print("Making all possible comparisons between the two bin sets...")

all_bin_pairs = {}
for bin_1 in good_bins_1:
    all_bin_pairs[bin_1] = {}
    for bin_2 in good_bins_2:
        match_1_length = 0
        match_2_length = 0
        mismatch_1_length = 0
        mismatch_2_length = 0

        for contig in bins_1[bin_1]:
            if contig in bins_2[bin_2]:
                match_1_length += bins_2[bin_2][contig]
            else:
                mismatch_1_length += bins_1[bin_1][contig]

        for contig in bins_2[bin_2]:
            if contig in bins_1[bin_1]:
                match_2_length += bins_1[bin_1][contig]
            else:
                mismatch_2_length += bins_2[bin_2][contig]

        # Use float division explicitly for Python 3 compatibility
        denom_1 = match_1_length + mismatch_1_length
        denom_2 = match_2_length + mismatch_2_length
        ratio_1 = 100.0 * match_1_length / denom_1 if denom_1 > 0 else 0.0
        ratio_2 = 100.0 * match_2_length / denom_2 if denom_2 > 0 else 0.0

        all_bin_pairs[bin_1][bin_2] = max([ratio_1, ratio_2])

################################################################################
# Load completion and contamination scores
################################################################################

print("Loading completion and contamination scores...")

bins_1_stats = {}
bins_2_stats = {}
bins_1_summary = {}
bins_2_summary = {}

for line in open(sys.argv[3]):
    if "completeness" in line:
        bins_1_summary["header"] = line
        continue
    cut = line.strip().split('\t')
    bins_1_stats[cut[0] + '.fa'] = (float(cut[1]), float(cut[2]))
    bins_1_summary[cut[0] + '.fa'] = line

for line in open(sys.argv[4]):
    if "completeness" in line:
        bins_2_summary["header"] = line
        continue
    cut = line.strip().split('\t')
    bins_2_stats[cut[0] + '.fa'] = (float(cut[1]), float(cut[2]))
    bins_2_summary[cut[0] + '.fa'] = line

################################################################################
# Score function — consistent with choose_best_bin.py
# SCORE = completeness + 5*(100-contamination)
################################################################################

def checkm1_score(completeness, contamination):
    return completeness + 5.0 * (100.0 - contamination)

################################################################################
# Go through all good bins and choose the best version
################################################################################

print("Choosing best bin version from each matched pair...")

os.system("mkdir " + sys.argv[5])
new_summary_file = bins_1_summary.get("header", "bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\n")
bins_2_matches = {}
bin_ct = 1

for bin_1 in all_bin_pairs:
    score_1 = checkm1_score(bins_1_stats[bin_1][0], bins_1_stats[bin_1][1])
    found_better = False

    for bin_2 in all_bin_pairs[bin_1]:
        # Require at least 80% contig overlap to consider bins the same genome
        if all_bin_pairs[bin_1][bin_2] < 80:
            continue
        bins_2_matches[bin_2] = None

        score_2 = checkm1_score(bins_2_stats[bin_2][0], bins_2_stats[bin_2][1])

        # Primary comparison: CheckM1 score
        if score_2 > score_1:
            use_bin_2 = True
        elif score_2 == score_1:
            # Tiebreaker 1: CheckM2 score if available
            c2_score_1 = checkm2_scores_1.get(bin_1, None)
            c2_score_2 = checkm2_scores_2.get(bin_2, None)
            if c2_score_1 is not None and c2_score_2 is not None:
                use_bin_2 = c2_score_2 > c2_score_1
            else:
                use_bin_2 = False
        else:
            use_bin_2 = False

        if use_bin_2:
            cmd = "cp {}/{} {}/bin.{}.fa".format(sys.argv[2], bin_2, sys.argv[5], bin_ct)
            new_summary_file += "bin.{}\t{}".format(
                bin_ct, "\t".join(bins_2_summary[bin_2].split("\t")[1:]))
            found_better = True

    if not found_better:
        new_summary_file += "bin.{}\t{}".format(
            bin_ct, "\t".join(bins_1_summary[bin_1].split("\t")[1:]))
        cmd = "cp {}/{} {}/bin.{}.fa".format(sys.argv[1], bin_1, sys.argv[5], bin_ct)

    os.system(cmd)
    bin_ct += 1

################################################################################
# Retrieve bins from second group not found in first group
################################################################################

print("Retrieving bins from second group not found in first group...")

for bin_2 in bins_2_stats:
    if bins_2_stats[bin_2][0] < c or bins_2_stats[bin_2][1] > x:
        continue
    if bin_2 in bins_2_matches:
        continue
    new_summary_file += "bin.{}\t{}".format(
        bin_ct, "\t".join(bins_2_summary[bin_2].split("\t")[1:]))
    cmd = "cp {}/{} {}/bin.{}.fa".format(sys.argv[2], bin_2, sys.argv[5], bin_ct)
    os.system(cmd)
    bin_ct += 1

################################################################################
# Write output stats file
################################################################################

with open(sys.argv[5] + ".stats", 'w') as f:
    f.write(new_summary_file)

print("There were {} bins cherry-picked from the original sets!".format(bin_ct))