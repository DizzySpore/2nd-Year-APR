#!/usr/bin/env python3

# Modified from metaWRAP's choose_best_bin.py
# Updated to Python 3 for use with the full metagenomics workflow
#
# This script takes in the reassembled_bins.stats file and chooses the best
# possible version for each bin. Each bin should have 3 versions: the original
# bin (.orig), the "strict" reassembled bin (.strict), and the "permissive"
# reassembled bin (.permissive). This program decides which one is better for
# each bin.
#
# Scoring: SCORE = completeness + 5*(100-contamination)
# Contamination is weighted more heavily as lower contamination is generally
# preferable to slightly higher completion.
# Tiebreaker 1: N50 score
#
# Usage:
# ./choose_best_bin.py stats_file min_completion max_contamination

import sys

best_bins = {}

for line in open(sys.argv[1]):
    if "contamination" in line:
        continue
    cut = line.strip().split("\t")

    # Skip bins below minimum completion or above maximum contamination
    if float(cut[1]) < float(sys.argv[2]) or float(cut[2]) > float(sys.argv[3]):
        continue

    # Bin name is everything before the last dot (strips .orig/.strict/.permissive)
    bin_name = ".".join(cut[0].split(".")[:-1])
    style = cut[0].split(".")[-1]
    score = float(cut[1]) + 5 * (100 - float(cut[2]))
    n50 = int(cut[5])

    if bin_name not in best_bins:
        best_bins[bin_name] = [style, score, n50]
    else:
        if score > best_bins[bin_name][1]:
            best_bins[bin_name] = [style, score, n50]
        elif score == best_bins[bin_name][1] and n50 > best_bins[bin_name][2]:
            best_bins[bin_name] = [style, score, n50]

for i in best_bins:
    print(i + '.' + best_bins[i][0])