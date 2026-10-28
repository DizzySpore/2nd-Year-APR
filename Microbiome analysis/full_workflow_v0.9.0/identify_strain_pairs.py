#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# identify_strain_pairs.py
#
# Identifies potential strain duplicate pairs from a CheckM stats file by
# finding bins that share the same lineage assignment and both meet quality
# thresholds. These pairs are candidates for strain separation by
# separate_strain_bins.py.
#
# Unlike contig-overlap based detection, lineage-based detection works for
# bins that represent the same organism but contain completely different
# contig sets -- the typical result when different binners independently
# partition the same organism's contigs in different ways.
#
# Lineage labels that are too broad to be meaningful for strain-level
# detection (root, Bacteria, Archaea) are excluded since these represent
# unresolved placements rather than true shared taxonomy.
#
# Usage:
#   identify_strain_pairs.py stats_file output_file
#                            min_completeness max_contamination
#
# Arguments:
#   stats_file        : summarize_checkm.py output
#                       columns: bin completeness contamination GC lineage N50 size
#   output_file       : path to write identified pairs (TSV)
#   min_completeness  : minimum completeness to consider a bin (e.g. 70)
#   max_contamination : maximum contamination to consider a bin (e.g. 10)

import sys
import os
from collections import defaultdict

# ---- Arguments --------------------------------------------------------------

stats_file    = sys.argv[1]
output_file   = sys.argv[2]
min_comp      = float(sys.argv[3])
max_cont      = float(sys.argv[4])

# Lineage labels too broad to be meaningful for strain detection
# These represent unresolved placements, not true shared taxonomy
EXCLUDE_LINEAGES = {
    'root',
    'Bacteria',
    'Archaea',
    'Eukaryota',
    ''
}

# ---- Load stats file --------------------------------------------------------

bins = {}   # bin_name -> (completeness, contamination, lineage)

with open(stats_file) as f:
    for line in f:
        if 'completeness' in line.lower():
            continue
        parts = line.strip().split('\t')
        if len(parts) < 5:
            continue
        name = parts[0]
        try:
            comp    = float(parts[1])
            cont    = float(parts[2])
            lineage = parts[4].strip()
        except (ValueError, IndexError):
            continue
        if comp >= min_comp and cont <= max_cont:
            bins[name] = (comp, cont, lineage)

print("Loaded {} bins meeting quality thresholds".format(len(bins)))

# ---- Group bins by lineage --------------------------------------------------

lineage_groups = defaultdict(list)
for bin_name, (comp, cont, lineage) in bins.items():
    if lineage not in EXCLUDE_LINEAGES:
        lineage_groups[lineage].append(bin_name)

# ---- Write pairs from bins sharing the same specific lineage ----------------

pairs_written = 0

with open(output_file, 'w') as f:
    f.write("# Strain pairs identified by shared CheckM lineage\n")
    f.write("# Bins meeting quality thresholds: {}\n".format(len(bins)))
    f.write("# Format: bin_a\tbin_b\tlineage\t"
            "comp_a\tcont_a\tcomp_b\tcont_b\n")

    for lineage, members in sorted(lineage_groups.items()):
        if len(members) < 2:
            continue

        # Write all pairwise combinations within this lineage group
        # Sort members by completeness descending so the better bin is bin_a
        members_sorted = sorted(
            members,
            key=lambda b: bins[b][0] - bins[b][1] * 5,
            reverse=True
        )

        for i, bin_a in enumerate(members_sorted):
            for bin_b in members_sorted[i+1:]:
                comp_a, cont_a, _ = bins[bin_a]
                comp_b, cont_b, _ = bins[bin_b]
                f.write("{}\t{}\t{}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\n".format(
                    bin_a, bin_b, lineage,
                    comp_a, cont_a,
                    comp_b, cont_b))
                pairs_written += 1
                print("  Pair: {} (comp={:.1f}% cont={:.1f}%) vs "
                      "{} (comp={:.1f}% cont={:.1f}%) -- {}".format(
                    bin_a, comp_a, cont_a,
                    bin_b, comp_b, cont_b,
                    lineage))

print("Found {} potential strain pairs by shared lineage".format(pairs_written))
print("Output written to {}".format(output_file))