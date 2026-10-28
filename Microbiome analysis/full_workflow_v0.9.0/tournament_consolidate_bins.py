#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tournament-Style Bin Consolidation with CheckM2 Tiebreaker
Outputs directly to binsO to integrate with downstream strain-pair resolution.
"""

import sys
import os
from collections import defaultdict

if len(sys.argv) < 5:
    print("Usage: python tournament_consolidate_bins.py <output_folder> <min_comp> <max_cont> <bin_set1> [bin_set2] ...")
    sys.exit(1)

output_dir = sys.argv[1]          # Should be "binsO"
min_comp = float(sys.argv[2])
max_cont = float(sys.argv[3])
bin_sets = sys.argv[4:]

print(f"Tournament consolidation started ? Output folder: {output_dir}")
print(f"Using {len(bin_sets)} input bin sets.")

# ====================== Load CheckM2 scores ======================
checkm2_scores = {}   # full bin path ? score

def load_checkm2(bin_set):
    checkm2_path = f"{bin_set}.checkm2/quality_report.tsv"
    if not os.path.isfile(checkm2_path):
        return
    with open(checkm2_path) as f:
        for line in f:
            if line.startswith("Name") or line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 3: continue
            try:
                bin_name = parts[0]
                if not bin_name.endswith('.fa'):
                    bin_name += '.fa'
                comp = float(parts[1])
                cont = float(parts[2]) if len(parts) > 2 else 0.0
                score = comp + 5.0 * (100.0 - cont)
                bin_path = os.path.join(bin_set, bin_name)
                checkm2_scores[bin_path] = score
            except:
                continue

for s in bin_sets:
    load_checkm2(s)

# ====================== Load all candidate bins ======================
all_candidates = []   # (checkm1_score, bin_path, set_name, bin_filename, contigs, stats_line, checkm2_score)

def load_set(bin_set):
    stats_file = f"{bin_set}.stats"
    if not os.path.isfile(stats_file):
        print(f"Warning: No stats file for {bin_set}")
        return

    with open(stats_file) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3: continue
            try:
                comp = float(parts[1])
                cont = float(parts[2])
                if comp < min_comp or cont > max_cont:
                    continue

                bin_name = parts[0]
                if not bin_name.endswith('.fa'):
                    bin_name += '.fa'

                bin_path = os.path.join(bin_set, bin_name)
                if not os.path.isfile(bin_path):
                    continue

                checkm1_score = comp + 5.0 * (100.0 - cont)
                checkm2_score = checkm2_scores.get(bin_path, checkm1_score)

                # Load contigs
                contigs = {}
                with open(bin_path) as bf:
                    current = ""
                    clen = 0
                    for bline in bf:
                        if bline.startswith('>'):
                            if current: contigs[current] = clen
                            current = bline[1:].strip()
                            clen = 0
                        else:
                            clen += len(bline.strip())
                    if current:
                        contigs[current] = clen

                all_candidates.append((checkm1_score, bin_path, bin_set, bin_name,
                                       contigs, line.strip(), checkm2_score))

            except:
                continue

for s in bin_sets:
    load_set(s)

print(f"Loaded {len(all_candidates)} qualified candidate bins.")

# ====================== Overlap & Grouping ======================
def overlap_score(c1, c2):
    if not c1 or not c2: return 0.0
    inter = sum(min(c1.get(k,0), c2.get(k,0)) for k in c1)
    len1 = sum(c1.values())
    len2 = sum(c2.values())
    return max(100.0 * inter / len1, 100.0 * inter / len2) if len1 and len2 else 0.0

genome_groups = defaultdict(list)
n = len(all_candidates)

for i in range(n):
    for j in range(i+1, n):
        if overlap_score(all_candidates[i][4], all_candidates[j][4]) >= 80.0:
            genome_groups[i].append(j)
            genome_groups[j].append(i)

# Select best per group
visited = set()
final_selection = []

for i in range(n):
    if i in visited: continue
    group = [i]
    visited.add(i)
    stack = genome_groups[i][:]
    while stack:
        node = stack.pop()
        if node not in visited:
            visited.add(node)
            group.append(node)
            stack.extend(genome_groups[node])

    # Best bin: CheckM1 primary, CheckM2 secondary
    best = max(group, key=lambda x: (all_candidates[x][0], all_candidates[x][6]))
    final_selection.append(all_candidates[best])

# ====================== Write Output ======================
os.makedirs(output_dir, exist_ok=True)

final_selection.sort(key=lambda x: x[0], reverse=True)

with open(f"{output_dir}.stats", 'w') as f:
    f.write("bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\tbinner\n")
    for idx, (_, path, set_name, _, _, stats_line, _) in enumerate(final_selection, 1):
        new_name = f"bin.{idx}.fa"
        os.system(f"cp '{path}' '{output_dir}/{new_name}'")
        f.write(f"bin.{idx}\t" + "\t".join(stats_line.split("\t")[1:]) + f"\t{set_name}\n")

print(f"\nTournament consolidation finished!")
print(f"? {len(final_selection)} best bins written to {output_dir}/")
print(f"? Stats file: {output_dir}.stats")