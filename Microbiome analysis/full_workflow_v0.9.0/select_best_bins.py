# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
select_best_bins.py - Global multi-set bin selector using a dual-layered uniform rule:
                     1. Native metaWRAP length-based contig overlap.
                     2. Uniform K-mer Jaccard similarity across all pairwise bins.
"""

import argparse
import csv
import shutil
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Scoring & Sequence Helpers
# ---------------------------------------------------------------------------

def checkm1_score(completeness, contamination):
    return completeness - 5.0 * contamination

def get_fast_kmers(fa_path, k=31, step=10):
    """
    Uniformly sketches a FASTA file into a downsampled set of K-mers.
    Using a step size keeps memory usage low while preserving a highly accurate 
    genomic similarity profile.
    """
    kmers = set()
    try:
        with open(fa_path, encoding="utf-8", errors="replace") as fh:
            current_seq = []
            for line in fh:
                if line.startswith(">"):
                    if current_seq:
                        seq_str = "".join(current_seq).upper()
                        for i in range(0, len(seq_str) - k + 1, step):
                            kmers.add(seq_str[i:i+k])
                        current_seq = []
                else:
                    current_seq.append(line.strip())
            if current_seq:
                seq_str = "".join(current_seq).upper()
                for i in range(0, len(seq_str) - k + 1, step):
                    kmers.add(seq_str[i:i+k])
    except OSError:
        pass
    return kmers

# ---------------------------------------------------------------------------
# Stats & Sequence Parsing
# ---------------------------------------------------------------------------

def load_stats_and_lengths(work_dir):
    stats_files = sorted(work_dir.glob("*.stats"))
    if not stats_files:
        sys.exit("ERROR: No *.stats files found in {}".format(work_dir))

    SKIP_STEMS = {"binsO", "binsM", "binsM1"}
    stats_files = [sf for sf in stats_files if sf.stem not in SKIP_STEMS]

    rows = []
    bin_contigs = {}

    print("Parsing bin FASTA sequences to map exact base-pair lengths...", file=sys.stderr)

    for sf in stats_files:
        file_stem = sf.stem   
        with open(sf, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                try:
                    comp = float(row["completeness"])
                    cont = float(row["contamination"])
                    gc   = float(row["GC"])
                    n50  = int(row["N50"])
                    size = int(row["size"])
                except (KeyError, ValueError):
                    continue

                bin_name = row["bin"].strip()
                binner_name = row.get("binner", file_stem).strip()

                bin_dir = work_dir / file_stem
                fa_path = bin_dir / "{}.fa".format(bin_name)
                if not fa_path.exists():
                    continue

                bin_id = "{}::{}".format(binner_name, bin_name)
                safe_name = "{}__{}".format(binner_name, bin_name)

                contigs_dict = {}
                current_name = ""
                current_len = 0
                
                try:
                    with open(fa_path, encoding="utf-8", errors="replace") as fa_fh:
                        for line in fa_fh:
                            if line.startswith(">"):
                                if current_name:
                                    contigs_dict[current_name] = current_len
                                current_name = line[1:].split()[0].strip()
                                current_len = 0
                            else:
                                current_len += len(line.strip())
                        if current_name:
                            contigs_dict[current_name] = current_len
                except OSError:
                    continue

                if not contigs_dict:
                    continue

                bin_contigs[bin_id] = contigs_dict
                total_calculated_size = sum(contigs_dict.values())

                entry = {
                    "bin_id":        bin_id,
                    "bin_name":      bin_name,
                    "safe_name":     safe_name,
                    "binner":        binner_name,
                    "completeness":  comp,
                    "contamination": cont,
                    "gc":            gc,
                    "lineage":       row.get("lineage", "root").strip(),
                    "n50":           n50,
                    "size":          total_calculated_size,
                    "score":         checkm1_score(comp, cont),
                    "fa_path":       fa_path,
                }
                rows.append(entry)

    print("Successfully mapped profile matrices for {} bins.".format(len(rows)), file=sys.stderr)
    return rows, bin_contigs

# ---------------------------------------------------------------------------
# Unified Multi-Layered Dereplication Engine
# ---------------------------------------------------------------------------

def deregplicate_metawrap_style(candidates, bin_contigs, overlap_threshold, ani_threshold):
    if not bin_contigs:
        return candidates

    # Global sort by composite CheckM score descending
    ranked_bins = sorted(candidates, key=lambda x: x["score"], reverse=True)
    accepted_bins = []
    
    # 1. Uniform Rule Step: Pre-compute K-mer profiles for all qualifying bins up front
    print("Pre-sketching K-mer sequences for uniform similarity checking...", file=sys.stderr)
    bin_kmers = {cand["bin_id"]: get_fast_kmers(cand["fa_path"]) for cand in ranked_bins}

    print("Evaluating uniform dereplication rules across groups...", file=sys.stderr)

    for cand in ranked_bins:
        cand_id = cand["bin_id"]
        cand_map = bin_contigs[cand_id]
        cand_total_bp = cand["size"]

        is_redundant = False
        
        for accepted in accepted_bins:
            acc_id = accepted["bin_id"]
            acc_map = bin_contigs[acc_id]
            
            # --- Rule A: Traditional Contig ID Overlap ---
            shared_bp = 0
            for contig, length in cand_map.items():
                if contig in acc_map:
                    shared_bp += length

            ratio_cand = (100.0 * shared_bp) / cand_total_bp
            ratio_acc  = (100.0 * shared_bp) / accepted["size"]
            max_overlap = max(ratio_cand, ratio_acc)

            if max_overlap >= overlap_threshold:
                is_redundant = True
                print("  [DROP] {} overlaps significantly with accepted bin {} ({:.1f}% contig overlap)".format(
                    cand_id, acc_id, max_overlap), file=sys.stderr)
                break

            # --- Rule B: Uniform K-mer Sequence Similarity ---
            kmers_cand = bin_kmers[cand_id]
            kmers_acc = bin_kmers[acc_id]
            
            if kmers_cand and kmers_acc:
                intersection = len(kmers_cand.intersection(kmers_acc))
                union = len(kmers_cand.union(kmers_acc))
                jaccard = intersection / union if union > 0 else 0
                
                if (jaccard * 100.0) >= ani_threshold:
                    is_redundant = True
                    print("  [DROP] {} is a genomic duplicate of accepted bin {} ({:.1f}% uniform sequence similarity)".format(
                        cand_id, acc_id, jaccard * 100.0), file=sys.stderr)
                    break

        if not is_redundant:
            accepted_bins.append(cand)

    return accepted_bins

# ---------------------------------------------------------------------------
# Execution Pipeline
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir",       required=True,  type=Path)
    parser.add_argument("--output-dir",     required=True,  type=Path)
    parser.add_argument("--stats-out",      required=True,  type=Path)
    parser.add_argument("--completeness",   default=50.0,   type=float)
    parser.add_argument("--contamination",  default=10.0,   type=float)
    parser.add_argument("--contig-overlap", default=80.0,   type=float, help="Min exact contig length match % (default: 80.0)")
    parser.add_argument("--ani-threshold",  default=75.0,   type=float, help="Min uniform K-mer Jaccard similarity % (default: 75.0)")
    parser.add_argument("--dry-run",        action="store_true")
    args = parser.parse_args()

    # -- Step 1: Load stats & parse structures --
    all_bins, bin_contigs = load_stats_and_lengths(args.work_dir)

    # -- Step 2: Lenient Pre-Filter --
    lenient_comp = args.completeness * 0.5
    lenient_cont = args.contamination * 3.0
    candidates = [
        b for b in all_bins
        if b["completeness"] >= lenient_comp and b["contamination"] <= lenient_cont
    ]
    print("  Candidates surviving lenient filtering pass: {}".format(len(candidates)), file=sys.stderr)

    if not candidates:
        sys.exit("ERROR: No bins met baseline filtering conditions.")

    # -- Step 3: Global Uniform Dereplication (Contigs + K-mers) --
    winners = deregplicate_metawrap_style(candidates, bin_contigs, args.contig_overlap, args.ani_threshold)

    # -- Step 4: Strict Quality Criteria Cut-offs --
    print("\nApplying final strict boundaries (comp >= {}%, cont <= {}%)...".format(
        args.completeness, args.contamination), file=sys.stderr)
    
    final_output_set = [
        w for w in winners 
        if w["completeness"] >= args.completeness and w["contamination"] <= args.contamination
    ]
    print("  Final selection count: {} bins pass verification.".format(len(final_output_set)), file=sys.stderr)

    if not final_output_set:
        sys.exit("ERROR: No bins passed strict evaluation boundaries.")

    # -- Step 5: Output Assembly & File Management --
    if not args.dry_run:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        
        with open(args.stats_out, "w", newline="", encoding="utf-8") as sf_out:
            writer = csv.writer(sf_out, delimiter="\t")
            writer.writerow(["bin", "completeness", "contamination", "GC", "lineage", "N50", "size", "binner"])
            for w in sorted(final_output_set, key=lambda x: x["score"], reverse=True):
                writer.writerow([
                    w["safe_name"],
                    "{:.4f}".format(w["completeness"]),
                    "{:.4f}".format(w["contamination"]),
                    "{:.3f}".format(w["gc"]),
                    w["lineage"],
                    w["n50"],
                    w["size"],
                    w["binner"]
                ])
        print("Wrote consolidated stats results out to: {}".format(args.stats_out), file=sys.stderr)

        for w in final_output_set:
            dest = args.output_dir / "{}.fa".format(w["safe_name"])
            shutil.copy2(w["fa_path"], dest)
        print("Completed moving files into output directory: {}".format(args.output_dir), file=sys.stderr)
    else:
        print("\n[DRY RUN RUNNING] Execution completed without disk modification flags.", file=sys.stderr)

if __name__ == "__main__":
    main()