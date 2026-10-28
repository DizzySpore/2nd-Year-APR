#!/bin/bash
################################################################################
# STAGE 14: Bin filtering - Quality filter high-confidence bins
################################################################################
#
# Input variables (from SLURM template):
#   - SAMPLE        : Sample name
#   - FINAL_BINS_DIR: Directory containing taxonomy-renamed bins (STAGE_08 output)
#   - BAKTA_OUT     : Directory containing per-bin Bakta annotations (STAGE_13 output)
#   - GOOD_BINS_DIR : Output directory for filtered high-quality bins
#   - STATUS_DIR    : Path to status directory
#
# Output:
#   - Filtered bins : ${GOOD_BINS_DIR}/*.fa
#   - Quality log   : ${GOOD_BINS_DIR}/${SAMPLE}_bin_quality_summary.txt
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_15
#
# Quality thresholds for a "good" bin:
#   - Completeness  >= 95%
#   - Contamination <= 5%
#   - At least one tRNA gene detected by Bakta
#   - At least one rRNA gene detected by Bakta
#
#   These are deliberately hardcoded (MIN_COMPLETENESS / MAX_CONTAMINATION
#   below) rather than reusing the global COMPLETENESS/CONTAMINATION
#   variables from 00_config.sh (50%/10%), since those are used earlier in
#   the pipeline for coarser bin refinement/selection and are a different,
#   more permissive bar than the "good bin" definition used here.
#
# rRNA/tRNA criterion:
#   - This checks for PRESENCE only (>=1 tRNA and >=1 rRNA gene called by
#     Bakta), not full MIMAG "high-quality genome" criteria (e.g. >=18
#     tRNAs and 5S/16S/23S rRNA specifically). Tighten the checks below if
#     you want the stricter MIMAG-style definition instead.
#
# Notes:
#   - Reads completeness/contamination from
#     ${FINAL_BINS_DIR}/reassembled_bins.stats. Bin names in that file must
#     match the taxonomy-renamed bins written by STAGE_08
#     (${SAMPLE}_<bin>_<species>.fa) - this requires the companion fix to
#     stage_08_gtdbtk.sh that actually wires up bin_rename_map /
#     update_stats_file (previously defined but never called/populated).
#   - rRNA/tRNA presence is read from the Bakta GFF3
#     (${BAKTA_OUT}/<bin_name>/<bin_name>.gff3), column 3 feature type.
#   - Bins missing a stats entry, a bin file, or a Bakta GFF3 are logged as
#     FAIL with a reason and skipped (missing GFF3 = treated as no
#     tRNA/rRNA rather than aborting the run).
#
################################################################################

echo "$(date): [STAGE 14/24] Starting bin quality filtering for ${SAMPLE}..."
echo "  Input     : ${FINAL_BINS_DIR}"
echo "  Bakta     : ${BAKTA_OUT}"
echo "  Output    : ${GOOD_BINS_DIR}"
echo "  Thresholds: completeness >= 95%, contamination <= 5%, tRNA present, rRNA present"

MIN_COMPLETENESS=95
MAX_CONTAMINATION=5

STATS_FILE="${FINAL_BINS_DIR}/reassembled_bins.stats"
QUALITY_LOG="${GOOD_BINS_DIR}/${SAMPLE}_bin_quality_summary.txt"

if [ ! -f "${STATS_FILE}" ]; then
    echo "$(date): ERROR - Stats file not found: ${STATS_FILE}" >&2
    exit 1
fi

echo "Bin_Name,Completeness,Contamination,tRNA_count,rRNA_count,Status,Reason" > "${QUALITY_LOG}"

BIN_COUNT=0
PASS_COUNT=0

while IFS=$'\t' read -r bin_name completeness contamination other; do
    [ "$bin_name" = "bin" ] && continue
    [ -z "$bin_name" ] && continue

    bin_file="${FINAL_BINS_DIR}/${bin_name}.fa"
    if [ ! -f "$bin_file" ]; then
        echo "$(date): WARNING - Bin file not found, skipping: ${bin_file}"
        echo "${bin_name},NA,NA,NA,NA,FAIL,bin_file_missing" >> "${QUALITY_LOG}"
        continue
    fi

    BIN_COUNT=$(( BIN_COUNT + 1 ))

    completeness_int=$(echo "$completeness" | awk '{print int($1)}')
    contamination_int=$(echo "$contamination" | awk '{print int($1)}')

    # -- rRNA/tRNA presence from Bakta ------------------------------------------
    gff_file="${BAKTA_OUT}/${bin_name}/${bin_name}.gff3"
    trna_count=0
    rrna_count=0
    if [ -f "$gff_file" ]; then
        trna_count=$(awk -F'\t' '$3=="tRNA"' "$gff_file" | wc -l)
        rrna_count=$(awk -F'\t' '$3=="rRNA"' "$gff_file" | wc -l)
    else
        echo "$(date): WARNING - Bakta GFF3 not found for ${bin_name}, treating tRNA/rRNA as absent: ${gff_file}"
    fi

    # -- Decision ----------------------------------------------------------------
    reason=""
    [ "$completeness_int" -lt "$MIN_COMPLETENESS" ] && reason="${reason}low_completeness;"
    [ "$contamination_int" -gt "$MAX_CONTAMINATION" ] && reason="${reason}high_contamination;"
    [ "$trna_count" -eq 0 ] && reason="${reason}no_tRNA;"
    [ "$rrna_count" -eq 0 ] && reason="${reason}no_rRNA;"

    if [ -z "$reason" ]; then
        echo "$(date): PASS - ${bin_name} (C:${completeness_int}% X:${contamination_int}% tRNA:${trna_count} rRNA:${rrna_count})"
        cp "$bin_file" "${GOOD_BINS_DIR}/${bin_name}.fa"
        echo "${bin_name},${completeness_int},${contamination_int},${trna_count},${rrna_count},PASS," >> "${QUALITY_LOG}"
        PASS_COUNT=$(( PASS_COUNT + 1 ))
    else
        echo "$(date): FAIL - ${bin_name} (C:${completeness_int}% X:${contamination_int}% tRNA:${trna_count} rRNA:${rrna_count}) [${reason%;}]"
        echo "${bin_name},${completeness_int},${contamination_int},${trna_count},${rrna_count},FAIL,${reason%;}" >> "${QUALITY_LOG}"
    fi

done < "$STATS_FILE"

echo "$(date): Filtering complete: ${PASS_COUNT}/${BIN_COUNT} bins passed thresholds."
if [ "${PASS_COUNT}" -eq 0 ]; then
    echo "$(date): WARNING - No bins passed quality thresholds for ${SAMPLE}."
fi

echo "STAGE_15" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 14/24] Bin filtering complete."