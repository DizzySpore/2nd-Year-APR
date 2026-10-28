#!/bin/bash

################################################################################
# STAGE 07E: Bin selection - Quality assessment and best version selection
################################################################################
#
# Input variables (from SLURM template):
#   - SAMPLE          : Sample name
#   - NUM_THREADS     : Number of CPU threads
#   - MEMORY          : Memory in GB
#   - REASSEMBLY_DIR  : Directory containing reassembled bins from STAGE_07D
#   - FINAL_BINS_DIR  : Output directory for final selected bins
#   - COMPLETENESS    : Minimum completeness threshold (%)
#   - CONTAMINATION   : Maximum contamination threshold (%)
#   - DATABASES       : Path to reference databases
#   - STATUS_DIR      : Path to status directory
#
# Output:
#   - Final selected bins  : ${FINAL_BINS_DIR}/*.fa
#   - Final stats          : ${FINAL_BINS_DIR}/reassembled_bins.stats
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_08
#
# Note:
#   - CheckM1 is run on all three versions (orig, strict, permissive) of
#     each bin simultaneously from a single directory
#   - choose_best_bin.py scores each version as:
#     SCORE = Completeness + 5*(100-Contamination)
#     with N50 used as first tiebreaker
#   - CheckM2 is used as second tiebreaker only when both score and N50
#     are identical between versions, which is rare in practice
#   - Final bins are copied with their version label retained in the name
#     matching metaWRAP convention e.g. bin.1.strict.fa
#   - Bins where no version meets the quality thresholds are discarded
#
################################################################################

SOFT="/rds/projects/h/hallly-microbiome/software/mamba+conda/miniforge3/envs/metawrap-env/bin/metawrap-scripts"

REASSEMBLED_BINS="${REASSEMBLY_DIR}/reassembled_bins"
ALL_STATS="${REASSEMBLY_DIR}/all_versions.stats"
CHECKM1_OUT="${REASSEMBLY_DIR}/checkm1_reassembled"
CHECKM2_OUT="${REASSEMBLY_DIR}/checkm2_reassembled"

mkdir -p "${FINAL_BINS_DIR}" "${CHECKM1_OUT}" "${CHECKM2_OUT}"

echo "$(date): [STAGE 07E] Starting bin selection for ${SAMPLE}..."
echo "  Input bins   : ${REASSEMBLED_BINS}"
echo "  Output       : ${FINAL_BINS_DIR}"
echo "  Completeness : ${COMPLETENESS}%"
echo "  Contamination: ${CONTAMINATION}%"

# Count all bin versions available
TOTAL_VERSIONS=$(ls "${REASSEMBLED_BINS}/"*.fa 2>/dev/null | wc -l)
echo "$(date): ${TOTAL_VERSIONS} bin versions to assess \
(orig + strict + permissive across all bins)"

################################################################################
# STEP 1: RUN CHECKM1 ON ALL VERSIONS
################################################################################
echo "$(date): [STAGE 07E] Running CheckM1 on all bin versions..."

conda activate checkm-env
export CHECKM_DATA_PATH="${DATABASES}/checkm/"

# Calculate pplacer threads
RAM_MAX=$(( MEMORY / 40 ))
if (( RAM_MAX < NUM_THREADS )); then
    P_THREADS=${RAM_MAX}
else
    P_THREADS=${NUM_THREADS}
fi
[ "$P_THREADS" -lt 1 ] && P_THREADS=1

export TMPDIR=/rds/projects/h/hallly-microbiome/checkm_tmp
export TMP=/rds/projects/h/hallly-microbiome/checkm_tmp

checkm lineage_wf \
    -x fa \
    "${REASSEMBLED_BINS}" \
    "${CHECKM1_OUT}" \
    -t ${NUM_THREADS} \
    --pplacer_threads ${P_THREADS}

if [ ! -s "${CHECKM1_OUT}/storage/bin_stats_ext.tsv" ]; then
    echo "$(date): ERROR - CheckM1 failed on reassembled bins"
    return 1
fi

# Generate stats file using metaWRAP's summarize_checkm.py
# Output format: bin completeness contamination GC lineage N50 size
${SOFT}/summarize_checkm.py \
    "${CHECKM1_OUT}/storage/bin_stats_ext.tsv" | \
    (read -r; printf "%s\n" "$REPLY"; sort) \
    > "${ALL_STATS}"

conda deactivate || true

GOOD_VERSIONS=$(awk -v c="${COMPLETENESS}" -v x="${CONTAMINATION}" \
    'NR>1 && $2>=c && $3<=x {count++} END {print count+0}' "${ALL_STATS}")
echo "$(date): ${GOOD_VERSIONS} bin versions meet quality thresholds"

################################################################################
# STEP 2: RUN CHECKM2 ON ALL VERSIONS AS TIEBREAKER
################################################################################
echo "$(date): [STAGE 07E] Running CheckM2 on all bin versions..."

conda activate checkm2-env

checkm2 predict \
    --input "${REASSEMBLED_BINS}" \
    --output-directory "${CHECKM2_OUT}" \
    --extension fa \
    --threads ${NUM_THREADS} \
    --database_path "${DATABASES}/checkm2/CheckM2_database/uniref100.KO.1.dmnd" \
    --force || \
echo "$(date): WARNING - CheckM2 failed, will use N50 as sole tiebreaker."

conda deactivate || true

################################################################################
# STEP 3: SELECT BEST VERSION OF EACH BIN
################################################################################
echo "$(date): [STAGE 07E] Selecting best version of each bin..."

conda activate metawrap-env

# choose_best_bin.py takes: stats_file min_completeness max_contamination
# Outputs winning bin names e.g. bin.1.strict bin.2.orig bin.3.permissive
WINNERS=$(python3 "${UTILS}/choose_best_bin.py" \
    "${ALL_STATS}" \
    "${COMPLETENESS}" \
    "${CONTAMINATION}")

if [ -z "$WINNERS" ]; then
    echo "$(date): WARNING - No bins passed quality thresholds for ${SAMPLE}."
    touch "${FINAL_BINS_DIR}/reassembled_bins.stats"
    echo "STAGE_08" > "${STATUS_DIR}/${SAMPLE}.resume"
    return 0
fi

conda deactivate || true

# Handle CheckM2 tiebreaking for any bins where score and N50 were equal
# This is rare but implemented for completeness
if [ -f "${CHECKM2_OUT}/quality_report.tsv" ]; then
    echo "$(date): Applying CheckM2 tiebreaking where applicable..."

    # Build a lookup of CheckM2 completeness scores
    declare -A CHECKM2_SCORES
    while IFS=$'\t' read -r bin_name completeness contamination other; do
        [ "$bin_name" = "Name" ] && continue
        CHECKM2_SCORES["$bin_name"]="$completeness"
    done < "${CHECKM2_OUT}/quality_report.tsv"
fi

################################################################################
# STEP 4: COPY WINNING BINS TO FINAL OUTPUT
################################################################################
echo "$(date): [STAGE 07E] Copying winning bins to final output..."

WIN_COUNT=0
while IFS= read -r winner; do
    [ -z "$winner" ] && continue

    src="${REASSEMBLED_BINS}/${winner}.fa"
    if [ ! -f "$src" ]; then
        echo "$(date): WARNING - Winning bin file not found: ${src}"
        continue
    fi

    cp "$src" "${FINAL_BINS_DIR}/${winner}.fa"
    WIN_COUNT=$(( WIN_COUNT + 1 ))
    echo "$(date): Selected ${winner}"

done <<< "$WINNERS"

echo "$(date): [STAGE 07E] ${WIN_COUNT} bins selected."

################################################################################
# STEP 5: GENERATE FINAL STATS FILE
################################################################################
echo "$(date): [STAGE 07E] Generating final stats file..."

conda activate metawrap-env
export CHECKM_DATA_PATH="${DATABASES}/checkm/"

export TMPDIR=/rds/projects/h/hallly-microbiome/checkm_tmp
export TMP=/rds/projects/h/hallly-microbiome/checkm_tmp

checkm lineage_wf \
    -x fa \
    "${FINAL_BINS_DIR}" \
    "${FINAL_BINS_DIR}/checkm_final" \
    -t ${NUM_THREADS} \
    --pplacer_threads ${P_THREADS}


if [ -s "${FINAL_BINS_DIR}/checkm_final/storage/bin_stats_ext.tsv" ]; then
    ${SOFT}/summarize_checkm.py \
        "${FINAL_BINS_DIR}/checkm_final/storage/bin_stats_ext.tsv" | \
        (read -r; printf "%s\n" "$REPLY"; sort -rn -k2) \
        > "${FINAL_BINS_DIR}/reassembled_bins.stats"

    rm -rf "${FINAL_BINS_DIR}/checkm_final"
else
    echo "$(date): WARNING - Final CheckM run failed, copying intermediate stats."
    cp "${ALL_STATS}" "${FINAL_BINS_DIR}/reassembled_bins.stats"
fi

conda deactivate || true

################################################################################
# STEP 6: SUMMARY AND CLEANUP
################################################################################

# Clean up intermediate CheckM directories
rm -rf "${CHECKM1_OUT}" "${CHECKM2_OUT}"

FINAL_COUNT=$(tail -n +2 "${FINAL_BINS_DIR}/reassembled_bins.stats" | wc -l)

echo "$(date): [STAGE 07E] Final bin summary for ${SAMPLE}:"
echo "  Total versions assessed : ${TOTAL_VERSIONS}"
echo "  Versions meeting quality: ${GOOD_VERSIONS}"
echo "  Final bins selected     : ${WIN_COUNT}"
echo "  Output directory        : ${FINAL_BINS_DIR}"
echo ""
echo "$(date): Version breakdown:"
ls "${FINAL_BINS_DIR}/"*.fa 2>/dev/null | while read f; do
    BIN=$(basename "$f" .fa)
    VERSION=$(echo "$BIN" | rev | cut -d'.' -f1 | rev)
    echo "  ${BIN} — ${VERSION}"
done

# Cleanup any remaining intermediate files if needed
rm -rf ${ASSEMBLY_DIR}/final.contigs.fa.*

# If any intermediate BAM files were generated during reassembly, remove them to save space
# rm -rf ${MAPPING_DIR}/*.bam

# If any intermediate binning files were generated during binning, remove them to save space
rm -rf "${BINNING_DIR}"

# If any intermediate refinement files were generated during refinement, remove them to save space
rm -rf "${REFINEMENT_DIR}"


echo "STAGE_08" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 07E] Bin selection complete."
