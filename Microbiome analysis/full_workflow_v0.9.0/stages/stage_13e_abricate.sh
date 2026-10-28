#!/bin/bash
################################################################################
# STAGE 13E: abricate - AMR/virulence gene screening
################################################################################
#
# Input variables (from SLURM template):
#   - SAMPLE          : Sample name
#   - NUM_THREADS     : Number of CPU threads
#   - FINAL_BINS_DIR  : Directory containing unfiltered nucleotide bins from STAGE_07E
#   - ABRICATE_OUT    : Output directory for abricate results
#   - STATUS_DIR      : Path to status directory
#
# Output:
#   - Per-bin, per-database reports : ${ABRICATE_OUT}/${BIN_NAME}/${BIN_NAME}_${DB}.tab
#   - Per-database sample summaries : ${ABRICATE_OUT}/${SAMPLE}_${DB}_summary.tab
#   - Completion marker             : ${ABRICATE_OUT}/loop_completed.txt
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_14
#
# Note:
#   - Unlike dbCAN/eggNOG-mapper (STAGE_13C/13D), abricate works directly on
#     assembled nucleotide contigs rather than Bakta's predicted proteins, so
#     it runs on the raw bin fasta files from FINAL_BINS_DIR (same input as
#     Bakta in STAGE_13).
#   - Runs against the core AMR/virulence database set: card, vfdb, resfinder,
#     plasmidfinder. Add/remove entries in DBS below to change the set.
#   - loop_completed.txt is written only after all bins/databases have been
#     processed and is used as the checkpoint for this stage.
#
################################################################################
echo "$(date): [STAGE 13E/24] Starting abricate AMR/virulence screening for ${SAMPLE}..."
echo "  Input bins : ${FINAL_BINS_DIR}/*.fa"
echo "  Output     : ${ABRICATE_OUT}"

DBS=(card vfdb resfinder plasmidfinder)
echo "  Databases  : ${DBS[*]}"

conda activate abricate-env
BIN_COUNT=0
for bin_file in "${FINAL_BINS_DIR}/"*.fa; do
    [ -f "$bin_file" ] || continue
    BIN_NAME=$(basename "$bin_file" .fa)
    BIN_COUNT=$(( BIN_COUNT + 1 ))
    echo "$(date): Screening ${BIN_NAME} (bin ${BIN_COUNT})..."

    mkdir -p "${ABRICATE_OUT}/${BIN_NAME}"

    for db in "${DBS[@]}"; do
        OUT_FILE="${ABRICATE_OUT}/${BIN_NAME}/${BIN_NAME}_${db}.tab"
        abricate \
            --db "$db" \
            --threads ${NUM_THREADS} \
            --nopath \
            "$bin_file" \
            > "$OUT_FILE"

        HIT_COUNT=$(($(wc -l < "$OUT_FILE") - 1))
        [ "$HIT_COUNT" -lt 0 ] && HIT_COUNT=0
        echo "$(date): ${BIN_NAME} - ${db}: ${HIT_COUNT} hits"
    done
done

# -- Per-database, per-sample summaries across all bins -------------------------
if [ "$BIN_COUNT" -gt 0 ]; then
    for db in "${DBS[@]}"; do
        SUMMARY_INPUTS=("${ABRICATE_OUT}/"*/"*_${db}.tab")
        if [ -e "${SUMMARY_INPUTS[0]}" ]; then
            abricate --summary "${SUMMARY_INPUTS[@]}" > "${ABRICATE_OUT}/${SAMPLE}_${db}_summary.tab"
            echo "$(date): Wrote ${db} summary for ${SAMPLE} (${#SUMMARY_INPUTS[@]} bins)"
        fi
    done
fi
conda deactivate || true

if [ "$BIN_COUNT" -eq 0 ]; then
    echo "$(date): WARNING - No bins found in ${FINAL_BINS_DIR} for ${SAMPLE}."
fi
echo "$(date): abricate screening complete for ${BIN_COUNT} bins."
# Write completion marker used by checkpoint scanner
touch "${ABRICATE_OUT}/loop_completed.txt"
echo "STAGE_14" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 13E/24] abricate screening complete."