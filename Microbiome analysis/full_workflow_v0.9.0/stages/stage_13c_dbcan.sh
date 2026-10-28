#!/bin/bash
################################################################################
# STAGE 13C: dbCAN - CAZyme annotation
################################################################################
#
# Input variables (from SLURM template):
#   - SAMPLE        : Sample name
#   - NUM_THREADS   : Number of CPU threads
#   - BAKTA_OUT     : Directory containing per-bin Bakta annotations from STAGE_13
#   - DBCAN_OUT     : Output directory for dbCAN annotations
#   - DATABASES     : Path to reference databases
#   - STATUS_DIR    : Path to status directory
#
# Output:
#   - Per-bin dbCAN directories : ${DBCAN_OUT}/${BIN_NAME}/
#   - Completion marker         : ${DBCAN_OUT}/loop_completed.txt
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_13D
#
# Note:
#   - Runs on the Bakta-predicted protein sequences (${BIN_NAME}.faa) rather
#     than re-predicting genes from the raw bin nucleotide fasta, since gene
#     calling has already been done by Bakta in STAGE_13.
#   - inputFile is run in 'protein' mode (dbCAN's protein/prok/meta distinction
#     only matters for its own internal gene calling, which is skipped here).
#   - loop_completed.txt is written only after all bins have been processed
#     and is used as the checkpoint for this stage.
#
################################################################################
echo "$(date): [STAGE 13C/24] Starting dbCAN CAZyme annotation for ${SAMPLE}..."
echo "  Input proteins : ${BAKTA_OUT}/*/*.faa"
echo "  Output         : ${DBCAN_OUT}"
echo "  Database       : ${DATABASES}/db_v5-2_9-13-2025"
conda activate dbcan5-env
BIN_COUNT=0
for bin_dir in "${BAKTA_OUT}/"*/; do
    [ -d "$bin_dir" ] || continue
    BIN_NAME=$(basename "$bin_dir")
    FAA_FILE="${bin_dir}${BIN_NAME}.faa"

    if [ ! -f "$FAA_FILE" ]; then
        echo "$(date): WARNING - No Bakta protein file found for ${BIN_NAME} (${FAA_FILE}), skipping."
        continue
    fi

    BIN_COUNT=$(( BIN_COUNT + 1 ))
    echo "$(date): Annotating ${BIN_NAME} (bin ${BIN_COUNT})..."

    run_dbcan CAZyme_annotation \
        --input_raw_data "$FAA_FILE" \
        --mode protein \
        --output_dir "${DBCAN_OUT}/${BIN_NAME}" \
        --db_dir "${DATABASES}/db_v5-2_9-13-2025" \
        --threads ${NUM_THREADS}

    # Log predicted CAZyme count from overview output
    OVERVIEW_FILE="${DBCAN_OUT}/${BIN_NAME}/overview.txt"
    if [ -f "$OVERVIEW_FILE" ]; then
        CAZYME_COUNT=$(($(wc -l < "$OVERVIEW_FILE") - 1))
        echo "$(date): ${BIN_NAME} - ${CAZYME_COUNT} CAZyme hits"
    fi
done
conda deactivate || true
if [ "$BIN_COUNT" -eq 0 ]; then
    echo "$(date): WARNING - No Bakta protein files found under ${BAKTA_OUT} for ${SAMPLE}."
fi
echo "$(date): dbCAN annotation complete for ${BIN_COUNT} bins."
# Write completion marker used by checkpoint scanner
touch "${DBCAN_OUT}/loop_completed.txt"
echo "STAGE_13D" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 13C/24] dbCAN annotation complete."
