#!/bin/bash
################################################################################
# STAGE 13: Bakta - Prokaryotic gene annotation
################################################################################
#
# Input variables (from SLURM template):
#   - SAMPLE        : Sample name
#   - NUM_THREADS   : Number of CPU threads
#   - FINAL_BINS_DIR: Directory containing unfiltered bins from STAGE_07E
#   - BAKTA_OUT     : Output directory for Bakta annotations
#   - DATABASES     : Path to reference databases
#   - STATUS_DIR    : Path to status directory
#
# Output:
#   - Per-bin annotation directories : ${BAKTA_OUT}/${BIN_NAME}/
#   - Completion marker              : ${BAKTA_OUT}/loop_completed.txt
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_13B
#
# Note:
#   - Bakta now runs on UNFILTERED bins straight from STAGE_07E, since bin
#     filtering (STAGE_14) runs after this stage in the current pipeline order.
#   - Bakta is run on each bin independently
#   - loop_completed.txt is written only after all bins have been
#     successfully annotated and is used as the checkpoint for this stage
#
################################################################################
echo "$(date): [STAGE 13/24] Starting Bakta annotation for ${SAMPLE}..."
echo "  Input bins : ${FINAL_BINS_DIR}/*.fa"
echo "  Output     : ${BAKTA_OUT}"
echo "  Database   : ${DATABASES}/bakta-db/db"
conda activate bakta-env
BIN_COUNT=0
for bin_file in "${FINAL_BINS_DIR}/"*.fa; do
    [ -f "$bin_file" ] || continue
    BIN_NAME=$(basename "$bin_file" .fa)
    BIN_COUNT=$(( BIN_COUNT + 1 ))
    echo "$(date): Annotating ${BIN_NAME} (bin ${BIN_COUNT})..."
    
    bakta \
        "$bin_file" \
        --meta \
        --output "${BAKTA_OUT}/${BIN_NAME}" \
        --prefix "${BIN_NAME}" \
        --threads ${NUM_THREADS} \
        --db "${DATABASES}/bakta-db/db" \
        --verbose
        
    # Log predicted protein count from GFF3 output
    GFF_FILE="${BAKTA_OUT}/${BIN_NAME}/${BIN_NAME}.gff3"
    if [ -f "$GFF_FILE" ]; then
        PROTEIN_COUNT=$(grep -c "product=" "$GFF_FILE" || true)
        echo "$(date): ${BIN_NAME} - ${PROTEIN_COUNT} predicted proteins"
    fi
done
conda deactivate || true
if [ "$BIN_COUNT" -eq 0 ]; then
    echo "$(date): WARNING - No bins found in ${FINAL_BINS_DIR} for ${SAMPLE}."
fi
echo "$(date): Bakta annotation complete for ${BIN_COUNT} bins."
# Write completion marker used by checkpoint scanner
touch "${BAKTA_OUT}/loop_completed.txt"
echo "STAGE_13B" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 13/24] Bakta annotation complete."