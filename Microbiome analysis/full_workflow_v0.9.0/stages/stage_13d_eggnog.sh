#!/bin/bash
################################################################################
# STAGE 13D: eggNOG-mapper - Functional/orthology annotation
################################################################################
#
# Input variables (from SLURM template):
#   - SAMPLE        : Sample name
#   - NUM_THREADS   : Number of CPU threads
#   - BAKTA_OUT     : Directory containing per-bin Bakta annotations from STAGE_13
#   - EGGNOG_OUT    : Output directory for eggNOG-mapper annotations
#   - DATABASES     : Path to reference databases
#   - STATUS_DIR    : Path to status directory
#
# Output:
#   - Per-bin eggNOG directories : ${EGGNOG_OUT}/${BIN_NAME}/
#   - Completion marker          : ${EGGNOG_OUT}/loop_completed.txt
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_13E
#
# Note:
#   - Runs on the Bakta-predicted protein sequences (${BIN_NAME}.faa) via
#     --itype proteins, so no gene calling is repeated here.
#   - --override is used because STAGE_13D always (re)runs from scratch for
#     custom-workflow resubmissions; loop_completed.txt is the checkpoint.
#
################################################################################
echo "$(date): [STAGE 13D/24] Starting eggNOG-mapper annotation for ${SAMPLE}..."
echo "  Input proteins : ${BAKTA_OUT}/*/*.faa"
echo "  Output         : ${EGGNOG_OUT}"
echo "  Database       : ${DATABASES}/eggnog-mapper-data"
conda activate eggnog-mapper-env
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

    mkdir -p "${EGGNOG_OUT}/${BIN_NAME}"

    emapper.py \
        -i "$FAA_FILE" \
        --itype proteins \
        -m diamond \
        --cpu ${NUM_THREADS} \
        --data_dir "${DATABASES}/eggnog-mapper-data" \
        --output "${BIN_NAME}" \
        --output_dir "${EGGNOG_OUT}/${BIN_NAME}" \
        --override

    # Log annotated gene count from annotations output
    ANNOT_FILE="${EGGNOG_OUT}/${BIN_NAME}/${BIN_NAME}.emapper.annotations"
    if [ -f "$ANNOT_FILE" ]; then
        ANNOT_COUNT=$(grep -vc "^#" "$ANNOT_FILE" || true)
        echo "$(date): ${BIN_NAME} - ${ANNOT_COUNT} annotated genes"
    fi
done
conda deactivate || true
if [ "$BIN_COUNT" -eq 0 ]; then
    echo "$(date): WARNING - No Bakta protein files found under ${BAKTA_OUT} for ${SAMPLE}."
fi
echo "$(date): eggNOG-mapper annotation complete for ${BIN_COUNT} bins."
# Write completion marker used by checkpoint scanner
touch "${EGGNOG_OUT}/loop_completed.txt"
echo "STAGE_13E" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 13D/24] eggNOG-mapper annotation complete."