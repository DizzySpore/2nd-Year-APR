#!/bin/bash
################################################################################
# STAGE 13B: QUAST Bakta - Annotation quality assessment
################################################################################
#
# Input variables (from SLURM template):
#   - SAMPLE          : Sample name
#   - NUM_THREADS     : Number of CPU threads
#   - BAKTA_OUT       : Directory containing Bakta annotation outputs
#   - BAKTA_QUAST_OUT : Output directory for QUAST results
#   - STATUS_DIR      : Path to status directory
#
# Output:
#   - QUAST report      : ${BAKTA_QUAST_OUT}/${SAMPLE}_quast_bakta/report.html
#   - Completion marker : ${BAKTA_QUAST_OUT}/quast_loop_completed.txt
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_14
#
# Note:
#   - QUAST is run on FFN nucleotide sequence files produced by Bakta
#   - If no FFN files are found the stage completes with a warning rather
#     than failing, since an empty bin set is a valid outcome
#
################################################################################
echo "$(date): [STAGE 13B/24] Starting QUAST on Bakta annotations for ${SAMPLE}..."
echo "  Input  : ${BAKTA_OUT}/*/*.ffn"
echo "  Output : ${BAKTA_QUAST_OUT}"
# Collect all FFN files produced by Bakta
FFN_FILES=()
for bakta_dir in "${BAKTA_OUT}"/*/; do
    [ -d "$bakta_dir" ] || continue
    BIN_NAME=$(basename "$bakta_dir")
    FFN_FILE="${bakta_dir}/${BIN_NAME}.ffn"
    [ -f "$FFN_FILE" ] && FFN_FILES+=("$FFN_FILE")
done
if [ "${#FFN_FILES[@]}" -eq 0 ]; then
    echo "$(date): WARNING - No Bakta FFN files found for ${SAMPLE}. Skipping QUAST."
else
    conda activate quast-env
    quast.py \
        --threads ${NUM_THREADS} \
        --output-dir "${BAKTA_QUAST_OUT}/${SAMPLE}_quast_bakta" \
        --silent \
        "${FFN_FILES[@]}"
    conda deactivate || true
    echo "$(date): QUAST assessment complete for ${#FFN_FILES[@]} annotated bins."
fi
# Write completion marker used by checkpoint scanner
touch "${BAKTA_QUAST_OUT}/quast_loop_completed.txt"
echo "STAGE_14" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 13B/24] QUAST Bakta assessment complete."