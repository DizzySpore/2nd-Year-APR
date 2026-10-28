#!/bin/bash

################################################################################
# STAGE 19: antiSMASH - Secondary metabolite detection
################################################################################
#
# Input variables (from SLURM template):
#   - SAMPLE        : Sample name
#   - NUM_THREADS   : Number of CPU threads
#   - BAKTA_OUT     : Directory containing Bakta annotation outputs (GBK files)
#   - ANTIMASH_OUT  : Output directory for antiSMASH results
#   - STATUS_DIR    : Path to status directory
#
# Output:
#   - Per-bin antiSMASH results : ${ANTIMASH_OUT}/${BIN_NAME}_antismash/
#   - Completion marker         : ${ANTIMASH_OUT}/antismash_loop_completed.txt
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_20
#
# Note:
#   - antiSMASH is run on GBK files produced by Bakta in STAGE_17
#   - Failures on individual bins are warnings not errors since not all
#     bins will contain biosynthetic gene clusters
#   - loop_completed.txt is written only after all bins have been
#     processed and is used as the checkpoint for this stage
#
################################################################################

echo "$(date): [STAGE 19/24] Starting antiSMASH secondary metabolite detection for ${SAMPLE}..."
echo "  Input  : ${BAKTA_OUT}/*/*.gbk"
echo "  Output : ${ANTISMASH_OUT}"

conda activate antismash-env

GENOME_COUNT=0
for bakta_dir in "${BAKTA_OUT}"/*/; do
    [ -d "$bakta_dir" ] || continue
    BIN_NAME=$(basename "$bakta_dir")
    GBK_FILE="${bakta_dir}/${BIN_NAME}.gbff"

    if [ ! -f "$GBK_FILE" ]; then
        echo "$(date): WARNING - GBK file not found for ${BIN_NAME}, skipping."
        continue
    fi

    GENOME_COUNT=$(( GENOME_COUNT + 1 ))
    echo "$(date): Running antiSMASH on ${BIN_NAME} (genome ${GENOME_COUNT})..."


    antismash \
        "$GBK_FILE" \
        --output-dir "${ANTISMASH_OUT}/${BIN_NAME}_antismash" \
        --output-basename "${BIN_NAME}" \
        --genefinding-tool prodigal \
        --cb-general \
        --cb-subclusters \
        --cb-knownclusters \
        --smcog-trees \
        --asf \
        --rre \
        --pfam2go \
        --clusterhmmer \
        --taxon bacteria \
        --cpus ${NUM_THREADS} \
        --allow-long-headers || \
    echo "$(date): WARNING - antiSMASH failed for ${BIN_NAME}, continuing."

    # Log detected cluster count from JSON output if available
    JSON_FILE="${ANTIMASH_OUT}/${BIN_NAME}_antismash/${BIN_NAME}.json"
    if [ -f "$JSON_FILE" ]; then
        CLUSTER_COUNT=$(grep -o '"type": "cluster"' "$JSON_FILE" | wc -l)
        echo "$(date): ${BIN_NAME} - ${CLUSTER_COUNT} biosynthetic gene clusters detected."
    fi

done

conda deactivate || true

if [ "$GENOME_COUNT" -eq 0 ]; then
    echo "$(date): WARNING - No Bakta GBK files found for ${SAMPLE}."
fi

echo "$(date): antiSMASH complete for ${GENOME_COUNT} genomes."

# Write completion marker used by checkpoint scanner
touch "${ANTISMASH_OUT}/antismash_loop_completed.txt"

echo "STAGE_20" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 19/24] antiSMASH complete."