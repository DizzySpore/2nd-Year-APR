#!/bin/bash

################################################################################
# STAGE 06: QUAST - Assembly quality assessment
################################################################################
#
# Input variables (from SLURM template):
#   - SAMPLE             : Sample name
#   - NUM_THREADS        : Number of CPU threads
#   - ASSEMBLY_FILE      : Path to MEGAHIT final contigs file
#   - QUAST_ASSEMBLY_OUT : Output directory for QUAST results
#   - STATUS_DIR         : Path to status directory
#
# Output:
#   - QUAST report : ${QUAST_ASSEMBLY_OUT}/report.html
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_07
#
################################################################################

echo "$(date): [STAGE 06/24] Starting QUAST assembly assessment for ${SAMPLE}..."
echo "  Input  : ${ASSEMBLY_FILE}"
echo "  Output : ${QUAST_ASSEMBLY_OUT}"

conda activate quast-env

quast "${ASSEMBLY_FILE}" \
    -o "${QUAST_ASSEMBLY_OUT}" \
    --threads ${NUM_THREADS}

conda deactivate || true

echo "STAGE_07" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 06/24] QUAST assembly assessment complete."