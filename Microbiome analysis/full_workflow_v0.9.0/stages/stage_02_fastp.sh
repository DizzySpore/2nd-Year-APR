#!/bin/bash

################################################################################
# STAGE 02: fastp - Read trimming and quality filtering
################################################################################
#
# Input variables (from SLURM template):
#   - SAMPLE      : Sample name
#   - NUM_THREADS : Number of CPU threads
#   - READ1       : Path to R1 fastq file
#   - READ2       : Path to R2 fastq file
#   - QC_READS    : Output directory for trimmed reads and reports
#   - STATUS_DIR  : Path to status directory
#
# Output:
#   - Trimmed reads : ${QC_READS}/out.${SAMPLE}.R1.fq.gz
#                     ${QC_READS}/out.${SAMPLE}.R2.fq.gz
#   - QC reports   : ${QC_READS}/${SAMPLE}_fastp.html
#                    ${QC_READS}/${SAMPLE}_fastp.json
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_03
#
################################################################################

echo "$(date): [STAGE 02/24] Starting fastp for ${SAMPLE}..."

module purge
module load bluebear
module load bear-apps/2023a
module load fastp/0.24.0-GCC-12.3.0

fastp --thread ${NUM_THREADS} -q 20 \
    -i "$READ1" -I "$READ2" \
    -o "${QC_READS}/out.${SAMPLE}.R1.fq.gz" \
    -O "${QC_READS}/out.${SAMPLE}.R2.fq.gz" \
    --html "${QC_READS}/${SAMPLE}_fastp.html" \
    --json "${QC_READS}/${SAMPLE}_fastp.json"

echo "STAGE_03" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 02/24] fastp complete."