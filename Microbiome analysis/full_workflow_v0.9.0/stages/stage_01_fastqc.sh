#!/bin/bash

################################################################################
# STAGE 01: FastQC - Quality check on raw reads
################################################################################
#
# Input variables (from SLURM template):
#   - SAMPLE        : Sample name
#   - NUM_THREADS   : Number of CPU threads
#   - READ1         : Path to R1 fastq file
#   - READ2         : Path to R2 fastq file
#   - FASTQC_RAWS_DIR : FastQC output directory
#   - STATUS_DIR    : Path to status directory
#
# Output:
#   - FastQC HTML and zip reports in FASTQC_RAWS_DIR
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_02
#
################################################################################

echo "$(date): [STAGE 01/24] Starting FastQC on raw reads for ${SAMPLE}..."

conda activate fastqc-env
fastqc -t ${NUM_THREADS} -o ${FASTQC_RAWS_DIR} "$READ1" "$READ2"
conda deactivate

echo "STAGE_02" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 01/24] FastQC complete."