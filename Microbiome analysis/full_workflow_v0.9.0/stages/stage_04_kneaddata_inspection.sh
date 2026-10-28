#!/bin/bash

################################################################################
# STAGE 04: KneadData inspection - FastQC on decontaminated reads
################################################################################
#
# Input variables (from SLURM template):
#   - SAMPLE          : Sample name
#   - NUM_THREADS     : Number of CPU threads
#   - KNEADDATA_OUT   : Directory containing decontaminated reads
#   - FASTQC_CLEAN_DIR: FastQC output directory for cleaned reads
#   - STATUS_DIR      : Path to status directory
#
# Output:
#   - FastQC HTML and zip reports in FASTQC_CLEAN_DIR
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_05
#
################################################################################

echo "$(date): [STAGE 04/24] Starting FastQC on decontaminated reads for ${SAMPLE}..."
echo "  Input R1 : ${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_1.fastq"
echo "  Input R2 : ${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_2.fastq"
echo "  Output   : ${FASTQC_CLEAN_DIR}"

conda activate fastqc-env

fastqc -t ${NUM_THREADS} -o "${FASTQC_CLEAN_DIR}" \
    "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_1.fastq" \
    "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_2.fastq"

conda deactivate || true

# Remove intermediate files to save space
rm -f "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_GRCh38_noalt_as_bowtie2_paired_contam_1.fastq" \
    "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_GRCh38_noalt_as_bowtie2_paired_contam_2.fastq" \
    "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_GRCh38_noalt_as_bowtie2_unmatched_1_contam.fastq" \
    "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_GRCh38_noalt_as_bowtie2_unmatched_2_contam.fastq" \
    "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_unmatched_1.fastq" \
    "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_unmatched_2.fastq"

echo "STAGE_05" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 04/24] KneadData inspection complete."