#!/bin/bash

################################################################################
# STAGE 03: KneadData - Human read decontamination
################################################################################
#
# Input variables (from SLURM template):
#   - SAMPLE        : Sample name
#   - NUM_THREADS   : Number of CPU threads
#   - QC_READS      : Directory containing fastp trimmed reads
#   - KNEADDATA_OUT : Output directory for decontaminated reads
#   - DATABASES     : Path to reference databases
#   - STATUS_DIR    : Path to status directory
#
# Output:
#   - Decontaminated reads : ${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_1.fastq
#                            ${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_2.fastq
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_04
#
# Note:
#   - --bypass-trim used as trimming was already performed by fastp in STAGE_02
#   - Trimmed reads from STAGE_02 are removed on success to save storage
#
################################################################################

echo "$(date): [STAGE 03/24] Starting KneadData for ${SAMPLE}..."
echo "  Input R1  : ${QC_READS}/out.${SAMPLE}.R1.fq.gz"
echo "  Input R2  : ${QC_READS}/out.${SAMPLE}.R2.fq.gz"
echo "  Reference : ${DATABASES}/human-GRCh38_no-alt/GRCh38_noalt_as/"
echo "  Output    : ${KNEADDATA_OUT}"

conda activate kneaddata-env

kneaddata --threads ${NUM_THREADS} --bypass-trim --bypass-trf --reorder \
    --input1 "${QC_READS}/out.${SAMPLE}.R1.fq.gz" \
    --input2 "${QC_READS}/out.${SAMPLE}.R2.fq.gz" \
    -db "${DATABASES}/human-GRCh38_no-alt/GRCh38_noalt_as/" \
    --output "${KNEADDATA_OUT}"

conda deactivate || true

# Clean up fastp trimmed reads now that kneaddata has consumed them
echo "$(date): Removing intermediate fastp reads to save storage..."
rm -f "${QC_READS}/out.${SAMPLE}.R1.fq.gz" \
      "${QC_READS}/out.${SAMPLE}.R2.fq.gz"

echo "STAGE_04" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 03/24] KneadData complete."