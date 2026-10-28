#!/bin/bash

################################################################################
# STAGE 05: MEGAHIT - Metagenomic assembly
################################################################################
#
# Input variables (from SLURM template):
#   - SAMPLE        : Sample name
#   - NUM_THREADS   : Number of CPU threads
#   - KNEADDATA_OUT : Directory containing decontaminated reads
#   - ASSEMBLY_DIR  : Output directory for assembly
#   - ASSEMBLY_FILE : Expected path to final contigs file
#   - STATUS_DIR    : Path to status directory
#
# Output:
#   - Assembled contigs : ${ASSEMBLY_FILE} (final.contigs.fa)
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_06
#
# Note:
#   - Intermediate MEGAHIT files and contamination-removed kneaddata
#     files are deleted on success to save storage
#   - Decontaminated paired reads are retained for binning in STAGE_07
#
################################################################################

echo "$(date): [STAGE 05/24] Starting MEGAHIT assembly for ${SAMPLE}..."
echo "  Input R1 : ${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_1.fastq"
echo "  Input R2 : ${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_2.fastq"
echo "  Output   : ${ASSEMBLY_DIR}"

conda activate megahit-env

megahit \
    -1 "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_1.fastq" \
    -2 "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_2.fastq" \
    --num-cpu-threads ${NUM_THREADS} \
    -o "${ASSEMBLY_DIR}" --force

conda deactivate || true

# Clean up intermediate MEGAHIT files to save storage
echo "$(date): Removing intermediate MEGAHIT files..."
rm -rf "${ASSEMBLY_DIR}/intermediate_contigs" \
       "${ASSEMBLY_DIR}/k*.gfa" \
       "${ASSEMBLY_DIR}/k*.fa" \
       "${ASSEMBLY_DIR}/tmp."* 2>/dev/null || true

# Clean up kneaddata contamination-removed files (not needed downstream)
# Paired reads are retained as they are required for binning in STAGE_07
echo "$(date): Removing intermediate KneadData files..."
rm -f "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata.repeats.removed."* \
      "${KNEADDATA_OUT}/out.${SAMPLE}_kneaddata_GRCh38_noalt_as_bowtie2_paired_contam_"* \
      "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_contam_"* 2>/dev/null || true

echo "STAGE_06" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 05/24] MEGAHIT assembly complete."