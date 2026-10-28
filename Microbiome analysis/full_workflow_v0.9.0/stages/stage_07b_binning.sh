#!/bin/bash

################################################################################
# STAGE 07B: Binning - CONCOCT, MetaBAT2, COMEbin, SemiBin2
################################################################################
#
# Input variables (from SLURM template):
#   - SAMPLE        : Sample name
#   - NUM_THREADS   : Number of CPU threads
#   - ASSEMBLY_FILE : Path to MEGAHIT final contigs
#   - MAPPING_DIR   : Directory containing BAM, depth and coverage files
#   - BINNING_DIR   : Output directory for all binner results
#   - STATUS_DIR    : Path to status directory
#
# Output:
#   - CONCOCT bins   : ${BINNING_DIR}/concoct_bins/*.fa
#   - MetaBAT2 bins  : ${BINNING_DIR}/metabat2_bins/*.fa
#   - COMEbin bins   : ${BINNING_DIR}/comebin_bins/*.fa
#   - SemiBin2 bins  : ${BINNING_DIR}/semibin2_bins/*.fa
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_07C
#
# Note:
#   - All four binners use the sorted BAM and coverage files from STAGE_07A
#   - Each binner outputs to its own subdirectory under BINNING_DIR
#   - Contigs are cut into chunks for CONCOCT as per its recommended workflow
#   - CONCOCT chunks are used only for binning — original contigs are used
#     for all downstream steps
#   - SemiBin2 is run in self-supervised mode with the global database
#     as a fallback if no reference is available
#
################################################################################

echo "$(date): [STAGE 07B] Starting binning for ${SAMPLE}..."
echo "  Assembly  : ${ASSEMBLY_FILE}"
echo "  BAM       : ${MAPPING_DIR}/${SAMPLE}.sorted.bam"
echo "  Depth     : ${MAPPING_DIR}/${SAMPLE}_depth.txt"
echo "  Coverage  : ${MAPPING_DIR}/${SAMPLE}_coverage.tsv"
echo "  Output    : ${BINNING_DIR}"

CONCOCT_DIR="${BINNING_DIR}/concoct_bins"
METABAT2_DIR="${BINNING_DIR}/metabat2_bins"
COMEBIN_DIR="${BINNING_DIR}/comebin_bins"
SEMIBIN2_DIR="${BINNING_DIR}/semibin2_bins"
CONCOCT_TMP="${BINNING_DIR}/concoct_tmp"

mkdir -p "${CONCOCT_DIR}" "${METABAT2_DIR}" "${COMEBIN_DIR}" \
         "${SEMIBIN2_DIR}" "${CONCOCT_TMP}"

################################################################################
# BINNER 1: CONCOCT
################################################################################
echo "$(date): [STAGE 07B] Starting CONCOCT binning..."

conda activate concoct-env

# Cut contigs into chunks as recommended by CONCOCT documentation
# 10000bp chunks balance composition signal with contig representation
echo "$(date): [STAGE 07B] Cutting contigs into 10kb chunks for CONCOCT..."
cut_up_fasta.py "${ASSEMBLY_FILE}" \
    -c 10000 \
    -o 0 \
    --merge_last \
    -b "${CONCOCT_TMP}/${SAMPLE}_chunks.bed" \
    > "${CONCOCT_TMP}/${SAMPLE}_chunks.fa"

# Generate coverage table for chunks
echo "$(date): [STAGE 07B] Generating chunk coverage table..."
concoct_coverage_table.py \
    "${CONCOCT_TMP}/${SAMPLE}_chunks.bed" \
    "${MAPPING_DIR}/${SAMPLE}.sorted.bam" \
    > "${CONCOCT_TMP}/${SAMPLE}_chunk_coverage.tsv"

# Run CONCOCT
echo "$(date): [STAGE 07B] Running CONCOCT..."
concoct \
    --composition_file "${CONCOCT_TMP}/${SAMPLE}_chunks.fa" \
    --coverage_file "${CONCOCT_TMP}/${SAMPLE}_chunk_coverage.tsv" \
    --threads ${NUM_THREADS} \
    -b "${CONCOCT_TMP}/"

# Merge clustering results back to original contigs
echo "$(date): [STAGE 07B] Merging CONCOCT chunk clusters back to contigs..."
merge_cutup_clustering.py \
    "${CONCOCT_TMP}/clustering_gt1000.csv" \
    > "${CONCOCT_TMP}/clustering_merged.csv"

# Extract bins as FASTA files
echo "$(date): [STAGE 07B] Extracting CONCOCT bins..."
extract_fasta_bins.py \
    "${ASSEMBLY_FILE}" \
    "${CONCOCT_TMP}/clustering_merged.csv" \
    --output_path "${CONCOCT_DIR}"

conda deactivate || true

# Clean up CONCOCT tmp directory
echo "$(date): [STAGE 07B] Removing CONCOCT intermediate files..."
rm -rf "${CONCOCT_TMP}"

echo "$(date): [STAGE 07B] CONCOCT binning complete."
echo "  Bins: $(ls ${CONCOCT_DIR}/*.fa 2>/dev/null | wc -l) produced"

################################################################################
# BINNER 2: MetaBAT2
################################################################################
echo "$(date): [STAGE 07B] Starting MetaBAT2 binning..."

conda activate metabat2-env

metabat2 \
    -i "${ASSEMBLY_FILE}" \
    -a "${MAPPING_DIR}/${SAMPLE}_depth.txt" \
    -o "${METABAT2_DIR}/bin" \
    -t ${NUM_THREADS} \
    -m 1500 \
    --unbinned \
    -v

conda deactivate || true

echo "$(date): [STAGE 07B] MetaBAT2 binning complete."
echo "  Bins: $(ls ${METABAT2_DIR}/*.fa 2>/dev/null | wc -l) produced"

################################################################################
# BINNER 3: COMEbin
################################################################################
echo "$(date): [STAGE 07B] Starting COMEbin binning..."

conda activate comebin-env

# COMEbin requires a directory of BAM files not a single file path
COMEBIN_BAM_DIR="${COMEBIN_DIR}/bam_input"
mkdir -p "${COMEBIN_BAM_DIR}"
ln -sf "${MAPPING_DIR}/${SAMPLE}.sorted.bam" \
    "${COMEBIN_BAM_DIR}/${SAMPLE}.sorted.bam"
ln -sf "${MAPPING_DIR}/${SAMPLE}.sorted.bam.bai" \
    "${COMEBIN_BAM_DIR}/${SAMPLE}.sorted.bam.bai"

# Remove any incomplete previous COMEbin run to avoid aug directory conflicts
rm -rf "${COMEBIN_DIR}/data_augmentation"

run_comebin.sh \
    -a "${ASSEMBLY_FILE}" \
    -o "${COMEBIN_DIR}" \
    -p "${COMEBIN_BAM_DIR}" \
    -t ${NUM_THREADS}

if [ -d "${COMEBIN_DIR}/comebin_res/comebin_res_bins" ]; then
    mv "${COMEBIN_DIR}/comebin_res/comebin_res_bins/"*.fa \
       "${COMEBIN_DIR}/" 2>/dev/null || true
fi

rm -rf "${COMEBIN_BAM_DIR}"

conda deactivate || true

echo "$(date): [STAGE 07B] COMEbin binning complete."
echo "  Bins: $(ls ${COMEBIN_DIR}/*.fa 2>/dev/null | wc -l) produced"

################################################################################
# BINNER 4: SemiBin2
################################################################################
echo "$(date): [STAGE 07B] Starting SemiBin2 binning..."

conda activate semibin2-env

SemiBin2 single_easy_bin \
    -i "${ASSEMBLY_FILE}" \
    -b "${MAPPING_DIR}/${SAMPLE}.sorted.bam" \
    -o "${SEMIBIN2_DIR}" \
    --threads ${NUM_THREADS} \
    --self-supervised \
    --compression none

# SemiBin2 outputs bins to output_bins subdirectory — move to expected location
if [ -d "${SEMIBIN2_DIR}/output_bins" ]; then
    mv "${SEMIBIN2_DIR}/output_bins/"*.fa "${SEMIBIN2_DIR}/" 2>/dev/null || true
fi

conda deactivate || true

echo "$(date): [STAGE 07B] SemiBin2 binning complete."
echo "  Bins: $(ls ${SEMIBIN2_DIR}/*.fa 2>/dev/null | wc -l) produced"

################################################################################
# SUMMARY
################################################################################
echo "$(date): [STAGE 07B] Binning summary for ${SAMPLE}:"
echo "  CONCOCT  : $(ls ${CONCOCT_DIR}/*.fa  2>/dev/null | wc -l) bins"
echo "  MetaBAT2 : $(ls ${METABAT2_DIR}/*.fa 2>/dev/null | wc -l) bins"
echo "  COMEbin  : $(ls ${COMEBIN_DIR}/*.fa  2>/dev/null | wc -l) bins"
echo "  SemiBin2 : $(ls ${SEMIBIN2_DIR}/*.fa 2>/dev/null | wc -l) bins"

echo "STAGE_07C" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 07B] Binning complete."