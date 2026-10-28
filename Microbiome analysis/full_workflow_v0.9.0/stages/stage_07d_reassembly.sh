#!/bin/bash

################################################################################
# STAGE 07D: SPAdes reassembly - Strict and permissive bin reassembly
################################################################################
#
# Input variables (from SLURM template):
#   - SAMPLE          : Sample name
#   - NUM_THREADS     : Number of CPU threads
#   - MEMORY          : Memory in GB
#   - MAPPING_DIR     : Directory containing sorted BAM file
#   - REFINEMENT_DIR  : Directory containing refined bins from STAGE_07C
#   - REASSEMBLY_DIR  : Output directory for reassembled bins
#   - COMPLETENESS    : Minimum completeness threshold (%)
#   - CONTAMINATION   : Maximum contamination threshold (%)
#   - DATABASES       : Path to reference databases
#   - STATUS_DIR      : Path to status directory
#
# Output:
#   - Per-bin reassembled FASTA : ${REASSEMBLY_DIR}/reassembled_bins/${BIN}.strict.fa
#                                  ${REASSEMBLY_DIR}/reassembled_bins/${BIN}.permissive.fa
#   - Original bins retained as : ${REASSEMBLY_DIR}/reassembled_bins/${BIN}.orig.fa
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_07E
#
# Note:
#   - Reads are extracted per bin using samtools by filtering the BAM file
#     for reads mapping to contigs belonging to that bin
#   - SPAdes is run without --meta since bin-extracted reads represent a
#     single genome or near single-genome dataset, not a complex mixture.
#     Standard SPAdes single-genome optimisations are more appropriate here.
#   - Permissive uses standard SPAdes settings for maximum completeness
#   - Strict adds --cov-cutoff auto to remove low coverage contigs that
#     are likely contaminants or assembly artefacts, trading completeness
#     for higher purity
#   - All three versions (orig, permissive, strict) are passed to STAGE_07F
#     for quality assessment and final selection
#   - Bins with fewer than 100 read pairs skip reassembly and retain orig only
#
################################################################################

SOFT="/rds/projects/h/hallly-microbiome/software/mamba+conda/miniforge3/envs/metawrap-env/bin/metawrap-scripts"

REFINED_BINS="${REFINEMENT_DIR}/work_files/binsO"
REASSEMBLED_BINS="${REASSEMBLY_DIR}/reassembled_bins"
BAM="${MAPPING_DIR}/${SAMPLE}.sorted.bam"
READS_TMP="${REASSEMBLY_DIR}/reads_tmp"
SPADES_TMP="${REASSEMBLY_DIR}/spades_tmp"

# Remove stale reassembly output
if [ -d "${REASSEMBLY_DIR}" ]; then
    echo "$(date): Removing stale reassembly directory..."
    rm -rf "${REASSEMBLY_DIR}"
fi
mkdir -p "${REASSEMBLY_DIR}"

echo "$(date): [STAGE 07D] Starting bin reassembly for ${SAMPLE}..."
echo "  Refined bins : ${REFINED_BINS}"
echo "  BAM file     : ${BAM}"
echo "  Output       : ${REASSEMBLED_BINS}"

echo "$(date): [STAGE 07D] Cleaning up any previous reassembly output..."
rm -rf "${REASSEMBLED_BINS}"
mkdir -p "${REASSEMBLED_BINS}" "${READS_TMP}" "${SPADES_TMP}"

# Verify BAM file exists
if [ ! -f "${BAM}" ]; then
    echo "$(date): ERROR - BAM file not found: ${BAM}"
    return 1
fi

################################################################################
# REASSEMBLE EACH BIN
################################################################################

BIN_COUNT=0
for bin_file in "${REFINED_BINS}/"*.fa; do
    [ -f "$bin_file" ] || continue
    BIN_NAME=$(basename "$bin_file" .fa)
    BIN_COUNT=$(( BIN_COUNT + 1 ))

    echo "$(date): [STAGE 07D] Processing bin ${BIN_NAME} (${BIN_COUNT})..."

    BIN_READS_DIR="${READS_TMP}/${BIN_NAME}"
    BIN_SPADES_STRICT="${SPADES_TMP}/${BIN_NAME}_strict"
    BIN_SPADES_PERMISSIVE="${SPADES_TMP}/${BIN_NAME}_permissive"
    mkdir -p "${BIN_READS_DIR}" "${BIN_SPADES_STRICT}" "${BIN_SPADES_PERMISSIVE}"

    # Copy original bin immediately — always retained regardless of reassembly
    cp "$bin_file" "${REASSEMBLED_BINS}/${BIN_NAME}.orig.fa"

    # -- Step 1: Extract contig names belonging to this bin --------------------
    echo "$(date): Extracting contig names for ${BIN_NAME}..."
    grep ">" "$bin_file" | sed 's/>//' | awk '{print $1}' \
        > "${BIN_READS_DIR}/contigs.txt"

    CONTIG_COUNT=$(wc -l < "${BIN_READS_DIR}/contigs.txt")
    echo "$(date): ${BIN_NAME} has ${CONTIG_COUNT} contigs"

    # -- Step 2: Extract reads mapping to bin contigs from BAM -----------------
    echo "$(date): Extracting reads mapping to ${BIN_NAME}..."

    conda activate samtools-env

    # Filter BAM to reads mapping to this bin's contigs
    samtools view \
        -@ ${NUM_THREADS} \
        -b "${BAM}" \
        $(cat "${BIN_READS_DIR}/contigs.txt" | tr '\n' ' ') \
        > "${BIN_READS_DIR}/${BIN_NAME}.bam"

    # Sort by read name for paired-end extraction
    samtools sort -n \
        -@ ${NUM_THREADS} \
        -o "${BIN_READS_DIR}/${BIN_NAME}.namesorted.bam" \
        "${BIN_READS_DIR}/${BIN_NAME}.bam"

    rm -f "${BIN_READS_DIR}/${BIN_NAME}.bam"

    # Extract paired reads to FASTQ
    samtools fastq \
        -@ ${NUM_THREADS} \
        -1 "${BIN_READS_DIR}/R1.fastq" \
        -2 "${BIN_READS_DIR}/R2.fastq" \
        -s "${BIN_READS_DIR}/singleton.fastq" \
        -0 /dev/null \
        -n \
        "${BIN_READS_DIR}/${BIN_NAME}.namesorted.bam"

    rm -f "${BIN_READS_DIR}/${BIN_NAME}.namesorted.bam"

    conda deactivate || true

    READ_COUNT=$(wc -l < "${BIN_READS_DIR}/R1.fastq")
    READ_COUNT=$(( READ_COUNT / 4 ))
    echo "$(date): Extracted ${READ_COUNT} read pairs for ${BIN_NAME}"

    # Skip reassembly if too few reads — retain orig only
    if [ "${READ_COUNT}" -lt 100 ]; then
        echo "$(date): WARNING - Too few reads (${READ_COUNT}) for ${BIN_NAME},"
        echo "$(date): Skipping reassembly, retaining original bin only."
        rm -rf "${BIN_READS_DIR}" "${BIN_SPADES_STRICT}" "${BIN_SPADES_PERMISSIVE}"
        continue
    fi

    conda activate spades-env

    # -- Step 3: SPAdes permissive reassembly ----------------------------------
    # Standard SPAdes — maximises completeness
    echo "$(date): Running SPAdes permissive reassembly for ${BIN_NAME}..."

    spades.py \
        -1 "${BIN_READS_DIR}/R1.fastq" \
        -2 "${BIN_READS_DIR}/R2.fastq" \
        -o "${BIN_SPADES_PERMISSIVE}" \
        --threads ${NUM_THREADS} \
        --memory ${MEMORY} || \
    echo "$(date): WARNING - SPAdes permissive failed for ${BIN_NAME}, continuing."

    if [ -f "${BIN_SPADES_PERMISSIVE}/scaffolds.fasta" ]; then
        cp "${BIN_SPADES_PERMISSIVE}/scaffolds.fasta" \
           "${REASSEMBLED_BINS}/${BIN_NAME}.permissive.fa"
        echo "$(date): Permissive reassembly complete for ${BIN_NAME}."
    else
        echo "$(date): WARNING - No permissive scaffolds produced for ${BIN_NAME}."
    fi

    # -- Step 4: SPAdes strict reassembly --------------------------------------
    # --cov-cutoff auto removes low coverage contigs likely to be contaminants
    echo "$(date): Running SPAdes strict reassembly for ${BIN_NAME}..."

    spades.py \
        -1 "${BIN_READS_DIR}/R1.fastq" \
        -2 "${BIN_READS_DIR}/R2.fastq" \
        -o "${BIN_SPADES_STRICT}" \
        --threads ${NUM_THREADS} \
        --memory ${MEMORY} \
        --cov-cutoff auto || \
    echo "$(date): WARNING - SPAdes strict failed for ${BIN_NAME}, continuing."

    if [ -f "${BIN_SPADES_STRICT}/scaffolds.fasta" ]; then
        cp "${BIN_SPADES_STRICT}/scaffolds.fasta" \
           "${REASSEMBLED_BINS}/${BIN_NAME}.strict.fa"
        echo "$(date): Strict reassembly complete for ${BIN_NAME}."
    else
        echo "$(date): WARNING - No strict scaffolds produced for ${BIN_NAME}."
    fi

    conda deactivate || true

    # Clean up per-bin temporary files immediately to save storage
    rm -rf "${BIN_READS_DIR}" \
           "${BIN_SPADES_STRICT}" \
           "${BIN_SPADES_PERMISSIVE}"

    echo "$(date): [STAGE 07D] Bin ${BIN_NAME} complete."
    echo "  orig       : ${REASSEMBLED_BINS}/${BIN_NAME}.orig.fa"
    echo "  permissive : ${REASSEMBLED_BINS}/${BIN_NAME}.permissive.fa"
    echo "  strict     : ${REASSEMBLED_BINS}/${BIN_NAME}.strict.fa"

done

# Clean up all temporary directories
rm -rf "${READS_TMP}" "${SPADES_TMP}"

echo "$(date): [STAGE 07D] Reassembly complete for ${BIN_COUNT} bins."
echo "$(date): All versions written to ${REASSEMBLED_BINS}"

echo "STAGE_07E" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 07D] SPAdes reassembly complete."