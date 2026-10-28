#!/bin/bash

################################################################################
# STAGE 12: MetaPhlAn profiling + per-bin MAG cross-check
################################################################################
#
# Part A - MetaPhlAn (whole sample):
#   Runs marker-gene taxonomic profiling on the decontaminated reads
#   (analysis type rel_ab_w_read_stats, which adds "coverage" and
#   "estimated_number_of_reads_from_the_clade" columns alongside the usual
#   clade_name / NCBI_tax_id / relative_abundance columns).
#
# Part B - Per-bin MAG cross-check:
#   For each MAG in final_bins/${SAMPLE}/, pulls the reads belonging to
#   that bin, runs MetaPhlAn on JUST those reads, and takes the top
#   species-level hit. That hit is compared directly against the GTDB-Tk
#   call embedded in the bin's filename, so you can see where the two
#   methods agree/disagree and what MetaPhlAn independently calls each bin.
#
#   HOW READS ARE OBTAINED:
#   Whichever version of a bin (orig / permissive / strict) ended up in
#   final_bins/, its reads are pulled using the contig names from the
#   matching ORIG bin in bin_reassembly/${SAMPLE}/reassembled_bins/
#   (${SAMPLE}_bin<N>.orig.fa). That file's contig headers are unchanged
#   from the original MEGAHIT assembly (e.g. k141_10486), so they match
#   the STAGE 07A assembly BAM directly -- a straightforward samtools
#   region query, no re-alignment needed. This works because permissive/
#   strict reassembly draws its reads from around the orig bin's own
#   contigs; if that ever pulls reads from beyond the orig contig set,
#   this will slightly under-count relative to what SPAdes actually used,
#   but still gives a direct read-based signal for that MAG.
#
# NOTE ON PIPELINE NUMBERING:
#   Part B assumes STAGE 07A has already run (bowtie2 mapping BAM exists)
#   and that binning/reassembly/GTDB-Tk have produced final_bins/ and
#   bin_reassembly/. It does not assume a fixed position among your 24
#   stages -- update the STAGE_13 marker written at the bottom to match
#   whatever your next stage actually is.
#
# Input variables (from SLURM template -- note FINAL_BINS_DIR and
# REASSEMBLY_DIR already resolve to ".../${SAMPLE}", per 02_slurm_template.sh):
#   - SAMPLE             : Sample name
#   - NUM_THREADS        : Number of CPU threads
#   - KNEADDATA_OUT      : Directory containing decontaminated paired reads
#   - METAPHLAN_OUT      : Output directory for MetaPhlAn results
#   - DATABASES          : Path to reference databases
#   - MAPPING_DIR        : ${DIRECTORY}/mapping/${SAMPLE}, containing
#                           ${SAMPLE}.sorted.bam
#   - FINAL_BINS_DIR     : ${DIRECTORY}/final_bins/${SAMPLE}, containing
#                           final MAG fastas + metaquast_references.txt
#   - REASSEMBLY_DIR     : ${DIRECTORY}/bin_reassembly/${SAMPLE}, containing
#                           reassembled_bins/${SAMPLE}_bin<N>.orig.fa
#   - CROSSCHECK_OUT     : ${DIRECTORY}/crosscheck_out (sample subdir is
#                           created by this script)
#   - STATUS_DIR         : Path to status directory
#
# Expected MAG filename pattern:
#   ${SAMPLE}_bin<N>.<orig|permissive|strict>_<Species_with_underscores>.fa
#   e.g. T1_bin1.orig_Enterobacter_flexneri_A.fa
#        T1_bin2.permissive_Campylobacter_D_upsaliensis.fa
#
# Output:
#   - MetaPhlAn profile (whole sample) :
#       ${METAPHLAN_OUT}/${SAMPLE}_metaphlan_profile.txt
#   - Bowtie2 alignments (whole sample) :
#       ${METAPHLAN_OUT}/${SAMPLE}_bowtie2out.bowtie2.bz2
#   - Per-bin MetaPhlAn profile :
#       ${CROSSCHECK_OUT}/${SAMPLE}/<bin>_metaphlan_profile.txt
#   - Cross-check summary table :
#       ${CROSSCHECK_OUT}/${SAMPLE}_MAG_metaphlan_crosscheck.tsv
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_13
#
################################################################################

set -uo pipefail

################################################################################
# PART A: MetaPhlAn profiling (whole sample)
################################################################################

echo "$(date): [STAGE 12/24] Starting MetaPhlAn profiling for ${SAMPLE}..."
echo "  Input R1 : ${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_1.fastq.gz"
echo "  Input R2 : ${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_2.fastq.gz"
echo "  Output   : ${METAPHLAN_OUT}/${SAMPLE}_metaphlan_profile.txt"

conda activate metaphlan-env

MAPOUT_FILE="${METAPHLAN_OUT}/${SAMPLE}_bowtie2out.bowtie2.bz2"
METAPHLAN_PROFILE="${METAPHLAN_OUT}/${SAMPLE}_metaphlan_profile.txt"
DB_DIR="${DATABASES}/metaphlandb/vJan26_202605/"

if [ -s "${MAPOUT_FILE}" ]; then
    echo "$(date): [STAGE 12/24] Found existing BowTie2 mapout file (${MAPOUT_FILE}). Skipping alignment..."

    metaphlan \
        "${MAPOUT_FILE}" \
        --input_type mapout \
        -t rel_ab_w_read_stats \
        -o "${METAPHLAN_PROFILE}" \
        --nproc ${NUM_THREADS} \
        --db_dir "${DB_DIR}"
else
    echo "$(date): [STAGE 12/24] No existing mapout file found. Running full BowTie2 alignment and profiling..."

    metaphlan \
        "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_1.fastq.gz","${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_2.fastq.gz" \
        --input_type fastq \
        -t rel_ab_w_read_stats \
        -o "${METAPHLAN_PROFILE}" \
        --nproc ${NUM_THREADS} \
        --db_dir "${DB_DIR}" \
        --mapout "${MAPOUT_FILE}"
fi

conda deactivate || true

echo "$(date): [STAGE 12/24] MetaPhlAn profiling complete."

################################################################################
# PART B: Per-bin MAG cross-check
################################################################################

# Calculate a staggered delay based on SLURM_ARRAY_TASK_ID
# Each array ID delays by (ID * 15) seconds plus a random 1-10 second jitter
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    DELAY=$(( (SLURM_ARRAY_TASK_ID % 10) * 15 + RANDOM % 10 ))
    echo "$(date): [STAGE 12/24] Staggering task start: sleeping for ${DELAY}s..."
    sleep "${DELAY}"
fi

echo "$(date): [STAGE 12/24] Starting per-bin MAG / MetaPhlAn cross-check for ${SAMPLE}..."

BIN_DIR="${FINAL_BINS_DIR}"
ORIG_BIN_DIR="${REASSEMBLY_DIR}/reassembled_bins"
BAM_FILE="${MAPPING_DIR}/${SAMPLE}.sorted.bam"
BIN_WORK_DIR="${CROSSCHECK_OUT}/${SAMPLE}"
OUT_TABLE="${CROSSCHECK_OUT}/${SAMPLE}_MAG_metaphlan_crosscheck.tsv"

mkdir -p "${BIN_WORK_DIR}"

# -- Sanity checks -----------------------------------------------------------
if [ ! -d "${BIN_DIR}" ]; then
    echo "$(date): [STAGE 12/24] ERROR: bin directory not found: ${BIN_DIR}"
    exit 1
fi
if [ ! -d "${ORIG_BIN_DIR}" ]; then
    echo "$(date): [STAGE 12/24] ERROR: orig-bin directory not found: ${ORIG_BIN_DIR}"
    exit 1
fi
if [ ! -s "${BAM_FILE}" ]; then
    echo "$(date): [STAGE 12/24] ERROR: sorted BAM not found: ${BAM_FILE}"
    exit 1
fi

shopt -s nullglob

# -- Step 1: pull each bin's reads from the assembly BAM using its orig -----
# -- bin's contig names, whichever version (orig/permissive/strict) won -----
conda activate samtools-env

BIN_JOBS=()   # "label|species_raw|r1_fastq|r2_fastq" (or NOREADS markers)

for BIN_FASTA in "${BIN_DIR}"/${SAMPLE}_bin*.fa; do
    BIN_BASENAME=$(basename "${BIN_FASTA}" .fa)

    # ${SAMPLE}_bin<N>.<status>_<Species_with_underscores>
    BIN_LABEL="${BIN_BASENAME%%.*}"                # ${SAMPLE}_bin<N>
    REMAINDER="${BIN_BASENAME#*.}"                 # <status>_<Species_...>
    SPECIES_RAW="${REMAINDER#*_}"                  # strip "orig_"/"permissive_"/"strict_"

    BIN_NUM="${BIN_LABEL#${SAMPLE}_bin}"
    ORIG_FASTA="${ORIG_BIN_DIR}/${SAMPLE}_bin${BIN_NUM}.orig.fa"

    if [ ! -s "${ORIG_FASTA}" ]; then
        echo "$(date): [STAGE 12/24] WARNING: no orig fasta found for ${BIN_LABEL} (${ORIG_FASTA}), skipping."
        continue
    fi

    echo "$(date): [STAGE 12/24] Pulling reads for ${BIN_LABEL} via orig contig names (GTDB-Tk: ${SPECIES_RAW//_/ })..."

    mapfile -t CONTIG_ARR < <(grep '^>' "${ORIG_FASTA}" | sed 's/^>//' | awk '{print $1}')

    if [ "${#CONTIG_ARR[@]}" -eq 0 ]; then
        echo "$(date): [STAGE 12/24] WARNING: no contigs found in ${ORIG_FASTA}, skipping."
        continue
    fi

    BIN_BAM="${BIN_WORK_DIR}/${BIN_LABEL}.bam"
    BIN_NSORT_BAM="${BIN_WORK_DIR}/${BIN_LABEL}.namesorted.bam"
    R1="${BIN_WORK_DIR}/${BIN_LABEL}_R1.fastq"
    R2="${BIN_WORK_DIR}/${BIN_LABEL}_R2.fastq"

    samtools view -b -@ ${NUM_THREADS} "${BAM_FILE}" "${CONTIG_ARR[@]}" > "${BIN_BAM}"
    samtools sort -n -@ ${NUM_THREADS} -o "${BIN_NSORT_BAM}" "${BIN_BAM}"
    samtools fastq -1 "${R1}" -2 "${R2}" -0 /dev/null -s /dev/null -n "${BIN_NSORT_BAM}" > /dev/null 2>&1
    rm -f "${BIN_BAM}" "${BIN_NSORT_BAM}"

    N_READS=$(( $(wc -l < "${R1}") / 4 ))
    echo "$(date): [STAGE 12/24]   ${N_READS} read pairs pulled for ${BIN_LABEL}."

    if [ "${N_READS}" -eq 0 ]; then
        echo "$(date): [STAGE 12/24]   WARNING: no reads found for ${BIN_LABEL}, will skip MetaPhlAn for this bin."
        rm -f "${R1}" "${R2}"
        BIN_JOBS+=("${BIN_LABEL}|${SPECIES_RAW}|NOREADS|NOREADS")
        continue
    fi

    BIN_JOBS+=("${BIN_LABEL}|${SPECIES_RAW}|${R1}|${R2}")
done

conda deactivate || true

# -- Step 2: run MetaPhlAn on each bin's reads and take the top hit ---------
conda activate metaphlan-env

echo -e "bin\tgtdbtk_species\tn_reads\tmetaphlan_top_hit_species\tmetaphlan_top_hit_pct\tgenus_match\tagreement_status" \
    > "${OUT_TABLE}"

for JOB in "${BIN_JOBS[@]}"; do
    IFS='|' read -r BIN_LABEL SPECIES_RAW R1 R2 <<< "${JOB}"
    GTDBTK_SPECIES="${SPECIES_RAW//_/ }"
    GTDBTK_KEY=$(echo "${GTDBTK_SPECIES}" | tr '[:upper:]' '[:lower:]' | xargs)

    if [ "${R1}" == "NOREADS" ]; then
        echo -e "${BIN_LABEL}\t${GTDBTK_SPECIES}\t0\tNA\tNA\tNA\tno_reads" >> "${OUT_TABLE}"
        continue
    fi

    N_READS=$(( $(wc -l < "${R1}") / 4 ))
    BIN_PROFILE="${BIN_WORK_DIR}/${BIN_LABEL}_metaphlan_profile.txt"
    BIN_MAPOUT="${BIN_WORK_DIR}/${BIN_LABEL}_bowtie2out.bowtie2.bz2"

    echo "$(date): [STAGE 12/24] Running MetaPhlAn on ${BIN_LABEL} (${N_READS} read pairs)..."

    # --mapout is required by MetaPhlAn whenever multiple comma-separated
    # FASTQ files are passed, even under --input_type fastq. Not needed
    # downstream here, so it's removed again once the profile is written.
    metaphlan \
        "${R1}","${R2}" \
        --input_type fastq \
        -t rel_ab_w_read_stats \
        -o "${BIN_PROFILE}" \
        --nproc ${NUM_THREADS} \
        --db_dir "${DB_DIR}" \
        --mapout "${BIN_MAPOUT}"

    rm -f "${BIN_MAPOUT}"

    # Top species-level hit = highest relative_abundance among s__ rows
    # (excluding strain-level t__ rows).
    TOP_LINE=$(awk -F'\t' '$1 !~ /^#/ && $1 ~ /s__/ && $1 !~ /t__/' "${BIN_PROFILE}" \
        | sort -t $'\t' -k3,3 -rn | head -n1)

    if [ -z "${TOP_LINE}" ]; then
        TOP_SPECIES=""
        TOP_PCT="NA"
        GENUS_MATCH="NA"
        STATUS="no_metaphlan_hit"
    else
        TOP_CLADE=$(echo "${TOP_LINE}" | awk -F'\t' '{print $1}')
        TOP_PCT=$(echo "${TOP_LINE}" | awk -F'\t' '{print $3}')
        TOP_SPECIES_RAW="${TOP_CLADE##*s__}"
        TOP_SPECIES=$(echo "${TOP_SPECIES_RAW}" | tr '_' ' ')
        TOP_KEY=$(echo "${TOP_SPECIES}" | tr '[:upper:]' '[:lower:]' | xargs)

        if [ "${TOP_KEY}" == "${GTDBTK_KEY}" ]; then
            STATUS="agree"
        else
            STATUS="disagree"
        fi

        # Genus = first whitespace-delimited word of each species name.
        # A simple string comparison (not a taxonomic judgment) so bins
        # that just differ on species/subspecies boundary or GTDB-style
        # placeholder genera (e.g. "Campylobacter D") can be told apart
        # at a glance from genuinely different-genus mismatches.
        GTDBTK_GENUS=$(echo "${GTDBTK_KEY}" | awk '{print $1}')
        TOP_GENUS=$(echo "${TOP_KEY}" | awk '{print $1}')

        if [ "${TOP_GENUS}" == "${GTDBTK_GENUS}" ]; then
            GENUS_MATCH="same_genus"
        else
            GENUS_MATCH="different_genus"
        fi
    fi

    echo -e "${BIN_LABEL}\t${GTDBTK_SPECIES}\t${N_READS}\t${TOP_SPECIES:-NA}\t${TOP_PCT}\t${GENUS_MATCH}\t${STATUS}" \
        >> "${OUT_TABLE}"

    # Drop the reads once profiled to save space; keep the per-bin MetaPhlAn profile.
    rm -f "${R1}" "${R2}"
done

conda deactivate || true

# -- Summary ------------------------------------------------------------
N_BINS=$(($(wc -l < "${OUT_TABLE}") - 1))
N_AGREE=$(awk -F'\t' 'NR>1 && $7=="agree" {c++} END {print c+0}' "${OUT_TABLE}")
N_DISAGREE=$(awk -F'\t' 'NR>1 && $7=="disagree" {c++} END {print c+0}' "${OUT_TABLE}")
N_SAME_GENUS=$(awk -F'\t' 'NR>1 && $6=="same_genus" {c++} END {print c+0}' "${OUT_TABLE}")

echo "$(date): [STAGE 12/24] Cross-check complete for ${SAMPLE}."
echo "  Bins checked : ${N_BINS}"
echo "  Agree        : ${N_AGREE}"
echo "  Disagree     : ${N_DISAGREE}"
echo "  Same genus   : ${N_SAME_GENUS}"
echo "  Output table : ${OUT_TABLE}"

# Adjust this to whatever the next stage in your pipeline actually is
echo "STAGE_13" > "${STATUS_DIR}/${SAMPLE}.resume"

echo "$(date): [STAGE 12/24] MetaPhlAn profiling and MAG cross-check complete."