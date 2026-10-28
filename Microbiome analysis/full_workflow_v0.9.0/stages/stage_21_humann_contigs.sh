#!/bin/bash

################################################################################
# STAGE 21: HUMAnN3 per-bin - Functional profiling on MAG-mapped reads
################################################################################
#
# CHANGE OF APPROACH:
#   This stage no longer runs HUMAnN3 on the whole assembly (final.contigs.fa).
#   Instead, for each MAG in final_bins/${SAMPLE}/ it pulls the reads that
#   belong to that bin -- using exactly the same contig-name / BAM-region
#   extraction method as STAGE 12 Part B -- and runs HUMAnN3 on just those
#   reads. This gives per-organism (per-MAG) functional profiles to sit
#   alongside the whole-sample, read-based profile from STAGE_20, instead
#   of one profile built by handing HUMAnN3 raw pooled-community contigs.
#
#   HOW READS ARE OBTAINED (identical logic to STAGE 12 Part B):
#   Whichever version of a bin (orig / permissive / strict) ended up in
#   final_bins/, its reads are pulled using the contig names from the
#   matching ORIG bin in bin_reassembly/${SAMPLE}/reassembled_bins/
#   (${SAMPLE}_bin<N>.orig.fa) against the STAGE 07A assembly BAM. See the
#   STAGE 12 header for the caveat on permissive/strict reassembly pulling
#   reads beyond the orig contig set.
#
#   TAXONOMIC PRESCREEN:
#   Each bin's reads are handed to HUMAnN3 with no --taxonomic-profile and
#   no --bypass-prescreen: HUMAnN3 runs its own bundled MetaPhlAn internally
#   for the prescreen step. STAGE_12's MetaPhlAn profile is deliberately
#   NOT fed in here, since STAGE_12 runs a newer, independently-installed
#   MetaPhlAn version than the one HUMAnN3 bundles/expects -- handing it a
#   taxonomic profile from a mismatched MetaPhlAn version risks a
#   clade-naming/database mismatch HUMAnN3 can't reconcile. The species
#   resolution gain here comes purely from splitting reads per MAG before
#   they ever reach HUMAnN3, not from cross-feeding MetaPhlAn's calls.
#
# NAMING NOTE:
#   Filename / STAGE_21 label / resume checkpoint are all kept as-is so
#   the rest of the pipeline (02_slurm_template.sh, 03_checkpoint.sh)
#   doesn't need touching. ASSEMBLY_FILE is no longer used by this stage.
#   HUMANN_CONTIGS_DIR now holds per-bin output rather than whole-assembly
#   contig output -- consider renaming it (and this file) to
#   HUMANN_BINS_DIR / stage_21_humann_bins.sh if you want the naming to
#   track the new behaviour; that would need matching edits in
#   02_slurm_template.sh and 03_checkpoint.sh.
#
# Input variables (from SLURM template):
#   - SAMPLE             : Sample name
#   - NUM_THREADS         : Number of CPU threads
#   - MAPPING_DIR         : ${DIRECTORY}/mapping/${SAMPLE}, containing
#                            ${SAMPLE}.sorted.bam
#   - FINAL_BINS_DIR      : ${DIRECTORY}/final_bins/${SAMPLE}
#   - REASSEMBLY_DIR      : ${DIRECTORY}/bin_reassembly/${SAMPLE}, containing
#                            reassembled_bins/${SAMPLE}_bin<N>.orig.fa
#   - HUMANN_CONTIGS_DIR  : ${DIRECTORY}/humann_contigs_output/${SAMPLE}
#                            (output dir; name kept from the old contigs
#                            approach, now holds per-bin output)
#   - STATUS_DIR          : Path to status directory
#
# Expected MAG filename pattern (same as STAGE 12):
#   ${SAMPLE}_bin<N>.<orig|permissive|strict>_<Species_with_underscores>.fa
#
# Output (per bin):
#   - ${HUMANN_CONTIGS_DIR}/<bin>_pathabundance.tsv (+ _normalised.tsv)
#   - ${HUMANN_CONTIGS_DIR}/<bin>_genefamilies.tsv  (+ _normalised.tsv)
#   - ${HUMANN_CONTIGS_DIR}/<bin>_pathcoverage.tsv
#   - ${HUMANN_CONTIGS_DIR}/${SAMPLE}_bin_humann_manifest.tsv : which bins
#       were run and their read counts -- useful input for a future
#       STAGE_22 merge
#   - Completion marker : ${HUMANN_CONTIGS_DIR}/humann_loop_completed.txt
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_22
#
################################################################################

set -uo pipefail

################################################################################
# MetaPhlAn database pinning for HUMAnN3's BUNDLED MetaPhlAn
################################################################################
# HUMAnN3 calls its own bundled MetaPhlAn for the prescreen. By default that
# MetaPhlAn checks online for a newer index and, if one exists, INTERACTIVELY
# prompts "[Y/N]" to download it. Under SLURM there is no stdin, so it dies
# with "EOFError: EOF when reading a line", taking the whole HUMAnN3 run for
# that bin with it.
#
# Fix: pin the index explicitly and run --offline so it never checks for or
# prompts about an upgrade. Note this is the OLDER database bundled with
# humann3-env -- deliberately NOT the newer vJan26 database STAGE_12 uses.
# The two are kept separate on purpose (see header note): HUMAnN3's
# ChocoPhlAn nucleotide database is version-matched to its own bundled
# MetaPhlAn index, and mixing them is what we're avoiding.
#
# To confirm the right values on your system:
#   conda activate humann3-env
#   python -c "import metaphlan, os; print(os.path.join(os.path.dirname(metaphlan.__file__),'metaphlan_databases'))"
#   ls "$(that path)"      # look for mpa_v*_CHOCOPhlAnSGB_*.pkl
# then set the two variables below to match. If HUMAnN3's MetaPhlAn database
# lives in the default in-package location you can leave METAPHLAN_DB_DIR
# empty and only set the index.
METAPHLAN_INDEX="mpa_vJun23_CHOCOPhlAnSGB_202403"
METAPHLAN_DB_DIR=""     # optional; leave empty to use the env's default location

MPA_OPTS="--index ${METAPHLAN_INDEX} --offline"
if [ -n "${METAPHLAN_DB_DIR}" ]; then
    MPA_OPTS="${MPA_OPTS} --bowtie2db ${METAPHLAN_DB_DIR}"
fi

echo "$(date): [STAGE 21/24] Starting per-bin HUMAnN3 for ${SAMPLE}..."
echo "  Bundled-MetaPhlAn options : ${MPA_OPTS}"

BIN_DIR="${FINAL_BINS_DIR}"
ORIG_BIN_DIR="${REASSEMBLY_DIR}/reassembled_bins"
BAM_FILE="${MAPPING_DIR}/${SAMPLE}.sorted.bam"
TMP_DIR="${HUMANN_CONTIGS_DIR}/tmp_reads"
MANIFEST="${HUMANN_CONTIGS_DIR}/${SAMPLE}_bin_humann_manifest.tsv"

mkdir -p "${HUMANN_CONTIGS_DIR}" "${TMP_DIR}"

# -- Sanity checks -----------------------------------------------------------
if [ ! -d "${BIN_DIR}" ]; then
    echo "$(date): [STAGE 21/24] ERROR: bin directory not found: ${BIN_DIR}"
    exit 1
fi
if [ ! -d "${ORIG_BIN_DIR}" ]; then
    echo "$(date): [STAGE 21/24] ERROR: orig-bin directory not found: ${ORIG_BIN_DIR}"
    exit 1
fi
if [ ! -s "${BAM_FILE}" ]; then
    echo "$(date): [STAGE 21/24] ERROR: sorted BAM not found: ${BAM_FILE}"
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
    REMAINDER="${BIN_BASENAME#*.}"                  # <status>_<Species_...>
    SPECIES_RAW="${REMAINDER#*_}"                   # strip "orig_"/"permissive_"/"strict_"

    BIN_NUM="${BIN_LABEL#${SAMPLE}_bin}"
    ORIG_FASTA="${ORIG_BIN_DIR}/${SAMPLE}_bin${BIN_NUM}.orig.fa"

    if [ ! -s "${ORIG_FASTA}" ]; then
        echo "$(date): [STAGE 21/24] WARNING: no orig fasta found for ${BIN_LABEL} (${ORIG_FASTA}), skipping."
        continue
    fi

    echo "$(date): [STAGE 21/24] Pulling reads for ${BIN_LABEL} via orig contig names (GTDB-Tk: ${SPECIES_RAW//_/ })..."

    mapfile -t CONTIG_ARR < <(grep '^>' "${ORIG_FASTA}" | sed 's/^>//' | awk '{print $1}')

    if [ "${#CONTIG_ARR[@]}" -eq 0 ]; then
        echo "$(date): [STAGE 21/24] WARNING: no contigs found in ${ORIG_FASTA}, skipping."
        continue
    fi

    BIN_BAM="${TMP_DIR}/${BIN_LABEL}.bam"
    BIN_NSORT_BAM="${TMP_DIR}/${BIN_LABEL}.namesorted.bam"
    R1="${TMP_DIR}/${BIN_LABEL}_R1.fastq"
    R2="${TMP_DIR}/${BIN_LABEL}_R2.fastq"

    samtools view -b -@ ${NUM_THREADS} "${BAM_FILE}" "${CONTIG_ARR[@]}" > "${BIN_BAM}"
    samtools sort -n -@ ${NUM_THREADS} -o "${BIN_NSORT_BAM}" "${BIN_BAM}"
    samtools fastq -1 "${R1}" -2 "${R2}" -0 /dev/null -s /dev/null -n "${BIN_NSORT_BAM}" > /dev/null 2>&1
    rm -f "${BIN_BAM}" "${BIN_NSORT_BAM}"

    N_READS=$(( $(wc -l < "${R1}") / 4 ))
    echo "$(date): [STAGE 21/24]   ${N_READS} read pairs pulled for ${BIN_LABEL}."

    if [ "${N_READS}" -eq 0 ]; then
        echo "$(date): [STAGE 21/24]   WARNING: no reads found for ${BIN_LABEL}, will skip HUMAnN3 for this bin."
        rm -f "${R1}" "${R2}"
        BIN_JOBS+=("${BIN_LABEL}|${SPECIES_RAW}|NOREADS|NOREADS")
        continue
    fi

    BIN_JOBS+=("${BIN_LABEL}|${SPECIES_RAW}|${R1}|${R2}")
done

conda deactivate || true

if [ ${#BIN_JOBS[@]} -eq 0 ]; then
    echo "$(date): [STAGE 21/24] WARNING: no bins matching ${SAMPLE}_bin*.fa found in ${BIN_DIR}."
fi

# -- Step 2: run HUMAnN3 on each bin's reads ---------------------------------
conda activate humann3-env

echo -e "bin\tgtdbtk_species\tn_reads\tstatus" > "${MANIFEST}"

for JOB in "${BIN_JOBS[@]}"; do
    IFS='|' read -r BIN_LABEL SPECIES_RAW R1 R2 <<< "${JOB}"
    GTDBTK_SPECIES="${SPECIES_RAW//_/ }"

    if [ "${R1}" == "NOREADS" ]; then
        echo -e "${BIN_LABEL}\t${GTDBTK_SPECIES}\t0\tno_reads" >> "${MANIFEST}"
        continue
    fi

    N_READS=$(( $(wc -l < "${R1}") / 4 ))
    READS_COMBINED="${TMP_DIR}/${BIN_LABEL}_combined.fastq"
    cat "${R1}" "${R2}" > "${READS_COMBINED}"

    echo "$(date): [STAGE 21/24] Running HUMAnN3 on ${BIN_LABEL} (${N_READS} read pairs)..."

    humann \
        --input "${READS_COMBINED}" \
        --output "${HUMANN_CONTIGS_DIR}" \
        --output-basename "${BIN_LABEL}" \
        --metaphlan-options "${MPA_OPTS}" \
        --threads ${NUM_THREADS} \
        --remove-temp-output \
        --verbose \
        --o-log "${HUMANN_CONTIGS_DIR}/${BIN_LABEL}_humann3.log"
    HUMANN_RC=$?

    PATHAB="${HUMANN_CONTIGS_DIR}/${BIN_LABEL}_pathabundance.tsv"
    GENEFAM="${HUMANN_CONTIGS_DIR}/${BIN_LABEL}_genefamilies.tsv"

    # Don't blindly renorm: if HUMAnN3 died (e.g. the bundled MetaPhlAn
    # failing), the tables won't exist and humann_renorm_table just prints
    # "Problem opening file" while the loop marches on marking bins "ok".
    # Record the failure in the manifest instead so it's visible later.
    if [ "${HUMANN_RC}" -ne 0 ] || [ ! -s "${PATHAB}" ] || [ ! -s "${GENEFAM}" ]; then
        echo "$(date): [STAGE 21/24]   ERROR: HUMAnN3 failed for ${BIN_LABEL} (exit ${HUMANN_RC}); see ${HUMANN_CONTIGS_DIR}/${BIN_LABEL}_humann3.log"
        echo -e "${BIN_LABEL}\t${GTDBTK_SPECIES}\t${N_READS}\thumann_failed" >> "${MANIFEST}"
        rm -f "${R1}" "${R2}" "${READS_COMBINED}"
        continue
    fi

    # Normalise pathway abundance and gene families tables
    humann_renorm_table \
        --input "${PATHAB}" \
        --output "${HUMANN_CONTIGS_DIR}/${BIN_LABEL}_pathabundance_normalised.tsv" \
        --units relab

    humann_renorm_table \
        --input "${GENEFAM}" \
        --output "${HUMANN_CONTIGS_DIR}/${BIN_LABEL}_genefamilies_normalised.tsv" \
        --units relab

    echo -e "${BIN_LABEL}\t${GTDBTK_SPECIES}\t${N_READS}\tok" >> "${MANIFEST}"

    # Drop the per-bin reads once profiled to save space; keep the per-bin
    # HUMAnN3 output tables.
    rm -f "${R1}" "${R2}" "${READS_COMBINED}"
done

conda deactivate || true

rmdir "${TMP_DIR}" 2>/dev/null || true

# -- Summary -------------------------------------------------------------
N_BINS=$(($(wc -l < "${MANIFEST}") - 1))
N_OK=$(awk -F'\t' 'NR>1 && $4=="ok" {c++} END {print c+0}' "${MANIFEST}")
N_NOREADS=$(awk -F'\t' 'NR>1 && $4=="no_reads" {c++} END {print c+0}' "${MANIFEST}")
N_FAILED=$(awk -F'\t' 'NR>1 && $4=="humann_failed" {c++} END {print c+0}' "${MANIFEST}")

echo "$(date): [STAGE 21/24] Per-bin HUMAnN3 complete for ${SAMPLE}."
echo "  Bins listed         : ${N_BINS}"
echo "  Bins profiled       : ${N_OK}"
echo "  Bins with no reads  : ${N_NOREADS}"
echo "  Bins failed         : ${N_FAILED}"
echo "  Manifest            : ${MANIFEST}"

# Don't mark the stage complete if nothing actually profiled -- otherwise the
# checkpoint scanner sees humann_loop_completed.txt and skips STAGE_21 on the
# next run, hiding a total failure like the MetaPhlAn database prompt.
if [ "${N_OK}" -eq 0 ]; then
    echo "$(date): [STAGE 21/24] ERROR: no bins profiled successfully for ${SAMPLE}; not writing completion marker."
    exit 1
fi

if [ "${N_FAILED}" -gt 0 ]; then
    echo "$(date): [STAGE 21/24] WARNING: ${N_FAILED} bin(s) failed; see manifest and per-bin logs."
fi

# Write completion marker used by checkpoint scanner
touch "${HUMANN_CONTIGS_DIR}/humann_loop_completed.txt"

echo "STAGE_22" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 21/24] HUMAnN3 per-bin mode complete."