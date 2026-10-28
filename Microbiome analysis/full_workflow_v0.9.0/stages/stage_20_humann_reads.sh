#!/bin/bash

################################################################################
# STAGE 20: HUMAnN3 reads - Functional profiling on reads
################################################################################
#
# Input variables (from SLURM template):
#   - SAMPLE          : Sample name
#   - NUM_THREADS     : Number of CPU threads
#   - KNEADDATA_OUT   : Directory containing decontaminated paired reads
#   - HUMANN_READS_DIR: Output directory for HUMAnN3 reads mode results
#   - STATUS_DIR      : Path to status directory
#
# Output:
#   - Pathway abundance : ${HUMANN_READS_DIR}/out.${SAMPLE}_combined_pathabundance.tsv
#   - Gene families     : ${HUMANN_READS_DIR}/out.${SAMPLE}_combined_genefamilies.tsv
#   - Normalised tables : *_pathabundance_normalised.tsv / *_genefamilies_normalised.tsv
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_21
#
# Note:
#   - Paired reads are concatenated into a single file for HUMAnN3 input
#   - Combined file is removed on completion to save storage
#   - No --bypass-prescreen / --bypass-nucleotide-index here: those flags
#     only make sense paired with a --taxonomic-profile (or a pre-built
#     custom nucleotide index) telling HUMAnN3 what to search against.
#     An earlier version of this script set them without either, and the
#     accompanying comment implied STAGE_12's MetaPhlAn output was being
#     reused -- it never actually was. Feeding STAGE_12's profile in here
#     was considered and rejected: STAGE_12 runs a newer, independently
#     installed MetaPhlAn version than the one HUMAnN3 bundles/expects,
#     and a --taxonomic-profile from a mismatched MetaPhlAn version risks
#     a clade-naming/database mismatch HUMAnN3 can't reconcile (same
#     reasoning as STAGE_21). So HUMAnN3 now just runs its own normal,
#     version-matched internal prescreen -- slower per sample than the
#     old (mis-configured) bypass, but the standard, correct workflow.
#   - --remove-temp-output added so HUMAnN3's own intermediate directory
#     (bowtie2/diamond output, often much larger than the input reads)
#     doesn't linger the way the combined-reads file used to, before it
#     was explicitly deleted below.
#
################################################################################

set -uo pipefail

echo "$(date): [STAGE 20/24] Starting HUMAnN3 reads mode for ${SAMPLE}..."
echo "  Input R1 : ${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_1.fastq.gz"
echo "  Input R2 : ${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_2.fastq.gz"
echo "  Output   : ${HUMANN_READS_DIR}"

if [ ! -s "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_1.fastq.gz" ] || \
   [ ! -s "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_2.fastq.gz" ]; then
    echo "$(date): [STAGE 20/24] ERROR: decontaminated paired reads not found in ${KNEADDATA_OUT}"
    exit 1
fi

################################################################################
# MetaPhlAn database pinning for HUMAnN3's BUNDLED MetaPhlAn
################################################################################
# Since this stage no longer bypasses the prescreen, HUMAnN3 now calls its own
# bundled MetaPhlAn. By default that checks online for a newer index and
# INTERACTIVELY prompts "[Y/N]" to download it -- which under SLURM (no stdin)
# dies with "EOFError: EOF when reading a line" and kills the HUMAnN3 run.
# Pin the index and run --offline so it never checks or prompts. Keep these
# values identical to STAGE_21's. This is the OLDER database bundled with
# humann3-env, deliberately not the newer one STAGE_12 uses.
METAPHLAN_INDEX="mpa_vJun23_CHOCOPhlAnSGB_202403"
METAPHLAN_DB_DIR=""     # optional; leave empty to use the env's default location

MPA_OPTS="--index ${METAPHLAN_INDEX} --offline"
if [ -n "${METAPHLAN_DB_DIR}" ]; then
    MPA_OPTS="${MPA_OPTS} --bowtie2db ${METAPHLAN_DB_DIR}"
fi
echo "  Bundled-MetaPhlAn options : ${MPA_OPTS}"

READS_COMBINED="${HUMANN_READS_DIR}/out.${SAMPLE}_combined.fastq.gz"

# Concatenate paired reads into a single input file for HUMAnN3. Cat-ing two
# gzip files together produces a valid multi-member gzip stream -- gunzip,
# zcat, and HUMAnN3's own reader all decompress multi-member gzip
# transparently -- so there's no need to decompress first and re-compress.
cat "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_1.fastq.gz" \
    "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_2.fastq.gz" \
    > "$READS_COMBINED"

conda activate humann3-env

humann \
    --input "$READS_COMBINED" \
    --output "${HUMANN_READS_DIR}" \
    --output-basename "out.${SAMPLE}_combined" \
    --metaphlan-options "${MPA_OPTS}" \
    --threads ${NUM_THREADS} \
    --remove-temp-output \
    --verbose \
    --o-log "${HUMANN_READS_DIR}/out.${SAMPLE}_combined_humann3.log"
HUMANN_RC=$?

PATHAB="${HUMANN_READS_DIR}/out.${SAMPLE}_combined_pathabundance.tsv"
GENEFAM="${HUMANN_READS_DIR}/out.${SAMPLE}_combined_genefamilies.tsv"

# Fail loudly rather than letting humann_renorm_table print "Problem opening
# file" and the stage report success on nonexistent tables.
if [ "${HUMANN_RC}" -ne 0 ] || [ ! -s "${PATHAB}" ] || [ ! -s "${GENEFAM}" ]; then
    echo "$(date): [STAGE 20/24] ERROR: HUMAnN3 failed for ${SAMPLE} (exit ${HUMANN_RC})."
    echo "  See ${HUMANN_READS_DIR}/out.${SAMPLE}_combined_humann3.log"
    rm -f "$READS_COMBINED"
    conda deactivate || true
    exit 1
fi

# Normalise pathway abundance and gene families tables
humann_renorm_table \
    --input "${PATHAB}" \
    --output "${HUMANN_READS_DIR}/out.${SAMPLE}_combined_pathabundance_normalised.tsv" \
    --units relab

humann_renorm_table \
    --input "${GENEFAM}" \
    --output "${HUMANN_READS_DIR}/out.${SAMPLE}_combined_genefamilies_normalised.tsv" \
    --units relab

conda deactivate || true

# Remove combined reads file to save storage
echo "$(date): Removing intermediate combined reads file..."
rm -f "$READS_COMBINED"

echo "STAGE_21" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 20/24] HUMAnN3 reads mode complete."