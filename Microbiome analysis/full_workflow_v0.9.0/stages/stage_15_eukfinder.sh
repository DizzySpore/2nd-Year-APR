#!/bin/bash
################################################################################
# STAGE 15: EukFinder - Eukaryotic sequence detection
################################################################################
#
# Input variables (from SLURM template):
#   - SAMPLE        : Sample name
#   - READ1, READ2  : Raw paired-end input fastq(.gz) files for the sample
#   - NUM_THREADS   : Number of CPU threads
#   - EUKFINDER_OUT : Output directory for EukFinder results
#   - DATABASES     : Path to reference databases
#   - STATUS_DIR    : Path to status directory
#   - DIRECTORY     : Working directory (returned to at end of stage)
#
# Output:
#   - EukFinder results: ${EUKFINDER_OUT}/gunzipped_reads/Eukfinder_results/
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_16
#
# Workflow:
#   EukFinder is run on raw reads (not assembled bins) using its two-step
#   short-read pipeline:
#     1. `eukfinder read_prep`  - adapter/quality trimming + human host
#                                  read removal + Centrifuge classification
#     2. `eukfinder short_seqs` - PLAST-based classification of the
#                                  read_prep output into
#                                  Euk/Unk/Bact/Arch/Misc/EUnk categories
#
#   This runs independently of KneadData/fastp QC done earlier in the
#   pipeline (STAGE_02/STAGE_03) - EukFinder does its own trimming and
#   host decontamination internally as part of `read_prep`.
#
#   Failures are logged to ${EUKFINDER_OUT}/error.log rather than treated
#   as fatal, since read_prep/short_seqs issues on a given sample
#   shouldn't necessarily block the rest of the workflow. Check that log
#   after the run to confirm this stage actually completed successfully
#   for each sample rather than silently no-oping.
#
################################################################################
echo "$(date): [STAGE 15/24] Starting EukFinder eukaryotic detection for ${SAMPLE}..."
echo "  Input reads : ${READ1}, ${READ2}"
echo "  Output      : ${EUKFINDER_OUT}"

mkdir -p "${EUKFINDER_OUT}"

set +e
conda activate eukfinder-env

# -- Stage reads into a scratch subfolder and gunzip -------------------------
mkdir -p "${EUKFINDER_OUT}/gunzipped_reads/"
cd "${EUKFINDER_OUT}/gunzipped_reads/" || {
    echo "$(date): ERROR - could not cd into ${EUKFINDER_OUT}/gunzipped_reads/" >> "${EUKFINDER_OUT}/error.log"
    conda deactivate
    set -e
    cd "${DIRECTORY}"
    return 1
}

cp "$READ1" .
cp "$READ2" .

echo "$(date): Gunzipping reads..."
gunzip "${SAMPLE}"*

GUNZIP1=""
GUNZIP2=""
for file in ${SAMPLE}*.fastq; do
    if [[ "$file" =~ _R1(_[0-9]{3})?\.fastq$ || "$file" =~ _1(_sequence)?\.fastq$ ]]; then
        GUNZIP1="$file"
    elif [[ "$file" =~ _R2(_[0-9]{3})?\.fastq$ || "$file" =~ _2(_sequence)?\.fastq$ ]]; then
        GUNZIP2="$file"
    fi
done

if [ -z "$GUNZIP1" ] || [ -z "$GUNZIP2" ]; then
    echo "$(date): ERROR - could not identify gunzipped R1/R2 for ${SAMPLE}" >> "${EUKFINDER_OUT}/error.log"
    conda deactivate
    set -e
    cd "${DIRECTORY}"
    return 1
fi

# -- Step 1: read_prep (trim + host decontamination + Centrifuge) ------------
eukfinder read_prep \
    --r1 "$GUNZIP1" \
    --r2 "$GUNZIP2" \
    -o "${SAMPLE}_eukfinder_read_prep" \
    -n "$NUM_THREADS" \
    --hcrop 10 -l 15 -t 15 --wsize 40 --qscore 25 --mlen 20 --mhlen 20 \
    -i "${DATABASES}/eukfinder-db/illumina_adapters/TrueSeq2_NexteraSE-PE.fa" \
    --hg "${DATABASES}/eukfinder-db/GRCh38.p13/GCF_000001405.39_GRCh38.p13_human_genome.fna" \
    --cdb "${DATABASES}/eukfinder-db/CentrifugeDB/Centrifuge_DB_2024"
READ_PREP_STATUS=$?
if [ $READ_PREP_STATUS -ne 0 ]; then
    echo "$(date): eukfinder read_prep failed for ${SAMPLE}" >> "${EUKFINDER_OUT}/error.log"
fi

rm -f "$GUNZIP1" "$GUNZIP2"

# -- Step 2: short_seqs (PLAST classification of read_prep output) ----------
if [ $READ_PREP_STATUS -eq 0 ]; then
    eukfinder short_seqs \
        --r1 "${SAMPLE}_eukfinder_read_prep_p.1.fastq" \
        --r2 "${SAMPLE}_eukfinder_read_prep_p.2.fastq" \
        --un "${SAMPLE}_eukfinder_read_prep_un.fastq" \
        -o "${SAMPLE}_eukfinder_short_seq_" \
        -t T --max_m 100 -e 0.01 --pid 60 --cov 30 --mhlen 50 -n "$NUM_THREADS" -z 10 \
        --pclass "${SAMPLE}_eukfinder_read_prep_centrifuge_P" \
        --uclass "${SAMPLE}_eukfinder_read_prep_centrifuge_UP" \
        --plast-database "${DATABASES}/eukfinder-db/PlastDB/PlastDB.fasta" \
        --plast-id-map "${DATABASES}/eukfinder-db/PlastDB/PlastDB_map.txt" \
        --cdb "${DATABASES}/eukfinder-db/CentrifugeDB/Centrifuge_DB_2024"
    SHORT_SEQS_STATUS=$?
    if [ $SHORT_SEQS_STATUS -ne 0 ]; then
        echo "$(date): eukfinder short_seqs failed for ${SAMPLE}" >> "${EUKFINDER_OUT}/error.log"
    fi
else
    echo "$(date): Skipping short_seqs for ${SAMPLE} - read_prep failed." >> "${EUKFINDER_OUT}/error.log"
fi

conda deactivate
set -e

# -- Cleanup: keep only the final Euk/Unk/EUnk contigs, discard the rest ----
# Eukfinder_short's intermediates (metaSPAdes assembly, per-category raw
# reads, Centrifuge intermediates, whole-sample trimmed reads) scale badly
# across hundreds of samples and are fully superseded once short_seqs has
# produced the final classified contig FASTAs. Set EUKFINDER_CLEANUP=false
# above the pipeline call (or comment this block out) to keep everything
# for debugging instead.
EUKFINDER_CLEANUP="${EUKFINDER_CLEANUP:-true}"

if [ "$EUKFINDER_CLEANUP" = "true" ] && [ "${SHORT_SEQS_STATUS:-1}" -eq 0 ]; then
    BEFORE_SIZE=$(du -sh "${EUKFINDER_OUT}/gunzipped_reads" 2>/dev/null | cut -f1)

    RESULTS_DIR="${EUKFINDER_OUT}/gunzipped_reads/Eukfinder_results"
    INTERMEDIATE_DIR="${EUKFINDER_OUT}/gunzipped_reads/Intermediate_data"

    # Non-eukaryotic, non-unknown final contig sets
    rm -f "${RESULTS_DIR}/scf_${SAMPLE}_eukfinder_short_seq_.Arch.fasta" \
          "${RESULTS_DIR}/scf_${SAMPLE}_eukfinder_short_seq_.Bact.fasta" \
          "${RESULTS_DIR}/scf_${SAMPLE}_eukfinder_short_seq_.Misc.fasta"

    # All intermediates: per-category raw reads, metaSPAdes working dir,
    # Centrifuge contig classification, scratch/tmp folders
    rm -rf "${INTERMEDIATE_DIR}"

    # Whole-sample read_prep outputs and their Centrifuge classification,
    # now fully consumed as input to short_seqs
    rm -f "${EUKFINDER_OUT}/gunzipped_reads/${SAMPLE}_eukfinder_read_prep_p.1.fastq" \
          "${EUKFINDER_OUT}/gunzipped_reads/${SAMPLE}_eukfinder_read_prep_p.2.fastq" \
          "${EUKFINDER_OUT}/gunzipped_reads/${SAMPLE}_eukfinder_read_prep_un.fastq" \
          "${EUKFINDER_OUT}/gunzipped_reads/${SAMPLE}_eukfinder_read_prep_centrifuge_P" \
          "${EUKFINDER_OUT}/gunzipped_reads/${SAMPLE}_eukfinder_read_prep_centrifuge_UP"

    AFTER_SIZE=$(du -sh "${EUKFINDER_OUT}/gunzipped_reads" 2>/dev/null | cut -f1)
    echo "$(date): Cleanup complete for ${SAMPLE} (${BEFORE_SIZE} -> ${AFTER_SIZE}). Kept Euk/Unk/EUnk fasta + summary_table.txt."
elif [ "$EUKFINDER_CLEANUP" = "true" ]; then
    echo "$(date): Skipping cleanup for ${SAMPLE} - short_seqs did not complete successfully, keeping intermediates for debugging."
fi

cd "${DIRECTORY}"

echo "$(date): EukFinder analysis complete for ${SAMPLE}."
echo "STAGE_16" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 15/24] EukFinder complete."