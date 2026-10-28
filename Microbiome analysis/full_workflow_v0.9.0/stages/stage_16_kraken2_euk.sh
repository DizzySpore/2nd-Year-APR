#!/bin/bash

################################################################################
# STAGE 16: Kraken2 Euk - Eukaryotic taxonomic profiling
################################################################################
#
# Input variables (from SLURM template):
#   - SAMPLE        : Sample name
#   - NUM_THREADS   : Number of CPU threads
#   - EUKFINDER_OUT : Directory containing EukFinder classified reads
#   - KRAKEN2_OUT   : Output directory for Kraken2/Bracken results
#   - DATABASES     : Path to reference databases
#   - STATUS_DIR    : Path to status directory
#
# Output:
#   - Kraken2 report  : ${KRAKEN2_OUT}/${SAMPLE}_euk-fin_kraken2_report.txt
#   - Bracken profile : ${KRAKEN2_OUT}/${SAMPLE}_euk-fin_bracken_abundance.txt
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_17
#
# Note:
#   - Runs on EukFinder classified reads rather than the full decontaminated
#     read set used in STAGE_10, providing eukaryote-specific profiling
#   - Bracken is skipped and an empty file created if no reads are classified
#
################################################################################

EUK_READS="${EUKFINDER_OUT}/gunzipped_reads/Intermediate_data/Classified_reads"

echo "$(date): [STAGE 16/24] Starting Kraken2 eukaryotic profiling for ${SAMPLE}..."
echo "  Input  : ${EUK_READS}"
echo "  Output : ${KRAKEN2_OUT}"

conda activate kraken2-env

kraken2 \
    --db "${DATABASES}/krakendb/Eukaryote/" \
    --paired --threads ${NUM_THREADS} --confidence 0.1 \
    --report "${KRAKEN2_OUT}/${SAMPLE}_euk-fin_kraken2_report.txt" \
    --output "${KRAKEN2_OUT}/${SAMPLE}_euk-fin_kraken2_output.txt" \
    "${EUK_READS}/${SAMPLE}_R1_classified.fastq" \
    "${EUK_READS}/${SAMPLE}_R2_classified.fastq"

UNCLASSIFIED_PCT=$(head -n1 "${KRAKEN2_OUT}/${SAMPLE}_euk-fin_kraken2_report.txt" | awk '{print $1}')

if (( $(echo "$UNCLASSIFIED_PCT < 100" | bc -l) )); then
    echo "$(date): Classified reads found (${UNCLASSIFIED_PCT}% unclassified). Running Bracken..."
    bracken \
        -d "${DATABASES}/krakendb/Eukaryote/" \
        -i "${KRAKEN2_OUT}/${SAMPLE}_euk-fin_kraken2_report.txt" \
        -o "${KRAKEN2_OUT}/${SAMPLE}_euk-fin_bracken_abundance.txt" \
        -r 100
else
    echo "$(date): No classified reads found. Skipping Bracken."
    touch "${KRAKEN2_OUT}/${SAMPLE}_euk-fin_bracken_abundance.txt"
fi

conda deactivate || true

echo "STAGE_17" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 16/24] Kraken2 eukaryotic profiling complete."