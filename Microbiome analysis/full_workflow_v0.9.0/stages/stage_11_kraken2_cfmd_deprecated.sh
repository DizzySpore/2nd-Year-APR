#!/bin/bash

################################################################################
# STAGE 11: Kraken2 cFMD - Food and custom database profiling
################################################################################
#
# Input variables (from SLURM template):
#   - SAMPLE        : Sample name
#   - NUM_THREADS   : Number of CPU threads
#   - KNEADDATA_OUT : Directory containing decontaminated paired reads
#   - KRAKEN2_OUT   : Output directory for Kraken2/Bracken results
#   - DATABASES     : Path to reference databases
#   - STATUS_DIR    : Path to status directory
#
# Output:
#   - cFMD Kraken2 report  : ${KRAKEN2_OUT}/${SAMPLE}.cFMD-kraken2_report.txt
#   - cFMD Bracken profile : ${KRAKEN2_OUT}/${SAMPLE}_cFMD-bracken_abundance.txt
#     (empty file created if no reads classified)
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_12
#
# Note:
#   - Bracken is only run if classified reads are detected in the Kraken2 report
#   - An empty abundance file is created if nothing is classified so downstream
#     checkpoint detection does not fail
#
################################################################################

echo "$(date): [STAGE 11/24] Starting Kraken2 cFMD database profiling for ${SAMPLE}..."

conda activate kraken2-env

kraken2 \
    --db "${DATABASES}/cFMD/kraken2-cFMD-db/cFMD/" \
    --paired --threads ${NUM_THREADS} --confidence 0.1 \
    --report "${KRAKEN2_OUT}/${SAMPLE}.cFMD-kraken2_report.txt" \
    --output "${KRAKEN2_OUT}/${SAMPLE}.cFMD-kraken2_output.txt" \
    "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_1.fastq" \
    "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_2.fastq"

# Only run Bracken if classified reads exist
UNCLASSIFIED_PCT=$(head -n1 "${KRAKEN2_OUT}/${SAMPLE}.cFMD-kraken2_report.txt" | awk '{print $1}')

if (( $(echo "$UNCLASSIFIED_PCT < 100" | bc -l) )); then
    echo "$(date): Classified reads found (${UNCLASSIFIED_PCT}% unclassified). Running Bracken..."
    bracken \
        -d "${DATABASES}/cFMD/kraken2-cFMD-db/cFMD/" \
        -i "${KRAKEN2_OUT}/${SAMPLE}.cFMD-kraken2_report.txt" \
        -o "${KRAKEN2_OUT}/${SAMPLE}_cFMD-bracken_abundance.txt" \
        -r 100
else
    echo "$(date): No classified reads found (100% unclassified). Skipping Bracken."
    touch "${KRAKEN2_OUT}/${SAMPLE}_cFMD-bracken_abundance.txt"
fi

conda deactivate || true

echo "STAGE_12" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 11/24] Kraken2 cFMD profiling complete."