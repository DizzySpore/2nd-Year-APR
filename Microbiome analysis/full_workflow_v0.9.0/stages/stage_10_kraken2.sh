#!/bin/bash

################################################################################
# STAGE 10: Kraken2/Bracken - Read-level taxonomic profiling
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
#   - Prokaryote profile : ${KRAKEN2_OUT}/${SAMPLE}_prok-bracken_abundance.txt
#   - Archaeal profile   : ${KRAKEN2_OUT}/${SAMPLE}_arc-test_bracken_abundance.txt
#   - Eukaryote profile  : ${KRAKEN2_OUT}/${SAMPLE}_euk-test_bracken_abundance.txt
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_11
#
# Note:
#   - Three databases are run: standard (prokaryotes), Archaea, Eukaryote
#   - Bracken read length set to 100bp
#   - Confidence threshold 0.1 applied to all Kraken2 runs
#
################################################################################

echo "$(date): [STAGE 10/24] Starting Kraken2/Bracken taxonomic profiling for ${SAMPLE}..."

conda activate kraken2-env

# Prokaryotes — standard database
echo "$(date): [STAGE 10/24] Running Kraken2/Bracken on prokaryotes..."
kraken2 \
    --db "${DATABASES}/krakendb/PlusPFP" \
    --paired --threads ${NUM_THREADS} --confidence 0.1 \
    --report "${KRAKEN2_OUT}/${SAMPLE}.prok-kraken2_report.txt" \
    --output "${KRAKEN2_OUT}/${SAMPLE}.prok-kraken2_output.txt" \
    "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_1.fastq.gz" \
    "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_2.fastq.gz"

bracken \
    -d "${DATABASES}/krakendb/PlusPFP" \
    -i "${KRAKEN2_OUT}/${SAMPLE}.prok-kraken2_report.txt" \
    -o "${KRAKEN2_OUT}/${SAMPLE}_prok-bracken_abundance.txt" \
    -r 100

# Archaea database
echo "$(date): [STAGE 10/24] Running Kraken2/Bracken on archaea..."
kraken2 \
    --db "${DATABASES}/krakendb/Archaea/" \
    --paired --threads ${NUM_THREADS} --confidence 0.1 \
    --report "${KRAKEN2_OUT}/${SAMPLE}.arc-test.kraken2_report.txt" \
    --output "${KRAKEN2_OUT}/${SAMPLE}.arc-test.kraken2_output.txt" \
    "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_1.fastq.gz" \
    "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_2.fastq.gz"

bracken \
    -d "${DATABASES}/krakendb/Archaea/" \
    -i "${KRAKEN2_OUT}/${SAMPLE}.arc-test.kraken2_report.txt" \
    -o "${KRAKEN2_OUT}/${SAMPLE}_arc-test_bracken_abundance.txt" \
    -r 100

# Eukaryotes database
echo "$(date): [STAGE 10/24] Running Kraken2/Bracken on eukaryotes..."
kraken2 \
    --db "${DATABASES}/krakendb/EuPathDB46/" \
    --paired --threads ${NUM_THREADS} --confidence 0.1 \
    --report "${KRAKEN2_OUT}/${SAMPLE}.euk-test.kraken2_report.txt" \
    --output "${KRAKEN2_OUT}/${SAMPLE}.euk-test.kraken2_output.txt" \
    "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_1.fastq.gz" \
    "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_2.fastq.gz"

bracken \
    -d "${DATABASES}/krakendb/EuPathDB46/" \
    -i "${KRAKEN2_OUT}/${SAMPLE}.euk-test.kraken2_report.txt" \
    -o "${KRAKEN2_OUT}/${SAMPLE}_euk-test_bracken_abundance.txt" \
    -r 100

conda deactivate || true

echo "STAGE_11" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 10/24] Kraken2/Bracken profiling complete."