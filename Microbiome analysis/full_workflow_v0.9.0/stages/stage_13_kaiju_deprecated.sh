#!/bin/bash

################################################################################
# STAGE 13: Kaiju - Protein-level taxonomic profiling
################################################################################
#
# Input variables (from SLURM template):
#   - SAMPLE        : Sample name
#   - NUM_THREADS   : Number of CPU threads
#   - KNEADDATA_OUT : Directory containing decontaminated paired reads
#   - KRAKEN2_OUT   : Directory containing Kraken2 outputs (for merge step)
#   - KAIJU_OUT     : Output directory for Kaiju results
#   - DATABASES     : Path to reference databases
#   - STATUS_DIR    : Path to status directory
#
# Output:
#   - Kaiju output       : ${KAIJU_OUT}/${SAMPLE}.kaiju.out
#   - Kaiju named output : ${KAIJU_OUT}/${SAMPLE}.kaiju.out.names
#   - Kaiju report       : ${KAIJU_OUT}/${SAMPLE}.kaiju.report.txt
#   - Krona plot         : ${KAIJU_OUT}/${SAMPLE}.kaiju.out.html
#   - Combined output    : ${KAIJU_OUT}/${SAMPLE}.kaiju.kraken2.combined.txt
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_14
#
# Note:
#   - Kaiju and Kraken2 prokaryote outputs are merged at the end of this stage
#     using kaiju-mergeOutputs at lowest common ancestor level
#   - Krona visualisation requires a separate krona-env conda environment
#
################################################################################

echo "$(date): [STAGE 13/24] Starting Kaiju profiling for ${SAMPLE}..."
echo "  Input R1 : ${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_1.fastq"
echo "  Input R2 : ${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_2.fastq"
echo "  Output   : ${KAIJU_OUT}"

conda activate kaiju-env

# Run Kaiju
kaiju \
    -t "${DATABASES}/kaijudb/nodes.dmp" \
    -f "${DATABASES}/kaijudb/nr_euk/kaiju_db_nr_euk.fmi" \
    -i "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_1.fastq" \
    -j "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_2.fastq" \
    -o "${KAIJU_OUT}/${SAMPLE}.kaiju.out" \
    -v -z ${NUM_THREADS}

# Convert to Krona format
kaiju2krona \
    -t "${DATABASES}/kaijudb/nodes.dmp" \
    -n "${DATABASES}/kaijudb/names.dmp" \
    -i "${KAIJU_OUT}/${SAMPLE}.kaiju.out" \
    -o "${KAIJU_OUT}/${SAMPLE}.kaiju.out.krona"

conda deactivate || true

# Generate Krona visualisation
conda activate krona-env
ktImportText \
    -o "${KAIJU_OUT}/${SAMPLE}.kaiju.out.html" \
    "${KAIJU_OUT}/${SAMPLE}.kaiju.out.krona"
conda deactivate || true

conda activate kaiju-env

# Add taxonomic names
kaiju-addTaxonNames \
    -t "${DATABASES}/kaijudb/nodes.dmp" \
    -n "${DATABASES}/kaijudb/names.dmp" \
    -i "${KAIJU_OUT}/${SAMPLE}.kaiju.out" \
    -o "${KAIJU_OUT}/${SAMPLE}.kaiju.out.names"

# Generate summary table
kaiju2table \
    -t "${DATABASES}/kaijudb/nodes.dmp" \
    -n "${DATABASES}/kaijudb/names.dmp" \
    -o "${KAIJU_OUT}/${SAMPLE}.kaiju.report.txt" \
    -r species -m 0.5 -v \
    "${KAIJU_OUT}/${SAMPLE}.kaiju.out.names"

# Merge Kaiju and Kraken2 prokaryote outputs at lowest common ancestor
echo "$(date): [STAGE 13/24] Merging Kaiju and Kraken2 outputs..."
sort -k2,2 "${KAIJU_OUT}/${SAMPLE}.kaiju.out" \
    > "${KAIJU_OUT}/${SAMPLE}.kaiju.out.sort"
sort -k2,2 "${KRAKEN2_OUT}/${SAMPLE}.prok-kraken2_output.txt" \
    > "${KRAKEN2_OUT}/${SAMPLE}.prok-kraken2_output.sort"

kaiju-mergeOutputs \
    -i "${KAIJU_OUT}/${SAMPLE}.kaiju.out.sort" \
    -j "${KRAKEN2_OUT}/${SAMPLE}.prok-kraken2_output.sort" \
    -o "${KAIJU_OUT}/${SAMPLE}.kaiju.kraken2.combined.txt" \
    -t "${DATABASES}/kaijudb/nodes.dmp" \
    -v -c lowest

conda deactivate || true

rm -f "${KAIJU_OUT}/${SAMPLE}.kaiju.out.sort" \
      "${KRAKEN2_OUT}/${SAMPLE}.prok-kraken2_output.sort"

echo "STAGE_14" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 13/24] Kaiju profiling and merge complete."