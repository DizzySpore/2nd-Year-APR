#!/bin/bash
################################################################################
# STAGE 11: MEDI - Dietary/food metagenomic quantification (Nextflow, batch)
################################################################################
#
# Unlike other stages, MEDI's quant.nf processes ALL samples in one Nextflow
# run. Since this script runs once per array task, only the task that arrives
# LAST (i.e. every other sample in this submission has already checked in)
# claims the run. Only samples in samples.txt (this submission) are quantified.
#
# Output:
#   - ${STATUS_DIR}/medi_complete.txt  - shared batch completion flag
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_12
#
################################################################################

echo "$(date): [STAGE 11/24] MEDI check-in for ${SAMPLE}..."

if [ -f "${STATUS_DIR}/medi_complete.txt" ]; then
    echo "$(date): MEDI already complete for this batch — skipping."
else
    touch "${STATUS_DIR}/stage11_arrived.${SLURM_ARRAY_JOB_ID}.${SAMPLE}"
    TOTAL_SAMPLES=$(wc -l < "${DIRECTORY}/samples.txt")
    ARRIVED=$(find "${STATUS_DIR}" -maxdepth 1 -name "stage11_arrived.${SLURM_ARRAY_JOB_ID}.*" | wc -l)
    echo "$(date): ${ARRIVED}/${TOTAL_SAMPLES} samples have reached STAGE_11."

    if [ "$ARRIVED" -ge "$TOTAL_SAMPLES" ] && mkdir "${STATUS_DIR}/stage11_medi.lock.${SLURM_ARRAY_JOB_ID}" 2>/dev/null; then
        echo "$(date): ${SAMPLE} is the last arrival — launching MEDI batch run..."

        MEDI_DIR="/rds/projects/h/hallly-microbiome/Leon/software/medi"   # personal dir, staying as-is for now
        MEDI_DB_PATH="${MEDI_DIR}/manual_medi_db"
        RUN_DIR="${DIRECTORY}/medi_run_output"

        mkdir -p "${RUN_DIR}/data/raw" "${MEDI_DIR}/data/raw"

        # Only symlink samples from THIS submission's samples.txt, not every
        # KneadData output that happens to exist in DEHUMAN_READS
        while IFS= read -r MEDI_SAMPLE; do
            [ -z "$MEDI_SAMPLE" ] && continue
            for f in "${DEHUMAN_READS}/kneaddata_output_${MEDI_SAMPLE}"/out.${MEDI_SAMPLE}.R1_kneaddata_paired_[12].fastq; do
                [ -e "$f" ] || continue
                ln -sf "$f" "${RUN_DIR}/data/raw/"
                ln -sf "$f" "${MEDI_DIR}/data/raw/"
            done
        done < "${DIRECTORY}/samples.txt"

        ln -sf "${MEDI_DIR}/bin"                       "${RUN_DIR}/bin"
        ln -sf "${MEDI_DIR}/quant.nf"                  "${RUN_DIR}/quant.nf"
        ln -sf "${MEDI_DIR}/data/dbs"                  "${RUN_DIR}/data/dbs"
        ln -sf "${MEDI_DIR}/data/sequences"            "${RUN_DIR}/data/sequences"
        ln -sf "${MEDI_DIR}/data/sketches"             "${RUN_DIR}/data/sketches"
        ln -sf "${MEDI_DIR}/data/uncompressed_sequences" "${RUN_DIR}/data/uncompressed_sequences"
        ln -sf "${MEDI_DIR}/data/manifest.csv"         "${RUN_DIR}/data/manifest.csv"
        ln -sf "${MEDI_DIR}/data/matches.csv"          "${RUN_DIR}/data/matches.csv"
        ln -sf "${MEDI_DIR}/data/mash_ani.csv"         "${RUN_DIR}/data/mash_ani.csv"

        pushd "${MEDI_DIR}" > /dev/null
        conda activate medi-lxm697

        nextflow run "${RUN_DIR}/quant.nf" \
            -resume \
            --db="${MEDI_DB_PATH}" \
            --executor.cpus="${NUM_THREADS}" \
            --executor.memory="${MEMORY} GB"
        NEXTFLOW_EXIT=$?

        conda deactivate || true
        popd > /dev/null

        if [ "$NEXTFLOW_EXIT" -ne 0 ]; then
            echo "Error: MEDI/Nextflow run failed (exit ${NEXTFLOW_EXIT})" >&2
            rmdir "${STATUS_DIR}/stage11_medi.lock.${SLURM_ARRAY_JOB_ID}" 2>/dev/null || true
            exit 1
        fi

        touch "${STATUS_DIR}/medi_complete.txt"
        echo "$(date): [STAGE 11/24] MEDI batch run complete."
    else
        echo "$(date): Not the last arrival (or lock already claimed) — deferring MEDI."
    fi
fi

echo "STAGE_12" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 11/24] Done for ${SAMPLE}."