#!/bin/bash

################################################################################
# STAGE 08: GTDB-Tk - Bin taxonomic classification and renaming
################################################################################
#
# Input variables (from SLURM template):
#   - SAMPLE            : Sample name
#   - NUM_THREADS       : Number of CPU threads
#   - REASSEMBLED_BINS_DIR: Directory containing reassembled bins
#   - GTDBTK_OUT        : Output directory for GTDB-Tk results
#   - DATABASES         : Path to reference databases
#   - STATUS_DIR        : Path to status directory
#
# Output:
#   - GTDB-Tk summary  : ${GTDBTK_OUT}/classify/gtdbtk.bac120.summary.tsv
#                        ${GTDBTK_OUT}/classify/gtdbtk.ar53.summary.tsv
#   - Renamed bins     : ${FINAL_BINS_DIR}/${SAMPLE}_${BIN}_${SPECIES}.fa
#   - MetaQUAST refs   : ${METAQUAST_REFS_FILE}
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_09
#
# Note:
#   - Bins are renamed in place with species information appended to filename
#   - MetaQUAST references file is generated for use in STAGE_09
#   - Bins with no species assignment are suffixed with genus_sp
#   - Bins with no genus or species (N/A) are labelled Unknown_organism
#
################################################################################

echo "$(date): [STAGE 08/24] Starting GTDB-Tk classification for ${SAMPLE}..."
echo "  Bins    : ${FINAL_BINS_DIR}"
echo "  Output  : ${GTDBTK_OUT}"
echo "  Database: ${DATABASES}/gtdbtk/release232"

conda activate gtdbtk-env

export GTDBTK_DATA_PATH="/rds/projects/h/hallly-microbiome/software/databases/gtdbtk/release232/release232"

# Remove stale GTDB-Tk output before reclassifying
if [ -d "${GTDBTK_OUT}" ]; then
    echo "$(date): Removing stale GTDB-Tk output directory..."
    rm -rf "${GTDBTK_OUT}"
fi
mkdir -p "${GTDBTK_OUT}"

conda activate gtdbtk-env

# -- Helper function: Strip taxonomy suffixes ------------------------------
# Strips any previously appended species/taxonomy suffix from a bin name/stem,
# ensuring we always work with the true original bin identifier.
get_pristine_name() {
    local name="$1"
    name="${name%.fa}"
    if [[ "$name" =~ \.([a-zA-Z0-9]+)_ ]]; then
        echo "$name" | sed -E 's/\.([a-zA-Z0-9]+)_.*/.\1/'
    elif [[ "$name" =~ (_bin[0-9]+)_ ]]; then
        echo "$name" | sed -E 's/(_bin[0-9]+)_.*/\1/'
    else
        echo "$name"
    fi
}

# -- Idempotent rename support ---------------------------------------------
PRISTINE_BINS_DIR="${FINAL_BINS_DIR}/.pristine_bins"

if [ ! -d "${PRISTINE_BINS_DIR}" ]; then
    echo "$(date): First run detected for ${SAMPLE} - archiving pristine bin names."
    mkdir -p "${PRISTINE_BINS_DIR}"
    for f in "${FINAL_BINS_DIR}"/*.fa; do
        [ -f "$f" ] || continue
        base=$(basename "$f")
        pristine_base=$(get_pristine_name "$base")
        cp "$f" "${PRISTINE_BINS_DIR}/${pristine_base}.fa"
    done
else
    # Sanitize existing archived files if an earlier run captured already-renamed bins
    echo "$(date): Sanitizing existing pristine bin archive..."
    for f in "${PRISTINE_BINS_DIR}"/*.fa; do
        [ -f "$f" ] || continue
        base=$(basename "$f")
        pristine_base=$(get_pristine_name "$base")
        if [ "$base" != "${pristine_base}.fa" ]; then
            mv "$f" "${PRISTINE_BINS_DIR}/${pristine_base}.fa"
        fi
    done
fi

echo "$(date): Restoring pristine bin names before classification..."
find "${FINAL_BINS_DIR}" -maxdepth 1 -name '*.fa' -delete
cp "${PRISTINE_BINS_DIR}"/*.fa "${FINAL_BINS_DIR}/"

# -- CPU cap ----------------------------------------------------------------
if [ "${NUM_THREADS}" -lt 8 ]; then
    GTDBTK_CPUS="${NUM_THREADS}"
else
    GTDBTK_CPUS=8
fi

# -- Thread Containment ----------------------------------------------------
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export HMMER_NCPU=1

echo "$(date): Using ${GTDBTK_CPUS} CPU(s) for GTDB-Tk (requested ${NUM_THREADS}, capped at 8)."

gtdbtk classify_wf \
    --cpus ${GTDBTK_CPUS} \
    --extension .fa \
    --genome_dir "${FINAL_BINS_DIR}" \
    --out_dir "${GTDBTK_OUT}"

conda deactivate || true

echo "$(date): [STAGE 08/24] GTDB-Tk classification complete."

# -- Bin renaming --------------------------------------------------------------
echo "$(date): [STAGE 08/24] Renaming bins with taxonomic information..."

echo "# References for MetaQUAST --references-list option"  > "${METAQUAST_REFS_FILE}"
echo "# Generated from GTDB-Tk classifications on $(date)" >> "${METAQUAST_REFS_FILE}"

declare -A bin_rename_map

process_gtdbtk_summary() {
    local summary_file="$1"
    local domain="$2"

    echo "$(date): Processing ${domain} taxonomy from ${summary_file}..."

    while IFS=$'\t' read -r user_genome classification other; do
        [ "$user_genome" = "user_genome" ] && continue

        species=$(echo "$classification" | rev | cut -d';' -f1 | rev | \
            sed 's/^s__//' | tr ' ' '_')
        genus=$(echo "$classification" | rev | cut -d';' -f2 | rev | \
            sed 's/^g__//' | tr ' ' '_')

        if [ "$genus" = "N/A" ] && [ "$species" = "N/A" ]; then
            genus="Unknown"
            species="organism"
        elif [ -z "$species" ]; then
            species="${genus}_sp"
        fi

        # Strip any species extension attached to user_genome to ensure single appending
        pristine_user_genome=$(get_pristine_name "$user_genome")
        old_path="${FINAL_BINS_DIR}/${user_genome}.fa"

        if [[ "$pristine_user_genome" == "${SAMPLE}_"* ]]; then
            new_bin_name="${pristine_user_genome}_${species}"
        else
            new_bin_name="${SAMPLE}_${pristine_user_genome}_${species}"
        fi
        new_path="${FINAL_BINS_DIR}/${new_bin_name}.fa"

        if [ -f "$old_path" ]; then
            mv "$old_path" "$new_path"
            echo "$(date): Renamed: $(basename "$old_path") -> $(basename "$new_path")"

            bin_rename_map["${pristine_user_genome}"]="${new_bin_name}"

            if [[ "$genus" != "Unknown" && "$species" != *"_sp" ]]; then
                species_clean="${species#${genus}_}"
                if [ "$species_clean" != "$species" ]; then
                    metaquast_name="${genus} ${species_clean}"
                else
                    metaquast_name="${genus} ${species}"
                fi
                echo "${metaquast_name}" | tr '_' ' ' >> "${METAQUAST_REFS_FILE}"
            fi
        else
            echo "$(date): WARNING - Bin file not found, skipping: ${old_path}"
        fi

    done < "$summary_file"
}

update_stats_file() {
    local stats_file="${FINAL_BINS_DIR}/reassembled_bins.stats"

    if [ ! -f "$stats_file" ]; then
        echo "$(date): WARNING - Stats file not found, skipping update: ${stats_file}"
        return
    fi

    if [ ${#bin_rename_map[@]} -eq 0 ]; then
        echo "$(date): WARNING - No bin renames recorded, skipping stats file update."
        return
    fi

    echo "$(date): Updating bin names in stats file: ${stats_file}"

    local tmp_file="${stats_file}.tmp"
    > "$tmp_file"

    while IFS=$'\t' read -r bin_id rest; do
        if [[ "$bin_id" == "#"* ]] || [[ "$bin_id" == "bin"* && "$rest" == *"length"* ]]; then
            if [ -n "$rest" ]; then
                echo -e "${bin_id}\t${rest}" >> "$tmp_file"
            else
                echo -e "${bin_id}" >> "$tmp_file"
            fi
            continue
        fi

        clean_key=$(get_pristine_name "$bin_id")

        if [ -n "${bin_rename_map[$clean_key]}" ]; then
            new_name="${bin_rename_map[$clean_key]}"
            if [ -n "$rest" ]; then
                echo -e "${new_name}\t${rest}" >> "$tmp_file"
            else
                echo -e "${new_name}" >> "$tmp_file"
            fi
        else
            if [ -n "$rest" ]; then
                echo -e "${bin_id}\t${rest}" >> "$tmp_file"
            else
                echo -e "${bin_id}" >> "$tmp_file"
            fi
        fi
    done < "$stats_file"

    mv "$tmp_file" "$stats_file"
    echo "$(date): Stats file updated: ${stats_file}"
}

if [ -f "${GTDBTK_OUT}/classify/gtdbtk.bac120.summary.tsv" ]; then
    process_gtdbtk_summary "${GTDBTK_OUT}/classify/gtdbtk.bac120.summary.tsv" "bacterial"
else
    echo "$(date): WARNING - No bacterial summary file found, skipping."
fi

if [ -f "${GTDBTK_OUT}/classify/gtdbtk.ar53.summary.tsv" ]; then
    process_gtdbtk_summary "${GTDBTK_OUT}/classify/gtdbtk.ar53.summary.tsv" "archaeal"
else
    echo "$(date): WARNING - No archaeal summary file found, skipping."
fi

update_stats_file

echo "$(date): Bin renaming complete."
echo "$(date): MetaQUAST references written to: ${METAQUAST_REFS_FILE}"

gzip -vf "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_1.fastq" || true
gzip -vf "${KNEADDATA_OUT}/out.${SAMPLE}.R1_kneaddata_paired_2.fastq" || true

echo "STAGE_09" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 08/24] GTDB-Tk complete."