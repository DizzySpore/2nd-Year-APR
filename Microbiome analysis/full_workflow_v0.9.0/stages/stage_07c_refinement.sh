#!/bin/bash

################################################################################
# STAGE 07C (Binette variant): Bin refinement using Binette instead of the
#             metaWRAP-style hybrid generation + CheckM1/CheckM2 pipeline
################################################################################
#
# Input variables (from SLURM template):
#   - SAMPLE        : Sample name
#   - NUM_THREADS   : Number of CPU threads
#   - MEMORY        : Memory in GB
#   - BINNING_DIR   : Directory containing per-binner bin subdirectories
#   - REFINEMENT_DIR: Output directory for refined bins
#   - ASSEMBLY_FILE : Path to MEGAHIT final contigs (same variable used by
#                      STAGE 07A / 07B). Binette requires this via --contigs;
#                      the old metaWRAP-style pipeline did not need it here.
#   - COMPLETENESS  : Minimum completeness threshold (%)
#   - CONTAMINATION : Maximum contamination threshold (%)
#   - DATABASES     : Path to reference databases
#   - STATUS_DIR    : Path to status directory
#
# Output:
#   - Refined bins     : ${REFINEMENT_DIR}/binette_${COMPLETENESS}_${CONTAMINATION}_bins/*.fa
#   - Bin stats        : ${REFINEMENT_DIR}/binette_${COMPLETENESS}_${CONTAMINATION}_bins.stats
#   - Work files       : ${REFINEMENT_DIR}/work_files/  (binsO/ populated so
#                         STAGE 07D needs no changes)
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_07E
#
# Note:
#   - Binette is run from the binette-env conda environment.
#   - Binette builds its own hybrid (intersection/difference/union) bins
#     from the input bin sets and scores everything with CheckM2 in one
#     pass, so the separate hybrid-generation, per-combo CheckM1/CheckM2,
#     and select_best_bins.py steps from the metaWRAP-style script are no
#     longer needed.
#   - Binette's own --min_completeness/--max_contamination filtering
#     replaces the old post-dereplication threshold filtering step.
#   - Output layout is deliberately kept compatible with STAGE 07D, which
#     only reads *.fa files from ${REFINEMENT_DIR}/work_files/binsO/.
#
################################################################################

echo "$(date): [STAGE 07C - Binette] Starting bin refinement for ${SAMPLE}..."
echo "  Binning dir  : ${BINNING_DIR}"
echo "  Contigs      : ${ASSEMBLY_FILE}"
echo "  Output       : ${REFINEMENT_DIR}"
echo "  Completeness : ${COMPLETENESS}%"
echo "  Contamination: ${CONTAMINATION}%"

if [ -z "${ASSEMBLY_FILE}" ] || [ ! -f "${ASSEMBLY_FILE}" ]; then
    echo "$(date): ERROR - ASSEMBLY_FILE assembly FASTA not found: ${ASSEMBLY_FILE}"
    echo "$(date): Binette requires the full assembly FASTA via --contigs."
    exit 1
fi

# Remove stale refinement output
if [ -d "${REFINEMENT_DIR}" ]; then
    echo "$(date): Removing stale refinement directory..."
    rm -rf "${REFINEMENT_DIR}"
fi
mkdir -p "${REFINEMENT_DIR}"
cd "${REFINEMENT_DIR}"

################################################################################
# STEP 1: COPY AND FILTER BINS INTO WORKING DIRECTORIES
################################################################################
echo "$(date): [STAGE 07C - Binette] Copying and filtering bins into working directories..."

declare -A BINNER_MAP=(
    ["binsA"]="${BINNING_DIR}/concoct_bins"
    ["binsB"]="${BINNING_DIR}/metabat2_bins"
    ["binsC"]="${BINNING_DIR}/comebin_bins"
    ["binsD"]="${BINNING_DIR}/semibin2_bins"
)

for BIN_SET in binsA binsB binsC binsD; do
    SRC="${BINNER_MAP[$BIN_SET]}"
    mkdir -p "${BIN_SET}"

    for bin_file in "${SRC}/"*.fa; do
        [ -f "$bin_file" ] || continue
        [[ "$(basename $bin_file)" == *"unbinned"* ]] && continue
        [[ "$(basename $bin_file)" == *"tooShort"* ]] && continue
        [[ "$(basename $bin_file)" == *"lowDepth"* ]] && continue

        SIZE=$(stat -c%s "$bin_file")
        if (( SIZE > 50000 )) && (( SIZE < 20000000 )); then
            BASE=$(basename "$bin_file")
            cp "$bin_file" "${BIN_SET}/${BASE%.*}.fa"
        else
            echo "$(date): Skipping $(basename $bin_file) — size outside 50kb-20Mb range"
        fi
    done

    COUNT=$(ls "${BIN_SET}/"*.fa 2>/dev/null | wc -l)
    echo "$(date): ${BIN_SET} ($(basename ${SRC})): ${COUNT} bins"
done

# Only pass bin-set directories that actually contain bins to Binette
BIN_DIR_ARGS=()
for BIN_SET in binsA binsB binsC binsD; do
    COUNT=$(ls "${BIN_SET}/"*.fa 2>/dev/null | wc -l)
    if [ "${COUNT}" -gt 0 ]; then
        BIN_DIR_ARGS+=("${BIN_SET}")
    else
        echo "$(date): ${BIN_SET} has no bins after filtering — excluding from Binette run"
    fi
done

if [ "${#BIN_DIR_ARGS[@]}" -lt 2 ]; then
    echo "$(date): ERROR - Binette needs at least two non-empty bin sets, found ${#BIN_DIR_ARGS[@]}"
    exit 1
fi

################################################################################
# STEP 2: RUN BINETTE
################################################################################
echo "$(date): [STAGE 07C - Binette] Running Binette on bin sets: ${BIN_DIR_ARGS[*]}"

conda activate binette-env

BINETTE_OUT="${REFINEMENT_DIR}/binette_out"

binette \
    --bin_dirs "${BIN_DIR_ARGS[@]}" \
    --contigs "${ASSEMBLY_FILE}" \
    --outdir "${BINETTE_OUT}" \
    --prefix "${SAMPLE}" \
    --threads "${NUM_THREADS}" \
    --min_completeness "${COMPLETENESS}" \
    --max_contamination "${CONTAMINATION}" \
    --checkm2_db "${DATABASES}/checkm2/CheckM2_database/uniref100.KO.1.dmnd"

BINETTE_STATUS=$?

conda deactivate || true

if [ ${BINETTE_STATUS} -ne 0 ] || [ ! -s "${BINETTE_OUT}/final_bins_quality_reports.tsv" ]; then
    echo "$(date): ERROR - Binette failed to produce output for ${SAMPLE}"
    exit 1
fi

FINAL_BIN_COUNT=$(ls "${BINETTE_OUT}/final_bins/"*.fa 2>/dev/null | wc -l)
echo "$(date): Binette produced ${FINAL_BIN_COUNT} final bins passing thresholds"

################################################################################
# STEP 3: ORGANISE OUTPUTS (kept compatible with STAGE 07D)
################################################################################
echo "$(date): [STAGE 07C - Binette] Organising output files..."

mkdir -p work_files
mv binsA binsB binsC binsD work_files/ 2>/dev/null || true
mv "${BINETTE_OUT}" work_files/ 2>/dev/null || true
BINETTE_OUT="work_files/$(basename "${BINETTE_OUT}")"

# STAGE 07D reads bins from ${REFINEMENT_DIR}/work_files/binsO — populate it
# from Binette's final_bins/ so no downstream changes are needed.
mkdir -p work_files/binsO
cp "${BINETTE_OUT}/final_bins/"*.fa work_files/binsO/ 2>/dev/null

FINAL_BINS="${REFINEMENT_DIR}/binette_${COMPLETENESS}_${CONTAMINATION}_bins"
mkdir -p "${FINAL_BINS}"

for bin_file in work_files/binsO/*.fa; do
    [ -f "$bin_file" ] || continue
    SIZE=$(stat -c%s "$bin_file")
    if [ "$SIZE" -gt 0 ]; then
        cp "$bin_file" "${FINAL_BINS}/"
    else
        echo "$(date): Skipping empty bin $(basename $bin_file)"
    fi
done

# Keep Binette's native quality report as the stats file (columns differ from
# the old metaWRAP-style stats: name, origin, is_original, original_name,
# completeness, contamination, checkm2_model, score, size, N50,
# coding_density, contig_count)
cp "${BINETTE_OUT}/final_bins_quality_reports.tsv" "${FINAL_BINS}.stats"

# Generate contig membership file for downstream use (same format as before)
for bin_file in "${FINAL_BINS}/"*.fa; do
    [ -f "$bin_file" ] || continue
    BIN_NAME=$(basename "$bin_file" .fa)
    grep ">" "$bin_file" | sed "s/>//" | \
        awk -v bin="$BIN_NAME" '{print $1"\t"bin}'
done > "${FINAL_BINS}.contigs"

echo "$(date): [STAGE 07C - Binette] Summary for ${SAMPLE}:"
echo "  Final bins : ${FINAL_BIN_COUNT}"
echo "  Output dir : ${FINAL_BINS}"
echo "  Stats file : ${FINAL_BINS}.stats"

echo "STAGE_07E" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 07C - Binette] Bin refinement complete."
