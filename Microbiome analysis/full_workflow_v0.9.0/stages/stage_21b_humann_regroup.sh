#!/bin/bash

################################################################################
# STAGE 21B: HUMAnN3 regroup - CPM renorm + regroup/rename to rxn, EC, KO
################################################################################
#
# Purpose:
#   Takes the raw (unnormalised) gene family tables already produced by
#   STAGE_20 (whole-sample) and STAGE_21 (per-bin), and for each:
#     1. Renormalises to CPM (humann_renorm_table --units cpm)
#     2. Regroups the CPM table from UniRef IDs to three separate feature
#        spaces: MetaCyc reactions (rxn), EC numbers (level4ec), KEGG
#        orthologs (ko)
#     3. Renames each regrouped table's IDs to human-readable descriptions
#
#   This is a distinct, CPM-based normalisation from the relative-abundance
#   (--units relab) one STAGE_20/21 already produce for the stacked-bar
#   visualizations -- both are kept; this stage doesn't touch or replace the
#   existing *_normalised.tsv files.
#
# UniRef version auto-detection:
#   humann_regroup_table needs to know whether your gene families are
#   UniRef50 or UniRef90 IDs (the --groups value is literally prefixed
#   uniref50_* or uniref90_*) -- this depends on which UniRef diamond
#   database HUMAnN3 searched against, which this script doesn't otherwise
#   know. Rather than hard-code one, each CPM table is sniffed for its own
#   ID prefix and the matching --groups value is picked automatically. If a
#   bin's gene families are ALL UNMAPPED/UNGROUPED (no UniRef IDs at all,
#   e.g. a very small or very divergent bin), regroup/rename is skipped for
#   that file with a warning -- the CPM table is still kept.
#
# --names / --groups VALUES:
#   Confirmed against `humann_rename_table --help` / `humann_regroup_table --help`
#   on this install: "ec", "kegg-orthology", and uniref{50,90}_level4ec are all
#   valid. No further verification needed.
#
# Input variables (from SLURM template):
#   - SAMPLE             : Sample name
#   - NUM_THREADS         : Number of CPU threads
#   - HUMANN_READS_DIR    : ${DIRECTORY}/humann_reads_output/${SAMPLE} (STAGE_20 output)
#   - HUMANN_CONTIGS_DIR  : ${DIRECTORY}/humann_contigs_output/${SAMPLE} (STAGE_21 output)
#   - STATUS_DIR          : Path to status directory
#
# Output (per input file: whole-sample, and each STAGE_21 bin with status=ok):
#   - <basename>_genefamilies-cpm.tsv
#   - <basename>_rxn-cpm.tsv          (+ <basename>_rxn-cpm-named.tsv)
#   - <basename>_ec-cpm.tsv           (+ <basename>_ec-cpm-named.tsv)
#   - <basename>_ko-cpm.tsv           (+ <basename>_ko-cpm-named.tsv)
#   - Completion marker : ${HUMANN_CONTIGS_DIR}/stage21b_regroup_completed.txt
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_22
#
################################################################################

set -uo pipefail

echo "$(date): [STAGE 21B/24] Starting HUMAnN3 CPM renorm + regroup/rename for ${SAMPLE}..."

declare -A REGROUP_SUFFIX=( [rxn]="rxn" [ec]="level4ec" [ko]="ko" )
declare -A RENAME_ARG=( [rxn]="metacyc-rxn" [ec]="ec" [ko]="kegg-orthology" )

# Renorm + regroup + rename one raw genefamilies.tsv. $1 = input path,
# $2 = output prefix (files written as ${2}-cpm.tsv, ${2}_rxn-cpm.tsv, etc).
process_genefamilies() {
    local input_file="$1"
    local out_prefix="$2"
    local label="$3"   # for logging only

    if [ ! -s "${input_file}" ]; then
        echo "$(date): [STAGE 21B/24]   WARNING: ${label} genefamilies file not found (${input_file}), skipping."
        return 0
    fi

    local cpm_file="${out_prefix}_genefamilies-cpm.tsv"

    humann_renorm_table \
        --input "${input_file}" \
        --output "${cpm_file}" \
        --units cpm \
        --update-snames

    # Sniff the UniRef build from the CPM table's own feature IDs -- skip
    # rows like UNMAPPED/UNGROUPED which never carry a UniRef prefix.
    local uniref_ver
    uniref_ver=$(grep -m1 -oE '^UniRef(50|90)_' "${cpm_file}" | grep -oE '(50|90)' || true)

    if [ -z "${uniref_ver}" ]; then
        echo "$(date): [STAGE 21B/24]   WARNING: no UniRef50/90 IDs found in ${cpm_file} (all UNMAPPED/UNGROUPED?); skipping regroup/rename for ${label}."
        return 0
    fi

    echo "$(date): [STAGE 21B/24]   ${label}: detected UniRef${uniref_ver}, regrouping to rxn/ec/ko..."

    for group in rxn ec ko; do
        local groups_arg="uniref${uniref_ver}_${REGROUP_SUFFIX[$group]}"
        local regrouped_file="${out_prefix}_${group}-cpm.tsv"
        local named_file="${out_prefix}_${group}-cpm-named.tsv"

        humann_regroup_table \
            --input "${cpm_file}" \
            --output "${regrouped_file}" \
            --groups "${groups_arg}"

        humann_rename_table \
            --input "${regrouped_file}" \
            --output "${named_file}" \
            --names "${RENAME_ARG[$group]}"
    done
}

conda activate humann3-env

# -- Whole-sample (STAGE_20) -------------------------------------------------
WHOLE_SAMPLE_INPUT="${HUMANN_READS_DIR}/out.${SAMPLE}_combined_genefamilies.tsv"
WHOLE_SAMPLE_PREFIX="${HUMANN_READS_DIR}/out.${SAMPLE}_combined"

echo "$(date): [STAGE 21B/24] Processing whole-sample gene families..."
process_genefamilies "${WHOLE_SAMPLE_INPUT}" "${WHOLE_SAMPLE_PREFIX}" "whole-sample"

# -- Per-bin (STAGE_21) -------------------------------------------------------
MANIFEST="${HUMANN_CONTIGS_DIR}/${SAMPLE}_bin_humann_manifest.tsv"

if [ ! -s "${MANIFEST}" ]; then
    echo "$(date): [STAGE 21B/24] WARNING: STAGE_21 manifest not found (${MANIFEST}); skipping per-bin regroup."
else
    echo "$(date): [STAGE 21B/24] Processing per-bin gene families..."
    # Read bins with status "ok" (skip header, NOREADS/failed bins)
    while IFS=$'\t' read -r BIN_LABEL GTDBTK_SPECIES N_READS STATUS; do
        [ "${BIN_LABEL}" == "bin" ] && continue   # header row
        [ "${STATUS}" != "ok" ] && continue

        BIN_INPUT="${HUMANN_CONTIGS_DIR}/${BIN_LABEL}_genefamilies.tsv"
        BIN_PREFIX="${HUMANN_CONTIGS_DIR}/${BIN_LABEL}"

        process_genefamilies "${BIN_INPUT}" "${BIN_PREFIX}" "${BIN_LABEL}"
    done < "${MANIFEST}"
fi

conda deactivate || true

# Write completion marker used by checkpoint scanner. Deliberately separate
# from STAGE_21's own humann_loop_completed.txt so the two stages' progress
# can be told apart on disk.
touch "${HUMANN_CONTIGS_DIR}/stage21b_regroup_completed.txt"

echo "STAGE_22" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 21B/24] HUMAnN3 CPM renorm + regroup/rename complete for ${SAMPLE}."