#!/bin/bash
################################################################################
# STAGE 18B: RAxML-NG - Phylogenetic tree per Panaroo/Snippy taxon group
################################################################################
#
# For every taxon group STAGE_18 (Snippy) or STAGE_17 (Panaroo, for
# genus-wide taxa) produced a usable multiple-sequence alignment for, this
# builds an ML tree with raxml-ng, batch-wide, using the same last-arrival
# pattern as STAGE_17/STAGE_18. Taxa are discovered by scanning
# ${PANAROO_OUT} (every taxon that got a pangenome has a folder there); each
# taxon's alignment then comes from one of two places depending on
# USE_PANAROO_ALIGNMENT below:
#   - Species-level taxa (the default): ${SNIPPY_OUT}/<taxon_key>/core.aln,
#     Snippy's core-SNP alignment.
#   - Genus-wide taxa in USE_PANAROO_ALIGNMENT: the core-gene alignment
#     Panaroo itself wrote (path recorded in
#     ${PANAROO_OUT}/<taxon_key>/core_alignment_path.txt by STAGE_17), since
#     Snippy's reference-coordinate mapping isn't appropriate once
#     divergence crosses species boundaries — see STAGE_17's header.
#
# MODEL CHOICE: Snippy's core.aln comes from `snp-sites -c` inside
# snippy-core, which strips invariant columns and keeps only variable
# sites. Feeding that straight to a plain substitution model overestimates
# branch lengths, since the model doesn't know how many invariant sites
# were removed. This is corrected with raxml-ng's Lewis ascertainment-bias
# correction (+ASC_LEWIS), which doesn't require knowing the exact
# invariant-site count (unlike +ASC_STAM).
#
# Panaroo's core-gene alignment is NOT SNP-only -- it's a normal multiple
# sequence alignment of orthologous genes, invariant sites included. Using
# +ASC_LEWIS on it would be wrong: that correction specifically assumes
# invariant sites were already removed, so applying it to an alignment that
# still has them double-corrects and miscalibrates branch lengths. MODEL is
# used for Snippy-derived (SNP-only) alignments; MODEL_GENE_ALIGNMENT (a
# plain model, no ASC correction) is used for Panaroo-derived (full gene)
# alignments instead.
#
# MIN_TAXA=3: a 2-taxon "tree" is just a single branch with no topology to
# speak of, so groups below this are skipped with a warning rather than
# handed to raxml-ng.
#
# Runs --all (ML search + bootstrapping) with autoMRE bootstopping so each
# tree comes with branch support without you having to guess a replicate
# count up front — this is slower than a plain --search and can be swapped
# out (see RAXML_MODE below) if you just want a fast topology first.
#
# Input variables (from SLURM template):
#   - SAMPLE        : Sample name (current array task)
#   - DIRECTORY     : Workflow working directory
#   - NUM_THREADS   : Number of CPU threads
#   - STATUS_DIR    : Path to status directory
#   - SNIPPY_OUT    : Snippy output directory (STAGE_18), source of core.aln
#   - RAXML_OUT     : Output directory for this stage
#
# Output:
#   - ${RAXML_OUT}/<taxon_key>/tree.raxml.bestTree   - ML tree (Newick)
#   - ${RAXML_OUT}/<taxon_key>/tree.raxml.support     - tree with bootstrap support
#   - ${RAXML_OUT}/<taxon_key>/tree.raxml.*           - remaining raxml-ng output
#   - ${RAXML_OUT}/<taxon_key>/itol/                  - iTOL annotation files (see
#       below); drag the tree file into iTOL's viewer, then drag these
#       alongside it — iTOL annotations are always separate files, never
#       embedded in the Newick itself:
#       * <taxon_key>_sample_colorstrip.txt   - leaf color by sample of origin
#       * <taxon_key>_reference_highlight.txt - bolds/highlights whichever
#           leaf is the type-strain/in-group reference row STAGE_18 labelled
#   - ${STATUS_DIR}/raxml_complete.txt                - shared batch completion flag
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_19
#
################################################################################

echo "$(date): [STAGE 18B/24] RAxML-NG check-in for ${SAMPLE}..."

# -- Configuration ----------------------------------------------------------
MODEL="GTR+G+ASC_LEWIS"    # see MODEL CHOICE note above (Snippy/SNP-only alignments)
MODEL_GENE_ALIGNMENT="GTR+G"   # see MODEL CHOICE note above (Panaroo/full gene alignments)
MIN_TAXA=3                 # skip groups with fewer sequences than this
RAXML_MODE="--all"         # --all = ML search + bootstrap; --search = topology only (fast)
BS_TREES="autoMRE{1000}"   # only used with --all
SEED=42                    # fixed for reproducibility; unset for a random seed each run

# Taxa in here were built directly from Panaroo's own core-gene alignment
# (STAGE_17) rather than Snippy (STAGE_18) -- see STAGE_17's header for why.
# Kept in sync with STAGE_17/STAGE_18's own USE_PANAROO_ALIGNMENT list.
declare -A USE_PANAROO_ALIGNMENT=(
    ["Segatella_genus"]=true
    ["Segatella_genus_relaxed"]=true
)

if [ -f "${STATUS_DIR}/raxml_complete.txt" ]; then
    echo "$(date): RAxML-NG already complete for this batch — skipping."
else
    touch "${STATUS_DIR}/stage18b_arrived.${SLURM_ARRAY_JOB_ID}.${SAMPLE}"
    TOTAL_SAMPLES=$(wc -l < "${DIRECTORY}/samples.txt")
    ARRIVED=$(find "${STATUS_DIR}" -maxdepth 1 -name "stage18b_arrived.${SLURM_ARRAY_JOB_ID}.*" | wc -l)
    echo "$(date): ${ARRIVED}/${TOTAL_SAMPLES} samples have reached STAGE_18B."

    if [ "$ARRIVED" -ge "$TOTAL_SAMPLES" ] && mkdir "${STATUS_DIR}/stage18b_raxml.lock.${SLURM_ARRAY_JOB_ID}" 2>/dev/null; then
        echo "$(date): ${SAMPLE} is the last arrival — launching RAxML-NG runs..."

        conda activate raxml-ng2-env
        FAILED_TAXA=()

        for taxon_dir in "${PANAROO_OUT}"/*/; do
            [ -d "$taxon_dir" ] || continue
            taxon_key=$(basename "$taxon_dir")

            if [ "${USE_PANAROO_ALIGNMENT[$taxon_key]:-false}" = true ]; then
                ALN_POINTER="${taxon_dir}core_alignment_path.txt"
                if [ ! -f "$ALN_POINTER" ]; then
                    echo "$(date): ${taxon_key} - uses Panaroo's own alignment but ${ALN_POINTER} doesn't exist (STAGE_17 may not have found one), skipping."
                    continue
                fi
                CORE_ALN=$(cat "$ALN_POINTER")
                if [ ! -f "$CORE_ALN" ]; then
                    echo "$(date): ${taxon_key} - alignment file recorded in ${ALN_POINTER} doesn't exist (${CORE_ALN}), skipping."
                    continue
                fi
                THIS_MODEL="$MODEL_GENE_ALIGNMENT"
                echo "$(date): ${taxon_key} - using Panaroo's own core-gene alignment: ${CORE_ALN} (model: ${THIS_MODEL})"
            else
                CORE_ALN="${SNIPPY_OUT}/${taxon_key}/core.aln"
                if [ ! -f "$CORE_ALN" ]; then
                    echo "$(date): ${taxon_key} - no core.aln (snippy-core didn't complete for this taxon), skipping."
                    continue
                fi
                THIS_MODEL="$MODEL"
            fi

            N_TAXA=$(grep -c '^>' "$CORE_ALN")
            if [ "$N_TAXA" -lt "$MIN_TAXA" ]; then
                echo "$(date): ${taxon_key} - only ${N_TAXA} sequence(s) in core.aln, below MIN_TAXA=${MIN_TAXA}; skipping."
                continue
            fi

            TAXON_RAXML_OUT="${RAXML_OUT}/${taxon_key}"
            mkdir -p "$TAXON_RAXML_OUT"

            RAXML_ARGS=(
                $RAXML_MODE
                --msa "$CORE_ALN"
                --model "$THIS_MODEL"
                --prefix "${TAXON_RAXML_OUT}/tree"
                --threads "${NUM_THREADS}"
                --force perf_threads
                --redo
            )
            [ "$RAXML_MODE" = "--all" ] && RAXML_ARGS+=(--bs-trees "$BS_TREES")
            [ -n "$SEED" ] && RAXML_ARGS+=(--seed "$SEED")

            echo "$(date): ${taxon_key} - building tree from ${N_TAXA} taxa..."
            if raxml-ng "${RAXML_ARGS[@]}"; then
                echo "$(date): ${taxon_key} - RAxML-NG complete."

                # -- iTOL annotation files --------------------------------
                # iTOL annotations are always separate plain-text files
                # dragged onto the tree in the viewer, never embedded in the
                # Newick itself, so these live alongside the tree output.
                ITOL_DIR="${TAXON_RAXML_OUT}/itol"
                mkdir -p "$ITOL_DIR"

                mapfile -t LEAVES < <(grep '^>' "$CORE_ALN" | sed 's/^>//')

                unset LEAF_CATEGORY CAT_COLOR
                declare -A LEAF_CATEGORY CAT_COLOR
                REF_LEAF=""

                for leaf in "${LEAVES[@]}"; do
                    if [[ "$leaf" == TYPE_STRAIN__* ]]; then
                        LEAF_CATEGORY["$leaf"]="Type strain"
                        REF_LEAF="$leaf"
                    elif [[ "$leaf" == INGROUP_REF__* ]]; then
                        LEAF_CATEGORY["$leaf"]="In-group reference"
                        REF_LEAF="$leaf"
                    elif [[ "$leaf" =~ ^([A-Za-z0-9]+)_bin[0-9]+\. ]]; then
                        LEAF_CATEGORY["$leaf"]="${BASH_REMATCH[1]}"
                    else
                        LEAF_CATEGORY["$leaf"]="Unknown"
                    fi
                done

                # Colorblind-friendly palette for actual samples; fixed
                # colors reserved for the non-sample categories below.
                PALETTE=(
                    "#e41a1c" "#377eb8" "#4daf4a" "#984ea3" "#ff7f00"
                    "#ffff33" "#a65628" "#f781bf" "#66c2a5" "#fc8d62"
                )
                mapfile -t ALL_CATS < <(printf '%s\n' "${LEAF_CATEGORY[@]}" | sort -u)
                pi=0
                for cat in "${ALL_CATS[@]}"; do
                    case "$cat" in
                        "Type strain")         CAT_COLOR["$cat"]="#000000" ;;
                        "In-group reference")  CAT_COLOR["$cat"]="#4d4d4d" ;;
                        "Unknown")              CAT_COLOR["$cat"]="#cccccc" ;;
                        *)
                            CAT_COLOR["$cat"]="${PALETTE[$(( pi % ${#PALETTE[@]} ))]}"
                            pi=$((pi + 1))
                            ;;
                    esac
                done

                COLORSTRIP_FILE="${ITOL_DIR}/${taxon_key}_sample_colorstrip.txt"
                {
                    printf 'DATASET_COLORSTRIP\n'
                    printf 'SEPARATOR TAB\n'
                    printf 'DATASET_LABEL\tSample of origin\n'
                    printf 'COLOR\t#000000\n'
                    printf 'LEGEND_TITLE\tSample / reference\n'
                    printf 'LEGEND_SHAPES'
                    for _ in "${ALL_CATS[@]}"; do printf '\t1'; done
                    printf '\n'
                    printf 'LEGEND_COLORS'
                    for cat in "${ALL_CATS[@]}"; do printf '\t%s' "${CAT_COLOR[$cat]}"; done
                    printf '\n'
                    printf 'LEGEND_LABELS'
                    for cat in "${ALL_CATS[@]}"; do printf '\t%s' "$cat"; done
                    printf '\n'
                    printf 'DATA\n'
                    for leaf in "${LEAVES[@]}"; do
                        printf '%s\t%s\t%s\n' "$leaf" "${CAT_COLOR[${LEAF_CATEGORY[$leaf]}]}" "${LEAF_CATEGORY[$leaf]}"
                    done
                } > "$COLORSTRIP_FILE"

                REFHIGHLIGHT_FILE="${ITOL_DIR}/${taxon_key}_reference_highlight.txt"
                {
                    printf 'TREE_COLORS\n'
                    printf 'SEPARATOR TAB\n'
                    printf 'DATA\n'
                    if [ -n "$REF_LEAF" ]; then
                        printf '%s\tlabel\t#000000\tbold\t1.5\n' "$REF_LEAF"
                        printf '%s\tlabel_background\t#ffff99\n' "$REF_LEAF"
                    fi
                } > "$REFHIGHLIGHT_FILE"

                echo "$(date): ${taxon_key} - iTOL annotation files written to ${ITOL_DIR}."
            else
                echo "$(date): WARNING - RAxML-NG failed for ${taxon_key}; skipping this taxon, continuing with the rest." >&2
                FAILED_TAXA+=("$taxon_key")
            fi
        done

        conda deactivate || true

        # A single taxon's failure no longer fails the whole batch — it's
        # logged and recorded here so the rest of the batch can complete
        # (and so the next submission doesn't redo everything just because
        # one taxon hit an edge case). Check this file after a run.
        if [ "${#FAILED_TAXA[@]}" -gt 0 ]; then
            printf '%s\n' "${FAILED_TAXA[@]}" > "${STATUS_DIR}/raxml_failed_taxa.txt"
            echo "$(date): WARNING - RAxML-NG failed for: ${FAILED_TAXA[*]} (see ${STATUS_DIR}/raxml_failed_taxa.txt)." >&2
        else
            rm -f "${STATUS_DIR}/raxml_failed_taxa.txt"
        fi

        touch "${STATUS_DIR}/raxml_complete.txt"
        echo "$(date): [STAGE 18B/24] RAxML-NG batch run complete."
    else
        echo "$(date): Not the last arrival (or lock already claimed) — deferring RAxML-NG."
    fi
fi

echo "STAGE_19" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 18B/24] Done for ${SAMPLE}."