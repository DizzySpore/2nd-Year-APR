#!/bin/bash
################################################################################
# STAGE 18: Snippy - within-group SNP calling on Panaroo taxon groups
################################################################################
#
# For every taxon group STAGE_17 (Panaroo) successfully built a pangenome
# for, this calls variants of every group member against a single reference
# using snippy's assembly mode (--ctgs — there's no per-bin read data to
# align, only per-sample raw reads from before binning), then merges the
# per-member snippy runs into a core-genome SNP alignment with snippy-core.
#
# ASSUMPTION: "SNP calling between the pangenome groups" is read as WITHIN
# each pangenome group (one within-group core-SNP alignment per taxon, same
# shape as the pangenome analysis) rather than a cross-species comparison —
# SNP/core-alignment calling across genomes as divergent as different
# genera/species isn't a standard or meaningful analysis, since positions
# don't correspond 1:1 across that much divergence. If a cross-group
# comparison was actually intended, this needs a different approach entirely
# (e.g. ANI/MLST-style distance rather than SNP calling) — flag it if so.
#
# Reference choice per group:
#   - Species-level taxa with a populated panaroo_references/<taxon_key>/
#     folder: the candidate with the fewest contigs is used (a rough
#     "most complete assembly" heuristic), preferring RefSeq (GCF_)
#     accessions over GenBank (GCA_) duplicates of the same assembly on a
#     tie. Every group member is called against it.
#   - Groups with no reference folder (currently just Veillonella_genus):
#     there's no external reference to use, so one of the group's OWN
#     members (again, fewest contigs) is picked as an ad-hoc internal
#     reference and excluded from being called against itself. This is a
#     real limitation, not a substitute for a proper reference — an
#     external reference should replace it as soon as one is available,
#     since an ad-hoc in-group choice biases the "core" toward whichever
#     genome happened to be picked.
#
# The type strain (or ad-hoc stand-in) is not just used as a comparison
# anchor — snippy-core includes it as a genuine row in core.aln/core.full.aln
# by design (labelled generically "Reference" by snippy-core itself), so it
# already appears as a leaf in the resulting tree. That generic label is
# rewritten here to TYPE_STRAIN__<reference filename> for species-level taxa
# with a genuine external reference, or INGROUP_REF__<bin name> for taxa
# where no external reference was available and one of the group's own bins
# had to stand in — so anyone reading the tree can immediately tell which
# leaf is the type strain and, for genuine type strains, which exact
# accession it was.
#
# Depends on STAGE_17 having written, for each successfully-pangenomed
# taxon, ${PANAROO_OUT}/<taxon_key>/group_members.txt — the exact good_bins
# fasta paths that went into that group's pangenome. SNP calling reads that
# file rather than re-deriving group membership here, so the gene-content
# and SNP analyses can never silently drift out of sync with each other.
#
# Input variables (from SLURM template):
#   - SAMPLE        : Sample name (current array task)
#   - DIRECTORY     : Workflow working directory
#   - NUM_THREADS   : Number of CPU threads
#   - STATUS_DIR    : Path to status directory
#   - PANAROO_OUT   : Panaroo output directory (STAGE_17)
#   - SNIPPY_OUT    : Output directory for this stage
#
# Output:
#   - ${SNIPPY_OUT}/<taxon_key>/<BIN_NAME>/   - per-genome snippy folder
#   - ${SNIPPY_OUT}/<taxon_key>/core.*        - snippy-core output
#       (core.aln, core.vcf, core.txt, core.tab, etc.)
#   - ${STATUS_DIR}/snippy_complete.txt       - shared batch completion flag
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_19
#
################################################################################

echo "$(date): [STAGE 18/24] Snippy check-in for ${SAMPLE}..."

# Redefined locally (not inherited from STAGE_17) since this stage can run
# in a submission where STAGE_17's own body never executes this run (e.g.
# resuming straight into STAGE_18 because STAGE_17 already checkpointed as
# complete in a prior submission).
REFERENCE_DIR="${DIRECTORY}/panaroo_references"

# Taxa in here skip Snippy entirely -- their tree is built directly from
# Panaroo's own core-gene alignment in STAGE_18B instead (see STAGE_17's
# header for why: Snippy's reference-coordinate mapping isn't appropriate
# across genus-level divergence). Keep this in sync with STAGE_17's own
# USE_PANAROO_ALIGNMENT list.
declare -A USE_PANAROO_ALIGNMENT=(
    ["Segatella_genus"]=true
    ["Segatella_genus_relaxed"]=true
)

if [ -f "${STATUS_DIR}/snippy_complete.txt" ]; then
    echo "$(date): Snippy already complete for this batch — skipping."
else
    touch "${STATUS_DIR}/stage18_arrived.${SLURM_ARRAY_JOB_ID}.${SAMPLE}"
    TOTAL_SAMPLES=$(wc -l < "${DIRECTORY}/samples.txt")
    ARRIVED=$(find "${STATUS_DIR}" -maxdepth 1 -name "stage18_arrived.${SLURM_ARRAY_JOB_ID}.*" | wc -l)
    echo "$(date): ${ARRIVED}/${TOTAL_SAMPLES} samples have reached STAGE_18."

    if [ "$ARRIVED" -ge "$TOTAL_SAMPLES" ] && mkdir "${STATUS_DIR}/stage18_snippy.lock.${SLURM_ARRAY_JOB_ID}" 2>/dev/null; then
        echo "$(date): ${SAMPLE} is the last arrival — launching Snippy runs..."

        conda activate snippy

        # Picks the fasta with the fewest contigs from a list of candidate
        # files (a rough "most complete assembly" proxy), preferring RefSeq
        # (GCF_) over GenBank (GCA_) duplicates on a tie.
        pick_best_fasta() {
            local best="" best_score=999999999 best_rank=9
            for f in "$@"; do
                [ -f "$f" ] || continue
                local n rank
                n=$(grep -c '^>' "$f")
                rank=1
                [[ "$(basename "$f")" == GCF_* ]] && rank=0
                if [ "$n" -lt "$best_score" ] || { [ "$n" -eq "$best_score" ] && [ "$rank" -lt "$best_rank" ]; }; then
                    best="$f"
                    best_score="$n"
                    best_rank="$rank"
                fi
            done
            echo "$best"
        }

        FAILED_TAXA=()

        for taxon_dir in "${PANAROO_OUT}"/*/; do
            [ -d "$taxon_dir" ] || continue
            taxon_key=$(basename "$taxon_dir")

            if [ "${USE_PANAROO_ALIGNMENT[$taxon_key]:-false}" = true ]; then
                echo "$(date): ${taxon_key} - uses Panaroo's own core-gene alignment (see STAGE_17); skipping Snippy for this taxon."
                continue
            fi

            MEMBERS_FILE="${taxon_dir}group_members.txt"

            if [ ! -f "$MEMBERS_FILE" ]; then
                echo "$(date): ${taxon_key} - no group_members.txt (no successful pangenome for this taxon), skipping."
                continue
            fi

            mapfile -t MEMBERS < "$MEMBERS_FILE"
            echo "$(date): ${taxon_key} - ${#MEMBERS[@]} member genome(s) from the pangenome run."

            REFSET_DIR="${REFERENCE_DIR}/${taxon_key}"
            REF_FASTA=""
            EXCLUDE_REF_MEMBER=""

            if [ -d "$REFSET_DIR" ]; then
                mapfile -t REF_CANDIDATES < <(find "$REFSET_DIR" -maxdepth 1 -type f \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" \))
                if [ "${#REF_CANDIDATES[@]}" -gt 0 ]; then
                    REF_FASTA=$(pick_best_fasta "${REF_CANDIDATES[@]}")
                fi
            fi

            if [ -z "$REF_FASTA" ]; then
                echo "$(date): ${taxon_key} - no external reference available; picking an in-group reference instead (see script header notes)."
                REF_FASTA=$(pick_best_fasta "${MEMBERS[@]}")
                EXCLUDE_REF_MEMBER="$REF_FASTA"
            fi

            if [ -z "$REF_FASTA" ]; then
                echo "$(date): WARNING - ${taxon_key} - could not determine any reference; skipping this taxon." >&2
                FAILED_TAXA+=("$taxon_key")
                continue
            fi

            echo "$(date): ${taxon_key} - using reference: ${REF_FASTA}"

            TAXON_SNIPPY_OUT="${SNIPPY_OUT}/${taxon_key}"
            mkdir -p "$TAXON_SNIPPY_OUT"

            for member_fasta in "${MEMBERS[@]}"; do
                [ -f "$member_fasta" ] || continue
                [ "$member_fasta" = "$EXCLUDE_REF_MEMBER" ] && continue

                BIN_NAME=$(basename "$member_fasta" .fa)
                MEMBER_OUT="${TAXON_SNIPPY_OUT}/${BIN_NAME}"

                if snippy \
                    --outdir "$MEMBER_OUT" \
                    --reference "$REF_FASTA" \
                    --ctgs "$member_fasta" \
                    --cpus "${NUM_THREADS}" \
                    --force; then
                    echo "$(date): ${taxon_key}/${BIN_NAME} - snippy complete."
                else
                    # A single member's snippy failure doesn't sink the
                    # whole taxon: snippy-core still runs below on whatever
                    # members did succeed.
                    echo "$(date): WARNING - snippy failed for ${taxon_key}/${BIN_NAME}, excluding it from this group's core alignment." >&2
                fi
            done

            mapfile -t SNIPPY_DIRS < <(find "$TAXON_SNIPPY_OUT" -mindepth 1 -maxdepth 1 -type d)
            if [ "${#SNIPPY_DIRS[@]}" -lt 1 ]; then
                echo "$(date): ${taxon_key} - no completed snippy runs, skipping snippy-core."
                FAILED_TAXA+=("$taxon_key")
                continue
            fi

            if snippy-core --ref "$REF_FASTA" --prefix "${TAXON_SNIPPY_OUT}/core" "${SNIPPY_DIRS[@]}"; then
                echo "$(date): ${taxon_key} - snippy-core complete (${#SNIPPY_DIRS[@]} genome(s) vs reference)."

                # snippy-core labels the reference row literally "Reference"
                # in both core.aln and core.full.aln, which is useless once
                # it ends up as an anonymous leaf in a tree. Rename it to
                # something traceable back to the actual file used, and
                # flag clearly whether it's a genuine type strain or just
                # one of the group's own bins pressed into service as a
                # stand-in reference (see reference-choice notes above).
                REF_LABEL=$(basename "$REF_FASTA")
                REF_LABEL="${REF_LABEL%.*}"
                if [ -n "$EXCLUDE_REF_MEMBER" ]; then
                    REF_LABEL="INGROUP_REF__${REF_LABEL}"
                else
                    REF_LABEL="TYPE_STRAIN__${REF_LABEL}"
                fi
                for aln in "${TAXON_SNIPPY_OUT}/core.aln" "${TAXON_SNIPPY_OUT}/core.full.aln"; do
                    [ -f "$aln" ] && sed -i "s/^>Reference\$/>${REF_LABEL}/" "$aln"
                done
                echo "$(date): ${taxon_key} - reference row relabelled as ${REF_LABEL} in core.aln/core.full.aln."
            else
                # Includes the "no SNPs detected" case (snp-sites has
                # nothing to output when every member is identical to the
                # reference at every core position) as well as genuine
                # snippy-core errors — either way, this taxon just doesn't
                # get a core alignment this run, not the whole batch.
                echo "$(date): WARNING - snippy-core failed for ${taxon_key} (see log above for the reason — may just be zero core SNPs); skipping this taxon." >&2
                FAILED_TAXA+=("$taxon_key")
            fi
        done

        conda deactivate || true

        # A single taxon's failure no longer fails the whole batch — it's
        # logged and recorded here so the rest of the batch can complete
        # (and so the next submission doesn't redo everything just because
        # one taxon hit an edge case). Check this file after a run.
        if [ "${#FAILED_TAXA[@]}" -gt 0 ]; then
            printf '%s\n' "${FAILED_TAXA[@]}" > "${STATUS_DIR}/snippy_failed_taxa.txt"
            echo "$(date): WARNING - Snippy failed for: ${FAILED_TAXA[*]} (see ${STATUS_DIR}/snippy_failed_taxa.txt)." >&2
        else
            rm -f "${STATUS_DIR}/snippy_failed_taxa.txt"
        fi

        touch "${STATUS_DIR}/snippy_complete.txt"
        echo "$(date): [STAGE 18/24] Snippy batch run complete."
    else
        echo "$(date): Not the last arrival (or lock already claimed) — deferring Snippy."
    fi
fi

echo "STAGE_19" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 18/24] Done for ${SAMPLE}."