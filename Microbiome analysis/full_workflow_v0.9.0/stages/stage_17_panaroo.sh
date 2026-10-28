#!/bin/bash
################################################################################
# STAGE 17: Panaroo - Targeted pangenome analysis (batch)
################################################################################
#
# Like STAGE_11 (MEDI), a pangenome has to be built from bins pooled across
# every sample in this submission, so this can only run once ALL array tasks
# have reached this stage. This script runs once per array task; only the
# task that arrives LAST (i.e. every other sample has already checked in)
# claims the run and launches Panaroo, once per target taxon group below.
#
# Bins are matched to a taxon by their name, using the labelling convention
# already present on disk (bakta_out/<SAMPLE>/<BIN_NAME>/), e.g.:
#   T1_bin1.orig_Enterobacter_flexneri_A
#   T1_bin5.strict_Streptococcus_lutetiensis
# i.e. <SAMPLE>_bin<N>.<version>_<Genus>[_<subgenus letter>]_<species...>
#
# INPUT FORMAT: panaroo is fed Bakta's .gbff (GenBank) files, collected into
# a manifest text file and passed as a single "-i manifest.txt" argument
# rather than listing paths directly on the command line.
#
# This specific combination was confirmed by hand on this cluster's data;
# nothing else tried worked: .gbff passed as multiple -i args failed with
# "Error reading prokka input!", and Bakta's .gff3 (even patched with a
# concatenated ##FASTA block, Prokka-style) failed the same way whether
# passed as multiple args or via a manifest file. Only .gbff read from a
# manifest file got through — apparently panaroo's multi-path CLI handling
# doesn't correctly detect GenBank format, but its manifest-list reader
# does. Each .gbff is sanity-checked for at least one CDS feature before use
# (Bakta has had versions where .gbff silently reports zero annotations).
#
# NOTE on Segatella/Prevotella naming:
# The GTDB-Tk database configured for this pipeline predates GTDB R214's
# split of Prevotella into Prevotella + Segatella, so bins that are
# taxonomically "Segatella copri" are still labelled "Prevotella_copri" in
# bakta_out/. The species-level pattern below matches on that actual label.
# A genus-wide "Segatella" pooling run is deliberately NOT included yet:
# Segatella is a subset of the old, much broader Prevotella genus, so a
# blanket "_Prevotella_" prefix match would sweep in many unrelated species.
# Once you can confirm the GTDB-Tk release in use (or just list the specific
# Prevotella-labelled species you want pooled as "Segatella"), add an entry
# here with an explicit species-name OR pattern, e.g.:
#   ["Segatella_genus"]="_Prevotella_(copri|stercorea|hominis)(_|$)"
#
# ASSUMPTIONS (easy to change — see the blocks flagged below):
#   1. Only bins that also passed STAGE_14 filtering (i.e. have a matching
#      fasta in good_bins/<SAMPLE>/) are pulled in, so low-quality/contaminated
#      bins don't get folded into the pangenome. To use ALL Bakta-annotated
#      bins regardless of STAGE_14 status instead, set REQUIRE_GOOD_BIN=false.
#   2. For species-level taxa, a bin is only included if it also passes a
#      fastANI check against reference/type-strain genomes you provide (see
#      "ACTION NEEDED" below). Genus-wide taxa skip this (REQUIRE_ANI=false)
#      since a handful of type strains can't represent whole-genus diversity.
#      Veillonella_genus was included early on purely as a two-sample smoke
#      test and has since been removed now that this is running on real
#      data — its snippy-core run hit a genuine zero-SNP edge case (its one
#      non-reference genome happened to match the ad-hoc in-group reference
#      at every core position), which isn't a bug so much as a sign that an
#      ad-hoc reference and small genus-wide group weren't a good fit for
#      SNP calling in the first place.
#
# >>> ACTION NEEDED (fastANI reference genomes):
#     Drop type-strain/reference fasta(s) into:
#       ${DIRECTORY}/panaroo_references/Bifidobacterium_longum/
#           (include BOTH the B. longum type strain and the B. longum
#           subsp. infantis type strain, as requested)
#       ${DIRECTORY}/panaroo_references/Escherichia_coli/
#       ${DIRECTORY}/panaroo_references/Prevotella_copri/
#       ${DIRECTORY}/panaroo_references/Bacteroides_fragilis/
#       ${DIRECTORY}/panaroo_references/Segatella_hominis/
#     Until a taxon's folder has at least one .fna/.fa/.fasta in it, ANI
#     filtering is skipped for that taxon (a warning is logged) and bins are
#     included on name-match alone, same as before this change.
# fastANI runs in the "fastani" conda env (FASTANI_ENV below).
#
# Input variables (from SLURM template):
#   - SAMPLE        : Sample name (current array task)
#   - DIRECTORY     : Workflow working directory
#   - NUM_THREADS   : Number of CPU threads
#   - STATUS_DIR    : Path to status directory
#   - PANAROO_OUT   : Output directory for panaroo results
#
# Output:
#   - ${PANAROO_OUT}/<taxon_key>/            - panaroo results per taxon group
#   - ${STATUS_DIR}/panaroo_complete.txt     - shared batch completion flag
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_19
#
################################################################################

echo "$(date): [STAGE 17/24] Panaroo check-in for ${SAMPLE}..."

# -- Configuration --------------------------------------------------------
# key : whether a bin must have passed STAGE_14 into good_bins/ to be
# included. Defaults to true for any taxon not listed here (see the
# unconditional default in the matching code below). Set a taxon to false
# for a "relaxed" variant that also includes bins STAGE_14 rejected,
# falling back to Bakta's own .fna when there's no good_bins copy -- useful
# for comparing a validated tree against a looser one on the same taxon.
declare -A REQUIRE_GOOD_BIN=(
    ["Bifidobacterium_longum"]=true
    ["Escherichia_coli"]=true
    ["Prevotella_copri"]=true
    ["Segatella_hominis"]=true
    ["Bacteroides_fragilis"]=true
    ["Segatella_genus"]=true
    ["Segatella_genus_relaxed"]=false
)
CLEAN_MODE="moderate"      # panaroo --clean-mode: strict|moderate|sensitive
                            # moderate is generally recommended when input
                            # comes from MAGs rather than curated isolates

# Taxa listed here skip Snippy/SNP-based tree building entirely (STAGE_18
# no-ops for them) and instead go straight from Panaroo's own core-gene
# alignment to RAxML-NG in STAGE_18B. Snippy assumes all genomes map onto
# one reference's coordinate system, which only holds for closely related
# strains -- across a whole genus, divergence is high enough that
# reference-coordinate mapping breaks down. Panaroo's own core-gene
# alignment (proper orthologous-cluster multiple sequence alignment, not
# reference mapping) tolerates that divergence, which is the whole point of
# a pangenome tool in the first place.
declare -A USE_PANAROO_ALIGNMENT=(
    ["Segatella_genus"]=true
    ["Segatella_genus_relaxed"]=true
)
CORE_ALIGNER="mafft"       # panaroo --aligner: mafft|prank|clustal|none
CORE_THRESHOLD=0.90        # panaroo --core_threshold: fraction of genomes a
                            # gene must be in to count as "core" (panaroo's
                            # own default is 0.95). Lowered here for
                            # genus-wide use since strict 95%-core across a
                            # whole genus can leave very few or zero
                            # universally-shared genes once you cross
                            # species boundaries -- watch the resulting
                            # alignment's gene/site count and adjust further
                            # if it comes out too small (or empty) to trust.

REFERENCE_DIR="${DIRECTORY}/panaroo_references"
FASTANI_ENV="fastani"
SPECIES_ANI_THRESHOLD=95   # fastANI %, ~95 is the conventional species boundary

# key : bash extended-regex pattern matched against BIN_NAME
declare -A TARGET_TAXA=(
    ["Bifidobacterium_longum"]="_Bifidobacterium_(longum|infantis)(_|$)"   # pools B. longum and B. longum subsp. infantis into one group
    ["Escherichia_coli"]="_Escherichia_coli(_|$)"
    ["Prevotella_copri"]="_Prevotella_copri(_|$)"           # aka "Segatella copri" - see note above
    ["Segatella_hominis"]="_Prevotella_hominis(_|$)"         # aka "Prevotella hominis" in this GTDB-Tk db; confirmed via bakta_out grep, distinct from Prevotella_intestinihominis
    ["Bacteroides_fragilis"]="_Bacteroides_fragilis(_|$)"
    ["Segatella_genus"]="\$^"   # sentinel: unreachable via GTDB-Tk matching by design (this taxon always
                                 # resolves via MetaPhlAn -- see METAPHLAN_OVERRIDE_GENERA). Kept non-empty
                                 # so a bin never matches by accident if METAPHLAN_OVERRIDE_GENERA is ever
                                 # cleared -- falls back to matching nothing rather than everything.
    ["Segatella_genus_relaxed"]="\$^"   # same sentinel, same reasoning
)

# key : whether to apply the fastANI reference check (see assumption 2)
declare -A REQUIRE_ANI=(
    ["Bifidobacterium_longum"]=true
    ["Escherichia_coli"]=true
    ["Prevotella_copri"]=true
    ["Segatella_hominis"]=true
    ["Bacteroides_fragilis"]=true
    ["Segatella_genus"]=false   # genus-wide pooling: a handful of species type strains can't represent whole-genus diversity
    ["Segatella_genus_relaxed"]=false
)

# key : the true/current genus for each taxon, used only to check against
# METAPHLAN_OVERRIDE_GENERA below. Independent of what GTDB-Tk's bin-folder
# label actually says (that mismatch is exactly the problem for Segatella).
declare -A TARGET_TAXA_GENUS=(
    ["Bifidobacterium_longum"]="Bifidobacterium"
    ["Escherichia_coli"]="Escherichia"
    ["Prevotella_copri"]="Segatella"
    ["Segatella_hominis"]="Segatella"
    ["Bacteroides_fragilis"]="Bacteroides"
    ["Segatella_genus"]="Segatella"
    ["Segatella_genus_relaxed"]="Segatella"
)

# Genera for which bin identity is resolved via the MetaPhlAn crosscheck
# (crosscheck_out/<SAMPLE>_MAG_metaphlan_crosscheck.tsv, column 4) instead
# of GTDB-Tk's bin-folder name. Add a genus here whenever MetaPhlAn's calls
# are trusted more than GTDB-Tk's for it. Segatella is the current case:
# this pipeline's GTDB-Tk database still uses the pre-split Prevotella name,
# and — more importantly — sometimes can't resolve species at all (GTDB
# placeholder names like "Prevotella sp934191715") for bins MetaPhlAn calls
# confidently (e.g. "Segatella copri" at 100%). For genera in this list,
# TARGET_TAXA's own pattern above is not consulted at all — matching goes
# entirely through the crosscheck lookup, so bins GTDB-Tk couldn't resolve
# are still caught.
METAPHLAN_OVERRIDE_GENERA=("Segatella")

# key : exact-match regex against column 4 (metaphlan_top_hit_species) of
# the crosscheck TSV. Only consulted for taxa whose TARGET_TAXA_GENUS is in
# METAPHLAN_OVERRIDE_GENERA above; every taxon in an override genus needs an
# entry here or that taxon is skipped with a warning (see below).
declare -A TARGET_TAXA_METAPHLAN_PATTERN=(
    ["Prevotella_copri"]="^Segatella copri$"
    ["Segatella_hominis"]="^Segatella hominis$"
    ["Segatella_genus"]="^(Candidatus )?Segatella "   # any species under the genus; deliberately
                                                         # overlaps with the two species-level entries
                                                         # above -- a copri/hominis bin is meant to land
                                                         # in both its species group AND this genus-wide one
    ["Segatella_genus_relaxed"]="^(Candidatus )?Segatella "   # same pattern; the only difference from
                                                                # Segatella_genus is REQUIRE_GOOD_BIN=false
)

if [ -f "${STATUS_DIR}/panaroo_complete.txt" ]; then
    echo "$(date): Panaroo already complete for this batch — skipping."
else
    touch "${STATUS_DIR}/stage17_arrived.${SLURM_ARRAY_JOB_ID}.${SAMPLE}"
    TOTAL_SAMPLES=$(wc -l < "${DIRECTORY}/samples.txt")
    ARRIVED=$(find "${STATUS_DIR}" -maxdepth 1 -name "stage17_arrived.${SLURM_ARRAY_JOB_ID}.*" | wc -l)
    echo "$(date): ${ARRIVED}/${TOTAL_SAMPLES} samples have reached STAGE_17."

    if [ "$ARRIVED" -ge "$TOTAL_SAMPLES" ] && mkdir "${STATUS_DIR}/stage17_panaroo.lock.${SLURM_ARRAY_JOB_ID}" 2>/dev/null; then
        echo "$(date): ${SAMPLE} is the last arrival — launching Panaroo runs..."

        mapfile -t BATCH_SAMPLES < "${DIRECTORY}/samples.txt"

        # -- MetaPhlAn crosscheck lookup (for METAPHLAN_OVERRIDE_GENERA) ------
        # Built once regardless of whether any taxon actually needs it; cheap,
        # and keeps the per-bin matching loop below simple. Key is
        # "SAMPLE|BIN_LABEL" (BIN_LABEL = bare "SAMPLE_binN", matching the
        # crosscheck TSV's own "bin" column) -> metaphlan_top_hit_species.
        CROSSCHECK_DIR="${DIRECTORY}/crosscheck_out"
        declare -A METAPHLAN_CALL
        if [ -d "$CROSSCHECK_DIR" ]; then
            for cc_file in "${CROSSCHECK_DIR}"/*_MAG_metaphlan_crosscheck.tsv; do
                [ -f "$cc_file" ] || continue
                cc_sample=$(basename "$cc_file" | sed 's/_MAG_metaphlan_crosscheck\.tsv$//')
                while IFS=$'\t' read -r bin_col gtdbtk_col nreads_col mpa_col _rest; do
                    [ "$bin_col" = "bin" ] && continue
                    METAPHLAN_CALL["${cc_sample}|${bin_col}"]="$mpa_col"
                done < "$cc_file"
            done
            echo "$(date): Loaded ${#METAPHLAN_CALL[@]} MetaPhlAn crosscheck call(s) from ${CROSSCHECK_DIR}."
        else
            echo "$(date): NOTE - ${CROSSCHECK_DIR} not found; any taxon in a METAPHLAN_OVERRIDE_GENERA genus will match nothing this run."
        fi

        WORKDIR="${STATUS_DIR}/panaroo_tmp.${SLURM_ARRAY_JOB_ID}"
        mkdir -p "$WORKDIR"

        # -- Phase 1: collect + fastANI-filter bins per taxon (into list files) --
        conda activate "$FASTANI_ENV" 2>/dev/null || echo "$(date): NOTE - could not activate ${FASTANI_ENV}; ANI checks will be skipped wherever a reference set is required."

        for taxon_key in "${!TARGET_TAXA[@]}"; do
            pattern="${TARGET_TAXA[$taxon_key]}"

            genus="${TARGET_TAXA_GENUS[$taxon_key]:-}"
            USE_METAPHLAN=false
            for g in "${METAPHLAN_OVERRIDE_GENERA[@]}"; do
                [ "$g" = "$genus" ] && USE_METAPHLAN=true && break
            done
            MPA_PATTERN="${TARGET_TAXA_METAPHLAN_PATTERN[$taxon_key]:-}"
            if [ "$USE_METAPHLAN" = true ] && [ -z "$MPA_PATTERN" ]; then
                echo "$(date): WARNING - ${taxon_key}'s genus (${genus}) is in METAPHLAN_OVERRIDE_GENERA but has no TARGET_TAXA_METAPHLAN_PATTERN entry; skipping this taxon entirely this run." >&2
                continue
            fi

            if [ "$USE_METAPHLAN" = true ]; then
                echo "$(date): Collecting bins for ${taxon_key} (via MetaPhlAn crosscheck, genus=${genus})..."
            else
                echo "$(date): Collecting bins for ${taxon_key} (via GTDB-Tk bin-folder name)..."
            fi

            REFSET_DIR="${REFERENCE_DIR}/${taxon_key}"
            REF_LIST_FILE="${WORKDIR}/${taxon_key}.refs.txt"
            HAVE_REFS=false
            if [ "${REQUIRE_ANI[$taxon_key]:-false}" = true ] && [ -d "$REFSET_DIR" ]; then
                find "$REFSET_DIR" -maxdepth 1 -type f \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" \) > "$REF_LIST_FILE"
                [ -s "$REF_LIST_FILE" ] && HAVE_REFS=true
            fi
            if [ "${REQUIRE_ANI[$taxon_key]:-false}" = true ] && [ "$HAVE_REFS" = false ]; then
                echo "$(date): WARNING - no reference genomes found in ${REFSET_DIR}; skipping ANI filtering for ${taxon_key} (name-match only)."
            fi

            LIST_FILE="${WORKDIR}/${taxon_key}.manifest.txt"
            MEMBERS_FILE="${WORKDIR}/${taxon_key}.members.txt"
            > "$LIST_FILE"
            > "$MEMBERS_FILE"

            for BSAMPLE in "${BATCH_SAMPLES[@]}"; do
                [ -z "$BSAMPLE" ] && continue
                SAMPLE_BAKTA_DIR="${DIRECTORY}/bakta_out/${BSAMPLE}"
                [ -d "$SAMPLE_BAKTA_DIR" ] || continue

                for bin_dir in "${SAMPLE_BAKTA_DIR}"/*/; do
                    [ -d "$bin_dir" ] || continue
                    BIN_NAME=$(basename "$bin_dir")

                    if [ "$USE_METAPHLAN" = true ]; then
                        BIN_LABEL="${BIN_NAME%%.*}"   # e.g. S003_bin22, matching crosscheck's "bin" column
                        MPA_SPECIES="${METAPHLAN_CALL["${BSAMPLE}|${BIN_LABEL}"]:-}"
                        if [ -z "$MPA_SPECIES" ]; then
                            continue   # no crosscheck row for this bin; not a match
                        fi
                        [[ "$MPA_SPECIES" =~ $MPA_PATTERN ]] || continue
                    else
                        [[ "$BIN_NAME" =~ $pattern ]] || continue
                    fi

                    GOOD_BIN_FASTA="${DIRECTORY}/good_bins/${BSAMPLE}/${BIN_NAME}.fa"
                    BAKTA_FNA="${bin_dir}${BIN_NAME}.fna"

                    if [ "${REQUIRE_GOOD_BIN[$taxon_key]:-true}" = true ]; then
                        [ -f "$GOOD_BIN_FASTA" ] || continue
                        MEMBER_FASTA="$GOOD_BIN_FASTA"
                    else
                        # Relaxed taxon: include bins that didn't clear
                        # STAGE_14's good_bins bar too, falling back to
                        # Bakta's own .fna (exists regardless of STAGE_14
                        # status) when there's no good_bins copy.
                        if [ -f "$GOOD_BIN_FASTA" ]; then
                            MEMBER_FASTA="$GOOD_BIN_FASTA"
                        elif [ -f "$BAKTA_FNA" ]; then
                            MEMBER_FASTA="$BAKTA_FNA"
                        else
                            continue   # no fasta available at all for this bin
                        fi
                    fi

                    gbff="${bin_dir}${BIN_NAME}.gbff"
                    [ -f "$gbff" ] || continue

                    # Skip silently-empty gbffs (known Bakta issue where the
                    # gbff reports zero annotations despite a populated gff3).
                    if ! grep -qE "^ +CDS" "$gbff"; then
                        echo "$(date): WARNING - ${BIN_NAME}.gbff has no CDS features, skipping."
                        continue
                    fi

                    # fastANI reference check (species-level taxa only)
                    if [ "$HAVE_REFS" = true ]; then
                        ANI_OUT="${WORKDIR}/${BIN_NAME}.ani.txt"
                        fastANI -q "$MEMBER_FASTA" --rl "$REF_LIST_FILE" -o "$ANI_OUT" >/dev/null 2>&1
                        BEST_ANI=$(cut -f3 "$ANI_OUT" 2>/dev/null | sort -rn | head -n1)
                        if [ -z "$BEST_ANI" ]; then
                            echo "$(date): ${BIN_NAME} - no fastANI hit against ${taxon_key} references (likely too divergent); excluding."
                            continue
                        fi
                        if ! awk -v a="$BEST_ANI" -v t="$SPECIES_ANI_THRESHOLD" 'BEGIN{exit !(a>=t)}'; then
                            echo "$(date): ${BIN_NAME} - best ANI ${BEST_ANI}% to ${taxon_key} references is below ${SPECIES_ANI_THRESHOLD}%; excluding."
                            continue
                        fi
                        echo "$(date): ${BIN_NAME} - best ANI ${BEST_ANI}% to ${taxon_key} references; included."
                    fi

                    echo "$gbff" >> "$LIST_FILE"
                    echo "$MEMBER_FASTA" >> "$MEMBERS_FILE"
                done
            done

            echo "$(date): ${taxon_key} - $(wc -l < "$LIST_FILE") matching bin(s) after filtering."
        done

        conda deactivate || true

        # -- Phase 2: run Panaroo per taxon on the filtered lists -----------------
        conda activate panaroo-env
        FAILED_TAXA=()

        for taxon_key in "${!TARGET_TAXA[@]}"; do
            LIST_FILE="${WORKDIR}/${taxon_key}.manifest.txt"
            N_GENOMES=$(wc -l < "$LIST_FILE")

            if [ "$N_GENOMES" -lt 2 ]; then
                echo "$(date): WARNING - fewer than 2 genomes for ${taxon_key}; skipping (panaroo needs >=2 genomes for a pangenome)."
                continue
            fi

            TAXON_OUT="${PANAROO_OUT}/${taxon_key}"
            mkdir -p "$TAXON_OUT"

            PANAROO_ARGS=(
                -i "$LIST_FILE"
                -o "$TAXON_OUT"
                --clean-mode "$CLEAN_MODE"
                --remove-invalid-genes
                -t "${NUM_THREADS}"
            )
            if [ "${USE_PANAROO_ALIGNMENT[$taxon_key]:-false}" = true ]; then
                PANAROO_ARGS+=(-a core --aligner "$CORE_ALIGNER" --core_threshold "$CORE_THRESHOLD")
            fi

            if panaroo "${PANAROO_ARGS[@]}"; then
                echo "$(date): Panaroo complete for ${taxon_key} (${N_GENOMES} genomes)."
                cp "${WORKDIR}/${taxon_key}.members.txt" "${TAXON_OUT}/group_members.txt"

                if [ "${USE_PANAROO_ALIGNMENT[$taxon_key]:-false}" = true ]; then
                    # Don't assume panaroo's exact output filename for -a core
                    # (it's varied across versions) -- find whatever
                    # alignment file it actually wrote and record its path
                    # for STAGE_18B, which skips Snippy for this taxon and
                    # reads straight from here instead.
                    CORE_ALN_FILE=$(find "$TAXON_OUT" -maxdepth 1 -iname "*core*align*.aln" | head -1)
                    if [ -z "$CORE_ALN_FILE" ]; then
                        CORE_ALN_FILE=$(find "$TAXON_OUT" -maxdepth 1 -iname "*.aln" | head -1)
                    fi
                    if [ -n "$CORE_ALN_FILE" ]; then
                        echo "${CORE_ALN_FILE}" > "${TAXON_OUT}/core_alignment_path.txt"
                        echo "$(date): ${taxon_key} - core alignment found at ${CORE_ALN_FILE}."
                    else
                        echo "$(date): WARNING - ${taxon_key} used USE_PANAROO_ALIGNMENT but no .aln file was found in ${TAXON_OUT}; STAGE_18B will have nothing to build a tree from." >&2
                        FAILED_TAXA+=("$taxon_key")
                    fi
                fi
            else
                echo "$(date): WARNING - Panaroo failed for ${taxon_key}; skipping this taxon, continuing with the rest." >&2
                FAILED_TAXA+=("$taxon_key")
            fi
        done

        conda deactivate || true
        rm -rf "$WORKDIR"

        # A single taxon's failure no longer fails the whole batch — it's
        # logged and recorded here so the rest of the batch can complete
        # (and so the next submission doesn't redo everything just because
        # one taxon hit an edge case). Check this file after a run.
        if [ "${#FAILED_TAXA[@]}" -gt 0 ]; then
            printf '%s\n' "${FAILED_TAXA[@]}" > "${STATUS_DIR}/panaroo_failed_taxa.txt"
            echo "$(date): WARNING - Panaroo failed for: ${FAILED_TAXA[*]} (see ${STATUS_DIR}/panaroo_failed_taxa.txt)." >&2
        else
            rm -f "${STATUS_DIR}/panaroo_failed_taxa.txt"
        fi

        touch "${STATUS_DIR}/panaroo_complete.txt"
        echo "$(date): [STAGE 17/24] Panaroo batch run complete."
    else
        echo "$(date): Not the last arrival (or lock already claimed) — deferring Panaroo."
    fi
fi

echo "STAGE_19" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 17/24] Done for ${SAMPLE}."