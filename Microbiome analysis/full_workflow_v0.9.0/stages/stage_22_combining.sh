#!/bin/bash

################################################################################
# STAGE 22: Combining - per-sample bin joins, cross-sample joins, barplots
################################################################################
#
# WHAT THIS STAGE DOES (4 levels of output)
#
#   LEVEL 1  Per-sample, per-bin (runs in EVERY array task, for its own sample)
#            - humann_join_tables over all of this sample's per-bin HUMAnN
#              tables, for EVERY table type found (or listed) -> one table per
#              sample per type, columns = bins
#            - humann_barplot on the FOCAL features only, NO metadata
#              (one bar per bin, raw per-bin values)
#            -> ${DIRECTORY}/combined_out/per_sample_bins/<SAMPLE>/
#
#   LEVEL 2  Whole-sample (STAGE_20) join across ALL samples + metadata
#            barplots (Part A). Also joins whole-sample counterparts of any
#            other bin table type, used for normalising Level 3.
#            -> ${DIRECTORY}/combined_out/whole_sample/
#
#   LEVEL 3  Bin-assigned aggregate (last-arrival only), ALL table types.
#            Sample-level table per type where each feature's species strata
#            are the species of the BIN the function came from. The per-bin
#            HUMAnN tables are depth-scaled (checked with bin_scale_diagnostic.py:
#            per-bin RPK tracks bin read count), so bin values are used at the
#            sample's own scale, NOT stretched to fill the whole-sample value:
#              - RPK-type tables (genefamilies, pathabundance): bin values are
#                added directly (they are already on the sample's scale)
#              - CPM / normalised tables: each bin value is multiplied by the
#                bin's share of the sample's RPK (bin total / whole-sample total
#                of genefamilies or pathabundance) to put it on sample scale
#              - coverage tables: the maximum over bins is taken (not additive)
#            What the bins do not account for is shown explicitly as the
#            "unclassified" stratum = max(0, whole-sample value - binned value),
#            so the binned fraction (typically well below 100%) is visible.
#            Where binned exceeds the whole-sample value (the two runs disagree)
#            the binned value is kept and there is no remainder.
#            Files: bin_aggregate_<type>_scaled.tsv (has remainder), _scaledbins
#            (no whole-sample table of that type: bins only), _maxbins (coverage).
#            Barplots (metadata grouped) are made for FOCAL features only.
#            -> ${DIRECTORY}/combined_out/bin_aggregate/
#
#   LEVEL 4  Per-species (STAGE_21 bins) cross-sample joins + barplots, for
#            B. infantis, S. copri, S. sinensis, E. coli and the whole
#            Segatella genus (Part B). Every species gets:
#              - the same focal pathways (FOCAL_FEATURES_GENERAL)
#              - its own species-specific features (e.g. B. infantis HMO ECs)
#              - "all pathways" variants (see LEVEL 5)
#            -> ${DIRECTORY}/combined_out/per_species/<species>/
#
#            SEGATELLA is assigned from MetaPhlAn, not GTDB-Tk: at the start of
#            the combine every bin whose MetaPhlAn top hit is a Segatella is
#            recorded in combined_out/segatella_bin_assignments.tsv. Those bins
#            take MetaPhlAn's species as their species for the bin-assigned
#            aggregate (Level 3) and define the segatella_copri /
#            segatella_sinensis / segatella_genus plot sets.
#
#   CUSTOM PATHWAYS  EC-based "pathways" (HMO, dietary-fibre depolymerisation
#            by substrate, resistant starch, pectin, and acetate / butyrate /
#            propionate production routes) are built by summing their member ECs
#            (rows AND species strata) from the ec-cpm-named tables, then
#            plotted like normal pathways: whole-sample + bin-assigned
#            (metadata barplots), per sample's bins (top-level groups only), and
#            per target species. Definitions: CUSTOM_PATHWAYS below.
#
#   LEVEL 5  "All pathways" barplots (one bar per sample / bin):
#              - whole-sample: all pathways summed, stacked by HUMAnN's own
#                species assignment, ordered/annotated by each metadata variant
#              - bin-assigned: same, but species = GTDB-Tk species of the
#                source bin (from the Level 3 aggregate)
#              - per species: stacked by HUMAnN taxon, and stacked by pathway
#            Each comes as _abs (absolute) and _rel (100% stacked) + a .tsv.
#            UNMAPPED / UNINTEGRATED / UNGROUPED are excluded (not pathways).
#
# ARCHITECTURE NOTE:
#   Same last-arrival pattern as STAGE_17. Level 1 runs in every task BEFORE
#   it checks in, so by the time the last task arrives every sample's per-bin
#   joined tables already exist. ${STATUS_DIR}/humann_combine_complete.txt is
#   the shared completion flag for Levels 2-4; Level 1 has its own per-sample
#   marker (${STATUS_DIR}/${SAMPLE}.stage22_bins_complete.txt). Delete the
#   relevant marker(s) to force a re-run.
#
# ============================================================================
# EDIT ME / CONFIRM BEFORE RUNNING:
#   1. Per-bin file naming is ASSUMED to be ${BIN_LABEL}_${SUFFIX}.tsv in
#      ${HUMANN_CONTIGS_DIR}. With BIN_TABLE_SUFFIXES=() all such suffixes are
#      auto-discovered; set BIN_TABLE_EXCLUDE_REGEX to drop non-HUMAnN files
#      or duplicate variants (e.g. "relab|stats").
#   2. Manifest columns assumed: 1=bin label, 2=GTDB-Tk species, 4=status.
#      If genus is held in a separate column, set MANIFEST_GENUS_COL.
#   3. Whole-sample counterparts are assumed to be
#      ${HUMANN_READS_DIR parent}/<SAMPLE>/*_${SUFFIX}.tsv (one per sample).
#   4. Verify focal feature ID strings against the joined tables.
# ============================================================================
#
# Input variables (from SLURM template):
#   - DIRECTORY          : project root
#   - SAMPLE             : this array task's sample
#   - HUMANN_READS_DIR   : per-sample STAGE_20 dir (parent derived below)
#   - HUMANN_CONTIGS_DIR : per-sample STAGE_21/21B dir (parent derived below)
#   - CROSSCHECK_OUT     : STAGE_12 crosscheck output (all-samples already)
#   - STATUS_DIR         : status directory
#
################################################################################

set -uo pipefail

################################################################################
# SHARED CONFIG
################################################################################

COMBINED_OUT="${DIRECTORY}/combined_out"
WHOLE_OUT="${COMBINED_OUT}/whole_sample"
SPECIES_OUT="${COMBINED_OUT}/per_species"
BINS_OUT="${COMBINED_OUT}/per_sample_bins"
AGG_OUT="${COMBINED_OUT}/bin_aggregate"
SEG_ASSIGN="${COMBINED_OUT}/segatella_bin_assignments.tsv"

# ---- Sample metadata (comma-delimited CSV, "sample_id" links to SAMPLE) ----
#   age_group : dietary-change based (milk_feeding / solid_feeding)
#   agegp     : flat 12-month intervals (0-12m / 13-24m)
METADATA_FILE="${DIRECTORY}/metadata_merged.csv"
METADATA_SAMPLE_COL="sample_id"
METADATA_VARIANTS=("age_group" "agegp")

# ---- Bin manifest columns (1-based) ----
MANIFEST_SPECIES_COL=2
MANIFEST_GENUS_COL=""     # EDIT ME if the manifest has a separate genus column

# ---- Segatella assignment (MetaPhlAn-based; see header) ----
# crosscheck_out/<SAMPLE>_MAG_metaphlan_crosscheck.tsv columns used (by name):
#   bin, gtdbtk_species, metaphlan_top_hit_species, metaphlan_top_hit_pct
SEGATELLA_METAPHLAN_REGEX="^segatella "   # case-insensitive, on metaphlan_top_hit_species
SEGATELLA_MIN_PCT=0                       # min metaphlan_top_hit_pct to accept (0 = any hit)
SEGATELLA_ALSO_GTDB=0                     # 1 = ALSO treat bins GTDB-Tk calls Segatella but
                                          #     MetaPhlAn doesn't as Segatella (keeps GTDB species)

# ---- Focal pathways (real MetaCyc IDs; verify against the first column of
# joined_pathabundance.tsv, before the ":" if present) ----
FOCAL_FEATURES_GENERAL=(
    "METSYN-PWY"
    "METHANOGENESIS-PWY"
    "PWY-6813"
    "PWY-6717"
    "PWY-6737"
)

# ---- Custom "pathways" (summed EC groups) ----
# For each group the member ECs' rows are summed (total AND per-species strata)
# from the ec-cpm-named tables and plotted like a normal pathway.
#
# EC NUMBERS WERE CHECKED AGAINST THEIR ENZYME NAMES (IUBMB/KEGG). Where an EC in
# the original lists is a different enzyme from the one named, the EC for the
# NAMED enzyme is used and the line says "listed as ...". ECs that are real
# enzymes with a relevant function stay in even if the label was off.
# Edit the lists freely; groups are unions built from the pieces below.

# ---- HMO ----
HMO_DIRECT_ECS="3.2.1.18 3.2.1.51 3.2.1.63 3.2.1.111 3.2.1.127 3.2.1.129 3.2.1.140 3.2.1.23 3.2.1.52 3.2.1.50"
HMO_DOWNSTREAM_ECS="3.1.1.53 4.2.2.15 2.4.1.211 3.2.1.85 3.2.1.38"

# ---- Fibre depolymerisation (by main substrate) ----
FIBRE_CELLULOSE_ECS="3.2.1.4 3.2.1.21 3.2.1.74"        # cellulase, beta-glucosidase, glucan 1,4-beta-glucosidase (listed as "endo-1,4-beta-xylanase")
FIBRE_XYLAN_ECS="3.2.1.8 3.2.1.37 3.2.1.55"            # 3.2.1.8 xylanase (listed as "mannanase"); 3.2.1.37 beta-xylosidase (listed as 3.2.1.39/3.2.1.78); 3.2.1.55 alpha-L-arabinofuranosidase (listed as 3.2.1.25/3.2.1.67)
FIBRE_MANNAN_ECS="3.2.1.78 3.2.1.25 3.2.1.22"          # 3.2.1.78 mannanase (listed as 3.2.1.8); 3.2.1.25 beta-mannosidase; 3.2.1.22 alpha-galactosidase
FIBRE_BETAGLUCAN_ECS="3.2.1.6 3.2.1.39"                # endo-1,3(4)-beta-glucanase; glucan endo-1,3-beta-glucosidase (3.2.1.6/3.2.1.39 listed as "alpha-L-galactosidase"/"beta-xylosidase")
FIBRE_FRUCTAN_ECS="3.2.1.26 3.2.1.80"                  # beta-fructofuranosidase; fructan beta-fructosidase (3.2.1.80 listed as "beta-glucanase")
FIBRE_GENERAL_ECS="${FIBRE_CELLULOSE_ECS} ${FIBRE_XYLAN_ECS} ${FIBRE_MANNAN_ECS} ${FIBRE_BETAGLUCAN_ECS} ${FIBRE_FRUCTAN_ECS} 3.2.1.51"   # 3.2.1.51 fucosidase (overlaps HMO-DIRECT)

# ---- Resistant starch ----
RESISTANT_STARCH_ECS="3.2.1.1 3.2.1.10 3.2.1.20 3.2.1.3"   # alpha-amylase, oligo-1,6-glucosidase, alpha-glucosidase, glucan 1,4-alpha-glucosidase

# ---- Pectin ----
PECTIN_BACKBONE_ECS="4.2.2.2 4.2.2.10 3.2.1.15 3.2.1.67 3.2.1.82"   # pectate lyase, pectin lyase, polygalacturonase; 3.2.1.67 galacturonan 1,4-alpha-galacturonidase and 3.2.1.82 exo-poly-alpha-galacturonosidase (listed as arabinofuranosidase / arabinanase)
PECTIN_ARABINAN_ECS="3.2.1.99 3.2.1.55"                              # 3.2.1.99 arabinan endo-1,5-alpha-L-arabinosidase (listed as 3.2.1.82); 3.2.1.55 arabinofuranosidase
PECTIN_RG_ECS="3.2.1.40 3.2.1.171 4.2.2.23 4.2.2.24"                 # alpha-L-rhamnosidase; rhamnogalacturonan HYDROLASE (3.2.1.171 listed as "lyase"); RG endo-/exo-lyase added
PECTIN_ECS="${PECTIN_BACKBONE_ECS} ${PECTIN_ARABINAN_ECS} ${PECTIN_RG_ECS}"

# ---- SCFA: acetate ----
SCFA_ACETATE_ECS="1.2.7.1 2.7.2.1 2.3.1.8"             # 1.2.7.1 pyruvate:ferredoxin oxidoreductase (listed as 1.2.7.4 = CO dehydrogenase); acetate kinase; phosphate acetyltransferase

# ---- SCFA: butyrate ----
BUTYRATE_UPPER_ECS="2.3.1.9 2.3.1.16 1.1.1.157 4.2.1.55 1.3.1.44 1.3.8.1"   # thiolases, 3-hydroxybutyryl-CoA DH, crotonase, trans-2-enoyl-CoA reductase; 1.3.8.1 butyryl-CoA dehydrogenase added
BUTYRATE_KINASE_ECS="2.3.1.19 2.7.2.7"                 # 2.3.1.19 phosphate butyryltransferase (listed as 2.7.2.7); 2.7.2.7 butyrate kinase (2.7.2.9 dropped: not butyrate kinase)
BUTYRATE_COAT_ECS="2.8.3.8"                            # butyryl-CoA:acetate CoA-transferase (key marker)

# ---- SCFA: propionate ----
PROP_SUCCINATE_ECS="6.4.1.1 1.1.1.37 4.2.1.2 1.3.1.6 1.3.5.4 5.4.99.2 2.1.3.1 7.2.4.3 4.1.1.41"   # 6.4.1.3 (propionyl-CoA carboxylase, listed as methylmalonyl-CoA carboxyltransferase) replaced by 2.1.3.1 transcarboxylase + 7.2.4.3 / 4.1.1.41 decarboxylase; 1.3.5.4 quinol fumarate reductase added
PROP_ACRYLATE_ECS="4.2.1.54 1.3.1.95 2.8.3.1"          # lactoyl-CoA dehydratase; 1.3.1.95 acryloyl-CoA reductase (listed as 1.3.1.31 = 2-enoate reductase); 2.8.3.1 propionate CoA-transferase (replaces 2.7.2.1)
PROP_PROPANEDIOL_ECS="1.2.1.87 1.1.1.2 4.2.1.28 2.3.1.222 2.7.2.15"   # 1.2.1.87 propanal dehydrogenase (listed as 1.2.1.27 = methylmalonate-semialdehyde DH); 1.1.1.2 ADH (non-specific); added: 4.2.1.28 diol dehydratase, 2.3.1.222 PduL, 2.7.2.15 propionate kinase
PROPIONATE_ECS="${PROP_SUCCINATE_ECS} ${PROP_ACRYLATE_ECS} ${PROP_PROPANEDIOL_ECS}"

# Each entry: ID::display name (no colons)::space-separated ECs
CUSTOM_PATHWAYS=(
    "HMO-DIRECT::Direct HMO degradation::${HMO_DIRECT_ECS}"
    "HMO-DOWNSTREAM::Downstream HMO utilisation::${HMO_DOWNSTREAM_ECS}"
    "FIBRE-ALL::Dietary fibre depolymerisation (fibre + resistant starch + pectin)::${FIBRE_GENERAL_ECS} ${RESISTANT_STARCH_ECS} ${PECTIN_ECS}"
    "FIBRE-GENERAL::Fibre-degrading enzymes (cellulose xylan mannan glucan fructan)::${FIBRE_GENERAL_ECS}"
    "FIBRE-CELLULOSE::Cellulose degradation::${FIBRE_CELLULOSE_ECS}"
    "FIBRE-XYLAN::Xylan and arabinoxylan degradation::${FIBRE_XYLAN_ECS}"
    "FIBRE-MANNAN::Mannan glucomannan galactomannan degradation::${FIBRE_MANNAN_ECS}"
    "FIBRE-BETAGLUCAN::Beta-glucan degradation::${FIBRE_BETAGLUCAN_ECS}"
    "FIBRE-FRUCTAN::Fructan and inulin degradation::${FIBRE_FRUCTAN_ECS}"
    "RESISTANT-STARCH::Resistant starch degradation::${RESISTANT_STARCH_ECS}"
    "PECTIN::Pectin degradation (all)::${PECTIN_ECS}"
    "PECTIN-ARABINAN::Arabinan degradation::${PECTIN_ARABINAN_ECS}"
    "PECTIN-RG::Rhamnogalacturonan degradation::${PECTIN_RG_ECS}"
    "SCFA-ACETATE::Acetate production::${SCFA_ACETATE_ECS}"
    "SCFA-BUTYRATE::Butyrate production (all routes)::${BUTYRATE_UPPER_ECS} ${BUTYRATE_KINASE_ECS} ${BUTYRATE_COAT_ECS}"
    "SCFA-BUTYRATE-KINASE::Butyrate - phosphotransbutyrylase and butyrate kinase route::${BUTYRATE_KINASE_ECS}"
    "SCFA-BUTYRATE-COAT::Butyrate - butyryl-CoA acetate CoA-transferase route::${BUTYRATE_COAT_ECS}"
    "SCFA-PROPIONATE::Propionate production (all routes)::${PROPIONATE_ECS}"
    "SCFA-PROP-SUCCINATE::Propionate - succinate pathway::${PROP_SUCCINATE_ECS}"
    "SCFA-PROP-ACRYLATE::Propionate - acrylate pathway::${PROP_ACRYLATE_ECS}"
    "SCFA-PROP-PROPANEDIOL::Propionate - propanediol pathway::${PROP_PROPANEDIOL_ECS}"
)
# IDs derived from the list above (kept in sync automatically)
CUSTOM_PATHWAY_IDS=()
for _cp in "${CUSTOM_PATHWAYS[@]}"; do
    CUSTOM_PATHWAY_IDS+=("${_cp%%::*}")
done
# Level 1 (one plot per sample per pathway) uses only the top-level groups, to
# keep the file count down; every pathway is plotted at the other levels.
CUSTOM_PATHWAY_IDS_LEVEL1=("HMO-DIRECT" "HMO-DOWNSTREAM" "FIBRE-ALL" "RESISTANT-STARCH" "PECTIN" "SCFA-ACETATE" "SCFA-BUTYRATE" "SCFA-PROPIONATE")

# Individual-EC focal list (all 15) used wherever single ECs are plotted
HMO_EC_FEATURES="${HMO_DIRECT_ECS} ${HMO_DOWNSTREAM_ECS}"

# ---- Per-bin table types ----
# Empty () = auto-discover every ${BIN_LABEL}_<suffix>.tsv in each sample's
# STAGE_21/21B directory. Or list explicitly, e.g. ("pathabundance" "ec-cpm-named").
# NOTE humann_join_tables' --file_name is a SUBSTRING match, which is why files
# are staged into a clean directory per type before joining.
BIN_TABLE_SUFFIXES=()
BIN_TABLE_EXCLUDE_REGEX=""

# Focal features per table type (space-separated). ONLY these types get
# barplots; every other type is joined and aggregated but not plotted.
declare -A BIN_TABLE_FOCAL=(
    ["pathabundance"]="${FOCAL_FEATURES_GENERAL[*]}"
    ["ec-cpm-named"]="${HMO_EC_FEATURES}"
)

# ---- "All pathways" plots ----
ALLPATH_TOP_N=15                                  # species / pathways shown individually
ALLPATH_EXCLUDE="UNMAPPED,UNINTEGRATED,UNGROUPED" # not real pathways

# ---- Table types made for each target species in Part B ----
SPECIES_TABLE_TYPES=("pathabundance" "ec-cpm-named")

################################################################################
# HELPERS
################################################################################

# make_custom_pathway_table <in.tsv> <out.tsv> [<ec_breakdown_out.tsv>]
# Builds the CUSTOM_PATHWAYS rows from an EC table (ec-cpm / ec-cpm-named style,
# optionally species-stratified): for each custom pathway, the total row is the
# sum of its member ECs' total rows, and each "|species" stratum row is the sum
# of the members' rows for that species. Returns non-zero if no member EC is
# present at all. If a third argument is given, also writes the individual
# member-EC rows (before summing) to that file as an EC composition table.
make_custom_pathway_table() {
    local in="$1" out="$2" breakdown_out="${3:-}"
    python3 - "${in}" "${out}" "${breakdown_out}" "${CUSTOM_PATHWAYS[@]}" << 'PYEOF'
import csv, sys
from collections import OrderedDict

in_path, out_path, breakdown_out = sys.argv[1:4]
defs = []
for d in sys.argv[4:]:
    pid, name, ecs = d.split("::", 2)
    defs.append((pid, name, set(ecs.split())))

with open(in_path) as fh:
    rd = csv.reader(fh, delimiter="\t")
    header = next(rd)
    ncol = len(header) - 1
    acc = OrderedDict()           # (pid, stratum) -> list of sums
    ec_rows = []                  # individual EC rows for the breakdown table
    members_seen = {pid: set() for pid, _, _ in defs}
    for r in rd:
        if not r:
            continue
        base, _, stratum = r[0].partition("|")
        ec = base.split(":", 1)[0].strip()
        for pid, _, ecs in defs:
            if ec in ecs:
                vals = []
                for x in r[1:ncol + 1]:
                    try:
                        vals.append(float(x))
                    except ValueError:
                        vals.append(0.0)
                vals += [0.0] * (ncol - len(vals))

                # Save for the detailed breakdown table
                if breakdown_out:
                    ec_rows.append([pid, r[0]] + vals)

                key = (pid, stratum)
                if key in acc:
                    acc[key] = [a + b for a, b in zip(acc[key], vals)]
                else:
                    acc[key] = vals
                if not stratum:
                    members_seen[pid].add(ec)

# Write the EC-level composition table if requested
if breakdown_out:
    with open(breakdown_out, "w", newline="") as bout:
        bw = csv.writer(bout, delimiter="\t", lineterminator="\n")
        bw.writerow(["# custom_pathway", "ec_feature"] + header[1:])
        for erow in ec_rows:
            bw.writerow([erow[0], erow[1]] + [repr(v) for v in erow[2:]])

n_written = 0
with open(out_path, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(header)
    for pid, name, ecs in defs:
        if (pid, "") not in acc:
            sys.stderr.write("NOTE: none of the {} ECs of {} found in {}\n".format(len(ecs), pid, in_path))
            continue
        label = "{}: {}".format(pid, name)
        w.writerow([label] + [repr(v) for v in acc[(pid, "")]])
        for (p2, stratum), vals in acc.items():
            if p2 == pid and stratum:
                w.writerow(["{}|{}".format(label, stratum)] + [repr(v) for v in vals])
        n_written += 1
        missing = sorted(ecs - members_seen[pid])
        sys.stderr.write("{}: {}/{} member ECs present{}\n".format(
            pid, len(ecs) - len(missing), len(ecs),
            "" if not missing else " (absent: {})".format(" ".join(missing))))
sys.exit(0 if n_written else 1)
PYEOF
}

# bin_unit_class <table type suffix>
# How a per-bin table's values combine into sample-level values:
#   additive            raw RPK-style (genefamilies, pathabundance): add directly
#   scale:genefamilies  CPM / normalised tables: multiply by the bin's share of
#                       the sample's gene-family RPK
#   scale:pathabundance normalised pathway tables: share of pathway abundance
#   coverage            not additive: take the maximum over bins
bin_unit_class() {
    case "$1" in
        *coverage*)                                                          echo "coverage" ;;
        pathabundance*cpm*|pathabundance*normalised*|pathabundance*relab*)   echo "scale:pathabundance" ;;
        *cpm*|*normalised*|*relab*)                                          echo "scale:genefamilies" ;;
        *)                                                                   echo "additive" ;;
    esac
}

# build_bin_scale_factors
# Writes ${BIN_SCALE_FACTORS}: for every (sample, bin) and each of genefamilies /
# pathabundance, the bin's table total, the whole-sample total and their ratio
# (the bin's share of the sample). Totals are sums of the unstratified rows
# (including UNMAPPED / UNGROUPED / UNINTEGRATED), i.e. what CPM was computed over.
build_bin_scale_factors() {
    python3 - "${DIRECTORY}/samples.txt" "${BINS_OUT}" "${CONTIGS_ROOT}" "${WHOLE_OUT}" "${BIN_SCALE_FACTORS}" << 'PYEOF'
import csv, os, re, sys
from collections import defaultdict
from statistics import median

samples_file, bins_root, contigs_root, whole_out, out_path = sys.argv[1:6]
with open(samples_file) as fh:
    samples = [l.strip() for l in fh if l.strip()]

def match_label(col, labels):
    best = None
    for lab in labels:
        if lab and re.search(r"(?<![A-Za-z0-9])" + re.escape(lab) + r"(?![A-Za-z0-9])", col):
            if best is None or len(lab) > len(best):
                best = lab
    return best

def column_totals(path, label_for_col):
    tot = defaultdict(float)
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        keys = [label_for_col(c) for c in header[1:]]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if "|" in f[0]:
                continue
            for k, x in zip(keys, f[1:]):
                if k is None:
                    continue
                try:
                    tot[k] += float(x)
                except ValueError:
                    pass
    return tot

kinds = ["genefamilies", "pathabundance"]
whole_tot = {}
for kind in kinds:
    pth = os.path.join(whole_out, "joined_{}.tsv".format(kind))
    if os.path.isfile(pth):
        sys.stderr.write("Scale factors: reading whole-sample {} totals...\n".format(kind))
        whole_tot[kind] = column_totals(pth, lambda c: match_label(c, samples))
    else:
        sys.stderr.write("WARNING: {} missing; no {} scale factors\n".format(pth, kind))

rows = []
for s in samples:
    manifest = os.path.join(contigs_root, s, "{}_bin_humann_manifest.tsv".format(s))
    if not os.path.isfile(manifest):
        continue
    bins = []
    with open(manifest) as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd, None)
        for r in rd:
            if len(r) > 3 and r[3] == "ok":
                bins.append(r[0])
    for kind in kinds:
        if kind not in whole_tot:
            continue
        j = os.path.join(bins_root, s, "joined_bins_{}.tsv".format(kind))
        if not os.path.isfile(j):
            continue
        wt = whole_tot[kind].get(s)
        if not wt:
            continue
        for b, t in column_totals(j, lambda c: match_label(c, bins)).items():
            rows.append([s, b, kind, t, wt, t / wt])

with open(out_path, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(["sample", "bin", "kind", "bin_total", "whole_total", "factor"])
    w.writerows(rows)

for kind in kinds:
    per = defaultdict(float)
    for r in rows:
        if r[2] == kind:
            per[r[0]] += r[5]
    if per:
        v = sorted(per.values())
        sys.stderr.write("Scale factors ({}): binned share of sample per sample: median {:.2f}, min {:.2f}, max {:.2f} ({} samples)\n".format(
            kind, median(v), v[0], v[-1], len(v)))
        if v[-1] > 1.0:
            sys.stderr.write("  WARNING: some samples have binned share > 1 -- bin tables are not on the sample's scale there.\n")
PYEOF
}

# build_segatella_assignments
# Run at the start of the combine. Writes ${SEG_ASSIGN} with one row per
# Segatella bin (manifest status=ok):
#   sample bin gtdbtk_species metaphlan_species metaphlan_pct source effective_species
# source = "metaphlan"    : MetaPhlAn top hit matches SEGATELLA_METAPHLAN_REGEX
#                           -> effective_species = MetaPhlAn's species (taken as true)
#        = "gtdbtk_only"  : GTDB-Tk says Segatella / Prevotella copri but MetaPhlAn
#                           doesn't -> effective_species = GTDB-Tk's species
build_segatella_assignments() {
    python3 - "${DIRECTORY}/samples.txt" "${CONTIGS_ROOT}" "${CROSSCHECK_OUT}" "${SEG_ASSIGN}" \
              "${SEGATELLA_METAPHLAN_REGEX}" "${SEGATELLA_MIN_PCT}" "${MANIFEST_SPECIES_COL}" << 'PYEOF'
import csv, os, re, sys
from collections import Counter

samples_file, contigs_root, crosscheck_out, out_path, mp_regex, min_pct, species_col = sys.argv[1:8]
mp_re = re.compile(mp_regex, re.I)
gt_re = re.compile(r"^(Segatella |Prevotella copri)", re.I)
min_pct = float(min_pct)
species_col = int(species_col) - 1

with open(samples_file) as fh:
    samples = [l.strip() for l in fh if l.strip()]

rows = []
n_no_cc = 0
n_cc_not_in_manifest = 0
n_manifest_bins = 0
n_matched_labels = 0

for sample in samples:
    manifest = os.path.join(contigs_root, sample, "{}_bin_humann_manifest.tsv".format(sample))
    cc_path = os.path.join(crosscheck_out, "{}_MAG_metaphlan_crosscheck.tsv".format(sample))
    if not os.path.isfile(manifest):
        continue

    bins = {}                                   # bin -> gtdb species (status ok only)
    with open(manifest) as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd, None)
        for r in rd:
            if len(r) > 3 and r[3] == "ok":
                bins[r[0]] = r[species_col] if len(r) > species_col else ""
    n_manifest_bins += len(bins)

    cc = {}                                     # bin -> (metaphlan species, pct)
    if os.path.isfile(cc_path):
        with open(cc_path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                try:
                    pct = float(r.get("metaphlan_top_hit_pct", "") or 0)
                except ValueError:
                    pct = 0.0
                cc[r.get("bin", "")] = ((r.get("metaphlan_top_hit_species", "") or "").strip(), pct)
    else:
        n_no_cc += 1

    n_cc_not_in_manifest += sum(1 for b in cc if b not in bins)
    n_matched_labels += sum(1 for b in bins if b in cc)

    for b, gt in bins.items():
        mp, pct = cc.get(b, ("", 0.0))
        if mp and mp_re.search(mp) and pct >= min_pct:
            rows.append([sample, b, gt, mp, "{:.2f}".format(pct), "metaphlan", mp])
        elif gt_re.search(gt or ""):
            rows.append([sample, b, gt, mp, "{:.2f}".format(pct), "gtdbtk_only", gt])

with open(out_path, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(["sample", "bin", "gtdbtk_species", "metaphlan_species", "metaphlan_pct",
                "source", "effective_species"])
    w.writerows(rows)

mpr = [r for r in rows if r[5] == "metaphlan"]
gto = [r for r in rows if r[5] == "gtdbtk_only"]
changed = [r for r in mpr if r[2] != r[6]]
sys.stderr.write("Segatella check: {} bin(s) via MetaPhlAn, {} GTDB-Tk-only bin(s) "
                 "(of {} manifest bins with status ok).\n".format(len(mpr), len(gto), n_manifest_bins))
sys.stderr.write("  MetaPhlAn species taken for those bins: {}\n".format(
    dict(Counter(r[6] for r in mpr).most_common(10))))
sys.stderr.write("  {} of the MetaPhlAn Segatella bins have a different GTDB-Tk species label "
                 "(the MetaPhlAn one is used).\n".format(len(changed)))
if gto:
    sys.stderr.write("  GTDB-Tk-only bins (excluded unless SEGATELLA_ALSO_GTDB=1): {}\n".format(
        [(r[0], r[1], r[2]) for r in gto[:5]]))
if n_no_cc:
    sys.stderr.write("WARNING: {} sample(s) have no crosscheck file.\n".format(n_no_cc))
if n_manifest_bins and n_matched_labels == 0:
    sys.stderr.write("WARNING: no crosscheck bin label matched any manifest bin label -- "
                     "MetaPhlAn assignment cannot work. Compare the 'bin' column of the "
                     "crosscheck file with column 1 of the manifest.\n")
elif n_cc_not_in_manifest:
    sys.stderr.write("NOTE: {} crosscheck bin(s) are not status=ok manifest bins (ignored).\n".format(
        n_cc_not_in_manifest))
PYEOF
}

# plot_all_pathways --abundance <tsv> --out-prefix <prefix> --stack-by species|pathway
#                   [--metadata <csv> --sample-col <col> --group-col <col>]
# Stacked barplot over ALL pathways (one bar per table column). Writes
# <prefix>_abs.png, <prefix>_rel.png, <prefix>.tsv. Uses numpy + matplotlib
# (already required by humann_barplot).
plot_all_pathways() {
    python3 - --top-n "${ALLPATH_TOP_N}" --exclude "${ALLPATH_EXCLUDE}" "$@" << 'PYEOF'
#!/usr/bin/env python3
"""
Stacked barplots summarising ALL pathways in a HUMAnN-style table (one bar per
table column, i.e. per sample or per bin).

  --stack-by species : each bar is the total abundance over all pathways,
        split by the species contributing it. Uses the stratified rows
        ("PWY-X: name|g__Genus.s__Genus_species"). Run on the whole-sample
        table this gives HUMAnN's own species assignment; run on the
        bin-aggregate table it gives species assigned by source bin.

  --stack-by pathway : each bar is split by pathway (top N + Other), using the
        unstratified (community-level) rows. Shows pathway composition.

Writes <prefix>_abs.png (absolute), <prefix>_rel.png (100% stacked) and
<prefix>.tsv (the plotted matrix). Optional metadata orders/annotates columns
by sample group. Needs numpy + matplotlib only.
"""
import argparse
import csv
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

OTHER = "Other"
UNCLASSIFIED = "unclassified"


def match_label(col, labels):
    """Longest label found in the header as a whole token (not flanked by
    letters/digits): 'out.S115_combined_Abundance' -> 'S115'; S11 != S115."""
    best = None
    for lab in labels:
        if lab and re.search(r"(?<![A-Za-z0-9])" + re.escape(lab) + r"(?![A-Za-z0-9])", col):
            if best is None or len(lab) > len(best):
                best = lab
    return best


def species_label(stratum):
    if stratum.lower() == UNCLASSIFIED:
        return UNCLASSIFIED
    m = re.search(r"s__(.+)$", stratum)
    if m:
        return m.group(1).replace("_", " ")
    m = re.search(r"g__(.+?)(\.|$)", stratum)
    if m:
        return m.group(1).replace("_", " ") + " (genus)"
    return stratum


def short_col(col):
    c = re.sub(r"^out\.", "", col)
    c = re.sub(r"_Abundance(-[A-Za-z0-9]+)?$", "", c)
    c = re.sub(r"_combined$", "", c)
    return c


def to_float(x):
    try:
        return float(x)
    except ValueError:
        return 0.0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--abundance", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--stack-by", choices=["species", "pathway"], default="species")
    ap.add_argument("--top-n", type=int, default=15)
    ap.add_argument("--exclude", default="UNMAPPED,UNINTEGRATED,UNGROUPED",
                    help="Comma-separated feature IDs to leave out (not real pathways)")
    ap.add_argument("--metadata", default="")
    ap.add_argument("--sample-col", default="sample_id")
    ap.add_argument("--group-col", default="")
    ap.add_argument("--group-order", default="")
    ap.add_argument("--title", default="")
    args = ap.parse_args()

    exclude = set(x.strip() for x in args.exclude.split(",") if x.strip())

    with open(args.abundance) as fh:
        rd = csv.reader(fh, delimiter="\t")
        header = next(rd)
        cols = header[1:]
        ncol = len(cols)
        mat = {}
        for r in rd:
            if not r:
                continue
            base, _, stratum = r[0].partition("|")
            fid = base.split(":", 1)[0].strip()
            if fid in exclude:
                continue
            if args.stack_by == "species":
                if not stratum:
                    continue
                key = species_label(stratum)
            else:
                if stratum:
                    continue
                key = fid if ":" not in base else "{}: {}".format(fid, base.split(":", 1)[1].strip())
                if len(key) > 48:
                    key = key[:45] + "..."
            vals = np.zeros(ncol)
            for i, x in enumerate(r[1:ncol + 1]):
                vals[i] = to_float(x)
            if key in mat:
                mat[key] += vals
            else:
                mat[key] = vals

    if not mat:
        sys.exit("ERROR: no usable rows in {} for --stack-by {} (stratified rows needed for 'species')".format(
            args.abundance, args.stack_by))

    # ---- collapse to top N + Other (+ unclassified kept separate) ----------
    named = [k for k in mat if k != UNCLASSIFIED]
    named.sort(key=lambda k: -mat[k].sum())
    keep, rest = named[:args.top_n], named[args.top_n:]
    final = [(k, mat[k]) for k in keep]
    if rest:
        final.append((OTHER, sum(mat[k] for k in rest)))
    if UNCLASSIFIED in mat:
        final.append((UNCLASSIFIED, mat[UNCLASSIFIED]))
    labels = [k for k, _ in final]
    data = np.vstack([v for _, v in final])
    totals = data.sum(axis=0)

    # ---- column order / groups ----------------------------------------------
    col_group = None
    if args.metadata and args.group_col:
        with open(args.metadata) as fh:
            meta = list(csv.DictReader(fh, delimiter=","))
        sid2g = {}
        for m in meta:
            g = (m.get(args.group_col) or "").strip()
            sid2g[m[args.sample_col]] = g if g and g.upper() != "NA" else "NA"
        col_group = []
        for c in cols:
            s = match_label(c, list(sid2g))
            col_group.append(sid2g.get(s, "NA") if s else "NA")
        uniq = sorted(set(col_group) - {"NA"})
        if args.group_order:
            uniq = [g.strip() for g in args.group_order.split(",") if g.strip() and g.strip() in uniq]
        gorder = uniq + (["NA"] if "NA" in col_group else [])
        rank = {g: i for i, g in enumerate(gorder)}
        order = sorted(range(ncol), key=lambda i: (rank.get(col_group[i], 99), -totals[i]))
    else:
        order = sorted(range(ncol), key=lambda i: -totals[i])

    data = data[:, order]
    names = [short_col(cols[i]) for i in order]
    groups_sorted = [col_group[i] for i in order] if col_group else None

    # ---- colours ---------------------------------------------------------------
    pal = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors)
    colours = {}
    j = 0
    for k in labels:
        if k == OTHER:
            colours[k] = (0.55, 0.55, 0.55)
        elif k == UNCLASSIFIED:
            colours[k] = (0.85, 0.85, 0.85)
        else:
            colours[k] = pal[j % len(pal)]
            j += 1

    # ---- matrix out ------------------------------------------------------------
    with open(args.out_prefix + ".tsv", "w") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["# " + ("species" if args.stack_by == "species" else "pathway")] + names)
        for k, row in zip(labels, data):
            w.writerow([k] + [repr(float(x)) for x in row])

    gcol = {}
    if groups_sorted:
        gp = plt.get_cmap("Set2").colors
        for i, g in enumerate(dict.fromkeys(groups_sorted)):
            gcol[g] = (0.8, 0.8, 0.8) if g == "NA" else gp[i % len(gp)]

    n = data.shape[1]
    x = np.arange(n)
    tick_fs = max(3.5, min(8, 520.0 / max(n, 1)))
    stack_name = "species" if args.stack_by == "species" else "pathway"

    for mode in ("abs", "rel"):
        plot = data.copy()
        if mode == "rel":
            tot = plot.sum(axis=0)
            tot[tot == 0] = 1.0
            plot = plot / tot * 100.0

        width = max(9.0, 0.17 * n + 5.5)
        if groups_sorted:
            fig, (axg, ax) = plt.subplots(2, 1, figsize=(width, 7.2), sharex=True,
                                          gridspec_kw={"height_ratios": [0.35, 10], "hspace": 0.02})
            for xi, g in zip(x, groups_sorted):
                axg.bar(xi, 1, width=1.0, color=gcol[g], linewidth=0)
            axg.set_ylim(0, 1)
            axg.axis("off")
            axg.legend([Patch(facecolor=gcol[g]) for g in gcol], list(gcol), title=args.group_col,
                       loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=8)
        else:
            fig, ax = plt.subplots(figsize=(width, 6.8))

        bottom = np.zeros(n)
        for k, row in zip(labels, plot):
            ax.bar(x, row, bottom=bottom, width=0.9, color=colours[k], linewidth=0)
            bottom += row

        ax.set_xlim(-0.7, n - 0.3)
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=90, fontsize=tick_fs)
        ax.set_ylabel("Relative contribution (%)" if mode == "rel"
                      else "Total abundance over all pathways (table units)")
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.legend([Patch(facecolor=colours[k]) for k in labels], labels,
                  title="Species" if stack_name == "species" else "Pathway",
                  loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=8)
        ttl = args.title or "All pathways, stacked by {}".format(stack_name)
        (axg if groups_sorted else ax).set_title(ttl, fontsize=11)
        out_png = "{}_{}.png".format(args.out_prefix, mode)
        fig.savefig(out_png, dpi=200, bbox_inches="tight")
        plt.close(fig)
        sys.stderr.write("Wrote {}\n".format(out_png))


if __name__ == "__main__":
    main()
PYEOF
}

# allpath_species_plots <table.tsv> <out_dir> <tag>
# Per-sample "all pathways" plots stacked by species, once per metadata variant
# (columns ordered/annotated by group), or once without metadata if none.
allpath_species_plots() {
    local table="$1" out_dir="$2" tag="$3" variant
    if [ ! -s "${table}" ]; then
        echo "$(date): [STAGE 22/24] ERROR: ${table} missing -- skipping all-pathways plots (${tag})." >&2
        return 1
    fi
    if [ -s "${METADATA_FILE}" ]; then
        for variant in "${METADATA_VARIANTS[@]}"; do
            echo "$(date): [STAGE 22/24] All-pathways plot (${tag}, grouped by ${variant})..."
            plot_all_pathways --abundance "${table}" \
                --out-prefix "${out_dir}/allpathways_${tag}_${variant}" \
                --stack-by species \
                --metadata "${METADATA_FILE}" --sample-col "${METADATA_SAMPLE_COL}" --group-col "${variant}" \
                --title "All pathways by species (${tag}), grouped by ${variant}" \
                || echo "$(date): [STAGE 22/24]   WARNING: all-pathways plot failed (${tag}, ${variant})"
        done
    else
        echo "$(date): [STAGE 22/24] All-pathways plot (${tag}, no metadata)..."
        plot_all_pathways --abundance "${table}" \
            --out-prefix "${out_dir}/allpathways_${tag}" \
            --stack-by species \
            --title "All pathways by species (${tag})" \
            || echo "$(date): [STAGE 22/24]   WARNING: all-pathways plot failed (${tag})"
    fi
}

# feature_present <table.tsv> <feature_id>
# True if first column (ignoring any "|taxon" stratum and ":" description)
# equals the feature ID. Avoids humann_barplot hard-failing on undetected features.
feature_present() {
    awk -F'\t' -v f="$2" '
        { split($1, a, "|"); split(a[1], b, ":"); if (b[1] == f) found = 1 }
        END { exit !found }' "$1"
}

# discover_bin_suffixes <sample_dir> <bin> [<bin> ...]
# Prints the distinct <suffix> values of ${BIN}_<suffix>.tsv files. Files that
# belong to a longer bin label (bin_1 vs bin_1_2) are not misattributed.
discover_bin_suffixes() {
    local sample_dir="$1"
    shift
    local bins=("$@") bin other f name suffix skip
    for bin in "${bins[@]}"; do
        for f in "${sample_dir}/${bin}_"*.tsv; do
            [ -e "${f}" ] || continue
            name=$(basename "${f}" .tsv)
            skip=0
            for other in "${bins[@]}"; do
                if [ "${#other}" -gt "${#bin}" ] && [[ "${name}" == "${other}_"* ]]; then
                    skip=1
                    break
                fi
            done
            [ "${skip}" -eq 1 ] && continue
            suffix="${name#"${bin}_"}"
            if [ -n "${BIN_TABLE_EXCLUDE_REGEX}" ] && [[ "${suffix}" =~ ${BIN_TABLE_EXCLUDE_REGEX} ]]; then
                continue
            fi
            echo "${suffix}"
        done
    done | sort -u
}

# build_pcl <abundance.tsv> <pcl_out> <metadata_col> [<metadata_col> ...]
# Prepends metadata rows to a HUMAnN table to make a PCL for humann_barplot.
build_pcl() {
    local abund="$1" pcl="$2"
    shift 2
    python3 - "${METADATA_FILE}" "${METADATA_SAMPLE_COL}" "${abund}" "${pcl}" "$@" << 'PYEOF'
import csv, sys

meta_path, sample_col, abund_path, out_path = sys.argv[1:5]
value_cols = sys.argv[5:]

# Metadata is a comma-delimited CSV, unlike every HUMAnN3 table (tab-delimited).
with open(meta_path) as fh:
    reader = csv.DictReader(fh, delimiter=",")
    meta_rows = list(reader)

# Duplicate sample IDs would silently overwrite (last wins) -- check first.
ids = [r[sample_col] for r in meta_rows]
dupes = sorted(set(i for i in ids if ids.count(i) > 1))
if dupes:
    sys.stderr.write(
        "WARNING: {} duplicate {} value(s) in {} (e.g. {}) -- only the LAST row for "
        "each is used. Verify this file is really one row per sample before trusting "
        "the metadata in the PCL.\n".format(len(dupes), sample_col, meta_path, dupes[:5])
    )

meta = {row[sample_col]: row for row in meta_rows}

with open(abund_path) as fh:
    reader = csv.reader(fh, delimiter="\t")
    header = next(reader)
    rows = list(reader)

sample_cols = header[1:]

# humann_join_tables suffixes headers (e.g. "S001_Abundance"), so fall back to
# prefix matching. Ambiguous prefixes are left unmatched rather than guessed.
import re

def match_label(col, labels):
    """Longest sample ID found in the header as a whole token, i.e. not
    flanked by letters/digits. Handles 'S115_Abundance' and
    'out.S115_combined_Abundance-RPKs' alike, and S11 never matches S115."""
    best = None
    for lab in labels:
        if lab and re.search(r"(?<![A-Za-z0-9])" + re.escape(lab) + r"(?![A-Za-z0-9])", col):
            if best is None or len(lab) > len(best):
                best = lab
    return best

col_to_meta_key = {}
for col in sample_cols:
    lab = match_label(col, list(meta))
    if lab is not None:
        col_to_meta_key[col] = lab

unmatched = [c for c in sample_cols if c not in col_to_meta_key]
if unmatched:
    sys.stderr.write(
        "WARNING: {} of {} sample columns have no metadata match even after prefix "
        "matching (e.g. {}). Check METADATA_SAMPLE_COL / sample naming.\n".format(
            len(unmatched), len(sample_cols), unmatched[:5]
        )
    )

def lookup(col, value_col):
    row = meta.get(col_to_meta_key.get(col, ""), {})
    val = row.get(value_col, "")
    return val if val not in ("", None) else "NA"

with open(out_path, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(header)
    for col in value_cols:
        writer.writerow([col] + [lookup(c, col) for c in sample_cols])
    for row in rows:
        writer.writerow(row)

sys.stderr.write("Wrote {} with {} metadata row(s) + {} feature row(s)\n".format(
    out_path, len(value_cols), len(rows)
))
PYEOF
}

# barplots_by_metadata <abundance.tsv> <out_dir> <tag> <feature> [<feature> ...]
# For each metadata variant: build a PCL, then two barplots per feature.
barplots_by_metadata() {
    local abund="$1" out_dir="$2" tag="$3"
    shift 3
    local features=("$@")
    local variant pcl feature

    if [ ! -s "${abund}" ]; then
        echo "$(date): [STAGE 22/24] ERROR: table ${abund} is missing or empty -- skipping barplots for ${tag}." >&2
        return 1
    fi

    if [ ! -s "${METADATA_FILE}" ]; then
        echo "$(date): [STAGE 22/24] METADATA_FILE not found (${METADATA_FILE}) -- skipping metadata barplots for ${tag}."
        return 0
    fi

    for variant in "${METADATA_VARIANTS[@]}"; do
        pcl="${out_dir}/hmp_${tag}_${variant}.pcl"
        echo "$(date): [STAGE 22/24] Building PCL for '${tag}' x '${variant}'..."
        build_pcl "${abund}" "${pcl}" "${variant}"

        for feature in "${features[@]}"; do
            if ! feature_present "${abund}" "${feature}"; then
                echo "$(date): [STAGE 22/24]   '${feature}' not in ${abund} -- skipping barplots."
                continue
            fi

            echo "$(date): [STAGE 22/24] humann_barplot '${tag}' x '${variant}' x '${feature}' (sort=sum)..."
            humann_barplot \
                --input "${pcl}" \
                --focal-metadata "${variant}" \
                --last-metadata "${variant}" \
                --output "${out_dir}/barplot_${tag}_${feature}_${variant}_sort-sum.png" \
                --focal-feature "${feature}" \
                --sort sum

            echo "$(date): [STAGE 22/24] humann_barplot '${tag}' x '${variant}' x '${feature}' (sort=sum+metadata, logstack)..."
            humann_barplot \
                --input "${pcl}" \
                --focal-metadata "${variant}" \
                --last-metadata "${variant}" \
                --output "${out_dir}/barplot_${tag}_${feature}_${variant}_sort-summeta-logstack.png" \
                --focal-feature "${feature}" \
                --sort sum metadata \
                --scaling logstack
        done
    done
}

################################################################################
# LEVEL 1: per-sample, per-bin join + focal barplots (no metadata). Every task.
################################################################################

run_stage22_per_sample_bins() {
    local marker="${STATUS_DIR}/${SAMPLE}.stage22_bins_complete.txt"
    local sample_dir="${HUMANN_CONTIGS_DIR}"
    local manifest="${sample_dir}/${SAMPLE}_bin_humann_manifest.tsv"

    # Trust the marker only if its outputs still exist (or the sample has no
    # manifest, in which case there is nothing to join). A marker whose outputs
    # were deleted is stale and is ignored.
    if [ -f "${marker}" ]; then
        if [ ! -s "${manifest}" ] || compgen -G "${BINS_OUT}/${SAMPLE}/joined_bins_*.tsv" > /dev/null; then
            echo "$(date): [STAGE 22/24] Per-bin combine already done for ${SAMPLE} -- skipping."
            return 0
        fi
        echo "$(date): [STAGE 22/24] Marker for ${SAMPLE} exists but its outputs are missing -- redoing per-bin combine."
        rm -f "${marker}"
    fi
    if [ ! -s "${manifest}" ]; then
        echo "$(date): [STAGE 22/24] No bin manifest for ${SAMPLE} (${manifest}) -- skipping per-bin combine."
        touch "${marker}"
        return 0
    fi

    echo "$(date): [STAGE 22/24] Per-bin combine for ${SAMPLE}..."
    conda activate humann3-env

    local out="${BINS_OUT}/${SAMPLE}"
    mkdir -p "${out}"

    # Bins that completed HUMAnN (status=ok), same filter as Part B
    local bins
    bins=$(awk -F'\t' 'NR>1 && $4=="ok" {print $1}' "${manifest}")

    local suffixes=()
    if [ "${#BIN_TABLE_SUFFIXES[@]}" -gt 0 ]; then
        suffixes=("${BIN_TABLE_SUFFIXES[@]}")
    else
        # shellcheck disable=SC2086
        mapfile -t suffixes < <(discover_bin_suffixes "${sample_dir}" ${bins})
    fi
    echo "$(date): [STAGE 22/24]   ${#suffixes[@]} bin table type(s) for ${SAMPLE}: ${suffixes[*]:-none}"

    local suffix stage n bin src joined feature
    local feats=()
    for suffix in "${suffixes[@]}"; do
        stage="${out}/staging_${suffix}"
        rm -rf "${stage}"
        mkdir -p "${stage}"

        n=0
        for bin in ${bins}; do
            src="${sample_dir}/${bin}_${suffix}.tsv"
            if [ -s "${src}" ]; then
                cp "${src}" "${stage}/${bin}_${suffix}.tsv"
                n=$((n + 1))
            fi
        done

        if [ "${n}" -eq 0 ]; then
            echo "$(date): [STAGE 22/24]   No ${suffix} bin tables for ${SAMPLE} -- skipping."
            continue
        fi

        joined="${out}/joined_bins_${suffix}.tsv"
        if ! humann_join_tables -i "${stage}" -o "${joined}" --file_name "${suffix}"; then
            echo "$(date): [STAGE 22/24]   WARNING: humann_join_tables failed for ${SAMPLE} ${suffix} (not a HUMAnN table?) -- skipping."
            rm -f "${joined}"
            continue
        fi
        echo "$(date): [STAGE 22/24]   Joined ${n} bin(s) (${suffix}) for ${SAMPLE} -> ${joined}"

        # Custom pathways (summed ECs) per bin, from the EC table
        if [ "${suffix}" = "ec-cpm-named" ]; then
            if make_custom_pathway_table "${joined}" "${out}/joined_bins_custom-pathways.tsv" "${out}/joined_bins_custom-pathways_ec_breakdown.tsv"; then
                for feature in "${CUSTOM_PATHWAY_IDS_LEVEL1[@]}"; do
                    feature_present "${out}/joined_bins_custom-pathways.tsv" "${feature}" || continue
                    humann_barplot \
                        --input "${out}/joined_bins_custom-pathways.tsv" \
                        --output "${out}/barplot_bins_custom-pathways_${feature}.png" \
                        --focal-feature "${feature}" \
                        --sort sum
                done
            fi
        fi

        # Barplots only where focal features are defined for this table type
        feats=()
        read -ra feats <<< "${BIN_TABLE_FOCAL[$suffix]:-}"
        if [ "${#feats[@]}" -eq 0 ]; then
            continue
        fi
        for feature in "${feats[@]}"; do
            if ! feature_present "${joined}" "${feature}"; then
                echo "$(date): [STAGE 22/24]   '${feature}' not in ${joined} -- skipping barplot."
                continue
            fi
            humann_barplot \
                --input "${joined}" \
                --output "${out}/barplot_bins_${suffix}_${feature}.png" \
                --focal-feature "${feature}" \
                --sort sum
        done
    done

    conda deactivate || true
    touch "${marker}"
}

################################################################################
# LEVELS 2-4: run only by the last-arriving array task
################################################################################

run_stage22_combine() {

echo "$(date): [STAGE 22/24] Starting cross-sample combining..."

conda activate humann3-env

READS_ROOT="$(dirname "${HUMANN_READS_DIR}")"
CONTIGS_ROOT="$(dirname "${HUMANN_CONTIGS_DIR}")"

mkdir -p "${WHOLE_OUT}" "${SPECIES_OUT}" "${AGG_OUT}"

# Segatella check, first thing: which bins does MetaPhlAn call Segatella?
echo "$(date): [STAGE 22/24] Segatella check (MetaPhlAn-based)..."
build_segatella_assignments || echo "$(date): [STAGE 22/24]   WARNING: Segatella check failed -- no Segatella bins will be assigned."
[ -s "${SEG_ASSIGN}" ] || printf 'sample\tbin\tgtdbtk_species\tmetaphlan_species\tmetaphlan_pct\tsource\teffective_species\n' > "${SEG_ASSIGN}"

# build_whole_sample_table <suffix>
# Joins the whole-sample (STAGE_20) tables of this type across samples into
# ${WHOLE_OUT}/joined_<suffix>.tsv, if one doesn't already exist. Returns
# non-zero if no whole-sample counterpart exists for this type.
build_whole_sample_table() {
    local suffix="$1"
    local joined="${WHOLE_OUT}/joined_${suffix}.tsv"
    [ -s "${joined}" ] && return 0

    local stage="${WHOLE_OUT}/staging_${suffix}"
    rm -rf "${stage}"
    mkdir -p "${stage}"

    local n=0 s hits
    while read -r s; do
        [ -z "${s}" ] && continue
        mapfile -t hits < <(find "${READS_ROOT}/${s}" -maxdepth 1 -name "*_${suffix}.tsv" 2>/dev/null)
        if [ "${#hits[@]}" -eq 1 ]; then
            cp "${hits[0]}" "${stage}/${s}_${suffix}.tsv"
            n=$((n + 1))
        elif [ "${#hits[@]}" -gt 1 ]; then
            echo "$(date): [STAGE 22/24]   WARNING: ${#hits[@]} whole-sample '${suffix}' files for ${s} -- skipping that sample."
        fi
    done < "${DIRECTORY}/samples.txt"

    [ "${n}" -eq 0 ] && return 1
    humann_join_tables -i "${stage}" -o "${joined}" --file_name "${suffix}" || { rm -f "${joined}"; return 1; }
    [ -s "${joined}" ]
}

################################################################################
# PART A (LEVEL 2): whole-sample (STAGE_20) joins + metadata barplots
################################################################################

echo "$(date): [STAGE 22/24] Part A: joining whole-sample tables from ${READS_ROOT}..."

# Raw (unnormalised), STRATIFIED pathabundance -- humann_barplot needs the
# taxon-stratified sub-rows to plot per-taxon contributions.
# Joined via per-type staging (flat folder of one file per sample), NOT by
# pointing humann_join_tables at the parent directory: the recursive form
# produced no output here, while the staged form works.
for WS_TYPE in pathabundance genefamilies; do
    if build_whole_sample_table "${WS_TYPE}"; then
        echo "$(date): [STAGE 22/24]   Whole-sample ${WS_TYPE} joined -> ${WHOLE_OUT}/joined_${WS_TYPE}.tsv"
    else
        FIRST_SAMPLE=$(head -1 "${DIRECTORY}/samples.txt")
        echo "$(date): [STAGE 22/24]   ERROR: could not build whole-sample '${WS_TYPE}' table." >&2
        echo "$(date): [STAGE 22/24]   Looked for ${READS_ROOT}/<SAMPLE>/*_${WS_TYPE}.tsv. Files in ${READS_ROOT}/${FIRST_SAMPLE}:" >&2
        ls "${READS_ROOT}/${FIRST_SAMPLE}" 2>&1 | head -30 >&2
    fi
done

if [ -s "${WHOLE_OUT}/joined_pathabundance.tsv" ]; then
    barplots_by_metadata \
        "${WHOLE_OUT}/joined_pathabundance.tsv" \
        "${WHOLE_OUT}" \
        "pathabund" \
        "${FOCAL_FEATURES_GENERAL[@]}"
else
    echo "$(date): [STAGE 22/24]   No joined_pathabundance.tsv -- skipping whole-sample metadata barplots."
fi

# All pathways, whole-sample, species = HUMAnN's own assignment
allpath_species_plots "${WHOLE_OUT}/joined_pathabundance.tsv" "${WHOLE_OUT}" "whole_sample_humann"

# Custom pathways, whole-sample (HUMAnN's own species assignment)
if build_whole_sample_table "ec-cpm-named" \
   && make_custom_pathway_table "${WHOLE_OUT}/joined_ec-cpm-named.tsv" "${WHOLE_OUT}/joined_custom-pathways.tsv" "${WHOLE_OUT}/joined_custom-pathways_ec_breakdown.tsv"; then
    barplots_by_metadata "${WHOLE_OUT}/joined_custom-pathways.tsv" "${WHOLE_OUT}" "custom_pathways" "${CUSTOM_PATHWAY_IDS[@]}"
else
    echo "$(date): [STAGE 22/24]   No whole-sample ec-cpm-named table (or no custom-pathway ECs in it) -- skipping whole-sample custom pathways."
fi


################################################################################
# PART A2 (LEVEL 3): bin-assigned aggregate, all table types, whole-sample
# normalised
################################################################################

echo "$(date): [STAGE 22/24] Part A2: bin-assigned aggregate..."

# Table types = whatever Level 1 produced across all samples
mapfile -t AGG_TYPES < <(find "${BINS_OUT}" -mindepth 2 -maxdepth 2 -name 'joined_bins_*.tsv' -printf '%f\n' \
    | sed -e 's/^joined_bins_//' -e 's/\.tsv$//' | sort -u)
echo "$(date): [STAGE 22/24]   ${#AGG_TYPES[@]} table type(s) to aggregate: ${AGG_TYPES[*]:-none}"

# The aggregates depend on the species assignment. If it changed (or these
# aggregates predate the MetaPhlAn Segatella step) discard them so they are
# rebuilt rather than silently reused.
STAMP_NOW="$(md5sum "${SEG_ASSIGN}" | cut -d' ' -f1)-sp${MANIFEST_SPECIES_COL}-g${MANIFEST_GENUS_COL:-none}-scaled2"
STAMP_FILE="${AGG_OUT}/.assignment_stamp"
if [ "$(cat "${STAMP_FILE}" 2>/dev/null)" != "${STAMP_NOW}" ]; then
    echo "$(date): [STAGE 22/24]   Species assignment inputs changed (or first run of this version) -- discarding old aggregates."
    rm -f "${AGG_OUT}"/bin_aggregate_*.tsv "${AGG_OUT}"/bin_scale_factors.tsv \
          "${AGG_OUT}"/hmp_*.pcl "${AGG_OUT}"/barplot_* "${AGG_OUT}"/allpathways_*
fi
echo "${STAMP_NOW}" > "${STAMP_FILE}"

# Each bin's share of its sample (for CPM / normalised tables). Cached; rebuilt
# when the stamp changes.
BIN_SCALE_FACTORS="${AGG_OUT}/bin_scale_factors.tsv"
if [ ! -s "${BIN_SCALE_FACTORS}" ]; then
    echo "$(date): [STAGE 22/24]   Building per-bin scale factors..."
    build_bin_scale_factors || echo "$(date): [STAGE 22/24]   WARNING: scale factor step failed; CPM-type aggregates will fail."
fi

for SUFFIX in "${AGG_TYPES[@]}"; do
    UNIT="$(bin_unit_class "${SUFFIX}")"
    WHOLE_TBL="-"
    if [ "${UNIT}" = "coverage" ]; then
        MODE="maxbins"
    elif build_whole_sample_table "${SUFFIX}"; then
        WHOLE_TBL="${WHOLE_OUT}/joined_${SUFFIX}.tsv"
        MODE="scaled"
    else
        MODE="scaledbins"
        echo "$(date): [STAGE 22/24]   No whole-sample '${SUFFIX}' table -- aggregate is bins only (no unbinned remainder) for this type."
    fi
    echo "$(date): [STAGE 22/24]   ${SUFFIX}: unit class ${UNIT}, mode ${MODE}"

    AGG_FILE="${AGG_OUT}/bin_aggregate_${SUFFIX}_${MODE}.tsv"

    # An existing aggregate is reused (the gene-family ones are slow). Delete
    # ${AGG_FILE} to force a rebuild, e.g. after changing the manifest columns.
    if [ -s "${AGG_FILE}" ]; then
        echo "$(date): [STAGE 22/24]   ${AGG_FILE} already exists -- reusing."
    else
    python3 - "${AGG_FILE}" "${BINS_OUT}" "${CONTIGS_ROOT}" "${SUFFIX}" "${DIRECTORY}/samples.txt" \
              "${WHOLE_TBL}" "${MANIFEST_SPECIES_COL}" "${MANIFEST_GENUS_COL}" "${SEG_ASSIGN}" \
              "${UNIT}" "${BIN_SCALE_FACTORS}" << 'PYEOF'
import csv, os, re, sys
from collections import OrderedDict, defaultdict

(out_path, bins_root, contigs_root, suffix, samples_file,
 whole_path, species_col, genus_col, seg_assign, unit, factors_path) = sys.argv[1:12]

species_col = int(species_col) - 1
genus_col = int(genus_col) - 1 if genus_col.strip() else None
have_whole = whole_path not in ("", "-") and os.path.isfile(whole_path)

# MetaPhlAn-identified Segatella bins: (sample, bin) -> MetaPhlAn species, which
# replaces the GTDB-Tk species for those bins.
override = {}
if seg_assign and os.path.isfile(seg_assign):
    with open(seg_assign) as fh:
        rd_o = csv.reader(fh, delimiter="\t")
        next(rd_o, None)
        for r in rd_o:
            if len(r) >= 7 and r[5] == "metaphlan":
                override[(r[0], r[1])] = r[6]
n_overridden = 0

# Per-bin share of the sample, for CPM / normalised tables
scale_kind = unit.split(":", 1)[1] if unit.startswith("scale:") else None
factor = {}
if scale_kind:
    if not os.path.isfile(factors_path):
        sys.exit("ERROR: scale factors file {} missing".format(factors_path))
    with open(factors_path) as fh:
        rd_f = csv.reader(fh, delimiter="\t")
        next(rd_f, None)
        for r in rd_f:
            if len(r) >= 6 and r[2] == scale_kind:
                try:
                    factor[(r[0], r[1])] = float(r[5])
                except ValueError:
                    pass
n_nofactor = 0

with open(samples_file) as fh:
    samples = [l.strip() for l in fh if l.strip()]

UNKNOWN = {"", "unclassified", "unknown", "na", "n/a", "none"}

def clean(s):
    return re.sub(r"[^A-Za-z0-9_]", "_", s.strip())

def taxon_stratum(species_val, genus_val=""):
    """GTDB-Tk call -> HUMAnN-style 'g__Genus.s__Genus_species'.
       - full species         -> g__Genus.s__Genus_species
       - genus only (species unknown) -> its own pseudo-species: g__Genus.s__Genus_sp
       - nothing known        -> g__Unclassified.s__Unclassified_bin
    Accepts 'Genus species', 's__Genus species', or a full GTDB lineage string."""
    sp = (species_val or "").strip()
    ge = (genus_val or "").strip()
    if ";" in sp or "g__" in sp:
        parts = {}
        for tok in sp.split(";"):
            tok = tok.strip()
            if tok[:3] in ("g__", "s__"):
                parts[tok[0]] = tok[3:].strip()
        ge = parts.get("g", ge)
        sp = parts.get("s", "")
    if sp.startswith("s__"):
        sp = sp[3:].strip()
    if ge.startswith("g__"):
        ge = ge[3:].strip()
    if sp.lower() in UNKNOWN:
        sp = ""
    if ge.lower() in UNKNOWN:
        ge = ""
    if sp and len(sp.split()) == 1:          # a lone word is a genus, not a species
        if not ge:
            ge = sp
        sp = ""
    if sp and not ge:
        ge = sp.split()[0]
    if not sp and not ge:
        return "g__Unclassified.s__Unclassified_bin"
    if not sp:
        return "g__{0}.s__{0}_sp".format(clean(ge))
    return "g__{}.s__{}".format(clean(ge), clean(sp))

def match_label(col, labels):
    """Longest label found in the column header as a whole token (not flanked
    by letters/digits). Handles 'S115_Abundance' and
    'out.S115_combined_Abundance-RPKs'; S11 never matches S115, bin.1 never
    matches bin.10."""
    best = None
    for lab in labels:
        if lab and re.search(r"(?<![A-Za-z0-9])" + re.escape(lab) + r"(?![A-Za-z0-9])", col):
            if best is None or len(lab) > len(best):
                best = lab
    return best

def feat_key(feat):
    return feat.split(":", 1)[0].strip()

names = {}                                    # feature key -> display row name

# ---- whole-sample (Level 2) values: sample -> key -> value ------------------
whole = defaultdict(dict)
if have_whole:
    with open(whole_path) as fh:
        rd = csv.reader(fh, delimiter="\t")
        header = next(rd)
        colmap = {}
        for ci, col in enumerate(header[1:], start=1):
            s = match_label(col, samples)
            if s is None:
                sys.stderr.write("WARNING: whole-sample column '{}' matches no sample ID\n".format(col))
            else:
                colmap[ci] = s
        for r in rd:
            if not r or "|" in r[0]:
                continue
            key = feat_key(r[0])
            names.setdefault(key, r[0])
            for ci, s in colmap.items():
                if ci < len(r):
                    try:
                        v = float(r[ci])
                    except ValueError:
                        continue
                    if v != 0.0:
                        whole[s][key] = v

if have_whole and not colmap:
    sys.exit("ERROR: none of the {} whole-sample columns in {} matched a sample ID in "
             "samples.txt -- refusing to write an empty/mislabelled aggregate.".format(
                 len(header) - 1, whole_path))

order = OrderedDict()                         # feature keys in first-seen order
totals = defaultdict(lambda: defaultdict(float))                  # key -> sample -> v
strata = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))  # key -> tag -> sample -> v
used_samples = []
n_not_in_whole = 0
n_over = 0
tot_binned = 0.0
tot_total = 0.0

for sample in samples:
    joined = os.path.join(bins_root, sample, "joined_bins_{}.tsv".format(suffix))
    manifest = os.path.join(contigs_root, sample, "{}_bin_humann_manifest.tsv".format(sample))
    if not (os.path.isfile(joined) and os.path.getsize(joined) > 0):
        sys.stderr.write("WARNING: no per-bin joined table for {} ({})\n".format(sample, joined))
        continue

    bin_taxon = {}
    if os.path.isfile(manifest):
        with open(manifest) as fh:
            rd = csv.reader(fh, delimiter="\t")
            next(rd, None)
            for r in rd:
                if len(r) > species_col:
                    g = r[genus_col] if (genus_col is not None and len(r) > genus_col) else ""
                    if (sample, r[0]) in override:
                        bin_taxon[r[0]] = taxon_stratum(override[(sample, r[0])], "")
                        n_overridden += 1
                    else:
                        bin_taxon[r[0]] = taxon_stratum(r[species_col], g)
    else:
        sys.stderr.write("WARNING: no manifest for {} -- all its bins -> unclassified\n".format(sample))

    bin_vals = defaultdict(lambda: defaultdict(float))            # key -> tag -> v
    with open(joined) as fh:
        rd = csv.reader(fh, delimiter="\t")
        header = next(rd)
        col_tag = {}
        col_fac = {}
        for ci, col in enumerate(header[1:], start=1):
            lab = match_label(col, list(bin_taxon))
            if lab is None:
                sys.stderr.write("WARNING: {}: column '{}' matches no manifest bin -> unclassified\n".format(sample, col))
                tag = "g__Unclassified.s__Unclassified_bin"
            else:
                tag = bin_taxon[lab]
            fac = 1.0
            if scale_kind:
                fac = factor.get((sample, lab)) if lab is not None else None
                if fac is None:
                    n_nofactor += 1            # cannot put this bin on sample scale
                    continue
            col_tag[ci] = tag
            col_fac[ci] = fac
        for r in rd:
            if not r or "|" in r[0]:       # ignore HUMAnN's own per-bin strata
                continue
            key = feat_key(r[0])
            names.setdefault(key, r[0])
            for ci, tag in col_tag.items():
                if ci >= len(r):
                    continue
                try:
                    v = float(r[ci])
                except ValueError:
                    continue
                if v != 0.0:
                    v *= col_fac[ci]
                    if unit == "coverage":
                        if v > bin_vals[key][tag]:
                            bin_vals[key][tag] = v
                    else:
                        bin_vals[key][tag] += v

    used_samples.append(sample)

    if unit == "coverage":
        # coverage is not additive: strata and total are the maximum over bins
        for key, tagd in bin_vals.items():
            order.setdefault(key, None)
            for tag, v in tagd.items():
                strata[key][tag][sample] = v
            totals[key][sample] = max(tagd.values())
        continue

    # Bin values are on the sample's scale (RPK-type: as is; CPM-type: scaled by
    # the bin's share). What the bins do not account for is the remainder
    # between them and the whole-sample value, reported as "unclassified".
    keys = set(bin_vals)
    if have_whole:
        keys |= set(whole[sample])
    for key in keys:
        tagd = bin_vals.get(key, {})
        binned = sum(tagd.values())
        total = binned
        if have_whole:
            W = whole[sample].get(key)
            if W is None:
                if binned > 0:
                    n_not_in_whole += 1        # in a bin, zero/absent in the whole sample
            elif binned > W:
                n_over += 1                    # binned exceeds whole: keep binned, no remainder
            else:
                rest = W - binned
                if rest > 0:
                    strata[key]["unclassified"][sample] += rest
                total = W
        if total <= 0:
            continue
        order.setdefault(key, None)
        for tag, v in tagd.items():
            strata[key][tag][sample] += v
        totals[key][sample] = total
        tot_binned += binned
        tot_total += total

if not used_samples or not order:
    sys.stderr.write("No data aggregated for {} -- nothing written.\n".format(suffix))
    sys.exit(1)

with open(out_path, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(["# Feature"] + used_samples)
    for key in order:
        name = names.get(key, key)
        w.writerow([name] + [repr(totals[key].get(s, 0.0)) for s in used_samples])
        for tag in sorted(strata[key]):
            w.writerow(["{}|{}".format(name, tag)] +
                       [repr(strata[key][tag].get(s, 0.0)) for s in used_samples])

sys.stderr.write("MetaPhlAn Segatella species applied to {} bin(s).\n".format(n_overridden))
sys.stderr.write("Wrote {} ({} samples, {} features; unit={}, whole-sample remainder={}).\n".format(
    out_path, len(used_samples), len(order), unit, "yes" if have_whole else "no"))
if have_whole:
    sys.stderr.write("  {} per-sample feature(s) in bins but absent from the whole-sample table (kept as binned only); "
                     "{} where binned exceeded the whole-sample value (kept, no remainder).\n".format(n_not_in_whole, n_over))
    if tot_total > 0:
        sys.stderr.write("  Binned share of whole-sample abundance (summed over features): {:.1%}\n".format(tot_binned / tot_total))
if n_nofactor:
    sys.stderr.write("WARNING: {} bin column(s) skipped: no scale factor.\n".format(n_nofactor))
PYEOF
    fi

    if [ ! -s "${AGG_FILE}" ]; then
        echo "$(date): [STAGE 22/24]   Aggregate for ${SUFFIX} not produced -- skipping."
        continue
    fi

    # Barplots only for table types with focal features defined
    AGG_FEATS=()
    read -ra AGG_FEATS <<< "${BIN_TABLE_FOCAL[$SUFFIX]:-}"
    if [ "${#AGG_FEATS[@]}" -gt 0 ]; then
        barplots_by_metadata "${AGG_FILE}" "${AGG_OUT}" "binagg_${SUFFIX}_${MODE}" "${AGG_FEATS[@]}"
    else
        echo "$(date): [STAGE 22/24]   No focal features for '${SUFFIX}' -- aggregate table written, no barplots."
    fi

    # All pathways, species = GTDB-Tk species of the source bin
    if [ "${SUFFIX}" = "pathabundance" ]; then
        allpath_species_plots "${AGG_FILE}" "${AGG_OUT}" "bin_assigned_${MODE}"
    fi

    # Custom pathways, species assigned by source bin
    if [ "${SUFFIX}" = "ec-cpm-named" ]; then
        CUSTOM_AGG="${AGG_OUT}/bin_aggregate_custom-pathways_${MODE}.tsv"
        CUSTOM_AGG_BREAKDOWN="${AGG_OUT}/bin_aggregate_custom-pathways_${MODE}_ec_breakdown.tsv"
        if make_custom_pathway_table "${AGG_FILE}" "${CUSTOM_AGG}" "${CUSTOM_AGG_BREAKDOWN}"; then
            barplots_by_metadata "${CUSTOM_AGG}" "${AGG_OUT}" "binagg_custom_pathways_${MODE}" "${CUSTOM_PATHWAY_IDS[@]}"
        else
            echo "$(date): [STAGE 22/24]   No custom-pathway ECs in ${AGG_FILE} -- skipping bin-assigned custom pathways."
        fi
    fi
done

################################################################################
# PART B (LEVEL 4): per-species (STAGE_21 bins) joins + barplots
################################################################################
#
# Target species / genus. Matching = GTDB-Tk species call from each sample's
# STAGE_21 manifest (awk ERE on column 2, status=ok), plus an optional MetaPhlAn
# fallback (regex, case-insensitive, on column 4 of STAGE_12's crosscheck
# summary) for bins GTDB-Tk doesn't label as wanted. Keep in sync with the R
# visualization script's TARGET_TAXA. Avoid backslashes in these patterns
# (awk -v would interpret them).
#
# Every species gets (1) the general focal pathways on the pathabundance table,
# (2) its own species-specific features (below) on the EC table, and (3) the
# "all pathways" variants. The tables are always joined, whether or not there
# is anything to plot.
#
# HMO / bifidobacterium_infantis is resolved at EC level, not CAZy-family level
# (no confirmed MetaCyc pathway for GH29/95/33/2/42/20/136/112 HMO degradation):
#   GH29, GH95 (alpha-L-fucosidase)      -> EC 3.2.1.51
#   GH33       (sialidase)               -> EC 3.2.1.18
#   GH2, GH42  (beta-galactosidase)      -> EC 3.2.1.23
#   GH20       (beta-hexosaminidase)     -> EC 3.2.1.52
#   GH136      (lacto-N-biosidase)       -> EC 3.2.1.140
#   GH112      (GNB/LNB phosphorylase)   -> EC 2.4.1.211
# A missing EC row means the activity wasn't detected via the UniRef/EC route
# in any bin (not necessarily absent biologically); dbCAN (STAGE_13C) is the
# more direct source.

# Where each target's bins come from:
#   gtdbtk    : regex on the GTDB-Tk species (manifest column 2), plus the
#               optional MetaPhlAn fallback below
#   segatella : rows of ${SEG_ASSIGN} (MetaPhlAn-based, built at the start); the
#               regex is applied to that file's effective_species column, so
#               species-level sets use MetaPhlAn's species name
declare -A TARGET_SPECIES_SOURCE=(
    ["bifidobacterium_infantis"]="gtdbtk"
    ["segatella_copri"]="segatella"
    ["segatella_sinensis"]="segatella"
    ["segatella_genus"]="segatella"
    ["escherichia_coli"]="gtdbtk"
)
declare -A TARGET_SPECIES_GTDBTK_REGEX=(
    ["bifidobacterium_infantis"]="^Bifidobacterium (infantis|longum_A)$"
    ["segatella_copri"]="^(Segatella|Prevotella) copri$"
    ["segatella_sinensis"]="^(Segatella|Prevotella) sinensis$"
    ["segatella_genus"]="^(Segatella .+|Prevotella copri.*)$"
    ["escherichia_coli"]="^Escherichia coli(_[A-Za-z])?$"
)
declare -A TARGET_SPECIES_METAPHLAN_FALLBACK=(
    ["escherichia_coli"]="^escherichia coli$"
)
# Species-specific EC features (space-separated); empty = none
declare -A TARGET_SPECIES_FOCAL_FEATURES=(
    ["bifidobacterium_infantis"]="${HMO_EC_FEATURES}"
    ["segatella_copri"]=""
    ["segatella_sinensis"]=""
    ["segatella_genus"]=""
    ["escherichia_coli"]=""
)

for SPECIES_KEY in "${!TARGET_SPECIES_GTDBTK_REGEX[@]}"; do
    echo "$(date): [STAGE 22/24] Part B: ${SPECIES_KEY}..."
    SPECIES_DIR="${SPECIES_OUT}/${SPECIES_KEY}"
    mkdir -p "${SPECIES_DIR}"

    GTDBTK_REGEX="${TARGET_SPECIES_GTDBTK_REGEX[$SPECIES_KEY]}"
    MP_FALLBACK="${TARGET_SPECIES_METAPHLAN_FALLBACK[$SPECIES_KEY]:-}"

    # ---- 1. which (sample, bin) pairs belong to this species ----
    MATCH_PAIRS=()
    if [ "${TARGET_SPECIES_SOURCE[$SPECIES_KEY]}" = "segatella" ]; then
        # MetaPhlAn-based Segatella assignment (built at the start of the combine)
        mapfile -t MATCH_PAIRS < <(awk -F'\t' -v pat="${GTDBTK_REGEX}" -v also="${SEGATELLA_ALSO_GTDB}" \
            'NR>1 && $7 ~ pat && ($6=="metaphlan" || (also==1 && $6=="gtdbtk_only")) {print $1 "\t" $2}' \
            "${SEG_ASSIGN}")
    else
        for SAMPLE_DIR in "${CONTIGS_ROOT}"/*/; do
            SAMPLE_NAME=$(basename "${SAMPLE_DIR}")
            MANIFEST="${SAMPLE_DIR}/${SAMPLE_NAME}_bin_humann_manifest.tsv"
            [ -s "${MANIFEST}" ] || continue

            MATCHED_BINS=$(awk -F'\t' -v pat="${GTDBTK_REGEX}" \
                'NR>1 && $4=="ok" && $2 ~ pat {print $1}' "${MANIFEST}")

            if [ -n "${MP_FALLBACK}" ]; then
                CROSSCHECK_FILE="${CROSSCHECK_OUT}/${SAMPLE_NAME}_MAG_metaphlan_crosscheck.tsv"
                if [ -s "${CROSSCHECK_FILE}" ]; then
                    FALLBACK_BINS=$(awk -F'\t' -v want="${MP_FALLBACK}" \
                        'NR>1 && tolower($4) ~ tolower(want) {print $1}' "${CROSSCHECK_FILE}")
                    for b in ${FALLBACK_BINS}; do
                        if ! grep -qx "${b}" <<< "${MATCHED_BINS}"; then
                            MATCHED_BINS="${MATCHED_BINS}"$'\n'"${b}"
                        fi
                    done
                fi
            fi

            for BIN_LABEL in ${MATCHED_BINS}; do
                [ -z "${BIN_LABEL}" ] && continue
                MATCH_PAIRS+=("${SAMPLE_NAME}"$'\t'"${BIN_LABEL}")
            done
        done
    fi

    if [ "${#MATCH_PAIRS[@]}" -eq 0 ]; then
        echo "$(date): [STAGE 22/24]   No bins matched ${SPECIES_KEY} -- skipping."
        continue
    fi
    printf '%s\n' "${MATCH_PAIRS[@]}" > "${SPECIES_DIR}/matched_bins.tsv"
    echo "$(date): [STAGE 22/24]   ${#MATCH_PAIRS[@]} bin(s) matched ${SPECIES_KEY} (list: ${SPECIES_DIR}/matched_bins.tsv)"

    # ---- 2. per table type: stage, join, plot ----
    for TABLE_SUFFIX in "${SPECIES_TABLE_TYPES[@]}"; do
        STAGE_DIR="${SPECIES_DIR}/staging_${TABLE_SUFFIX}"
        rm -rf "${STAGE_DIR}"
        mkdir -p "${STAGE_DIR}"

        N_MATCHED=0
        for PAIR in "${MATCH_PAIRS[@]}"; do
            SAMPLE_NAME="${PAIR%%$'\t'*}"
            BIN_LABEL="${PAIR#*$'\t'}"
            SRC="${CONTIGS_ROOT}/${SAMPLE_NAME}/${BIN_LABEL}_${TABLE_SUFFIX}.tsv"
            if [ -s "${SRC}" ]; then
                cp "${SRC}" "${STAGE_DIR}/${SAMPLE_NAME}_${BIN_LABEL}_${TABLE_SUFFIX}.tsv"
                N_MATCHED=$((N_MATCHED + 1))
            fi
        done

        if [ "${N_MATCHED}" -eq 0 ]; then
            echo "$(date): [STAGE 22/24]   No ${TABLE_SUFFIX} tables found for ${SPECIES_KEY} -- skipping this type."
            continue
        fi

        JOINED_FILE="${SPECIES_DIR}/joined_${SPECIES_KEY}_${TABLE_SUFFIX}.tsv"
        if ! humann_join_tables -i "${STAGE_DIR}" -o "${JOINED_FILE}" --file_name "${TABLE_SUFFIX}"; then
            echo "$(date): [STAGE 22/24]   WARNING: humann_join_tables failed for ${SPECIES_KEY} ${TABLE_SUFFIX} -- skipping."
            continue
        fi
        echo "$(date): [STAGE 22/24]   Joined ${N_MATCHED} bin(s) (${TABLE_SUFFIX}) for ${SPECIES_KEY} -> ${JOINED_FILE}"

        # Focal features: general pathways on pathabundance; species-specific
        # features on the EC table
        FOCAL_FEATURE_LIST=()
        if [ "${TABLE_SUFFIX}" = "pathabundance" ]; then
            FOCAL_FEATURE_LIST=("${FOCAL_FEATURES_GENERAL[@]}")
        else
            read -ra FOCAL_FEATURE_LIST <<< "${TARGET_SPECIES_FOCAL_FEATURES[$SPECIES_KEY]:-}"
        fi

        for FOCAL_FEATURE in "${FOCAL_FEATURE_LIST[@]}"; do
            if ! feature_present "${JOINED_FILE}" "${FOCAL_FEATURE}"; then
                echo "$(date): [STAGE 22/24]   '${FOCAL_FEATURE}' not in ${JOINED_FILE} -- skipping barplot."
                continue
            fi
            echo "$(date): [STAGE 22/24]   humann_barplot ${SPECIES_KEY} x ${FOCAL_FEATURE}..."
            humann_barplot \
                --input "${JOINED_FILE}" \
                --output "${SPECIES_DIR}/barplot_${SPECIES_KEY}_${FOCAL_FEATURE}.png" \
                --focal-feature "${FOCAL_FEATURE}" \
                --sort sum
        done

        # Custom pathways (summed ECs) for this species' bins
        if [ "${TABLE_SUFFIX}" = "ec-cpm-named" ]; then
            CUSTOM_JOINED="${SPECIES_DIR}/joined_${SPECIES_KEY}_custom-pathways.tsv"
            CUSTOM_BREAKDOWN="${SPECIES_DIR}/joined_${SPECIES_KEY}_custom-pathways_ec_breakdown.tsv"
            if make_custom_pathway_table "${JOINED_FILE}" "${CUSTOM_JOINED}" "${CUSTOM_BREAKDOWN}"; then
                for CP in "${CUSTOM_PATHWAY_IDS[@]}"; do
                    feature_present "${CUSTOM_JOINED}" "${CP}" || continue
                    echo "$(date): [STAGE 22/24]   humann_barplot ${SPECIES_KEY} x ${CP}..."
                    humann_barplot \
                        --input "${CUSTOM_JOINED}" \
                        --output "${SPECIES_DIR}/barplot_${SPECIES_KEY}_${CP}.png" \
                        --focal-feature "${CP}" \
                        --sort sum
                done
            else
                echo "$(date): [STAGE 22/24]   No custom-pathway ECs in ${JOINED_FILE} -- no custom pathway plots for ${SPECIES_KEY}."
            fi
        fi

        # All pathways (pathway table only): stacked by HUMAnN taxon, and by pathway
        if [ "${TABLE_SUFFIX}" = "pathabundance" ]; then
            for STACK in species pathway; do
                echo "$(date): [STAGE 22/24]   All-pathways plot ${SPECIES_KEY} (stacked by ${STACK})..."
                plot_all_pathways --abundance "${JOINED_FILE}" \
                    --out-prefix "${SPECIES_DIR}/allpathways_${SPECIES_KEY}_by-${STACK}" \
                    --stack-by "${STACK}" \
                    --title "${SPECIES_KEY}: all pathways, stacked by ${STACK}" \
                    || echo "$(date): [STAGE 22/24]   WARNING: all-pathways plot failed (${SPECIES_KEY}, ${STACK})"
            done
        fi
    done
done

conda deactivate || true

}   # end run_stage22_combine

################################################################################
# DRIVER
################################################################################

mkdir -p "${BINS_OUT}"

# LEVEL 1 runs for every sample, BEFORE check-in, so that the last arrival can
# rely on every sample's per-bin joined tables already existing.
run_stage22_per_sample_bins

# LEVELS 2-4: last-arrival check-in (same pattern as STAGE_17).
# ${SLURM_ARRAY_JOB_ID:-manual} so this also works outside an array.
ARRAY_ID="${SLURM_ARRAY_JOB_ID:-manual}"

echo "$(date): [STAGE 22/24] Combine check-in for ${SAMPLE}..."

# Same idea for the shared flag: a flag without its outputs is stale.
if [ -f "${STATUS_DIR}/humann_combine_complete.txt" ] && [ ! -s "${WHOLE_OUT}/joined_pathabundance.tsv" ]; then
    echo "$(date): [STAGE 22/24] Completion flag exists but combined outputs are missing -- ignoring stale flag."
    rm -f "${STATUS_DIR}/humann_combine_complete.txt"
fi

if [ -f "${STATUS_DIR}/humann_combine_complete.txt" ]; then
    echo "$(date): [STAGE 22/24] Combine already complete for this batch -- skipping."
else
    touch "${STATUS_DIR}/stage22_arrived.${ARRAY_ID}.${SAMPLE}"
    TOTAL_SAMPLES=$(wc -l < "${DIRECTORY}/samples.txt")
    ARRIVED=$(find "${STATUS_DIR}" -maxdepth 1 -name "stage22_arrived.${ARRAY_ID}.*" | wc -l)
    echo "$(date): [STAGE 22/24] ${ARRIVED}/${TOTAL_SAMPLES} samples have reached STAGE_22."

    if [ "$ARRIVED" -ge "$TOTAL_SAMPLES" ] && mkdir "${STATUS_DIR}/stage22_combine.lock.${ARRAY_ID}" 2>/dev/null; then
        echo "$(date): [STAGE 22/24] ${SAMPLE} is the last arrival -- running combine..."
        run_stage22_combine
        touch "${STATUS_DIR}/humann_combine_complete.txt"
        echo "$(date): [STAGE 22/24] Combining complete."
    else
        echo "$(date): [STAGE 22/24] Not the last arrival (or lock already claimed) -- deferring combine."
    fi
fi

echo "STAGE_23" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 22/24] Done for ${SAMPLE}."