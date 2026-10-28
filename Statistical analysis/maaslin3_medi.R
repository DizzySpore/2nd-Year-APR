#!/usr/bin/env Rscript
# =============================================================================
# Maaslin3 Analysis: MetaPhlAn bacteria + Bracken archaea + HUMAnN + MEDI
# =============================================================================
# Taxonomic / functional sources:
#   Bacteria (species + genus) : MetaPhlAn *_metaphlan_profile_wreads.txt
#   Archaea                    : Bracken   *_archaea_bracken_abundance.c01.txt
#   Functional (HUMAnN3, CPM)  : humann_dir/SXXX/out.SXXX_combined_pathabundance_normalised.tsv
#                                out.SXXX_combined_genefamilies-cpm.tsv
#                                out.SXXX_combined_ec-cpm-named.tsv
#                                out.SXXX_combined_rxn-cpm-named.tsv
#                                out.SXXX_combined_ko-cpm-named.tsv
#   Food abundance (MEDI)      : food_abundance.csv  (organisms detected in the diet)
#   NB: all MEDI tables keep ONLY sample IDs matching S<number>; other datasets are dropped.
#   Sample exclusion           : WHO anthropometric survey CSV; every child with a
#                                non-empty 'Flag' is dropped from all analyses
#
# Dietary sources (MEDI):
#   food_content.csv  -> nutrients / compounds
#
# Analyses (all optional, controlled by --analyses):
#   bacteria_species_metadata   bacterial species  ~ metadata
#   bacteria_genus_metadata     bacterial genera   ~ metadata
#   archaea_metadata            archaeal species   ~ metadata
#   pathways_metadata           HUMAnN pathways    ~ metadata
#   genefamilies_metadata       HUMAnN gene fams   ~ metadata
#   reactions_metadata          HUMAnN reactions   ~ metadata
#   ko_metadata                 HUMAnN KO groups   ~ metadata
#   ec_metadata                 HUMAnN EC numbers  ~ metadata
#   pathways_strat_metadata     species-assigned pathways  ~ metadata
#   ec_strat_metadata / reactions_strat_metadata / ko_strat_metadata   (same, other layers)
#   nutrient_vs_pathways_strat  species-assigned pathways  ~ nutrients (univariate)
#   nutrient_vs_ec_strat / nutrient_vs_reactions_strat / nutrient_vs_ko_strat
#   food_metadata               food abundance    ~ metadata
#   food_vs_bacteria_species/genus, food_vs_archaea   taxa ~ each food (univariate)
#   food_vs_pathways/ec/reactions/ko (+ _strat)       functions ~ each food (univariate)
#   nutrient_vs_ec              EC abundances      ~ nutrients (univariate)
#   ec_vs_bacteria_species / ec_vs_bacteria_genus / ec_vs_archaea
#   nutrient_metadata           nutrients          ~ metadata
#   nutrient_bacteria_species   bacterial species  ~ nutrients + covariates
#   nutrient_bacteria_genus     bacterial genera   ~ nutrients + covariates
#   nutrient_archaea            archaeal species   ~ nutrients + covariates
#   nutrient_vs_pathways        pathway abundances ~ nutrients (univariate)
#   nutrient_vs_reactions       reactions          ~ nutrients (univariate)
#   nutrient_vs_ko              KO abundances      ~ nutrients (univariate)
#   pathways_vs_bacteria_species bacterial species ~ pathways  (univariate)
#   pathways_vs_bacteria_genus  bacterial genera   ~ pathways  (univariate)
#   pathways_vs_archaea         archaeal species   ~ pathways  (univariate)
#   reactions_vs_bacteria_species bacterial species ~ reactions (univariate)
#   reactions_vs_bacteria_genus  bacterial genera  ~ reactions (univariate)
#   reactions_vs_archaea         archaeal species  ~ reactions (univariate)
#   ko_vs_bacteria_species      bacterial species  ~ KOs       (univariate)
#   ko_vs_bacteria_genus        bacterial genera   ~ KOs       (univariate)
#   ko_vs_archaea               archaeal species   ~ KOs       (univariate)
#
# Usage:
#   Rscript maaslin3_medi.R \
#     --metaphlan_dir   metaphlan_wreads_output/ \
#     --archaea_dir     kraken2_bracken_out/ \
#     --humann_dir      humann_reads_output/ \
#     --food_content    medi/data/food_content.csv \
#     --food_content_type Nutrient \
#     --metadata        metadata_merged.csv \
#     --output_dir      maaslin3_output/ \
#     --fixed_effects   age,gender,excl_brstfed,keep_animals \
#     --random_effects  lab_id
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(reformulas)
  library(lme4)
  library(lmerTest)
  library(optimx)
  library(maaslin3)
})

# =============================================================================
# 1. Build bacterial abundance table from MetaPhlAn *_wreads.txt profiles
#    - Extracts species-level rows (last clade token starts with "s__")
#    - Uses relative_abundance column, divided by 100 to give 0-1 fractions
#    - Sample ID derived from filename: "S001_metaphlan_profile_wreads.txt" -> "S001"
# =============================================================================
# Shared MetaPhlAn file parser - returns a long data frame with columns:
# sample_id, taxon, abundance, rank ("species" or "genus")
.parse_metaphlan_files <- function(metaphlan_dir) {
  files <- Sys.glob(file.path(metaphlan_dir, "*_metaphlan_profile_wreads.txt"))
  if (length(files) == 0)
    files <- Sys.glob(file.path(metaphlan_dir, "*_metaphlan_profile.txt"))
  if (length(files) == 0)
    stop("No MetaPhlAn files found in: ", metaphlan_dir,
         "\nExpected: *_metaphlan_profile_wreads.txt or *_metaphlan_profile.txt")

  parse_one <- function(f) {
    lines      <- readLines(f, warn = FALSE)
    data_start <- which(grepl("^k__", lines))[1]
    if (is.na(data_start)) { warning("No k__ data in: ", basename(f)); return(NULL) }

    df <- read.table(text = lines[data_start:length(lines)],
                     sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                     quote = "", fill = TRUE, comment.char = "")
    base_cols <- c("clade_name", "ncbi_tax_id", "relative_abundance", "coverage")
    colnames(df)[seq_len(min(ncol(df), length(base_cols)))] <-
      base_cols[seq_len(min(ncol(df), length(base_cols)))]
    if (!"clade_name" %in% colnames(df) || !"relative_abundance" %in% colnames(df))
      return(NULL)

    sample_id <- basename(f) %>%
      str_remove("_metaphlan_profile_wreads\\.txt$") %>%
      str_remove("_metaphlan_profile\\.txt$")

    # Rows containing |s__ have species-level assignment.
    # Extract both the s__ and g__ tokens from those rows so we can aggregate
    # to either level without re-reading the files.
    sp_rows <- df %>%
      filter(grepl("\\|s__[^|]+", clade_name)) %>%
      mutate(
        kingdom_raw  = sub("\\|.*", "", clade_name),
        is_bacteria  = grepl("Bacteria", kingdom_raw, ignore.case = TRUE),
        abundance    = suppressWarnings(as.numeric(relative_abundance)) / 100,
        species_name = sub(".*\\|s__([^|]+).*", "\\1", clade_name),
        genus_name   = sub(".*\\|g__([^|]+)\\|.*", "\\1",
                           sub(".*\\|g__([^|]+)$",   "\\1", clade_name))
      ) %>%
      filter(is_bacteria, !is.na(abundance), abundance > 0)

    if (nrow(sp_rows) == 0) return(NULL)

    # Species-level: aggregate SGBs sharing the same s__ token
    sp <- sp_rows %>%
      group_by(taxon = make.names(species_name)) %>%
      summarise(abundance = sum(abundance), .groups = "drop") %>%
      mutate(sample_id = sample_id, rank = "species")

    # Genus-level: aggregate all s__ rows sharing the same g__ token
    gn <- sp_rows %>%
      group_by(taxon = make.names(genus_name)) %>%
      summarise(abundance = sum(abundance), .groups = "drop") %>%
      mutate(sample_id = sample_id, rank = "genus")

    bind_rows(sp, gn)
  }

  bind_rows(lapply(files, parse_one))
}

# Build wide abundance table from parsed MetaPhlAn data for a given rank
.metaphlan_wide <- function(long_data, rank, label) {
  dat <- long_data %>% filter(rank == !!rank)
  if (nrow(dat) == 0)
    stop("No ", rank, "-level rows found in MetaPhlAn data.")
  wide <- dat %>%
    group_by(sample_id, taxon) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    pivot_wider(names_from = taxon, values_from = abundance, values_fill = 0) %>%
    column_to_rownames("sample_id")
  mat <- as.matrix(wide); storage.mode(mat) <- "double"
  wide <- as.data.frame(mat)
  cat(sprintf("  MetaPhlAn %s: %d samples, %d %s.\n",
              label, nrow(wide), ncol(wide), rank))
  wide
}

build_metaphlan_table <- function(metaphlan_dir, rank = "species") {
  long <- .parse_metaphlan_files(metaphlan_dir)
  if (nrow(long) == 0)
    stop("No bacterial rows extracted from MetaPhlAn files in: ", metaphlan_dir)
  .metaphlan_wide(long, rank, "bacteria")
}

# =============================================================================
# 2. Build archaeal abundance table from Bracken *_archaea_bracken_abundance.c01.txt
#    - Already species-level (taxonomy_lvl == "S" expected but not enforced)
#    - Uses fraction_total_reads column directly
#    - Sample ID: "S001_archaea_bracken_abundance.c01.txt" -> "S001"
# =============================================================================
build_archaea_table <- function(archaea_dir) {
  files <- Sys.glob(file.path(archaea_dir, "*_archaea_bracken_abundance.c01.txt"))
  if (length(files) == 0)
    stop("No archaea Bracken files found in: ", archaea_dir,
         "\nExpected pattern: *_archaea_bracken_abundance.c01.txt")

  parse_one <- function(f) {
    df <- read_tsv(f, show_col_types = FALSE)
    required <- c("name", "fraction_total_reads")
    missing  <- setdiff(required, colnames(df))
    if (length(missing) > 0) {
      warning("Skipping ", basename(f), ": missing column(s) ",
              paste(missing, collapse = ", "))
      return(NULL)
    }
    sample_id <- basename(f) %>%
      str_remove("_archaea_bracken_abundance\\.c01\\.txt$")
    df %>%
      transmute(
        sample_id = sample_id,
        taxon     = make.names(name),
        abundance = as.numeric(fraction_total_reads)
      )
  }

  all_data <- bind_rows(lapply(files, parse_one))
  if (nrow(all_data) == 0)
    stop("No archaea rows extracted from Bracken files in: ", archaea_dir)

  wide <- all_data %>%
    group_by(sample_id, taxon) %>%
    summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = taxon, values_from = abundance, values_fill = 0) %>%
    column_to_rownames("sample_id")

  mat <- as.matrix(wide); storage.mode(mat) <- "double"
  wide <- as.data.frame(mat)

  cat(sprintf("  Bracken archaea: %d samples, %d species.\n",
              nrow(wide), ncol(wide)))
  wide
}

# =============================================================================
# 3. Build HUMAnN feature tables (pathways, gene families, ECs, reactions, KOs)
#
# File layout (per-sample sub-directories, CPM-normalised, 'named' where given):
#   Pathways    : humann_dir/SXXX/out.SXXX_combined_pathabundance_normalised.tsv
#                 (or ..._pathabundance-cpm.tsv). The raw ..._pathabundance.tsv is
#                 deliberately NOT used because HUMAnN tables are analysed with
#                 normalisation NONE and so must already be normalised.
#   Gene fams   : humann_dir/SXXX/out.SXXX_combined_genefamilies-cpm.tsv
#   ECs         : humann_dir/SXXX/out.SXXX_combined_ec-cpm-named.tsv
#   Reactions   : humann_dir/SXXX/out.SXXX_combined_rxn-cpm-named.tsv
#   KO groups   : humann_dir/SXXX/out.SXXX_combined_ko-cpm-named.tsv
#
# All files:
#   - First line is a header ('# Gene Family<TAB>out.SXXX_..._Abundance-CPM')
#   - UNMAPPED / UNINTEGRATED / UNGROUPED rows are dropped (also when stratified)
#   - Community-total rows (no '|')          -> "total" tables  (get_*())
#   - Species-stratified rows 'FEAT|g__X.s__Y' -> "stratified" tables
#       (get_*_strat()); columns are named  <FEATURE>__IN__<Species>
#       e.g.  PWY.7219__IN__Escherichia_coli
#     '|unclassified' strata are dropped -- only species-assigned rows kept.
# =============================================================================

.HUMANN_DROP    <- c("UNMAPPED", "UNINTEGRATED", "UNGROUPED")
.STRAT_SEP      <- "__IN__"

# Candidate file-name patterns per layer, tried in order (first with hits wins).
.HUMANN_PATTERNS <- list(
  pathways    = c("out.*_combined_pathabundance_normalised.tsv",
                  "out.*_combined_pathabundance-cpm.tsv"),
  genefamilies= c("out.*_combined_genefamilies-cpm.tsv"),
  ec          = c("out.*_combined_ec-cpm-named.tsv"),
  reactions   = c("out.*_combined_rxn-cpm-named.tsv"),
  ko          = c("out.*_combined_ko-cpm-named.tsv")
)

# Locate files for one layer (per-sample sub-directories; also tolerates a flat
# layout with the same file names directly in humann_dir).
.humann_files <- function(humann_dir, layer) {
  for (pat in .HUMANN_PATTERNS[[layer]]) {
    f <- unique(c(Sys.glob(file.path(humann_dir, "*", pat)),
                  Sys.glob(file.path(humann_dir, pat))))
    if (length(f) > 0) return(f)
  }
  character()
}

# Sample ID from "out.S001_combined_....tsv" -> "S001"
.humann_sample_id <- function(f) {
  bn <- basename(f)
  m  <- regmatches(bn, regexpr("(?<=^out\\.)(S\\d+)", bn, perl = TRUE))
  if (length(m) == 1 && nchar(m) > 0) return(m)
  sub("_.*$", "", bn)
}

# Parse one HUMAnN TSV into a long data frame:
#   sample_id, func (raw feature name), stratum (NA for totals), abundance
# keep_stratified = FALSE skips stratified rows entirely (used for gene families
# where stratified rows would be enormous).
.parse_humann_file <- function(f, sample_id, keep_stratified = TRUE) {
  df <- tryCatch(
    read_tsv(f, skip = 1, col_names = c("feature", "abundance"),
             col_types = cols(feature = col_character(), abundance = col_double(),
                              .default = col_skip()),
             quote = "", progress = FALSE, show_col_types = FALSE),
    error = function(e) { warning("Cannot parse ", f, ": ", e$message); NULL })
  if (is.null(df) || nrow(df) == 0) return(NULL)

  df <- df[!is.na(df$abundance) & !is.na(df$feature) & nzchar(df$feature) &
             !grepl("^#", df$feature), , drop = FALSE]

  is_strat <- grepl("|", df$feature, fixed = TRUE)
  if (!keep_stratified) df <- df[!is_strat, , drop = FALSE]
  if (nrow(df) == 0) return(NULL)

  parts   <- strsplit(df$feature, "|", fixed = TRUE)
  func    <- vapply(parts, `[`, character(1), 1)
  stratum <- vapply(parts, function(x) if (length(x) > 1) x[2] else NA_character_,
                    character(1))

  out <- data.frame(sample_id = sample_id, func = func, stratum = stratum,
                    abundance = df$abundance, stringsAsFactors = FALSE)
  out[!out$func %in% .HUMANN_DROP, , drop = FALSE]
}

# Read all files of one layer -> long data frame (cached by the caller).
.humann_long <- function(humann_dir, layer, keep_stratified = TRUE) {
  files <- .humann_files(humann_dir, layer)
  if (length(files) == 0)
    stop("No HUMAnN '", layer, "' files found under ", humann_dir,
         " (tried: ", paste(.HUMANN_PATTERNS[[layer]], collapse = ", "), ").")
  cat(sprintf("  Reading %d HUMAnN %s file(s) ...\n", length(files), layer))
  long <- bind_rows(lapply(files, function(f)
    .parse_humann_file(f, .humann_sample_id(f), keep_stratified)))
  if (is.null(long) || nrow(long) == 0)
    stop("No usable rows extracted for HUMAnN ", layer, ".")

  # Some HUMAnN outputs (e.g. '*_normalised.tsv' from humann_renorm_table
  # --units relab) hold fractions, not CPM. The abundance filters, and the
  # CPM tables for the other layers, use CPM scale, so convert fractions to CPM.
  tot_sum <- tapply(long$abundance[is.na(long$stratum)],
                    long$sample_id[is.na(long$stratum)], sum)
  if (length(tot_sum) > 0 && max(tot_sum, na.rm = TRUE) <= 1.0001) {
    long$abundance <- long$abundance * 1e6
    cat(sprintf("  HUMAnN %s values are relative abundance (max per-sample total %.3f): ",
                layer, max(tot_sum)), "rescaled x1e6 to CPM scale.\n", sep = "")
  }
  long
}

# Species label from a stratum string such as "g__Escherichia.s__Escherichia_coli"
.stratum_species <- function(stratum) {
  sp <- ifelse(grepl("s__", stratum, fixed = TRUE),
               sub("^.*s__", "", stratum),
               sub("^g__", "", stratum))
  make.names(sp)
}

.long_to_wide <- function(df, label) {
  wide <- df %>%
    group_by(sample_id, feature) %>%
    summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = feature, values_from = abundance, values_fill = 0) %>%
    column_to_rownames("sample_id")
  mat <- as.matrix(wide); storage.mode(mat) <- "double"
  wide <- as.data.frame(mat)
  cat(sprintf("  HUMAnN %s: %d samples, %d features.\n", label, nrow(wide), ncol(wide)))
  wide
}

# Community-total table (unstratified rows only).
humann_total_table <- function(long, label) {
  tot <- long[is.na(long$stratum), , drop = FALSE]
  if (nrow(tot) == 0) stop("No unstratified rows for HUMAnN ", label, ".")
  # Unique, syntactic names; distinct raw names never get merged by make.names().
  fn  <- unique(tot$func)
  map <- setNames(make.names(fn, unique = TRUE), fn)
  tot$feature <- unname(map[tot$func])
  .long_to_wide(tot[, c("sample_id", "feature", "abundance")], label)
}

# Species-assigned table: one column per (function, species) pair.
humann_strat_table <- function(long, label) {
  st <- long[!is.na(long$stratum) & !grepl("unclassified", long$stratum,
                                           ignore.case = TRUE), , drop = FALSE]
  if (nrow(st) == 0) stop("No species-stratified rows for HUMAnN ", label, ".")
  fn  <- unique(st$func)
  fmap <- setNames(make.names(fn, unique = TRUE), fn)
  st$feature <- paste0(unname(fmap[st$func]), .STRAT_SEP, .stratum_species(st$stratum))
  .long_to_wide(st[, c("sample_id", "feature", "abundance")],
                paste0(label, " (species-stratified)"))
}

# =============================================================================
# 4. Build food abundance table from MEDI food_abundance.csv
# =============================================================================
build_food_abundance_table <- function(food_abundance_file,
                                       level = c("taxon", "food_group", "food")) {
  level <- match.arg(level)
  fa <- read_csv(food_abundance_file, show_col_types = FALSE,
                 col_select = any_of(c("sample_id", "relative_abundance", "taxon",
                                       "wikipedia_id", "food_group")))
  need <- c("sample_id", "relative_abundance",
            switch(level, taxon = "taxon", food_group = c("food_group", "taxon"),
                          food = "wikipedia_id"))
  miss <- setdiff(need, colnames(fa))
  if (length(miss) > 0)
    stop("food_abundance file missing column(s): ", paste(miss, collapse = ", "))

  fa <- fa[!is.na(fa$sample_id) & !is.na(fa$relative_abundance), , drop = FALSE]

  # MEDI lists one row per FOOD, so different foods that map to the same
  # organism (e.g. Roselle and Hibiscus tea -> Hibiscus sabdariffa) repeat the
  # SAME reads. 'taxon' / 'food_group' levels de-duplicate these so one organism
  # is counted once; 'food' keeps the old one-column-per-food behaviour.
  if (level == "taxon") {
    long <- fa %>% filter(!is.na(taxon)) %>%
      group_by(sample_id, feature = taxon) %>%
      summarise(value = max(relative_abundance), .groups = "drop")
  } else if (level == "food_group") {
    long <- fa %>% filter(!is.na(food_group), !is.na(taxon)) %>%
      group_by(sample_id, food_group, taxon) %>%
      summarise(value = max(relative_abundance), .groups = "drop") %>%
      group_by(sample_id, feature = food_group) %>%
      summarise(value = sum(value), .groups = "drop")
  } else {
    long <- fa %>% filter(!is.na(wikipedia_id)) %>%
      group_by(sample_id, feature = wikipedia_id) %>%
      summarise(value = sum(relative_abundance), .groups = "drop")
  }

  fn  <- unique(long$feature)
  map <- setNames(make.names(fn, unique = TRUE), fn)
  long$feature <- unname(map[long$feature])

  wide <- long %>%
    pivot_wider(names_from = feature, values_from = value, values_fill = 0) %>%
    column_to_rownames("sample_id")
  mat <- as.matrix(wide); storage.mode(mat) <- "double"
  wide <- as.data.frame(mat)
  wide <- keep_study_samples(wide, "food abundance")
  cat(sprintf("  Food abundance (%s level): %d samples, %d features.\n",
              level, nrow(wide), ncol(wide)))
  wide
}

# MEDI outputs mix samples from several datasets. Only IDs of the form S<number>
# (e.g. S001, S123) belong to this study; everything else (e.g.
# '1808-run7_E010BM4') is dropped from every MEDI-derived table.
.STUDY_SAMPLE_PATTERN <- "^S[0-9]+$"

keep_study_samples <- function(tab, label = "MEDI table") {
  ids <- trimws(rownames(tab))
  ok  <- grepl(.STUDY_SAMPLE_PATTERN, ids)
  if (!any(ok))
    stop(label, ": no sample IDs match ", .STUDY_SAMPLE_PATTERN,
         ". First IDs: ", paste(head(rownames(tab), 5), collapse = ", "))
  if (any(!ok))
    cat(sprintf("  %s: dropped %d/%d non-study sample(s) (not S<number>), e.g. %s.\n",
                label, sum(!ok), length(ok),
                paste(head(rownames(tab)[!ok], 3), collapse = ", ")))
  tab <- tab[ok, , drop = FALSE]
  rownames(tab) <- ids[ok]
  tab
}

# =============================================================================
# 4. Build food content table from MEDI food_content.csv
# =============================================================================
build_food_content_table <- function(food_content_file, source_type = NULL) {
  fc <- read_csv(food_content_file, show_col_types = FALSE)

  if (!is.null(source_type) && "source_type" %in% colnames(fc)) {
    fc <- fc %>% filter(source_type %in% !!source_type)
    if (nrow(fc) == 0)
      stop("No rows remain after source_type filter: ", paste(source_type, collapse = ", "))
  }

  # Identify the compound name column - MEDI uses 'name'
  name_col <- intersect(c("name", "compound_name", "compound"), colnames(fc))[1]
  if (is.na(name_col))
    stop("food_content file: cannot find a compound name column ",
         "(looked for: name, compound_name, compound). ",
         "Available columns: ", paste(colnames(fc), collapse = ", "))

  required <- c("sample_id", "abundance")
  missing  <- setdiff(required, colnames(fc))
  if (length(missing) > 0)
    stop("food_content file missing column(s): ", paste(missing, collapse = ", "))

  cat(sprintf("  Using '%s' as nutrient label column.\n", name_col))

  wide <- fc %>%
    mutate(.feature = make.names(.data[[name_col]])) %>%
    group_by(sample_id, .feature) %>%
    summarise(value = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = .feature, values_from = value, values_fill = 0) %>%
    column_to_rownames("sample_id")

  mat <- as.matrix(wide); storage.mode(mat) <- "double"
  wide <- as.data.frame(mat)
  wide <- keep_study_samples(wide, "food content")
  cat(sprintf("  Food content: %d samples, %d nutrients/compounds.\n",
              nrow(wide), ncol(wide)))
  wide
}

# =============================================================================
# 4b. Exclude samples flagged in the WHO anthropometric survey table
#     Any child with a non-empty 'Flag' value (e.g. HAZ, WAZ, WHZ, BAZ -- biologically
#     implausible z-scores) is removed from EVERY analysis. The exclusion list is
#     stored in .SAMPLE_FILTER and applied in filter_features() (so prevalence /
#     abundance filters are computed without these samples) and load_metadata().
# =============================================================================
.SAMPLE_FILTER <- new.env()
.SAMPLE_FILTER$exclude <- character()

drop_excluded_rows <- function(x) {
  ex <- .SAMPLE_FILTER$exclude
  if (length(ex) == 0) return(x)
  x[!rownames(x) %in% ex, , drop = FALSE]
}

build_excluded_samples <- function(anthro_file, metadata_file,
                                   anthro_id_col    = "lab_id",
                                   metadata_id_col  = "lab_id",
                                   flag_col         = "Flag",
                                   output_dir       = NULL) {
  an <- read_csv(anthro_file, col_types = cols(.default = col_character()),
                 show_col_types = FALSE)
  for (cc in c(anthro_id_col, flag_col))
    if (!cc %in% colnames(an))
      stop("Anthropometry file lacks column '", cc, "'. Available: ",
           paste(colnames(an), collapse = ", "))

  flag <- trimws(an[[flag_col]])
  flagged <- an[!is.na(flag) & nzchar(flag), , drop = FALSE]
  flagged_ids <- unique(trimws(flagged[[anthro_id_col]]))
  flagged_ids <- flagged_ids[!is.na(flagged_ids) & nzchar(flagged_ids)]
  cat(sprintf("  Anthropometry file: %d rows, %d flagged (%d unique %s).\n",
              nrow(an), nrow(flagged), length(flagged_ids), anthro_id_col))
  if (nrow(flagged) > 0) {
    tb <- table(trimws(flagged[[flag_col]]))
    cat("  Flag breakdown: ",
        paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = ", "), "\n", sep = "")
  }

  sep  <- if (grepl("\\.csv$", metadata_file, ignore.case = TRUE)) "," else "\t"
  meta <- as.data.frame(read_delim(metadata_file, delim = sep, show_col_types = FALSE,
                                   col_types = cols(.default = col_character())))
  if (!metadata_id_col %in% colnames(meta))
    stop("Metadata lacks column '", metadata_id_col, "' needed to link to the ",
         "anthropometry file (use --metadata_id_col). Available: ",
         paste(head(colnames(meta), 40), collapse = ", "))

  sample_ids <- meta[[1]]
  hit <- trimws(meta[[metadata_id_col]]) %in% flagged_ids
  excl <- unique(sample_ids[hit])
  cat(sprintf("  Matched %d/%d flagged IDs to metadata -> %d sample(s) excluded.\n",
              length(intersect(flagged_ids, trimws(meta[[metadata_id_col]]))),
              length(flagged_ids), length(excl)))
  unmatched <- setdiff(flagged_ids, trimws(meta[[metadata_id_col]]))
  if (length(unmatched) > 0)
    cat(sprintf("  Note: %d flagged ID(s) not present in metadata (nothing to exclude).\n",
                length(unmatched)))

  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    write_tsv(data.frame(sample_id = excl,
                         lab_id = trimws(meta[[metadata_id_col]])[hit]),
              file.path(output_dir, "excluded_flagged_samples.tsv"))
  }
  excl
}

# =============================================================================
# 5. Filter features by prevalence and mean value
# =============================================================================
filter_features <- function(tab, min_prevalence = 0.10, min_abundance = 0.01,
                             label = "features") {
  tab <- drop_excluded_rows(tab)   # flagged samples never contribute to filters
  prev_mask  <- colMeans(tab > 0) >= min_prevalence
  abund_mask <- colMeans(tab)     >= min_abundance
  filtered   <- tab[, prev_mask & abund_mask, drop = FALSE]
  cat(sprintf("  %s: %d retained (prevalence >= %.2f, mean >= %.4f).\n",
              label, ncol(filtered), min_prevalence, min_abundance))
  filtered
}

# =============================================================================
# 5b. Load and merge alpha diversity statistics
#     Source: diversity_analysis_outputs/diversity_stats_per_sample.csv
#     Columns used: Observed, shannon, inverse_simpson, chao1
#     Sample IDs in that file have a "_metaphlan_wreads" suffix which is
#     stripped before merging so they align with the main metadata IDs.
# =============================================================================
load_diversity <- function(diversity_file,
                            metrics = c("Observed", "shannon",
                                        "inverse_simpson", "chao1")) {
  div <- read_csv(diversity_file, show_col_types = FALSE)

  # First column is the sample ID (may be unnamed, shown as "...1" or "")
  id_col <- colnames(div)[1]
  div <- as.data.frame(div)
  rownames(div) <- sub("_metaphlan_wreads$", "", div[[id_col]])

  # Keep only the requested alpha diversity columns
  keep <- intersect(metrics, colnames(div))
  missing <- setdiff(metrics, colnames(div))
  if (length(missing) > 0)
    warning("Diversity metric(s) not found in file: ",
            paste(missing, collapse = ", "))
  if (length(keep) == 0)
    stop("No requested diversity metrics found in: ", diversity_file)

  div <- div[, keep, drop = FALSE]

  # Ensure numeric
  for (col in colnames(div))
    div[[col]] <- suppressWarnings(as.numeric(div[[col]]))

  cat(sprintf("  Diversity file: %d samples, metrics: %s.
",
              nrow(div), paste(keep, collapse = ", ")))
  div
}

# =============================================================================
# 6. Load and validate metadata
# =============================================================================
load_metadata <- function(metadata_file, fixed_effects = NULL,
                           random_effects = NULL) {
  sep  <- if (grepl("\\.csv$", metadata_file, ignore.case = TRUE)) "," else "\t"
  meta <- read_delim(metadata_file, delim = sep, show_col_types = FALSE)
  id_col <- colnames(meta)[1]
  meta   <- as.data.frame(meta)
  rownames(meta) <- meta[[id_col]]

  n_before <- nrow(meta)
  meta <- drop_excluded_rows(meta)
  if (nrow(meta) < n_before)
    cat(sprintf("  Excluded %d flagged sample(s) from metadata.\n", n_before - nrow(meta)))

  cat(sprintf("  Metadata: %d samples, %d columns.\n", nrow(meta), ncol(meta)))

  present <- intersect(c(fixed_effects, random_effects), colnames(meta))
  for (col in present) {
    cat(sprintf("    %-25s  %d unique, %d NAs\n",
                col, length(unique(na.omit(meta[[col]]))), sum(is.na(meta[[col]]))))
  }
  missing <- setdiff(c(fixed_effects, random_effects), colnames(meta))
  if (length(missing) > 0)
    warning("Effect column(s) not in metadata: ", paste(missing, collapse = ", "))
  meta
}

# =============================================================================
# 7. Preprocess metadata (numeric coercion + factor conversion)
# =============================================================================
preprocess_metadata <- function(meta, fixed_effects = NULL, random_effects = NULL) {
  numeric_cols <- c("age", "age_days", "age_months", "muac", "muac_z",
                    "wfh_zscore", "hfa_zscore", "wfa_zscore",
                    "height", "weight", "hemocue_result", "temperature")
  factor_cols  <- c("gender", "age_group", "diarrhoea", "excl_brstfed",
                    "keep_animals", "lab_id", "cluster", "village",
                    "malaria_prev", "hiv_test")

  for (col in intersect(numeric_cols, colnames(meta))) {
    orig <- meta[[col]]
    conv <- suppressWarnings(as.numeric(as.character(orig)))
    n    <- sum(!is.na(orig) & is.na(conv))
    if (n > 0) warning(sprintf("%d non-numeric value(s) in '%s' coerced to NA.", n, col))
    meta[[col]] <- conv
  }
  for (col in intersect(factor_cols, colnames(meta)))
    meta[[col]] <- as.factor(meta[[col]])

  all_effects <- unique(c(fixed_effects, random_effects))
  na_counts   <- colSums(is.na(meta[, intersect(all_effects, colnames(meta)),
                                     drop = FALSE]))
  if (any(na_counts > 0))
    warning("NAs in effect columns: ",
            paste(sprintf("%s(%d)", names(na_counts[na_counts > 0]),
                          na_counts[na_counts > 0]), collapse = ", "))
  meta
}

# =============================================================================
# 8. Align samples between a feature table and metadata
# =============================================================================
align_samples <- function(tab, meta, label = "feature table") {
  common <- intersect(rownames(tab), rownames(meta))
  if (length(common) == 0)
    stop("No common samples between ", label, " and metadata.\n",
         "  Feature table IDs (first 5): ",
         paste(head(rownames(tab), 5), collapse = ", "), "\n",
         "  Metadata IDs (first 5): ",
         paste(head(rownames(meta), 5), collapse = ", "))
  if (length(common) < nrow(tab))
    cat(sprintf("  Note: %d %s samples not in metadata (dropped).\n",
                nrow(tab) - length(common), label))
  if (length(common) < nrow(meta))
    cat(sprintf("  Note: %d metadata samples not in %s (dropped).\n",
                nrow(meta) - length(common), label))
  list(tab = tab[common, , drop = FALSE], meta = meta[common, , drop = FALSE])
}

# =============================================================================
# 9. Core Maaslin3 wrapper
# =============================================================================
`%||%` <- function(a, b) if (!is.null(a)) a else b

run_maaslin3 <- function(input_tab, meta, output_dir,
                          fixed_effects  = character(),
                          random_effects = NULL,
                          normalization  = "TSS",
                          transform      = "LOG",
                          cores          = 4) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  fe <- intersect(fixed_effects,                   colnames(meta))
  re <- intersect(random_effects %||% character(), colnames(meta))

  if (length(fe) == 0)
    stop("None of the requested fixed effects found in metadata.")

  # Subset metadata to ONLY the columns Maaslin3 will use.
  # Passing a wide data frame (e.g. 351 columns) causes Maaslin3 to process
  # all columns and crash on NAs in unrelated fields.
  meta <- meta[, unique(c(fe, re)), drop = FALSE]

  # Final NA drop restricted to effect columns only.
  complete_mask <- complete.cases(meta)
  if (!all(complete_mask)) {
    n_drop    <- sum(!complete_mask)
    keep      <- rownames(meta)[complete_mask]
    meta      <- meta[keep,      , drop = FALSE]
    input_tab <- input_tab[keep, , drop = FALSE]
    cat(sprintf("  Dropped %d sample(s) with NA in effect columns. Remaining: %d\n",
                n_drop, nrow(meta)))
  }

  cat(sprintf("  Fixed : %s\n", paste(fe, collapse = ", ")))
  if (length(re) > 0)
    cat(sprintf("  Random: %s\n", paste(re, collapse = ", ")))
  cat(sprintf("  Norm: %s | Transform: %s | Samples: %d | Features: %d\n",
              normalization, transform, nrow(input_tab), ncol(input_tab)))

  if (ncol(input_tab) == 0 || nrow(input_tab) < 10) {
    cat(sprintf("  SKIPPING: only %d feature(s) / %d sample(s) after filtering.\n",
                ncol(input_tab), nrow(input_tab)))
    return(invisible(NULL))
  }

  fit <- maaslin3(
    input_data     = t(input_tab),
    input_metadata = meta,
    output         = output_dir,
    fixed_effects  = fe,
    random_effects = if (length(re) > 0) re else NULL,
    normalization  = normalization,
    transform      = transform,
    min_abundance  = 0,
    min_prevalence = 0,
    cores          = cores
  )

  if (!is.null(fit$results) && nrow(fit$results) > 0) {
    sig <- fit$results %>% filter(qval_joint < 0.25) %>% arrange(qval_joint)
    write_tsv(fit$results, file.path(output_dir, "all_results.tsv"))
    write_tsv(sig,         file.path(output_dir, "significant_results.tsv"))
    cat(sprintf("  Significant (q < 0.25): %d\n", nrow(sig)))
  } else {
    cat("  No results returned.\n")
  }
  invisible(fit)
}

# =============================================================================
# 10. Helpers for univariate diet-vs-taxa pooled results
# =============================================================================
.pval_col <- function(df) {
  for (col in c("pval_joint", "pval_individual", "pval")) {
    if (col %in% colnames(df)) return(col)
  }
  stop("No p-value column found. Available: ", paste(colnames(df), collapse = ", "))
}

.extract_diet_rows <- function(results_df) {
  if (is.null(results_df) || nrow(results_df) == 0) return(NULL)
  for (col in c("metadata", "value", "variable", "term")) {
    if (col %in% colnames(results_df) && "diet_value" %in% results_df[[col]])
      return(results_df[results_df[[col]] == "diet_value", , drop = FALSE])
  }
  results_df
}

# =============================================================================
# 11. Diet vs taxa: univariate (one diet feature at a time) pipeline
# =============================================================================
run_diet_taxa_pipeline <- function(taxa_tab, diet_tab, output_dir,
                                    covariates_meta  = NULL,
                                    fixed_covariates = NULL,
                                    random_effects   = NULL,
                                    normalization    = "TSS",
                                    transform        = "LOG",
                                    cores            = 4) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  common <- Reduce(intersect, Filter(Negate(is.null),
    list(rownames(taxa_tab), rownames(diet_tab),
         if (!is.null(covariates_meta)) rownames(covariates_meta) else NULL)))

  if (length(common) == 0)
    stop("No common samples between taxa table, diet table",
         if (!is.null(covariates_meta)) ", and covariates" else "", ".")

  taxa_tab <- taxa_tab[common, , drop = FALSE]
  diet_tab <- diet_tab[common, , drop = FALSE]

  cat(sprintf("  %d common samples | %d taxa | %d diet features.\n",
              length(common), ncol(taxa_tab), ncol(diet_tab)))

  results_list <- list()

  for (feat in colnames(diet_tab)) {
    feat_dir <- file.path(output_dir, "per_feature", make.names(feat))
    meta <- data.frame(diet_value = diet_tab[[feat]],
                       row.names  = common, check.names = FALSE)
    fe <- "diet_value"
    re <- NULL

    if (!is.null(covariates_meta) && !is.null(fixed_covariates)) {
      fc <- intersect(fixed_covariates, colnames(covariates_meta))
      if (length(fc) > 0) {
        meta <- cbind(meta, covariates_meta[common, fc, drop = FALSE])
        fe   <- c(fe, fc)
      }
    }
    if (!is.null(covariates_meta) && !is.null(random_effects)) {
      rc <- intersect(random_effects, colnames(covariates_meta))
      if (length(rc) > 0) {
        meta <- cbind(meta, covariates_meta[common, rc, drop = FALSE])
        re   <- rc
      }
    }

    # Drop samples with NA in any fixed effect column before fitting
    complete <- complete.cases(meta)
    if (!all(complete)) {
      n_drop   <- sum(!complete)
      meta     <- meta[complete,                    , drop = FALSE]
      taxa_sub <- taxa_tab[rownames(meta), , drop = FALSE]
      cat(sprintf("    '%s': dropped %d sample(s) with NA in covariates.\n",
                  feat, n_drop))
    } else {
      taxa_sub <- taxa_tab[rownames(meta), , drop = FALSE]
    }

    fit <- tryCatch(
      run_maaslin3(taxa_sub, meta, feat_dir,
                   fixed_effects = fe, random_effects = re,
                   normalization = normalization, transform = transform,
                   cores = cores),
      error = function(e) {
        warning(sprintf("Maaslin3 failed for '%s': %s", feat, conditionMessage(e)))
        NULL
      }
    )

    if (!is.null(fit) && !is.null(fit$results) && nrow(fit$results) > 0) {
      rows <- .extract_diet_rows(fit$results)
      if (!is.null(rows) && nrow(rows) > 0) {
        rows$diet_feature <- feat
        results_list[[feat]] <- rows
      }
    }
  }

  if (length(results_list) == 0) {
    cat("  No diet-taxa associations returned.\n")
    return(invisible(NULL))
  }

  combined <- bind_rows(results_list)
  pcol     <- .pval_col(combined)
  combined <- combined %>%
    mutate(qval_global_BH = p.adjust(.data[[pcol]], method = "BH")) %>%
    arrange(qval_global_BH)

  sig <- combined %>% filter(qval_global_BH < 0.25)
  write_tsv(combined, file.path(output_dir, "diet_taxa_all_results.tsv"))
  write_tsv(sig,      file.path(output_dir, "diet_taxa_significant_results.tsv"))
  cat(sprintf("  Total associations: %d | Significant (global BH q < 0.25): %d\n",
              nrow(combined), nrow(sig)))
  invisible(combined)
}

# =============================================================================
# 12. Diet vs metadata pipeline
# =============================================================================
run_diet_metadata_pipeline <- function(diet_tab, metadata_file, output_dir,
                                        fixed_effects  = character(),
                                        random_effects = NULL,
                                        normalization  = "TSS",
                                        transform      = "LOG",
                                        cores          = 4) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  meta    <- load_metadata(metadata_file, fixed_effects, random_effects)
  meta    <- preprocess_metadata(meta, fixed_effects, random_effects)
  aligned <- align_samples(diet_tab, meta, "diet table")
  diet_tab <- aligned$tab
  meta     <- aligned$meta

  fe <- intersect(fixed_effects, colnames(meta))
  re <- intersect(random_effects %||% character(), colnames(meta))

  # Drop samples with NA in any fixed effect column
  if (length(fe) > 0) {
    complete  <- complete.cases(meta[, fe, drop = FALSE])
    n_dropped <- sum(!complete)
    meta     <- meta[complete,     , drop = FALSE]
    diet_tab <- diet_tab[complete, , drop = FALSE]
    if (n_dropped > 0)
      cat(sprintf("  Dropped %d sample(s) with NA in fixed effect column(s). Remaining: %d\n",
                  n_dropped, nrow(meta)))
  }

  run_maaslin3(diet_tab, meta, output_dir,
               fixed_effects  = fe,
               random_effects = if (length(re) > 0) re else NULL,
               normalization = normalization, transform = transform, cores = cores)
}

# =============================================================================
# 13. Taxa vs metadata pipeline (bacteria or archaea)
# =============================================================================
run_taxa_metadata_pipeline <- function(taxa_tab, metadata_file, output_dir,
                                        fixed_effects  = character(),
                                        random_effects = NULL,
                                        min_prevalence = 0.10,
                                        min_abundance  = 0.01,
                                        label          = "taxa",
                                        extra_meta_tab = NULL,
                                        normalization  = "TSS",
                                        transform      = "LOG",
                                        cores          = 4) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  taxa_tab <- filter_features(taxa_tab, min_prevalence, min_abundance, label)

  meta    <- load_metadata(metadata_file, fixed_effects, random_effects)
  meta    <- preprocess_metadata(meta, fixed_effects, random_effects)

  # Merge extra predictor columns (e.g. taxa abundances) into metadata.
  # Restricts to samples present in all three tables.
  if (!is.null(extra_meta_tab)) {
    common <- Reduce(intersect, list(rownames(meta),
                                     rownames(extra_meta_tab),
                                     rownames(taxa_tab)))
    if (length(common) == 0)
      stop("No common samples between feature table, metadata, and extra_meta_tab.")
    meta     <- cbind(meta[common, , drop = FALSE],
                      extra_meta_tab[common, , drop = FALSE])
    taxa_tab <- taxa_tab[common, , drop = FALSE]
    cat(sprintf("  Merged %d extra predictor column(s) into metadata. %d common samples.\n",
                ncol(extra_meta_tab), length(common)))
  }

  aligned  <- align_samples(taxa_tab, meta, label)
  taxa_tab <- aligned$tab
  meta     <- aligned$meta

  cat(sprintf("  Final: %d samples, %d %s.\n", nrow(taxa_tab), ncol(taxa_tab), label))

  fe <- intersect(fixed_effects, colnames(meta))
  re <- intersect(random_effects %||% character(), colnames(meta))

  missing_fe <- setdiff(fixed_effects, colnames(meta))
  if (length(missing_fe) > 0)
    cat(sprintf("  Skipping absent fixed effect(s): %s\n",
                paste(missing_fe, collapse = ", ")))

  # Drop samples with NA in any fixed effect or z-score column.
  # Maaslin3 cannot handle NA values in fixed effects and will crash.
  zscore_cols    <- c("wfh_zscore", "hfa_zscore", "wfa_zscore")
  drop_cols      <- unique(c(intersect(zscore_cols, colnames(meta)),
                              intersect(fe, colnames(meta))))
  if (length(drop_cols) > 0) {
    complete  <- complete.cases(meta[, drop_cols, drop = FALSE])
    n_dropped <- sum(!complete)
    meta     <- meta[complete,     , drop = FALSE]
    taxa_tab <- taxa_tab[complete, , drop = FALSE]
    if (n_dropped > 0)
      cat(sprintf("  Dropped %d sample(s) with NA in fixed effect/z-score column(s). Remaining: %d\n",
                  n_dropped, nrow(meta)))
  }

  # Save aligned inputs
  write_tsv(as.data.frame(t(taxa_tab)) %>% rownames_to_column("taxon"),
            file.path(output_dir, "taxa_input.tsv"))
  meta_save <- meta[, colnames(meta) != "sample_id", drop = FALSE]
  write_tsv(meta_save %>% rownames_to_column("sample_id"),
            file.path(output_dir, "metadata_input.tsv"))

  run_maaslin3(taxa_tab, meta, file.path(output_dir, "maaslin3_results"),
               fixed_effects = fe, random_effects = if (length(re) > 0) re else NULL,
               normalization = normalization, transform = transform, cores = cores)
}


# =============================================================================
# 14. Univariate metadata screening pipeline
#
# Runs each metadata variable independently as a single fixed effect against
# a feature table (bacteria, pathways, KOs etc.), then pools all results and
# applies global Benjamini-Hochberg correction. This avoids the underdetermined
# model problem that arises when fitting 40+ variables simultaneously with
# ~109 samples, and gives each variable a fair, independent test.
#
# Per-variable Maaslin3 outputs go to output_dir/per_variable/VARNAME/.
# Pooled results go to output_dir/univariate_all_results.tsv and
# output_dir/univariate_significant_results.tsv (global BH q < 0.25).
# =============================================================================
run_univariate_pipeline <- function(
    taxa_tab,
    metadata_file,
    output_dir,
    fixed_effects,          # vector of variable names to test, one at a time
    random_effects  = NULL,
    min_prevalence  = 0.10,
    min_abundance   = 0.01,
    label           = "features",
    extra_meta_tab  = NULL, # merged into metadata alongside each variable
    normalization   = "TSS",
    transform       = "LOG",
    cores           = 4
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Filter feature table once up front
  taxa_tab <- filter_features(taxa_tab, min_prevalence, min_abundance, label)

  # Load metadata once
  meta_full <- load_metadata(metadata_file, fixed_effects, random_effects)
  meta_full <- preprocess_metadata(meta_full, fixed_effects, random_effects)

  # Merge extra predictors (e.g. nutrient table) if supplied
  if (!is.null(extra_meta_tab)) {
    common <- Reduce(intersect, list(rownames(meta_full),
                                     rownames(extra_meta_tab),
                                     rownames(taxa_tab)))
    if (length(common) == 0)
      stop("No common samples between feature table, metadata, and extra_meta_tab.")
    meta_full <- cbind(meta_full[common, , drop = FALSE],
                       extra_meta_tab[common, , drop = FALSE])
    taxa_tab  <- taxa_tab[common, , drop = FALSE]
    cat(sprintf("  Merged %d extra predictor column(s). %d common samples.\n",
                ncol(extra_meta_tab), length(common)))
  }

  # Align to common samples
  aligned  <- align_samples(taxa_tab, meta_full, label)
  taxa_tab <- aligned$tab
  meta_full <- aligned$meta

  cat(sprintf("  Univariate screen: %d samples | %d %s | %d variables to test.\n",
              nrow(taxa_tab), ncol(taxa_tab), label, length(fixed_effects)))

  if (ncol(taxa_tab) == 0 || nrow(taxa_tab) < 10) {
    cat("  SKIPPING univariate screen: no features or too few samples after filtering.\n")
    return(invisible(NULL))
  }

  re <- intersect(random_effects %||% character(), colnames(meta_full))

  results_list <- list()

  for (var in fixed_effects) {
    if (!var %in% colnames(meta_full)) {
      cat(sprintf("  [SKIP] '%s' not in metadata.\n", var))
      next
    }

    var_dir <- file.path(output_dir, "per_variable", make.names(var))

    # Build per-variable metadata: just this variable + random effects
    meta_var <- meta_full[, unique(c(var, re)), drop = FALSE]

    # Drop samples with NA in this variable
    complete  <- complete.cases(meta_var[, var, drop = FALSE])
    n_drop    <- sum(!complete)
    if (n_drop > 0)
      cat(sprintf("  [%s] Dropping %d NA sample(s). Remaining: %d\n",
                  var, n_drop, sum(complete)))
    meta_var  <- meta_var[complete,      , drop = FALSE]
    taxa_sub  <- taxa_tab[rownames(meta_var), , drop = FALSE]

    if (nrow(meta_var) < 10) {
      cat(sprintf("  [SKIP] '%s': fewer than 10 complete samples after NA removal.\n", var))
      next
    }

    # Skip variables with zero variance after NA removal
    col_vals <- meta_var[[var]]
    if (is.numeric(col_vals) && var(col_vals, na.rm = TRUE) == 0) {
      cat(sprintf("  [SKIP] '%s': zero variance.\n", var))
      next
    }
    if (is.factor(col_vals) && length(levels(droplevels(col_vals))) < 2) {
      cat(sprintf("  [SKIP] '%s': only one level after NA removal.\n", var))
      next
    }

    fit <- tryCatch(
      run_maaslin3(taxa_sub, meta_var, var_dir,
                   fixed_effects  = var,
                   random_effects = if (length(re) > 0) re else NULL,
                   normalization  = normalization,
                   transform      = transform,
                   cores          = cores),
      error = function(e) {
        warning(sprintf("Maaslin3 failed for '%s': %s", var, conditionMessage(e)))
        NULL
      }
    )

    if (!is.null(fit) && !is.null(fit$results) && nrow(fit$results) > 0) {
      # Extract rows for this variable
      rows <- fit$results
      for (col in c("metadata", "value", "variable", "term")) {
        if (col %in% colnames(rows) && var %in% rows[[col]]) {
          rows <- rows[rows[[col]] == var, , drop = FALSE]
          break
        }
      }
      if (nrow(rows) > 0) {
        rows$metadata_variable <- var
        results_list[[var]] <- rows
      }
    }
  }

  if (length(results_list) == 0) {
    cat("  No univariate associations returned.\n")
    return(invisible(NULL))
  }

  combined <- bind_rows(results_list)
  pcol     <- .pval_col(combined)
  combined <- combined %>%
    mutate(qval_global_BH = p.adjust(.data[[pcol]], method = "BH")) %>%
    arrange(qval_global_BH)

  sig <- combined %>% filter(qval_global_BH < 0.25)

  write_tsv(combined, file.path(output_dir, "univariate_all_results.tsv"))
  write_tsv(sig,      file.path(output_dir, "univariate_significant_results.tsv"))

  cat(sprintf("  Tests: %d | Significant (global BH q < 0.25): %d\n",
              nrow(combined), nrow(sig)))
  cat(sprintf("  Results: %s\n", file.path(output_dir, "univariate_all_results.tsv")))

  invisible(combined)
}

# =============================================================================
# 15. Master pipeline
# =============================================================================
# Full variable list for univariate screening -- one per biological concept,
# D1 enrolment timepoint only, redundant measures collapsed to one representative.
SCREENING_VARIABLES <- c(
  # Alpha diversity metrics (from diversity_stats_per_sample.csv).
  # Ensure these columns are present in your metadata_merged.csv, or merge
  # diversity_stats_per_sample.csv into it first (strip the _metaphlan_wreads
  # suffix from sample IDs before merging).
  "Observed", "shannon", "inverse_simpson", "chao1",
  # Anthropometry
  "age_months", "wfh_zscore", "hfa_zscore", "wfa_zscore", "muac", "hemocue_result",
  # Socioeconomic / household
  "education", "enough_food", "economic_step", "housing", "electricity",
  "financial_account", "drinking_water_src", "toilet_type", "treat_dwater",
  "num_of_people",
  # Animal contact
  "keep_animals", "livestock", "companion_pet",
  "chicken_count", "cattle_count", "goat_count",
  # Clinical at enrolment
  "gender", "diarrhoea", "fever", "temperature",
  "malnutrition", "anaemia", "sickle_cell",
  "malariatest_result", "malaria_episodes", "malaria_doses",
  "excl_brstfed", "hiv_exposed", "hivtest_result",
  # Vaccination
  "bcg", "opv0", "opv1", "opv2", "opv3",
  "measlesrubella1", "pcoccal",
  "rtss1", "rtss2", "rtss3", "rtss4",
  # Location / clustering
  "location", "cluster"
)

run_medi_maaslin3_pipeline <- function(
    metaphlan_dir         = NULL,
    archaea_dir           = NULL,
    humann_dir            = NULL,
    food_content_file     = NULL,
    food_abundance_file   = NULL,   # MEDI food_abundance.csv
    food_abundance_level  = "taxon", # taxon | food_group | food
    food_normalization    = "NONE", # foods are already fractions of total reads
    diversity_file        = NULL,   # diversity_stats_per_sample.csv
    metadata_file,
    output_dir             = "maaslin3_output",
    analyses               = "auto",
    run_univariate         = TRUE,    # run univariate screening for all feature tables
    screening_variables    = SCREENING_VARIABLES,
    fixed_effects          = c("age_months", "gender", "excl_brstfed", "keep_animals"),
    random_effects         = NULL,
    food_content_type      = NULL,
    min_prevalence         = 0.10,
    min_abundance          = 0.01,
    diet_min_prevalence  = 0.10,
    diet_min_abundance   = 0,
    # HUMAnN-specific filters applied before Maaslin3.
    # Gene families (UniRef90) can have 500k+ features; these defaults cut to
    # a tractable number. Raise humann_min_prevalence to 0.5+ if still slow.
    humann_min_prevalence = 0.30,
    humann_min_abundance  = 0.001,
    # Filters for species-stratified (function x species) tables. Looser than
    # the community-total filters because each species only carries part of a
    # function's abundance.
    strat_min_prevalence  = 0.10,
    strat_min_abundance   = 0.001,
    # HUMAnN tables are already CPM-normalised, so re-normalising (TSS) after
    # prevalence filtering would distort them. Used whenever a HUMAnN table is
    # the Maaslin3 response; taxa/nutrient responses keep `normalization`.
    humann_normalization  = "NONE",
    # Flagged-sample exclusion (WHO anthropometric survey table)
    anthro_file           = NULL,
    anthro_id_col         = "lab_id",   # ID column in the anthropometry CSV
    metadata_id_col       = "lab_id",   # matching column in the metadata
    anthro_flag_col       = "Flag",
    normalization        = "TSS",
    transform            = "LOG",
    cores                = 4
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  valid_analyses <- c(
    "bacteria_species_metadata", "bacteria_genus_metadata",
    "archaea_metadata",
    "pathways_metadata",         "genefamilies_metadata",
    "reactions_metadata",        "ko_metadata",        "ec_metadata",
    "pathways_strat_metadata",   "ec_strat_metadata",
    "reactions_strat_metadata",  "ko_strat_metadata",
    "nutrient_metadata",
    "nutrient_bacteria_species", "nutrient_bacteria_genus",
    "nutrient_archaea",
    # Functional vs diet
    "nutrient_vs_pathways",      "nutrient_vs_reactions",     "nutrient_vs_ko",
    "nutrient_vs_ec",
    "nutrient_vs_pathways_strat","nutrient_vs_ec_strat",
    "nutrient_vs_reactions_strat","nutrient_vs_ko_strat",
    # Functional vs taxonomy
    "pathways_vs_bacteria_species", "pathways_vs_bacteria_genus", "pathways_vs_archaea",
    "reactions_vs_bacteria_species","reactions_vs_bacteria_genus","reactions_vs_archaea",
    "ko_vs_bacteria_species",    "ko_vs_bacteria_genus",      "ko_vs_archaea",
    "ec_vs_bacteria_species",    "ec_vs_bacteria_genus",      "ec_vs_archaea",
    # MEDI food abundance (each food/taxon tested univariately)
    "food_metadata",
    "food_vs_bacteria_species",  "food_vs_bacteria_genus",    "food_vs_archaea",
    "food_vs_pathways",          "food_vs_ec",                "food_vs_reactions",
    "food_vs_ko",
    "food_vs_pathways_strat",    "food_vs_ec_strat",
    "food_vs_reactions_strat",   "food_vs_ko_strat"
  )

  if (identical(analyses, "auto")) {
    sel <- character()
    if (!is.null(metaphlan_dir)) {
      sel <- c(sel, "bacteria_species_metadata", "bacteria_genus_metadata")
      if (!is.null(food_content_file))
        sel <- c(sel, "nutrient_bacteria_species", "nutrient_bacteria_genus")
    }
    if (!is.null(archaea_dir)) {
      sel <- c(sel, "archaea_metadata")
      if (!is.null(food_content_file))
        sel <- c(sel, "nutrient_archaea")
    }
    if (!is.null(humann_dir)) {
      # genefamilies_metadata excluded from auto (UniRef90 too large).
      sel <- c(sel, "pathways_metadata", "reactions_metadata", "ko_metadata",
                    "ec_metadata", "pathways_strat_metadata")
      if (!is.null(food_content_file))
        sel <- c(sel, "nutrient_vs_pathways", "nutrient_vs_reactions", "nutrient_vs_ko",
                      "nutrient_vs_ec", "nutrient_vs_pathways_strat")
      if (!is.null(metaphlan_dir))
        sel <- c(sel, "pathways_vs_bacteria_species", "pathways_vs_bacteria_genus",
                       "reactions_vs_bacteria_species", "reactions_vs_bacteria_genus",
                       "ko_vs_bacteria_species", "ko_vs_bacteria_genus",
                       "ec_vs_bacteria_species", "ec_vs_bacteria_genus")
      if (!is.null(archaea_dir))
        sel <- c(sel, "pathways_vs_archaea", "reactions_vs_archaea", "ko_vs_archaea",
                      "ec_vs_archaea")
    }
    if (!is.null(food_content_file))
      sel <- c(sel, "nutrient_metadata")
    if (!is.null(food_abundance_file)) {
      sel <- c(sel, "food_metadata")
      if (!is.null(metaphlan_dir))
        sel <- c(sel, "food_vs_bacteria_species", "food_vs_bacteria_genus")
      if (!is.null(archaea_dir))
        sel <- c(sel, "food_vs_archaea")
      if (!is.null(humann_dir))
        sel <- c(sel, "food_vs_pathways", "food_vs_ec", "food_vs_reactions",
                      "food_vs_ko", "food_vs_pathways_strat")
    }
  } else if (identical(analyses, "all")) {
    sel <- valid_analyses
  } else {
    sel     <- trimws(strsplit(analyses, ",")[[1]])
    unknown <- setdiff(sel, valid_analyses)
    if (length(unknown) > 0)
      stop("Unknown --analyses value(s): ", paste(unknown, collapse = ", "),
           "\nValid: ", paste(valid_analyses, collapse = ", "))
  }

  cat(sprintf("\nAnalyses to run: %s\n", paste(sel, collapse = ", ")))

  # ---------------------------------------------------------------------------
  # Exclude samples flagged in the WHO anthropometric survey table.
  # Stored globally; applied in filter_features() and load_metadata().
  # ---------------------------------------------------------------------------
  .SAMPLE_FILTER$exclude <- character()
  if (!is.null(anthro_file)) {
    cat("\nExcluding flagged samples (WHO anthropometry) ...\n")
    .SAMPLE_FILTER$exclude <- build_excluded_samples(
      anthro_file, metadata_file,
      anthro_id_col = anthro_id_col, metadata_id_col = metadata_id_col,
      flag_col = anthro_flag_col, output_dir = output_dir)
  }

  # ---------------------------------------------------------------------------
  # Merge alpha diversity metrics into metadata file if supplied.
  # This is done by writing a temporary merged CSV that all downstream
  # load_metadata() calls will read, ensuring diversity columns are available
  # as screening variables without modifying the original metadata file.
  # ---------------------------------------------------------------------------
  .metadata_file_to_use <- metadata_file
  if (!is.null(diversity_file)) {
    cat("\nMerging alpha diversity metrics into metadata ...\n")
    div_tab <- load_diversity(diversity_file)

    # Read current metadata, merge diversity, write to temp file
    sep  <- if (grepl("\\.csv$", metadata_file, ignore.case = TRUE)) "," else "\t"
    meta_tmp <- read_delim(metadata_file, delim = sep, show_col_types = FALSE)
    id_col   <- colnames(meta_tmp)[1]
    meta_tmp <- as.data.frame(meta_tmp)
    rownames(meta_tmp) <- meta_tmp[[id_col]]

    common <- intersect(rownames(meta_tmp), rownames(div_tab))
    if (length(common) == 0)
      warning("No common samples between metadata and diversity file. ",
              "Check that sample IDs match after stripping _metaphlan_wreads suffix.")
    else {
      for (col in colnames(div_tab))
        meta_tmp[common, col] <- div_tab[common, col]
      cat(sprintf("  Added diversity metrics for %d/%d samples.\n",
                  length(common), nrow(meta_tmp)))
    }

    tmp_path <- file.path(output_dir, ".metadata_with_diversity.csv")
    meta_save <- meta_tmp[, colnames(meta_tmp) != id_col, drop = FALSE]
    write_csv(meta_save %>% rownames_to_column(id_col), tmp_path)
    .metadata_file_to_use <- tmp_path
    cat(sprintf("  Temporary merged metadata: %s\n", tmp_path))
  }

  # Replace metadata_file with the (possibly diversity-enriched) version
  # for all downstream calls
  metadata_file <- .metadata_file_to_use

  # ---------------------------------------------------------------------------
  # Lazy loaders -- each table built at most once
  # ---------------------------------------------------------------------------
  .bact_sp_tab   <- NULL
  .bact_gn_tab   <- NULL
  .arch_tab      <- NULL
  .pathway_tab   <- NULL
  .genefam_tab   <- NULL
  .reaction_tab  <- NULL
  .ko_tab        <- NULL
  .nutrient_tab  <- NULL
  .metaphlan_raw <- NULL   # shared parsed long-format data

  get_metaphlan_raw <- function() {
    if (is.null(.metaphlan_raw)) {
      cat("\nParsing MetaPhlAn files ...\n")
      .metaphlan_raw <<- .parse_metaphlan_files(metaphlan_dir)
    }
    .metaphlan_raw
  }
  get_bacteria_species <- function() {
    if (is.null(.bact_sp_tab)) {
      cat("\nBuilding bacterial species table (MetaPhlAn) ...\n")
      .bact_sp_tab <<- .metaphlan_wide(get_metaphlan_raw(), "species", "bacteria")
    }
    .bact_sp_tab
  }
  get_bacteria_genus <- function() {
    if (is.null(.bact_gn_tab)) {
      cat("\nBuilding bacterial genus table (MetaPhlAn) ...\n")
      .bact_gn_tab <<- .metaphlan_wide(get_metaphlan_raw(), "genus", "bacteria")
    }
    .bact_gn_tab
  }
  get_archaea <- function() {
    if (is.null(.arch_tab)) {
      cat("\nBuilding archaeal taxa table (Bracken) ...\n")
      .arch_tab <<- build_archaea_table(archaea_dir)
    }
    .arch_tab
  }
  # ---- HUMAnN getters -------------------------------------------------------
  .tab_cache <- list()
  .cached <- function(key, builder) {
    if (is.null(.tab_cache[[key]])) .tab_cache[[key]] <<- builder()
    .tab_cache[[key]]
  }
  .hum_long_cache <- list()
  .layer_label <- c(pathways     = "pathway abundances",
                    genefamilies = "gene families",
                    ec           = "EC abundances",
                    reactions    = "reactions",
                    ko           = "KO abundances")

  # Long-format parse of one layer, shared by its total and stratified tables.
  # Flagged samples are dropped immediately to save memory.
  get_humann_long <- function(layer) {
    if (layer == "genefamilies") {           # never cached / never stratified
      lg <- .humann_long(humann_dir, layer, keep_stratified = FALSE)
    } else {
      if (is.null(.hum_long_cache[[layer]])) {
        .hum_long_cache[[layer]] <<- .humann_long(humann_dir, layer, keep_stratified = TRUE)
      }
      lg <- .hum_long_cache[[layer]]
    }
    lg[!lg$sample_id %in% .SAMPLE_FILTER$exclude, , drop = FALSE]
  }

  get_humann <- function(layer) .cached(layer, function() {
    cat(sprintf("\nBuilding HUMAnN %s table (community totals) ...\n", .layer_label[[layer]]))
    if (layer == "genefamilies") {
      cat("  WARNING: UniRef90 gene families can have 500k+ features.\n")
      cat(sprintf("  Applying strict filters: prevalence >= %.2f, mean >= %.4f\n",
                  humann_min_prevalence, humann_min_abundance))
    }
    tab <- humann_total_table(get_humann_long(layer), .layer_label[[layer]])
    cat(sprintf("  Median per-sample sum of %s: %.4g  (CPM-scale expected ~1e5-1e6; relative abundance <= 1)\n",
                .layer_label[[layer]], median(rowSums(tab))))
    filter_features(tab, humann_min_prevalence, humann_min_abundance,
                    .layer_label[[layer]])
  })

  get_humann_strat <- function(layer) .cached(paste0(layer, "_strat"), function() {
    lab <- paste0(.layer_label[[layer]], " (species-stratified)")
    cat(sprintf("\nBuilding HUMAnN %s table ...\n", lab))
    tab <- humann_strat_table(get_humann_long(layer), .layer_label[[layer]])
    filter_features(tab, strat_min_prevalence, strat_min_abundance, lab)
  })

  get_pathways     <- function() get_humann("pathways")
  get_genefamilies <- function() get_humann("genefamilies")
  get_ec           <- function() get_humann("ec")
  get_reactions    <- function() get_humann("reactions")
  get_ko           <- function() get_humann("ko")

  get_pathways_strat  <- function() get_humann_strat("pathways")
  get_ec_strat        <- function() get_humann_strat("ec")
  get_reactions_strat <- function() get_humann_strat("reactions")
  get_ko_strat        <- function() get_humann_strat("ko")

  get_food <- function() .cached("food", function() {
    cat("\nBuilding MEDI food abundance table ...\n")
    ft <- build_food_abundance_table(food_abundance_file, level = food_abundance_level)
    filter_features(ft, diet_min_prevalence, diet_min_abundance, "food features")
  })

  get_nutrients <- function() {
    if (is.null(.nutrient_tab)) {
      cat("\nBuilding nutrient content table ...\n")
      nt <- build_food_content_table(food_content_file, food_content_type)
      .nutrient_tab <<- filter_features(nt, diet_min_prevalence,
                                         diet_min_abundance, "nutrients")
    }
    .nutrient_tab
  }

  step <- 0L
  run_step <- function(label, expr) {
    step <<- step + 1L
    cat(sprintf("\n=== [%d] %s ===\n", step, label))
    expr
  }

  # Helper: run univariate screening for a pre-built feature table
  screen <- function(tab_getter, tab_label, out_subdir,
                     filt_prev = min_prevalence, filt_abund = min_abundance,
                     extra = NULL, norm = normalization) {
    if (!run_univariate || length(screening_variables) == 0) return(invisible(NULL))
    step <<- step + 1L
    cat(sprintf("\n=== [%d] Univariate screen: %s vs %d metadata variables ===\n",
                step, tab_label, length(screening_variables)))
    run_univariate_pipeline(
      taxa_tab        = tab_getter(),
      metadata_file   = metadata_file,
      output_dir      = file.path(output_dir, out_subdir),
      fixed_effects   = screening_variables,
      random_effects  = random_effects,
      min_prevalence  = filt_prev,
      min_abundance   = filt_abund,
      label           = tab_label,
      extra_meta_tab  = extra,
      normalization   = norm,
      transform       = transform,
      cores           = cores
    )
  }

  # Helper: run nutrient-vs-taxa with filtered taxa merged as predictors
  # Flipped design: taxa are the FEATURES, nutrients are the PREDICTORS.
  # This uses all samples with bacterial profiles (not just the overlap with
  # MEDI food-detected samples), avoids the massively underdetermined model
  # that results from putting 30+ species simultaneously as fixed effects,
  # and answers the more natural question: which bacteria associate with
  # estimated nutrient intake?
  # Nutrients are merged into the metadata as fixed effects alongside the
  # standard covariates (age, gender).
  run_nutrient_vs_taxa <- function(taxa_tab, taxa_label, out_subdir) {
    nutrient_filt <- filter_features(get_nutrients(), diet_min_prevalence,
                                      diet_min_abundance,
                                      paste("nutrients (for", taxa_label, "model)"))
    if (ncol(nutrient_filt) == 0) {
      cat(sprintf("  No nutrients passed filters -- skipping %s.\n", taxa_label))
      return(invisible(NULL))
    }
    run_taxa_metadata_pipeline(
      taxa_tab      = taxa_tab,
      metadata_file = metadata_file,
      output_dir    = file.path(output_dir, out_subdir),
      fixed_effects = c(colnames(nutrient_filt),
                        intersect(fixed_effects, c("age", "gender"))),
      random_effects = random_effects,
      min_prevalence = min_prevalence, min_abundance = min_abundance,
      label          = taxa_label,
      extra_meta_tab = nutrient_filt,
      normalization  = normalization, transform = transform, cores = cores)
  }

  # --- 1. Bacterial species vs metadata ------------------------------------
  if ("bacteria_species_metadata" %in% sel) {
    run_step("Bacterial species vs metadata",
      run_taxa_metadata_pipeline(
        taxa_tab = get_bacteria_species(), metadata_file = metadata_file,
        output_dir = file.path(output_dir, "bacteria_species_vs_metadata"),
        fixed_effects = fixed_effects, random_effects = random_effects,
        min_prevalence = min_prevalence, min_abundance = min_abundance,
        label = "bacterial species",
        normalization = normalization, transform = transform, cores = cores))
    screen(get_bacteria_species, "bacterial species",
           "bacteria_species_univariate")
  }

  # --- 2. Bacterial genus vs metadata --------------------------------------
  if ("bacteria_genus_metadata" %in% sel) {
    run_step("Bacterial genus vs metadata",
      run_taxa_metadata_pipeline(
        taxa_tab = get_bacteria_genus(), metadata_file = metadata_file,
        output_dir = file.path(output_dir, "bacteria_genus_vs_metadata"),
        fixed_effects = fixed_effects, random_effects = random_effects,
        min_prevalence = min_prevalence, min_abundance = min_abundance,
        label = "bacterial genera",
        normalization = normalization, transform = transform, cores = cores))
    screen(get_bacteria_genus, "bacterial genera",
           "bacteria_genus_univariate")
  }

  # --- 3. Archaeal species vs metadata -------------------------------------
  if ("archaea_metadata" %in% sel) {
    run_step("Archaeal species vs metadata",
      run_taxa_metadata_pipeline(
        taxa_tab = get_archaea(), metadata_file = metadata_file,
        output_dir = file.path(output_dir, "archaea_vs_metadata"),
        fixed_effects = fixed_effects, random_effects = random_effects,
        min_prevalence = min_prevalence, min_abundance = min_abundance,
        label = "archaeal species",
        normalization = normalization, transform = transform, cores = cores))
    screen(get_archaea, "archaeal species",
           "archaea_univariate")
  }

  # --- 4. Pathway abundances vs metadata -----------------------------------
  if ("pathways_metadata" %in% sel) {
    run_step("HUMAnN pathway abundances vs metadata",
      run_taxa_metadata_pipeline(
        taxa_tab = get_pathways(), metadata_file = metadata_file,
        output_dir = file.path(output_dir, "pathways_vs_metadata"),
        fixed_effects = fixed_effects, random_effects = random_effects,
        min_prevalence = min_prevalence, min_abundance = min_abundance,
        label = "pathway abundances",
        normalization = humann_normalization, transform = transform, cores = cores))
    screen(get_pathways, "pathway abundances",
           "pathways_univariate",
           filt_prev = humann_min_prevalence, filt_abund = humann_min_abundance,
           norm = humann_normalization)
  }

  # --- 5. Gene families vs metadata ----------------------------------------
  if ("genefamilies_metadata" %in% sel)
    run_step("HUMAnN gene families vs metadata",
      run_taxa_metadata_pipeline(
        taxa_tab = get_genefamilies(), metadata_file = metadata_file,
        output_dir = file.path(output_dir, "genefamilies_vs_metadata"),
        fixed_effects = fixed_effects, random_effects = random_effects,
        min_prevalence = min_prevalence, min_abundance = min_abundance,
        label = "gene families",
        normalization = humann_normalization, transform = transform, cores = cores))

  # --- 6. Reaction abundances vs metadata ----------------------------------
  if ("reactions_metadata" %in% sel) {
    run_step("HUMAnN reaction abundances vs metadata",
      run_taxa_metadata_pipeline(
        taxa_tab = get_reactions(), metadata_file = metadata_file,
        output_dir = file.path(output_dir, "reactions_vs_metadata"),
        fixed_effects = fixed_effects, random_effects = random_effects,
        min_prevalence = min_prevalence, min_abundance = min_abundance,
        label = "reactions",
        normalization = humann_normalization, transform = transform, cores = cores))
    screen(get_reactions, "reactions",
           "reactions_univariate",
           filt_prev = humann_min_prevalence, filt_abund = humann_min_abundance,
           norm = humann_normalization)
  }

  # --- 7. KO abundances vs metadata ----------------------------------------
  if ("ko_metadata" %in% sel) {
    run_step("HUMAnN KO abundances vs metadata",
      run_taxa_metadata_pipeline(
        taxa_tab = get_ko(), metadata_file = metadata_file,
        output_dir = file.path(output_dir, "ko_vs_metadata"),
        fixed_effects = fixed_effects, random_effects = random_effects,
        min_prevalence = min_prevalence, min_abundance = min_abundance,
        label = "KO abundances",
        normalization = humann_normalization, transform = transform, cores = cores))
    screen(get_ko, "KO abundances",
           "ko_univariate",
           filt_prev = humann_min_prevalence, filt_abund = humann_min_abundance,
           norm = humann_normalization)
  }

  # --- 7b. EC abundances vs metadata ---------------------------------------
  if ("ec_metadata" %in% sel) {
    run_step("HUMAnN EC abundances vs metadata",
      run_taxa_metadata_pipeline(
        taxa_tab = get_ec(), metadata_file = metadata_file,
        output_dir = file.path(output_dir, "ec_vs_metadata"),
        fixed_effects = fixed_effects, random_effects = random_effects,
        min_prevalence = min_prevalence, min_abundance = min_abundance,
        label = "EC abundances",
        normalization = humann_normalization, transform = transform, cores = cores))
    screen(get_ec, "EC abundances", "ec_univariate",
           filt_prev = humann_min_prevalence, filt_abund = humann_min_abundance,
           norm = humann_normalization)
  }

  # --- 7c. Species-assigned (stratified) functions vs metadata -------------
  # Columns are <FEATURE>__IN__<Species>; each is a function's abundance
  # contributed by one species. Multivariate model + univariate screen.
  .run_strat_vs_metadata <- function(getter, label, out_subdir, univ_subdir) {
    run_step(sprintf("HUMAnN %s (species-stratified) vs metadata", label),
      run_taxa_metadata_pipeline(
        taxa_tab = getter(), metadata_file = metadata_file,
        output_dir = file.path(output_dir, out_subdir),
        fixed_effects = fixed_effects, random_effects = random_effects,
        min_prevalence = strat_min_prevalence, min_abundance = strat_min_abundance,
        label = paste(label, "(species-stratified)"),
        normalization = humann_normalization, transform = transform, cores = cores))
    screen(getter, paste(label, "(species-stratified)"), univ_subdir,
           filt_prev = strat_min_prevalence, filt_abund = strat_min_abundance,
           norm = humann_normalization)
  }

  if ("pathways_strat_metadata" %in% sel)
    .run_strat_vs_metadata(get_pathways_strat, "pathway abundances",
                           "pathways_strat_vs_metadata", "pathways_strat_univariate")
  if ("ec_strat_metadata" %in% sel)
    .run_strat_vs_metadata(get_ec_strat, "EC abundances",
                           "ec_strat_vs_metadata", "ec_strat_univariate")
  if ("reactions_strat_metadata" %in% sel)
    .run_strat_vs_metadata(get_reactions_strat, "reactions",
                           "reactions_strat_vs_metadata", "reactions_strat_univariate")
  if ("ko_strat_metadata" %in% sel)
    .run_strat_vs_metadata(get_ko_strat, "KO abundances",
                           "ko_strat_vs_metadata", "ko_strat_univariate")

  # --- 8. Nutrients vs metadata --------------------------------------------
  if ("nutrient_metadata" %in% sel) {
    run_step("Nutrients vs metadata",
      run_taxa_metadata_pipeline(
        taxa_tab = get_nutrients(), metadata_file = metadata_file,
        output_dir = file.path(output_dir, "nutrient_vs_metadata"),
        fixed_effects = fixed_effects, random_effects = random_effects,
        min_prevalence = diet_min_prevalence, min_abundance = diet_min_abundance,
        label = "nutrients",
        normalization = normalization, transform = transform, cores = cores))
    screen(get_nutrients, "nutrients",
           "nutrients_univariate",
           filt_prev = diet_min_prevalence, filt_abund = diet_min_abundance)
  }

  # --- 9. Nutrients vs bacterial species -----------------------------------
  if ("nutrient_bacteria_species" %in% sel)
    run_step("Nutrients vs bacterial species",
      run_nutrient_vs_taxa(get_bacteria_species(),
                            "bacterial species", "nutrient_vs_bacteria_species"))

  # --- 10. Nutrients vs bacterial genus -------------------------------------
  if ("nutrient_bacteria_genus" %in% sel)
    run_step("Nutrients vs bacterial genus",
      run_nutrient_vs_taxa(get_bacteria_genus(),
                            "bacterial genera", "nutrient_vs_bacteria_genus"))

  # --- 11. Nutrients vs archaeal species ------------------------------------
  if ("nutrient_archaea" %in% sel)
    run_step("Nutrients vs archaeal species",
      run_nutrient_vs_taxa(get_archaea(),
                            "archaeal species", "nutrient_vs_archaea"))

  # -----------------------------------------------------------------------
  # Functional layer <-> diet associations
  # -----------------------------------------------------------------------
  # Each nutrient is tested univariately against the functional feature table.
  # Uses run_univariate_pipeline with nutrients as the screening variables.

  .run_nutrient_vs_functional <- function(feat_getter, feat_label, out_subdir,
                                           filt_prev, filt_abund) {
    if (is.null(food_content_file)) return(invisible(NULL))
    nutrient_tab <- filter_features(get_nutrients(), diet_min_prevalence,
                                     diet_min_abundance, "nutrients")
    if (ncol(nutrient_tab) == 0) {
      cat("  No nutrients passed filters -- skipping.
"); return(invisible(NULL))
    }
    run_univariate_pipeline(
      taxa_tab       = feat_getter(),
      metadata_file  = metadata_file,
      output_dir     = file.path(output_dir, out_subdir),
      fixed_effects  = colnames(nutrient_tab),
      random_effects = random_effects,
      min_prevalence = filt_prev, min_abundance = filt_abund,
      label          = feat_label,
      extra_meta_tab = nutrient_tab,
      normalization  = humann_normalization, transform = transform, cores = cores
    )
  }

  if ("nutrient_vs_pathways" %in% sel)
    run_step("Nutrients vs pathway abundances (univariate)",
      .run_nutrient_vs_functional(get_pathways, "pathway abundances",
        "nutrient_vs_pathways", humann_min_prevalence, humann_min_abundance))

  if ("nutrient_vs_reactions" %in% sel)
    run_step("Nutrients vs reaction abundances (univariate)",
      .run_nutrient_vs_functional(get_reactions, "reactions",
        "nutrient_vs_reactions", humann_min_prevalence, humann_min_abundance))

  if ("nutrient_vs_ko" %in% sel)
    run_step("Nutrients vs KO abundances (univariate)",
      .run_nutrient_vs_functional(get_ko, "KO abundances",
        "nutrient_vs_ko", humann_min_prevalence, humann_min_abundance))

  if ("nutrient_vs_ec" %in% sel)
    run_step("Nutrients vs EC abundances (univariate)",
      .run_nutrient_vs_functional(get_ec, "EC abundances",
        "nutrient_vs_ec", humann_min_prevalence, humann_min_abundance))

  # Species-assigned functions vs nutrients
  if ("nutrient_vs_pathways_strat" %in% sel)
    run_step("Nutrients vs species-stratified pathway abundances (univariate)",
      .run_nutrient_vs_functional(get_pathways_strat, "pathway abundances (species-stratified)",
        "nutrient_vs_pathways_strat", strat_min_prevalence, strat_min_abundance))
  if ("nutrient_vs_ec_strat" %in% sel)
    run_step("Nutrients vs species-stratified EC abundances (univariate)",
      .run_nutrient_vs_functional(get_ec_strat, "EC abundances (species-stratified)",
        "nutrient_vs_ec_strat", strat_min_prevalence, strat_min_abundance))
  if ("nutrient_vs_reactions_strat" %in% sel)
    run_step("Nutrients vs species-stratified reactions (univariate)",
      .run_nutrient_vs_functional(get_reactions_strat, "reactions (species-stratified)",
        "nutrient_vs_reactions_strat", strat_min_prevalence, strat_min_abundance))
  if ("nutrient_vs_ko_strat" %in% sel)
    run_step("Nutrients vs species-stratified KO abundances (univariate)",
      .run_nutrient_vs_functional(get_ko_strat, "KO abundances (species-stratified)",
        "nutrient_vs_ko_strat", strat_min_prevalence, strat_min_abundance))

  # -----------------------------------------------------------------------
  # Functional layer <-> taxonomy associations
  # Each pathway/reaction/KO is tested univariately against the taxa table.
  # -----------------------------------------------------------------------
  .run_functional_vs_taxa <- function(taxa_getter, taxa_label, taxa_prev, taxa_abund,
                                       func_getter, func_label, func_prev, func_abund,
                                       out_subdir) {
    func_tab <- filter_features(func_getter(), func_prev, func_abund, func_label)
    if (ncol(func_tab) == 0) {
      cat(sprintf("  No %s passed filters -- skipping.
", func_label))
      return(invisible(NULL))
    }
    run_univariate_pipeline(
      taxa_tab       = taxa_getter(),
      metadata_file  = metadata_file,
      output_dir     = file.path(output_dir, out_subdir),
      fixed_effects  = colnames(func_tab),
      random_effects = random_effects,
      min_prevalence = taxa_prev, min_abundance = taxa_abund,
      label          = taxa_label,
      extra_meta_tab = func_tab,
      normalization  = normalization, transform = transform, cores = cores
    )
  }

  if ("pathways_vs_bacteria_species" %in% sel)
    run_step("Pathway abundances vs bacterial species (univariate)",
      .run_functional_vs_taxa(get_bacteria_species, "bacterial species",
        min_prevalence, min_abundance,
        get_pathways, "pathway abundances",
        humann_min_prevalence, humann_min_abundance,
        "pathways_vs_bacteria_species"))

  if ("pathways_vs_bacteria_genus" %in% sel)
    run_step("Pathway abundances vs bacterial genera (univariate)",
      .run_functional_vs_taxa(get_bacteria_genus, "bacterial genera",
        min_prevalence, min_abundance,
        get_pathways, "pathway abundances",
        humann_min_prevalence, humann_min_abundance,
        "pathways_vs_bacteria_genus"))

  if ("pathways_vs_archaea" %in% sel)
    run_step("Pathway abundances vs archaeal species (univariate)",
      .run_functional_vs_taxa(get_archaea, "archaeal species",
        min_prevalence, min_abundance,
        get_pathways, "pathway abundances",
        humann_min_prevalence, humann_min_abundance,
        "pathways_vs_archaea"))

  if ("reactions_vs_bacteria_species" %in% sel)
    run_step("Reaction abundances vs bacterial species (univariate)",
      .run_functional_vs_taxa(get_bacteria_species, "bacterial species",
        min_prevalence, min_abundance,
        get_reactions, "reactions",
        humann_min_prevalence, humann_min_abundance,
        "reactions_vs_bacteria_species"))

  if ("reactions_vs_bacteria_genus" %in% sel)
    run_step("Reaction abundances vs bacterial genera (univariate)",
      .run_functional_vs_taxa(get_bacteria_genus, "bacterial genera",
        min_prevalence, min_abundance,
        get_reactions, "reactions",
        humann_min_prevalence, humann_min_abundance,
        "reactions_vs_bacteria_genus"))

  if ("reactions_vs_archaea" %in% sel)
    run_step("Reaction abundances vs archaeal species (univariate)",
      .run_functional_vs_taxa(get_archaea, "archaeal species",
        min_prevalence, min_abundance,
        get_reactions, "reactions",
        humann_min_prevalence, humann_min_abundance,
        "reactions_vs_archaea"))

  if ("ko_vs_bacteria_species" %in% sel)
    run_step("KO abundances vs bacterial species (univariate)",
      .run_functional_vs_taxa(get_bacteria_species, "bacterial species",
        min_prevalence, min_abundance,
        get_ko, "KO abundances",
        humann_min_prevalence, humann_min_abundance,
        "ko_vs_bacteria_species"))

  if ("ko_vs_bacteria_genus" %in% sel)
    run_step("KO abundances vs bacterial genera (univariate)",
      .run_functional_vs_taxa(get_bacteria_genus, "bacterial genera",
        min_prevalence, min_abundance,
        get_ko, "KO abundances",
        humann_min_prevalence, humann_min_abundance,
        "ko_vs_bacteria_genus"))

  if ("ko_vs_archaea" %in% sel)
    run_step("KO abundances vs archaeal species (univariate)",
      .run_functional_vs_taxa(get_archaea, "archaeal species",
        min_prevalence, min_abundance,
        get_ko, "KO abundances",
        humann_min_prevalence, humann_min_abundance,
        "ko_vs_archaea"))

  if ("ec_vs_bacteria_species" %in% sel)
    run_step("EC abundances vs bacterial species (univariate)",
      .run_functional_vs_taxa(get_bacteria_species, "bacterial species",
        min_prevalence, min_abundance,
        get_ec, "EC abundances",
        humann_min_prevalence, humann_min_abundance,
        "ec_vs_bacteria_species"))

  if ("ec_vs_bacteria_genus" %in% sel)
    run_step("EC abundances vs bacterial genera (univariate)",
      .run_functional_vs_taxa(get_bacteria_genus, "bacterial genera",
        min_prevalence, min_abundance,
        get_ec, "EC abundances",
        humann_min_prevalence, humann_min_abundance,
        "ec_vs_bacteria_genus"))

  if ("ec_vs_archaea" %in% sel)
    run_step("EC abundances vs archaeal species (univariate)",
      .run_functional_vs_taxa(get_archaea, "archaeal species",
        min_prevalence, min_abundance,
        get_ec, "EC abundances",
        humann_min_prevalence, humann_min_abundance,
        "ec_vs_archaea"))

  # -----------------------------------------------------------------------
  # MEDI food abundance
  #   food_metadata : food features (response) ~ metadata
  #   food_vs_*     : each food feature (predictor, one at a time) vs every
  #                   taxon / function (response), pooled + global BH.
  #   Univariate by design: a model with 100+ foods as simultaneous fixed
  #   effects is unidentifiable at this sample size.
  # -----------------------------------------------------------------------
  .run_food_univariate <- function(resp_getter, resp_label, out_subdir,
                                    resp_prev, resp_abund, resp_norm) {
    if (is.null(food_abundance_file)) return(invisible(NULL))
    food_tab <- get_food()
    if (ncol(food_tab) == 0) {
      cat("  No food features passed filters -- skipping.\n"); return(invisible(NULL))
    }
    run_univariate_pipeline(
      taxa_tab       = resp_getter(),
      metadata_file  = metadata_file,
      output_dir     = file.path(output_dir, out_subdir),
      fixed_effects  = colnames(food_tab),
      random_effects = random_effects,
      min_prevalence = resp_prev, min_abundance = resp_abund,
      label          = resp_label,
      extra_meta_tab = food_tab,
      normalization  = resp_norm, transform = transform, cores = cores
    )
  }
  .food_vs <- function(key, step_label, getter, resp_label, out_subdir,
                       resp_prev, resp_abund, resp_norm)
    if (key %in% sel)
      run_step(step_label,
        .run_food_univariate(getter, resp_label, out_subdir,
                             resp_prev, resp_abund, resp_norm))

  if ("food_metadata" %in% sel) {
    run_step("Food abundance vs metadata",
      run_taxa_metadata_pipeline(
        taxa_tab = get_food(), metadata_file = metadata_file,
        output_dir = file.path(output_dir, "food_vs_metadata"),
        fixed_effects = fixed_effects, random_effects = random_effects,
        min_prevalence = diet_min_prevalence, min_abundance = diet_min_abundance,
        label = "food features",
        normalization = food_normalization, transform = transform, cores = cores))
    screen(get_food, "food features", "food_univariate",
           filt_prev = diet_min_prevalence, filt_abund = diet_min_abundance,
           norm = food_normalization)
  }

  # Taxa as response (taxa keep the standard normalisation)
  .food_vs("food_vs_bacteria_species", "Food vs bacterial species (univariate)",
           get_bacteria_species, "bacterial species", "food_vs_bacteria_species",
           min_prevalence, min_abundance, normalization)
  .food_vs("food_vs_bacteria_genus", "Food vs bacterial genera (univariate)",
           get_bacteria_genus, "bacterial genera", "food_vs_bacteria_genus",
           min_prevalence, min_abundance, normalization)
  .food_vs("food_vs_archaea", "Food vs archaeal species (univariate)",
           get_archaea, "archaeal species", "food_vs_archaea",
           min_prevalence, min_abundance, normalization)

  # HUMAnN functions as response (already CPM -> humann_normalization)
  .food_vs("food_vs_pathways", "Food vs pathway abundances (univariate)",
           get_pathways, "pathway abundances", "food_vs_pathways",
           humann_min_prevalence, humann_min_abundance, humann_normalization)
  .food_vs("food_vs_ec", "Food vs EC abundances (univariate)",
           get_ec, "EC abundances", "food_vs_ec",
           humann_min_prevalence, humann_min_abundance, humann_normalization)
  .food_vs("food_vs_reactions", "Food vs reactions (univariate)",
           get_reactions, "reactions", "food_vs_reactions",
           humann_min_prevalence, humann_min_abundance, humann_normalization)
  .food_vs("food_vs_ko", "Food vs KO abundances (univariate)",
           get_ko, "KO abundances", "food_vs_ko",
           humann_min_prevalence, humann_min_abundance, humann_normalization)

  # Species-assigned functions as response
  .food_vs("food_vs_pathways_strat", "Food vs species-stratified pathways (univariate)",
           get_pathways_strat, "pathway abundances (species-stratified)",
           "food_vs_pathways_strat",
           strat_min_prevalence, strat_min_abundance, humann_normalization)
  .food_vs("food_vs_ec_strat", "Food vs species-stratified ECs (univariate)",
           get_ec_strat, "EC abundances (species-stratified)", "food_vs_ec_strat",
           strat_min_prevalence, strat_min_abundance, humann_normalization)
  .food_vs("food_vs_reactions_strat", "Food vs species-stratified reactions (univariate)",
           get_reactions_strat, "reactions (species-stratified)", "food_vs_reactions_strat",
           strat_min_prevalence, strat_min_abundance, humann_normalization)
  .food_vs("food_vs_ko_strat", "Food vs species-stratified KOs (univariate)",
           get_ko_strat, "KO abundances (species-stratified)", "food_vs_ko_strat",
           strat_min_prevalence, strat_min_abundance, humann_normalization)

  cat(sprintf("\nAll analyses complete. Outputs under: %s/\n\n", output_dir))
}


# =============================================================================
# 15. CLI entry point
# =============================================================================
if (!interactive()) {
  option_list <- list(
    make_option("--metaphlan_dir",  type = "character", default = NULL,
                help = "Directory of *_metaphlan_profile_wreads.txt files (bacteria species + genus)"),
    make_option("--archaea_dir",    type = "character", default = NULL,
                help = "Directory of *_archaea_bracken_abundance.c01.txt files"),
    make_option("--humann_dir",     type = "character", default = NULL,
                help = paste0("HUMAnN directory with per-sample sub-directories SXXX/ containing ",
                              "out.SXXX_combined_pathabundance_normalised.tsv, out.SXXX_combined_genefamilies-cpm.tsv and ",
                              "out.SXXX_combined_{ec,rxn,ko}-cpm-named.tsv")),
    make_option("--food_content",   type = "character", default = NULL,
                help = "MEDI food_content.csv"),
    make_option("--food_abundance", type = "character", default = NULL,
                help = paste0("MEDI food_abundance.csv (relative abundance of dietary taxa). ",
                              "Each food feature is tested against taxa and HUMAnN functions.")),
    make_option("--food_abundance_level", type = "character", default = "taxon",
                help = paste0("Aggregation of food_abundance: 'taxon' (de-duplicated organisms), ",
                              "'food_group', or 'food' (one column per food item) [default: taxon]")),
    make_option("--food_normalization", type = "character", default = "NONE",
                help = "Maaslin3 normalisation when food features are the response [default: NONE]"),
    make_option("--diversity_file", type = "character", default = NULL,
                help = paste0("diversity_stats_per_sample.csv from diversity_analysis_outputs/. ",
                              "Alpha diversity metrics (Observed, shannon, inverse_simpson, chao1) ",
                              "are merged into metadata and added to univariate screening. ",
                              "Sample IDs with _metaphlan_wreads suffix are stripped automatically.")),
    make_option("--food_content_type", type = "character", default = NULL,
                help = "Comma-separated source_type filter for food_content, e.g. Nutrient"),
    make_option("--metadata",       type = "character",
                help = "TSV/CSV metadata file (first column = sample_id)"),
    make_option("--output_dir",     type = "character", default = "maaslin3_output",
                help = "Output directory [default: maaslin3_output]"),
    make_option("--analyses",       type = "character", default = "auto",
                help = paste0("'auto', 'all', or comma-separated from:\n",
                              "  bacteria_species_metadata, bacteria_genus_metadata,\n",
                              "  archaea_metadata,\n",
                              "  pathways_metadata, genefamilies_metadata,\n",
                              "  reactions_metadata, ko_metadata, ec_metadata,\n",
                              "  pathways_strat_metadata, ec_strat_metadata,\n",
                              "  reactions_strat_metadata, ko_strat_metadata,\n",
                              "  nutrient_metadata,\n",
                              "  nutrient_bacteria_species, nutrient_bacteria_genus,\n",
                              "  nutrient_archaea,\n",
                              "  nutrient_vs_pathways, nutrient_vs_reactions, nutrient_vs_ko,\n",
                              "  nutrient_vs_ec, nutrient_vs_pathways_strat, nutrient_vs_ec_strat,\n",
                              "  nutrient_vs_reactions_strat, nutrient_vs_ko_strat,\n",
                              "  pathways_vs_bacteria_species, pathways_vs_bacteria_genus,\n",
                              "  pathways_vs_archaea,\n",
                              "  reactions_vs_bacteria_species, reactions_vs_bacteria_genus,\n",
                              "  reactions_vs_archaea,\n",
                              "  ko_vs_bacteria_species, ko_vs_bacteria_genus,\n",
                              "  ko_vs_archaea,\n",
                              "  ec_vs_bacteria_species, ec_vs_bacteria_genus, ec_vs_archaea,\n",
                              "  food_metadata, food_vs_bacteria_species, food_vs_bacteria_genus,\n",
                              "  food_vs_archaea, food_vs_pathways, food_vs_ec, food_vs_reactions,\n",
                              "  food_vs_ko, food_vs_{pathways,ec,reactions,ko}_strat")),
    make_option("--run_univariate", action = "store_true", default = FALSE,
                help = paste0("Run univariate metadata screening for every feature table. ",
                              "Each variable in --screening_variables is tested ",
                              "independently; results are pooled and BH-corrected globally. ",
                              "[default: FALSE]")),
    make_option("--screening_variables", type = "character", default = NULL,
                help = paste0("Comma-separated metadata variables for univariate screening. ",
                              "Defaults to the built-in SCREENING_VARIABLES list (44 variables) ",
                              "when --run_univariate is set and this is not supplied.")),
    make_option("--fixed_effects",  type = "character",
                default = "wfh_zscore,hfa_zscore,wfa_zscore,age,gender,excl_brstfed,keep_animals,diarrhoea",
                help = "Comma-separated fixed effects for *_metadata analyses"),
    make_option("--random_effects", type = "character", default = NULL,
                help = "Comma-separated random effects [default: none]"),
    make_option("--min_prevalence", type = "double",  default = 0.10,
                help = "Min prevalence for taxa/pathways/genes [default: 0.10]"),
    make_option("--min_abundance",  type = "double",  default = 0.01,
                help = "Min mean abundance for taxa/pathways/genes [default: 0.01]"),
    make_option("--diet_min_prevalence", type = "double", default = 0.10,
                help = "Min prevalence for nutrient features [default: 0.10]"),
    make_option("--diet_min_abundance",  type = "double", default = 0,
                help = "Min mean for nutrient features [default: 0]"),
    make_option("--humann_min_prevalence", type = "double", default = 0.30,
                help = paste0("Min prevalence for HUMAnN features (pathways, gene families, ",
                              "reactions, KOs). Higher values reduce feature count and run ",
                              "time. Gene families may need 0.5+ [default: 0.30]")),
    make_option("--humann_min_abundance",  type = "double", default = 0.001,
                help = "Min mean abundance for HUMAnN features [default: 0.001]"),
    make_option("--strat_min_prevalence", type = "double", default = 0.10,
                help = "Min prevalence for species-stratified HUMAnN features [default: 0.10]"),
    make_option("--strat_min_abundance",  type = "double", default = 0.001,
                help = "Min mean abundance for species-stratified HUMAnN features [default: 0.001]"),
    make_option("--humann_normalization", type = "character", default = "NONE",
                help = paste0("Maaslin3 normalisation for HUMAnN tables used as the response ",
                              "(already CPM, so NONE avoids re-normalising) [default: NONE]")),
    make_option("--anthro_file",    type = "character", default = NULL,
                help = paste0("WHO anthropometric survey CSV. Any row with a non-empty 'Flag' ",
                              "has its sample excluded from ALL analyses.")),
    make_option("--anthro_id_col",  type = "character", default = "lab_id",
                help = "ID column in --anthro_file [default: lab_id]"),
    make_option("--metadata_id_col", type = "character", default = "lab_id",
                help = "Metadata column matching --anthro_id_col [default: lab_id]"),
    make_option("--anthro_flag_col", type = "character", default = "Flag",
                help = "Flag column in --anthro_file [default: Flag]"),
    make_option("--normalization",  type = "character", default = "TSS",
                help = "Maaslin3 normalisation (TSS, CLR, CSS, NONE) [default: TSS]"),
    make_option("--transform",      type = "character", default = "LOG",
                help = "Maaslin3 transform (LOG, LOGIT, AST, NONE) [default: LOG]"),
    make_option("--cores",          type = "integer",  default = 4,
                help = "CPU cores [default: 4]")
  )

  opt <- parse_args(OptionParser(
    usage       = "%prog [options]",
    description = paste0("\nMaaslin3 analysis: MetaPhlAn bacteria (species + genus) + ",
                         "Bracken archaea + HUMAnN pathways/genes + MEDI nutrients.\n"),
    option_list = option_list
  ))

  if (is.null(opt$metadata))
    stop("--metadata is required.", call. = FALSE)
  if (is.null(opt$metaphlan_dir) && is.null(opt$archaea_dir) &&
      is.null(opt$humann_dir)    && is.null(opt$food_content) &&
      is.null(opt$food_abundance))
    stop("At least one of --metaphlan_dir, --archaea_dir, --humann_dir, ",
         "--food_content must be supplied.", call. = FALSE)

  run_medi_maaslin3_pipeline(
    metaphlan_dir         = opt$metaphlan_dir,
    archaea_dir           = opt$archaea_dir,
    humann_dir            = opt$humann_dir,
    food_content_file     = opt$food_content,
    food_abundance_file   = opt$food_abundance,
    food_abundance_level  = opt$food_abundance_level,
    food_normalization    = opt$food_normalization,
    diversity_file        = opt$diversity_file,
    metadata_file         = opt$metadata,
    output_dir            = opt$output_dir,
    analyses              = opt$analyses,
    run_univariate        = opt$run_univariate,
    screening_variables   = if (!is.null(opt$screening_variables))
                              strsplit(opt$screening_variables, ",")[[1]]
                            else SCREENING_VARIABLES,
    fixed_effects         = strsplit(opt$fixed_effects, ",")[[1]],
    random_effects        = if (!is.null(opt$random_effects))
                              strsplit(opt$random_effects, ",")[[1]] else NULL,
    food_content_type     = if (!is.null(opt$food_content_type))
                              strsplit(opt$food_content_type, ",")[[1]] else NULL,
    min_prevalence        = opt$min_prevalence,
    min_abundance         = opt$min_abundance,
    diet_min_prevalence   = opt$diet_min_prevalence,
    diet_min_abundance    = opt$diet_min_abundance,
    humann_min_prevalence = opt$humann_min_prevalence,
    humann_min_abundance  = opt$humann_min_abundance,
    strat_min_prevalence  = opt$strat_min_prevalence,
    strat_min_abundance   = opt$strat_min_abundance,
    humann_normalization  = opt$humann_normalization,
    anthro_file           = opt$anthro_file,
    anthro_id_col         = opt$anthro_id_col,
    metadata_id_col       = opt$metadata_id_col,
    anthro_flag_col       = opt$anthro_flag_col,
    normalization         = opt$normalization,
    transform             = opt$transform,
    cores                 = opt$cores
  )
}

# =============================================================================
# 16. Build food organism abundance table from MEDI S_counts.csv
#     These are the eukaryotic/dietary species detected by Kraken2 -- the
#     complement of the bacterial/archaeal taxa table. We keep everything that
#     is NOT in the microbial kingdoms so we're working with the diet signal.
#     Sample IDs are stripped of their level prefix (S_S003 -> S003).
# =============================================================================
build_food_organism_table <- function(taxa_counts_file,
                                       microbial_kingdoms = c("K__Bacteria",
                                                              "K__Archaea")) {
  counts <- read_csv(taxa_counts_file, show_col_types = FALSE)

  required <- c("sample_id", "name", "fraction_total_reads", "lineage")
  missing  <- setdiff(required, colnames(counts))
  if (length(missing) > 0)
    stop("taxa_counts file missing column(s): ", paste(missing, collapse = ", "))

  food_rows <- counts %>%
    mutate(
      kingdom_token = str_extract(lineage, "^[^;]+"),
      sample_id     = sub("^[^_]+_", "", sample_id),   # S_S003 -> S003
      taxon         = make.names(name)
    ) %>%
    filter(!kingdom_token %in% microbial_kingdoms)

  if (nrow(food_rows) == 0)
    stop("No non-microbial rows found in taxa_counts file. ",
         "Check that it is a MEDI S_counts/G_counts file containing dietary reads.")

  wide <- food_rows %>%
    group_by(sample_id, taxon) %>%
    summarise(abundance = sum(fraction_total_reads, na.rm = TRUE),
              .groups = "drop") %>%
    pivot_wider(names_from = taxon, values_from = abundance,
                values_fill = 0) %>%
    column_to_rownames("sample_id")

  mat <- as.matrix(wide); storage.mode(mat) <- "double"
  wide <- as.data.frame(mat)

  wide <- keep_study_samples(wide, "food organisms")
  cat(sprintf("  Food organism table: %d samples, %d dietary species.\n",
              nrow(wide), ncol(wide)))
  wide
}

# =============================================================================
# 17. Food organisms vs nutrient content pipeline
#
#     Tests whether the relative abundance of each detected dietary species
#     (e.g. Zea mays, Bos taurus, Thunnus thynnus) predicts the estimated
#     nutrient / compound content of that sample. Each food organism is tested
#     univariately against the nutrient table (same approach as diet-vs-taxa).
#     This validates biological plausibility: maize reads -> carbohydrates,
#     bovine/fish reads -> protein/fatty acids.
# =============================================================================
run_food_organism_nutrient_pipeline <- function(
    food_org_tab,
    nutrient_tab,
    output_dir,
    normalization = "TSS",
    transform     = "LOG",
    cores         = 4
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  common <- intersect(rownames(food_org_tab), rownames(nutrient_tab))
  if (length(common) == 0)
    stop("No common samples between food organism table and nutrient table.")

  food_org_tab <- food_org_tab[common, , drop = FALSE]
  nutrient_tab <- nutrient_tab[common, , drop = FALSE]

  cat(sprintf("  Food organisms vs nutrients: %d samples | %d organisms | %d nutrients.\n",
              length(common), ncol(food_org_tab), ncol(nutrient_tab)))

  results_list <- list()

  for (org in colnames(food_org_tab)) {
    org_dir <- file.path(output_dir, "per_organism", make.names(org))
    meta <- data.frame(organism_abundance = food_org_tab[[org]],
                       row.names = common, check.names = FALSE)

    fit <- tryCatch(
      run_maaslin3(nutrient_tab, meta, org_dir,
                   fixed_effects = "organism_abundance",
                   normalization = normalization,
                   transform     = transform,
                   cores         = cores),
      error = function(e) {
        warning(sprintf("Maaslin3 failed for organism '%s': %s",
                        org, conditionMessage(e)))
        NULL
      }
    )

    if (!is.null(fit) && !is.null(fit$results) && nrow(fit$results) > 0) {
      # Extract rows for the organism predictor
      rows <- fit$results
      for (col in c("metadata", "value", "variable", "term")) {
        if (col %in% colnames(rows) &&
            "organism_abundance" %in% rows[[col]]) {
          rows <- rows[rows[[col]] == "organism_abundance", , drop = FALSE]
          break
        }
      }
      if (nrow(rows) > 0) {
        rows$food_organism <- org
        results_list[[org]] <- rows
      }
    }
  }

  if (length(results_list) == 0) {
    cat("  No food-organism-nutrient associations returned.\n")
    return(invisible(NULL))
  }

  combined <- bind_rows(results_list)
  pcol     <- .pval_col(combined)
  combined <- combined %>%
    mutate(qval_global_BH = p.adjust(.data[[pcol]], method = "BH")) %>%
    arrange(qval_global_BH)

  sig <- combined %>% filter(qval_global_BH < 0.25)
  write_tsv(combined, file.path(output_dir, "food_org_nutrient_all_results.tsv"))
  write_tsv(sig,      file.path(output_dir, "food_org_nutrient_significant_results.tsv"))
  cat(sprintf("  Associations tested: %d | Significant (BH q < 0.25): %d\n",
              nrow(combined), nrow(sig)))
  invisible(combined)
}