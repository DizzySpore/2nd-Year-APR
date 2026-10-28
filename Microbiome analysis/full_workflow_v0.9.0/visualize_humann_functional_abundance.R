#!/usr/bin/env Rscript

################################################################################
# visualize_humann_functional_abundance.R
################################################################################
#
# Purpose:
#   Stacked bar charts of HUMAnN3 functional (pathway) abundance:
#     1. Whole-sample, read-based profiles from STAGE_20
#        (humann_reads_output/) -- one bar per sample.
#     2. Per-MAG profiles from STAGE_21 (humann_contigs_output/) -- one bar
#        per bin, faceted by sample, labelled with the GTDB-Tk species call
#        pulled from STAGE_21's manifest file.
#
#   NOTE ON NAMING: STAGE_21 was reworked to run HUMAnN3 on reads pulled
#   per-MAG (bin), not on assembly contigs -- "per contig" in the original
#   ask is treated here as "per bin", since that's what STAGE_21 now
#   actually produces. If you genuinely want per-contig functional calls
#   (rather than per-bin), that's a different HUMAnN3 input and this
#   script would need a different STAGE_21 reader.
#
# Usage:
#   Edit the CONFIG block below to match your directory layout, then:
#     Rscript visualize_humann_functional_abundance.R
#
# Requires: dplyr, tidyr, readr, stringr, ggplot2, scales, purrr
#   install.packages(c("dplyr","tidyr","readr","stringr","ggplot2","scales","purrr"))
#
################################################################################

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
    library(stringr)
    library(ggplot2)
    library(scales)
    library(purrr)
})

################################################################################
# CONFIG -- edit these to match your run
################################################################################

# Matches ${DIRECTORY} in the pipeline
BASE_DIR <- "/rds/projects/h/hallly-microbiome/Leon/malawi_data/attempt2/work"

# NOTE: these are the top-level, ALL-SAMPLES directories -- i.e. one level
# up from the pipeline's own per-sample HUMANN_READS_DIR / HUMANN_CONTIGS_DIR
# variables (which point at .../humann_reads_output/${SAMPLE}). This script
# scans across every sample subdirectory underneath, so point it at the
# parent, not a single sample's folder.
HUMANN_READS_ROOT   <- file.path(BASE_DIR, "humann_reads_output")     # STAGE_20 output (per sample subdir)
HUMANN_CONTIGS_ROOT <- file.path(BASE_DIR, "humann_contigs_output")   # STAGE_21 output (per sample subdir, per-bin files)
CROSSCHECK_ROOT      <- file.path(BASE_DIR, "crosscheck_out")         # STAGE_12 output (used only for taxon-focus plots)
OUTPUT_DIR          <- file.path(BASE_DIR, "visualizations", "humann")

# "pathabundance" (interpretable pathway names) or "genefamilies" (UniRef IDs).
# Pathway abundance is the sane default for a readable stacked bar chart.
PROFILE_TYPE <- "pathabundance"

# Pathway-level UNMAPPED/UNINTEGRATED currently dominate these samples'
# abundance so heavily that they crowd out everything else in a stacked bar.
# Excluded from every plot below for now. This is a MetaCyc-pathway-level
# problem specifically -- gene families (UniRef, or translated to PFAM/EC/
# etc, which HUMAnN3 can do via its regroup utility) generally have a much
# lower unclassified fraction, since building a full pathway requires several
# consecutive steps to all be identified, and PFAM/EC-level classification
# doesn't require that. Worth revisiting with a gene-family-based view later;
# for now, set to character(0) to bring UNMAPPED/UNINTEGRATED back if needed.
EXCLUDE_FEATURES <- c("UNMAPPED", "UNINTEGRATED")

# How many top (non-excluded) features to show individually before lumping
# the remainder into "Other".
TOP_N <- 15

# Truncate long pathway descriptions in the legend to this many characters
LEGEND_LABEL_WIDTH <- 45

# STAGE_21's per-bin plot only: as more samples finish, the number of
# sample x bin panels grows and the facet grid gets unreadable, and most of
# those bins are one-off genomes that don't recur elsewhere anyway. Keep only
# species (by GTDB-Tk call) present in at least this fraction of samples
# CURRENTLY in the dataset -- i.e. this is relative to however many samples
# have STAGE_21 output at the time you run the script, not a fixed sample
# target, so which species clear the bar can shift between runs while the
# pipeline is still going. Set to 0 to disable and show every bin again.
MAG_PREVALENCE_THRESHOLD <- 0.30

################################################################################
# TAXON-FOCUS TARGETS
################################################################################
# One stacked-bar plot per entry, subsetting STAGE_21's per-bin output to
# only the bin(s) identified as that taxon, per sample.
#
# Matching is against GTDB-Tk species names as they appear in STAGE_21's
# manifest (bin_humann_manifest.tsv), which come straight from the bin
# filenames -- i.e. whatever your GTDB-Tk release actually calls things.
# GTDB SPLITS/RENAMES SPECIES BETWEEN RELEASES, AND WE CANNOT KNOW WHAT
# EXACT STRING YOUR INSTALL PRODUCES FROM HERE. Before trusting these plots,
# check what's actually in your manifests, e.g.:
#   grep -h "Bifidobacterium" humann_contigs_output/*/*_bin_humann_manifest.tsv
#   grep -h "Escherichia"     humann_contigs_output/*/*_bin_humann_manifest.tsv
# and adjust the gtdbtk_regex values below to match. The Bifidobacterium
# infantis pattern below is a best guess (either a validly-published
# "Bifidobacterium infantis" species name, or a GTDB placeholder suffix like
# "Bifidobacterium longum_A") and should be verified against your own data --
# don't take it on faith.
#
# metaphlan_fallback / metaphlan_species: for Segatella copri specifically,
# some GTDB-Tk releases still call this genus "Prevotella" and/or lag the
# Segatella split entirely. For bins GTDB-Tk does NOT already match against
# gtdbtk_regex, this also checks STAGE_12's MetaPhlAn crosscheck summary
# (${CROSSCHECK_ROOT}/<sample>_MAG_metaphlan_crosscheck.tsv) for a bin whose
# metaphlan_top_hit_species equals metaphlan_species, and includes it too --
# so a bin only shows up via the MetaPhlAn fallback when GTDB-Tk didn't
# already identify it directly. Requires STAGE_12 to have been run.
TARGET_TAXA <- list(
    "Bifidobacterium longum" = list(
        gtdbtk_regex       = "(?i)^Bifidobacterium longum$",
        metaphlan_fallback = FALSE
    ),
    "Bifidobacterium infantis" = list(
        gtdbtk_regex       = "(?i)^Bifidobacterium (infantis|longum_A)$",
        metaphlan_fallback = FALSE
    ),
    "Escherichia coli" = list(
        gtdbtk_regex       = "(?i)^Escherichia coli$",
        metaphlan_fallback = FALSE
    ),
    "Segatella copri" = list(
        gtdbtk_regex       = "(?i)^(Segatella|Prevotella) copri$",
        metaphlan_fallback = TRUE,
        metaphlan_species  = "Segatella copri"
    )
)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Shown on every plot so the exclusion is visible on the image itself, not
# just in code -- bars will not sum to 100% while this is non-empty.
EXCLUSION_CAPTION <- if (length(EXCLUDE_FEATURES) > 0) {
    paste0("Excludes: ", paste(EXCLUDE_FEATURES, collapse = ", "), " -- bars do not sum to 100%")
} else {
    NULL
}

################################################################################
# HELPERS
################################################################################

# Read one HUMAnN3 *_normalised.tsv table (2 columns: feature, abundance),
# keep only community-level (unstratified) rows -- i.e. drop any row whose
# feature name contains "|", which marks a species-stratified contribution
# to that pathway/gene family rather than the pathway/gene family itself --
# and drop anything in EXCLUDE_FEATURES (see CONFIG).
read_humann_table <- function(path) {
    df <- suppressMessages(read_tsv(path, col_names = TRUE, show_col_types = FALSE))
    names(df) <- c("feature", "abundance")
    df %>%
        filter(!str_detect(feature, "\\|")) %>%
        filter(!feature %in% EXCLUDE_FEATURES) %>%
        mutate(abundance = as.numeric(abundance))
}

# Shorten a HUMAnN pathway id like "PWY-5484: glycolysis II ..." to a
# legend-friendly label; falls back to the raw feature for UNMAPPED etc.
short_label <- function(feature, width = LEGEND_LABEL_WIDTH) {
    str_trunc(feature, width, ellipsis = "...")
}

# Collapse a long-format (group, feature, abundance) table to "top N + Other".
top_n_plus_other <- function(df, top_n = TOP_N) {
    ranking <- df %>%
        group_by(feature) %>%
        summarise(mean_abundance = mean(abundance), .groups = "drop") %>%
        arrange(desc(mean_abundance))

    keep_features <- head(ranking$feature, top_n)

    df %>%
        mutate(
            plot_feature = ifelse(feature %in% keep_features, feature, "Other (below top N)")
        ) %>%
        group_by(across(-c(feature, abundance)), plot_feature) %>%
        summarise(abundance = sum(abundance), .groups = "drop")
}

# Consistent fill scale: always-there categories get fixed grey/black tones,
# "Other" gets a fixed grey, everything else gets a qualitative palette.
build_fill_scale <- function(levels_present) {
    special <- c(
        "Other (below top N)" = "grey35"
    )
    named <- setdiff(levels_present, names(special))
    pal <- setNames(
        scales::hue_pal()(length(named)),
        named
    )
    all_colors <- c(pal, special[intersect(names(special), levels_present)])
    all_colors[levels_present]
}

################################################################################
# STAGE 20 -- whole-sample, read-based functional profile
################################################################################

collect_stage20 <- function(humann_reads_dir, profile_type) {
    pattern <- paste0("_", profile_type, "_normalised\\.tsv$")
    files <- list.files(humann_reads_dir, pattern = pattern, full.names = TRUE, recursive = TRUE)

    if (length(files) == 0) {
        warning("No STAGE_20 ", profile_type, " tables found under ", humann_reads_dir)
        return(tibble(sample = character(), feature = character(), abundance = double()))
    }

    map_dfr(files, function(f) {
        # out.<SAMPLE>_combined_pathabundance_normalised.tsv -> <SAMPLE>
        sample <- basename(f) %>%
            str_remove("^out\\.") %>%
            str_remove(paste0("_combined_", profile_type, "_normalised\\.tsv$"))

        read_humann_table(f) %>% mutate(sample = sample, .before = 1)
    })
}

plot_stage20 <- function(df, profile_type) {
    if (nrow(df) == 0) return(invisible(NULL))

    plot_df <- top_n_plus_other(df, TOP_N)
    plot_df$label <- short_label(plot_df$plot_feature)
    plot_df$label <- factor(plot_df$label, levels = short_label(unique(plot_df$plot_feature)))

    fill_colors <- build_fill_scale(unique(plot_df$plot_feature))
    names(fill_colors) <- short_label(names(fill_colors))

    p <- ggplot(plot_df, aes(x = sample, y = abundance, fill = label)) +
        geom_col(width = 0.7, colour = "white", linewidth = 0.15) +
        scale_fill_manual(values = fill_colors, name = str_to_title(profile_type)) +
        scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
        labs(
            title = "Whole-sample functional abundance (STAGE_20, reads mode)",
            x = "Sample", y = "Relative abundance",
            caption = EXCLUSION_CAPTION
        ) +
        theme_minimal(base_size = 12) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.text = element_text(size = 8),
            legend.key.size = unit(0.4, "cm"),
            panel.grid.major.x = element_blank()
        )

    out_file <- file.path(OUTPUT_DIR, paste0("stage20_", profile_type, "_stacked_bar.png"))
    ggsave(out_file, p, width = max(6, 0.6 * length(unique(plot_df$sample)) + 4), height = 7, dpi = 300)
    message("Wrote ", out_file)
    invisible(p)
}

################################################################################
# STAGE 21 -- per-MAG (bin) functional profile
################################################################################

collect_stage21 <- function(humann_contigs_dir, profile_type) {
    sample_dirs <- list.dirs(humann_contigs_dir, full.names = TRUE, recursive = FALSE)

    if (length(sample_dirs) == 0) {
        warning("No sample subdirectories found under ", humann_contigs_dir)
        return(tibble(sample = character(), bin = character(), feature = character(),
                      abundance = double(), gtdbtk_species = character(), bin_label = character()))
    }

    map_dfr(sample_dirs, function(sdir) {
        sample <- basename(sdir)

        manifest_file <- list.files(sdir, pattern = "_bin_humann_manifest\\.tsv$", full.names = TRUE)
        manifest <- if (length(manifest_file) == 1) {
            suppressMessages(read_tsv(
                manifest_file,
                col_types = cols(
                    bin            = col_character(),
                    gtdbtk_species = col_character(),
                    n_reads        = col_double(),
                    status         = col_character()
                )
            )) %>%
                filter(status == "ok") %>%
                select(bin, gtdbtk_species)
        } else {
            tibble(bin = character(), gtdbtk_species = character())
        }

        pattern <- paste0("_", profile_type, "_normalised\\.tsv$")
        files <- list.files(sdir, pattern = pattern, full.names = TRUE)
        # The manifest file itself never matches this pattern, but guard anyway
        files <- files[!str_detect(files, "_bin_humann_manifest")]

        if (length(files) == 0) {
            return(tibble(sample = character(), bin = character(), feature = character(),
                          abundance = double(), gtdbtk_species = character(), bin_label = character()))
        }

        map_dfr(files, function(f) {
            bin_label <- basename(f) %>% str_remove(paste0("_", profile_type, "_normalised\\.tsv$"))
            read_humann_table(f) %>%
                mutate(sample = sample, bin = bin_label, .before = 1)
        }) %>%
            left_join(manifest, by = "bin") %>%
            mutate(
                gtdbtk_species = ifelse(is.na(gtdbtk_species), bin, gtdbtk_species),
                bin_label = paste0(bin, "\n(", gtdbtk_species, ")")
            )
    })
}

# Keep only bins whose GTDB-Tk species appears in at least `threshold` of the
# samples PRESENT IN df -- i.e. samples with STAGE_21 output already, not
# necessarily every sample the run will eventually cover. Prevalence is
# computed per species (one count per sample a species appears in at all,
# regardless of how many bins of it that sample has), not per bin.
filter_prevalent_species <- function(df, threshold = MAG_PREVALENCE_THRESHOLD) {
    if (threshold <= 0) return(df)

    total_samples <- n_distinct(df$sample)

    prevalence <- df %>%
        distinct(sample, gtdbtk_species) %>%
        count(gtdbtk_species, name = "n_samples") %>%
        mutate(prevalence = n_samples / total_samples)

    keep_species <- prevalence %>% filter(prevalence >= threshold) %>% pull(gtdbtk_species)

    message("  MAG prevalence filter: keeping ", length(keep_species), "/", nrow(prevalence),
            " species present in >=", scales::percent(threshold), " of ", total_samples,
            " sample(s) currently collected (>= ", ceiling(threshold * total_samples), " sample(s)).")

    df %>% filter(gtdbtk_species %in% keep_species)
}

plot_stage21 <- function(df, profile_type) {
    if (nrow(df) == 0) return(invisible(NULL))

    df <- filter_prevalent_species(df, MAG_PREVALENCE_THRESHOLD)

    if (nrow(df) == 0) {
        message("  No species cleared the MAG_PREVALENCE_THRESHOLD -- skipping STAGE_21 plot. ",
                "Lower MAG_PREVALENCE_THRESHOLD if this happens with real data (e.g. early in the run, ",
                "with few samples collected so far).")
        return(invisible(NULL))
    }

    # Rank top-N pathways ONCE across the whole pooled dataset (all bins,
    # all samples) so the same feature always maps to the same fill colour
    # in every facet -- ranking per-bin would give each facet a different
    # top-N set and make the legend meaningless across bins.
    plot_df <- top_n_plus_other(df, TOP_N)

    plot_df$label <- short_label(plot_df$plot_feature)
    plot_df$label <- factor(plot_df$label, levels = short_label(unique(plot_df$plot_feature)))

    fill_colors <- build_fill_scale(unique(plot_df$plot_feature))
    names(fill_colors) <- short_label(names(fill_colors))

    n_samples <- length(unique(plot_df$sample))

    prevalence_caption <- if (MAG_PREVALENCE_THRESHOLD > 0) {
        paste0("Species present in <", scales::percent(MAG_PREVALENCE_THRESHOLD),
               " of collected samples are omitted")
    } else {
        NULL
    }
    combined_caption <- paste(c(EXCLUSION_CAPTION, prevalence_caption), collapse = "\n")
    if (combined_caption == "") combined_caption <- NULL

    p <- ggplot(plot_df, aes(x = bin_label, y = abundance, fill = label)) +
        geom_col(width = 0.7, colour = "white", linewidth = 0.15) +
        scale_fill_manual(values = fill_colors, name = str_to_title(profile_type)) +
        scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
        facet_wrap(~ sample, scales = "free_x", ncol = min(3, n_samples)) +
        labs(
            title = "Per-MAG functional abundance (STAGE_21, bin-mapped reads)",
            x = "Bin (GTDB-Tk species)", y = "Relative abundance",
            caption = combined_caption
        ) +
        theme_minimal(base_size = 11) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
            legend.text = element_text(size = 8),
            legend.key.size = unit(0.4, "cm"),
            panel.grid.major.x = element_blank(),
            strip.text = element_text(face = "bold")
        )

    out_file <- file.path(OUTPUT_DIR, paste0("stage21_", profile_type, "_stacked_bar_per_bin.png"))
    n_bins <- length(unique(plot_df$bin_label))
    ggsave(out_file, p, width = max(8, 0.5 * n_bins + 4), height = 4 * ceiling(n_samples / min(3, n_samples)) + 2, dpi = 300, limitsize = FALSE)
    message("Wrote ", out_file)
    invisible(p)
}

################################################################################
# TAXON FOCUS -- subset STAGE_21's per-bin output to specific species
################################################################################

# Load STAGE_12's whole-sample cross-check summary tables
# (${CROSSCHECK_ROOT}/<sample>_MAG_metaphlan_crosscheck.tsv), one row per bin,
# with both the GTDB-Tk call and MetaPhlAn's independent top hit for that
# bin's reads side by side. Only needed for the MetaPhlAn-fallback taxa.
load_crosscheck_summary <- function(crosscheck_root) {
    files <- list.files(crosscheck_root, pattern = "_MAG_metaphlan_crosscheck\\.tsv$", full.names = TRUE)

    if (length(files) == 0) {
        warning("No STAGE_12 crosscheck summary tables found under ", crosscheck_root,
                " -- MetaPhlAn-fallback taxon matching will be skipped.")
        return(tibble(sample = character(), bin = character(), gtdbtk_species = character(),
                      n_reads = integer(), metaphlan_top_hit_species = character(),
                      metaphlan_top_hit_pct = double(), genus_match = character(),
                      agreement_status = character()))
    }

    map_dfr(files, function(f) {
        sample <- basename(f) %>% str_remove("_MAG_metaphlan_crosscheck\\.tsv$")
        # Explicit col_types, not show_col_types=FALSE-with-guessing: readr infers
        # types per file from a sample of that file's own rows, and a file where
        # every row happens to look like an integer (or is header-only) can get
        # guessed differently than one with more varied values. Across ~dozens of
        # per-sample files that inconsistency is exactly what broke bind_rows here
        # -- pin every column explicitly so every file parses identically.
        suppressMessages(read_tsv(
            f,
            col_types = cols(
                bin                       = col_character(),
                gtdbtk_species            = col_character(),
                n_reads                   = col_double(),
                metaphlan_top_hit_species = col_character(),
                metaphlan_top_hit_pct     = col_double(),
                genus_match               = col_character(),
                agreement_status          = col_character()
            )
        )) %>%
            mutate(sample = sample, .before = 1)
    })
}

# For each target taxon: match bins by GTDB-Tk species (from STAGE_21's own
# manifest, embedded in stage21_df via bin_label/gtdbtk_species), and -- only
# for taxa with metaphlan_fallback = TRUE -- also pull in any bin GTDB-Tk did
# NOT already match, but that STAGE_12's MetaPhlAn crosscheck called as
# metaphlan_species for that taxon. GTDB-Tk match always takes priority where
# both would apply.
identify_target_bins <- function(stage21_df, crosscheck_df, target_taxa) {
    bin_species <- stage21_df %>% distinct(sample, bin, gtdbtk_species)

    map_dfr(names(target_taxa), function(taxon) {
        cfg <- target_taxa[[taxon]]

        gtdb_hits <- bin_species %>%
            filter(str_detect(str_squish(gtdbtk_species), regex(cfg$gtdbtk_regex))) %>%
            transmute(sample, bin, target_taxon = taxon, match_source = "gtdbtk",
                      matched_label = gtdbtk_species)

        mp_hits <- if (isTRUE(cfg$metaphlan_fallback) && nrow(crosscheck_df) > 0) {
            crosscheck_df %>%
                filter(str_to_lower(str_squish(metaphlan_top_hit_species)) == str_to_lower(cfg$metaphlan_species)) %>%
                transmute(sample, bin, target_taxon = taxon, match_source = "metaphlan_fallback",
                          matched_label = metaphlan_top_hit_species)
        } else {
            tibble(sample = character(), bin = character(), target_taxon = character(),
                   match_source = character(), matched_label = character())
        }

        # gtdb_hits first so distinct() keeps the GTDB-Tk match over the
        # fallback whenever a bin was identified both ways.
        bind_rows(gtdb_hits, mp_hits) %>%
            distinct(sample, bin, target_taxon, .keep_all = TRUE)
    })
}

plot_taxon_focus <- function(stage21_df, taxon_matches, taxon, profile_type) {
    bins_for_taxon <- taxon_matches %>% filter(target_taxon == taxon)

    if (nrow(bins_for_taxon) == 0) {
        message("  No bins matched '", taxon, "' in any sample -- skipping.")
        return(invisible(NULL))
    }

    fallback_samples <- bins_for_taxon %>% filter(match_source == "metaphlan_fallback") %>% pull(sample)
    if (length(fallback_samples) > 0) {
        message("  '", taxon, "' identified via STAGE_12 MetaPhlAn fallback (not GTDB-Tk) in: ",
                paste(fallback_samples, collapse = ", "))
    }

    df <- stage21_df %>%
        inner_join(bins_for_taxon %>% select(sample, bin, match_source), by = c("sample", "bin"))

    if (nrow(df) == 0) {
        message("  '", taxon, "' matched ", nrow(bins_for_taxon), " bin(s) via crosscheck, but none have ",
                "STAGE_21 pathway data yet (job likely still running/pending) -- skipping for now.")
        return(invisible(NULL))
    }

    plot_df <- top_n_plus_other(df %>% select(-gtdbtk_species, -bin_label), TOP_N)
    plot_df$label <- short_label(plot_df$plot_feature)
    plot_df$label <- factor(plot_df$label, levels = short_label(unique(plot_df$plot_feature)))

    fill_colors <- build_fill_scale(unique(plot_df$plot_feature))
    names(fill_colors) <- short_label(names(fill_colors))

    # Flag which samples used the MetaPhlAn fallback directly on the x-axis
    plot_df <- plot_df %>%
        mutate(sample_label = ifelse(match_source == "metaphlan_fallback",
                                      paste0(sample, "*"), sample))

    p <- ggplot(plot_df, aes(x = sample_label, y = abundance, fill = label)) +
        geom_col(width = 0.6, colour = "white", linewidth = 0.15) +
        scale_fill_manual(values = fill_colors, name = str_to_title(profile_type)) +
        scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
        labs(
            title = paste0("Functional abundance: ", taxon, " (per-sample MAG)"),
            subtitle = if (length(fallback_samples) > 0) {
                "* = identified via STAGE_12 MetaPhlAn (not GTDB-Tk) for this sample"
            } else NULL,
            x = "Sample", y = "Relative abundance",
            caption = EXCLUSION_CAPTION
        ) +
        theme_minimal(base_size = 12) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.text = element_text(size = 8),
            legend.key.size = unit(0.4, "cm"),
            panel.grid.major.x = element_blank()
        )

    slug <- taxon %>% str_replace_all("[^A-Za-z0-9]+", "_") %>% str_to_lower()
    out_file <- file.path(OUTPUT_DIR, paste0("taxon_", slug, "_", profile_type, "_stacked_bar.png"))
    ggsave(out_file, p, width = max(6, 0.6 * length(unique(plot_df$sample_label)) + 4), height = 7, dpi = 300)
    message("  Wrote ", out_file)
    invisible(p)
}

################################################################################
# RUN
################################################################################

message("Collecting STAGE_20 (whole-sample) ", PROFILE_TYPE, " tables...")
stage20_df <- collect_stage20(HUMANN_READS_ROOT, PROFILE_TYPE)
plot_stage20(stage20_df, PROFILE_TYPE)

message("Collecting STAGE_21 (per-MAG) ", PROFILE_TYPE, " tables...")
stage21_df <- collect_stage21(HUMANN_CONTIGS_ROOT, PROFILE_TYPE)
plot_stage21(stage21_df, PROFILE_TYPE)

message("Matching bins to taxon-focus targets...")
crosscheck_df <- load_crosscheck_summary(CROSSCHECK_ROOT)
taxon_matches <- identify_target_bins(stage21_df, crosscheck_df, TARGET_TAXA)

message("Building taxon-focus plots...")
for (taxon in names(TARGET_TAXA)) {
    plot_taxon_focus(stage21_df, taxon_matches, taxon, PROFILE_TYPE)
}

if (nrow(taxon_matches) > 0) {
    matches_file <- file.path(OUTPUT_DIR, "taxon_focus_bin_matches.tsv")
    write_tsv(taxon_matches, matches_file)
    message("Wrote ", matches_file, " (audit trail of which bin matched which taxon, and how)")
}

message("Done. Plots written to ", OUTPUT_DIR)