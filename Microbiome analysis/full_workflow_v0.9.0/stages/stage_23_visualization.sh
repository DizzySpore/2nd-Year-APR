#!/bin/bash

################################################################################
# STAGE 23: Visualization - Kraken2/Bracken and abundance plots
################################################################################
#
# Input variables (from SLURM template):
#   - SAMPLE           : Sample name
#   - KRAKEN2_OUT      : Directory containing Kraken2/Bracken outputs
#   - METAPHLAN_OUT    : Directory containing MetaPhlAn outputs
#   - HUMANN_READS_DIR : Directory containing HUMAnN3 reads mode outputs
#   - VISUALIZATION_OUT: Output directory for plots and visualizations
#   - STATUS_DIR       : Path to status directory
#
# Output:
#   - Krona plots      : ${KRAKEN2_OUT}/${SAMPLE}_*_krona.html
#   - Kraken2 plot     : ${VISUALIZATION_OUT}/${SAMPLE}_kraken2_toptaxa.png
#   - MetaPhlAn plot   : ${VISUALIZATION_OUT}/${SAMPLE}_metaphlan_species.png
#   - HUMAnN3 plot     : ${VISUALIZATION_OUT}/${SAMPLE}_humann3_pathways.png
#   - Completion marker: ${STATUS_DIR}/abundance_visualization_complete.txt
#   - ${STATUS_DIR}/${SAMPLE}.resume updated to STAGE_24
#
# Note:
#   - Krona plots are generated from each Bracken abundance file using
#     ktImportTaxonomy. The -t and -m column flags should be verified
#     against the actual Bracken output format on your system.
#   - R visualizations use ggplot2 and are warnings not errors if they
#     fail since missing input files are a valid outcome for some samples
#
################################################################################

echo "$(date): [STAGE 23/24] Starting visualization for ${SAMPLE}..."
echo "  Krona input  : ${KRAKEN2_OUT}/${SAMPLE}*_bracken_abundance.txt"
echo "  R plots out  : ${VISUALIZATION_OUT}"

# -- Krona interactive plots ---------------------------------------------------
conda activate krona-env

for abundance_file in "${KRAKEN2_OUT}/${SAMPLE}"*_bracken_abundance.txt; do
    [ -f "$abundance_file" ] || continue
    BASENAME=$(basename "$abundance_file" _bracken_abundance.txt)

    echo "$(date): Generating Krona plot for ${BASENAME}..."
    ktImportTaxonomy \
        -t 5 -m 3 \
        -o "${KRAKEN2_OUT}/${BASENAME}_krona.html" \
        "$abundance_file" || \
    echo "$(date): WARNING - Krona plot failed for ${BASENAME}, continuing."
done

conda deactivate || true

# -- R-based abundance plots ---------------------------------------------------
conda activate r-viz-env

Rscript - << REOF
suppressPackageStartupMessages({
    library(tidyverse)
    library(ggplot2)
})

sample_name     <- "${SAMPLE}"
kraken_dir      <- "${KRAKEN2_OUT}"
metaphlan_dir   <- "${METAPHLAN_OUT}"
humann_dir      <- "${HUMANN_READS_DIR}"
output_dir      <- "${VISUALIZATION_OUT}"

cat(sprintf("[%s] Generating R visualizations for %s\n", Sys.time(), sample_name))

# -- Kraken2 top taxa bar chart ---------------------------------------------
kraken_report <- file.path(kraken_dir,
    sprintf("%s.prok-kraken2_report.txt", sample_name))

if (file.exists(kraken_report)) {
    kraken_df <- read_tsv(kraken_report,
        col_names = c("abundance", "reads", "reads_assigned", "rank", "taxid", "taxonomy"),
        show_col_types = FALSE)

    top_kraken <- kraken_df %>%
        filter(rank == "S") %>%
        arrange(desc(reads)) %>%
        head(15) %>%
        mutate(taxonomy = str_trim(taxonomy))

    if (nrow(top_kraken) > 0) {
        p_kraken <- ggplot(top_kraken,
                aes(x = reorder(taxonomy, reads), y = reads)) +
            geom_bar(stat = "identity", fill = "steelblue") +
            coord_flip() +
            labs(
                title = sprintf("Top 15 Taxa - %s (Kraken2)", sample_name),
                x = "Taxonomy", y = "Read Count") +
            theme_minimal() +
            theme(plot.title = element_text(size = 12, face = "bold"))

        ggsave(file.path(output_dir,
            sprintf("%s_kraken2_toptaxa.png", sample_name)),
            p_kraken, width = 8, height = 6, dpi = 300)
        cat(sprintf("[%s] Kraken2 top taxa plot saved.\n", Sys.time()))
    }
} else {
    cat(sprintf("[%s] WARNING: Kraken2 report not found, skipping.\n", Sys.time()))
}

# -- MetaPhlAn species bar chart --------------------------------------------
metaphlan_prof <- file.path(metaphlan_dir,
    sprintf("%s_metaphlan_profile.txt", sample_name))

if (file.exists(metaphlan_prof)) {
    # Skip comment lines beginning with # before reading
    metaphlan_df <- read_tsv(metaphlan_prof,
        comment = "#", show_col_types = FALSE) %>%
        filter(str_detect(clade_name, "s__"),
               !str_detect(clade_name, "t__")) %>%
        arrange(desc(relative_abundance)) %>%
        head(15) %>%
        mutate(species = str_extract(clade_name, "s__[^|]+") %>%
                   str_remove("s__") %>%
                   str_replace_all("_", " "))

    if (nrow(metaphlan_df) > 0) {
        p_metaphlan <- ggplot(metaphlan_df,
                aes(x = reorder(species, relative_abundance),
                    y = relative_abundance)) +
            geom_bar(stat = "identity", fill = "coral") +
            coord_flip() +
            labs(
                title = sprintf("Top 15 Species - %s (MetaPhlAn)", sample_name),
                x = "Species", y = "Relative Abundance (%)") +
            theme_minimal() +
            theme(plot.title = element_text(size = 12, face = "bold"),
                  axis.text.y = element_text(face = "italic"))

        ggsave(file.path(output_dir,
            sprintf("%s_metaphlan_species.png", sample_name)),
            p_metaphlan, width = 8, height = 6, dpi = 300)
        cat(sprintf("[%s] MetaPhlAn species plot saved.\n", Sys.time()))
    }
} else {
    cat(sprintf("[%s] WARNING: MetaPhlAn profile not found, skipping.\n", Sys.time()))
}

# -- HUMAnN3 pathway bar chart ----------------------------------------------
humann_path <- file.path(humann_dir,
    sprintf("out.%s_combined_pathabundance_normalised.tsv", sample_name))

if (file.exists(humann_path)) {
    humann_df <- read_tsv(humann_path, show_col_types = FALSE) %>%
        rename(Pathway = 1, Abundance = 2) %>%
        filter(!str_detect(Pathway, "\\|")) %>%
        arrange(desc(Abundance)) %>%
        head(20)

    if (nrow(humann_df) > 0) {
        p_humann <- ggplot(humann_df,
                aes(x = reorder(Pathway, Abundance), y = Abundance)) +
            geom_bar(stat = "identity", fill = "mediumseagreen") +
            coord_flip() +
            labs(
                title = sprintf("Top 20 Pathways - %s (HUMAnN3)", sample_name),
                x = "Pathway", y = "Relative Abundance") +
            theme_minimal() +
            theme(plot.title = element_text(size = 12, face = "bold"),
                  axis.text.y = element_text(size = 7))

        ggsave(file.path(output_dir,
            sprintf("%s_humann3_pathways.png", sample_name)),
            p_humann, width = 10, height = 8, dpi = 300)
        cat(sprintf("[%s] HUMAnN3 pathway plot saved.\n", Sys.time()))
    }
} else {
    cat(sprintf("[%s] WARNING: HUMAnN3 pathway file not found, skipping.\n", Sys.time()))
}

cat(sprintf("[%s] R visualization complete for %s.\n", Sys.time(), sample_name))
REOF

conda deactivate || true

# Write shared completion marker used by checkpoint scanner
touch "${STATUS_DIR}/abundance_visualization_complete.txt"

echo "STAGE_24" > "${STATUS_DIR}/${SAMPLE}.resume"
echo "$(date): [STAGE 23/24] Visualization complete."