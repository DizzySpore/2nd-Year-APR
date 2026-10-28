#!/bin/bash

################################################################################
# STAGE 24: R analyses - MaAsLin3 statistical analysis
################################################################################
#
# Input variables (from SLURM template):
#   - SAMPLE      : Sample name
#   - NUM_THREADS : Number of CPU threads
#   - DIRECTORY   : Working directory (contains all sample outputs)
#   - STATUS_DIR  : Path to status directory
#
# Output:
#   - MaAsLin3 Bracken results  : ${DIRECTORY}/maaslin3_output_bracken/
#   - MaAsLin3 pathway results  : ${DIRECTORY}/maaslin3_output_pathways/
#   - ${STATUS_DIR}/${SAMPLE}.resume cleared and .done written on completion
#
# NOTE: This stage is a placeholder and will be updated to reflect changes
#       in study design, metadata structure, and sample ID patterns.
#       Current implementation is based on MaAsLin3_Analysis_Scriptv0_5.r
#
# Dependencies:
#   - R with maaslin3 installed at:
#     <R_LIB_PATH>
#   - Metadata file at: ${DIRECTORY}/metadata.csv
#   - Bracken outputs in: ${DIRECTORY}/kraken2_bracken_out/
#   - HUMAnN3 outputs in: ${DIRECTORY}/humann_reads_output/
#
################################################################################

R_LIB="${R_LIB:-<R_LIB_PATH>}"
METADATA="${DIRECTORY}/metadata.csv"
BRACKEN_DIR="${DIRECTORY}/kraken2_bracken_out"
HUMANN_DIR="${DIRECTORY}/humann_reads_output"
MAASLIN_BRACKEN_OUT="${DIRECTORY}/maaslin3_output_bracken"
MAASLIN_PATHWAY_OUT="${DIRECTORY}/maaslin3_output_pathways"

echo "$(date): [STAGE 24/24] Starting MaAsLin3 statistical analysis..."
echo "  Metadata    : ${METADATA}"
echo "  Bracken dir : ${BRACKEN_DIR}"
echo "  HUMAnN3 dir : ${HUMANN_DIR}"
echo "  Output      : ${DIRECTORY}"

echo "  R library   : ${R_LIB}"

# Verify metadata file exists
if [ ! -f "${METADATA}" ]; then
    echo "$(date): ERROR - Metadata file not found: ${METADATA}"
    echo "$(date): Please provide metadata.csv in the working directory."
    return 1
fi

mkdir -p "${MAASLIN_BRACKEN_OUT}" "${MAASLIN_PATHWAY_OUT}"

conda activate r-env

Rscript - << REOF
# -- Libraries -----------------------------------------------------------------
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
    library(maaslin3, lib.loc = "${R_LIB}")
})

# -- Paths ---------------------------------------------------------------------
directory     <- "${DIRECTORY}"
bracken_dir   <- "${BRACKEN_DIR}"
humann_dir    <- "${HUMANN_DIR}"
metadata_file <- "${METADATA}"
bracken_out   <- "${MAASLIN_BRACKEN_OUT}"
pathway_out   <- "${MAASLIN_PATHWAY_OUT}"
num_threads   <- as.integer("${NUM_THREADS}")

cat(sprintf("[%s] Loading metadata...\n", Sys.time()))

# -- Step 1: Load metadata -----------------------------------------------------
# TODO: update metagenome_id column name if metadata structure changes
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

# -- Step 2: Load and combine Bracken abundance files -------------------------
cat(sprintf("[%s] Loading Bracken abundance files...\n", Sys.time()))

bracken_files <- list.files(
    path       = bracken_dir,
    pattern    = "_prok-bracken_abundance.txt",
    full.names = TRUE
)

if (length(bracken_files) == 0) {
    stop("No Bracken abundance files found in: ", bracken_dir)
}

bracken_data_list <- lapply(bracken_files, function(file) {
    # TODO: update sample ID extraction pattern to match actual filenames
    sample_id <- sub(".*/(.*?)_prok-bracken_abundance.txt", "\\1", file)
    tryCatch({
        df <- read_tsv(file,
            col_types = cols(
                name                 = col_character(),
                taxonomy_id          = col_character(),
                taxonomy_lvl         = col_character(),
                kraken_assigned_reads = col_double(),
                added_reads          = col_double(),
                new_est_reads        = col_double(),
                fraction_total_reads = col_double()
            ),
            quote = "")

        if (!all(c("name", "fraction_total_reads") %in% colnames(df))) {
            warning(paste("File", file, "missing required columns"))
            return(NULL)
        }

        df %>%
            select(name, fraction_total_reads) %>%
            mutate(
                fraction_total_reads = as.numeric(fraction_total_reads),
                sample_id = sample_id
            )
    }, error = function(e) {
        warning(paste("Error reading file", file, ":", e$message))
        return(NULL)
    })
})

bracken_data_list <- bracken_data_list[!sapply(bracken_data_list, is.null)]

if (length(bracken_data_list) == 0) {
    stop("No Bracken files were successfully read.")
}

bracken_data <- bind_rows(bracken_data_list)

# Handle duplicates
duplicate_check <- bracken_data %>%
    group_by(name, sample_id) %>%
    summarise(count = n(), .groups = "drop") %>%
    filter(count > 1)

if (nrow(duplicate_check) > 0) {
    warning("Duplicate name/sample_id combinations found. Aggregating by sum.")
    bracken_data <- bracken_data %>%
        group_by(name, sample_id) %>%
        summarise(fraction_total_reads = sum(fraction_total_reads, na.rm = TRUE),
                  .groups = "drop")
}

# Pivot to abundance matrix
bracken_abundance <- bracken_data %>%
    pivot_wider(
        names_from  = sample_id,
        values_from = fraction_total_reads,
        values_fill = 0
    ) %>%
    as.data.frame()

rownames(bracken_abundance) <- bracken_abundance$name
bracken_abundance$name <- NULL

cat(sprintf("[%s] Bracken: %d taxa across %d samples\n",
    Sys.time(), nrow(bracken_abundance), ncol(bracken_abundance)))

# -- Step 3: Load HUMAnN3 pathway abundance files ------------------------------
cat(sprintf("[%s] Loading HUMAnN3 pathway files...\n", Sys.time()))

pathabundance_files <- list.files(
    path      = humann_dir,
    pattern   = "_combined_pathabundance.tsv",
    full.names = TRUE,
    recursive  = TRUE
)

pathabundance_data_list <- lapply(pathabundance_files, function(file) {
    sample_id <- sub(".*out\\.(.*?)_combined_pathabundance\\.tsv", "\\1", file)
    tryCatch({
        df <- read_tsv(file,
            col_names = c("pathway", "abundance"),
            col_types = cols(.default = col_character()),
            quote     = "")
        df$sample_id <- sample_id
        return(df)
    }, error = function(e) {
        warning(paste("Error reading file", file, ":", e$message))
        return(NULL)
    })
})

pathabundance_data_list <- pathabundance_data_list[
    !sapply(pathabundance_data_list, is.null)]

pathabundance_available <- length(pathabundance_data_list) > 0

if (pathabundance_available) {
    pathabundance_abundance <- bind_rows(pathabundance_data_list) %>%
        mutate(abundance = as.numeric(abundance)) %>%
        pivot_wider(
            names_from  = sample_id,
            values_from = abundance,
            values_fill = 0
        ) %>%
        as.data.frame()

    rownames(pathabundance_abundance) <- pathabundance_abundance$pathway
    pathabundance_abundance$pathway <- NULL

    cat(sprintf("[%s] HUMAnN3: %d pathways across %d samples\n",
        Sys.time(), nrow(pathabundance_abundance),
        ncol(pathabundance_abundance)))
} else {
    warning("No HUMAnN3 pathway files found. Skipping pathway analysis.")
}

# -- Step 4: Align metadata and abundance data ---------------------------------
cat(sprintf("[%s] Aligning metadata and abundance data...\n", Sys.time()))

# TODO: update metagenome_id to match the sample ID column in your metadata
metadata <- metadata %>%
    filter(metagenome_id %in% colnames(bracken_abundance)) %>%
    as.data.frame() %>%
    `rownames<-`(.$metagenome_id) %>%
    select(-metagenome_id)

common_samples <- intersect(colnames(bracken_abundance), rownames(metadata))

if (pathabundance_available) {
    common_samples <- intersect(common_samples, colnames(pathabundance_abundance))
}

if (length(common_samples) == 0) {
    stop("No common samples across datasets. Check sample IDs.")
}

cat(sprintf("[%s] Common samples: %d\n", Sys.time(), length(common_samples)))

bracken_abundance <- bracken_abundance[, common_samples, drop = FALSE]
metadata          <- metadata[common_samples, , drop = FALSE]
if (pathabundance_available) {
    pathabundance_abundance <- pathabundance_abundance[,
        common_samples, drop = FALSE]
}

# -- Step 5: Preprocess metadata columns ---------------------------------------
# TODO: update column list to match actual metadata structure
numeric_cols     <- c("muac", "age", "hemocue", "height", "weight")
categorical_cols <- c("gender", "mrdt", "diarrhoea", "excl_brstfed")

for (col in numeric_cols) {
    if (col %in% colnames(metadata)) {
        metadata[[col]] <- as.numeric(as.character(metadata[[col]]))
    } else {
        warning(paste(col, "column not found in metadata."))
    }
}

for (col in categorical_cols) {
    if (col %in% colnames(metadata)) {
        invalid <- sum(!is.na(metadata[[col]]) &
                       !metadata[[col]] %in% c("0", "1"))
        if (invalid > 0) {
            warning(paste(invalid, "invalid values in", col, "- converting to NA"))
            metadata[[col]][!metadata[[col]] %in% c("0", "1")] <- NA
        }
    } else {
        warning(paste(col, "column not found in metadata."))
    }
}

if ("lab_id" %in% colnames(metadata)) {
    metadata$lab_id <- as.factor(metadata$lab_id)
}

# Remove samples with missing key variables
key_columns <- intersect(c("muac", "age", "gender", "hemocue",
                            "height", "weight", "mrdt",
                            "diarrhoea", "excl_brstfed"),
                          colnames(metadata))
metadata <- metadata[complete.cases(metadata[, key_columns]), , drop = FALSE]
bracken_abundance <- bracken_abundance[, rownames(metadata), drop = FALSE]
if (pathabundance_available) {
    pathabundance_abundance <- pathabundance_abundance[,
        rownames(metadata), drop = FALSE]
}

cat(sprintf("[%s] Samples after metadata filtering: %d\n",
    Sys.time(), nrow(metadata)))

# -- Step 6: Filter low-abundance and low-prevalence features ------------------
min_prevalence <- 0.1 * ncol(bracken_abundance)
min_abundance  <- 0.001

bracken_abundance <- bracken_abundance[
    rowSums(bracken_abundance > 0) >= min_prevalence &
    rowMeans(bracken_abundance) >= min_abundance, , drop = FALSE]

cat(sprintf("[%s] Filtered Bracken taxa: %d\n",
    Sys.time(), nrow(bracken_abundance)))

# -- Step 7: Build fixed effects from available metadata columns ---------------
# TODO: update fixed_effects to match study design
fixed_effects <- intersect(
    c("muac", "age", "gender", "hemocue", "height",
      "weight", "mrdt", "diarrhoea", "excl_brstfed"),
    colnames(metadata)
)

cat(sprintf("[%s] Fixed effects: %s\n",
    Sys.time(), paste(fixed_effects, collapse = ", ")))

# -- Step 8: Run MaAsLin3 on Bracken data -------------------------------------
cat(sprintf("[%s] Running MaAsLin3 on Bracken abundance...\n", Sys.time()))

maaslin_results_bracken <- maaslin3(
    input_data     = bracken_abundance,
    input_metadata = metadata,
    output         = bracken_out,
    fixed_effects  = fixed_effects,
    normalization  = "TSS",
    transform      = "LOG",
    min_abundance  = 0,
    min_prevalence = 0,
    cores          = num_threads
)

write.csv(maaslin_results_bracken$results,
    file.path(bracken_out, "significant_results.csv"))

cat(sprintf("[%s] MaAsLin3 Bracken analysis complete.\n", Sys.time()))

# -- Step 9: Run MaAsLin3 on HUMAnN3 pathway data -----------------------------
if (pathabundance_available && ncol(pathabundance_abundance) > 0) {
    pathabundance_abundance <- pathabundance_abundance[
        rowSums(pathabundance_abundance > 0) >= min_prevalence &
        rowMeans(pathabundance_abundance) >= min_abundance, , drop = FALSE]

    cat(sprintf("[%s] Filtered pathways: %d\n",
        Sys.time(), nrow(pathabundance_abundance)))
    cat(sprintf("[%s] Running MaAsLin3 on pathway abundance...\n", Sys.time()))

    maaslin_results_pathways <- maaslin3(
        input_data     = pathabundance_abundance,
        input_metadata = metadata,
        output         = pathway_out,
        fixed_effects  = fixed_effects,
        normalization  = "TSS",
        transform      = "LOG",
        min_abundance  = 0,
        min_prevalence = 0,
        cores          = num_threads
    )

    write.csv(maaslin_results_pathways$results,
        file.path(pathway_out, "significant_results.csv"))

    cat(sprintf("[%s] MaAsLin3 pathway analysis complete.\n", Sys.time()))
} else {
    cat(sprintf("[%s] Skipping pathway analysis - no data available.\n",
        Sys.time()))
}

cat(sprintf("[%s] All R analyses complete.\n", Sys.time()))
REOF

conda deactivate || true

# Mark sample as fully complete
touch "${STATUS_DIR}/${SAMPLE}.done"
rm -f "${STATUS_DIR}/${SAMPLE}.resume"

echo "$(date): [STAGE 24/24] R analyses complete. Pipeline finished for ${SAMPLE}."
