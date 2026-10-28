# ==============================================================================
# UNIFIED MICROBIOME ENGINE - COMPLETE CLINICAL, ECONOMIC & DEMOGRAPHIC ENGINE
# ==============================================================================
options(bitmapType = "cairo")

library(mia)
library(scater)
library(ggplot2)
library(ggrepel)
library(scales)
library(scatterplot3d)
library(vegan)
library(dplyr)   # FIX: was used via dplyr::case_when() in plot_permanova_results
                  #      but never loaded explicitly anywhere in the script.

merged_file   <- "metaphlan_wreads_output/merged_abundance_wreads.tsv"
metadata_file <- "metadata_merged.csv"
output_dir    <- "diversity_analysis_outputs"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================
find_col <- function(patterns, df) {
  for (p in patterns) {
    m <- colnames(df)[grep(p, colnames(df), ignore.case = TRUE)]
    if (length(m) > 0) return(m[1])
  }
  return(NULL)
}

recode_factor <- function(raw_num, mapping, na_codes = c(88, 99)) {
  raw_num[raw_num %in% na_codes] <- NA
  out <- rep(NA_character_, length(raw_num))
  for (val in names(mapping)) out[raw_num == as.numeric(val)] <- mapping[[val]]
  factor(out, levels = unname(mapping))
}

recode_yesno <- function(raw_num, na_codes = c(88, 99, 9)) {
  raw_num[raw_num %in% na_codes] <- NA
  factor(ifelse(raw_num == 1, "Yes", ifelse(raw_num == 0, "No", NA)),
         levels = c("No", "Yes"))
}

clean_numeric <- function(x, min_val = -Inf, max_val = Inf,
                          na_codes = c(88, 99, 999, 9999)) {
  x <- as.numeric(x)
  x[x %in% na_codes] <- NA
  x[x < min_val | x > max_val] <- NA
  x
}

stratify_z <- function(col_name, df) {
  if (is.null(col_name)) return(data.frame(raw    = rep(NA_real_,      nrow(df)),
                                           factor = rep(NA_character_, nrow(df))))
  raw <- as.numeric(df[[col_name]])
  raw[raw > 5 | raw < -5] <- NA
  grp <- rep(NA_character_, length(raw))
  grp[raw < -2]             <- "Malnourished"
  grp[raw >= -2 & raw <= 2] <- "Normal"
  grp[raw > 2]              <- "Overweight"
  data.frame(raw = raw, factor = grp, stringsAsFactors = FALSE)
}

recode_water <- function(raw_num) {
  grp <- rep(NA_character_, length(raw_num))
  grp[raw_num %in% c(1,2,3)]  <- "Piped/Tap"
  grp[raw_num %in% c(4,5,6)]  <- "Borehole"
  grp[raw_num == 7]            <- "Bottled"
  grp[raw_num == 8]            <- "Protected well"
  grp[raw_num %in% c(9,10)]   <- "Unprotected/Surface"
  grp[raw_num == 99]           <- NA
  factor(grp, levels = c("Piped/Tap","Borehole","Bottled",
                         "Protected well","Unprotected/Surface"))
}

sig_label <- function(p) {
  ifelse(p <= 0.001, "***",
         ifelse(p <= 0.01,  "**",
                ifelse(p <= 0.05, "*", "ns")))
}

# Qualitative magnitude labels for non-parametric effect sizes.
# Thresholds:
#  - Rank-biserial correlation |r| (Wilcoxon):  Cohen-style convention for r
#  - Epsilon-squared (Kruskal-Wallis):          Rea & Parker (2014) convention
# Qualitative magnitude labels for non-parametric effect sizes.
# Thresholds:
#  - Rank-biserial correlation |r| (Wilcoxon):  Cohen-style convention for r
#  - Epsilon-squared (Kruskal-Wallis):          Rea & Parker (2014) convention
#
# NOTE: effect_size and method are vectors (one entry per row of test_alpha()'s
# output, i.e. one per alpha metric). Different metrics for the SAME category
# can legitimately end up using different tests (NAs differ metric-to-metric,
# which can change how many groups survive the min-count-2 filter), so this
# function is written to be fully vectorised rather than branching with a
# scalar if() on `method`, which errors ("the condition has length > 1") the
# moment a category mixes Wilcoxon and Kruskal-Wallis rows.
effect_label <- function(effect_size, method) {
  a   <- abs(effect_size)
  out <- rep(NA_character_, length(effect_size))

  is_rb <- !is.na(method) & method == "Rank-biserial correlation"
  out[is_rb] <- ifelse(is.na(a[is_rb]), NA_character_,
                       ifelse(a[is_rb] < 0.10, "negligible",
                              ifelse(a[is_rb] < 0.30, "small",
                                     ifelse(a[is_rb] < 0.50, "medium", "large"))))

  is_eps <- !is.na(method) & method == "Epsilon-squared"
  out[is_eps] <- ifelse(is.na(a[is_eps]), NA_character_,
                        ifelse(a[is_eps] < 0.01, "negligible",
                               ifelse(a[is_eps] < 0.08, "small",
                                      ifelse(a[is_eps] < 0.26, "medium", "large"))))

  out
}

# Compute alpha diversity only for indices not already present in colData
ensure_alpha <- function(tse_obj) {
  alpha_indices <- c("Observed","shannon","inverse_simpson","gini_simpson","chao1")
  existing      <- colnames(colData(tse_obj))
  for (idx in alpha_indices) {
    if (!idx %in% existing) {
      tse_obj <- addAlpha(tse_obj, assay.type = "counts", index = idx, name = idx)
    }
  }
  tse_obj
}

# ==============================================================================
# ALPHA DIVERSITY STATISTICAL TESTING (Non-parametric)
# ------------------------------------------------------------------------------
# 2-level categorical variables  -> Wilcoxon rank-sum test (wilcox.test)
# 3+-level categorical variables -> Kruskal-Wallis test (kruskal.test)
# Both are distribution-free tests appropriate for alpha diversity indices,
# which are rarely normally distributed and often have unequal variances
# across small/uneven group sizes.
# ==============================================================================
test_alpha <- function(tse_obj, category) {
  metrics <- c("Observed", "shannon", "inverse_simpson", "gini_simpson", "chao1")
  results <- list()

  for (metric in metrics) {
    if (!metric %in% colnames(colData(tse_obj))) next

    alpha_vals <- colData(tse_obj)[[metric]]
    group_var  <- colData(tse_obj)[[category]]

    valid <- !is.na(alpha_vals) & !is.na(group_var) &
      as.character(group_var) != "" & as.character(group_var) != "NA"
    if (sum(valid) < 6) next

    alpha_vals <- as.numeric(alpha_vals[valid])
    group_var  <- droplevels(factor(as.character(group_var[valid])))

    n_levels <- nlevels(group_var)
    if (n_levels < 2) next

    # Drop levels with fewer than 2 observations - neither test can use them
    tab <- table(group_var)
    keep_levels <- names(tab)[tab >= 2]
    if (length(keep_levels) < 2) next
    keep_mask  <- group_var %in% keep_levels
    alpha_vals <- alpha_vals[keep_mask]
    group_var  <- droplevels(factor(as.character(group_var[keep_mask])))
    n_levels   <- nlevels(group_var)
    if (n_levels < 2) next

    if (n_levels == 2) {
      lvls <- levels(group_var)
      grp1 <- alpha_vals[group_var == lvls[1]]
      grp2 <- alpha_vals[group_var == lvls[2]]
      n1   <- length(grp1)
      n2   <- length(grp2)

      test_res <- tryCatch(wilcox.test(grp1, grp2, exact = FALSE),
                           error = function(e) NULL)
      if (is.null(test_res)) next

      # Effect size: rank-biserial correlation, derived directly from the
      # Mann-Whitney U statistic (W) returned by wilcox.test().
      #   r_rb = 1 - 2W / (n1 * n2),  range [-1, 1]
      # Positive => level 1 tends to have higher values than level 2.
      w_stat      <- unname(test_res$statistic)
      effect_size <- 1 - (2 * w_stat) / (n1 * n2)
      effect_method <- "Rank-biserial correlation"

      res_df <- data.frame(
        category        = category,
        metric          = metric,
        n_levels        = 2,
        n_total         = length(alpha_vals),
        test            = "Wilcoxon rank-sum",
        statistic       = round(w_stat, 4),
        p_value         = test_res$p.value,
        effect_size     = round(effect_size, 4),
        effect_method   = effect_method,
        effect_direction = ifelse(effect_size > 0,
                                  paste0(lvls[1], " > ", lvls[2]),
                                  ifelse(effect_size < 0,
                                         paste0(lvls[2], " > ", lvls[1]), "no difference")),
        stringsAsFactors = FALSE
      )
    } else {
      test_res <- tryCatch(kruskal.test(alpha_vals ~ group_var),
                           error = function(e) NULL)
      if (is.null(test_res)) next

      # Effect size: epsilon-squared (H-based, rank ANOVA analogue of eta-squared)
      #   eps^2 = H / (n - 1),  range [0, 1]
      h_stat        <- unname(test_res$statistic)
      n_tot         <- length(alpha_vals)
      effect_size   <- h_stat / (n_tot - 1)
      effect_method <- "Epsilon-squared"

      res_df <- data.frame(
        category        = category,
        metric          = metric,
        n_levels        = n_levels,
        n_total         = n_tot,
        test            = "Kruskal-Wallis",
        statistic       = round(h_stat, 4),
        p_value         = test_res$p.value,
        effect_size     = round(effect_size, 4),
        effect_method   = effect_method,
        effect_direction = NA_character_,
        stringsAsFactors = FALSE
      )
    }
    results[[paste(metric, category, sep = "_")]] <- res_df
  }
  if (length(results) == 0) return(NULL)
  out <- do.call(rbind, results)
  out$effect_magnitude <- effect_label(out$effect_size, out$effect_method)
  out
}

# ==============================================================================
# AGE-ADJUSTED ALPHA DIVERSITY TESTING
# ------------------------------------------------------------------------------
# Companion to test_alpha() above. Does NOT replace or modify it - the
# unadjusted Wilcoxon/Kruskal-Wallis results from test_alpha() are kept
# exactly as-is; this adds a second, independent significance test per
# metric/category that asks whether the association survives adjustment
# for age.
#
# For each alpha metric, fits a rank-based linear model:
#     rank(alpha_diversity) ~ age_months + group_var
# and uses drop1(model, test = "F") to test the categorical term after
# accounting for age_months - a rank-based ANCOVA. rank() is used instead of
# the raw metric values so this stays in the same non-parametric spirit as
# the Wilcoxon/Kruskal-Wallis tests above, rather than assuming normality.
#
# This mirrors the age-adjusted ("partial") PERMANOVA already used for beta
# diversity (run_beta_diversity()'s condition_on argument), just applied to
# the univariate alpha diversity metrics instead of the community-level
# distance matrices.
#
# age_group and agegp are themselves age proxies by construction - adjusting
# an age variable for age is nonsensical - so NA is returned for those
# categories without attempting to fit a model.
# ==============================================================================
test_alpha_age_adjusted <- function(tse_obj, category) {
  metrics <- c("Observed", "shannon", "inverse_simpson", "gini_simpson", "chao1")
  
  # Categories for which "age-adjusted" has no sensible meaning.
  age_native_vars <- c("age_group", "agegp")
  
  age_months_vec <- colData(tse_obj)$age_months
  
  results <- list()
  
  for (metric in metrics) {
    if (!metric %in% colnames(colData(tse_obj))) next
    
    # Age variables can't sensibly be adjusted for age - record NA and move on.
    if (category %in% age_native_vars) {
      results[[paste(metric, category, sep = "_")]] <- data.frame(
        category             = category,
        metric               = metric,
        p_value_age_adjusted = NA_real_,
        model_n              = NA_integer_,
        stringsAsFactors     = FALSE
      )
      next
    }
    
    if (is.null(age_months_vec)) next  # nothing to adjust for
    
    alpha_vals <- colData(tse_obj)[[metric]]
    group_var  <- colData(tse_obj)[[category]]
    
    valid <- !is.na(alpha_vals) & !is.na(group_var) & !is.na(age_months_vec) &
      as.character(group_var) != "" & as.character(group_var) != "NA"
    if (sum(valid) < 8) next  # a bit higher than test_alpha()'s n<6 cutoff:
                              # this model has more parameters (intercept +
                              # age_months + k-1 group dummies), so needs
                              # more residual degrees of freedom to be stable
    
    alpha_vals     <- as.numeric(alpha_vals[valid])
    group_var      <- droplevels(factor(as.character(group_var[valid])))
    age_months_sub <- as.numeric(age_months_vec[valid])
    
    # Drop levels with fewer than 2 observations, exactly as test_alpha() does.
    tab <- table(group_var)
    keep_levels <- names(tab)[tab >= 2]
    if (length(keep_levels) < 2) next
    keep_mask      <- group_var %in% keep_levels
    alpha_vals     <- alpha_vals[keep_mask]
    group_var      <- droplevels(factor(as.character(group_var[keep_mask])))
    age_months_sub <- age_months_sub[keep_mask]
    
    if (nlevels(group_var) < 2) next
    if (length(alpha_vals) < nlevels(group_var) + 3) next  # keep some residual df
    
    model_df <- data.frame(
      alpha_vals = alpha_vals,
      age_months = age_months_sub,
      group_var  = group_var
    )
    
    p_adj <- tryCatch({
      fit <- lm(rank(alpha_vals) ~ age_months + group_var, data = model_df)
      dr1 <- drop1(fit, test = "F")
      row_idx <- which(rownames(dr1) == "group_var")
      if (length(row_idx) == 0) NA_real_ else dr1[["Pr(>F)"]][row_idx]
    }, error = function(e) {
      message(sprintf("    Age-adjusted model failed for '%s' (%s): %s",
                      category, metric, conditionMessage(e)))
      NA_real_
    })
    
    results[[paste(metric, category, sep = "_")]] <- data.frame(
      category             = category,
      metric               = metric,
      p_value_age_adjusted = p_adj,
      model_n              = length(alpha_vals),
      stringsAsFactors     = FALSE
    )
  }
  
  if (length(results) == 0) return(NULL)
  do.call(rbind, results)
}

# ==============================================================================
# ALPHA DIVERSITY PLOT HELPER
# ==============================================================================
plot_alpha <- function(tse_obj, category, color_vec, out_path_fn, subtitle = NULL) {
  for (metric in c("Observed","shannon","inverse_simpson","gini_simpson","chao1")) {
    vals <- colData(tse_obj)[[metric]]
    if (is.null(vals) || all(is.na(vals))) next
    alpha_df <- data.frame(
      x     = colData(tse_obj)[[category]],
      y     = vals,
      group = colData(tse_obj)[[category]]
    )
    p <- ggplot(alpha_df, aes(x=x, y=y, colour=group, fill=group)) +
      geom_violin(alpha=0.3, trim=FALSE) +
      geom_jitter(width=0.15, size=1.5, alpha=0.8) +
      scale_colour_manual(values=color_vec, name=category) +
      scale_fill_manual(values=color_vec, name=category) +
      theme_bw() +
      labs(title    = paste(metric, "by", category),
           subtitle = subtitle,
           x = category, y = metric)
    ggsave(out_path_fn(metric), plot=p, width=8, height=5, dpi=300)
  }
}

# ==============================================================================
# BETA DIVERSITY HELPER
# ==============================================================================
run_beta_diversity <- function(tse_obj, category, color_palette, file_prefix,
                               condition_on = NULL) {
  
  current_factors <- droplevels(factor(colData(tse_obj)[[category]]))
  num_levels      <- nlevels(current_factors)
  current_colors  <- color_palette[seq_len(num_levels)]
  sample_colors   <- current_colors[as.numeric(current_factors)]
  
  tse_obj <- transformAssay(tse_obj, assay.type="counts",
                            method="relabundance", name="rel")
  tse_obj <- transformAssay(tse_obj, assay.type="counts", method="clr",
                            name="clr_assay", pseudocount=1)
  
  meta_df <- as.data.frame(colData(tse_obj))
  results <- list()
  
  for (metric_label in c("bray","jaccard","aitchison")) {
    
    if (metric_label == "bray") {
      dist_matrix <- vegan::vegdist(t(assay(tse_obj,"rel")),  method="bray")
    } else if (metric_label == "jaccard") {
      dist_matrix <- vegan::vegdist(t(assay(tse_obj,"rel")),  method="jaccard", binary=TRUE)
    } else {
      dist_matrix <- dist(t(assay(tse_obj,"clr_assay")), method="euclidean")
    }
    
    message(paste("   PERMANOVA:", toupper(metric_label), "|", category))
    
    # Standard PERMANOVA
    form_std <- as.formula(paste("dist_matrix ~", category))
    perm_std <- vegan::adonis2(form_std, data=meta_df, permutations=999)
    r2_std   <- round(perm_std$R2[1], 4)
    f_std    <- round(perm_std$F[1],  4)
    p_std    <- perm_std$`Pr(>F)`[1]
    
    # Partial PERMANOVA - condition out age proxies using by = "margin"
    r2_partial           <- NA_real_
    f_partial            <- NA_real_
    p_partial            <- NA_real_
    conditioned_on_label <- NA_character_
    
    if (!is.null(condition_on)) {
      valid_conditions <- condition_on[sapply(condition_on, function(cv) {
        !is.null(meta_df[[cv]]) &&
          sum(!is.na(meta_df[[cv]])) > 0 &&
          length(unique(na.omit(as.character(meta_df[[cv]])))) > 1 &&
          cv != category
      })]
      
      if (length(valid_conditions) > 0) {
        
        # Subset to rows complete across ALL formula terms
        partial_cols  <- c(valid_conditions, category)
        complete_mask <- complete.cases(meta_df[, partial_cols, drop=FALSE])
        
        if (sum(complete_mask) >= 6 &&
            length(unique(na.omit(
              as.character(meta_df[[category]][complete_mask])))) >= 2) {
          
          dist_matrix_partial <- as.dist(
            as.matrix(dist_matrix)[complete_mask, complete_mask])
          meta_partial <- meta_df[complete_mask, , drop=FALSE]
          
          for (cv in valid_conditions) {
            meta_partial[[cv]] <- droplevels(
              factor(as.character(meta_partial[[cv]])))
          }
          
          cond_terms   <- paste(valid_conditions, collapse=" + ")
          form_partial <- as.formula(
            paste("dist_matrix_partial ~", cond_terms, "+", category))
          
          tryCatch({
            perm_partial <- vegan::adonis2(form_partial, data=meta_partial,
                                           by="margin", permutations=999)
            cat_row <- which(rownames(perm_partial) == category)
            if (length(cat_row) > 0) {
              r2_partial           <- round(perm_partial$R2[cat_row], 4)
              f_partial            <- round(perm_partial$F[cat_row],  4)
              p_partial            <- perm_partial$`Pr(>F)`[cat_row]
              conditioned_on_label <- paste(valid_conditions, collapse=" + ")
            }
          }, error = function(e) {
            message(paste("    Partial PERMANOVA failed:", e$message))
          })
          
        } else {
          message("    Partial PERMANOVA skipped: insufficient complete cases.")
        }
      }
    }
    
    results[[metric_label]] <- data.frame(
      category             = category,
      metric               = metric_label,
      R2                   = r2_std,
      F_stat               = f_std,
      p_value              = p_std,
      R2_partial           = r2_partial,
      F_stat_partial       = f_partial,
      p_value_partial      = p_partial,
      conditioned_on       = conditioned_on_label,
      context              = NA_character_,
      label                = NA_character_,
      stringsAsFactors     = FALSE
    )
    
    # 3D PCoA
    pcoa_res    <- cmdscale(dist_matrix, k=3, eig=TRUE)
    pcoa_coords <- pcoa_res$points
    
    png(filename=paste0(file_prefix,"_pcoa_",metric_label,".png"),
        width=8, height=7, units="in", res=300)
    s3d <- scatterplot3d(
      x=as.numeric(pcoa_coords[,1]),
      y=as.numeric(pcoa_coords[,2]),
      z=as.numeric(pcoa_coords[,3]),
      color=sample_colors, pch=19, cex.symbols=1.2, grid=TRUE, box=FALSE,
      xlab="PC1", ylab="PC2", zlab="PC3",
      main=paste0("3D PCoA (",toupper(metric_label),")\n",
                  gsub("_"," ",basename(file_prefix))," (p=",p_std,")")
    )
    for (g in levels(current_factors)) {
      idx <- which(current_factors == g)
      if (length(idx) > 1) {
        cen <- colMeans(pcoa_coords[idx,1:3,drop=FALSE])
        for (i in idx)
          s3d$points3d(
            x=c(pcoa_coords[i,1],cen[1]),
            y=c(pcoa_coords[i,2],cen[2]),
            z=c(pcoa_coords[i,3],cen[3]),
            type="l",
            col=adjustcolor(
              current_colors[which(levels(current_factors)==g)], alpha.f=0.3))
      }
    }
    legend(s3d$xyz.convert(max(pcoa_coords[,1])*0.7,
                           max(pcoa_coords[,2]),
                           min(pcoa_coords[,3])),
           legend=levels(current_factors), col=current_colors,
           pch=19, bty="n", cex=0.8)
    dev.off()
    
    # 2D NMDS
    nmds_res    <- vegan::metaMDS(dist_matrix, k=2, trace=FALSE, trymax=100)
    nmds_coords <- nmds_res$points
    
    nmds_df  <- data.frame(NMDS1=nmds_coords[,1], NMDS2=nmds_coords[,2],
                           group=colData(tse_obj)[[category]])
    ell_grps <- names(table(nmds_df$group)[table(nmds_df$group) >= 4])
    nmds_ell <- nmds_df[nmds_df$group %in% ell_grps, ]
    
    nmds_plot <- ggplot(nmds_df, aes(x=NMDS1, y=NMDS2, colour=group)) +
      geom_point(size=2) +
      { if (nrow(nmds_ell) > 0)
        suppressWarnings(
          stat_ellipse(data=nmds_ell, type="t", level=0.95, linewidth=0.8))
        else NULL } +
      scale_colour_manual(values=current_colors, name=category) +
      theme_bw() +
      labs(title    = paste0("2D NMDS (",toupper(metric_label)," Dissimilarity)"),
           subtitle = paste0(gsub("_"," ",basename(file_prefix)),
                             " | Stress:",round(nmds_res$stress,4),
                             " | p=",p_std,
                             if (!is.na(p_partial))
                               paste0(" | p(age-adj)=",p_partial) else ""),
           x="NMDS 1", y="NMDS 2")
    ggsave(paste0(file_prefix,"_nmds_",metric_label,".png"),
           plot=nmds_plot, width=8, height=6, dpi=300)
  }
  
  do.call(rbind, results)
}

# ==============================================================================
# STRATIFIED ANALYSIS HELPER
# ------------------------------------------------------------------------------
# UPDATED: now runs test_alpha() alongside run_beta_diversity() and returns
# BOTH as a list(beta = ..., alpha = ...) instead of returning beta rows only.
# Every call site below has been updated to use $beta / $alpha accordingly.
# ==============================================================================
run_stratified <- function(tse_subset, strat_var, strat_colors,
                           color_palette, out_dir, label, subtitle,
                           condition_on = NULL) {
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)
  
  colData(tse_subset)[[strat_var]] <-
    droplevels(factor(as.character(colData(tse_subset)[[strat_var]])))
  
  if (nlevels(colData(tse_subset)[[strat_var]]) < 2 || ncol(tse_subset) < 6)
    return(list(beta = NULL, alpha = NULL))
  
  tse_subset <- ensure_alpha(tse_subset)
  
  lvls   <- levels(colData(tse_subset)[[strat_var]])
  colors <- if (!is.null(strat_colors)) strat_colors else color_palette[seq_along(lvls)]
  
  plot_alpha(
    tse_obj     = tse_subset,
    category    = strat_var,
    color_vec   = colors,
    out_path_fn = function(m) file.path(out_dir, paste0("alpha_",m,"_",label,".png")),
    subtitle    = subtitle
  )
  
  beta_rows <- tryCatch(
    run_beta_diversity(
      tse_obj       = tse_subset,
      category      = strat_var,
      color_palette = colors,
      file_prefix   = file.path(out_dir, paste0("beta_",label)),
      condition_on  = condition_on
    ),
    error = function(e) { message(paste("    Beta diversity failed:", e$message)); NULL }
  )
  
  alpha_rows <- tryCatch(test_alpha(tse_subset, strat_var),
                        error = function(e) { message(paste("    Alpha test failed:", e$message)); NULL })
  
  if (!is.null(beta_rows)) {
    beta_rows$label   <- label
    beta_rows$context <- subtitle
  }
  if (!is.null(alpha_rows)) {
    alpha_rows$label   <- label
    alpha_rows$context <- subtitle
  }
  
  list(beta = beta_rows, alpha = alpha_rows)
}

# ==============================================================================
# PERMANOVA VISUALISATION FUNCTION
# ==============================================================================
plot_permanova_results <- function(results_df, out_dir,
                                   q_thresholds   = c(0.001, 0.01, 0.05),
                                   context_filter = "main") {
  
  plot_dir <- file.path(out_dir, "permanova_plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive=TRUE)
  
  df <- results_df[!is.na(results_df$context) &
                     results_df$context == context_filter, ]
  if (nrow(df) == 0) {
    message("No rows match context_filter = '", context_filter,
            "'. Skipping permanova plots.")
    return(invisible(NULL))
  }
  
  df$fdr_tier <- cut(
    df$q_value,
    breaks = c(-Inf, q_thresholds[1], q_thresholds[2], q_thresholds[3], Inf),
    labels = c(
      paste0("q \u2264 ", q_thresholds[1]),
      paste0("q \u2264 ", q_thresholds[2]),
      paste0("q \u2264 ", q_thresholds[3]),
      "Not significant"
    ),
    right = TRUE
  )
  
  fdr_colours <- setNames(
    c("#B2182B", "#E6614A", "#F4A460", "#CCCCCC"),
    c(paste0("q \u2264 ", q_thresholds[1]),
      paste0("q \u2264 ", q_thresholds[2]),
      paste0("q \u2264 ", q_thresholds[3]),
      "Not significant")
  )
  
  df$age_proxy <- !is.na(df$p_value_partial) &
    df$q_value         <  0.05 &
    df$q_value_partial >= 0.05
  
  # ----------------------------------------------------------
  # PLOT 1: Bubble plot
  # ----------------------------------------------------------
  message("  Plotting: main effect bubble plot")
  
  var_order     <- tapply(df$R2, df$category, mean, na.rm=TRUE)
  var_order     <- names(sort(var_order, decreasing=FALSE))
  df$category_f <- factor(df$category, levels=var_order)
  
  f_cap      <- quantile(df$F_stat, 0.95, na.rm=TRUE)
  df$F_capped <- pmin(df$F_stat, f_cap)
  
  p1 <- ggplot(df, aes(x=R2, y=category_f,
                       size=F_capped, colour=fdr_tier, shape=age_proxy)) +
    geom_vline(xintercept=0, colour="grey80", linewidth=0.3) +
    geom_point(alpha=0.85) +
    scale_colour_manual(values=fdr_colours, name="FDR tier", drop=FALSE) +
    scale_size_continuous(name  = paste0("Pseudo-F\n(capped at ",round(f_cap,1),")"),
                          range = c(1.5, 8)) +
    scale_shape_manual(values = c("FALSE"=16, "TRUE"=8),
                       labels = c("FALSE"="Independent effect",
                                  "TRUE" ="Likely age proxy \u26a0"),
                       name   = "") +
    facet_wrap(~metric, nrow=1) +
    theme_bw(base_size=11) +
    theme(
      panel.grid.major.y = element_line(colour="grey93"),
      panel.grid.major.x = element_line(colour="grey93"),
      axis.text.y        = element_text(size=8),
      strip.background   = element_rect(fill="#2C3E50"),
      strip.text         = element_text(colour="white", face="bold"),
      legend.position    = "right"
    ) +
    labs(
      title    = "PERMANOVA Effect Sizes with FDR Correction",
      subtitle = paste0("Bubble size = Pseudo-F statistic | ",
                        "\u26a0 = significant unadjusted but not after age conditioning"),
      x        = expression(R^2~"(proportion of community variance explained)"),
      y        = NULL
    )
  
  plot_height <- max(6, length(var_order) * 0.28 + 2.5)
  ggsave(file.path(plot_dir,"permanova_bubble_plot.png"),
         plot=p1, width=14, height=plot_height, dpi=300, limitsize=FALSE)
  
  # ----------------------------------------------------------
  # PLOT 2: Standard vs age-adjusted R2 scatter
  # ----------------------------------------------------------
  message("  Plotting: standard vs age-adjusted R2 comparison")
  
  df_partial <- df[!is.na(df$R2_partial), ]
  
  if (nrow(df_partial) > 0) {
    
    df_partial$sig_change <- factor(
      dplyr::case_when(
        df_partial$q_value <  0.05 & df_partial$q_value_partial <  0.05 ~
          "Significant both",
        df_partial$q_value <  0.05 & df_partial$q_value_partial >= 0.05 ~
          "Lost after age adjustment",
        df_partial$q_value >= 0.05 & df_partial$q_value_partial <  0.05 ~
          "Gained after age adjustment",
        TRUE ~ "Not significant"
      ),
      levels = c("Significant both","Lost after age adjustment",
                 "Gained after age adjustment","Not significant")
    )
    
    change_colours <- setNames(
      c("#2C7BB6", "#D7191C", "#1A9641", "#CCCCCC"),
      c("Significant both","Lost after age adjustment",
        "Gained after age adjustment","Not significant")
    )
    
    r2_max <- max(c(df_partial$R2, df_partial$R2_partial), na.rm=TRUE) * 1.05
    
    p2 <- ggplot(df_partial,
                 aes(x=R2, y=R2_partial, colour=sig_change, label=category)) +
      geom_abline(slope=1, intercept=0,
                  linetype="dashed", colour="grey50", linewidth=0.5) +
      geom_point(size=3, alpha=0.85) +
      geom_text_repel(
        data          = df_partial[df_partial$sig_change != "Not significant", ],
        size          = 2.8, max.overlaps=20,
        segment.colour="grey70", segment.size=0.3
      ) +
      scale_colour_manual(values=change_colours, name="FDR outcome") +
      scale_x_continuous(limits=c(0, r2_max)) +
      scale_y_continuous(limits=c(0, r2_max)) +
      coord_equal() +
      facet_wrap(~metric, nrow=1) +
      theme_bw(base_size=11) +
      theme(
        strip.background = element_rect(fill="#2C3E50"),
        strip.text       = element_text(colour="white", face="bold"),
        legend.position  = "bottom"
      ) +
      labs(
        title    = "Effect Size Before vs After Age Conditioning",
        subtitle = paste0("Points below the diagonal lost variance explained after ",
                          "accounting for age - likely age proxies\n",
                          "One point per variable per beta metric"),
        x        = expression(R^2~"(standard PERMANOVA)"),
        y        = expression(R^2~"(age-adjusted partial PERMANOVA)")
      )
    
    ggsave(file.path(plot_dir,"permanova_age_adjustment_scatter.png"),
           plot=p2, width=14, height=6, dpi=300)
    
  } else {
    message("  Skipping age-adjustment scatter: no partial PERMANOVA results.")
  }
  
  # ----------------------------------------------------------
  # PLOT 3: FDR significance heatmap
  # ----------------------------------------------------------
  message("  Plotting: FDR significance heatmap")
  
  df$log10q    <- pmin(-log10(df$q_value + 1e-10), 4)
  df$stars     <- sig_label(df$q_value)
  df$stars_adj <- ifelse(!is.na(df$q_value_partial),
                         sig_label(df$q_value_partial), "")
  
  df$annotation <- paste0(
    df$stars,
    ifelse(df$stars_adj != "" & df$stars_adj != df$stars,
           paste0("\n(", df$stars_adj, ")"), "")
  )
  
  heat_order    <- tapply(df$log10q, df$category, mean, na.rm=TRUE)
  heat_order    <- names(sort(heat_order, decreasing=TRUE))
  df$category_h <- factor(df$category, levels=rev(heat_order))
  df$metric_f   <- factor(toupper(df$metric),
                          levels=c("BRAY","JACCARD","AITCHISON"))
  
  p3 <- ggplot(df, aes(x=metric_f, y=category_h, fill=log10q)) +
    geom_tile(colour="white", linewidth=0.5) +
    geom_text(aes(label=annotation),
              size=2.8, lineheight=0.9, colour="white", fontface="bold") +
    scale_fill_gradientn(
      colours = c("#F7F7F7","#FDD49E","#FC8D59","#D7301F","#7F0000"),
      values  = scales::rescale(c(0,1,2,3,4)),
      limits  = c(0,4),
      name    = expression(-log[10](q)),
      guide   = guide_colourbar(barwidth=8, barheight=0.8,
                                title.position="top", title.hjust=0.5)
    ) +
    scale_x_discrete(position="top") +
    theme_bw(base_size=11) +
    theme(
      axis.text.x     = element_text(face="bold", size=10),
      axis.text.y     = element_text(size=8),
      panel.grid      = element_blank(),
      legend.position = "bottom",
      plot.title      = element_text(face="bold"),
      plot.subtitle   = element_text(size=9, colour="grey40")
    ) +
    labs(
      title    = "PERMANOVA FDR Significance Heatmap",
      subtitle = paste0("Fill = -log10(FDR q-value), capped at 4\n",
                        "Stars = FDR tier | Bracketed stars = age-adjusted significance"),
      x        = NULL, y = NULL
    )
  
  heat_height <- max(5, length(heat_order) * 0.22 + 2.5)
  ggsave(file.path(plot_dir,"permanova_fdr_heatmap.png"),
         plot=p3, width=7, height=heat_height, dpi=300, limitsize=FALSE)
  
  message(paste("PERMANOVA plots saved to:", plot_dir))
  invisible(list(bubble=p1, heatmap=p3))
}

# ==============================================================================
# ALPHA DIVERSITY EFFECT SIZE / P-VALUE / FDR VISUALISATION
# ------------------------------------------------------------------------------
# Mirrors plot_permanova_results() above, but for the non-parametric alpha
# diversity tests. Wilcoxon (rank-biserial correlation, range [-1,1]) and
# Kruskal-Wallis (epsilon-squared, range [0,1]) results are plotted separately
# since their effect-size scales aren't directly comparable.
# ==============================================================================
plot_alpha_effect_results <- function(alpha_df, out_dir,
                                      q_thresholds   = c(0.001, 0.01, 0.05),
                                      context_filter = "main") {
  
  plot_dir <- file.path(out_dir, "alpha_diversity_plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive=TRUE)
  
  df <- alpha_df[!is.na(alpha_df$context) & alpha_df$context == context_filter, ]
  if (nrow(df) == 0) {
    message("No rows match context_filter = '", context_filter,
            "'. Skipping alpha diversity effect plots.")
    return(invisible(NULL))
  }
  
  df$fdr_tier <- cut(
    df$q_value,
    breaks = c(-Inf, q_thresholds[1], q_thresholds[2], q_thresholds[3], Inf),
    labels = c(
      paste0("q \u2264 ", q_thresholds[1]),
      paste0("q \u2264 ", q_thresholds[2]),
      paste0("q \u2264 ", q_thresholds[3]),
      "Not significant"
    ),
    right = TRUE
  )
  
  fdr_colours <- setNames(
    c("#B2182B", "#E6614A", "#F4A460", "#CCCCCC"),
    c(paste0("q \u2264 ", q_thresholds[1]),
      paste0("q \u2264 ", q_thresholds[2]),
      paste0("q \u2264 ", q_thresholds[3]),
      "Not significant")
  )
  
  df$neglog10p <- pmin(-log10(df$p_value + 1e-300), 10)
  
  # NULL-safe: age_proxy may be absent if this is called on alpha results
  # that were never age-adjusted (e.g. an older results object).
  if (!"age_proxy" %in% colnames(df)) df$age_proxy <- FALSE
  df$age_proxy[is.na(df$age_proxy)] <- FALSE
  
  make_bubble <- function(sub_df, effect_label_text, x_limits, file_name, title_suffix) {
    if (nrow(sub_df) == 0) return(invisible(NULL))
    
    var_order       <- tapply(abs(sub_df$effect_size), sub_df$category, mean, na.rm=TRUE)
    var_order       <- names(sort(var_order, decreasing=FALSE))
    sub_df$category_f <- factor(sub_df$category, levels=var_order)
    
    p <- ggplot(sub_df, aes(x=effect_size, y=category_f,
                            size=neglog10p, colour=fdr_tier, shape=age_proxy)) +
      geom_vline(xintercept=0, colour="grey80", linewidth=0.3) +
      geom_point(alpha=0.85) +
      scale_colour_manual(values=fdr_colours, name="FDR tier", drop=FALSE) +
      scale_size_continuous(name = expression(-log[10](p)~"(capped at 10)"),
                            range = c(1.5, 8)) +
      scale_shape_manual(values = c("FALSE"=16, "TRUE"=8),
                         labels = c("FALSE"="Independent effect",
                                    "TRUE" ="Likely age proxy \u26a0"),
                         name   = "") +
      scale_x_continuous(limits = x_limits) +
      facet_wrap(~metric, nrow=1) +
      theme_bw(base_size=11) +
      theme(
        panel.grid.major.y = element_line(colour="grey93"),
        panel.grid.major.x = element_line(colour="grey93"),
        axis.text.y        = element_text(size=8),
        strip.background   = element_rect(fill="#2C3E50"),
        strip.text         = element_text(colour="white", face="bold"),
        legend.position    = "right"
      ) +
      labs(
        title    = paste0("Alpha Diversity Effect Sizes with FDR Correction (", title_suffix, ")"),
        subtitle = paste0("Bubble size = -log10(p) | Colour = FDR tier | ",
                          "\u26a0 = significant unadjusted but not after age conditioning"),
        x        = effect_label_text,
        y        = NULL
      )
    
    plot_height <- max(6, length(var_order) * 0.28 + 2.5)
    ggsave(file.path(plot_dir, file_name),
           plot=p, width=12, height=plot_height, dpi=300, limitsize=FALSE)
    invisible(p)
  }
  
  # Wilcoxon (2-level variables): rank-biserial correlation, signed, [-1, 1]
  wilcox_df <- df[df$test == "Wilcoxon rank-sum", ]
  make_bubble(wilcox_df,
             expression("Rank-biserial correlation "*italic(r)),
             c(-1, 1),
             "alpha_effect_bubble_wilcoxon.png",
             "Wilcoxon, 2-level variables")
  
  # Kruskal-Wallis (3+ level variables): epsilon-squared, unsigned, [0, 1]
  kw_df <- df[df$test == "Kruskal-Wallis", ]
  make_bubble(kw_df,
             expression("Epsilon-squared "*(epsilon^2)),
             c(0, max(1, max(kw_df$effect_size, na.rm=TRUE) * 1.05)),
             "alpha_effect_bubble_kruskal.png",
             "Kruskal-Wallis, 3+ level variables")
  
  # Combined FDR significance heatmap across both test types
  message("  Plotting: alpha diversity FDR significance heatmap")
  
  df$log10q    <- pmin(-log10(df$q_value + 1e-10), 4)
  df$stars     <- sig_label(df$q_value)
  df$stars_adj <- ifelse(!is.na(df$q_value_age_adjusted),
                         sig_label(df$q_value_age_adjusted), "")
  
  df$annotation <- paste0(
    df$stars,
    ifelse(df$stars_adj != "" & df$stars_adj != df$stars,
           paste0("\n(", df$stars_adj, ")"), "")
  )
  
  heat_order      <- tapply(df$log10q, df$category, mean, na.rm=TRUE)
  heat_order      <- names(sort(heat_order, decreasing=TRUE))
  df$category_h   <- factor(df$category, levels=rev(heat_order))
  df$metric_f     <- factor(df$metric,
                            levels=c("Observed","shannon","inverse_simpson",
                                    "gini_simpson","chao1"))
  
  p_heat <- ggplot(df, aes(x=metric_f, y=category_h, fill=log10q)) +
    geom_tile(colour="white", linewidth=0.5) +
    geom_text(aes(label=annotation),
              size=2.8, lineheight=0.9, colour="white", fontface="bold") +
    scale_fill_gradientn(
      colours = c("#F7F7F7","#FDD49E","#FC8D59","#D7301F","#7F0000"),
      values  = scales::rescale(c(0,1,2,3,4)),
      limits  = c(0,4),
      name    = expression(-log[10](q)),
      guide   = guide_colourbar(barwidth=8, barheight=0.8,
                                title.position="top", title.hjust=0.5)
    ) +
    scale_x_discrete(position="top") +
    theme_bw(base_size=11) +
    theme(
      axis.text.x     = element_text(face="bold", size=9, angle=20, hjust=0),
      axis.text.y     = element_text(size=8),
      panel.grid      = element_blank(),
      legend.position = "bottom",
      plot.title      = element_text(face="bold"),
      plot.subtitle   = element_text(size=9, colour="grey40")
    ) +
    labs(
      title    = "Alpha Diversity FDR Significance Heatmap",
      subtitle = paste0("Fill = -log10(FDR q-value), capped at 4\n",
                        "Stars = FDR tier | Bracketed stars = age-adjusted significance"),
      x        = NULL, y = NULL
    )
  
  heat_height <- max(5, length(heat_order) * 0.22 + 2.5)
  ggsave(file.path(plot_dir,"alpha_diversity_fdr_heatmap.png"),
         plot=p_heat, width=8, height=heat_height, dpi=300, limitsize=FALSE)
  
  message(paste("Alpha diversity effect/FDR plots saved to:", plot_dir))
  invisible(list(heatmap=p_heat))
}

# ==============================================================================
# 1. PARSE & INTEGERIZE MATRIX
# ==============================================================================
message("\n--- Step 1: Parsing MetaPhlAn Counts Matrix ---")
raw_data      <- read.delim(merged_file, skip=1, check.names=FALSE)
filtered_data <- raw_data[raw_data$clade_name != "UNCLASSIFIED", ]
temp_file     <- tempfile(fileext=".tsv")
write.table(filtered_data, file=temp_file, sep="\t", row.names=FALSE, quote=FALSE)

tse_base     <- importMetaPhlAn(file=temp_file)
assay_matrix <- assay(tse_base, "metaphlan")
assay(tse_base, "counts") <- round(as.matrix(assay_matrix))

# ==============================================================================
# 2. LOAD METADATA
# ==============================================================================
message("\n--- Step 2: Loading metadata and applying WHO Anthro flagged-sample exclusion ---")
metadata <- read.csv(metadata_file, stringsAsFactors=FALSE, check.names=FALSE,
                     na.strings=c("", "NA"))

# Exclude samples only when their lab_id matches a WHO Anthro row with a
# non-empty Flag. Blank/missing flags and samples without a WHO match are kept.
who_file <- "WHO_Anthro_Nutritional_Survey_Malawi.csv"
if (!file.exists(who_file)) {
  stop("WHO Anthro flag file not found: ", who_file,
       ". Place it in the working directory or update who_file.")
}
who <- read.csv(who_file, stringsAsFactors=FALSE, check.names=FALSE,
                na.strings=c("", "NA"))

required_metadata_cols <- c("sample_id", "lab_id")
missing_metadata_cols <- setdiff(required_metadata_cols, colnames(metadata))
if (length(missing_metadata_cols)) {
  stop("Metadata file missing required column(s): ",
       paste(missing_metadata_cols, collapse=", "))
}
required_who_cols <- c("lab_id", "Flag")
missing_who_cols <- setdiff(required_who_cols, colnames(who))
if (length(missing_who_cols)) {
  stop("WHO Anthro file missing required column(s): ",
       paste(missing_who_cols, collapse=", "))
}

who_lab_ids <- trimws(as.character(who$lab_id))
who_flags   <- trimws(as.character(who$Flag))
flagged_lab_ids <- unique(who_lab_ids[
  !is.na(who_lab_ids) & nzchar(who_lab_ids) &
  !is.na(who_flags) & nzchar(who_flags)
])
message(sprintf("WHO Anthro: %d unique lab_id(s) have a non-empty Flag.",
                length(flagged_lab_ids)))

metadata$lab_id   <- trimws(as.character(metadata$lab_id))
metadata$sample_id <- trimws(as.character(metadata$sample_id))
flagged_metadata <- !is.na(metadata$lab_id) &
                    metadata$lab_id %in% flagged_lab_ids
flagged_sample_ids <- unique(metadata$sample_id[flagged_metadata &
                                                 !is.na(metadata$sample_id) &
                                                 nzchar(metadata$sample_id)])

sample_names_all <- colnames(tse_base)
base_ids_all <- gsub("_metaphlan_wreads$", "", sample_names_all)
exclude_data_cols <- base_ids_all %in% flagged_sample_ids
if (any(exclude_data_cols)) {
  message(sprintf("WHO Anthro exclusion: removing %d MetaPhlAn sample(s) matched by lab_id.",
                  sum(exclude_data_cols)))
  tse_base <- tse_base[, !exclude_data_cols]
} else {
  message("WHO Anthro exclusion: no flagged samples matched the MetaPhlAn sample IDs.")
}

# Save an audit trail of flagged metadata rows and sample IDs.
if (any(flagged_metadata)) {
  write.csv(metadata[flagged_metadata, c("sample_id", "lab_id"), drop=FALSE],
            file=file.path(output_dir, "who_anthro_excluded_samples.csv"),
            row.names=FALSE)
}

sample_names_in_data <- colnames(tse_base)
mapping_df <- data.frame(
  data_colname = sample_names_in_data,
  base_id      = gsub("_metaphlan_wreads$", "", sample_names_in_data),
  stringsAsFactors = FALSE
)
merged_meta <- merge(mapping_df, metadata,
                     by.x="base_id", by.y="sample_id", all.x=TRUE, sort=FALSE)
rownames(merged_meta) <- merged_meta$data_colname
merged_meta <- merged_meta[sample_names_in_data, , drop=FALSE]

# ==============================================================================
# 3. MAP METADATA
# ==============================================================================
message("\n--- Step 3: Mapping metadata variables ---")
cd <- colData(tse_base)

## Continuous
cd$age_days       <- clean_numeric(merged_meta$age_days,       0,   1825)
cd$age_months     <- clean_numeric(merged_meta$age_months,     0,   60)
cd$muac           <- clean_numeric(merged_meta$muac,           50,  250)
cd$height         <- clean_numeric(merged_meta$height,         30,  150)
cd$weight         <- clean_numeric(merged_meta$weight,         1,   30)
cd$temperature    <- clean_numeric(merged_meta$temperature,    30,  45)
cd$hemocue_result <- clean_numeric(merged_meta$hemocue_result, 2,   25)
cd$sample_volume  <- clean_numeric(merged_meta$sample_volume,  0,   Inf)
cd$mother_age     <- clean_numeric(merged_meta$mother_age,     12,  60)

for (v in c("malaria_prev","averagempoconcentration","mrdt_percent",
            "mrdt_percent_fup","mrdt_percent_comb")) {
  col <- find_col(v, merged_meta)
  if (!is.null(col)) cd[[v]] <- clean_numeric(merged_meta[[col]], 0)
}

for (vs in list(c("heart_rate",         40,  250),
                c("resp_rate",          5,   100),
                c("oxy_satu",           50,  100),
                c("bld_pressure",       40,  200),
                c("systolic_pressure",  40,  200),
                c("d90_v2_temperature", 30,  45))) {
  col <- find_col(vs[1], merged_meta)
  if (!is.null(col)) cd[[vs[1]]] <- clean_numeric(merged_meta[[col]],
                                                  as.numeric(vs[2]),
                                                  as.numeric(vs[3]))
}

for (gv in c("gps_latitude","gps_longitude")) {
  col <- find_col(gv, merged_meta)
  if (!is.null(col)) cd[[gv]] <- as.numeric(merged_meta[[col]])
}

## Core factors
cd$age_group <- as.character(merged_meta$age_group)

agegp_col <- find_col("^agegp$", merged_meta)
if (!is.null(agegp_col)) cd$agegp <- as.character(merged_meta[[agegp_col]])

gender_src <- find_col(c("^gender$","^sex$"), merged_meta)
raw_gender <- as.numeric(merged_meta[[gender_src]])
cd$gender  <- factor(ifelse(raw_gender==1,"Male",ifelse(raw_gender==0,"Female",NA)),
                     levels=c("Female","Male"))

## Socioeconomic
econ_col <- find_col(c("economic_step","econ_step","ladder"), merged_meta)
if (!is.null(econ_col))
  cd$economic_step <- recode_factor(as.numeric(merged_meta[[econ_col]]),
                                    c("1"="Step 1 (Poorest)","2"="Step 2","3"="Step 3",
                                      "4"="Step 4","5"="Step 5","6"="Step 6 (Richest)"))

house_col <- find_col(c("housing","house_adequacy","house_adeq"), merged_meta)
if (!is.null(house_col))
  cd$housing <- recode_factor(as.numeric(merged_meta[[house_col]]),
                              c("1"="Less than adequate","2"="Just adequate","3"="More than adequate"))

clothing_col <- find_col("^clothing$", merged_meta)
if (!is.null(clothing_col))
  cd$clothing <- recode_factor(as.numeric(merged_meta[[clothing_col]]),
                               c("1"="Less than adequate","2"="Just adequate","3"="More than adequate"))

edu_col <- find_col(c("^education$","educ"), merged_meta)
if (!is.null(edu_col))
  cd$education <- recode_factor(as.numeric(merged_meta[[edu_col]]),
                                c("1"="Never been to school","2"="Primary",
                                  "3"="Secondary","4"="Tertiary/College"))

hh_col <- find_col(c("num_of_people","hh_size","hhsize","household_size","num_people"),
                   merged_meta)
if (!is.null(hh_col)) {
  raw_hh <- as.numeric(merged_meta[[hh_col]])
  sz     <- rep(NA_character_, length(raw_hh))
  sz[raw_hh>=1 & raw_hh<=2] <- "1-2"
  sz[raw_hh>=3 & raw_hh<=6] <- "3-6"
  sz[raw_hh>=7]              <- "7+"
  cd$household_size <- factor(sz, levels=c("1-2","3-6","7+"))
}

food_col <- find_col(c("enough_food","food_security","enoughfood"), merged_meta)
if (!is.null(food_col))
  cd$enough_food <- recode_yesno(as.numeric(merged_meta[[food_col]]))

rw_col <- find_col("read_write", merged_meta)
if (!is.null(rw_col)) {
  raw_rw <- as.numeric(merged_meta[[rw_col]]); raw_rw[raw_rw==2] <- NA
  cd$read_write <- recode_yesno(raw_rw)
}

lighting_col <- find_col("^lighting$", merged_meta)
if (!is.null(lighting_col))
  cd$lighting <- recode_factor(as.numeric(merged_meta[[lighting_col]]),
                               c("1"="Electric","2"="Solar","3"="Paraffin/Oil","4"="Battery",
                                 "5"="Candles","6"="No light","7"="Other","8"="Unknown"), na_codes=c())

toilet_col <- find_col("toilet_type", merged_meta)
if (!is.null(toilet_col))
  cd$toilet_type <- recode_factor(as.numeric(merged_meta[[toilet_col]]),
                                  c("1"="Flush/pour flush","2"="Pit latrine (slab)",
                                    "3"="Pit latrine (earth)","4"="Open defaecation"), na_codes=c(5,99))

for (wv in c("cooking_water_src","washing_water_src",
             "bathing_water_src","drinking_water_src")) {
  col <- find_col(wv, merged_meta)
  if (!is.null(col)) cd[[wv]] <- recode_water(as.numeric(merged_meta[[col]]))
}

treat_col <- find_col("treat_dwater", merged_meta)
if (!is.null(treat_col))
  cd$treat_dwater <- recode_yesno(as.numeric(merged_meta[[treat_col]]))

for (av in c("own_bed","own_table","financial_account","electricity",
             "chair_sofa","keep_animals","companion_pet","livestock")) {
  col <- find_col(paste0("^",av,"$"), merged_meta)
  if (!is.null(col)) cd[[av]] <- recode_yesno(as.numeric(merged_meta[[col]]))
}

sleep_col <- find_col("sleeping_material", merged_meta)
if (!is.null(sleep_col))
  cd$sleeping_material <- recode_factor(as.numeric(merged_meta[[sleep_col]]),
                                        c("1"="Bed frame + Mattress","2"="Bed frame + Mat","3"="Bed frame alone",
                                          "4"="Mattress on floor","5"="Mat on floor","6"="Cloth/Sack on floor",
                                          "7"="Floor (nothing)","8"="Other"), na_codes=c(99))

loc_col <- find_col("^location$", merged_meta)
if (!is.null(loc_col))
  cd$location <- recode_factor(as.numeric(merged_meta[[loc_col]]),
                               c("1"="Participant's Home","2"="Study Tent","3"="Other"), na_codes=c())

clus_col <- find_col("^cluster$", merged_meta)
if (!is.null(clus_col))
  cd$cluster <- recode_factor(as.numeric(merged_meta[[clus_col]]),
                              c("1"="Chikwawa","2"="Montfort"), na_codes=c())

for (iv in c("intensive_area","intensive_route")) {
  col <- find_col(iv, merged_meta)
  if (!is.null(col)) cd[[iv]] <- recode_yesno(as.numeric(merged_meta[[col]]))
}

## Clinical
# FIX: previously used merged_meta$diarrhoea / merged_meta$excl_brstfed directly.
# Every other variable in this script uses find_col() so a slightly different
# column name (case, underscore, etc.) degrades gracefully instead of throwing
# an error or silently producing a zero-length vector. Made consistent here.
diarrhoea_col <- find_col("^diarrhoea$", merged_meta)
if (!is.null(diarrhoea_col))
  cd$diarrhoea <- recode_yesno(as.numeric(merged_meta[[diarrhoea_col]]), c(88,99))

brstfed_col <- find_col("excl_brstfed", merged_meta)
if (!is.null(brstfed_col))
  cd$excl_brstfed <- recode_yesno(as.numeric(merged_meta[[brstfed_col]]), c(9,88,99))

oedema_col <- find_col("^oedema$", merged_meta)
if (!is.null(oedema_col))
  cd$oedema <- recode_yesno(as.numeric(merged_meta[[oedema_col]]), na_codes=c())

for (yn in c("fever","tested_malaria","receive_malariatreatment","hiv_exposed",
             "art_treatment","anaemia","malnutrition","lost_weight","eating_less",
             "weight_gain","underweight","chronic_illness_yn","sickle_cell","mrdt",
             "d90_v2_fever","d90_v2_anaemia","d90_v2_malnutrition")) {
  col <- find_col(paste0("^",yn,"$"), merged_meta)
  if (!is.null(col)) cd[[yn]] <- recode_yesno(as.numeric(merged_meta[[col]]))
}

for (pn in c("malariatest_result","hivtest_result","mrdt_result",
             "d90_fup_mrdt_result","d90_stl_mrdt_result","d90_v2_malariatest_result")) {
  col <- find_col(paste0("^",pn,"$"), merged_meta)
  if (!is.null(col)) {
    raw <- as.numeric(merged_meta[[col]])
    cd[[pn]] <- factor(ifelse(raw==1,"Positive",ifelse(raw==0,"Negative",NA)),
                       levels=c("Negative","Positive"))
  }
}

mrdt_master_col <- find_col("mrdt_result_master", merged_meta)
if (!is.null(mrdt_master_col)) {
  raw <- as.numeric(merged_meta[[mrdt_master_col]])
  cd$mrdt_result_master <- factor(ifelse(raw==1,"Positive",ifelse(raw==0,"Negative",NA)),
                                  levels=c("Negative","Positive"))
}

sickle_col <- find_col("sicklecell_test_result", merged_meta)
if (!is.null(sickle_col))
  cd$sicklecell_test_result <- recode_factor(as.numeric(merged_meta[[sickle_col]]),
                                             c("1"="HbAA (Normal)","2"="HbAS (Trait)",
                                               "3"="HbSC (Disease)","4"="HbSS (Disease)"), na_codes=c())

sc_col <- find_col("^sample_collected$", merged_meta)
if (!is.null(sc_col))
  cd$sample_collected <- recode_factor(as.numeric(merged_meta[[sc_col]]),
                                       c("1"="Yes","2"="Yes (inadequate)","0"="No"), na_codes=c())

for (ep in c("malaria_episodes","d90_v2_malaria_episodes")) {
  col <- find_col(paste0("^",ep,"$"), merged_meta)
  if (!is.null(col)) {
    ep_map <- c("6"="None","1"="One","2"="Two","3"="Three",
                "4"="Four","5"="More than four")
    cd[[ep]] <- factor(
      recode_factor(as.numeric(merged_meta[[col]]), ep_map, na_codes=c(99)),
      levels=c("None","One","Two","Three","Four","More than four"))
  }
}

doses_col <- find_col("^malaria_doses$", merged_meta)
if (!is.null(doses_col))
  cd$malaria_doses <- recode_factor(as.numeric(merged_meta[[doses_col]]),
                                    c("1"="One","2"="Two","3"="Three","4"="Four"), na_codes=c(99))

for (vacc in c("bcg","opv0","opv1","opv2","opv3","measlesrubella1","measlesrubella2",
               "pcoccal","rtss1","rtss2","rtss3","rtss4")) {
  col <- find_col(paste0("^",vacc,"$"), merged_meta)
  if (!is.null(col)) cd[[vacc]] <- recode_yesno(as.numeric(merged_meta[[col]]))
}

conc_col <- find_col("^concordant$", merged_meta)
if (!is.null(conc_col)) {
  raw <- as.numeric(merged_meta[[conc_col]])
  cd$concordant <- factor(ifelse(raw==1,"Yes",ifelse(raw==0,"No",NA)),
                          levels=c("No","Yes"))
}

for (dc in c("malaria_prevgp","malaria_gp_2017","mrdt_ea_cat")) {
  col <- find_col(dc, merged_meta)
  if (!is.null(col)) cd[[dc]] <- as.character(merged_meta[[col]])
}

## Anthropometrics
muac_col <- find_col(c("muac_zscore","muac_z","^muac_z$"), merged_meta)
whz_col  <- find_col(c("wfh_zscore","whz","whf","weight_for_height"), merged_meta)
haz_col  <- find_col(c("hfa_zscore","haz","hfa","height_for_age"), merged_meta)
waz_col  <- find_col(c("wfa_zscore","waz","wfa","weight_for_age"), merged_meta)

muac_data <- stratify_z(muac_col, merged_meta)
whz_data  <- stratify_z(whz_col,  merged_meta)
haz_data  <- stratify_z(haz_col,  merged_meta)
waz_data  <- stratify_z(waz_col,  merged_meta)

for (zg in list(list("muac_group", muac_data$factor),
                list("whz_group",  whz_data$factor),
                list("haz_group",  haz_data$factor),
                list("waz_group",  waz_data$factor))) {
  cd[[zg[[1]]]] <- factor(zg[[2]], levels=c("Malnourished","Normal","Overweight"))
}

comp <- rep("Normal", ncol(tse_base))
comp[(!is.na(muac_data$raw)&muac_data$raw< -2)|
       (!is.na(whz_data$raw) &whz_data$raw < -2)|
       (!is.na(haz_data$raw) &haz_data$raw < -2)|
       (!is.na(waz_data$raw) &waz_data$raw < -2)] <- "Malnourished"
comp[is.na(muac_data$raw)&is.na(whz_data$raw)&
       is.na(haz_data$raw)&is.na(waz_data$raw)] <- NA
cd$composite_malnutrition <- factor(comp, levels=c("Normal","Malnourished"))

colData(tse_base) <- cd
message("Metadata mapping complete.")

# ==============================================================================
# MICROBIAL ABUNDANCE CO-OCCURRENCE ANALYSIS
# ------------------------------------------------------------------------------
# Estimates pairwise Spearman associations between species-level abundances
# across samples. Positive correlations indicate taxa that tend to be abundant
# together; negative correlations indicate opposing abundance patterns. These
# are associations, not proof of direct ecological interactions, and relative
# abundance data are compositional, so results should be interpreted cautiously.
# ==============================================================================
run_cooccurrence_analysis <- function(tse_obj, out_dir, raw_metaphlan = NULL,
                                      min_prevalence = 0.10,
                                      max_taxa = 40,
                                      min_abs_rho = 0.60,
                                      max_q = 0.05) {
  message("\n--- Co-occurrence: pairwise species abundance associations ---")

  # Prefer the original merged MetaPhlAn table, whose clade_name column
  # preserves full lineage labels (including s__ species rows). The imported
  # TreeSE can expose terminal t__ SGB labels and lose the parent s__ label.
  abund <- NULL
  if (!is.null(raw_metaphlan) && "clade_name" %in% colnames(raw_metaphlan)) {
    sample_cols <- intersect(colnames(tse_obj), colnames(raw_metaphlan))
    clade <- as.character(raw_metaphlan$clade_name)
    # Select rows whose terminal clade is species-level, not strain/SGB-level.
    species_rows <- grepl("(^|\\|)s__[^|]+$", clade)
    if (length(sample_cols) >= 4 && any(species_rows)) {
      abund <- as.matrix(raw_metaphlan[species_rows, sample_cols, drop = FALSE])
      storage.mode(abund) <- "numeric"
      rownames(abund) <- sub(".*\\|s__", "s__", clade[species_rows])
    }
  }

  if (is.null(abund)) {
    message("Co-occurrence skipped: no terminal s__ species rows were found in the original MetaPhlAn table. The script will not substitute t__ SGBs for species.")
    return(invisible(NULL))
  }
  if (nrow(abund) < 2 || ncol(abund) < 4) {
    message("Co-occurrence skipped: species abundance matrix is too small.")
    return(invisible(NULL))
  }

  abund[!is.finite(abund) | abund < 0] <- 0
  sample_totals <- colSums(abund, na.rm = TRUE)
  valid_samples <- is.finite(sample_totals) & sample_totals > 0
  abund <- abund[, valid_samples, drop = FALSE]
  if (ncol(abund) < 4) {
    message("Co-occurrence skipped: fewer than four non-empty samples remain.")
    return(invisible(NULL))
  }

  # Normalize each sample to proportions, independent of the scale of the
  # source abundance assay.
  rel <- sweep(abund, 2, colSums(abund), "/")
  prevalence <- rowMeans(rel > 0, na.rm = TRUE)
  mean_abundance <- rowMeans(rel, na.rm = TRUE)
  keep <- is.finite(prevalence) & prevalence >= min_prevalence &
          is.finite(mean_abundance) & mean_abundance > 0
  rel <- rel[keep, , drop = FALSE]
  mean_abundance <- mean_abundance[keep]
  prevalence <- prevalence[keep]

  if (nrow(rel) < 2) {
    message(sprintf("Co-occurrence skipped: fewer than two taxa passed %.0f%% prevalence.",
                    min_prevalence * 100))
    return(invisible(NULL))
  }

  # Restrict to the most abundant taxa among those passing prevalence so the
  # pairwise analysis remains tractable and the heatmap remains interpretable.
  ord <- order(mean_abundance, decreasing = TRUE)
  chosen <- ord[seq_len(min(length(ord), max_taxa))]
  rel <- rel[chosen, , drop = FALSE]
  mean_abundance <- mean_abundance[chosen]
  prevalence <- prevalence[chosen]

  taxa_names <- rownames(rel)
  if (is.null(taxa_names)) taxa_names <- paste0("Taxon_", seq_len(nrow(rel)))
  # Simplify labels for output and plots while retaining full clade names in CSV.
  display_names <- sub(".*\\|s__", "", taxa_names)
  display_names <- sub("^s__", "", display_names)
  display_names <- gsub("_", " ", display_names, fixed = TRUE)
  display_names[duplicated(display_names)] <- paste0(display_names[duplicated(display_names)], " [", which(duplicated(display_names)), "]")

  pairs <- utils::combn(seq_len(nrow(rel)), 2)
  pair_results <- vector("list", ncol(pairs))
  for (i in seq_len(ncol(pairs))) {
    a <- pairs[1, i]
    b <- pairs[2, i]
    x <- as.numeric(rel[a, ])
    y <- as.numeric(rel[b, ])
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 4 || length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2) {
      rho <- pval <- NA_real_
    } else {
      test <- suppressWarnings(tryCatch(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE),
                                        error = function(e) NULL))
      if (is.null(test)) {
        rho <- pval <- NA_real_
      } else {
        rho <- unname(test$estimate)
        pval <- test$p.value
      }
    }
    pair_results[[i]] <- data.frame(
      taxon_1 = taxa_names[a], taxon_2 = taxa_names[b],
      taxon_1_label = display_names[a], taxon_2_label = display_names[b],
      prevalence_1 = prevalence[a], prevalence_2 = prevalence[b],
      mean_relative_abundance_1 = mean_abundance[a],
      mean_relative_abundance_2 = mean_abundance[b],
      spearman_rho = rho, p_value = pval,
      stringsAsFactors = FALSE
    )
  }
  pairs_df <- do.call(rbind, pair_results)
  pairs_df$q_value <- p.adjust(pairs_df$p_value, method = "BH")
  pairs_df$association <- ifelse(is.na(pairs_df$spearman_rho), NA_character_,
                                 ifelse(pairs_df$spearman_rho > 0, "positive", "negative"))
  pairs_df$significant_edge <- !is.na(pairs_df$q_value) &
    pairs_df$q_value <= max_q & abs(pairs_df$spearman_rho) >= min_abs_rho

  write.csv(pairs_df,
            file = file.path(out_dir, "microbial_cooccurrence_all_pairs.csv"),
            row.names = FALSE)
  edges <- pairs_df[pairs_df$significant_edge, , drop = FALSE]
  edges <- edges[order(edges$q_value, -abs(edges$spearman_rho)), , drop = FALSE]
  write.csv(edges,
            file = file.path(out_dir, "microbial_cooccurrence_significant_edges.csv"),
            row.names = FALSE)

  # Correlation heatmap for the selected taxa; non-significant correlations are
  # greyed out to make the strongest FDR-supported associations easier to see.
  rho_mat <- suppressWarnings(cor(t(rel), method = "spearman", use = "pairwise.complete.obs"))
  p_mat <- matrix(NA_real_, nrow = nrow(rel), ncol = nrow(rel),
                  dimnames = list(taxa_names, taxa_names))
  diag(p_mat) <- 0
  for (i in seq_len(nrow(rel) - 1)) {
    for (j in (i + 1):nrow(rel)) {
      x <- as.numeric(rel[i, ]); y <- as.numeric(rel[j, ])
      test <- suppressWarnings(tryCatch(cor.test(x, y, method = "spearman", exact = FALSE),
                                        error = function(e) NULL))
      if (!is.null(test)) p_mat[i, j] <- p_mat[j, i] <- test$p.value
    }
  }
  q_mat <- matrix(NA_real_, nrow = nrow(rel), ncol = nrow(rel),
                  dimnames = list(taxa_names, taxa_names))
  diag(q_mat) <- 0
  upper <- upper.tri(p_mat)
  q_mat[upper] <- p.adjust(p_mat[upper], method = "BH")
  q_mat[lower.tri(q_mat)] <- t(q_mat)[lower.tri(q_mat)]

  heat_df <- expand.grid(taxon_1 = seq_len(nrow(rel)), taxon_2 = seq_len(nrow(rel)))
  heat_df$rho <- rho_mat[cbind(heat_df$taxon_1, heat_df$taxon_2)]
  heat_df$q <- q_mat[cbind(heat_df$taxon_1, heat_df$taxon_2)]
  heat_df$label_1 <- factor(display_names[heat_df$taxon_1], levels = display_names)
  heat_df$label_2 <- factor(display_names[heat_df$taxon_2], levels = rev(display_names))
  heat_df$plot_rho <- heat_df$rho
  heat_df$plot_rho[is.na(heat_df$q) | heat_df$q > max_q] <- NA_real_

  p_co <- ggplot(heat_df, aes(x = label_1, y = label_2, fill = plot_rho)) +
    geom_tile(colour = "white", linewidth = 0.15) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1), na.value = "grey90",
                         name = "Spearman\nrho") +
    theme_minimal(base_size = 8) +
    theme(axis.title = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 7),
          panel.grid = element_blank(),
          plot.title = element_text(face = "bold")) +
    labs(title = "Microbial abundance co-occurrence",
         subtitle = paste0("Species-level s__ clades passing ", min_prevalence * 100,
                           "% prevalence; top ", nrow(rel), " by mean abundance. ",
                           "Colour shown only for BH-FDR q <= ", max_q,
                           "; |rho| >= ", min_abs_rho,
                           " defines significant edges in the edge-list CSV."))
  ggsave(file.path(out_dir, "microbial_cooccurrence_correlation_heatmap.png"),
         plot = p_co, width = max(8, min(18, 0.32 * nrow(rel) + 5)),
         height = max(7, min(18, 0.32 * nrow(rel) + 5)), dpi = 300,
         limitsize = FALSE)

  message(sprintf("Co-occurrence complete: %d taxa, %d pairs, %d significant edge(s).",
                  nrow(rel), nrow(pairs_df), nrow(edges)))
  message("Outputs: microbial_cooccurrence_all_pairs.csv, microbial_cooccurrence_significant_edges.csv, microbial_cooccurrence_correlation_heatmap.png")
  invisible(list(all_pairs = pairs_df, significant_edges = edges, correlation = rho_mat))
}

run_cooccurrence_analysis(tse_base, output_dir, raw_metaphlan = raw_data)

# ==============================================================================
# 4. RAW ALPHA DIVERSITY TABLE FOR ALL SAMPLES
# ==============================================================================
message("\n--- Step 4: Computing alpha diversity for all samples ---")

tse_base <- addAlpha(tse_base, assay.type="counts", index="Observed",        name="Observed")
tse_base <- addAlpha(tse_base, assay.type="counts", index="shannon",         name="shannon")
tse_base <- addAlpha(tse_base, assay.type="counts", index="inverse_simpson", name="inverse_simpson")
tse_base <- addAlpha(tse_base, assay.type="counts", index="gini_simpson",    name="gini_simpson")
tse_base <- addAlpha(tse_base, assay.type="counts", index="chao1",           name="chao1")

write.csv(as.data.frame(colData(tse_base)),
          file = file.path(output_dir,"alpha_diversity_all_samples.csv"),
          row.names=TRUE)
message("Raw alpha diversity table written.")

# ==============================================================================
# 4b. COMPREHENSIVE PER-SAMPLE DIVERSITY TABLE
# ==============================================================================
message("\n--- Step 4b: Computing per-sample diversity stats table ---")

tse_full <- transformAssay(tse_base, assay.type="counts",
                           method="relabundance", name="rel")
tse_full <- transformAssay(tse_full, assay.type="counts", method="clr",
                           name="clr_assay", pseudocount=1)

alpha_metrics <- c("Observed","shannon","inverse_simpson","gini_simpson","chao1")

diversity_table <- data.frame(
  sample_id    = colnames(tse_full),
  row.names    = colnames(tse_full),
  stringsAsFactors = FALSE
)

for (m in alpha_metrics) {
  vals <- colData(tse_full)[[m]]
  diversity_table[[m]] <- if (!is.null(vals)) as.numeric(vals) else NA_real_
}

beta_specs <- list(
  bray      = list(mat    = t(assay(tse_full,"rel")),
                   dist_fn = function(m) vegan::vegdist(m, method="bray")),
  jaccard   = list(mat    = t(assay(tse_full,"rel")),
                   dist_fn = function(m) vegan::vegdist(m, method="jaccard", binary=TRUE)),
  aitchison = list(mat    = t(assay(tse_full,"clr_assay")),
                   dist_fn = function(m) dist(m, method="euclidean"))
)

for (metric_name in names(beta_specs)) {
  message(paste("  Computing beta coordinates:", metric_name))
  spec        <- beta_specs[[metric_name]]
  dist_matrix <- spec$dist_fn(spec$mat)
  
  # PCoA
  pcoa_res    <- cmdscale(dist_matrix, k=3, eig=TRUE)
  pcoa_coords <- pcoa_res$points
  eig_vals    <- pmax(pcoa_res$eig, 0)
  pct_var     <- round(eig_vals / sum(eig_vals) * 100, 2)
  
  diversity_table[[paste0("PCoA1_",metric_name)]] <- as.numeric(pcoa_coords[,1])
  diversity_table[[paste0("PCoA2_",metric_name)]] <- as.numeric(pcoa_coords[,2])
  diversity_table[[paste0("PCoA3_",metric_name)]] <- as.numeric(pcoa_coords[,3])
  attr(diversity_table, paste0("PCoA_pctvar_",metric_name)) <-
    c(PC1=pct_var[1], PC2=pct_var[2], PC3=pct_var[3])
  
  # NMDS
  nmds_res    <- suppressWarnings(
    vegan::metaMDS(dist_matrix, k=2, trace=FALSE, trymax=100))
  nmds_coords <- nmds_res$points
  
  diversity_table[[paste0("NMDS1_",metric_name)]]       <- as.numeric(nmds_coords[,1])
  diversity_table[[paste0("NMDS2_",metric_name)]]       <- as.numeric(nmds_coords[,2])
  diversity_table[[paste0("NMDS_stress_",metric_name)]] <- round(nmds_res$stress, 6)
}

# Append key metadata columns
meta_cols_to_include <- c(
  "age_group","agegp","gender","age_days","age_months",
  "diarrhoea","excl_brstfed","fever",
  "muac","height","weight","temperature","hemocue_result",
  "muac_group","whz_group","haz_group","waz_group","composite_malnutrition",
  "mrdt_result","mrdt_result_master","malariatest_result",
  "malaria_episodes","malaria_prev",
  "education","economic_step","housing","enough_food","household_size",
  "cluster","location","intensive_area"
)

for (col in meta_cols_to_include) {
  val <- colData(tse_full)[[col]]
  if (!is.null(val)) diversity_table[[col]] <- as.character(val)
}

write.csv(diversity_table,
          file = file.path(output_dir,"diversity_stats_per_sample.csv"),
          row.names=TRUE)
message(paste("Per-sample diversity table written:",
              nrow(diversity_table),"samples x",ncol(diversity_table),"columns"))

# PCoA variance explained companion file
pcoa_var_rows <- list()
for (metric_name in names(beta_specs)) {
  pv <- attr(diversity_table, paste0("PCoA_pctvar_",metric_name))
  if (!is.null(pv))
    pcoa_var_rows[[metric_name]] <- data.frame(
      metric  = metric_name,
      axis    = c("PC1","PC2","PC3"),
      pct_var = as.numeric(pv),
      stringsAsFactors = FALSE
    )
}
write.csv(do.call(rbind, pcoa_var_rows),
          file = file.path(output_dir,"pcoa_variance_explained.csv"),
          row.names=FALSE)
message("PCoA variance explained table written.")

# ==============================================================================
# STEP 4c: DYNAMICALLY EXTRACT ALL CATEGORIES MATCHING THE SINGLE-THREAD SCRIPT
# ==============================================================================
message("\n--- Step 4c: Dynamically reading all categorical targets ---")

full_meta_df <- as.data.frame(colData(tse_base))

all_possible_targets <- colnames(full_meta_df)[sapply(full_meta_df, function(x) is.factor(x) || is.character(x))]

excluded_cols <- c("sample_id", "SampleID", "age_days", "age_months")
targets_categorical <- setdiff(all_possible_targets, excluded_cols)

color_palette <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
  "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a"
)

message(paste("Successfully loaded", length(targets_categorical), "categorical variables for parallel analysis."))

# ==============================================================================
# STEP 5: MAIN ANALYSIS LOOP (SEQUENTIAL)
# ------------------------------------------------------------------------------
# Parallelisation via parallel::makeCluster()/parLapply() has been removed.
# The PSOCK cluster setup was crashing with heap corruption ("malloc(): invalid
# next size", "unsorted double linked list corrupted") when the 16 worker
# processes rendered Cairo-based plots concurrently - this persisted even
# after pre-warming the fontconfig cache, so the concurrency itself (not just
# the cache race) appears to be the trigger. Running everything sequentially
# in a single R process removes that concurrency entirely and trades
# wall-clock time for stability.
#
# Each category still gets: alpha diversity plots, non-parametric alpha
# testing (Wilcoxon/Kruskal-Wallis), beta diversity + PERMANOVA, and
# age-group-stratified sub-analyses - exactly the same work the parallel
# workers did, just one variable at a time.
# ==============================================================================
message("\n--- Step 5: Main analysis loop (Sequential) ---")

if (exists("categorical_targets") && !exists("targets_categorical")) {
  targets_categorical <- categorical_targets
} else if (exists("targets_categorical") && !exists("categorical_targets")) {
  categorical_targets <- targets_categorical
} else if (!exists("targets_categorical") && !exists("categorical_targets")) {
  stop("CRITICAL ERROR: Could not find your target variables list. Check earlier script steps.")
}

if (exists("palette_base") && !exists("color_palette")) {
  color_palette <- palette_base
} else if (!exists("color_palette")) {
  stop("CRITICAL ERROR: Could not find your color palette object (palette_base).")
}

message(paste("Processing", length(targets_categorical),
              "categorical variables sequentially..."))

run_category_analysis <- function(category) {
  age_proxy_vars <- c("age_group", "agegp", "age_days", "age_months")
  
  if (is.null(colData(tse_base)[[category]])) return(NULL)
  
  raw_vec    <- colData(tse_base)[[category]]
  valid_mask <- !is.na(raw_vec) & as.character(raw_vec) != "" & as.character(raw_vec) != "NA"
  
  if (sum(valid_mask) < 6) return(NULL)
  
  tse_clean <- tse_base[, valid_mask]
  colData(tse_clean)[[category]] <- droplevels(factor(colData(tse_clean)[[category]]))
  
  if (nlevels(colData(tse_clean)[[category]]) < 2) return(NULL)
  
  tse_clean <- ensure_alpha(tse_clean)
  
  # A. Generate Alpha Diversity Plots
  current_factors <- colData(tse_clean)[[category]]
  current_colors  <- color_palette[seq_len(nlevels(current_factors))]
  
  tryCatch({
    plot_alpha(
      tse_obj     = tse_clean,
      category    = category,
      color_vec   = current_colors,
      out_path_fn = function(m) file.path(output_dir, paste0("alpha_", m, "_", category, ".png"))
    )
  }, error = function(e) message(sprintf("    plot_alpha failed for '%s': %s",
                                         category, conditionMessage(e))))
  
  # A2. Non-parametric alpha diversity testing
  #     Wilcoxon rank-sum for 2-level variables, Kruskal-Wallis for 3+ levels.
  alpha_rows_main <- tryCatch(test_alpha(tse_clean, category), error = function(e) NULL)
  if (!is.null(alpha_rows_main)) {
    alpha_rows_main$context <- "main"
    alpha_rows_main$label   <- category
  }
  
  # A3. Age-adjusted alpha diversity testing (mirrors the partial/age-adjusted
  #     PERMANOVA used for beta diversity). Adds p_value_age_adjusted per
  #     metric, merged in by category + metric so the unadjusted Wilcoxon/
  #     Kruskal-Wallis results above are left completely untouched.
  if (!is.null(alpha_rows_main)) {
    alpha_rows_adj <- tryCatch(test_alpha_age_adjusted(tse_clean, category),
                               error = function(e) {
                                 message(sprintf("    Age-adjusted alpha test failed for '%s': %s",
                                                 category, conditionMessage(e)))
                                 NULL
                               })
    if (!is.null(alpha_rows_adj)) {
      alpha_rows_main <- merge(alpha_rows_main,
                               alpha_rows_adj[, c("category", "metric", "p_value_age_adjusted")],
                               by = c("category", "metric"), all.x = TRUE, sort = FALSE)
    } else {
      alpha_rows_main$p_value_age_adjusted <- NA_real_
    }
  }
  
  # B. Generate Beta Diversity Calculations & Plots
  cond_vars <- age_proxy_vars[age_proxy_vars != category]
  
  beta_rows <- tryCatch({
    run_beta_diversity(
      tse_obj       = tse_clean,
      category      = category,
      color_palette = color_palette,
      file_prefix   = file.path(output_dir, paste0("beta_", category)),
      condition_on  = cond_vars
    )
  }, error = function(e) {
    message(sprintf("    run_beta_diversity failed for '%s': %s",
                    category, conditionMessage(e)))
    NULL
  })
  
  if (!is.null(beta_rows)) {
    beta_rows$context <- "main"
    beta_rows$label   <- category
  }
  
  # C. Generate Stratified Analysis & Plots (by age_group)
  strat_beta_list  <- list()
  strat_alpha_list <- list()
  
  if (!is.null(colData(tse_clean)$age_group)) {
    for (ag in unique(na.omit(colData(tse_clean)$age_group))) {
      tse_sub_age <- tse_clean[, colData(tse_clean)$age_group == ag & !is.na(colData(tse_clean)$age_group)]
      if (ncol(tse_sub_age) >= 6 && nlevels(droplevels(factor(colData(tse_sub_age)[[category]]))) >= 2) {
        strat_res <- tryCatch({
          run_stratified(
            tse_subset    = tse_sub_age,
            strat_var     = category,
            strat_colors  = current_colors,
            color_palette = color_palette,
            out_dir       = file.path(output_dir, paste0("stratified_age_", ag)),
            label         = paste0(category, "_age_", ag),
            subtitle      = paste("Age Group:", ag)
          )
        }, error = function(e) {
          message(sprintf("    run_stratified failed for '%s' age '%s': %s",
                          category, ag, conditionMessage(e)))
          NULL
        })
        
        if (!is.null(strat_res)) {
          if (!is.null(strat_res$beta)) {
            strat_res$beta$context <- paste0("age_group_stratified:", ag)
            strat_beta_list[[paste0(category, "_age_", ag)]] <- strat_res$beta
          }
          if (!is.null(strat_res$alpha)) {
            strat_res$alpha$context <- paste0("age_group_stratified:", ag)
            strat_alpha_list[[paste0(category, "_age_", ag)]] <- strat_res$alpha
          }
        }
      }
    }
  }
  
  list(
    beta        = beta_rows,
    alpha       = alpha_rows_main,
    strat_beta  = strat_beta_list,
    strat_alpha = strat_alpha_list
  )
}

start_serial_time <- Sys.time()

main_loop_results <- vector("list", length(targets_categorical))
names(main_loop_results) <- targets_categorical

# Per-category timeout (seconds), as a safety net against a category that
# hangs rather than crashes (e.g. a pathological metaMDS convergence loop).
# Typical per-category time observed so far has been under a minute, so this
# leaves a very large margin.
category_timeout_secs <- 1200
crashed_categories <- character(0)

for (i in seq_along(targets_categorical)) {
  category  <- targets_categorical[i]
  cat_start <- Sys.time()
  message(sprintf("\n[%d/%d] Processing category: %s", i, length(targets_categorical), category))
  
  # Run this category's analysis in an isolated forked child process rather
  # than directly in this session. A previous run crashed with heap
  # corruption ("malloc(): invalid next size", SIGABRT/exit 134) during a
  # specific category's PERMANOVA/plotting step even with parallelisation
  # removed entirely - i.e. a native (C/Fortran-level) memory bug somewhere
  # in the compiled graphics/numerical stack, not something an R-level
  # tryCatch() can trap. Forking means that if a child's heap corrupts and
  # it aborts, only that child process dies: this parent process (and the
  # rest of the 83-variable run) is unaffected, since fork() gives the
  # child its own independent copy of memory. Only one child ever runs at a
  # time (we wait for it immediately below), so this does NOT reintroduce
  # the concurrent-rendering issue from the original PSOCK cluster - it is
  # still fully sequential, just crash-isolated.
  job <- tryCatch(
    parallel::mcparallel(run_category_analysis(category), silent = FALSE),
    error = function(e) NULL
  )
  
  cat_result <- NULL
  
  if (is.null(job)) {
    message(sprintf("    Could not fork a worker process for '%s' - skipping.", category))
    crashed_categories <- c(crashed_categories, category)
  } else {
    collected <- tryCatch(
      parallel::mccollect(job, wait = TRUE, timeout = category_timeout_secs),
      error = function(e) {
        message(sprintf("    mccollect error for '%s': %s", category, conditionMessage(e)))
        NULL
      }
    )
    
    if (is.null(collected) || length(collected) == 0) {
      message(sprintf("    Child process for '%s' crashed, timed out, or produced no result - skipping.",
                      category))
      crashed_categories <- c(crashed_categories, category)
      # Make sure a hung/zombie child doesn't linger if it timed out rather
      # than actually crashed.
      tryCatch(tools::pskill(job$pid, signal = tools::SIGKILL), error = function(e) NULL)
    } else {
      cat_result <- collected[[1]]
      if (inherits(cat_result, "try-error")) {
        message(sprintf("    Child process for '%s' returned an error - skipping.", category))
        crashed_categories <- c(crashed_categories, category)
        cat_result <- NULL
      }
    }
  }
  
  main_loop_results[[i]] <- cat_result
  
  cat_elapsed <- round(as.numeric(difftime(Sys.time(), cat_start, units = "secs")), 1)
  message(sprintf("  [%d/%d] Done in %s sec.", i, length(targets_categorical), cat_elapsed))
}

if (length(crashed_categories) > 0) {
  message(paste0("\nWARNING: ", length(crashed_categories),
                 " categor", ifelse(length(crashed_categories) == 1, "y", "ies"),
                 " could not be processed (crashed / timed out / errored) and were skipped:\n  ",
                 paste(crashed_categories, collapse = ", "),
                 "\nThese are candidates for individual follow-up (e.g. running them one at a time",
                 " interactively) since something specific to their data or the underlying compiled",
                 " libraries appears to be triggering the failure."))
  writeLines(crashed_categories, file.path(output_dir, "crashed_categories.txt"))
}

# Reconstruct the output structure - captures beta, alpha, AND the
# age-stratified beta/alpha results from every category.
all_permanova_rows <- list()
all_alpha_results  <- list()

for (i in seq_along(targets_categorical)) {
  cat_name <- targets_categorical[i]
  res      <- main_loop_results[[i]]
  if (is.null(res)) next
  
  if (!is.null(res$beta))  all_permanova_rows[[paste0(cat_name, "_main")]] <- res$beta
  if (!is.null(res$alpha)) all_alpha_results[[paste0(cat_name, "_main")]]  <- res$alpha
  
  if (length(res$strat_beta) > 0) {
    for (nm in names(res$strat_beta))
      all_permanova_rows[[paste0(nm, "_stratbeta")]] <- res$strat_beta[[nm]]
  }
  if (length(res$strat_alpha) > 0) {
    for (nm in names(res$strat_alpha))
      all_alpha_results[[paste0(nm, "_stratalpha")]] <- res$strat_alpha[[nm]]
  }
}

elapsed_total <- round(as.numeric(difftime(Sys.time(), start_serial_time, units="mins")), 1)
message(paste("Sequential (crash-isolated) main loop completed all math and graphics in",
              elapsed_total, "minutes."))

# ==============================================================================
# 6. AGE-GROUP x GENDER AND AGEGP x GENDER STRATIFIED
# ==============================================================================
message("\n\n>>> AGE-GROUP x GENDER STRATIFIED <<<")

gender_colors <- setNames(c("#1B9E77","#D95F02"), c("Female","Male"))
age_groups_all   <- unique(na.omit(as.character(colData(tse_base)$age_group)))
agegp_groups_all <- unique(na.omit(as.character(colData(tse_base)$agegp)))

strat_dir_ag <- file.path(output_dir,"age_gender_stratified")
for (ag in age_groups_all) {
  mask <- !is.na(colData(tse_base)$age_group) &
    colData(tse_base)$age_group == ag &
    !is.na(colData(tse_base)$gender)
  
  if (sum(mask) < 6 ||
      length(unique(as.character(colData(tse_base)$gender[mask]))) < 2) next
  
  ag_label <- gsub("[^A-Za-z0-9]","_",ag)
  message(paste("  age_group:", ag, "| n =", sum(mask)))
  
  tse_ag <- tse_base[, mask]
  colData(tse_ag)$gender <- factor(as.character(colData(tse_ag)$gender),
                                   levels=c("Female","Male"))
  strat_res <- run_stratified(
    tse_subset    = tse_ag,
    strat_var     = "gender",
    strat_colors  = unname(gender_colors),
    color_palette = color_palette,
    out_dir       = strat_dir_ag,
    label         = paste0("agegroup_",ag_label,"_by_gender"),
    subtitle      = paste("age_group:",ag,"| n=",sum(mask))
  )
  if (!is.null(strat_res$beta)) {
    strat_res$beta$context <- paste0("age_group_stratified:",ag)
    strat_res$beta$label   <- paste0("agegroup_",ag_label,"_by_gender")
    all_permanova_rows[[paste0("gender_agegroup_",ag_label)]] <- strat_res$beta
  }
  if (!is.null(strat_res$alpha)) {
    strat_res$alpha$context <- paste0("age_group_stratified:",ag)
    strat_res$alpha$label   <- paste0("agegroup_",ag_label,"_by_gender")
    all_alpha_results[[paste0("gender_agegroup_",ag_label)]] <- strat_res$alpha
  }
}

message("\n\n>>> AGEGP x GENDER STRATIFIED <<<")
strat_dir_agegp <- file.path(output_dir,"agegp_gender_stratified")
for (ag in agegp_groups_all) {
  mask <- !is.na(colData(tse_base)$agegp) &
    colData(tse_base)$agegp == ag &
    !is.na(colData(tse_base)$gender)
  
  if (sum(mask) < 6 ||
      length(unique(as.character(colData(tse_base)$gender[mask]))) < 2) next
  
  ag_label <- gsub("[^A-Za-z0-9]","_",ag)
  message(paste("  agegp:", ag, "| n =", sum(mask)))
  
  tse_ag <- tse_base[, mask]
  colData(tse_ag)$gender <- factor(as.character(colData(tse_ag)$gender),
                                   levels=c("Female","Male"))
  strat_res <- run_stratified(
    tse_subset    = tse_ag,
    strat_var     = "gender",
    strat_colors  = unname(gender_colors),
    color_palette = color_palette,
    out_dir       = strat_dir_agegp,
    label         = paste0("agegp_",ag_label,"_by_gender"),
    subtitle      = paste("agegp:",ag,"| n=",sum(mask))
  )
  if (!is.null(strat_res$beta)) {
    strat_res$beta$context <- paste0("agegp_stratified:",ag)
    strat_res$beta$label   <- paste0("agegp_",ag_label,"_by_gender")
    all_permanova_rows[[paste0("gender_agegp_",ag_label)]] <- strat_res$beta
  }
  if (!is.null(strat_res$alpha)) {
    strat_res$alpha$context <- paste0("agegp_stratified:",ag)
    strat_res$alpha$label   <- paste0("agegp_",ag_label,"_by_gender")
    all_alpha_results[[paste0("gender_agegp_",ag_label)]] <- strat_res$alpha
  }
}

# ==============================================================================
# 7. SIGNIFICANT VARIABLE x AGE_GROUP AND x AGEGP STRATIFIED
# ==============================================================================
message("\n\n>>> SIGNIFICANT-VARIABLE x AGE STRATIFIED <<<")

exclude_from_strat <- c("age_group","agegp","gender")

main_rows <- do.call(rbind, Filter(function(x)
  !is.null(x) && nrow(x) > 0 && "context" %in% names(x) && x$context[1]=="main",
  all_permanova_rows))

if (!is.null(main_rows) && nrow(main_rows) > 0) {
  sig_categories <- unique(main_rows$category[
    main_rows$p_value < 0.05 &
      !main_rows$category %in% exclude_from_strat
  ])
} else {
  sig_categories <- character(0)
}

message(paste("Significant categories:", paste(sig_categories, collapse=", ")))

for (sig_cat in sig_categories) {
  
  # by age_group
  strat_dir_sig_ag <- file.path(output_dir, paste0("agegroup_stratified_",sig_cat))
  for (ag in age_groups_all) {
    age_mask <- !is.na(colData(tse_base)$age_group) &
      colData(tse_base)$age_group == ag
    cat_vec  <- colData(tse_base)[[sig_cat]]
    cat_mask <- !is.na(cat_vec) &
      as.character(cat_vec) != "" &
      as.character(cat_vec) != "NA"
    mask     <- age_mask & cat_mask
    
    if (sum(mask) < 6) next
    if (length(unique(na.omit(as.character(cat_vec[mask])))) < 2) next
    
    ag_label  <- gsub("[^A-Za-z0-9]","_",ag)
    label_str <- paste0(sig_cat,"_agegroup_",ag_label)
    message(paste(" ",sig_cat,"| age_group:",ag,"| n=",sum(mask)))
    
    tse_sub <- tse_base[, mask]
    lvls    <- unique(na.omit(as.character(colData(tse_sub)[[sig_cat]])))
    colors  <- color_palette[seq_along(lvls)]
    
    strat_res <- run_stratified(
      tse_subset    = tse_sub,
      strat_var     = sig_cat,
      strat_colors  = colors,
      color_palette = color_palette,
      out_dir       = strat_dir_sig_ag,
      label         = label_str,
      subtitle      = paste("age_group:",ag,"| n=",sum(mask)),
      condition_on  = "agegp"
    )
    if (!is.null(strat_res$beta)) {
      strat_res$beta$context <- paste0("agegroup_stratified:",ag)
      strat_res$beta$label   <- label_str
      all_permanova_rows[[paste0(sig_cat,"_agegroup_",ag_label)]] <- strat_res$beta
    }
    if (!is.null(strat_res$alpha)) {
      strat_res$alpha$context <- paste0("agegroup_stratified:",ag)
      strat_res$alpha$label   <- label_str
      all_alpha_results[[paste0(sig_cat,"_agegroup_",ag_label)]] <- strat_res$alpha
    }
  }
  
  # by agegp
  strat_dir_sig_agegp <- file.path(output_dir, paste0("agegp_stratified_",sig_cat))
  for (ag in agegp_groups_all) {
    age_mask <- !is.na(colData(tse_base)$agegp) &
      colData(tse_base)$agegp == ag
    cat_vec  <- colData(tse_base)[[sig_cat]]
    cat_mask <- !is.na(cat_vec) &
      as.character(cat_vec) != "" &
      as.character(cat_vec) != "NA"
    mask     <- age_mask & cat_mask
    
    if (sum(mask) < 6) next
    if (length(unique(na.omit(as.character(cat_vec[mask])))) < 2) next
    
    ag_label  <- gsub("[^A-Za-z0-9]","_",ag)
    label_str <- paste0(sig_cat,"_agegp_",ag_label)
    message(paste(" ",sig_cat,"| agegp:",ag,"| n=",sum(mask)))
    
    tse_sub <- tse_base[, mask]
    lvls    <- unique(na.omit(as.character(colData(tse_sub)[[sig_cat]])))
    colors  <- color_palette[seq_along(lvls)]
    
    strat_res <- run_stratified(
      tse_subset    = tse_sub,
      strat_var     = sig_cat,
      strat_colors  = colors,
      color_palette = color_palette,
      out_dir       = strat_dir_sig_agegp,
      label         = label_str,
      subtitle      = paste("agegp:",ag,"| n=",sum(mask)),
      condition_on  = "age_group"
    )
    if (!is.null(strat_res$beta)) {
      strat_res$beta$context <- paste0("agegp_stratified:",ag)
      strat_res$beta$label   <- label_str
      all_permanova_rows[[paste0(sig_cat,"_agegp_",ag_label)]] <- strat_res$beta
    }
    if (!is.null(strat_res$alpha)) {
      strat_res$alpha$context <- paste0("agegp_stratified:",ag)
      strat_res$alpha$label   <- label_str
      all_alpha_results[[paste0(sig_cat,"_agegp_",ag_label)]] <- strat_res$alpha
    }
  }
}

# ==============================================================================
# 8. FDR CORRECTION & FINAL REPORT
# ==============================================================================
message("\n--- Step 8: FDR correction and final report (Beta / PERMANOVA) ---")

all_results <- do.call(rbind, Filter(function(x) !is.null(x) && nrow(x) > 0,
                                     all_permanova_rows))
rownames(all_results) <- NULL

# BH-FDR across all standard p-values jointly
all_results$q_value <- p.adjust(all_results$p_value, method="BH")

# BH-FDR across partial p-values separately
all_results$q_value_partial <- NA_real_
partial_mask <- !is.na(all_results$p_value_partial)
if (any(partial_mask))
  all_results$q_value_partial[partial_mask] <-
  p.adjust(all_results$p_value_partial[partial_mask], method="BH")

all_results$p_sig         <- sig_label(all_results$p_value)
all_results$q_sig         <- sig_label(all_results$q_value)
all_results$p_sig_partial <- sig_label(all_results$p_value_partial)
all_results$q_sig_partial <- sig_label(all_results$q_value_partial)

write.csv(all_results,
          file = file.path(output_dir,"permanova_all_results_with_FDR.csv"),
          row.names=FALSE)

# Human-readable PERMANOVA report
report_file <- file.path(output_dir,"permanova_results_summary.txt")
cat("==========================================================================\n",
    file=report_file)
cat("              MASTER COMPREHENSIVE PERMANOVA REPORT\n",
    file=report_file, append=TRUE)
cat("  FDR: Benjamini-Hochberg | Partial PERMANOVA conditioned on age proxies\n",
    file=report_file, append=TRUE)
cat("==========================================================================\n\n",
    file=report_file, append=TRUE)

for (ctx in unique(all_results$context)) {
  sub <- all_results[all_results$context == ctx, ]
  cat(paste0("\n--- Context: ", ctx, " ---\n"), file=report_file, append=TRUE)
  cat(sprintf("  %-30s %-10s %-6s %-7s %-8s %-4s %-8s %-4s | %-6s %-7s %-8s %-4s %-8s %-4s %-25s\n",
              "Category","Metric","R2","F","p","sig","q","sig",
              "R2adj","Fadj","padj","sig","qadj","sig","Conditioned on"),
      file=report_file, append=TRUE)
  cat(paste(rep("-",145),collapse=""),"\n", file=report_file, append=TRUE)
  
  for (i in seq_len(nrow(sub))) {
    r <- sub[i,]
    cat(sprintf(
      "  %-30s %-10s %-6s %-7s %-8s %-4s %-8s %-4s | %-6s %-7s %-8s %-4s %-8s %-4s %-25s\n",
      r$category, r$metric,
      format(r$R2,      digits=3, nsmall=3),
      format(r$F_stat,  digits=3, nsmall=3),
      format(r$p_value, digits=3, nsmall=3), r$p_sig,
      format(r$q_value, digits=3, nsmall=3), r$q_sig,
      ifelse(is.na(r$R2_partial),      "-", format(r$R2_partial,      digits=3,nsmall=3)),
      ifelse(is.na(r$F_stat_partial),  "-", format(r$F_stat_partial,  digits=3,nsmall=3)),
      ifelse(is.na(r$p_value_partial), "-", format(r$p_value_partial, digits=3,nsmall=3)),
      ifelse(is.na(r$p_sig_partial),   "-", r$p_sig_partial),
      ifelse(is.na(r$q_value_partial), "-", format(r$q_value_partial, digits=3,nsmall=3)),
      ifelse(is.na(r$q_sig_partial),   "-", r$q_sig_partial),
      ifelse(is.na(r$conditioned_on),  "-", r$conditioned_on)),
      file=report_file, append=TRUE)
  }
  cat("\n", file=report_file, append=TRUE)
}

cat("\nLEGEND\n", file=report_file, append=TRUE)
cat("R2           : variance explained by variable (unadjusted)\n",
    file=report_file, append=TRUE)
cat("R2adj        : variance explained AFTER conditioning out age proxies\n",
    file=report_file, append=TRUE)
cat("Pseudo-F     : between-group to within-group variance ratio\n",
    file=report_file, append=TRUE)
cat("p / q        : raw and BH-FDR-adjusted p-values (standard PERMANOVA)\n",
    file=report_file, append=TRUE)
cat("padj / qadj  : raw and BH-FDR-adjusted p-values (age-conditioned partial PERMANOVA)\n",
    file=report_file, append=TRUE)
cat("sig symbols  : *** <=0.001 | ** <=0.01 | * <=0.05 | ns = not significant\n",
    file=report_file, append=TRUE)
cat("\nINTERPRETATION NOTE\n", file=report_file, append=TRUE)
cat("A variable significant in standard PERMANOVA (p/q *) but NOT in partial\n",
    file=report_file, append=TRUE)
cat("PERMANOVA (padj/qadj ns) is likely acting as a proxy for age rather than\n",
    file=report_file, append=TRUE)
cat("having an independent effect on microbiome composition.\n",
    file=report_file, append=TRUE)

# ------------------------------------------------------------------------------
# NEW: Alpha diversity FDR correction and report
# ------------------------------------------------------------------------------
message("\n--- Step 8b: FDR correction and final report (Alpha diversity) ---")

if (length(all_alpha_results) > 0) {
  # Only "main" (unstratified) rows ever get p_value_age_adjusted, since
  # age-adjustment is only run once per category in run_category_analysis()
  # (see requirement 6) - stratified/gender-split alpha rows never attempt
  # it. Backfill the column as NA where absent so every element has an
  # identical column set before row-binding.
  all_alpha_results <- lapply(all_alpha_results, function(df) {
    if (!"p_value_age_adjusted" %in% colnames(df)) df$p_value_age_adjusted <- NA_real_
    df
  })

  all_alpha_df <- do.call(rbind, Filter(function(x) !is.null(x) && nrow(x) > 0,
                                        all_alpha_results))
  rownames(all_alpha_df) <- NULL

  all_alpha_df$q_value <- p.adjust(all_alpha_df$p_value, method = "BH")
  all_alpha_df$p_sig   <- sig_label(all_alpha_df$p_value)
  all_alpha_df$q_sig   <- sig_label(all_alpha_df$q_value)

  # Age-adjusted FDR correction - BH applied only across rows where an
  # age-adjusted p-value was actually computed, mirroring how q_value_partial
  # is handled for the beta-diversity/PERMANOVA results above.
  all_alpha_df$q_value_age_adjusted <- NA_real_
  age_adj_mask <- !is.na(all_alpha_df$p_value_age_adjusted)
  if (any(age_adj_mask)) {
    all_alpha_df$q_value_age_adjusted[age_adj_mask] <-
      p.adjust(all_alpha_df$p_value_age_adjusted[age_adj_mask], method = "BH")
  }
  all_alpha_df$p_sig_age_adjusted <- sig_label(all_alpha_df$p_value_age_adjusted)
  all_alpha_df$q_sig_age_adjusted <- sig_label(all_alpha_df$q_value_age_adjusted)

  # Flag variables whose apparent alpha-diversity association disappears
  # after adjusting for age - i.e. likely age proxies rather than
  # independent effects. The !is.na() guard forces FALSE (rather than NA)
  # when age-adjustment was never attempted for that row (age_group/agegp,
  # stratified contexts, or insufficient data for the model), matching the
  # same pattern already used for the beta-diversity age_proxy flag in
  # plot_permanova_results() above.
  all_alpha_df$age_proxy <- !is.na(all_alpha_df$q_value_age_adjusted) &
    all_alpha_df$q_value < 0.05 &
    all_alpha_df$q_value_age_adjusted >= 0.05

  write.csv(all_alpha_df,
            file = file.path(output_dir, "alpha_diversity_test_results_with_FDR.csv"),
            row.names = FALSE)

  alpha_report_file <- file.path(output_dir, "alpha_diversity_results_summary.txt")
  cat("==========================================================================\n",
      file = alpha_report_file)
  cat("           ALPHA DIVERSITY NON-PARAMETRIC TEST REPORT\n",
      file = alpha_report_file, append = TRUE)
  cat("  2-level variables : Wilcoxon rank-sum test\n",
      file = alpha_report_file, append = TRUE)
  cat("  3+ level variables: Kruskal-Wallis test\n",
      file = alpha_report_file, append = TRUE)
  cat("  Age-adjusted      : rank(alpha_diversity) ~ age_months + variable,\n",
      file = alpha_report_file, append = TRUE)
  cat("                      drop1(model, test='F') on the variable term\n",
      file = alpha_report_file, append = TRUE)
  cat("  FDR: Benjamini-Hochberg, corrected jointly across all tests\n",
      file = alpha_report_file, append = TRUE)
  cat("        (unadjusted and age-adjusted p-values corrected separately)\n",
      file = alpha_report_file, append = TRUE)
  cat("==========================================================================\n\n",
      file = alpha_report_file, append = TRUE)

  for (ctx in unique(all_alpha_df$context)) {
    sub <- all_alpha_df[all_alpha_df$context == ctx, ]
    cat(paste0("\n--- Context: ", ctx, " ---\n"), file = alpha_report_file, append = TRUE)
    cat(sprintf("  %-30s %-16s %-4s %-20s %-9s %-8s %-4s %-8s %-4s | %-8s %-22s %-10s | %-8s %-8s %-4s %-10s\n",
                "Category","Metric","n","Test","Statistic","p","sig","q","sig",
                "Effect","Effect method","Magnitude","p_age","q_age","sig","age_proxy"),
        file = alpha_report_file, append = TRUE)
    cat(paste(rep("-", 195), collapse=""), "\n", file = alpha_report_file, append = TRUE)

    for (i in seq_len(nrow(sub))) {
      r <- sub[i, ]
      cat(sprintf(
        "  %-30s %-16s %-4s %-20s %-9s %-8s %-4s %-8s %-4s | %-8s %-22s %-10s | %-8s %-8s %-4s %-10s\n",
        r$category, r$metric, r$n_total, r$test,
        format(r$statistic, digits = 3, nsmall = 3),
        format(r$p_value,   digits = 3, nsmall = 3), r$p_sig,
        format(r$q_value,   digits = 3, nsmall = 3), r$q_sig,
        format(r$effect_size, digits = 3, nsmall = 3),
        r$effect_method, r$effect_magnitude,
        ifelse(is.na(r$p_value_age_adjusted), "-", format(r$p_value_age_adjusted, digits=3, nsmall=3)),
        ifelse(is.na(r$q_value_age_adjusted), "-", format(r$q_value_age_adjusted, digits=3, nsmall=3)),
        ifelse(is.na(r$p_sig_age_adjusted),   "-", r$p_sig_age_adjusted),
        ifelse(isTRUE(r$age_proxy), "YES", "-")),
        file = alpha_report_file, append = TRUE)
    }
    cat("\n", file = alpha_report_file, append = TRUE)
  }

  cat("\nLEGEND\n", file = alpha_report_file, append = TRUE)
  cat("n_levels     : number of groups compared (2 -> Wilcoxon, 3+ -> Kruskal-Wallis)\n",
      file = alpha_report_file, append = TRUE)
  cat("n            : total samples included in the test (after removing NAs / singleton groups)\n",
      file = alpha_report_file, append = TRUE)
  cat("p / q        : raw and BH-FDR-adjusted p-values (unadjusted test)\n",
      file = alpha_report_file, append = TRUE)
  cat("sig symbols  : *** <=0.001 | ** <=0.01 | * <=0.05 | ns = not significant\n",
      file = alpha_report_file, append = TRUE)
  cat("Effect       : rank-biserial correlation r (Wilcoxon, range -1 to 1, sign shows\n",
      file = alpha_report_file, append = TRUE)
  cat("               which level is higher - see effect_direction column in the CSV)\n",
      file = alpha_report_file, append = TRUE)
  cat("               OR epsilon-squared (Kruskal-Wallis, range 0 to 1)\n",
      file = alpha_report_file, append = TRUE)
  cat("Magnitude    : rank-biserial r -> negligible <0.10 | small <0.30 | medium <0.50 | large >=0.50\n",
      file = alpha_report_file, append = TRUE)
  cat("               epsilon-squared -> negligible <0.01 | small <0.08 | medium <0.26 | large >=0.26\n",
      file = alpha_report_file, append = TRUE)
  cat("p_age / q_age: raw and BH-FDR-adjusted p-values from the age-adjusted model\n",
      file = alpha_report_file, append = TRUE)
  cat("               ('-' = not computed: age_group/agegp, or insufficient data)\n",
      file = alpha_report_file, append = TRUE)
  cat("age_proxy    : YES if significant before age adjustment (q < 0.05) but NOT\n",
      file = alpha_report_file, append = TRUE)
  cat("               after (q_age >= 0.05)\n",
      file = alpha_report_file, append = TRUE)

  cat("\nINTERPRETATION NOTE - AGE-ADJUSTED ALPHA DIVERSITY\n",
      file = alpha_report_file, append = TRUE)
  cat("Age-adjusted p-values are derived from rank-based linear models:\n",
      file = alpha_report_file, append = TRUE)
  cat("\n    rank(alpha_diversity) ~ age_months + variable\n\n",
      file = alpha_report_file, append = TRUE)
  cat("Variables significant before adjustment but not after age adjustment are\n",
      file = alpha_report_file, append = TRUE)
  cat("likely acting as age proxies rather than independently influencing alpha\n",
      file = alpha_report_file, append = TRUE)
  cat("diversity. This mirrors the age-adjusted partial PERMANOVA already used\n",
      file = alpha_report_file, append = TRUE)
  cat("for beta diversity, applied here to the univariate alpha diversity metrics.\n",
      file = alpha_report_file, append = TRUE)
  cat("age_group and agegp are excluded from age-adjustment (p_age/q_age = '-')\n",
      file = alpha_report_file, append = TRUE)
  cat("since adjusting an age variable for age is not meaningful.\n",
      file = alpha_report_file, append = TRUE)

  n_age_proxies <- sum(all_alpha_df$age_proxy, na.rm = TRUE)
  message(paste("Alpha diversity statistical tests (Kruskal-Wallis / Wilcoxon) completed:",
                nrow(all_alpha_df), "rows across", length(unique(all_alpha_df$context)), "contexts."))
  message(paste("Age-adjusted alpha diversity testing flagged", n_age_proxies,
                "metric/category combination(s) as likely age proxies."))

  message("\n--- Generating alpha diversity effect/FDR visualisation plots ---")
  plot_alpha_effect_results(all_alpha_df, out_dir=output_dir, context_filter="main")
} else {
  message("No alpha diversity test results were generated (check group sizes / NA filtering).")
}

# PERMANOVA visualisations
message("\n--- Generating PERMANOVA visualisation plots ---")
plot_permanova_results(all_results, out_dir=output_dir, context_filter="main")

message("\n======================================================================")
message("COMPLETE!")
message(paste("Alpha table (raw)      :", file.path(output_dir,"alpha_diversity_all_samples.csv")))
message(paste("Per-sample table       :", file.path(output_dir,"diversity_stats_per_sample.csv")))
message(paste("PCoA variance table    :", file.path(output_dir,"pcoa_variance_explained.csv")))
message(paste("Beta/PERMANOVA CSV     :", file.path(output_dir,"permanova_all_results_with_FDR.csv")))
message(paste("Beta/PERMANOVA report  :", file.path(output_dir,"permanova_results_summary.txt")))
message(paste("Alpha test CSV         :", file.path(output_dir,"alpha_diversity_test_results_with_FDR.csv")))
message(paste("Alpha test report      :", file.path(output_dir,"alpha_diversity_results_summary.txt")))
message(paste("Alpha effect/FDR plots :", file.path(output_dir,"alpha_diversity_plots/")))
message(paste("PERMANOVA plots        :", file.path(output_dir,"permanova_plots/")))
message("======================================================================")