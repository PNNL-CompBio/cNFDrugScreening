# ---------------------------------------------------------------------------
# 02_analyze_modality_correlations.R
# ---------------------------------------------------------------------------
# Purpose
# - Given (1) drug response fits and (2) a long-format omics table for one modality,
#   this script builds sample x drug and sample x feature matrices, makes a few summary
#   plots (drug efficacy/variability + optional heatmap), and computes Spearman
#   correlations between drug response and molecular features.
#
# Main entry (This does both of the individual calls)
# - analyze_modality(fits, df_long, sample_col, feature_col, value_col, ...)
#
# Individual Calls
# - analyze_drug_response(fits, metric, outdir, heatmap_filename, ...)
#     * Drug-only: builds drug_mat and writes 3 plots by default
#       - most_efficacious.pdf
#       - most_variable.pdf
#       - drug heatmap (heatmap_filename)
# - analyze_modality_correlations(df_long, sample_col, feature_col, value_col, drug_mat, ...)
#     * Modality-only: builds feat_mat and computes correlations + summary plot
#       - cor_features_by_drug.pdf
#
# Inputs
# - fits: long drug response table with improve_sample_id, improve_drug_id,
#         dose_response_metric, dose_response_value
# - df_long: long omics table with sample IDs + feature IDs + values
# - sample_col / feature_col / value_col: column names in df_long that identify
#         the sample, the molecular feature, and the measurement to analyze
#
# Outputs (written to outdir)
# - most_efficacious.pdf, most_variable.pdf
# - drug_heatmap_large.pdf (optional; only for drugs measured in all samples)
# - cor_features_by_drug.pdf (counts of significant correlated features per drug)
#
# Returns (as a list)
# - drug_mat, feat_mat: wide matrices used for analysis
# - cor_tbl: per drug-feature correlations (Spearman) + p-values + FDR
# - cor_summary / cor_plot: summary of significant correlations
# - drug_summary: per-drug mean response, # measured, and variability
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
})

dir.create("figs", showWarnings = FALSE)


# Helpers

make_feature_matrix <- function(df_long, shared_ids, sample_col, feature_col, value_col) {
  #   Build a sample × feature matrix from a long-format omics table, restricted to a set
  #   of shared sample IDs. Cleans IDs, drops blank IDs, pivots to wide, fills missing with 0,
  #   and averages duplicates (mean) per sample-feature pair.
  # Inputs:
  #   df_long: long omics data.frame
  #   shared_ids: character vector of sample IDs to keep
  #   sample_col: column name in df_long for sample ID
  #   feature_col: column name in df_long for feature ID
  #   value_col: column name in df_long for numeric value
  # Output:
  #   data.frame (wide) with rownames = samples and columns = features
  df <- df_long %>%
    ungroup() %>%
    dplyr::filter(.data[[sample_col]] %in% shared_ids) %>%
    mutate(
      !!sample_col  := trimws(as.character(.data[[sample_col]])),
      !!feature_col := trimws(as.character(.data[[feature_col]]))
    )

  # Drop rows with NA/blank sample/feature IDs
  bad_sample  <- is.na(df[[sample_col]])  | df[[sample_col]]  == ""
  bad_feature <- is.na(df[[feature_col]]) | df[[feature_col]] == ""
  df <- df[!(bad_sample | bad_feature), , drop = FALSE]

  if (!nrow(df)) {
    warning("[make_feature_matrix] No rows left after cleaning; returning empty frame.")
    out <- data.frame(check.names = FALSE)
    return(out)
  }

  # Pivot to wide
  wide <- df %>%
    dplyr::select(all_of(c(sample_col, feature_col, value_col))) %>%
    tidyr::pivot_wider(
      names_from  = all_of(feature_col),
      values_from = all_of(value_col),
      values_fill = 0,
      values_fn   = mean
    ) %>%
    as.data.frame(check.names = FALSE)

  rn <- wide[[sample_col]]
  bad_rn <- is.na(rn) | rn == ""
  if (any(bad_rn)) {
    message("[make_feature_matrix] Removing ", sum(bad_rn), " rows with NA/blank rownames after pivot.")
    wide <- wide[!bad_rn, , drop = FALSE]
    rn   <- rn[!bad_rn]
  }

  if (!nrow(wide)) {
    warning("[make_feature_matrix] Wide table is empty after removing bad rownames; returning empty frame.")
    out <- data.frame(check.names = FALSE)
    return(out)
  }

  # Finalize rownames
  rownames(wide) <- make.unique(as.character(rn), sep = "_dup")
  wide[[sample_col]] <- NULL
  wide
}

make_drug_matrix <- function(
    #   Build a sample × drug response matrix from the long drug fits table for a chosen metric
  #   (e.g., uM_viability). Pivots to wide and averages duplicates (mean).
  # Inputs:
  #   fits: long drug response data.frame
  #   metric: metric value to select from metric_col (default "uM_viability")
  #   sample_col: sample ID column name in fits (default "improve_sample_id")
  #   drug_col: drug ID column name in fits (default "improve_drug_id")
  #   value_col: response value column name in fits (default "dose_response_value")
  #   metric_col: metric label column name in fits (default "dose_response_metric")
  # Output:
  #   data.frame (wide) with rownames = samples and columns = drugs
  fits, metric = "uM_viability",
  sample_col = "improve_sample_id",
  drug_col   = "improve_drug_id",
  value_col  = "dose_response_value",
  metric_col = "dose_response_metric"
) {
  fits %>%
    dplyr::filter(.data[[metric_col]] == metric) %>%
    dplyr::select(all_of(c(sample_col, drug_col, value_col))) %>%
    tidyr::pivot_wider(
      names_from  = all_of(drug_col),
      values_from = all_of(value_col),
      values_fn   = mean
    ) %>%
    tibble::column_to_rownames(sample_col)
}

summarize_drugs <- function(
    #   Summarize per-drug response for a selected metric and write two PDF scatter plots:
  #   (1) "most_efficacious" (low mean viability) and (2) "most_variable" (high SD).
  # Inputs:
  #   fits: long drug response data.frame
  #   metric: metric to analyze (default "uM_viability")
  #   metric_col: column holding metric labels (default "dose_response_metric")
  #   outdir: directory to write PDFs (default "figs")
  #   rotate_x: x-axis label rotation angle for readability
  # Output:
  #   list with:
  #     summary: data.frame of per-drug meanResponse, nMeasured, variability
  #     p_eff: ggplot object for efficacious drugs
  #     p_var: ggplot object for variable drugs
  fits, metric = "uM_viability", metric_col = "dose_response_metric",
  outdir = "figs", rotate_x = 45
) {
  ds <- fits %>%
    dplyr::filter(.data[[metric_col]] == metric) %>%
    group_by(.data$improve_drug_id) %>%
    distinct() %>%
    summarize(
      meanResponse = mean(.data$dose_response_value, na.rm = TRUE),
      nMeasured    = n_distinct(.data$improve_sample_id),
      variability  = sd(.data$dose_response_value, na.rm = TRUE),
      .groups = "drop"
    )

  p_eff <- ds %>%
    arrange(desc(.data$meanResponse)) %>%
    dplyr::filter(.data$meanResponse < 0.5) %>%
    ggplot(aes(y = .data$meanResponse, x = .data$improve_drug_id,
               colour = .data$nMeasured, size = .data$variability)) +
    geom_point() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = rotate_x, hjust = 1)) +
    labs(title = "Most efficacious drugs",
         y = "Mean cell viability (fraction)", x = "Drug")

  p_var <- ds %>%
    arrange(desc(.data$variability)) %>%
    dplyr::filter(.data$variability > 0.15) %>%
    ggplot(aes(y = .data$meanResponse, x = .data$improve_drug_id,
               colour = .data$nMeasured, size = .data$variability)) +
    geom_point() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = rotate_x, hjust = 1)) +
    labs(title = "Most variable drugs",
         y = "Mean cell viability (fraction)", x = "Drug")

  ggsave(file.path(outdir, "most_efficacious.pdf"), p_eff, width = 12, height = 8, dpi = 300)
  ggsave(file.path(outdir, "most_variable.pdf"),   p_var, width = 12, height = 8, dpi = 300)

  list(summary = ds, p_eff = p_eff, p_var = p_var)
}

compute_cors <- function(drug_mat, feat_mat, shared_samples = NULL) {
  #   Compute Spearman correlations between each drug response column and each feature column
  #   across shared samples. Also computes per-pair p-values (cor.test) when enough data exists
  #   and applies BH FDR correction.
  # Inputs:
  #   drug_mat: numeric matrix/data.frame (samples × drugs), rownames = sample IDs
  #   feat_mat: numeric matrix/data.frame (samples × features), rownames = sample IDs
  #   shared_samples: optional character vector of sample IDs to use; if NULL, uses rowname intersection
  # Output:
  #   tibble/data.frame with columns: drug, feature, cor, pval, fdr, direction
  if (is.null(shared_samples)) {
    shared_samples <- base::intersect(rownames(drug_mat), rownames(feat_mat))
  }
  if (length(shared_samples) == 0L) {
    return(tibble(
      drug = character(), feature = character(),
      cor = numeric(), pval = numeric(), fdr = numeric(),
      direction = character()
    ))
  }

  drug_mat <- drug_mat[shared_samples, , drop = FALSE]
  feat_mat <- feat_mat[shared_samples, , drop = FALSE]

  cres <- suppressWarnings(
    stats::cor(drug_mat, feat_mat, use = "pairwise.complete.obs", method = "spearman")
  ) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("drug") %>%
    tidyr::pivot_longer(cols = - "drug", names_to = "feature", values_to = "cor")

  csig <- do.call(rbind, lapply(colnames(drug_mat), function(d) {
    do.call(rbind, lapply(colnames(feat_mat), function(f) {
      dv <- drug_mat[, d]; fv <- feat_mat[, f]
      p <- NA_real_
      if (sum(is.finite(dv) & is.finite(fv)) >= 3) {
        p <- tryCatch(
          stats::cor.test(dv, fv, method = "spearman", use = "pairwise.complete.obs")$p.value,
          error = function(e) NA_real_
        )
      }
      c(drug = d, feature = f, pval = p)
    })) %>%
      as.data.frame()
  })) %>%
    as.data.frame() %>%
    mutate(pval = as.numeric(.data$pval)) %>%
    mutate(fdr  = p.adjust(.data$pval, method = "BH"))

  left_join(cres, csig, by = c("drug","feature")) %>%
    mutate(direction = ifelse(.data$cor < 0, "neg", "pos"))
}

summarize_correlated_features <- function(cor_tbl, fdr_thresh = 0.25, outdir = "figs") {
  #   Summarize significant drug-feature associations by counting how many features are
  #   significantly correlated with each drug (split by positive/negative direction),
  #   and write a bar plot PDF.
  # Inputs:
  #   cor_tbl: correlation table from compute_cors()
  #   fdr_thresh: significance threshold on FDR (default 0.25)
  #   outdir: directory to write the PDF (default "figs")
  # Output:
  #   list with:
  #     summary: tibble of per-drug counts and mean correlation by direction
  #     plot: ggplot object (or NULL if no significant results)
  if (nrow(cor_tbl) == 0L) return(list(summary = tibble(), plot = NULL))
  corsummary <- cor_tbl %>%
    dplyr::filter(is.finite(.data$fdr), !is.na(.data$fdr), .data$fdr < fdr_thresh) %>%
    mutate(direction = ifelse(.data$cor > 0, "pos", "neg")) %>%
    group_by(.data$drug, .data$direction) %>%
    summarize(features = n(), meanCor = mean(.data$cor), .groups = "drop")

  if (nrow(corsummary) == 0L) return(list(summary = corsummary, plot = NULL))

  p <- corsummary %>%
    dplyr::filter(.data$features > 1) %>%
    ggplot(aes(x = .data$drug, y = .data$features, fill = .data$direction)) +
    geom_col(position = "dodge") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = paste0("Significant feature counts per drug (FDR < ", fdr_thresh, ")"),
         x = "Drug", y = "# Features")

  ggsave(file.path(outdir, "cor_features_by_drug.pdf"), p, width = 12, height = 6, units = "in")
  list(summary = corsummary, plot = p)
}


# ---------------------------
# Drug-only analysis stage
# ---------------------------
analyze_drug_response <- function(
    #   Drug-only workflow:
  #   - builds sample x drug matrix for a metric
  #   - writes drug summary plots (most_efficacious, most_variable)
  #   - writes drug heatmap (by default) for drugs measured in all samples
  # Inputs:
  #   fits: long drug response table
  #   metric: drug response metric to analyze (default "uM_viability")
  #   outdir: output directory for plots (default "figs")
  #   heatmap_filename: filename for drug heatmap PDF; set NULL to skip
  # Output:
  #   list containing:
  #     drug_mat, drug_summary, p_eff, p_var
  fits,
  metric = "uM_viability",
  outdir = "figs",
  heatmap_filename = "drug_heatmap_large.pdf"
) {

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # Drug matrix for the metric
  drug_mat <- make_drug_matrix(
    fits        = fits,
    metric      = metric,
    sample_col  = "improve_sample_id",
    drug_col    = "improve_drug_id",
    value_col   = "dose_response_value",
    metric_col  = "dose_response_metric"
  )

  # Summaries (writes most_efficacious.pdf + most_variable.pdf)
  dsum <- summarize_drugs(
    fits, metric = metric, metric_col = "dose_response_metric", outdir = outdir
  )

  # Heatmap (by default)
  if (!is.null(heatmap_filename) && nrow(drug_mat) > 0 && ncol(drug_mat) > 0) {
    fulldrugs <- dsum$summary %>%
      dplyr::filter(.data$nMeasured == nrow(drug_mat)) %>%
      pull(.data$improve_drug_id)

    subm <- drug_mat[, colnames(drug_mat) %in% fulldrugs, drop = FALSE]
    if (nrow(subm) > 1 && ncol(subm) > 0) {
      pheatmap::pheatmap(
        as.matrix(subm),
        filename     = file.path(outdir, heatmap_filename),
        width        = 28, height = 16,
        angle_col    = 45, fontsize_col = 6,
        cluster_rows = TRUE, cluster_cols = TRUE,
        show_rownames = TRUE, show_colnames = TRUE
      )
    }
  }

  list(
    drug_mat     = drug_mat,
    drug_summary = dsum$summary,
    p_eff        = dsum$p_eff,
    p_var        = dsum$p_var
  )
}


# ---------------------------
# Modality-only analysis stage
# ---------------------------
analyze_modality_correlations <- function(
    #   Modality-only workflow:
  #   - aligns samples shared between drug_mat and omics table
  #   - builds sample x feature matrix
  #   - computes drug-feature Spearman correlations + p-values + FDR
  #   - summarizes significant features per drug and writes a summary plot
  # Inputs:
  #   df_long: long omics table for one modality
  #   sample_col: sample ID column name in df_long (e.g., "Specimen")
  #   feature_col: feature ID column name in df_long (e.g., "feature_id" or "Gene")
  #   value_col: numeric value column name in df_long (e.g., "correctedAbundance")
  #   drug_mat: sample x drug matrix (rownames = sample IDs)
  #   outdir: output directory for plots (default "figs")
  #   fdr_thresh: FDR cutoff used for correlation summary (default 0.25)
  # Output:
  #   list containing feat_mat, cor_tbl, cor_summary, cor_plot, shared_ids
  df_long,
  sample_col,
  feature_col,
  value_col,
  drug_mat,
  outdir = "figs",
  fdr_thresh = 0.25
) {

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  shared_ids <- base::intersect(rownames(drug_mat), unique(df_long[[sample_col]]))

  feat_mat <- make_feature_matrix(
    df_long = df_long,
    shared_ids = shared_ids,
    sample_col = sample_col,
    feature_col = feature_col,
    value_col = value_col
  )

  # Correlations
  shared_after <- base::intersect(rownames(drug_mat), rownames(feat_mat))
  cor_tbl <- if (length(shared_after) > 0L) {
    compute_cors(drug_mat, feat_mat, shared_samples = shared_after)
  } else {
    tibble(drug = character(), feature = character(), cor = numeric(),
           pval = numeric(), fdr = numeric(), direction = character())
  }
  cor_res <- summarize_correlated_features(cor_tbl, fdr_thresh = fdr_thresh, outdir = outdir)

  list(
    feat_mat    = feat_mat,
    shared_ids  = shared_ids,
    cor_tbl     = cor_tbl,
    cor_summary = cor_res$summary,
    cor_plot    = cor_res$plot
  )
}


# ---------------------------
# Main wrapper (backwards compatible)
# ---------------------------
analyze_modality <- function(
    #   End-to-end wrapper for running one omics modality:
  #   - runs drug-only analysis (drug_mat + drug plots + heatmap)
  #   - runs modality-only correlations (feat_mat + cor_tbl + summary plot)
  # Inputs:
  #   fits: long drug response table (must include improve_sample_id, improve_drug_id, dose_response_metric, dose_response_value)
  #   df_long: long omics table for one modality
  #   sample_col: sample ID column name in df_long (e.g., "Specimen")
  #   feature_col: feature ID column name in df_long (e.g., "feature_id" or "Gene")
  #   value_col: numeric value column name in df_long (e.g., "correctedAbundance")
  #   metric: drug response metric to analyze (default "uM_viability")
  #   outdir: output directory for plots (default "figs")
  #   heatmap_filename: filename for drug heatmap PDF; set NULL to skip
  #   fdr_thresh: FDR cutoff used for correlation summary (default 0.25)
  # Output:
  #   list containing matrices, correlation results, summaries, and (optionally) plot objects:
  #     drug_mat, feat_mat, shared_ids, cor_tbl, cor_summary, cor_plot, drug_summary
  fits,
  df_long,
  sample_col,        # e.g., "Specimen"
  feature_col,       # e.g., "feature_id" | "Gene" | "site"
  value_col,         # e.g., "correctedAbundance"
  metric = "uM_viability",    # Or fit_auc
  outdir = "figs",
  heatmap_filename = "drug_heatmap_large.pdf",
  fdr_thresh = 0.25
) {

  drug_res <- analyze_drug_response(
    fits = fits,
    metric = metric,
    outdir = outdir,
    heatmap_filename = heatmap_filename
  )

  mod_res <- analyze_modality_correlations(
    df_long = df_long,
    sample_col = sample_col,
    feature_col = feature_col,
    value_col = value_col,
    drug_mat = drug_res$drug_mat,
    outdir = outdir,
    fdr_thresh = fdr_thresh
  )

  list(
    drug_mat     = drug_res$drug_mat,
    feat_mat     = mod_res$feat_mat,
    shared_ids   = mod_res$shared_ids,
    cor_tbl      = mod_res$cor_tbl,
    cor_summary  = mod_res$cor_summary,
    cor_plot     = mod_res$cor_plot,
    drug_summary = drug_res$drug_summary
  )
}
