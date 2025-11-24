# analyze_modality.R
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
  # Strict cleaning: trim to character, drop NA/blank IDs, then pivot
  df <- df_long %>%
    ungroup() %>%
    # keep only shared sample IDs first (as before)
    dplyr::filter(.data[[sample_col]] %in% shared_ids) %>%
    # coerce and trim IDs
    mutate(
      !!sample_col  := trimws(as.character(.data[[sample_col]])),
      !!feature_col := trimws(as.character(.data[[feature_col]]))
    )

  # Drop rows with NA/blank sample/feature IDs
  bad_sample  <- is.na(df[[sample_col]])  | df[[sample_col]]  == ""
  bad_feature <- is.na(df[[feature_col]]) | df[[feature_col]] == ""
  if (any(bad_sample))  message("[make_feature_matrix] Dropping ", sum(bad_sample),  " rows with NA/blank ", sample_col)
  if (any(bad_feature)) message("[make_feature_matrix] Dropping ", sum(bad_feature), " rows with NA/blank ", feature_col)
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

  # Guard against NA/blank rownames after pivot
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

  # Finalize rownames (must be unique, non-empty)
  rownames(wide) <- make.unique(as.character(rn), sep = "_dup")
  wide[[sample_col]] <- NULL
  wide
}

make_drug_matrix <- function(
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
# Main wrapper
# ---------------------------
# Returns: list(drug_mat, feat_mat, shared_ids, cor_tbl, cor_summary, cor_plot, drug_summary)
analyze_modality <- function(
    fits,
    df_long,
    sample_col,        # e.g., "Specimen"
    feature_col,       # e.g., "feature_id" | "Gene" | "site"
    value_col,         # e.g., "correctedAbundance"
    metric = "uM_viability",
    outdir = "figs",
    heatmap_filename = "drug_heatmap_large.pdf",
    fdr_thresh = 0.25
) {
  # Pre-intersection like original (all metrics)
  shared_ids <- base::intersect(unique(fits$improve_sample_id), unique(df_long[[sample_col]]))

  # Feature matrix over shared IDs
  feat_mat <- make_feature_matrix(
    df_long = df_long,
    shared_ids = shared_ids,
    sample_col = sample_col,
    feature_col = feature_col,
    value_col = value_col
  )

  # Drug matrix for the exact metric (no normalization/fallbacks)
  drug_mat <- make_drug_matrix(
    fits        = fits,
    metric      = metric,
    sample_col  = "improve_sample_id",
    drug_col    = "improve_drug_id",
    value_col   = "dose_response_value",
    metric_col  = "dose_response_metric"
  )

  # Summaries & heatmap (original style)
  dsum <- summarize_drugs(
    fits, metric = metric, metric_col = "dose_response_metric", outdir = outdir
  )

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

  # Correlations (only if overlap)
  shared_after <- base::intersect(rownames(drug_mat), rownames(feat_mat))
  cor_tbl <- if (length(shared_after) > 0L) {
    compute_cors(drug_mat, feat_mat, shared_samples = shared_after)
  } else {
    tibble(drug = character(), feature = character(), cor = numeric(),
           pval = numeric(), fdr = numeric(), direction = character())
  }
  cor_res <- summarize_correlated_features(cor_tbl, fdr_thresh = fdr_thresh, outdir = outdir)

  list(
    drug_mat    = drug_mat,
    feat_mat    = feat_mat,
    shared_ids  = shared_ids,
    cor_tbl     = cor_tbl,
    cor_summary = cor_res$summary,
    cor_plot    = cor_res$plot,
    drug_summary = dsum$summary
  )
}
