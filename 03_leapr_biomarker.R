# ---------------------------------------------------------------------------
# 03_leapr_biomarker.R
# ---------------------------------------------------------------------------
# Purpose
# - For each drug, correlate drug response (uM_viability) with omics features across
#   samples, then use leapR to find enriched pathways/genesets among:
#     * TOP    = features positively correlated with viability (more resistant)
#     * BOTTOM = features negatively correlated with viability (more sensitive)
#
# Main entry
# - run_leapr_directional_one_cached(drugs, df_long, sample_col, feature_col, value_col, omic_label, cache_path, ...)
#
# Inputs
# - drugs: long drug-response table with improve_drug_id, improve_sample_id,
#          dose_response_metric, dose_response_value
# - df_long: long omics table (sample/feature/value columns)
# - sample_col / feature_col / value_col: column names in df_long that identify
#          the sample ID, feature ID, and numeric measurement
# - omic_label: label used for reporting/assay naming (e.g., "global", "rna", "phospho")
# - cache_path: .RData file path used to save/load results (skips recompute unless always_rerun=TRUE)
#
# Options
# - geneset_name / geneset_object: choose the leapR geneset DB (defaults depend on omic_label)
# - min_features: minimum number of correlated features required to run leapR for TOP/BOTTOM
# - write_csvs: write per-drug leapR tables to CSV
# - test_one: run only the first drug - make sure things are actually working.
#
# Outputs
# - Cached results saved to cache_path (if provided)
# - Optional CSVs: leapR_top_paths/dir_split/*_{TOP|BOTTOM}.csv (if write_csvs=TRUE)
# - Plots: save_leapr_plots() writes pathway barplots to figs/pathways_*.pdf
#
# Returns
# - Named list by drug: res_list[[drug]]$top and res_list[[drug]]$bottom (leapR result tables)
# ---------------------------------------------------------------------------


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(leapR)
  library(ggplot2)
  library(grDevices)
})

# -----------------------------
# Helpers
# -----------------------------
long_to_matrix <- function(df_long, sample_col, feature_col, value_col) {
  #   Convert a long-format omics table into a numeric matrix (rows = samples, cols = features).
  #   Cleans sample/feature IDs, drops NA/blank IDs, pivots to wide, fills missing with 0,
  #   and averages duplicates (mean) per sample-feature pair.
  # Inputs:
  #   df_long: long omics data.frame
  #   sample_col: column name in df_long for sample IDs
  #   feature_col: column name in df_long for feature IDs
  #   value_col: column name in df_long for numeric values
  # Output:
  #   numeric matrix with rownames = sample IDs and colnames = feature IDs, or NULL if empty
  if (is.null(df_long) || !nrow(df_long)) return(NULL)

  df <- df_long |>
    dplyr::mutate(
      !!sample_col  := trimws(as.character(.data[[sample_col]])),
      !!feature_col := trimws(as.character(.data[[feature_col]]))
    )

  # Drop NA/blank sample or feature IDs with messages
  bad_sample  <- is.na(df[[sample_col]])  | df[[sample_col]]  == ""
  bad_feature <- is.na(df[[feature_col]]) | df[[feature_col]] == ""
  n_bad_s     <- sum(bad_sample,  na.rm = TRUE)
  n_bad_f     <- sum(bad_feature, na.rm = TRUE)
  if (n_bad_s > 0) message("[long_to_matrix] Dropping ", n_bad_s, " rows with NA/blank ", sample_col)
  if (n_bad_f > 0) message("[long_to_matrix] Dropping ", n_bad_f, " rows with NA/blank ", feature_col)
  df <- df[!(bad_sample | bad_feature), , drop = FALSE]

  if (!nrow(df)) {
    warning("[long_to_matrix] No rows left after removing NA/blank sample/feature IDs.")
    return(NULL)
  }

  # Pivot to wide
  wide <- df %>%
    dplyr::select(
      !!rlang::sym(sample_col),
      !!rlang::sym(feature_col),
      !!rlang::sym(value_col)
    ) %>%
    tidyr::pivot_wider(
      names_from  = !!rlang::sym(feature_col),
      values_from = !!rlang::sym(value_col),
      values_fill = 0,
      values_fn   = mean
    ) %>%
    as.data.frame(check.names = FALSE)

  rn <- wide[[sample_col]]
  bad_rn <- is.na(rn) | rn == ""
  if (any(bad_rn)) {
    message("[long_to_matrix] Removing ", sum(bad_rn), " rows with NA/blank rownames after pivot.")
    wide <- wide[!bad_rn, , drop = FALSE]
    rn   <- rn[!bad_rn]
  }
  if (!nrow(wide)) {
    warning("[long_to_matrix] Wide table is empty after cleaning.")
    return(NULL)
  }

  rownames(wide) <- make.unique(as.character(rn), sep = "_dup")
  wide[[sample_col]] <- NULL
  as.matrix(wide)
}



# ---- PHOSPHO
.extract_gene_from_site <- function(site_id) {
  #   Extract a gene symbol from a phosphosite/site identifier string. Primarily uses the
  #   substring before the first '-' (e.g., "AAAS-S495s" -> "AAAS"); falls back to splitting
  #   on common delimiters if needed.
  # Inputs:
  #   site_id: character scalar (site/feature ID)
  # Output:
  #   character scalar gene symbol (uppercase) or NA if not parseable
  if (is.na(site_id) || site_id == "") return(NA_character_)
  x <- as.character(site_id)

  # Get chars before first '-' (e.g., "AAAS-S495s" -> "AAAS")
  gene <- sub("^([^\\-]+)-.*$", "\\1", x, perl = TRUE)

  # If no '-' present, fall back to splitting on common delimiters
  if (identical(gene, x)) {
    parts <- strsplit(x, "[|:_\\-\\.]", fixed = FALSE)[[1]]
    gene  <- parts[1]
  }
  gene <- sub("^([A-Za-z0-9]+).*", "\\1", gene)
  gene <- toupper(gene)
  if (nchar(gene) == 0) return(NA_character_)
  gene
}

# Build phospho site gene map from long table
.build_phospho_gene_map_from_long <- function(df_long, feature_col) {
  #   Build a mapping from phosphosite IDs to gene symbols using columns present in the
  #   long table (preferred). If no suitable gene column exists, falls back to parsing gene
  #   symbols from the site IDs themselves.
  # Inputs:
  #   df_long: long omics data.frame
  #   feature_col: column name in df_long containing phosphosite IDs
  # Output:
  #   named character vector mapping site -> gene, or NULL if no sites found
  gene_cols <- c("Gene","gene","hgnc_id","hgnc_symbol","protein","Protein","Symbol","symbol")
  has <- gene_cols[gene_cols %in% colnames(df_long)]
  if (length(has)) {
    gcol <- has[[1]]
    mp <- df_long %>%
      dplyr::select(!!rlang::sym(feature_col), !!rlang::sym(gcol)) %>%
      dplyr::rename(site = !!rlang::sym(feature_col), gene = !!rlang::sym(gcol)) %>%
      dplyr::mutate(site = trimws(as.character(site)),
                    gene = toupper(trimws(as.character(gene)))) %>%
      dplyr::filter(!is.na(site), site != "", !is.na(gene), gene != "") %>%
      dplyr::distinct(site, gene)
    if (nrow(mp)) return(setNames(mp$gene, mp$site))
  }
  sites <- unique(trimws(as.character(df_long[[feature_col]])))
  sites <- sites[!is.na(sites) & sites != ""]
  if (!length(sites)) return(NULL)
  genes <- vapply(sites, .extract_gene_from_site, FUN.VALUE = character(1))
  genes[genes == ""] <- NA_character_
  setNames(genes, sites)
}


.collapse_sites_to_genes <- function(cor_named_vec, map_site2gene, agg = c("mean","maxabs")) {
  #   Collapse a named vector of site-level correlation values to gene-level values using a
  #   site->gene mapping. Supports aggregation by mean or by the max absolute correlation per gene.
  # Inputs:
  #   cor_named_vec: named numeric vector (names = site IDs, values = correlations)
  #   map_site2gene: named character vector mapping site -> gene
  #   agg: "mean" or "maxabs" (how to aggregate multiple sites per gene)
  # Output:
  #   named numeric vector (names = genes, values = aggregated correlations)
  agg <- match.arg(agg)
  if (is.null(map_site2gene) || !length(cor_named_vec)) return(cor_named_vec)

  # Align and drop unmapped
  genes <- map_site2gene[names(cor_named_vec)]
  keep  <- !is.na(genes) & genes != ""
  v     <- cor_named_vec[keep]
  g     <- genes[keep]
  if (!length(v)) return(setNames(numeric(0), character(0)))

  if (agg == "mean") {
    # mean per gene
    df <- tibble(gene = g, val = as.numeric(v)) %>%
      group_by(gene) %>% summarise(val = mean(val, na.rm = TRUE), .groups = "drop")
    out <- stats::setNames(df$val, df$gene)
  } else {
    # max by absolute value, keep sign
    df <- tibble(gene = g, val = as.numeric(v)) %>%
      mutate(ord = order(-abs(val))) %>%
      group_by(gene) %>%
      slice_max(order_by = abs(val), n = 1, with_ties = FALSE) %>%
      ungroup()
    out <- stats::setNames(df$val, df$gene)
  }
  out
}

# ---- Normalize phosphosite IDs to match kinasesubstrates (e.g. "AAAS-S495s" -> "AAAS-S495")
.normalize_kinase_site_id <- function(x) {
  #   Normalize phosphosite IDs to better match leapR kinasesubstrates formatting by trimming
  #   trailing lowercase letters (e.g., "AAAS-S495s" -> "AAAS-S495").
  # Inputs:
  #   x: character vector of phosphosite IDs
  # Output:
  #   character vector of normalized site IDs
  x <- as.character(x)
  x <- trimws(x)
  # Drop trailing lowercase letters
  sub("[a-z]+$", "", x)
}

# Spearman correlations
.col_spearman <- function(vec, mat) {
  #   Compute Spearman correlation between a drug response vector and each feature column in a
  #   sample × feature matrix. Uses sample ID intersection and pairwise complete observations.
  # Inputs:
  #   vec: named numeric vector of responses (names = sample IDs)
  #   mat: numeric matrix (rownames = sample IDs, cols = features)
  # Output:
  #   named numeric vector of correlations (one per feature column; NA where not computable)
  shared <- intersect(names(vec), rownames(mat))
  if (length(shared) < 3) return(setNames(rep(NA_real_, ncol(mat)), colnames(mat)))
  v <- vec[shared]
  m <- as.matrix(mat[shared, , drop = FALSE])
  apply(m, 2, function(col) {
    if (all(is.na(col))) return(NA_real_)
    if (sd(col, na.rm = TRUE) == 0 || sd(v, na.rm = TRUE) == 0) return(NA_real_)
    suppressWarnings(cor(v, col, method = "spearman", use = "pairwise.complete.obs"))
  })
}

# Build SummarizedExperiment to feed into leapR
.build_se_from_corvec <- function(cor_named_vec, features_all, col_label,
                                  map_to_gene = NULL, assay_label = "proteomics") {
  #   Build a single-column SummarizedExperiment containing correlation scores for a set of features,
  #   suitable as input to leapR enrichment functions. Optionally stores a mapped gene ID in rowData.
  # Inputs:
  #   cor_named_vec: named numeric vector of scores (names = feature IDs)
  #   features_all: character vector of features to include (sets row order and rownames)
  #   col_label: column/sample label to assign in the SE (e.g., "<drug>_TOP")
  #   map_to_gene: optional named vector mapping feature -> gene ID/symbol (stored as hgnc_id)
  #   assay_label: assay name label to assign (e.g., "proteomics", "phospho", "rna")
  # Output:
  #   SummarizedExperiment with 1 assay column holding the scores
  v <- rep(NA_real_, length(features_all)); names(v) <- features_all
  common <- intersect(names(cor_named_vec), features_all)
  v[common] <- cor_named_vec[common]
  mat <- matrix(v, nrow = length(v), ncol = 1, dimnames = list(features_all, col_label))
  rd  <- S4Vectors::DataFrame(feature_id = features_all)
  rd$hgnc_id <- if (is.null(map_to_gene)) features_all else map_to_gene[features_all]
  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(values = mat),
    rowData = rd,
    colData = S4Vectors::DataFrame(sample = col_label)
  )
  SummarizedExperiment::assayNames(se) <- assay_label
  se
}

.safe_leapr <- function(...) {
  #   Run leapR::leapR() safely: catches errors, prints a readable message, and returns NULL
  #   instead of stopping the whole pipeline. This is used for Debugging.
  # Inputs:
  #   ...: arguments passed directly to leapR::leapR()
  # Output:
  #   leapR result object/table, or NULL on error
  tryCatch(leapR::leapR(...),
           error = function(e) { message("[leapR] ", conditionMessage(e)); NULL })
}

# Load a leapR built-in geneset by name
.load_leapr_geneset_by_name <- function(name) {
  #   Load a built-in leapR geneset dataset by name (e.g., "kinasesubstrates", "krbpaths").
  #   Validates the name and errors if the dataset cannot be loaded.
  # Inputs:
  #   name: character scalar geneset name
  # Output:
  #   geneset object loaded from the leapR package
  valid <- c("kinasesubstrates", "ncipid", "krbpaths", "longlist", "shortlist")
  if (!(name %in% valid)) {
    stop("Unknown geneset name: '", name, "'. Valid: ", paste(valid, collapse = ", "))
  }
  suppressWarnings(utils::data(list = name, package = "leapR", envir = environment()))
  if (!exists(name, inherits = FALSE)) {
    stop("leapR dataset '", name, "' not found in the installed {leapR}.")
  }
  get(name, inherits = FALSE)
}

# Decide default geneset from omic label when no override is provided
.default_geneset_for_omic <- function(omic_label) {
  #   Choose a default geneset database based on the omics label:
  #   - phospho-like labels -> kinasesubstrates
  #   - otherwise -> krbpaths
  # Inputs:
  #   omic_label: character scalar describing the modality (e.g., "phospho", "rna", "global")
  # Output:
  #   geneset object to use with leapR
  ol <- tolower(omic_label)
  if (ol %in% c("phospho","phosphoproteomics","phosphoprotein","phosphoproteome")) {
    .load_leapr_geneset_by_name("kinasesubstrates")
  } else {
    .load_leapr_geneset_by_name("krbpaths")
  }
}

# -----------------------------
# Main
# -----------------------------
run_leapr_directional_one_cached <- function(
  #   For each drug, correlate uM_viability with each omics feature across samples, split features
  #   into TOP (positive; more resistant) and BOTTOM (negative; more sensitive), then run leapR
  #   enrichment separately on each direction. Supports phospho-specific site->gene handling,
  #   optional site-normalization for kinasesubstrates, CSV writing, and caching to .RData.
  # Inputs:
  #   drugs: long drug-response data.frame (must include improve_drug_id, improve_sample_id, dose_response_metric, dose_response_value)
  #   df_long: long omics data.frame (sample/feature/value columns)
  #   sample_col: column name in df_long for sample IDs
  #   feature_col: column name in df_long for feature IDs (gene/site)
  #   value_col: column name in df_long for numeric measurement
  #   omic_label: modality label used in assay naming and output filenames (e.g., "rna", "global", "phospho")
  #   cache_path: file path to save/load cached results (.RData); skipped if always_rerun=TRUE
  #   write_csvs: TRUE/FALSE; write per-drug TOP/BOTTOM leapR tables to CSV
  #   always_rerun: TRUE/FALSE; ignore cache and recompute
  #   min_features: minimum features required to run leapR for TOP/BOTTOM
  #   test_one: TRUE/FALSE; only run the first drug (debug)
  #   geneset_name: optional built-in leapR geneset name
  #   geneset_object: optional geneset object to use directly (overrides geneset_name/default)
  # Output:
  #   named list by drug: res_list[[drug]]$top and res_list[[drug]]$bottom (leapR results or NULL)
    drugs,                 # Character vector of drugs to test (IDs/names used by your fits/model)
    df_long,               # Long-format omics table (one row per sample x feature)
    sample_col,            # Column name in df_long containing sample IDs
    feature_col,           # Column name in df_long containing feature IDs (e.g., gene/site)
    value_col,             # Column name in df_long containing numeric values to analyze
    omic_label,            # Short label for this modality (used in logs/output names), e.g. "RNA"
    cache_path,            # File path to cache (read/write) computed results
    write_csvs     = FALSE,# If TRUE, write result/intermediate CSVs to disk
    always_rerun   = FALSE,# If TRUE, ignore cache and recompute even if cache exists
    min_features   = 5,    # Minimum # of features required to run; otherwise skip/return early
    test_one       = FALSE,# If TRUE, run a single test case (e.g., first drug) for debugging
    geneset_name   = NULL, # Optional geneset label (used for naming outputs/plot titles)
    geneset_object = NULL  # Optional geneset definition (e.g., character vector) to filter features
  ) {

  # cache check! If the cached value exists, stop there.
  if (!always_rerun && is.character(cache_path) && nzchar(cache_path) && file.exists(cache_path)) {
    load(cache_path) # loads res_list
    if (exists("res_list")) return(res_list)
  }

  # pivot long to matrix
  feat_mat <- long_to_matrix(df_long, sample_col, feature_col, value_col)
  if (is.null(feat_mat) || !nrow(feat_mat) || !ncol(feat_mat)) {
    warning("[run_leapr_directional_one_cached] Empty feature matrix after pivot; returning empty list.")
    return(list())
  }

  # pick geneset
  if (!is.null(geneset_object)) {
    geneset_db <- geneset_object
  } else if (!is.null(geneset_name)) {
    geneset_db <- .load_leapr_geneset_by_name(geneset_name)
  } else {
    geneset_db <- .default_geneset_for_omic(omic_label)
  }

  # optional phospho site gene mapping
  map_site2gene <- NULL
  is_phospho <- tolower(omic_label) %in% c("phospho","phosphoproteomics","phosphoprotein","phosphoproteome")
  if (is_phospho) {
    map_site2gene <- .build_phospho_gene_map_from_long(df_long, feature_col)
    if (is.null(map_site2gene) || !length(map_site2gene)) {
      phos_features <- colnames(feat_mat)
      map_site2gene <- setNames(
        vapply(phos_features, .extract_gene_from_site, FUN.VALUE = character(1)),
        phos_features
      )
    }
    # Print statements
    feats <- colnames(feat_mat)
    mapped <- map_site2gene[feats]
    n_mapped <- sum(!is.na(mapped) & mapped != "")
    message(sprintf("[phospho mapping] %d/%d sites mapped to gene symbols (%.1f%%)",
                    n_mapped, length(feats), 100 * n_mapped / max(1, length(feats))))
    if (n_mapped < length(feats)) {
      unm <- feats[is.na(mapped) | mapped == ""]
      if (length(unm)) {
        show_n <- min(5L, length(unm))
        message("[phospho mapping] Unmapped examples: ",
                paste(utils::head(unm, show_n), collapse = ", "),
                if (length(unm) > show_n) paste0(" ... +", length(unm) - show_n, " more") else "")
      }
    }
  }

  # For phospho, detect when we are using site-level kinase substrates
  uses_kinase_sites <- is_phospho && {
    if (!is.null(geneset_name)) {
      identical(geneset_name, "kinasesubstrates")
    } else {
      identical(geneset_db, .load_leapr_geneset_by_name("kinasesubstrates"))
    }
  }

  # If using kinasesubstrates, normalize column names so they match site IDs
  if (uses_kinase_sites) {
    old_sites  <- colnames(feat_mat)
    norm_sites <- .normalize_kinase_site_id(old_sites)
    if (!identical(old_sites, norm_sites)) {
      message("[kinasesubstrates] Normalizing phosphosite IDs (e.g. 'AAAS-S495s' -> 'AAAS-S495')")
      colnames(feat_mat) <- make.unique(norm_sites)
    }
  }

  res_list <- list()
  out_csv_dir <- file.path("leapR_top_paths", "dir_split")
  if (write_csvs && !dir.exists(out_csv_dir)) dir.create(out_csv_dir, recursive = TRUE)

  all_drugs <- unique(drugs$improve_drug_id)
  if (test_one && length(all_drugs) > 0) {
    message("[run_leapr_directional_one_cached] test_one=TRUE then running only the first drug: ", all_drugs[[1]])
    all_drugs <- all_drugs[[1]]
  } else {
    all_drugs <- sort(all_drugs)
  }

  total <- length(all_drugs)
  for (i in seq_along(all_drugs)) {
    drug <- all_drugs[[i]]
    message(sprintf("[%-3d/%-3d] %s", i, total, drug))

    # mean response per sample for uM_viability
    dv <- drugs %>%
      dplyr::filter(.data$improve_drug_id == !!drug,
                    .data$dose_response_metric == "uM_viability") %>%
      dplyr::group_by(.data$improve_sample_id) %>%
      dplyr::summarise(resp = mean(.data$dose_response_value, na.rm = TRUE),
                       .groups = "drop")

    if (!nrow(dv)) {
      message(" No response rows for metric 'uM_viability'; skipping.")
      next
    }
    dv_vec <- stats::setNames(dv$resp, dv$improve_sample_id)

    # correlations at site-level (or normalized site-level for kinasesubstrates)
    cors <- .col_spearman(dv_vec, feat_mat)
    pos  <- cors[!is.na(cors) & cors > 0]   # resistant (TOP)
    neg  <- cors[!is.na(cors) & cors < 0]   # sensitive (BOTTOM; flip)
    message(sprintf(" Features (site-level): pos=%d, neg=%d (min_features=%d)",
                    length(pos), length(neg), min_features))

    if (uses_kinase_sites && i == 1) {
      ks_sites <- unique(unlist(geneset_db[["matrix"]]))
      ov <- intersect(names(cors), ks_sites)
      message("[kinasesubstrates] Overlapping sites with geneset (first drug): ", length(ov))
      if (length(ov)) {
        message(" Example overlaps: ", paste(utils::head(ov, 5), collapse = ", "))
      }
    }

    # For phospho + GENE-LEVEL sets, collapse site to gene before SE
    if (is_phospho && !uses_kinase_sites) {
      pos <- .collapse_sites_to_genes(pos, map_site2gene, agg = "mean")
      neg <- .collapse_sites_to_genes(neg, map_site2gene, agg = "mean")
      message(sprintf(" Gene-level: pos=%d, neg=%d", length(pos), length(neg)))
    }

    res_list[[drug]] <- list(top = NULL, bottom = NULL)

    # TOP (resistant)
    if (length(pos) >= min_features) {
      feats_top <- names(pos)
      se_top <- .build_se_from_corvec(
        cor_named_vec = pos,
        features_all  = feats_top,
        col_label     = paste0(drug, "_TOP"),
        map_to_gene   = if (is_phospho) NULL else NULL,   # not needed when features are genes/sites
        assay_label   = omic_label
      )
      top_res <- .safe_leapr(
        geneset            = geneset_db,
        enrichment_method  = "enrichment_in_order",
        eset               = se_top,
        assay_name         = omic_label,
        primary_columns    = paste0(drug, "_TOP"),
        id_column          = NULL
      )
      res_list[[drug]]$top <- top_res
      message(" TOP  (resistant): ", if (is.null(top_res)) "no result" else "OK")
      if (write_csvs && !is.null(top_res)) {
        utils::write.csv(as.data.frame(top_res),
                         file = file.path(out_csv_dir, paste0(drug, "_", omic_label, "_TOP.csv")),
                         row.names = FALSE)
      }
    } else {
      message(" TOP  (resistant): skipped (too few positive features)")
    }

    # BOTTOM (sensitive)
    if (length(neg) >= min_features) {
      # flip sign so strong negatives rank to top
      neg_flip <- -neg
      feats_bot <- names(neg_flip)
      se_bot <- .build_se_from_corvec(
        cor_named_vec = neg_flip,
        features_all  = feats_bot,
        col_label     = paste0(drug, "_BOTTOM"),
        map_to_gene   = if (is_phospho) NULL else NULL,
        assay_label   = omic_label
      )
      bot_res <- .safe_leapr(
        geneset            = geneset_db,
        enrichment_method  = "enrichment_in_order",
        eset               = se_bot,
        assay_name         = omic_label,
        primary_columns    = paste0(drug, "_BOTTOM"),
        id_column          = NULL
      )
      res_list[[drug]]$bottom <- bot_res
      message(" BOTTOM(sensitive): ", if (is.null(bot_res)) "no result" else "OK")
      if (write_csvs && !is.null(bot_res)) {
        utils::write.csv(as.data.frame(bot_res),
                         file = file.path(out_csv_dir, paste0(drug, "_", omic_label, "_BOTTOM.csv")),
                         row.names = FALSE)
      }
    } else {
      message(" BOTTOM(sensitive): skipped (too few negative features)")
    }

    if (isTRUE(test_one)) break
  }

  if (is.character(cache_path) && nzchar(cache_path)) {
    save(res_list, file = cache_path)
  }
  res_list
}

# -----------------------------
# Plot and save using leapR builtin plotter
# -----------------------------
save_leapr_plots <- function(
  #   Save leapR pathway barplots for TOP (resistant) and BOTTOM (sensitive) results for each drug.
  #   Supports plotting all drugs or a requested subset (case-insensitive matching).
  # Inputs:
  #   res_list: named list returned by run_leapr_directional_one_cached()
  #   omic_label: modality label used in plot titles and filenames
  #   top_n: number of top pathways to plot per direction
  #   drugs: NULL to plot all; otherwise character vector of drug IDs/names to plot (case-insensitive)
  #   outdir: output directory for PDF files
  # Output:
  #   invisible(NULL); side-effect is writing PDF plots to outdir
    res_list,
    omic_label,
    top_n = 15,
    drugs = NULL,            # NULL = plot all drugs in res_list; otherwise character vector of drug IDs/names (case-insensitive)
    outdir = "figs"          # output directory for PDFs
) {
  if (!length(res_list)) return(invisible(NULL))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  safelabel <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)

  all_drugs <- names(res_list)
  if (!length(all_drugs)) {
    message("[save_leapr_plots] res_list has no named drug entries.")
    return(invisible(NULL))
  }

  # Case-insensitive lookup table: UPPER(drug) -> original name in res_list
  key_upper <- toupper(all_drugs)
  drug_map  <- stats::setNames(all_drugs, key_upper)

  # Decide which drugs to plot
  if (is.null(drugs)) {
    plot_drugs <- all_drugs
  } else {
    req <- as.character(drugs)
    req_upper <- toupper(req)

    found_upper <- intersect(req_upper, names(drug_map))
    plot_drugs  <- unname(drug_map[found_upper])

    missing_upper <- setdiff(req_upper, names(drug_map))
    if (length(missing_upper)) {
      missing_original <- req[req_upper %in% missing_upper]
      message("[save_leapr_plots] Skipping ", length(missing_original),
              " requested drug(s) not present in res_list (case-insensitive match): ",
              paste(utils::head(missing_original, 10), collapse = ", "),
              if (length(missing_original) > 10) paste0(" ... +", length(missing_original) - 10, " more") else "")
    }
  }

  if (!length(plot_drugs)) {
    message("[save_leapr_plots] No matching drugs to plot.")
    return(invisible(NULL))
  }

  for (drug in plot_drugs) {
    two <- res_list[[drug]]
    if (is.null(two)) next

    # TOP (resistant)
    if (!is.null(two$top)) {
      p_top <- leapR::plot_leapr_bar(
        two$top,
        title = paste0(drug, " — ", omic_label, " (Resistant)"),
        top_n = top_n
      )
      if (!is.null(p_top)) {
        fn <- file.path(outdir, paste0(
          "pathways_", safelabel(drug), "_", omic_label, "_resistant_top", top_n, ".pdf"
        ))
        ggplot2::ggsave(fn, p_top, width = 7, height = 5, device = grDevices::cairo_pdf)
      }
    }

    # BOTTOM (Sensitive)
    if (!is.null(two$bottom)) {
      p_bot <- leapR::plot_leapr_bar(
        two$bottom,
        title = paste0(drug, " — ", omic_label, " (Sensitive)"),
        top_n = top_n
      )
      if (!is.null(p_bot)) {
        fn <- file.path(outdir, paste0(
          "pathways_", safelabel(drug), "_", omic_label, "_sensitive_top", top_n, ".pdf"
        ))
        ggplot2::ggsave(fn, p_bot, width = 7, height = 5, device = grDevices::cairo_pdf)
      }
    }
  }

  invisible(NULL)
}
