# =============================================================================
# normalize_omics_pipeline.R
#
# Main entry point:
#   run_modality(modality, batches, meta, syn, ...)
#
# Key inputs (with examples):
#
# Example batches value:
# batches <- list(
#   list(syn_id = "syn69963552", cohort = 1, value_start_col = 5, fname_aliquot_index = 8),
#   list(syn_id = "syn69947351", cohort = 2, value_start_col = 5, fname_aliquot_index = 9)
# )
#
# Arguments:
#  - modality: Which data type to run: "phospho", "global", or "rna"
#  - batches: List of batch configs; each element should include at least:
#       * syn_id: Synapse file ID for the wide feature×sample table
#       * cohort: Cohort/batch label used for joining meta and for ComBat batching
#     Optional per-batch fields:
#       * value_start_col: Column index where sample measurement columns begin (auto-detected if NULL)
#       * fname_aliquot_index: Token index (split on "_") used to parse aliquot number from sample filenames
#  - meta: Sample metadata table used to join batch sample IDs to Patient/Tumor/Specimen
#          (cnF_helper_code.R creates this.)
#  - syn: Synapse client object used for Synapse reading/upload (synapser)
#  - drop_name_substrings: Regex pattern(s); sample columns matching any pattern will be removed (QC/blank runs)
#  - out_dir: Output directory for generated CSV/PDF files
#  - out_prefix: Default base name for output files; if NULL (recommended), derived from modality (sanitized)
#  - upload_parent_id: Synapse project/folder ID to upload outputs into
#                      (ignored if NULL or write_outputs=FALSE)
#  - pcols: Color vector for PCA plotting (names should match Patient IDs)
#           (cnF_helper_code.R creates this.)
#  - write_outputs: Master toggle to write CSV/PDF outputs and perform uploads
#  - save_basename: Override base name used in output filenames (supersedes out_prefix)
#  - do_batch_correct: If FALSE, skip ComBat and use combined matrix instead
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(SummarizedExperiment)
  library(ggplot2)
  library(readr)
  library(rlang)
})

# Small helpers
modified_zscore <- function(x, na.rm = TRUE) {
  #   Robust z-score for a numeric vector using median and MAD (less sensitive to
  #   outliers than mean/SD). If MAD is 0 or NA, returns all zeros.
  # Inputs:
  #   x: numeric vector
  #   na.rm: TRUE/FALSE; whether to ignore NA when computing median/MAD
  # Output:
  #   numeric vector (same length as x)
  m  <- suppressWarnings(stats::median(x, na.rm = na.rm))
  md <- suppressWarnings(stats::mad(x, constant = 1, na.rm = na.rm))
  if (is.na(md) || md == 0) return(rep(0, length(x)))
  0.6745 * (x - m) / md
}

filter_by_missingness <- function(mat) {
  #   Filter features (rows) by missingness: keep rows where <= 50% of values are NA.
  # Inputs:
  #   mat: numeric matrix (features x samples)
  # Output:
  #   numeric matrix with a subset of rows retained
  keep <- apply(mat, 1, function(r) mean(is.na(r)) <= 0.5)
  mat[keep, , drop = FALSE]
}

union_rows_fill_NA <- function(mats) {
  #   Align a list of matrices to the union of all rownames (features), filling
  #   missing feature rows in each matrix with NA.
  # Inputs:
  #   mats: list of numeric matrices with rownames
  # Output:
  #   list of numeric matrices, each reindexed to the same union row set
  all_feats <- Reduce(union, lapply(mats, rownames))
  lapply(mats, function(m) {
    mm <- matrix(NA_real_, nrow = length(all_feats), ncol = ncol(m),
                 dimnames = list(all_feats, colnames(m)))
    mm[rownames(m), colnames(m)] <- m
    mm
  })
}

collapse_duplicate_features <- function(mat) {
  #   Collapse duplicate feature IDs (duplicate rownames) by summing values across
  #   duplicates for each sample (NA treated as 0 for summation).
  # Inputs:
  #   mat: numeric matrix (features x samples) with rownames as feature IDs
  # Output:
  #   numeric matrix with unique rownames (duplicates collapsed)
  if (!any(duplicated(rownames(mat)))) return(mat)
  grp <- split(seq_len(nrow(mat)), rownames(mat))
  collapsed <- do.call(rbind, lapply(grp, function(ix) colSums(mat[ix, , drop = FALSE], na.rm = TRUE)))
  rownames(collapsed) <- names(grp)
  collapsed
}

make_dropper <- function(substrings) {
  #   Build a function that flags sample column names to drop based on one or more
  #   regex patterns. Useed for removing protocol optimization runs.
  # Inputs:
  #   substrings: NULL or character vector of regex patterns
  # Output:
  #   function(x): logical vector; TRUE means "drop this name"
  if (is.null(substrings) || length(substrings) == 0) return(function(x) rep(FALSE, length(x)))
  pattern <- paste0(substrings, collapse = "|")
  function(x) grepl(pattern, x, fixed = FALSE)
}

# Functions to clean up irregular names
#   Extract the filename portion from a path or portion of path (drops directories).
# Inputs:
#   x: character vector of file paths
# Output:
#   character vector of basenames
basename_only   <- function(x) sub("^.*[\\\\/]", "", x)
basename_no_ext <- function(x) sub("\\.[^.]+$", "", basename_only(x))

normalize_specimen_like <- function(x) {
  #   Normalize specimen strings to a consistent form for joining (lowercase,
  #   remove whitespace, unify separators, normalize organoid/tissue/skin labels).
  # Inputs:
  #   x: character vector
  # Output:
  #   character vector of normalized specimen-like strings
  y <- tolower(x)
  y <- gsub("\\s+", "", y)
  y <- gsub("\\.", "-", y)
  y <- gsub("_", "-", y)
  y <- gsub("organoids?$", "organoid", y)
  y <- gsub("-organoids?-", "-organoid-", y)
  y <- gsub("skin$", "skin", y)
  y <- gsub("tissues?$", "tissue", y)
  y <- gsub("--+", "-", y)
  y <- gsub("^-|-$", "", y)
  y
}

parse_rna_header_triplet <- function(fnames) {
  #   Parse RNA sample headers expected to look like "sample.T1.condition" (dot-delimited).
  #   Extracts sample_id, optional tumor (T#), and condition; also builds a normalized
  #   specimen key to help join against metadata.
  # Inputs:
  #   fnames: character vector of RNA sample column names
  # Output:
  #   data.frame with columns: fname, sample_id, tumor, condition_raw, condition_norm, specimen_norm
  toks_list <- strsplit(fnames, "\\.")
  out <- lapply(seq_along(toks_list), function(i) {
    toks <- toks_list[[i]]
    sample_id <- if (length(toks) >= 1) toks[[1]] else NA_character_

    tumor <- NA_character_
    condition_raw <- NA_character_

    if (length(toks) >= 2) {
      if (grepl("^T\\d+$", toks[[2]], ignore.case = TRUE)) {
        tumor <- toupper(toks[[2]])
        cond_tokens <- toks[-c(1,2)]
        condition_raw <- if (length(cond_tokens)) paste(cond_tokens, collapse = ".") else NA_character_
      } else {
        cond_tokens <- toks[-1]
        condition_raw <- if (length(cond_tokens)) paste(cond_tokens, collapse = ".") else NA_character_
      }
    }

    condition_norm <- condition_raw
    if (!is.na(condition_norm)) {
      low <- tolower(condition_norm)
      if (grepl("^organoids?$", low)) condition_norm <- "organoid"
      else if (grepl("^tissues?$", low)) condition_norm <- "tissue"
      else if (grepl("^skin$", low)) condition_norm <- "skin"
    }

    data.frame(
      fname         = fnames[i],
      sample_id     = sample_id,
      tumor         = ifelse(is.na(tumor), NA_character_, toupper(tumor)),
      condition_raw = condition_raw,
      condition_norm= condition_norm,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, out)
  df$specimen_norm <- with(df, {
    sid <- sample_id
    tmr <- ifelse(is.na(tumor) | tumor == "", "", paste0("_", tumor))
    cnd <- ifelse(is.na(condition_norm) | condition_norm == "", "", paste0("_", condition_norm))
    paste0(sid, tmr, cnd)
  })
  df
}

# Functions to get data from Synapse
read_wide_from_synapse <- function(syn, syn_id) {
  #   Download a Synapse file and read it as a wide tab-delimited table (features x samples).
  # Inputs:
  #   syn: Synapse client object (synapser)
  #   syn_id: Synapse file ID (e.g., "syn69963552")
  # Output:
  #   data.frame containing the wide table (annotation columns + sample columns)
  message(" Reading Synapse file: ", syn_id)
  df <- read.table(
    syn$get(syn_id)$path,
    sep = "\t", header = TRUE, quote = '"',
    fill = TRUE, check.names = FALSE
  )
  message(sprintf("    - Read %d rows × %d cols; first cols: %s",
                  nrow(df), ncol(df), paste(head(colnames(df), 8), collapse = ", ")))
  df
}

detect_value_start_col <- function(wide_df, fallback = 5) {
  #   Auto-detect where sample measurement columns start in a wide table by looking for
  #   headers that resemble file paths or RAW/mzML names. Falls back if not found.
  # Inputs:
  #   wide_df: data.frame (wide)
  #   fallback: integer column index to use if auto-detection fails
  # Output:
  #   integer column index for the first sample column
  nms <- colnames(wide_df)
  is_pathy <- grepl("\\.(raw|mzml)$", nms, ignore.case = TRUE) |
    grepl("[/\\\\]", nms) |
    grepl("^[A-Za-z]:\\\\", nms)
  if (any(is_pathy)) {
    i <- which(is_pathy)[1]
    message(" Auto-detected first sample column at index ", i, " to '", nms[i], "'")
    return(i)
  }
  message(" Did not detect path/RAW headers; using fallback value_start_col=", fallback)
  fallback
}

parse_fnames <- function(fnames, aliquot_field_index, cohort) {
  #   Parse sample column names into (fname, aliquot, cohort). Attempts to extract aliquot
  #   from a specific underscore token index; otherwise tries the last numeric token.
  # Inputs:
  #   fnames: character vector of sample column names
  #   aliquot_field_index: integer token index (split on "_") or NULL
  #   cohort: cohort label to attach to all parsed samples
  # Output:
  #   data.frame with columns: fname, aliquot (numeric or NA), cohort
  message("Parsing filenames to (fname, aliquot, cohort)")
  rows <- lapply(fnames, function(fname) {
    toks <- strsplit(fname, "_", fixed = TRUE)[[1]]
    aliq <- NA_real_

    if (!is.null(aliquot_field_index) &&
        aliquot_field_index >= 1 &&
        aliquot_field_index <= length(toks)) {
      aliq_try <- suppressWarnings(as.double(toks[[aliquot_field_index]]))
      if (!is.na(aliq_try)) aliq <- aliq_try
    }

    if (is.na(aliq)) {
      num_tokens <- suppressWarnings(as.double(toks))
      if (any(!is.na(num_tokens))) aliq <- tail(num_tokens[!is.na(num_tokens)], 1)
    }

    data.frame(
      fname   = fname,
      aliquot = aliq,
      cohort  = cohort,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  message(sprintf("    - Parsed %d samples (aliquot NA: %d)",
                  nrow(out), sum(is.na(out$aliquot))))
  out
}

#####
#Feature ID builders
#####
build_phospho_ids <- function(df) {
  #   Build unique phosphosite feature IDs from phospho annotation columns.
  # Inputs:
  #   df: data.frame containing at least Gene.Names, Residue, Site
  # Output:
  #   character vector of feature IDs (one per row)
  lsite <- tolower(df$Residue)
  paste0(df$`Gene.Names`, "-", df$Residue, df$Site, lsite)
}

#   Build global proteomics feature IDs (gene symbols).
# Inputs:
#   df: data.frame containing a Genes column
# Output:
#   character vector of feature IDs (one per row)
build_global_ids <- function(df) as.character(df$Genes)

build_rna_ids <- function(df) {
  #   Build RNA feature IDs by selecting a gene identifier column (tries common names
  #   like gene_id, gene_name, Symbol, Ensembl). Errors if none found.
  # Inputs:
  #   df: data.frame containing a recognized gene ID column
  # Output:
  #   character vector of feature IDs (one per row)
  cand <- c("gene_id","Gene","gene","gene_name","Symbol","symbol","ENSEMBL","Ensembl","ensembl_gene_id")
  hit  <- cand[cand %in% names(df)]
  if (length(hit) == 0) stop("RNA feature-id column not found.")
  if ("gene_id" %in% hit) return(as.character(df[["gene_id"]]))
  if ("gene_name" %in% hit) return(as.character(df[["gene_name"]]))
  as.character(df[[hit[[1]]]])
}
pick_builder <- function(modality) {
  #   Choose the correct feature-ID builder function based on modality.
  # Inputs:
  #   modality: "phospho", "global", or "rna" (case-insensitive)
  # Output:
  #   function(df) -> character vector of feature IDs
  m <- tolower(modality)
  if (m == "phospho") return(build_phospho_ids)
  if (m == "global")  return(build_global_ids)
  if (m == "rna")     return(build_rna_ids)
  stop("Unknown modality: ", modality)
}

# Functions to Normalize Data. Uses SummarizedExperiment

coldata_tbl <- function(se) {
  #   Convert SummarizedExperiment colData into a clean data.frame with consistent
  #   filename fields (fname, basename, stem). Avoids name collisions.
  # Inputs:
  #   se: SummarizedExperiment
  # Output:
  #   data.frame of sample metadata; includes fname, fname_base, fname_stem
  cd <- as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors = FALSE)
  if ("fname" %in% names(cd)) names(cd)[names(cd) == "fname"] <- ".coldata_fname"
  names(cd) <- make.unique(names(cd), sep = "_")
  cd$fname       <- rownames(cd)
  cd$fname_base  <- basename_only(cd$fname)
  cd$fname_stem  <- basename_no_ext(cd$fname)
  cd <- cd[, c("fname","fname_base","fname_stem", setdiff(names(cd), c("fname","fname_base","fname_stem"))), drop = FALSE]
  cd
}

looks_like_sample_header <- function(x) {
  #   Heuristic test for whether a column name looks like a raw file/sample path
  #   (e.g., contains slashes or ends in .raw/.mzml).
  # Inputs:
  #   x: character vector of column names
  # Output:
  #   logical vector; TRUE indicates "looks like a sample header"
  # Example file:
  # "I:\UserData\LeDay\Piehowski_orgonoids_Feb25\RawData\1338241_cNF_organoid_DIA_P_01_29Jan25_Ned_BEHCoA-25-01-02.raw"
  grepl("\\.(raw|mzml)$", x, ignore.case = TRUE) |
    grepl("[/\\\\]", x) |
    grepl("^[A-Za-z]:\\\\", x)
}

make_se <- function(wide_df, value_start_col, feature_ids, fnames_df, meta, drop_name, modality) {
  #   Build a SummarizedExperiment from a wide feature x sample table:
  #   - selects sample columns
  #   - converts values to numeric
  #   - drops unwanted sample columns by name pattern
  #   - attaches sample metadata by joining on (aliquot, cohort)
  #   - applies extra RNA-specific parsing and metadata reconciliation
  # Inputs:
  #   wide_df: wide data.frame (features + sample columns)
  #   value_start_col: integer index of first sample column
  #   feature_ids: character vector of feature IDs (length = nrow(wide_df))
  #   fnames_df: data.frame mapping fname->aliquot/cohort (from parse_fnames)
  #   meta: metadata table used to map aliquot/cohort to Patient/Tumor/Specimen
  #   drop_name: function(x)->logical; TRUE means drop that sample column
  #   modality: "phospho", "global", or "rna"
  # Output:
  #   SummarizedExperiment with assay "values" and populated colData/rowData
  all_candidate <- colnames(wide_df)[value_start_col:ncol(wide_df)]
  has_pathy <- any(looks_like_sample_header(all_candidate))
  if (has_pathy) {
    sample_cols <- all_candidate[looks_like_sample_header(all_candidate)]
    sample_cols <- setdiff(sample_cols, c("Site", "Sequence"))
  } else {
    sample_cols <- setdiff(all_candidate, c("Site", "Sequence"))
  }

  # message(" Candidate sample columns (first 6):")
  # print(utils::head(sample_cols, 6))

  message(" Casting measurement block to numeric")
  raw_block <- wide_df[, sample_cols, drop = FALSE]
  clean_block <- as.data.frame(
    lapply(raw_block, function(col) {
      if (is.factor(col)) col <- as.character(col)
      col[col %in% c("", "NA", "NaN", "na", "n/a", "NULL")] <- NA
      suppressWarnings(as.numeric(col))
    }),
    check.names = FALSE
  )

  keep <- !drop_name(colnames(clean_block))
  if (any(!keep)) {
    message(" Dropping unwanted sample columns by pattern: ", sum(!keep))
    clean_block <- clean_block[, keep, drop = FALSE]
    sample_cols <- sample_cols[keep]
  }

  mat <- as.matrix(clean_block)
  rownames(mat) <- feature_ids
  colnames(mat) <- sample_cols

  fnames_df <- fnames_df %>% dplyr::semi_join(data.frame(fname = sample_cols), by = "fname")

  parsed_df <- if (tolower(modality) == "rna") parse_rna_header_triplet(fnames_df$fname) else
    data.frame(fname = fnames_df$fname, stringsAsFactors = FALSE)

  message(" Joining sample meta (by aliquot & cohort)")
  meta_join <- meta
  if ("fname" %in% names(meta_join)) {
    message(" ! Debug: Dropping 'fname' column from 'meta' to prevent duplication")
    meta_join <- dplyr::select(meta_join, -fname)
  }

  cdata <- fnames_df %>%
    dplyr::left_join(meta_join, by = c("aliquot","cohort")) %>%
    dplyr::left_join(parsed_df, by = "fname") %>%
    dplyr::mutate(cohort = as.factor(cohort),
                  cohort_key = as.character(cohort))

  if (!"aliquot" %in% names(cdata)) cdata$aliquot <- NA_real_
  if ("aliquot.x" %in% names(cdata) || "aliquot.y" %in% names(cdata)) {
    cdata$aliquot <- dplyr::coalesce(cdata$aliquot, cdata$aliquot.x, cdata$aliquot.y)
    cdata <- dplyr::select(cdata, -dplyr::any_of(c("aliquot.x","aliquot.y")))
  }

  if (tolower(modality) == "rna") {
    meta_norm <- meta %>%
      dplyr::mutate(
        Specimen_norm = normalize_specimen_like(Specimen),
        cohort_key    = as.character(cohort)
      ) %>%
      dplyr::select(Specimen, Specimen_norm, Patient, Tumor, cohort_key)

    cdata2 <- cdata %>%
      dplyr::mutate(fname_norm = normalize_specimen_like(ifelse(
        is.na(specimen_norm) | specimen_norm == "", fname, specimen_norm
      ))) %>%
      dplyr::left_join(meta_norm, by = c("cohort_key", "fname_norm" = "Specimen_norm"))

    for (nm in c("Specimen","Patient","Tumor")) {
      if (!nm %in% names(cdata)) cdata[[nm]] <- NA
      cdata[[nm]] <- dplyr::coalesce(cdata[[nm]], cdata2[[nm]])
    }

    if (!"Specimen" %in% names(cdata)) cdata$Specimen <- NA_character_
    if (!"Patient"  %in% names(cdata)) cdata$Patient  <- NA_character_
    if (!"Tumor"    %in% names(cdata)) cdata$Tumor    <- NA_character_

    cdata$Patient <- dplyr::coalesce(cdata$Patient, cdata$sample_id)
    cdata$Tumor   <- dplyr::coalesce(cdata$Tumor, cdata$tumor)

    make_specimen <- function(pid, tmr, cond) {
      pid_clean  <- pid
      tmr_clean  <- ifelse(is.na(tmr)  | tmr  == "", "", paste0("_", tmr))
      cond_clean <- ifelse(is.na(cond) | cond == "", "", paste0("_", cond))
      paste0(pid_clean, tmr_clean, cond_clean)
    }
    need_spec <- is.na(cdata$Specimen) | cdata$Specimen == ""
    if (any(need_spec)) {
      cdata$Specimen[need_spec] <- make_specimen(
        pid  = cdata$Patient[need_spec],
        tmr  = cdata$Tumor[need_spec],
        cond = cdata$condition_norm[need_spec]
      )
    }
  }

  message(" - Example cdata rows:")
  # print(utils::head(cdata[, intersect(c(
  #   "fname","aliquot","cohort","Specimen","Patient","Tumor",
  #   "sample_id","tumor","condition_raw","condition_norm","specimen_norm"
  # ), names(cdata)), drop = FALSE], 10))

  rn <- cdata$fname
  cdata_nofname <- dplyr::select(cdata, -fname)

  rd <- S4Vectors::DataFrame(feature_id = rownames(mat)); rownames(rd) <- rd$feature_id

  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = S4Vectors::SimpleList(values = mat),
    rowData = rd,
    colData = S4Vectors::DataFrame(cdata_nofname, row.names = rn)
  )

  message(sprintf(" - SE assay dims: %d feats × %d samples", nrow(se), ncol(se)))
  message(sprintf(" - colData names: %s", paste(names(SummarizedExperiment::colData(se)), collapse = ", ")))
  se
}

scale_columns_modified_z <- function(m) {
  #   Apply modified_zscore() to each column of a matrix (per-sample robust scaling).
  # Inputs:
  #   m: numeric matrix (features x samples)
  # Output:
  #   numeric matrix of same dimensions (column-wise robust z-scored)
  out <- m
  for (j in seq_len(ncol(m))) out[, j] <- modified_zscore(m[, j])
  out
}

normalize_by_modality <- function(se, modality) {
  #   Normalize a SummarizedExperiment assay using modality-specific transforms:
  #   phospho: 0->NA, filter missingness, log2(x+0.01), robust zscore
  #   global: log2(x), robust zscore
  #   rna: filter missingness, log2(x+1), robust zscore
  #   Also collapses duplicated feature IDs.
  # Inputs:
  #   se: SummarizedExperiment with assay "values"
  #   modality: "phospho", "global", or "rna"
  # Output:
  #   SummarizedExperiment with normalized assay "values"
  mtype <- tolower(modality)
  mat0  <- as.matrix(SummarizedExperiment::assay(se, "values"))
  c0    <- colnames(mat0)

  message(" Normalizing modality = ", modality)
  if (mtype == "phospho") {
    mat0[mat0 == 0] <- NA_real_
    mat1 <- filter_by_missingness(mat0)
    mlog <- log2(mat1 + 0.01)
    mat2 <- scale_columns_modified_z(mlog)
  } else if (mtype == "global") {
    mlog <- log2(mat0)
    mat2 <- scale_columns_modified_z(mlog)
  } else if (mtype == "rna") {
    mat1 <- filter_by_missingness(mat0)
    mlog <- log2(mat1 + 1)
    mat2 <- scale_columns_modified_z(mlog)
  } else {
    abort(paste0("Unknown modality: ", modality))
  }

  mat3 <- collapse_duplicate_features(mat2)
  stopifnot(identical(colnames(mat3), c0))

  rd <- S4Vectors::DataFrame(feature_id = rownames(mat3)); rownames(rd) <- rd$feature_id
  se_out <- SummarizedExperiment::SummarizedExperiment(
    assays  = S4Vectors::SimpleList(values = mat3),
    rowData = rd,
    colData = SummarizedExperiment::colData(se)[colnames(mat3), , drop = FALSE]
  )
  message(sprintf(" - Normalized assay dims: %d feats × %d samples", nrow(se_out), ncol(se_out)))
  se_out
}

# Function to Combine Data (batches)

combine_batches_intersection <- function(se_list) {
  #   Combine multiple normalized batches by intersecting shared features (rownames),
  #   then concatenating samples (cbind). Also stacks colData.
  # Inputs:
  #   se_list: list of SummarizedExperiment objects (normalized)
  # Output:
  #   SummarizedExperiment containing combined matrix and combined colData
  message("Combining batches (intersection of features, then cbind samples)")
  mats  <- lapply(se_list, function(se) as.matrix(SummarizedExperiment::assay(se, "values")))
  feats <- Reduce(intersect, lapply(mats, rownames))
  feats <- feats[!is.na(feats) & feats != ""]
  message("    - Intersection feature count: ", length(feats))
  if (length(feats) == 0) stop("No common features across batches after cleaning feature IDs.")
  matsI <- lapply(seq_along(mats), function(i) {
    m <- mats[[i]][feats, , drop = FALSE]
    message(sprintf("      Batch %d: %d feats × %d samples after intersect", i, nrow(m), ncol(m)))
    m
  })
  matCB <- do.call(cbind, matsI)

  cd <- do.call(S4Vectors::rbind, lapply(se_list, SummarizedExperiment::colData))
  rd <- S4Vectors::DataFrame(feature_id = rownames(matCB)); rownames(rd) <- rd$feature_id

  se <- SummarizedExperiment(
    assays  = S4Vectors::SimpleList(values = matCB),
    rowData = rd,
    colData = cd
  )
  message(sprintf(" - Combined assay dims: %d feats × %d samples", nrow(se), ncol(se)))
  message(" - Head(sample names) in combined assay:")
  # print(utils::head(colnames(SummarizedExperiment::assay(se, "values")), 6))
  message(" - Head(rownames) in combined colData (should match):")
  # print(utils::head(rownames(SummarizedExperiment::colData(se)), 6))
  invisible(se)
}

# Combat Function. (Lots of messages to help debug)

combat_by_cohort <- function(se) {
  #   Batch-correct the combined matrix using ComBat (sva) with colData$cohort as the
  #   batch variable. Replaces non-finite values with 0 before correction.
  # Inputs:
  #   se: SummarizedExperiment with assay "values" and colData column "cohort"
  # Output:
  #   SummarizedExperiment with batch-corrected assay "values"
  message("Running ComBat by cohort (batch-only; mean.only = FALSE)")
  suppressPackageStartupMessages(library(sva))

  mat <- as.matrix(SummarizedExperiment::assay(se, "values"))
  message(sprintf("Matrix dims before ComBat: %d features × %d samples", nrow(mat), ncol(mat)))

  n_bad <- sum(!is.finite(mat))
  if (n_bad > 0) message(" Replacing ", n_bad, " non-finite values with 0.")
  mat[!is.finite(mat)] <- 0

  cd <- as.data.frame(SummarizedExperiment::colData(se))
  cd <- cd[colnames(mat), , drop = FALSE]
  if (!"cohort" %in% names(cd)) stop("colData must contain 'cohort' for ComBat batching.")

  batch <- droplevels(as.factor(cd$cohort))
  # message(" Batch table (pre-drop):"); print(table(batch, useNA = "ifany"))

  keep <- !is.na(batch)
  if (any(!keep)) {
    message(" Dropping ", sum(!keep), " samples with NA cohort before ComBat.")
    mat   <- mat[, keep, drop = FALSE]
    batch <- droplevels(batch[keep])
    cd    <- cd[keep, , drop = FALSE]
  }

  # message(" Final check — ncol(mat)=", ncol(mat), "; length(batch)=", length(batch))
  # message(" Batch table (final):"); print(table(batch, useNA = "ifany"))

  pre_by_cohort <- tapply(colMeans(mat), batch, sd)
  # message(" Pre-ComBat: SD of column means by cohort:"); print(pre_by_cohort)

  cb <- sva::ComBat(dat = mat, batch = batch, mean.only = FALSE, par.prior = TRUE)

  post_by_cohort <- tapply(colMeans(cb), batch, sd)
  # message(" Post-ComBat: SD of column means by cohort:"); print(post_by_cohort)

  SummarizedExperiment::assay(se, "values") <- cb
  # message(" Matrix dims after ComBat:  ", nrow(cb), " × ", ncol(cb))
  invisible(se)
}

# Plot Functions (PCA) + more debug messages

se_to_long <- function(se, modality) {
  #   Convert a SummarizedExperiment matrix into long format (one row per feature-sample
  #   pair) and join sample metadata from colData.
  # Inputs:
  #   se: SummarizedExperiment with assay "values"
  #   modality: "phospho", "global", or "rna" (controls feature column name)
  # Output:
  #   data.frame in long format with correctedAbundance + sample metadata columns
  feature_col <- if (tolower(modality) == "global") "Gene" else "feature_id"

  avals <- as.data.frame(SummarizedExperiment::assay(se, "values"), check.names = FALSE)
  avals[[feature_col]] <- rownames(avals)

  long <- avals |>
    tidyr::pivot_longer(cols = -all_of(feature_col), names_to = "fname", values_to = "correctedAbundance")

  cd <- coldata_tbl(se)

  if (!"aliquot" %in% names(cd)) cd$aliquot <- NA_real_
  if ("aliquot.x" %in% names(cd) || "aliquot.y" %in% names(cd)) {
    cd$aliquot <- dplyr::coalesce(cd$aliquot, cd$aliquot.x, cd$aliquot.y)
  }

  want <- c("fname","aliquot","cohort","Specimen","Patient","Tumor",
            "fname_base","fname_stem",
            "sample_id","tumor","condition_raw","condition_norm","specimen_norm")
  have <- intersect(want, names(cd))
  long1 <- dplyr::left_join(long, cd[, have, drop = FALSE], by = "fname")

  message("    - se_to_long(): non-NA counts to Patient=",
          sum(!is.na(long1$Patient)), "; Tumor=", sum(!is.na(long1$Tumor)),
          "; Specimen=", sum(!is.na(long1$Specimen)), "; cohort=", sum(!is.na(long1$cohort)))
  distinct(long1)
}

pca_df_present_in_all <- function(se) {
  #   Prepare a PCA data.frame using only features that are complete (finite) across
  #   all samples. Joins PCA scores with sample metadata for plotting.
  # Inputs:
  #   se: SummarizedExperiment with assay "values"
  # Output:
  #   data.frame with PC1/PC2 and metadata columns (Patient/Tumor/Specimen/cohort, etc.)
  message("Preparing PCA (features present in ALL samples)")
  mat <- as.matrix(SummarizedExperiment::assay(se, "values"))

  keep_rows <- apply(mat, 1, function(r) all(is.finite(r)))
  n_keep <- sum(keep_rows); n_all <- nrow(mat)
  message(" - Kept ", n_keep, " / ", n_all, " features with complete data for PCA")
  if (n_keep < 2) stop("Too few complete features for PCA after intersection filter.")

  pcs <- prcomp(t(mat[keep_rows, , drop = FALSE]))

  cd <- coldata_tbl(se)
  df <- as.data.frame(pcs$x[, 1:2, drop = FALSE])
  df$fname <- rownames(df)

  df1 <- dplyr::left_join(df, cd, by = "fname")

  message(" - colData columns present: ", paste(setdiff(names(cd), c("fname","fname_base","fname_stem")), collapse = ", "))
  message(" - Non-NA counts in colData: Patient=", sum(!is.na(cd$Patient)),
          "; Tumor=", sum(!is.na(cd$Tumor)), "; Specimen=", sum(!is.na(cd$Specimen)),
          "; cohort=", sum(!is.na(cd$cohort)))
  message(" - After join: n rows = ", nrow(df1))
  message(" - Non-NA counts after join: Patient=", sum(!is.na(df1$Patient)),
          "; Tumor=", sum(!is.na(df1$Tumor)), "; Specimen=", sum(!is.na(df1$Specimen)),
          "; cohort=", sum(!is.na(df1$cohort)))

  if (sum(!is.na(df1$Patient)) == 0) {
    message(" ! Warning: Patient is NA for all samples after join. Will color/shape by cohort.")
    df1$Patient_fallback <- as.character(df1$cohort)
    df1$Tumor_fallback   <- as.character(df1$cohort)
    # message(" - DEBUG: head(df1$fname):"); print(utils::head(df1$fname, 6))
    # message(" - DEBUG: head(cd$fname):");   print(utils::head(cd$fname, 6))
  }

  if ("Specimen" %in% names(df1)) {
    df1$Tumor <- stringr::str_remove(stringr::str_extract(df1$Specimen, "_T\\d+"), "^_")
  }

  df1
}

plot_pca <- function(pc_df, title_text, pcols = NULL) {
  #   Create a PCA scatter plot (PC1 vs PC2), choosing a sensible color/shape mapping
  #   based on available metadata (prefers Patient/Tumor; falls back to cohort/condition).
  # Inputs:
  #   pc_df: data.frame returned by pca_df_present_in_all()
  #   title_text: plot title string
  #   pcols: optional named vector of colors for Patient values
  # Output:
  #   ggplot object (PCA scatter)
  color_col <- if ("Patient" %in% names(pc_df) && any(!is.na(pc_df$Patient))) {
    "Patient"
  } else if ("condition_norm" %in% names(pc_df) && any(!is.na(pc_df$condition_norm))) {
    "condition_norm"
  } else if ("tumor" %in% names(pc_df) && any(!is.na(pc_df$tumor))) {
    "tumor"
  } else {
    "Patient_fallback"
  }

  shape_col <- if ("Tumor" %in% names(pc_df) && any(!is.na(pc_df$Tumor))) {
    "Tumor"
  } else if ("tumor" %in% names(pc_df) && any(!is.na(pc_df$tumor))) {
    "tumor"
  } else {
    "Tumor_fallback"
  }

  g <- ggplot(pc_df, aes(PC1, PC2, col = .data[[color_col]])) +
    geom_point(aes(shape = .data[[shape_col]]), size = 3) +
    labs(title = title_text, color = color_col, shape = shape_col) +
    theme_bw()
  if (!is.null(pcols) && color_col == "Patient") g <- g + scale_color_manual(values = pcols)
  print(g)
  g
}

plot_hist <- function(se, title_text) {
  #   Plot a histogram of all assay values, filled by cohort, to visualize distributions
  #   (e.g., pre- vs post-ComBat).
  # Inputs:
  #   se: SummarizedExperiment with assay "values"
  #   title_text: plot title string
  # Output:
  #   ggplot object (histogram)
  cd <- coldata_tbl(se)
  df <- as.data.frame(SummarizedExperiment::assay(se, "values")) |>
    tidyr::pivot_longer(everything(), names_to = "fname", values_to = "val") |>
    dplyr::left_join(cd[, c("fname","cohort")], by = "fname")
  g <- ggplot(df, aes(x = val, fill = as.factor(cohort))) +
    geom_histogram(bins = 60, alpha = 0.9) +
    labs(title = title_text, x = "Value", fill = "Cohort") +
    theme_bw()
  print(g)
  g
}

# Upload function

perform_uploads <- function(paths, syn, parent_id) {
  #   Upload a set of local output files to a Synapse folder/project using syn$store().
  # Inputs:
  #   paths: character vector of local file paths
  #   syn: Synapse client object (synapser)
  #   parent_id: Synapse folder/project ID to store into
  # Output:
  #   invisible(NULL); side-effect is file uploads to Synapse
  if (is.null(parent_id) || length(paths) == 0) return(invisible(NULL))
  message("All steps succeeded — uploading ", length(paths), " file(s) to Synapse…")
  for (p in paths) {
    fullp <- normalizePath(p, winslash = "/", mustWork = FALSE)
    message(" Uploading: ", basename(p), " (", fullp, ")")
    f <- syn$store(synapser::File(p, parentId = parent_id))
    message(" Uploaded: ", basename(p), " (local: ", fullp, ") to Synapse ID: ", f$properties$id)
  }
  message("Uploads complete.")
}

#####
# Main entry
#####
# This is how we call the function / pipeline

# Example batches value:
# batches <- list(
#   list(syn_id = "syn69963552", cohort = 1, value_start_col = 5, fname_aliquot_index = 8),
#   list(syn_id = "syn69947351", cohort = 2, value_start_col = 5, fname_aliquot_index = 9)
# )


run_modality <- function(
  #   End-to-end normalization pipeline for one modality across one or more batches:
  #   - read wide tables from Synapse
  #   - build feature IDs + construct SummarizedExperiment per batch
  #   - modality-specific normalization per batch
  #   - combine batches on shared feature intersection
  #   - QC plots (PCA + hist) pre and post
  #   - optional ComBat batch correction by cohort
  #   - export long CSVs and optional Synapse uploads
  # Inputs:
  #   modality: "phospho", "global", or "rna"
  #   batches: list of batch configs (syn_id, cohort, optional parsing hints)
  #   meta: sample metadata table for joining
  #   syn: Synapse client object
  #   drop_name_substrings: optional regex patterns to drop sample columns
  #   out_dir: directory for outputs
  #   out_prefix: base name for outputs (defaults from modality)
  #   upload_parent_id: Synapse folder/project ID for uploads (optional)
  #   pcols: optional named color vector for Patient PCA coloring
  #   write_outputs: TRUE/FALSE to write CSV/PDF and upload
  #   save_basename: override base output stem
  #   do_batch_correct: TRUE/FALSE to run ComBat
  # Output:
  #   list containing SE objects, long tables, PCA data, plot objects, and written file paths
    modality,                # Which data type to run: "phospho", "global", or "rna"
    batches,                 # List of batch configs (each element should include at least: syn_id, cohort; optionally: value_start_col, fname_aliquot_index)
    meta,                    # Sample metadata table used to join batch sample IDs to Patient/Tumor/Specimen - cnF_helper_code.R creates this.
    syn,                     # Synapse client object used for synapse reading/upload (synapser)
    drop_name_substrings = NULL, # Regex pattern - sample columns matching any pattern will be removed (QC/blank runs)
    out_dir = ".",           # Output directory for generated CSV/PDF files
    out_prefix = NULL,       # Default base name for output files; if NULL (recommended), derived from modality (sanitized)
    upload_parent_id = NULL, # Synapse project ID to upload outputs into (ignored if NULL or write_outputs=FALSE)
    pcols = NULL,            # Color vector for PCA plotting (names should match Patient IDs) - cnF_helper_code.R creates this.
    write_outputs = TRUE,    # master toggle to write CSV/PDF & upload
    save_basename = NULL,    # override base name used in output files
    do_batch_correct = TRUE  # if FALSE, skip ComBat and use combined matrix
) {
  message("==================================================")
  message("Starting run_modality(): ", modality)
  message("Output directory: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  if (is.null(out_prefix)) out_prefix <- gsub("[^A-Za-z0-9]+", "_", tolower(modality))
  file_stem <- if (!is.null(save_basename) && nzchar(save_basename)) save_basename else out_prefix

  upload_queue <- character(0)
  results <- NULL

  tryCatch({
    drop_name <- make_dropper(drop_name_substrings)
    build_ids <- pick_builder(modality)

    # ---- Per-batch normalization ------------------------------------------------
    se_list <- vector("list", length(batches))
    for (i in seq_along(batches)) {
      message("--------------------------------------------------")
      message("Batch ", i, " of ", length(batches))
      b <- batches[[i]]

      wide <- read_wide_from_synapse(syn, b$syn_id)

      if (tolower(modality) == "phospho") {
        message(" Preprocessing phospho table: drop blank Gene.Names; build site ids")
        wide <- wide %>% dplyr::filter(!is.na(.data$Gene.Names), .data$Gene.Names != "")
      } else if (tolower(modality) == "global") {
        message(" Preprocessing global table: split multi-symbol rows (Genes by ';')")
        wide <- tidyr::separate_rows(wide, Genes, sep = ";") |>
          dplyr::mutate(Genes = trimws(Genes)) |>
          dplyr::filter(!is.na(Genes) & Genes != "")
      } else if (tolower(modality) == "rna") {
        message(" Preprocessing RNA table: (using merged Salmon matrix; samples start after gene_id/gene_name)")
      }

      feats <- build_ids(wide)
      ok <- !is.na(feats) & feats != ""
      if (!all(ok)) {
        message(" Dropping ", sum(!ok), " empty/NA feature IDs before SE construction.")
        wide  <- wide[ok, , drop = FALSE]
        feats <- feats[ok]
      }

      fallback_col <- if (tolower(modality) == "rna") 3 else 5
      first_col <- if (!is.null(b$value_start_col)) b$value_start_col else detect_value_start_col(wide, fallback = fallback_col)

      fnmap <- parse_fnames(colnames(wide)[first_col:ncol(wide)], b$fname_aliquot_index, b$cohort)
      if (tolower(modality) == "rna" && any(is.na(fnmap$aliquot))) {
        message("  ! Note: ", sum(is.na(fnmap$aliquot)), " RNA sample(s) with NA aliquot after parsing (expected for RNA headers).")
      }

      se0  <- make_se(wide, value_start_col = first_col, feature_ids = feats,
                      fnames_df = fnmap, meta = meta, drop_name = drop_name, modality = modality)
      se_n <- normalize_by_modality(se0, modality)
      se_list[[i]] <- se_n

      if (write_outputs) {
        message(" Writing per-batch normalized long table (pre-ComBat)")
        batch_long <- se_to_long(se_n, modality)
        batch_tag  <- paste0("batch", b$cohort)
        batch_path <- file.path(out_dir, paste0(file_stem, "_", batch_tag, "_normalized_long.csv"))
        readr::write_csv(batch_long, batch_path)
        upload_queue <- c(upload_queue, batch_path)
      }
    }

    # Combine & pre-QC
    message("--------------------------------------------------")
    se_combined  <- combine_batches_intersection(se_list)

    message(" Pre-QC plots (PCA & histogram) on combined (pre-ComBat)")
    pre_pc_df  <- pca_df_present_in_all(se_combined)
    pre_pca    <- plot_pca(pre_pc_df, paste0(modality, " samples"), pcols = pcols)
    pre_hist   <- plot_hist(se_combined, paste0(modality, ": value distribution (pre-ComBat)"))
    if (write_outputs) {
      pre_pca_pdf  <- file.path(out_dir, paste0(file_stem, "_preComBat_PCA.pdf"))
      pre_hist_pdf <- file.path(out_dir, paste0(file_stem, "_preComBat_Hist.pdf"))
      ggsave(pre_pca_pdf,  pre_pca,  width = 7, height = 4.5, device = cairo_pdf)
      ggsave(pre_hist_pdf, pre_hist, width = 7, height = 4.5, device = cairo_pdf)
      upload_queue <- c(upload_queue, pre_pca_pdf, pre_hist_pdf)
    }

    # ComBat
    if (isTRUE(do_batch_correct)) {
      message("--------------------------------------------------")
      se_post <- combat_by_cohort(se_combined)
      post_suffix <- "_batchCorrected"
      post_title  <- paste0("Batch-corrected ", modality, " samples")
    } else {
      message("--------------------------------------------------")
      message("Skipping ComBat per do_batch_correct=FALSE; using combined matrix as 'post'.")
      se_post <- se_combined
      post_suffix <- "_noBatchCorrect"
      post_title  <- paste0("Combined ", modality, " samples (no ComBat)")
    }

    # Exports
    message(" Building long tables")
    long_pre  <- se_to_long(se_combined,  modality) |>
      dplyr::filter(is.finite(correctedAbundance))
    long_post <- se_to_long(se_post, modality)

    if (write_outputs) {
      path_pre  <- file.path(out_dir, paste0(file_stem, "_preComBat_long.csv"))
      path_post <- file.path(out_dir, paste0(file_stem, post_suffix, ".csv"))
      write_csv(long_pre,  path_pre)
      write_csv(long_post, path_post)
      upload_queue <- c(upload_queue, path_pre, path_post)
    }

    # Post ComBat/QC
    message(" Post-QC plots (PCA & histogram)")
    pc_df  <- pca_df_present_in_all(se_post)
    gpca   <- plot_pca(pc_df, post_title, pcols = pcols)
    ghist  <- plot_hist(se_post, paste0(modality, ": value distribution", ifelse(isTRUE(do_batch_correct), "", " (no ComBat)")))
    if (write_outputs) {
      post_pca_pdf  <- file.path(out_dir, paste0(file_stem, ifelse(isTRUE(do_batch_correct), "_PCA.pdf", "_PCA_noComBat.pdf")))
      post_hist_pdf <- file.path(out_dir, paste0(file_stem, ifelse(isTRUE(do_batch_correct), "_Hist.pdf", "_Hist_noComBat.pdf")))
      ggsave(post_pca_pdf,  gpca,  width = 7, height = 4.5, device = cairo_pdf)
      ggsave(post_hist_pdf, ghist, width = 7, height = 4.5, device = cairo_pdf)
      upload_queue <- c(upload_queue, post_pca_pdf, post_hist_pdf)
    }

    # Pack results for return
    results <- list(
      se_batches     = se_list,
      se_combined    = se_combined,
      se_corrected   = if (isTRUE(do_batch_correct)) se_post else NULL,
      se_post        = se_post,                # always populated (corrected or not)
      did_combat     = isTRUE(do_batch_correct),
      long_pre       = long_pre,
      long_post      = long_post,
      pca_df_pre     = pre_pc_df,
      pca_df_post    = pc_df,
      plots          = list(pre_pca = pre_pca, pre_hist = pre_hist, pca = gpca, hist = ghist),
      files          = if (write_outputs) list(queued = upload_queue) else list()
    )

    # FINAL STEP: Uploads
    if (write_outputs && !is.null(upload_parent_id)) {
      perform_uploads(upload_queue, syn, upload_parent_id)
    } else if (!write_outputs) {
      message("write_outputs=FALSE — skipping all writes/uploads.")
    } else {
      message("No upload_parent_id provided — skipping uploads.")
    }

  }, error = function(e) {
    message("ERROR: run_modality() failed for ", modality, ". No uploads were attempted.")
    message("       ", conditionMessage(e))
    message("------- DEBUG SNAPSHOT -------")
    message(" traceback:"); print(sys.calls())
    message(" sessionInfo():"); print(utils::sessionInfo())
    message("------- END DEBUG  -----------")
    stop(e)
  })

  message("run_modality() finished successfully for: ", modality)
  results
}
