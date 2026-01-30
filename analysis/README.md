# cNF multi-omics analysis pipeline overview
Purpose:  cNF batch merging with optional batch correction, drug/omics correlation assessment, pathway enrichment analysis
The notebooks (`.Rmd`) are the primary analysis entry points; the `.R` scripts are sourced helpers.
The `.html` files result from the latest run/knit of the R Markdown files.

## Quick run order (and what each file sources)

1) **[01_run_normalize_omics.Rmd](01_run_normalize_omics.Rmd)**  
   - *Sources:* [../source/00_cNF_helper_code.R](../source/00_cNF_helper_code.R), [../source/01_normalize_batchcorrect_omics.R](../source/01_normalize_batchcorrect_omics.R)  
   - *Goal:* Per-modality preprocessing/normalization and batch correction (when needed).  
   - *Outputs (examples):*
     - Batch-corrected long tables written to synapse.
     - PCA/QC plots for before and after ComBat batch correction (e.g., `globalCorrectedPCA.pdf`, `phosphoCorrectedPCA.pdf`)

2) **[02_analyze_modality.Rmd](02_analyze_modality.Rmd)**  
   - *Sources:* [../source/00_cNF_helper_code.R](../source/00_cNF_helper_code.R), [../source/02_analyze_modality_correlations.R](../source/02_analyze_modality_correlations.R)  
   - *Goal:* Builds the drug response matrix once, then runs per-modality correlations using the long-format omics tables.  
   - *Outputs (written to `outdir`, examples):*
     - Drug-only plots (written once, based on drug response table):
       - `most_efficacious.pdf`
       - `most_variable.pdf`
       - `drug_heatmap_large_viability.pdf` (heatmap of drugs measured in all samples)
     - Modality correlation summary plots (written for each omics type):
       - `<modality>_cor_features_by_drug.pdf`
   - *Notes:*
     - RNA typically runs the end-to-end wrapper (`analyze_modality()`), producing drug-only + modality outputs.
     - Global/phospho typically run modality-only (`analyze_modality_correlations()`), reusing the RNA-derived `drug_mat`.

3) **[03_pathway_enrichment.Rmd](03_pathway_enrichment.Rmd)**  
   - *Sources:* [../source/00_cNF_helper_code.R](../source/00_cNF_helper_code.R), [../source/03_leapr_biomarker.R](../source/03_leapr_biomarker.R)  
   - *Goal:* Direction-aware pathway enrichment (leapR) using correlated features (resistant & sensitive).  
   - *Outputs (examples):*
     - Per-drug pathway barplots (e.g., `pathways_<Drug>_<omic>_<direction>_top15.pdf`)
     - Summary plots across drugs (e.g., top recurrent pathways)

## Helper scripts (sourced by notebooks)
- **[00_cNF_helper_code.R](../source/00_cNF_helper_code.R)**: Shared utilities (Synapse helpers, plotting helpers, common metadata).
- **[01_normalize_batchcorrect_omics.R](../source/01_normalize_batchcorrect_omics.R)**: Normalization / Joining Batches / ComBat batch correction code.
- **[02_analyze_modality_correlations.R](../source/02_analyze_modality_correlations.R)**: Drug summary plots, correlations (includes separated drug-only + modality-only functions).
- **[03_leapr_biomarker.R](../source/03_leapr_biomarker.R)**: Directional feature ranking + leapR enrichment + pathway plotting.
