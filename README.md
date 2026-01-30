# cNFDrugScreening
This repository manages the code to analyze the results of drug screens in cNF organoids together with omics measurements of the tumor data. 

## Data repository
The data for this code is hosted on Synapse at http://synapse.org/cnfDrugResponse

## Current analysis
We are collecting omics measurements from cNF organoid samples together with drug response data to identify potential biomarkers of drug response.

The current analysis workflow (see `analysis/`) is organized as a small set of notebooks that:
- Join batches and (when needed) batch-correct RNA, global proteomics, and phosphoproteomics (via ComBat), and generate basic QC plots (e.g., PCA).
- Summarize drug response behavior across the cohort (most efficacious / most variable drugs + a cohort-wide drug heatmap).
- Correlate drug response with molecular features per modality (Spearman correlations + FDR), producing per-drug summaries of correlated features.
- Perform pathway enrichment on correlated features (direction-aware enrichment with leapR) to interpret putative mechanisms and biomarkers.
