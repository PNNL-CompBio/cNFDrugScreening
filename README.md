# cNFDrugScreening
This repository manages the code to analyze the results of drug screens in cNF organoids together with omics measurements of the tumor data. 

## Data repository
The data for this code is hosted on Synapse at http://synapse.org/cnfDrugResponse

## Current analysis
We are collecting omics measurements from cNF organoid samples together with drug response data to identify potential biomarkers of drug response. 

## Quick run order (and what each file sources)
1) **01_harmonize_drug_data.Rmd**  
   - *Sources:* `cNF_helper_code.R`  
   - Download & quick QC of drug fits.

2) **02_run_normalize_omics.Rmd**  
   - *Sources:* `cNF_helper_code.R`, `02_normalize_batchcorrect_omics.R`  
   - Runs per-modality normalization and ComBat; writes long tables/plots if enabled.

3) **04_analyze_modality_and_pathway_enrich.Rmd**  
   - *Sources:* `cNF_helper_code.R`, `03_analyze_modality_correlations.R`, `04_leapr_biomarker.R`  
   - Builds drug/feature matrices, modality correlations, leapR pathway enrichment and plots for all.
