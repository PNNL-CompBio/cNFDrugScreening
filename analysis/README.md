


# cNF multi-omics analysis pipeline overview
Purpose: Reference for running the cNF analysis end-to-end and understanding what each stage produces.  
  
- Run order: Lists the recommended notebook/script sequence and what each file sources (dependencies).  
  
- Data processing: Summarizes cohort-wise preprocessing, normalization, and ComBat batch correction for RNA, global proteomics and phosphoproteomics.  
  
- Exploratory outputs: Describes the PCA figures generated from batch-corrected data (global + phospho).  
  
- Biomarker discovery: Outlines how drug response matrices and omics feature matrices are constructed, filtered to shared samples, and correlated.  
  
- Pathway interpretation: Summarizes direction-aware enrichment (leapR / GSEA) and the key plots produced across drugs, modalities, and metrics.  
  
## Quick run order (and what each file sources)

1) **01_run_normalize_omics.Rmd**  
   - *Sources:* `cNF_helper_code.R`, `01_normalize_batchcorrect_omics.R`  
   - Runs per-modality normalization and ComBat; writes long tables/plots if enabled.

2) **02_analyze_modality.Rmd**  
   - *Sources:* `cNF_helper_code.R`, `02_analyze_modality_correlations.R`  
   - Builds drug/feature matrices, modality correlations.

3) **03_pathway_enrichment.Rmd**  
   - *Sources:* `cNF_helper_code.R`, `03_leapr_biomarker.R`  
   - LeapR pathway enrichment and plots for all.
  

## Normalized Global / Phosphoproteomics Data for Clustering via PCA

The datasets analyzed in this study included both global proteomics and phosphoproteomics
measurements collected from two experimental cohorts of cNF (cutaneous neurofibroma) organoid
samples. Each cohort comprised multiple tumor specimens derived from different patients.

For each dataset (global proteins and phosphosites), raw intensity values were preprocessed to
ensure comparability across cohorts:

- **Missing values and zeros:** Zero measurements were replaced with missing values (NA) to prevent skew during normalization.
- **Feature filtering:** Proteins or phosphosites absent in more than 50% of samples were removed to reduce noise and improve robustness.
- **Normalization:** Remaining values were log2-transformed with a small offset (to prevent log(0)/undefined) and scaled using a modified z-score (median centering and median absolute deviation scaling). This approach places samples on a common scale while down-weighting outliers.
- **Cohort separation:** Each cohort was processed independently to account for technical differences, then merged into a single combined dataset.
- **Batch correction:** To address systematic differences between cohorts, ComBat (sva R package) was used, producing batch-adjusted abundance matrices for both global proteomics and phosphoproteomics.

The resulting normalized and corrected data provided a harmonized view of protein and phosphosite
abundance across all patients and cohorts, enabling joint exploratory analysis such as PCA (Figures:
phosphoCorrectedPCA.pdf and globalCorrectedPCA.pdf).

#### Figure: `phosphoCorrectedPCA.pdf`

This figure shows the results of a principal component analysis (PCA) performed on the
phosphoproteomics data after normalization and batch correction using ComBat. Prior to correction,
samples clustered primarily by cohort, reflecting a strong batch effect. After adjustment, the major
sources of variation correspond more closely to biological factors, such as patient identity and
tumor replicate. Each point represents an individual sample, colored by patient and shaped by tumor
designation, enabling visualization of patient-specific clustering patterns. The correction
substantially reduced separation by cohort, suggesting that technical batch effects were effectively
mitigated and that the remaining signal more likely reflects underlying biological differences in
phosphoproteome profiles.

#### Figure: `globalCorrectedPCA.pdf`

This figure presents a PCA of the global proteomics dataset following the same normalization and
ComBat batch correction procedure. Similar to the phosphoproteomics analysis, initial clustering was
dominated by cohort effects; however, correction aligned the data so that patient identity, rather
than cohort, explained the primary axes of variation. Each sample is again colored by patient and
shaped by tumor replicate, allowing assessment of within-patient reproducibility and between-patient
divergence. The improved alignment indicates that batch effects were successfully minimized, and the
corrected data provide a more reliable basis for downstream analyses of patient- and tumor-specific
proteomic signatures.

### Biomarker Evaluation using Batch-Corrected Global / Phosphoproteomics Data

For this biomarker evaluation analysis, we brought together four different data types collected on
the same patient-derived samples:

- **Drug sensitivity measurements:** single-dose viability values (and some full dose–response curves) across a large panel of compounds. For drugs with multiple doses enabling curve fitting, dose-response fit_auc (area under the fitted viability curve; lower AUC = higher sensitivity) was also analyzed.
- **RNA sequencing:** transcript abundance data (to be harmonized into long format in subsequent steps).
- **Global proteomics:** normalized protein abundance values, batch-corrected across cohorts.
- **Phosphoproteomics:** normalized and batch-corrected site-specific phosphorylation abundances.

To enable cross-modality analyses, all datasets were filtered to retain only those specimens present
in both the drug response data and at least one molecular modality. This produced a shared set of
samples suitable for correlation and biomarker discovery analyses.

- **Proteomics:** Both global and phosphoproteomic data were restructured into sample-by-feature matrices, where rows represent specimens and columns represent either proteins (global) or phosphosites (phospho). Abundance values correspond to batch-corrected, log-scaled measurements.
- **Drug sensitivity:** Viability data (expressed in relative units of surviving fraction at given doses) were reformatted into a specimen-by-drug matrix, with average values across replicates used where multiple measurements existed. Separately, a fit_auc specimen×drug matrix was constructed by averaging per-drug fit_auc across replicates, retaining only drugs present in syn69947322.
- **Filtering:** Only drugs with complete measurements across all available specimens were retained for certain analyses, yielding a consistent comparison set ("full drugs"). For fit_auc, additional coverage summaries (mean, SD, n specimens per drug) were calculated to support ranking and visualization.
- **Drug Count:** 238

<details>
<summary>Drug list</summary>

- Abemaciclib
- ABT-737
- Adagrasib
- Adavosertib
- Afatinib
- Alectinib
- Alisertib
- Alpelisib
- Alvespimycin
- Alvocidib
- Anlotinib
- Apitolisib
- ARS-1620
- Avapritinib
- AVL-292
- Avutometinib
- Axitinib
- Barasertib
- Batoprotafib
- Belinostat
- Belzutifan
- BI-3406
- BI-847325
- BI-D1870
- Bimiralisib
- Binimetinib
- Birinapant
- BLU9931
- BMS-536924
- Bosutinib
- Brigatinib
- Brivanib
- Brukinsa
- Brustaol
- Cabozantinib
- Calquence
- Capivasertib
- Capmatinib
- CBL0137
- CCT-018159
- Cedazuridine
- Ceralasertib
- Cerdulatinib
- Ceritinib
- CFT8634
- Cobimetinib
- Copanlisib
- CPI-613
- Crenolanib
- Crizotinib
- Dabrafenib
- Dacomitinib
- Danusertib
- Daporinad
- Daprodustat
- Dasatinib
- Defactinib
- Derazantinib
- Digoxin
- Dinaciclib
- Dovitinib
- Doxorubicin
- Duvelisib
- Eganelisib
- Elimusertib
- Enasidenib
- Encorafenib
- Enitociclib
- Ensartinib
- Entinostat
- Entrectinib
- Enzastaurin
- Everolimus
- Fadraciclib
- Famotidine
- Fedratinib
- Fimepinostat
- Fisogatinib
- FRAX597
- Futibatinib
- Galunisertib
- Ganetespib
- Gefitinib
- Geldanamycin
- Gilteritinib
- Glasdegib
- GNE-617
- GSK2256098
- H 89 2HCl
- Ibrutinib
- Idasanutlin
- IKK-16
- Imatinib
- INCB188053
- INCB191856
- Infigratinib
- Ipatasertib
- Ivosidenib
- JNJ-7706621
- Ketotifen
- KO-947
- KW-2478
- Lapatinib
- Larotrectinib
- Lenvatinib
- Lestaurtinib
- Letrozole
- Lifirafenib
- Linsitinib
- Lorlatinib
- Losartan
- Losmapimod
- Lovastatin
- Luminespib
- LY231514
- LY3023414
- Manumycin A
- Merestinib
- Metformin
- Methotrexate
- Midostaurin
- Mirdametinib
- MK-2206
- ML 210
- MLN2480
- Molibresib
- Napabucasin
- Navitoclax
- Nedisertib
- Nedometinib
- Neratinib
- Nilotinib
- Nintedanib
- Nirogacestat
- NT219
- NVP-BEZ235
- Octreotide
- Olaparib
- Olutasidenib
- Omipalisib
- Onalespib
- Onametostat
- Onvansertib
- Osimertinib
- Pacritinib
- Palbociclib
- Pan-RAS-IN-1
- Panobinostat
- Pazopanib
- PCNA-I1
- Pelabresib
- Pemigatinib
- Pemrametostat
- PF-562271
- Pimozide
- Pirtobrutinib
- PLX-3397
- Ponatinib
- Pralsetinib
- Prexasertib
- PT2385
- Quizartinib
- Ranitidine
- Rapamycin
- Ravoxertinib
- Regorafenib
- Repotrectinib
- Retaspimycin
- Ribociclib
- Ripretinib
- RMC-6236
- RO4929097
- Rogaratinib
- Romidepsin
- Roxadustat
- Rucaparib
- Ruxolitinib
- Sapanisertib
- SAR405838
- Seclidemstat
- Selinexor
- Selpercatinib
- Selumetinib
- Sertraline
- Siremadlin
- Sitravatinib
- SNX-2112
- SNX-5422
- Sonidegib
- Sorafenib
- Sotatercept
- Sotorasib
- Subasumstat
- Sulfasalazine
- Sunitinib
- Surufatinib
- SW106065
- Tadalafil
- TAE226
- TAK-243
- Tanespimycin
- Tazemetostat
- TED-347
- Tegavivint
- Telaglenastat
- Temozolomide
- Temsirolimus
- Tepotinib
- THZ1
- Ticlopidine
- Tipifarnib
- Tivantinib
- Tivozanib
- TK216
- Tofacitinib
- Tomivosertib
- Tovorafenib
- Tozasertib
- Trametinib
- Tretinoin
- Trilaciclib
- Triptolide
- Ulixertinib
- Umbralisib
- UNC2025
- Unesbulin
- Vactosertib
- Vandetanib
- Vemurafenib
- Venetoclax
- Vismodegib
- Vistusertib
- Vociprotafib
- Volasertib
- Vorinostat
- VS-6766
- Y-27632
- Zanzalintinib
</details>

### Exploratory Analysis of Drug Responses

We performed initial exploration of the drug dataset to understand variability and efficacy across
the compound panel.

### Correlation of molecular features with drug response

Spearman correlations were computed between every drug (response profile) and every molecular
feature (protein or phosphosite).

Significance values were estimated via permutation-based correlation tests, with multiple testing
correction applied (FDR).

Features were separated into positive correlations (higher abundance/phosphorylation associated with
higher viability, i.e. resistance) and negative correlations (higher abundance/phosphorylation
associated with lower viability, i.e. sensitivity).

#### Figure: `most_efficacious_<metric>.png/pdf`

This scatterplot highlights the subset of drugs that were most efficacious across the patient
samples, defined as compounds with an average cell viability below 0.5 (i.e., less than 50%
survival). Each point represents a drug, positioned by its mean viability (y-axis) and labeled along
the x-axis. Point size reflects the variability in response across specimens, while point color
indicates the number of samples tested. Several compounds demonstrate both strong overall activity
and consistent performance across patients, suggesting broad-spectrum effectiveness.

#### Figure: `most_variable_<metric>.png/pdf`

This scatterplot highlights drugs with the highest variability in response across specimens, defined
as those with a standard deviation greater than 0.15. Here, each point again represents a drug, with
mean viability on the y-axis and drug identity on the x-axis. Point size encodes variability, and
color denotes the number of samples measured. Unlike the most efficacious plot, these compounds are
not necessarily the most potent but instead show strong heterogeneity between patients. Such drugs
may provide the greatest opportunity for biomarker discovery, as differences in response are more
likely to be explained by underlying molecular features.

#### Figure: `<omics>_drug_heatmap_large_viability.pdf`

This heatmap displays drug response values (viability) for the subset of drugs measured consistently
across all specimens ("full drugs"). Rows correspond to patient samples and columns to drugs.
Hierarchical clustering of both rows and columns highlights patterns of similarity, revealing groups
of patients with shared sensitivity profiles as well as clusters of compounds with correlated
activity. The visualization provides a global overview of drug response heterogeneity across the
cohort, serving as a baseline reference for linking molecular features to therapeutic sensitivity.

#### Figure: `<metric>_cor_features_by_drug.pdf`

Bar charts summarizing, for each drug, the number of molecular features (proteins or phosphosites)
that correlate with response at FDR < 0.25, separated by direction (features associated with
resistance vs sensitivity) and faceted by data modality. These counts provide a quick sense of which
drugs show the richest correlational signal.

**Additional fit_auc figures:**
- Figure: most_sensitive_auc.pdf — scatter plot of drugs ranked by mean fit_auc (y-axis); point size = SD of fit_auc across specimens; color = number of specimens.
- Figure: most_variable_auc.pdf — scatter plot of drugs ranked by SD of fit_auc (y-axis); point size = number of specimens.

### Functional enrichment of correlated features

Using the leapR package, correlated features were ranked and tested for enrichment against gene set
collections..

Enrichment was conducted in a direction-aware manner:

"Top" (positively correlated features) highlight pathways enriched among resistance-associated
features.

"Bottom" (negatively correlated features, flipped in rank) highlight pathways enriched among
sensitivity-associated features.

Results were summarized at multiple levels:

- **Pathways per drug:** Number of significant pathways identified for each drug, separated by resistant vs sensitive associations.
- **Pathways across drugs:** Top 15 most recurrent pathways across the full drug panel, faceted by proteome vs phosphoproteome.
- **Drug-specific profiles:** Bar charts of the top 15 pathways per drug and modality, annotated with significance stars (* <0.05, ** <0.01, *** <0.001)

#### Figure: `pathways_across_drugs_top15.pdf`

Pathway-level summary of enrichment results across the full drug panel. Bars show the top 15
pathways (per modality) most frequently enriched among significant feature–response correlations
(FDR < 0.05), with direction indicating whether pathways are associated with relative resistance
(features positively correlated with viability) or sensitivity (features negatively correlated with
viability). This highlights recurrent biological programs linked to drug response.

Generalized figures:

Generalized pathway enrichment figures were generated for the two most efficacious drugs, the two
most variable drugs, and Onalespib.

#### Figure: `pathways_<Drug>_<omics>_<efficacy>_<drug_response>_top15.pdf`

Each pathway figure (pathways_<Drug>_<omics>_<efficacy>_<drug_response>_top15.pdf) displays the top
15 enriched pathways whose molecular features are most strongly associated with the drug response
profile for <Drug>. The <omics> field indicates whether enrichment was performed on global
proteomics (global) or phosphoproteomics (phospho) data. The <efficacy> field specifies whether the
pathways are associated with relative resistance (resistant, features positively correlated with
drug response) or relative sensitivity (sensitive, features negatively correlated with drug
response). The <drug_response> field denotes whether correlations were computed using single-dose
viability (viability) or dose-response fit_auc (fit_auc). In all cases, pathways are ranked by
−log10(FDR), and significance is annotated using standard thresholds (* <0.05, ** <0.01, ***
<0.001). These figures collectively summarize the biological programs whose protein abundance or
phosphorylation levels track with either enhanced or diminished drug effect across patient samples.

Note on Onalespib (and other uniformly potent drugs):
Onalespib was one of the most efficacious compounds in the panel, but because nearly all samples
were highly sensitive, the variability across specimens was minimal. Correlation-based enrichment
relies on inter-patient heterogeneity; without it, even strong biological effects cannot be linked
robustly to specific pathways. This explains why Onalespib yields few significant pathways despite
its overall potency.

### RNA-Seq / Differential Expression Analysis / Gene Set Enrichment

Samples and preprocessing.
RNA-seq quantifications (Salmon gene-level quant.genes.sf) were retrieved for two conditions from
the same organoid line (NF0021-T1):

1. untreated organoids
2. Onalespib (1 µM) treated organoids

Raw read counts were assembled into a gene × sample matrix and filtered to retain genes with total
counts ≥ 10 across the two samples. Size-factor normalization was performed with DESeq2 to obtain
normalized counts.

### Effect-size estimation (no replication)

Because we have one sample per condition, formal differential testing is not performed. Instead, we
report exploratory effect sizes: per-gene log2 fold-change (log2FC) computed from normalized counts
as log2FC = log2((Onalespib + 1) / (Untreated + 1)). These values support ranking and visualization,
and hypothesis generation but should be interpreted as descriptive (no p-values).

### Ranked list for enrichment

Symbols were converted to ENTREZ ID; when multiple rows mapped to the same ENTREZ, we used the
median log2FC per ENTREZ for enrichment.

Gene set enrichment analysis (GSEA).
We ran GSEA using clusterProfiler using:

MSigDB Hallmark (H) via msigdbr (concise, non-redundant "meta-pathways"),

GO Biological Process (BP) via gseGO

Both dotplots display GeneRatio on the x-axis (fraction of the pathway's genes that are enriched)
and count as point size.

#### Figure: `heatmap_top30_all_log2fc.pdf`

Top 30 transcripts by |log2FC| (including unmapped IDs).
A clustered heatmap of the strongest responders (by absolute log2FC) using normalized counts
centered per gene. Labels use symbols when available, otherwise RefSeq base IDs. This heatmap
preserves completeness and has all of the true top changes.

#### Figure: `heatmap_top30_mapped_log2fc.pdf`

Same selection logic but restricted to symbol-mapped genes for the heatmap. This is much more
interpretable and cleaner to share but is technically missing the top differences.

#### Figure: `gsea_dot_GSEA-MSigDB_Hallmark.pdf`

Dotplot of MSigDB Hallmark GSEA results
This figure displays enriched Hallmark gene sets ranked by normalized enrichment score, with points
colored by enrichment direction, sized by the number of leading-edge genes, and positioned by
GeneRatio.

#### Figure: `gsea_dot_GSEA-GO_BP_.pdf`

Dotplot of GO Biological Processes results
The figure shows enriched GO BP terms from the ranked gene list, with points colored by enrichment
direction, sized by leading-edge gene count, and positioned by GeneRatio.

