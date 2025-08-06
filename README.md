NeuroPhenol: Neurodegenerative DEG & Phenolic CTD Workflow
A fully automated R-based pipeline for

identifying differentially expressed genes (DEGs) in Huntington’s disease vs. non‐demented controls

comparing HD DEGs with Alzheimer’s DEGs (volcano plots, Venn diagrams, master tables)

visualizing top shared DEGs via ComplexHeatmap

extracting and summarizing phenolic compound–gene interactions from the Comparative Toxicogenomics Database (CTD)

Repository Structure
.
├── data/
│   ├── CTD_gene_cgixns_1754297452048.csv       # Raw CTD chemical–gene interactions
│   ├── neurophenol_analysis_session.RData      # R workspace with `expr` and `pheno_mat`
│   └── (other raw expression & phenotype files)
│
├── results/
│   ├── HD_vs_Control/                          # HD differential expression outputs
│   └── phenolic_ctd/                           # CTD phenolic compound processing
│
├── scripts/
│   ├── analysis_hd_vs_control_and_compare.R    # DE analysis, volcano, Venn, master tables
│   ├── heatmap_top20_common.R                  # ComplexHeatmap of top 20 shared DEGs
│   └── phenolic_ctd_processing.R               # Filter CTD for phenolics & export summaries
│
└── README.md                                   # This file
Prerequisites
R ≥ 4.0.0

CRAN packages: limma, ggplot2, ggrepel, ggvenn, VennDiagram, dplyr, readr, stringr, tidyr, pheatmap

Bioconductor packages: org.Hs.eg.db, clusterProfiler, msigdbr, enrichplot, ReactomePA, ComplexHeatmap, circlize, RColorBrewer

Install dependencies
r
# CRAN
install.packages(c(
  "limma","ggplot2","ggrepel","ggvenn","VennDiagram",
  "dplyr","readr","stringr","tidyr","pheatmap"
))

# Bioconductor
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "org.Hs.eg.db","clusterProfiler","msigdbr","enrichplot",
  "ReactomePA","ComplexHeatmap","circlize","RColorBrewer"
))
Workflow & Usage
1. Differential Expression & Comparison
Run the main DE pipeline:

bash
Rscript scripts/analysis_hd_vs_control_and_compare.R
Inputs: data/neurophenol_analysis_session.RData (must contain expr, pheno_mat) Alzheimer’s DEG CSV (from previous AD analysis) in results/HD_vs_Control

Outputs (in results/HD_vs_Control/):

deg_hd_all.csv (full limma results)

deg_hd_significant.csv (FDR < 0.0005 & |log2FC| > 1)

volcano_hd.png

venn_alz_hd.png

master_deg.csv, significant_master_deg.csv

sessionInfo_HD.txt

2. Heatmap of Top 20 Shared DEGs
Generate a ComplexHeatmap of the top 20 Common DEGs:

bash
Rscript scripts/heatmap_top20_common.R
Inputs:

Z-scored expression matrix (mat_top20_z) and annotation objects (annot_row, annot_col) created by the DE script

Output: Top20_CommonGenes_ComplexHeatmap.pdf in your working directory

3. Phenolic CTD Processing
Filter CTD interactions for phenolic compounds and export summaries:

bash
Rscript scripts/phenolic_ctd_processing.R
Input: data/CTD_gene_cgixns_1754297452048.csv

Outputs (in results/phenolic_ctd/):

phenolic_chemical_gene_interactions.csv

compound_gene_map.csv

gene_chemical_counts.csv

compound_gene_edges.tsv

sessionInfo_phenolic_ctd.txt

Results Overview
DE Analysis (results/HD_vs_Control/): Full and filtered DEG tables, volcano & Venn plots, master DEG list

Heatmap: Publication‐ready ComplexHeatmap of top 20 shared DEGs

CTD Phenolics (results/phenolic_ctd/): Interaction tables, summary counts, network‐ready edge list

Reproducibility Tips
Pin package versions with renv

Encapsulate workflow steps using drake or targets

Use continuous integration (GitHub Actions) to automatically install dependencies and run key scripts

License & Contact
This project is released under the MIT License. For questions or contributions, please contact Preetam Banerjee (<banerjeepreetam8759@gmail.com>).
