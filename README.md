# SKCM Data Fusion analysis, Amaro et al. 2022
script used for Amaro et al. 2022 [1]

This repository stores all scripts for the analysis used in [1], to make all analysis you should create 3 folders (src, data and results).

## src

### Download_data_SKCM.R
Script to download all data used to make plot and Data Fusion (DF) analysis:
  * Download RNA-seq data with RTCGAToolbox [2] and normalization with DESeq2 [3].
  * Methylation data preparation
  * Download maf files for supplementary figure 1 with TCGAbiolinks [4].

### SKCM_data_fusion.R
Script to make DF with R packages and produce all plot, table and analysis of the paper [1]

### JSVD_Pfeffer_et_al_2019.py
This script computes a python version of JSVD as reported in [5] and implemented in [6], it requires pymanopt [7] and other libraries to be installed in order to work (check the first part of the script for more information).
  * Check the path to output and input files in the script and change it or create the appropriate folder structure before running the script. 
  * Just use the command "python3 JSVD_Pfeffer_et_al_2019.py 4" to run the script on an Ubuntu terminal.



## REFERENCES

[1] Amaro, A.; Pfeffer, M.; Pfeffer, U.; Reggiani, F. Evaluation and Comparison of Multi-Omics Data Integration Methods for Subtyping of Cutaneous Melanoma. Biomedicines 2022, 10, 3240. https://doi.org/10.3390/biomedicines10123240

[2] Samur, M.K. RTCGAToolbox: A New Tool for Exporting TCGA Firehose Data. PLoS ONE 2014, 9, e106397, doi:10.1371/journal.pone.0106397.

[3] Love, M.I.; Huber, W.; Anders, S. Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2. Genome Biol. 2014, 15, 550, doi:10.1186/s13059-014-0550-8.

[4] Colaprico, A.; Silva, T.C.; Olsen, C.; Garofano, L.; Cava, C.; Garolini, D.; Sabedot, T.S.; Malta, T.M.; Pagnotta, S.M.; Castiglioni, I.; et al. TCGAbiolinks: An R/Bioconductor Package for Integrative Analysis of TCGA Data. Nucleic Acids Res. 2016, 44, e71, doi:10.1093/nar/gkv1507.

[5] Sato, H. Joint Singular Value Decomposition Algorithm Based on the Riemannian Trust-Region Method. JSIAM Lett. 2015, 7, 13–16, doi:10.14495/jsiaml.7.13.

[6] Pfeffer, M.; Uschmajew, A.; Amaro, A.; Pfeffer, U. Data Fusion Techniques for the Integration of Multi-Domain Genomic Data from Uveal Melanoma. Cancers 2019, 11, 1434, doi:10.3390/cancers11101434.

[7] Townsend, J.; Koep, N.; Weinchwald, S. Pymanopt: A Python Toolbox for Optimization on Manifolds Using Automatic Differentiation. J. Mach. Learn. Res. 17, 1–5.
