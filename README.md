README: Analysis Pipeline for Human Lymph Node-on-a-Chip

This repository contains the custom R and MATLAB scripts used for the bioinformatic analysis, cellular subset quantification, and statistical modeling presented in the manuscript: "A Human Lymph node-on-a-chip for Personalized Evaluation of Vaccine Immunogenicity".

1. Code Availability 
Custom Code: All custom scripts used in this study (R for scRNA-seq, MATLAB for statistical modeling) are provided.

Repository: The code is deposited in this public GitHub repository: [https://github.com/fs2698-web/LN-on-a-chip-manuscript-code]

Versioning: This release (v1.0) corresponds to the version used during the peer-review process.

2. Software & Environment 
Operating System: Windows 11 (Tested).

R Environment (>= 4.2.0):

CRAN: Seurat, ggplot2, dplyr, rstudioapi, Cairo, pheatmap, patchwork, scales, data.table.

Bioconductor: CellChat, slingshot, fgsea, msigdbr, edgeR.

MATLAB Environment: Required Statistics and Machine Learning Toolbox.

Installation: R scripts are equipped with automatic dependency checks to prompt the installation of missing CRAN packages.

3. Hardware Requirements & Runtime 
Hardware: Minimum 16GB RAM (32GB recommended for large RDS objects); Multi-core CPU.

Expected Runtime:

R scripts (UMAP, DE plots): 2–10 minutes.
MATLAB scripts: < 1 minute.


4. Data Preparation
Input data objects (~1.15 GB) required to run the code are available via GitHub Releases:

Download: [Access Data Assets (v1.0)](https://github.com/fs2698-web/LN-on-a-chip-manuscript-code/releases/tag/v1.0)

Setup: Place all .rds and .xlsx files into the data_demo/ folder as specified above.

5. How to Run
For R: Open RStudio -> Open a script in code-for-LN/ -> Session -> Set Working Directory -> To Source File Location -> Run.

For MATLAB: Set code-for-LN/ as the current folder -> Run the .m scripts.

Expected Output: Results (PDF/PNG/Excel) will be saved in the output_example/ folder.

6. Reproducibility & Variations
UMAP Coordinates: UMAP visualization is inherently stochastic. Due to differences in operating systems and underlying libraries (e.g., uwot), re-running the scripts may result in minor coordinate shifts compared to the manuscript. However, the biological clustering and cell-type relationships remain consistent.

7. License 
This code package is distributed under the MIT License.
