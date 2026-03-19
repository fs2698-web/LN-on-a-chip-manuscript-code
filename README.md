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

Download: Access Data Assets (v1.0)

Setup: Place all .rds and .xlsx files into the data_demo/ folder as specified above.

6. How to Run (运行指引)
For R: Open RStudio -> Open a script in code-for-LN/ -> Session -> Set Working Directory -> To Source File Location -> Run.

For MATLAB: Set code-for-LN/ as the current folder -> Run the .m scripts.

Expected Output: Results (PDF/PNG/Excel) will be saved in the output_example/ folder.

7. Reproducibility & Variations (重现性说明)
UMAP Coordinates: UMAP visualization is inherently stochastic. Due to differences in operating systems and underlying libraries (e.g., uwot), re-running the scripts may result in minor coordinate shifts compared to the manuscript. However, the biological clustering and cell-type relationships remain consistent.

8. License (许可协议)
This code package is distributed under the MIT License.























## Overview
The codebase is designed to integrate multi-modal data from the Lymph Node-on-a-chip (LN-oC) platform, including:

Single-cell RNA sequencing (scRNA-seq): Profiling the transcriptomic landscape of B cells, T cells, and DCs within the biomimetic environment.

On-chip Effluent Analysis: Quantifying immune cell activation kinetics and subset differentiation (e.g., Plasma cells) using automated thresholding.

System Modeling: Correlating on-chip immune metrics with clinical vaccine responses (e.g., serum antibody titers) to evaluate donor-specific immunogenicity.

## Folder structure
LN_chip_code/
- README.md
- data_demo/             
- scripts_R/            
- scripts_MATLAB/       
- output_example/        
- manuscript_mapping/    

## Installation & Dependencies
Operating System: Windows 11
Expected Runtime: Most R scripts: 2–10 minutes.
                  MATLAB scripts: < 1 minute.
R Environment (>= 4.2.0)
Packages: Seurat, CellChat, slingshot, fgsea, ggplot2, dplyr, rstudioapi, Cairo, pheatmap.

Note: Scripts will automatically check and prompt for missing CRAN packages.

MATLAB Environment
Toolboxes: Statistics and Machine Learning Toolbox.

## Data Preparation
The processed data objects (**~1.15 GB** in total) are essential for running the analysis and are hosted as a GitHub Release.

1.  **Direct Download**: Access the data via this link:  
     **[Download Processed RDS Data (v1.0)](https://github.com/fs2698-web/LN-on-a-chip-manuscript-code/releases/tag/v1.0)**
2.  **Files to Download**: Please download all `.rds` files and `.xlsx` files from the **Assets** section.
3.  **Local Setup**: Move all downloaded files into a folder named `data_demo` within your project root (as shown in the Directory Structure above).

## How to Run
For R Scripts:
Open RStudio.

Open any script from the code-for-LN/ folder.

Set Working Directory: Go to Session -> Set Working Directory -> To Source File Location.

Run the script. Results will be saved in ../output_example/.

For MATLAB Scripts:
Open MATLAB.

Navigate to the code-for-LN/ folder.

Run the .m file. The script will automatically look for data in ../data_demo/.

## Expected output
Example output files may be written to the `output_example/` folder, including processed result tables, summary statistics, heatmaps, UMAP plots, pseudotime plots, and B-cell subset percentage tables.


## Data availability note
Large raw data files and intermediate objects may not be included in this package. 

## License
This code package is distributed under the MIT License unless otherwise noted.
