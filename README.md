# HackBio-Internship---Stage-2-Bioinformatics-Data-Visualization

Project Overview: 
This repository contains a comprehensive R-based analysis focusing on differential gene expression and clinical data exploration. The project was completed as part of the HackBio Internship (Stage 2), demonstrating the ability to process raw genomic data, perform statistical analyses, and generate publication-quality visualisations.


📊 Key Analysis Components1. 

1. Gene Expression Profiling (HBR vs. UHR)
Clustered Heatmaps: Visualising the top differentially expressed genes using pheatmap to identify expression patterns across biological replicates.
Volcano Plots: Identifying statistically significant genes by plotting $Log_2$ Fold Change against statistical significance ($-log10$ Adjusted P-value).

2. Clinical Data Exploration (Breast Cancer Dataset)
Feature Correlation: A correlation matrix and heatmap of tumour features (radius, texture, area, etc.) to understand multi-collinearity in diagnostic data.
Distribution Analysis: Using Kernel Density Estimate (KDE) plots and scatter plots to differentiate between Malignant and Benign diagnoses.

   Advanced Bioinformatics Panels
Kinetics & Half-life: Reproducing complex research figures involving log-transformed kinetic regimes and cell-type ratio distributions.
Network Analysis: Visualising cell-cell interaction networks using directed graphs with edge weights proportional to interaction strength.

🛠️ Tech Stack & Libraries
1. Language: R

2. Visualization: ggplot2, pheatmap, patchwork, RColorBrewer

3. Data Manipulation: dplyr, tidyr, reshape2

4. Network Analysis: igraph

5. File I/O: readxl

📂 Repository Structure/data:
1. /data: Contains the raw CSV and Excel datasets (HBR/UHR counts, Breast Cancer metrics).

2. /scripts: The complete R script containing logic for all 8 tasks.

3. /outputs: High-resolution (300 DPI) PNG exports of the final multi-panel figures.

🚀 How to RunClone the repository. 
Ensure you have the required R packages installed: install.packages(c("ggplot2", "pheatmap", "dplyr", "patchwork", "igraph")). 
Set your working directory to the project folder. Run the main script to generate the Final_Figure_Part1_2.png and Final_Figure_Tasks.png.

