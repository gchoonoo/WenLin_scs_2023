# Single Cell Analysis - Humanized StRG T cells

## Package versions
* R version 3.6.3
* Seurat_3.1.4
* SeuratWrappers_0.3.0
* scProportionTest_0.0.0.9000
* monocle3_1.0.0
* ggplot2_3.3.5 
* dplyr_1.0.7

## Reproduction Instructions
1. Download .RData files from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223026
2. Run Rscript WenLin_scs_2023.R in R to reproduce single cell figures from the paper

## Code Functionality
The Rscript functions to reproduce the following from Figure 5 and Figure S4:
* Single cell clustering UMAPs
* Feature UMAPs
* Feature violin plots
* Proportion fold change plot
* Trajectory UMAP
* Trajectory boxplot
* TCR clonotype barplot
* Shannon diversity plot
