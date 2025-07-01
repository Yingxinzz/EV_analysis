# SVAtlas: A Comprehensive Single Extracellular Vesicle Omics Resource

## Overiew

This  includes an automated single-EV analysis pipeline for clustering and subpopulation visualization based on input expression matrices. One could obtain the data from the Download of our website. Outputs include a t-SNE plot of EV subclusters and clustering outputs in JSON format.


## Workflow for Single Extracellular Vesicle (sEV) Analysis

The project includes the following key modules:

### 1.Data Organization & Matrix Transformation

(1) Convert long-format data (Sample, sEV, Protein, Value) into wide-format expression matrix: proteins as rows, sEVs as columns.

(2) Create Seurat objects for each sample based on the wide-format matrices.

### 2.Merge Seurat Objects

(1) If there are two or more Seurat objects, merge them using the merge() function to create a single combined Seurat object.

(2) Retain sample IDs as metadata (orig.ident) to facilitate downstream batch correction.

### 3.Normalization, Dimensionality Reduction, Batch Correction

(1) Apply SCTransform normalization to the merged Seurat object, retaining all genes rather than only variable features.

(2) Perform batch effect correction using Harmony based on the orig.ident metadata.

(3) Save convergence plots of the Harmony algorithm for quality assessment

### 4.Dimensionality Reduction, Visualization, Clustering

(1) Run dimensionality reduction methods (e.g., PCA, UMAP, t-SNE) on batch-corrected data.

(2) Perform clustering (e.g., FindNeighbors(), FindClusters()) to identify groups of sEVs with similar profiles.

(3) Visualize clusters and sample distributions using appropriate plotting functions (e.g., DimPlot()).


## Dependencies

data.table v1.17.0

dplyr v1.1.4

tidyr v1.3.1

stringr v1.5.1

Seurat v5.2.1

readr v2.1.5

purrr v1.0.4

harmony v1.2.3

ggplot2 v3.5.2

scales v1.3.0