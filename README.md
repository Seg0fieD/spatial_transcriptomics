# Spatial Transcriptomics Analysis 

This project focuses on the analysis of spatial transcriptomics data using a 10X Visium dataset of a mouse brain, utilizing the Seurat package in R.

## Key Tasks:

1. **Spatial Transcriptomics**:
   - Properties of the slides were analyzed, including the size of the point, distance between points, and the number of points used in the technology.
   - The resolution of the spatial transcriptomics technology was compared to the size of an average eukaryotic cell, considering its implications for data analysis.
   - An initial overview of the data was performed, including visualization of the sample image, coordinates of the spots, and one gene expression matrix.

2. **Spatial Transcriptomics Data in Seurat**:
   - The data was loaded using Seurat's `Load10X_Spatial` function.
   - The Seurat object was inspected to identify where the gene expression data and tissue image are stored.
   - Gene expression for two randomly chosen genes was visualized in the tissue.

3. **Data Preprocessing**:
   - The data was filtered with reasonable cut-off values, justified with appropriate plots. The filtering thresholds were compared to those used in scRNA-seq from Project 1.
   - SCTransform was used for preprocessing, replacing some steps from Project 1.

4. **Dimensionality Reduction, Clustering, and Visualization**:
   - Dimensionality reduction was performed using PCA and UMAP, and a plot was used to explain the choice of dimensions.
   - Clustering was done based on PCA results and visualized in both 2D UMAP space and on the tissue slide.

5. **Differential Expression Analysis**:
   - Differentially expressed genes (DEGs) were analyzed based on clustering and spatial patterning.
   - The top three spatially variable features were identified and visualized on the tissue slide. Comparisons were made with the DEG results from clustering.

6. **Merging the Data**:
   - Two datasets were merged without batch correction, followed by pre-processing, dimensionality reduction, and clustering.
   - Merging with batch correction was performed, and the results were displayed in UMAP space.
   - Batch effects were assessed and resolved, if necessary, for further processing.

7. **Cell-type Identification**:
   - Automatic annotation was performed using Data Integration with scRNA-seq data, and the labels were transferred to the spatial transcriptomics dataset. The annotation was visualized in UMAP space.
   - Manual annotation was performed by identifying marker genes for selected cell types, and gene expression for these markers was visualized on UMAP plots and tissue slides.

8. **Deconvolution**:
   - Deconvolution was performed using the SCDC package. The methodology behind deconvolution was explained, including its necessity and limitations.
   - A reference dataset was prepared, and genes for deconvolution were selected based on DEG analysis.
   - An ExpressionSet object was created for the reference and query data, and deconvolution was performed. The distribution of cell types on the tissue slide was visualized and compared to previous results.


## R Script:
The R script used for analysis is named **Spatial-Transcriptomics-script**.

## Results:
The resulting plots and images can be found in the `png/` directory of this repository.

## References:
For more details, refer to the official Seurat documentation for Spatial Transcriptomics: [Seurat Spatial Vignette](https://satijalab.org/seurat/articles/spatial_vignette.html).
