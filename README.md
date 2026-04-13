# 6110_assignment_04
**Single Cell RNA seq analysis of nasal influenza infection in Mouse**  
Author: Maryanne Ogakwu  
Dataset:   Primary nasal influenza infection rewires tissue-scale memory response dynamics by Kazer et al. 2025

----

## Introduction
#### Justifications
LogNormalize was selected over SCTransform due to computational constraints imposed by the dataset size of 156,572 cells. SCTransform required approximately 29.3 GB of contiguous RAM during residual matrix computation, exceeding available memory. LogNormalize is a well-established alternative that produces comparable clustering and differential expression results, particularly in large datasets where increased cell count provides sufficient statistical power. Mitochondrial gene expression was additionally regressed out during scaling to compensate for technical variation that SCTransform would otherwise account for.

----
## Methods
### 1. Data Acqusition
Data used in this study was collected from the 2025 study by Kazer et al on  "Primary nasal influenza infection rewires tissue-scale memory response dynamic". The dataset was transformed  into a seurat object containing the metadata and the data, and was provided to the class by the course instrutor.  
Data can be downloaded via the link below.  
```
https://aacgenomicspublic.blob.core.windows.net/public/seurat_ass4.rds
```
### 2. Quality Control 
Qualitity Control checks were performed on the seurat object by calculating the mitochondrial percentage present in the cells and filtering to remove cells with very low numbers of unique genes detected as they could be low quality droplets or empty cells. Cells with very high gene counts will also be removed as these consititute the characteristics of doublets.  
Caclculating mitochondrial contant. Note that the data set contains cells from *Mus musculus* (house mouse), so mitochondrial is denoted as `mt`. 
~~~
seurat_ass4[["percent.mt"]] <- PercentageFeatureSet(seurat_ass4, pattern = "^mt-")
~~~
The low quality cells were removed using the code below. `nFeature_RNA` > 200 used to keep cells with more that 200 genes detected.  This removes empty droplets, cellular debris and dying cells with low gene detection.
`percent.mt < 20` used to keep cells with less than 20% mitochondrial reads. This removed damaged cells that leak cytoplasmic RNA out leaving high concentrations of mitochondrial RNA
~~~
seurat_ass4 <- subset(seurat_ass4, subset = nFeature_RNA > 200 & percent.mt < 20)
~~~
Before and fter quality control, the data was visualised using single cell violin plots to determine the thresholds to be used to filter the data and a regresion line graph to view the qc metrics.  
Single cell violin plot:
~~~
VlnPlot(seurat_ass4, features = c ("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)
~~~
Regression line graph:
~~~
FeatureScatter(seurat_ass4, feature1= "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_smooth(method = 'lm')
~~~

### 3. Normalisation and Scaling
Raw data counts were normalized using `LogNormalize` method in Seurat which divides each gene count by the total counts per cell, multiplies by a scale factor of 10,000 and applies a log transformation. This makes the cells more comparable to each other as each cell captures a different total number of RNA molecules due to technical variation.
~~~
seurat_ass4 <- NormalizeData(seurat_ass4, normalization.method = "LogNormalize")
~~~
The top 2000 variable genes were extracted using `Variance Stabilizing Transformation (VST)` method for the analysis because not all 25,083 genes present in the dataset are informative for distinguishing cell types, for example, housekeeping genes that are expressed equally in all cells add noise without adding information. VST accounts for the relationship between mean expression and variance.
~~~
seurat_ass4 <- FindVariableFeatures(seurat_ass4, selection.method = "vst", nfeatures = 2000)
~~~
Data was scaled using `ScaleData` to prepare the dataset for PCA and Clustering. This step prevents highly expressed genes from dominating PCA simply due to their magnitude. `vars.to.regress = "percent.mt"` was used to remove the influence of mitochondrial gene expression from the data so it doesn't drive clustering. Scaling all 25,083 genes required 29.3 GB of RAM which exceeded my system's available memory causing crashing. To remedy this, scaling was carried out on only variable features in the dataset, as downsteam PCA and UMAP require only the variable features.
~~~
seurat_ass4 <- ScaleData(seurat_ass4, vars.to.regress = "percent.mt")
~~~

### 4. PCA
Principal component analysis (PCA) was performed on the 2,000 most variable genes using `VariableFeatures` to reduce dimensionality prior to clustering.  
~~~
seurat_ass4 <- RunPCA(seurat_ass4, features = VariableFeatures(object = seurat_ass4))
~~~
Visualizations were carried out using a heat map, a scatter plot and an elbow plot. Heatmaps of top gene loadings were examined to confirm that principal components captured biologically meaningful sources of variation. The scatter plot was used to visualise the cells in the PCA space to determine whether the cells are separating into distinct groups, identify batch effects and confirm that prior normalization and scaling worked. The Elbow plot was plotted to determine the number of significant principal components and identify the point at which additional components explained minimal additional variance. 
~~~
DimHeatmap(seurat_ass4, dims = 1, cells = 500, balanced = TRUE)
DimPlot(seurat_ass4, reduction = 'pca') + NoLegend() +
  ggtitle("PCA after Normalization and Scaling")
ElbowPlot(seurat_ass4)
~~~

### 5. UMAP

### 6. Feature Plots

### 7. Automated Annotation

### 8. Manual Annotation

### 9. Differential Expression Analysis using MAST

### 10. ORA Enrichment of top genes in cluster 0

### 11. Cell Trajectory using Slingshot

----
## Results

----
## Discussion


----
## Conclusion

----


## Dependencies and Packages
| Tool | Version | Source | Purpose |
|------|---------|--------|---------|
| R | 4.5.1 | CRAN | Statistical analysis |
|Seurat | 5.4.0 | CRAN | Tool for Single Cell Genomics |
| SeuratObject | 5.3.0 | Data Structures for Single Cell data |
| SingleR | 2.10.0 | Refernce based ScRNA sequence annotation |
| TidyR | 1.3.2 | CRAN | Tidy Messy Data |
| Tidyverse | 2.0.0 | CRAN | Tidyverse |
| ggplot2 | 4.0.2 | CRAN | Data Visualization |
| Clusterprofiler | 4.16.0 | Bio C 3.21 |Universal Erichment tool for interpreting omics data |
| Org.Mm.eg.db | 3.21.0 | Bioconductor | Genome Wide Annotation for Mouse |
| Scrapper | 1.2.1 | Bio C 3.21 | Binding to C++ libraries for SC analysis, dependency for single R |
| Celldex | 1.18.0 | Bio  3.21 | Index of Refernce Cell Type datasets |
| dplyr | 1.2.0 | Bio C 3.21 | Data Manipulation |
| MAST | 1.33.0 | Bio C 3.21 | Model Based Analysis of Single Cell data |
| SingleCellExperiment | 1.30.1 |Bio C 3.21 | S4 classes for Single Cell data |
| Slingshot | 2.16.0 | Bio C 3.21 | Tools for ordering single cell sequences |

## References

