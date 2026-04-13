# 6110_assignment_04
**Single Cell RNA seq analysis of nasal influenza infection in Mouse**  
Author: Maryanne Ogakwu  
Dataset:   Primary nasal influenza infection rewires tissue-scale memory response dynamics by Kazer et al. 2025

----

## Introduction

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
After quality control, the data was visualised using single cell violin plots to determine the thresholds to be used to filter the data and a regresion line graph to view the qc metrics.  
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

### 4. PCA

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

