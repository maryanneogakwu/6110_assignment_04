# 6110_assignment_04
**Single Cell RNA seq analysis of nasal influenza infection in Mouse**  
Author: Maryanne Ogakwu  
Dataset:   Primary nasal influenza infection rewires tissue-scale memory response dynamics by Kazer et al. 2025

----

## Introduction
#### Justifications
##### Normalization and Scaling
LogNormalize was selected over SCTransform due to computational constraints imposed by the dataset size of 156,572 cells. SCTransform required approximately 29.3 GB of contiguous RAM during residual matrix computation, exceeding available memory. LogNormalize is a well-established alternative that produces comparable clustering and differential expression results, particularly in large datasets where increased cell count provides sufficient statistical power. Mitochondrial gene expression was additionally regressed out during scaling to compensate for technical variation that SCTransform would otherwise account for.  

##### Differential Expression Analysis  
MAST was selected over DESeq2 as it was specifically developed for single-cell RNA sequencing data. Unlike DESeq2, which was originally designed for bulk RNA sequencing, MAST employs a hurdle model that simultaneously accounts for both the proportion of cells expressing a gene and the expression level among expressing cells, making it well suited to handle the zero-inflated nature of scRNA-seq data where many genes are detected in only a subset of cells. MAST additionally accounts for cellular detection rate as a covariate, reducing confounding technical variation. While DESeq2 has been adapted for single-cell use, its negative binomial model does not account for dropout events as effectively as MAST, making MAST the more statistically appropriate method for this dataset.  

##### Cell Trajectory Analysis  
Slingshot was selected over Monocle3 as it integrates directly with the SingleCellExperiment framework compatible with our existing Seurat workflow, avoiding the additional object conversion required by Monocle3. While Monocle3 supports complex trajectory topologies, its graph learning algorithm is considerably slower on datasets exceeding 100,000 cells, making it less practical for this 156,572 cell dataset. Additionally, Slingshot's ability to incorporate predefined cell type labels from our SingleR annotation as trajectory nodes allowed biologically informed trajectory inference consistent with our annotation results

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
The ElbowPlot showed standard deviation of explained variance levels off around PC 15, meaning PCs beyond this point explain minimal additional biological variation and would add noise to the clustering rather than information.

### 5. Clustering 
A K-nearest neighbor graph was first constructed using the first 15 principal components, selected based on the ElbowPlot, followed by shared nearest neighbor (SNN) graph refinement using `FindNeighbors`. 
~~~
seurat_ass4 <- FindNeighbors(seurat_ass4, dims = 1:15)
~~~
Cells were then clustered using the Louvain algorithm implemented in FindClusters, applied to the SNN graph at a resolution of 0.8, yielding 40 distinct clusters. Cell metadata including cluster assignments were extracted and saved for downstream analysis
~~~
seurat_ass4 <- FindClusters(seurat_ass4, resolution = 0.8)
~~~
### 6. Uniform Manifold Approximation and Projection (UMAP)
UMAP dimensionality reduction was performed using the first 15 principal components consistent with clustering, generating a two-dimensional embedding for visualization of the 40 identified clusters. It preserves both local and global structures of he data making it better for visualization of distinct cell populations, unlike PCA which is linear. The same 15 PCs used for clustering are used here to ensure consistency between your clusters and your visualization.
~~~
DimPlot(seurat_ass4, reduction = "umap", label = TRUE) + 
  ggtitle("UMAP of clusters")
~~~
Cluster marker genes were identified using FindAllMarkers with the Wilcoxon rank sum test, retaining only positive markers expressed in at least 25% of cells within a cluster and exceeding a log2 fold change threshold of 0.25. Each cluster was downsampled to a maximum of 500 cells to manage computational memory requirements while maintaining statistical reliability. The markers identified in this step was used to carry out Manual annotation.
~~~
all.markers <- FindAllMarkers(seurat_ass4, 
                              only.pos = TRUE,     
                              min.pct = 0.25,      
                              logfc.threshold = 0.25,
                              max.cells.per.ident = 500) 
~~~
### 7. Feature Plots
Violin plots and feature plots were generated for genes selected based on their top contributions to the first five principal components, representing the major sources of transcriptional variation in the dataset. Feature plots were overlaid onto the UMAP embedding with a minimum expression cutoff at the 10th percentile to reduce background noise. These genes were selected to represent the diversity of cell types present in the nasal mucosa including myeloid, epithelial, fibroblast, endothelial, and macrophage populations.
~~~
Violin plot
VlnPlot(seurat_ass4, features = c("Tyrobp", "Sparc", "Krt18", "Cst3", "S100a5", "Prdx6", "Col1a2", "Flt1", "C1qc"))
Feature plot
FeaturePlot(seurat_ass4,
            features = c(
              "Tyrobp",
              "Sparc",
              "Krt18",
              "Cst3",
              "S100a5",
              "Prdx6",
              "Col1a2",
              "Flt1",
              "C1qc"),
            min.cutoff = "q10",
            ncol = 3) &
  theme(plot.title = element_text(size = 10, face = "bold"))
~~~
### 8. Automated Annotation
Automated cell type annotation was performed using SingleR with two complementary mouse reference datasets: ImmGenData for immune cell populations and MouseRNAseqData for broader cell type coverage including non-immune populations. To manage computational complexity given the 156,572 cell dataset, annotation was performed at the cluster level rather than the individual cell level, assigning a single label to each of the 40 clusters based on transcriptional similarity to reference profiles. Labels were subsequently mapped back to individual cells.
~~~
singler_results <- SingleR(test = seurat_counts,
                            ref = list(ImmGen = mouse_immgen,
                                       MouseRNAseq = mouse_rnaseq),
                            labels = list(mouse_immgen$label.main,
                                          mouse_rnaseq$label.main),
                            clusters = seurat_ass4@meta.data$seurat_clusters)
~~~
~~~
seurat_ass4@meta.data$singler_labels <- singler_results$labels[
  seurat_ass4@meta.data$seurat_clusters]
~~~
Annotation confidence was assessed using SingleR score heatmaps, and results were validated by comparing cell type assignments with known marker gene expression from the principal component analysis.

### 9. Manual Annotation

### 10. Differential Expression Analysis using MAST
Differential expression analysis between Naive and Day 8 post infection cells within cluster 0 was performed using MAST via Seurat's FindMarkers function. Genes were filtered to those expressed in at least 25% of cells with a minimum log2 fold change threshold of 0.25. Each group was downsampled to 500 cells to manage computational memory requirements.
~~~
cluster0_DE_MAST <- FindMarkers(seurat_ass4,
                                 ident.1 = "Naive",
                                 ident.2 = "D08",
                                 group.by = "orig.ident",
                                 subset.ident = 0,
                                 test.use = "MAST",
                                 min.pct = 0.25,
                                 logfc.threshold = 0.25,
                                 max.cells.per.ident = 500)
~~~
Results were visualized using volcano plots generated with ggplot2, with genes classified as significant if they exceeded an adjusted p-value threshold of 0.05 and an absolute log2 fold change greater than 0.5. The top 20 differentially expressed genes by fold change were labelled using ggrepel. This was done to identify the specific genes driving differences between Naive and D08 which feeds directly into your ORA/GSEA analysis.
~~~
ggplot(cluster0_DE_MAST, aes(x = avg_log2FC,
                             y = -log10(p_val_adj),
                             color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Significant" = "red",
                                "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(data = top_genes,
                  aes(label = gene),
                  size = 3,
                  max.overlaps = 20) +
  ggtitle("Volcano Plot - Cluster 0 Naive vs D08 (MAST)") +
  xlab("log2 Fold Change") +
  ylab("-log10 Adjusted P-value") +
  theme_classic()
~~~

### 11. Over-representation Analysis of top genes in cluster 0
Over-representation analysis (ORA) was performed on significant differentially expressed genes from the MAST analysis using the `enrichGO` function from clusterProfiler, with genes filtered to an adjusted p-value below 0.05 and absolute log2 fold change exceeding 0.5.
~~~
sig_genes_MAST <- rownames(cluster0_DE_MAST[
  cluster0_DE_MAST$p_val_adj < 0.05 & 
    abs(cluster0_DE_MAST$avg_log2FC) > 0.5, ])
~~~
Gene Ontology Biological Process terms were tested against the mouse genome annotation database `org.Mm.eg.db.` Enriched pathways were filtered at a p-value and q-value cutoff of 0.05. 
~~~
ora_results <- enrichGO(gene = sig_genes_MAST,
                         OrgDb = org.Mm.eg.db,
                         keyType = "SYMBOL",
                         ont = "BP",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05,
                         readable = TRUE)
~~~
Results were visualized using dot plots and bar plots displaying the top 20 enriched biological processes, with dot size representing gene count and colour representing statistical significance.
~~~
dotplot(ora_results, showCategory = 20) +
  ggtitle("ORA - Cluster 0 Naive vs D08")

barplot(ora_results, showCategory = 20) +
  ggtitle("ORA - Cluster 0 Naive vs D08")
~~~
### 12. Cell Trajectory using Slingshot
Cell trajectory analysis was performed using Slingshot to infer developmental and activation trajectories across cell types in the nasal mucosa during IAV infection. The Seurat object was converted to a SingleCellExperiment object and Slingshot was applied to the UMAP embedding using SingleR cell type annotations as cluster labels, with epithelial cells defined as the trajectory root consistent with their role as the primary target of IAV infection. 
~~~
sce <- slingshot(sce,
                  clusterLabels = "singler_labels",
                  reducedDim = "UMAP",
                  start.clus = "Epithelial cells")
~~~
Trajectories were visualized by overlaying inferred lineage curves onto the UMAP embedding colored by cell type.
~~~
plot(reducedDims(sce)$UMAP,
     col = colors[sce$singler_labels],
     pch = 16,
     cex = 0.3)
lines(SlingshotDataSet(sce),
      lwd = 2,
      col = "black")
title("Slingshot Trajectory")
legend("bottomright",
        legend = cell_types,
        col = colors,
        pch = 16,
        cex = 0.7)
~~~
Pseudotime values were extracted for each cell representing their position along the inferred trajectory and saved for downstream analysis.
~~~
pseudotime <- slingPseudotime(sce)
~~~
----
## Results

<img width="1200" height="500" alt="Regression QC plot" src="https://github.com/user-attachments/assets/92543a5b-2283-478b-b19b-fd284c6d1d6c" />

**Figure 1:** **Regression QC plot**  
A scatter plot of total UMI counts against the number of detected genes per cell revealed a strong positive correlation of 0.83, indicating that cells with greater sequencing depth detected proportionally more unique genes as expected in high quality scRNA-seq data. Cells from all five timepoints showed consistent distributions, confirming the absence of timepoint-specific technical bias. A small population of cells with disproportionately high UMI counts relative to gene detection was observed, potentially representing doublets or low quality cells that passed the mitochondrial and gene count filtering thresholds applied during quality control.

<img width="1000" height="500" alt="PC Elbow plot" src="https://github.com/user-attachments/assets/b22318c0-4921-4081-9d19-f3ac931ff628" />

**Figure 2:** **Elbow plot of Principal Components**  
The elbow plot revealed a steep decline in standard deviation from PC1 through PC4, capturing the major sources of transcriptional variation including myeloid versus fibroblast (PC1), epithelial versus myeloid (PC2), and macrophage versus neutrophil (PC5) cell type differences. The curve plateaued at approximately PC15, beyond which additional components explained minimal additional variance. Consequently, the first 15 principal components were selected for downstream clustering and UMAP embedding, balancing biological information capture with noise reduction.





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

