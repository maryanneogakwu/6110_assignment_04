# ----------- Loading Packages ----------------
library(tidyverse)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(ggplot2)
library(SeuratData)
library(SeuratWrappers)
library(clusterProfiler)
library(org.Mm.eg.db)
library(SingleR)
library(celldex)
library(MAST)
library(slingshot)
library(SingleCellExperiment)

#saving seurat object
saveRDS(seurat_ass4, file = 'seurat_ass4.RDS')
# ------------- Package installation -----------
remotes::install_github("bnprks/BPCells/r")
remotes::install_github("mojaveazure/seurat-disk")
devtools::install_github('satijalab/seurat-data')
BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("scrapper")
BiocManager::install("MAST")
BiocManager::install("slingshot")
getwd()
setwd("C:/Users/marya/Downloads/Assignment_04/Data")

#Loading seurat object
seurat_ass4 <- readRDS("C:/Users/marya/Downloads/Assignment_04/Data/seurat_ass4.rds")
rm(seurat_ass4.1)
# ----------- Part 1: Quality Control Metrics ---------
#Calculate mtRNA percentage
seurat_ass4[["percent.mt"]] <- PercentageFeatureSet(seurat_ass4, pattern = "^mt-")
View(seurat_ass4@meta.data)

#Visualize metadata columns
#This is done to decide which thresholds to use to filter the data
VlnPlot(seurat_ass4, features = c ("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)
#Regression line to view the qc metrics
FeatureScatter(seurat_ass4, feature1= "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_smooth(method = 'lm')
#Subsetting
seurat_ass4 <- subset(seurat_ass4, subset = nFeature_RNA > 200 & percent.mt < 20)
VlnPlot(seurat_ass4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Regression line using subset
FeatureScatter(seurat_ass4, feature1= "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_smooth(method = 'lm')

# ------------ Part 2: Normalisation and Scaling ---------
# Basic log normalization of the data and scaling
seurat_ass4 <- NormalizeData(seurat_ass4, normalization.method = "LogNormalize")
#Idenitfy highly variable features
seurat_ass4 <- FindVariableFeatures(seurat_ass4, selection.method = "vst", nfeatures = 2000)
#Top 50 highly variable genes
top50 <- head(VariableFeatures(seurat_ass4), 50)
top50_plot <- VariableFeaturePlot(seurat_ass4)
LabelPoints(plot = top50_plot, points = top50, repel = TRUE) +
  ggtitle("Top 50 variable genes")

# By default Seurat only scales variable features. Here, we're instead scaling all features (better downstream visualization)
# The scaling phase is also where we would regress out unwanted sources of variation, e.g. cell cycle stage
all.genes <- rownames(seurat_ass4)
#seurat_ass4 <- ScaleData(seurat_ass4, features = all.genes)

# Scale only variable features instead of all genes
seurat_ass4 <- ScaleData(seurat_ass4, vars.to.regress = "percent.mt")
# without specifying features= it defaults to variable features only

# -------------- Part 3: PCA ----------------
# We run a PCA to produce principal components that can be used to cluster our cells
seurat_ass4 <- RunPCA(seurat_ass4, features = VariableFeatures(object = seurat_ass4))
print(seurat_ass4[['pca']], dims =1:5, nfeatures = 5)
#Loadings for PC1 in 500 cells along with the top genes contributing to the PC1
DimHeatmap(seurat_ass4, dims = 1, cells =500, balanced = TRUE)
DimPlot(seurat_ass4, reduction ='pca') + NoLegend() +
  ggtitle("PCA after Normalization and Scaling")

#Elbow plot to determine the dimensionality of the data
ElbowPlot(seurat_ass4)

# ---------- Part 4: Clustering -----------
seurat_ass4 <- FindNeighbors(seurat_ass4, dims = 1:15)
seurat_ass4 <- FindClusters(seurat_ass4, resolution = 0.8)
view(seurat_ass4@meta.data)

#Saving the metadata to a file with mitorchondrial.pt, RNA snn resolution and cluster information
# Extract metadata
metadata <- seurat_ass4@meta.data
head(metadata)
dim(metadata)
#Saving as a CSV 
write.csv(metadata, "metadata.csv", row.names = TRUE)

#Saving as a RDS for reloading in R
saveRDS(metadata, "metadata.rds")

# --------------------- Part 5: UMAP ------------------------
# UMAP with clusters displayed
seurat_ass4 <- RunUMAP(seurat_ass4, dims = 1:15)
DimPlot(seurat_ass4, reduction = "umap", label = TRUE) + ggtitle("UMAP of clusters")


#Finding all markers in all 40 clusters at once
#Find markers for all clusters at once
all.markers <- FindAllMarkers(seurat_ass4, 
                              only.pos = TRUE,     
                              min.pct = 0.25,      
                              logfc.threshold = 0.25,
                              max.cells.per.ident = 500) 

#Save as CSV for easy viewing
write.csv(all.markers, "all_markers.csv", row.names = TRUE)
#Save RDS to avoid rerunning in R
saveRDS(all.markers, "all_markers.rds")
#Reload with
#all.markers <- readRDS("all_markers.rds")
#The parameters were set to only.pos to return only positive markers, min.pct to check genes that are found in 25% of cells and the minimum logfold change was set to 0.25. Max.cells.per.ident is used to downsample each cluster to 500 cells for statistical testing cutting down the RAM usage while still giving reliable markers

# View top 5 markers per cluster
top_5_markers <- all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)

head(top_5_markers)
#Save as CSV for easy viewing
write.csv(top_5_markers, "top_5_markers.csv", row.names = TRUE)
# ----------------------Part 5: Feature Plots ----------------------
VlnPlot(seurat_ass4, features = c("Tyrobp", "Sparc", "Krt18", "Cst3", "S100a5", "Prdx6", "Col1a2", "Flt1", "C1qc"))
#Display feature on the UMAP plot
#Based on your PCs, these are the most informative genes to highlight
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
#min.cutoff = "q10" used to remove background noise
#max.cutoff = "q95" used to prevent outliers skewing colour scale


# --------------------------- Part 6: Automated Annotation ------------------
#Load both mouse references
mouse_immgen <- celldex::ImmGenData()
mouse_rnaseq <- celldex::MouseRNAseqData()
#Save these as files
#Extract normalised counts
seurat_counts <- GetAssayData(seurat_ass4, 
                              assay = "RNA", 
                              layer = "data")

#Run SingleR with both refernces on cluster level instead of cell level to reduce computational complexity
singler_results <- SingleR(test = seurat_counts,
                           ref = list(ImmGen = mouse_immgen,
                                      MouseRNAseq = mouse_rnaseq),
                           labels = list(mouse_immgen$label.main,
                                         mouse_rnaseq$label.main),
                           clusters = seurat_ass4@meta.data$seurat_clusters)

#Map cluster labels back to cells
#Gives the number of clusters per cell type 
seurat_ass4@meta.data$singler_labels <- singler_results$labels[
  match(seurat_ass4@meta.data$seurat_clusters, 
        rownames(singler_results))]
#Gives the number of cells per cell type out of 156,545 cells after mapping 
seurat_ass4@meta.data$singler_labels <- singler_results$labels[
  seurat_ass4@meta.data$seurat_clusters]

#Seeing exactly which cluster equals which cell type
data.frame(
  Cluster = rownames(singler_results),
  CellType = singler_results$labels
)

#Show cell type for every cell
head(seurat_ass4@meta.data$singler_labels)

#Sums to 156545
table(seurat_ass4@meta.data$singler_labels)

#Check predicted cell types
table(singler_results$labels)

#Add labels to Seurat object
seurat_ass4@meta.data$singler_labels <- singler_results$labels

#Visualize on UMAP
DimPlot(seurat_ass4, 
        group.by = "singler_labels",
        label = TRUE,
        repel = TRUE,
        label.size = 3) +
  ggtitle("Automated Cell Type Annotation") +
  theme(legend.position = "right")

#Check annotation quality to show how confident SingleR was in each annotation.
plotScoreHeatmap(singler_results)

#Cluster to Cell Type Table to check how cluster numbers map to cell types
table(Cluster = seurat_ass4@meta.data$seurat_clusters,
      CellType = seurat_ass4@meta.data$singler_labels)


#Save just the metadata that contains cluster assignments and SingleR labels,
metadata.singler <- seurat_ass4@meta.data
write.csv(metadata, "metadata_annotated.csv", row.names = TRUE)

#Save the SingleR results object , contains full SingleR scores showing confidence per cell type
saveRDS(singler_results, "singler_results.rds")

#Save the cluster to cell type mapping
cluster_celltype_map <- data.frame(
  Cluster = rownames(singler_results),
  CellType = singler_results$labels
)
write.csv(cluster_celltype_map, "cluster_celltype_mapping.csv", row.names = FALSE)

# Reload annotations without resaving full object
#metadata <- read.csv("metadata_annotated.csv", row.names = 1)
#seurat_ass_4@meta.data$singler_labels <- metadata$singler_labels
# --------------- Part 7: Manual Annotation ---------

# ------------------- Part 8: Differential Expression Analysis with MAST ---------------------

# Run DE with MAST for one cluster comparing Naive vs D08
cluster0_DE_MAST <- FindMarkers(seurat_ass4,
                                ident.1 = "Naive",
                                ide0nt.2 = "D08",
                                group.by = "orig.ident",
                                subset.ident = 0,
                                test.use = "MAST",
                                min.pct = 0.25,
                                logfc.threshold = 0.25,
                                max.cells.per.ident = 500)

#View results
head(cluster0_DE_MAST)

# Save results
write.csv(cluster0_DE_MAST, "MAST_DE_cluster0_Naive_vs_D08.csv", row.names = TRUE)
saveRDS(cluster0_DE_MAST, "MAST_DE_cluster0_Naive_vs_D08.rds")

#Visualization of top genes using Volcano Plot
cluster0_DE_MAST$gene <- rownames(cluster0_DE_MAST)
cluster0_DE_MAST$significant <- ifelse(cluster0_DE_MAST$p_val_adj < 0.05 & 
                                         abs(cluster0_DE_MAST$avg_log2FC) > 0.5,
                                       "Significant", "Not Significant")

ggplot(cluster0_DE_MAST, aes(x = avg_log2FC, 
                             y = -log10(p_val_adj),
                             color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Significant" = "red", 
                                "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ggtitle("Volcano Plot - Cluster 0 Naive vs D08 (MAST)") +
  xlab("log2 Fold Change") +
  ylab("-log10 Adjusted P-value") +
  theme_classic()
#
top_genes <- cluster0_DE_MAST %>%
  filter(p_val_adj < 0.05) %>%
  top_n(n = 20, wt = abs(avg_log2FC))

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
# ----------------------- ORA Enrichment of top genes ----------------------
#Extracting significant genes from your MAST results
sig_genes_MAST <- rownames(cluster0_DE_MAST[
  cluster0_DE_MAST$p_val_adj < 0.05 & 
    abs(cluster0_DE_MAST$avg_log2FC) > 0.5, ])

cat("Number of significant genes:", length(sig_genes_MAST))

#Running ORA
ora_results <- enrichGO(gene = sig_genes_MAST,
                        OrgDb = org.Mm.eg.db,
                        keyType = "SYMBOL",
                        ont = "BP",             
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        readable = TRUE)

#Checking results
head(as.data.frame(ora_results))

#Visualization of results
dotplot(ora_results, showCategory = 20) +
  ggtitle("ORA - Cluster 0 Naive vs D08")

barplot(ora_results, showCategory = 20) +
  ggtitle("ORA - Cluster 0 Naive vs D08")

#Saving ORA results to file
write.csv(as.data.frame(ora_results), 
          "ORA_cluster0_Naive_vs_D08.csv", 
          row.names = TRUE)

# --------------------  Part 10: Cell Trajectory using Slingshot -----------------------
#Convert Seurat to SingleCellExperiment
sce <- as.SingleCellExperiment(seurat_ass4)

#Running Slingshot using UMAP and cluster labels
sce <- slingshot(sce,
                 clusterLabels = "singler_labels",  # use your annotations
                 reducedDim = "UMAP",
                 start.clus = "Epithelial cells")   # starting cell type
# Create color palette for cell types
cell_types <- unique(sce$singler_labels)
colors <- setNames(rainbow(length(cell_types)), cell_types)

# Plot with correct colors
plot(reducedDims(sce)$UMAP,
     col = colors[sce$singler_labels],
     pch = 16,
     cex = 0.3)
lines(SlingshotDataSet(sce),
      lwd = 2,
      col = "black")
title("Slingshot Trajectory")

# Add legend
legend("bottomright",
       legend = cell_types,
       col = colors,
       pch = 16,
       cex = 0.7)

# Step 4 - Extract pseudotime values
pseudotime <- slingPseudotime(sce)
head(pseudotime)
