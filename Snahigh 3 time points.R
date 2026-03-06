##SNAhigh in Omentum at 3 different time points

##Start after line 110 to avoid re- extractign, and re pre processing! Just import and do downstream analysis
library(Seurat)
##Extarcting SNAhigh from different time points
### -------------------------
### Day 14 – SNAhigh (Sample 855)
### -------------------------

day14_snahigh <- subset(
  integrated,
  subset = sample == "855"
)

### -------------------------
### Day 25 – SNAhigh
### -------------------------

day25_snahigh <- subset(
  sna,
  subset = SNA_level_D25 == "SNAhigh"
)

### -------------------------
### Day 40 – SNAhigh Omentum
### -------------------------

day40_snahigh_omentum <- subset(
  snahigh,
  subset = organ == "omentum"
)

### -------------------------
### Quick sanity checks
### -------------------------

table(day14_snahigh$sample)
table(day25_snahigh$SNA_level_D25)
table(day40_snahigh_omentum$organ)

### Cell counts
ncol(day14_snahigh)
ncol(day25_snahigh)
ncol(day40_snahigh_omentum)

### Add timepoint metadata
day14_snahigh$timepoint <- "D14"
day25_snahigh$timepoint <- "D25"
day40_snahigh_omentum$timepoint <- "D40"

### Make cell names unique
day14_snahigh <- RenameCells(day14_snahigh, add.cell.id = "D14")
day25_snahigh <- RenameCells(day25_snahigh, add.cell.id = "D25")
day40_snahigh_omentum <- RenameCells(day40_snahigh_omentum, add.cell.id = "D40")

##Merge the SNAhigh together
sna_timecourse <- merge(
  x = day14_snahigh,
  y = list(day25_snahigh, day40_snahigh_omentum)
)

table(sna_timecourse$timepoint)

DefaultAssay(sna_timecourse) <- "RNA"

sna_timecourse <- NormalizeData(sna_timecourse)

sna_timecourse <- FindVariableFeatures(
  sna_timecourse,
  selection.method = "vst",
  nfeatures = 2000
)

sna_timecourse <- ScaleData(sna_timecourse)
sna_timecourse <- RunPCA(sna_timecourse)
ElbowPlot(sna_timecourse)

sna_timecourse <- FindNeighbors(sna_timecourse, dims = 1:15)
sna_timecourse <- FindClusters(sna_timecourse, resolution = 0.5)

sna_timecourse <- RunUMAP(sna_timecourse, dims = 1:15)

DimPlot(sna_timecourse, label = TRUE, split.by="timepoint")
DimPlot(sna_timecourse, group.by = "timepoint")

saveRDS(sna_timecourse, file = "sna_timecourse_seurat.rds")

table(sna_timecourse$seurat_clusters, sna_timecourse$timepoint)
library(ggplot2)

prop <- prop.table(
  table(sna_timecourse$seurat_clusters, sna_timecourse$timepoint),
  margin = 2
)

prop_df <- as.data.frame(prop)

colnames(prop_df) <- c("cluster", "timepoint", "proportion")

ggplot(prop_df, aes(x = timepoint, y = proportion, fill = cluster)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  ylab("Cell proportion") +
  xlab("Timepoint") +
  labs(fill = "Cluster")




##Run from here to import pre processed file
sna_timecourse <- readRDS("sna_timecourse_seurat.rds")

##FInd Marker Genes to identify what clusters are
sna_timecourse <- JoinLayers(sna_timecourse)
Idents(sna_timecourse) <- "seurat_clusters"
cluster_markers <- FindAllMarkers(
  sna_timecourse,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

#Extract Top 10 marker genes
library(dplyr)

top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)

write.csv(cluster_markers, "SNA_timecourse_cluster_markers.csv")
write.csv(top_markers, "SNA_timecourse_top10_markers.csv")

library(dplyr)

##Dotplot for Top 5 genes per cluster
top5 <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)
top5_genes <- unique(top5$gene)
DotPlot(
  sna_timecourse,
  features = top5_genes,
  group.by = "seurat_clusters"
) +
  RotatedAxis()


#Check for canonical markers
## Second layer of annotation with canonical markers

DotPlot(
  sna_timecourse,
  features = list(
    Macrophages = c("Adgre1","Cd68","Cd14","Lyz2","Mrc1","Marco","Retnla"),
    DCs = c("Itgax","Itgam","Xcr1","Clec9a","Cd209a","Clec10a","Flt3","Zbtb46"),
    T_cells = c("Cd3d","Cd3e","Cd3g","Cd4","Cd8a","Cd8b1","Trbc1","Trbc2"),
    NK = c("Nkg7","Klrc1","Klrc2","Klrb1c","Gzmb","Prf1"),
    B_cells = c("Cd79a","Cd79b","Ms4a1","Cd19","Ighm","Jchain"),
    Fibroblast = c("Col1a1","Col1a2","Dpt","Lum","Pdgfra"),
    Epithelial = c("Krt8","Krt18","Krt14","Krt5"),
    MDSC = c("S100a8","S100a9","Il1f9","Cxcr2")
  ),
  group.by = "seurat_clusters"
) +
  RotatedAxis()

Idents(sna_timecourse) <- "seurat_clusters"
##Rename clusters
new_cluster_ids <- c(
  "0" = "B_cells_1",
  "1" = "Macrophage_1",
  "2" = "NK_cells",
  "3" = "Macrophage_2",
  "4" = "Fibroblasts",
  "5" = "T_cells_1",
  "6" = "B_cells_2",
  "7" = "B_cells_3",
  "8" = "DC_1",
  "9" = "Stromal_cells",
  "10" = "Epithelial",
  "11" = "MDSC",
  "12" = "DC_2",
  "13" = "B_cells_4",
  "14" = "DC_3",
  "15" = "T_cells_2"
)

sna_timecourse <- RenameIdents(sna_timecourse, new_cluster_ids)
DimPlot(sna_timecourse, label = TRUE, repel = TRUE)
DimPlot(sna_timecourse, label = TRUE, split.by = "timepoint")
