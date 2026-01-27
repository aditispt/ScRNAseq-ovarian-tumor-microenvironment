##Comparing Fibroblasts to CAFs day 25
caf_modules <- list(
  iCAF = c("Cxcl1", "Cxcl2", "Cxcl12", "Has1", "Has2", "Il6"),
  myCAF = c("Acta2", "Myl9", "Tagln", "Postn"),
  apCAF = c("H2-Aa", "H2-Ab1", "Cd74", "Ccl21"),
  ecmCAF = c("Postn", "Cthrc1", "Col6a3", "Fn1", "Mmp14", "Mmp11")
)

# Keep only genes present in Day 25 dataset
caf_modules <- lapply(
  caf_modules,
  function(g) g[g %in% rownames(sna)]
)

fibro_d25 <- subset(
  sna,
  idents = "Fibroblasts"
)

fibro_d25_list <- SplitObject(
  fibro_d25,
  split.by = "SNA_level_D25"
)

names(fibro_d25_list)
# should be "SNAhigh" and "SNAlow"

calc_caf_presence_matrix <- function(seurat_obj, caf_modules, pct_thresh = 10) {
  
  DefaultAssay(seurat_obj) <- "RNA"
  
  clusters <- levels(Idents(seurat_obj))
  
  mat <- matrix(
    0,
    nrow = length(clusters),
    ncol = length(caf_modules)
  )
  
  rownames(mat) <- paste0("Fibroblast_", clusters)
  colnames(mat) <- names(caf_modules)
  
  for (cl in clusters) {
    cl_cells <- WhichCells(seurat_obj, idents = cl)
    
    expr <- FetchData(
      seurat_obj,
      vars = unique(unlist(caf_modules)),
      cells = cl_cells,
      layer = "counts"
    )
    
    for (caf in names(caf_modules)) {
      genes <- caf_modules[[caf]]
      
      present_genes <- sapply(genes, function(g) {
        mean(expr[[g]] > 0) * 100 >= pct_thresh
      })
      
      mat[paste0("Fibroblast_", cl), caf] <-
        sum(present_genes) / length(genes) * 100
    }
  }
  
  as.data.frame(mat)
}

##Compute CAF overload
caf_presence_d25_high <- calc_caf_presence_matrix(
  fibro_d25_list[["SNAhigh"]],
  caf_modules
)

caf_presence_d25_low <- calc_caf_presence_matrix(
  fibro_d25_list[["SNAlow"]],
  caf_modules
)

#Prepare Heatmap
prep_heatmap_df <- function(df) {
  df %>%
    rownames_to_column("Cluster") %>%
    pivot_longer(
      cols = -Cluster,
      names_to = "CAF_State",
      values_to = "Percent"
    ) %>%
    mutate(
      Cluster = factor(Cluster, levels = rev(unique(Cluster))),
      CAF_State = factor(
        CAF_State,
        levels = c("iCAF", "myCAF", "apCAF", "ecmCAF")
      )
    )
}

hm_d25_high <- prep_heatmap_df(caf_presence_d25_high)
hm_d25_low  <- prep_heatmap_df(caf_presence_d25_low)

plot_caf_heatmap <- function(df, title) {
  ggplot(df, aes(x = CAF_State, y = Cluster, fill = Percent)) +
    geom_tile(color = "grey70", size = 0.6) +
    geom_text(aes(label = sprintf("%.1f", Percent)), size = 4) +
    scale_fill_gradient(
      low = "white",
      high = "limegreen",
      limits = c(0, 100),
      name = "%"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    ggtitle(title)
}

plot_caf_heatmap(
  hm_d25_high,
  "CAF Marker Presence in Fibroblasts – Day 25 SNAhigh"
)

plot_caf_heatmap(
  hm_d25_low,
  "CAF Marker Presence in Fibroblasts – Day 25 SNAlow"
)

##Assign module scores
fibro_d25 <- subset(
  sna,
  idents = "Fibroblasts"
)

fibro_d25 <- AddModuleScore(
  fibro_d25,
  features = caf_modules,
  name = "CAF_Score",
  nbin = 10
)

fibro_d25@meta.data <- fibro_d25@meta.data %>%
  dplyr::rename(
    iCAF_score   = CAF_Score1,
    myCAF_score  = CAF_Score2,
    apCAF_score  = CAF_Score3,
    ecmCAF_score = CAF_Score4
  )

VlnPlot(
  fibro_d25,
  features = c("iCAF_score", "myCAF_score", "apCAF_score", "ecmCAF_score"),
  group.by = "SNA_level_D25",
  pt.size = 0,
  ncol = 2
)

library(ggplot2)

fibro_d25@meta.data %>%
  dplyr::select(
    SNA_level_D25,
    iCAF_score,
    myCAF_score,
    apCAF_score,
    ecmCAF_score
  ) %>%
  tidyr::pivot_longer(
    -SNA_level_D25,
    names_to = "CAF_Type",
    values_to = "Score"
  ) %>%
  ggplot(aes(x = SNA_level_D25, y = Score, fill = SNA_level_D25)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ CAF_Type, scales = "free_y") +
  theme_minimal(base_size = 14)

