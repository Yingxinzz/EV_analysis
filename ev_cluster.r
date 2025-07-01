library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(Seurat)
library(readr)
library(purrr)
library(harmony)
library(ggplot2)
library(scales)
library(jsonlite)

rm(list = ls())
gc()
file_path <- "ARZ_human/datasets"
output_matrix_path <- "ARZ_human/all_results"
output_seurat_path <- "ARZ_human/SeuratObjects"
output_merged_path <- "ARZ_human/merge"
if (!dir.exists(output_seurat_path)) {
  dir.create(output_seurat_path, recursive = TRUE)
}

file_names <- list.files(path = file_path, pattern = "*.total_ev_protein.csv", full.names = TRUE)

seurat_list <- list()

for (file in file_names) {
  message("Processing: ", basename(file))
  
  sc_data <- fread(file, data.table = FALSE)
  colnames(sc_data) <- c("symbol", "ev", "variable", "value")
  sc_data$value <- as.numeric(sc_data$value)
  
  sample_name <- unique(sc_data$symbol)
  if (length(sample_name) != 1) {
    warning("Multiple or missing sample names in file: ", basename(file))
    next
  }
  
  dat_wide <- pivot_wider(sc_data,
                          id_cols = ev,
                          names_from = variable,
                          values_from = value,
                          values_fill = 0)
  rownames(dat_wide) <- dat_wide$ev
  dat_wide <- dat_wide[, -1, drop = FALSE]
  
  valid_ev_idx <- apply(dat_wide, 1, function(x) sum(x >= 3) >= 3)
  dat_filtered <- dat_wide[valid_ev_idx, , drop = FALSE]
  n_ev <- nrow(dat_filtered)
  
  if (n_ev >= 50) {
    group_size <- 30
  } else if (n_ev >= 30) {
    group_size <- 20
  } else if (n_ev >= 10) {
    group_size <- 10
  } else {
    message("Too few EVs in ", sample_name, ", skip.")
    next
  }
  
  num_groups <- floor(n_ev / group_size)
  if (num_groups == 0) {
    message("EVs < group size in ", sample_name, ", skip.")
    next
  }
  
  set.seed(123)
  sampled_ev_idx <- sample(1:n_ev, num_groups * group_size, replace = FALSE)
  dat_sub <- dat_filtered[sampled_ev_idx, , drop = FALSE]
  groups <- rep(1:num_groups, each = group_size)
  
  pseudo_ev_list <- lapply(unique(groups), function(g) {
    rows <- which(groups == g)
    colSums(dat_sub[rows, , drop = FALSE])
  })
  names(pseudo_ev_list) <- paste0(sample_name, "_sEV", seq_along(pseudo_ev_list))
  
  pseudo_matrix <- do.call(cbind, pseudo_ev_list)
  pseudo_matrix <- as.data.frame(pseudo_matrix)
  pseudo_matrix <- pseudo_matrix[order(rownames(pseudo_matrix)), ]
  
  seu <- CreateSeuratObject(counts = pseudo_matrix,
                            project = sample_name,
                            min.features = 2,
                            min.cells = 1)
  
  seu_file <- file.path(output_seurat_path, paste0(sample_name, ".pseudo.rds"))
  write_rds(seu, seu_file)
  seurat_list[[sample_name]] <- seu
  
  message("Saved SeuratObject: ", seu_file)
  
  rm(sc_data, dat_wide, dat_filtered, dat_sub, pseudo_matrix, seu)
  gc()
}

if (length(seurat_list) >= 2) {
  message("Merging ", length(seurat_list), " Seurat objects...")
  merged_seu <- merge(seurat_list[[1]], y = seurat_list[-1], project = "MergedPseudoEVs")
  merged_seu$orig.ident <- sapply(strsplit(colnames(merged_seu), "_"), `[`, 1)
  saveRDS(merged_seu, file = file.path(output_merged_path, "MergedPseudoEVs.rds"))
  message("Merged Seurat object saved to: MergedPseudoEVs.rds")
}


merged_seu <- readRDS(file.path(output_merged_path, "MergedPseudoEVs.rds"))
merged_seu <- SCTransform(merged_seu, return.only.var.genes = FALSE, verbose = TRUE)
merged_seu <- RunPCA(merged_seu, npcs = 30, assay = "SCT", verbose = FALSE)
merged_seu <- RunHarmony(merged_seu, group.by.vars = "orig.ident", plot_convergence = TRUE)
ggsave(file.path(output_matrix_path, "Harmony_Convergence.png"), width = 15, height = 8)

merged_seu <- merged_seu %>%
  RunUMAP(reduction = "harmony", dims = 1:30, n.neighbors = 30, min.dist = 0.2, spread = 1) %>%
  RunTSNE(reduction = "harmony", dims = 1:30, perplexity = 50, learning.rate = 1000, theta = 0.5, check_duplicates = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.4)
Idents(merged_seu) <- "seurat_clusters"
saveRDS(merged_seu, file = file.path(output_matrix_path, "MergedPseudoCells_SCT_Harmony_UMAP_arzh.rds"))

merged_seu<-readRDS(file = file.path(output_matrix_path, "MergedPseudoCells_SCT_Harmony_UMAP_arzh.rds"))
min_cells <- 50
cluster_counts <- table(Idents(merged_seu))
valid_clusters <- names(cluster_counts[cluster_counts >= min_cells])
merged_seu$cluster_simplified <- ifelse(
  Idents(merged_seu) %in% valid_clusters,
  as.character(Idents(merged_seu)),
  "Other"
)
merged_seu$cluster_simplified <- factor(merged_seu$cluster_simplified)
cluster_table <- table(merged_seu$cluster_simplified)
cluster_prop <- prop.table(cluster_table) * 100
cluster_labels <- paste0(names(cluster_prop), " (", round(cluster_prop, 1), "%)")
names(cluster_labels) <- names(cluster_prop)
merged_seu$cluster_with_pct <- plyr::mapvalues(
  merged_seu$cluster_simplified,
  from = names(cluster_labels),
  to = cluster_labels
)
tsne_coords <- Embeddings(merged_seu, "tsne")
tsne_df <- data.frame(
  cell = rownames(tsne_coords),
  tSNE_1 = tsne_coords[, 1],
  tSNE_2 = tsne_coords[, 2],
  cluster = merged_seu$cluster_with_pct
)
write_json(tsne_df, path = file.path(output_matrix_path, "tsne_results.json"), pretty = TRUE)

p_batch <- DimPlot(merged_seu, reduction = "umap", group.by = "orig.ident") +
  ggtitle("UMAP by Sample (orig.ident)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18),
    legend.text = element_text(size = 12)
  )
ggsave(file.path(output_matrix_path, "UMAP_by_orig.ident.pdf"),
       p_batch, width = 10, height = 7)
print(p_batch)


plot_umap_tsne <- function(seu, group = "cluster_with_pct", output_path = ".") {
  n_cluster <- length(unique(seu[[group]][, 1]))
  palette <- hue_pal()(n_cluster)
  
  p_umap <- DimPlot(seu, reduction = "umap", group.by = group, pt.size = 0.5,
                    raster = TRUE, cols = palette) +
    ggtitle("UMAP") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "gray95", color = NA),
      plot.title = element_text(hjust = 0.5, size = 18),
      legend.text = element_text(size = 12)
    )
  ggsave(file.path(output_path, "UMAP.pdf"), p_umap, width = 12, height = 8, dpi = 300)
  print(p_umap)
  
  p_tsne <- DimPlot(seu, reduction = "tsne", group.by = group, pt.size = 0.5,
                    raster = TRUE, cols = palette) +
    ggtitle("t-SNE") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "gray95", color = NA),
      plot.title = element_text(hjust = 0.5, size = 18),
      legend.text = element_text(size = 12)
    )
  ggsave(file.path(output_path, "tSNE.pdf"), p_tsne, width = 12, height = 8, dpi = 300)
  print(p_tsne)
  
  message("UMAP and t-SNE plots saved to: ", output_path)
  return(seu)
}

merged_seu <- plot_umap_tsne(merged_seu, group = "cluster_with_pct", output_path = output_matrix_path)



