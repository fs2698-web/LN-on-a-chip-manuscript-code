
library(Seurat)
library(ggplot2)
library(patchwork)
library(scales)
library(cowplot) 

options(stringsAsFactors = FALSE)
set.seed(42)

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  this_file <- rstudioapi::getActiveDocumentContext()$path
  base_dir <- dirname(dirname(this_file))
} else {
  base_dir <- getwd()
}

base_path <- file.path(base_dir, "data_demo")
seurat_obj <- readRDS(file.path(base_path, "seurat_obj.rds"))

output_dir <- file.path(base_dir, "output_example", "UMAP_outputs")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 拆分对象
seurat_h5 <- subset(seurat_obj, subset = hash.ID == "Hashtag-5-Antibody")
seurat_h6 <- subset(seurat_obj, subset = hash.ID == "Hashtag-6-Antibody")
seurat_h7 <- subset(seurat_obj, subset = hash.ID == "Hashtag-7-Antibody")

rename_cells <- function(obj) {
    obj$mapped_type <- obj$final_type
    obj$mapped_type[obj$mapped_type == "Endothelial_cells"] <- "ECs"
    return(obj)
}

seurat_h5 <- rename_cells(seurat_h5)
seurat_h6 <- rename_cells(seurat_h6)
seurat_h7 <- rename_cells(seurat_h7)

DefaultAssay(seurat_h7) <- "RNA"
seurat_h7 <- NormalizeData(seurat_h7, verbose = FALSE)
seurat_h7 <- FindVariableFeatures(seurat_h7, nfeatures = 2000, verbose = FALSE)
seurat_h7 <- ScaleData(seurat_h7, verbose = FALSE)
seurat_h7 <- RunPCA(seurat_h7, npcs = 30, verbose = FALSE)
seurat_h7 <- RunUMAP(seurat_h7, reduction = "pca", dims = 1:30, 
                     reduction.name = "ref.umap", return.model = TRUE, verbose = FALSE)

map_to_h7 <- function(query, ref) {
    DefaultAssay(query) <- "RNA"
    query <- NormalizeData(query, verbose = FALSE)
    query <- FindVariableFeatures(query, verbose = FALSE)
    anchors <- FindTransferAnchors(reference = ref, query = query, reference.reduction = "pca", dims = 1:30)
    MapQuery(anchorset = anchors, query = query, reference = ref, 
             reference.reduction = "pca", reduction.model = "ref.umap")
}

seurat_h5_mapped <- map_to_h7(seurat_h5, seurat_h7)
seurat_h6_mapped <- map_to_h7(seurat_h6, seurat_h7)

message("Processing cell types and filtering rare populations...")

target_types <- c("Fibroblasts", "ECs", "DC", "Monocyte", "CD4 T", 
                  "Treg", "CD8 T", "Tfh", "Tfr", "NK", 
                  "Naive B", "Memory B", "Plasmablast")

custom_colors <- c(
  "Fibroblasts" = "#E15E8A", "ECs" = "#9E76A3", "DC" = "#00A4D1",
  "Monocyte" = "#4397D3", "CD4 T" = "#D87000", "Treg" = "#7B9625",
  "CD8 T" = "#E67F6B", "Tfh" = "#A88714", "Tfr" = "#8B8B21",
  "NK" = "#C86E9E", "Naive B" = "#20A352", "Memory B" = "#00B08B",
  "Plasmablast" = "#14ACB3", "Other" = "#E15E8C"
)

prepare_for_fig <- function(obj, is_pbmc = FALSE) {
  if ("predicted.id" %in% colnames(obj@meta.data)) {
    obj$mapped_type <- as.character(obj$predicted.id)
  } else {
    obj$mapped_type <- as.character(obj$final_type_updated)
  }
  obj$mapped_type[obj$mapped_type == "Endothelial_cells"] <- "ECs"
  stem_types <- c("CMP", "GMP", "HSC_CD34+", "MEP", "MSC", "Tissue_stem_cells")
  obj$mapped_type[obj$mapped_type %in% stem_types] <- "Other"
  

  counts <- table(obj$mapped_type)
  rare_types <- names(counts[counts < 11])
  obj$mapped_type[obj$mapped_type %in% rare_types | !(obj$mapped_type %in% target_types)] <- "Other"
  

  if (is_pbmc) {
    obj <- subset(obj, subset = mapped_type != "Other")
  }
  

  obj$mapped_type <- factor(obj$mapped_type, levels = c(target_types, "Other"))
  return(obj)
}

h5_final <- prepare_for_fig(seurat_h5_mapped)
h6_final <- prepare_for_fig(seurat_h6_mapped, is_pbmc = TRUE)
h7_final <- prepare_for_fig(seurat_h7)

message("Generating final plots...")

plot_custom <- function(obj, title) {
  DimPlot(obj, reduction = "ref.umap", group.by = "mapped_type", 
          label = FALSE, pt.size = 0.5) +
    scale_color_manual(values = custom_colors, breaks = target_types) +
    ggtitle(title) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      legend.position = "none", 
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

p1 <- plot_custom(h5_final, "Lymph node mononuclear cells")
p2 <- plot_custom(h6_final, "PBMC")
p3 <- plot_custom(h7_final, "On-chip lymph node niche")

dummy_data <- data.frame(x = 1, y = 1, type = factor(target_types, levels = target_types))
shared_legend_plot <- ggplot(dummy_data, aes(x, y, color = type)) +
  geom_point() +
  scale_color_manual(values = custom_colors, name = "Cell Types") +
  theme_void() +
  guides(color = guide_legend(override.aes = list(size = 4)))

unique_legend <- cowplot::get_legend(shared_legend_plot)

final_combined_plot <- (p1 | p2 | p3 | wrap_plots(unique_legend)) + 
  plot_layout(widths = c(1, 1, 1, 0.2))

print(final_combined_plot)

output_path <- file.path(output_dir, "Combined_UMAP.png")
ggsave(output_path, final_combined_plot, width = 20, height = 5, dpi = 300)

message("Success! The combined UMAP has been saved.")