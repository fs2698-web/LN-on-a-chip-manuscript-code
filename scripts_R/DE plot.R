
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  script_dir <- getwd()
}
project_root <- dirname(script_dir) 

library(Seurat)
library(ggplot2)
library(scales)
library(dplyr)

rds_path <- file.path(project_root, "data_demo", "seurat_obj.rds")
outdir   <- file.path(project_root, "output_example", "DE plot")

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)


seurat_obj <- readRDS(rds_path)

# ==========================================
# 2. Metadata 
# ==========================================
seurat_obj$condition4 <- dplyr::case_when(
    seurat_obj$HTO_maxID == "Hashtag-5-Antibody" ~ "LN",
    seurat_obj$HTO_maxID == "Hashtag-6-Antibody" ~ "PBMC",
    seurat_obj$HTO_maxID == "Hashtag-8-Antibody" ~ "Control",
    seurat_obj$HTO_maxID == "Hashtag-7-Antibody" ~ "Vaccine",
    TRUE ~ NA_character_
)
seurat_obj$condition4 <- factor(
    seurat_obj$condition4,
    levels = c("LN", "PBMC", "Control", "Vaccine")
)

# ==========================================
# 3. DATA
# ==========================================
genes_order <- c(
    "MS4A1", "CD19", "BCL6", "CD27", "CD38", "MZB1", 
    "CD3D", "CD4", "CD69", "ICOS", "ITGAX", "HLA-DRA", 
    "CXCL13", "PDPN", "PECAM1", "MKI67", "BCL2"
)

genes_present <- genes_order[genes_order %in% rownames(seurat_obj)]

p_temp <- DotPlot(
    object = seurat_obj,
    features = genes_present,
    group.by = "condition4"
)

plot_df2 <- p_temp$data %>%
    filter(id %in% c("Control", "Vaccine")) %>%
    mutate(
        id = factor(id, levels = c("Control", "Vaccine")),
        features.plot = factor(features.plot, levels = rev(genes_present)) # rev 保证从上往下排序
    )

# ==========================================
# 4. ggplot2 
# ==========================================
p_two_from_four <- ggplot(plot_df2, aes(x = id, y = features.plot)) +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
    scale_color_gradientn(
        colors = c("#f2eef7", "#dadaeb", "#bcbddc", "#807dba", "#3f51a2"),
        limits = c(-1, 1),
        oob = squish,
        name = "Average\nExpression"
    ) +
    scale_size(
        range = c(0.5, 8),
        breaks = c(20, 40, 60, 80),
        name = "Percent\nExpressed"
    ) +
    theme_classic(base_size = 16) +
    theme(
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        plot.margin = margin(10, 10, 10, 10)
    )

print(p_two_from_four)

# ==========================================
# 5. SAVE
# ==========================================
file_base <- file.path(outdir, "Control_vs_Vaccine_DE_plot")

# PNG
ggsave(
    filename = paste0(file_base, ".png"),
    plot = p_two_from_four,
    width = 5,
    height = 7,
    dpi = 300,
    bg = "white"
)

# PDF
ggsave(
    filename = paste0(file_base, ".pdf"),
    plot = p_two_from_four,
    width = 5,
    height = 7,
    device = cairo_pdf, 
    bg = "white"
)

message("Done! Results saved in: ", outdir)
