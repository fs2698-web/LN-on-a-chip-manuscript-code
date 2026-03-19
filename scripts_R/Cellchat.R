
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {

  script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  
  script_dir <- getwd() 
}

project_root <- dirname(script_dir)
message("location: ", project_root)

rds_path <- file.path(project_root, "data_demo", "b_seurat_object.rds")
out_dir  <- file.path(project_root, "output_example", "CellChat_Results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


if (!file.exists(rds_path)) {
  stop("File not found")
}


library(Seurat)
library(CellChat)
library(dplyr)

seurat_obj <- readRDS(rds_path)


seurat_obj$condition <- ifelse(
  seurat_obj$hash.ID == "Hashtag-7-Antibody", "Vaccine",
  ifelse(seurat_obj$hash.ID == "Hashtag-8-Antibody", "Control", NA)
)

seurat_obj$group_DT <- case_when(
  seurat_obj$cell_type %in% c("CD4 T", "CD8 T", "Treg", "Tfh") ~ "T",
  seurat_obj$cell_type == "DC" ~ "DC",
  TRUE ~ NA_character_
)

Idents(seurat_obj) <- seurat_obj$group_DT
sub_ctrl_DT <- subset(seurat_obj, 
                      idents = c("DC", "T"), 
                      subset = condition == "Control")

# ============================================================
# CellChat 
# ============================================================
data_ctrl <- GetAssayData(sub_ctrl_DT, assay = "RNA", slot = "data")
meta_ctrl <- data.frame(cell_type = Idents(sub_ctrl_DT),
                        row.names  = colnames(sub_ctrl_DT))

cc_ctrl_DT <- createCellChat(object = data_ctrl,
                             meta   = meta_ctrl,
                             group.by = "cell_type")

cc_ctrl_DT@DB <- CellChatDB.human
cc_ctrl_DT <- subsetData(cc_ctrl_DT)
cc_ctrl_DT <- identifyOverExpressedGenes(cc_ctrl_DT)
cc_ctrl_DT <- identifyOverExpressedInteractions(cc_ctrl_DT)
cc_ctrl_DT <- computeCommunProb(cc_ctrl_DT, raw.use = TRUE)
cc_ctrl_DT <- computeCommunProbPathway(cc_ctrl_DT)
cc_ctrl_DT <- aggregateNet(cc_ctrl_DT)

# ============================================================
# save
# ============================================================
plot_path <- file.path(out_dir, "Control_DC_vs_T_chord.pdf")

pdf(plot_path, width = 6, height = 6)
par(mar = c(2, 2, 2, 2))
netVisual_chord_gene(
  object      = cc_ctrl_DT,
  slot.name   = "net",
  sources.use = "DC",
  targets.use = "T",
  lab.cex     = 0.8,
  small.gap   = 2,
  big.gap     = 15,
  scale       = TRUE,
  directional = 1,
  title.name  = "Control: DC <-> T"
)
dev.off()

message("complete: ", plot_path)