# ============================================================
# Violin Plots (Tfh & AP-1)
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(Cairo)
})
# ---------------------------
# 0) User settings
# ---------------------------

if (!requireNamespace("rstudioapi", quietly = TRUE)) install.packages("rstudioapi")


script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
project_root <- dirname(script_dir) 


rds_path <- file.path(project_root, "data_demo", "Tcells_annotated.rds")
out_dir  <- file.path(project_root, "output_example", "AP1 and Tfh score")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---------------------------
# 1) Load object
# ---------------------------
obj <- readRDS(rds_path)
meta_cols <- colnames(obj[[]])
donor_col_candidates <- c("hash.ID", "hashtag", "donor")
donor_col <- donor_col_candidates[donor_col_candidates %in% meta_cols][1]
current_vals <- as.character(obj[[donor_col, drop = TRUE]])
obj$donor_new <- case_when(
  grepl("5", current_vals) ~ "D3",
  grepl("8", current_vals) ~ "D7",
  grepl("6", current_vals) ~ "D9",
  grepl("7", current_vals) ~ "D10",
  TRUE ~ NA_character_
)
obj <- subset(obj, subset = donor_new %in% c("D3", "D7", "D9", "D10"))
obj$donor_new <- factor(obj$donor_new, levels = c("D3", "D7", "D9", "D10"))
# color
my_cols <- c("D3" = "#E69F00", "D7" = "#F04444", "D9" = "#56B4E9", "D10" = "#009E73")

# ---------------------------
# 2) Define gene sets & module scores
# ---------------------------
tfh_genes <- c("BCL6","CXCR5","ICOS","PDCD1","IL21")
ap1_genes <- c("FOS","FOSB","JUN","JUNB","JUND","ATF3")

tfh_present <- tfh_genes[tfh_genes %in% rownames(obj)]
ap1_present <- ap1_genes[ap1_genes %in% rownames(obj)]

obj <- AddModuleScore(obj, features = list(tfh_present), name = "TfhScore")
obj <- AddModuleScore(obj, features = list(ap1_present), name = "AP1Score")

# ---------------------------
# 3) Plotting
# ---------------------------

# Tfh Score
pC1 <- VlnPlot(obj, features = "TfhScore1", group.by = "donor_new", pt.size = 0, cols = my_cols) +
  theme_classic() +
  NoLegend() +
  labs(title = "Tfh score", x = "")

# AP1 Score
pD1 <- VlnPlot(obj, features = "AP1Score1", group.by = "donor_new", pt.size = 0, cols = my_cols) +
  theme_classic() +
  NoLegend() +
  labs(title = "AP-1 score", x = "")

# save
ggsave(file.path(out_dir, "Tfh_Violin.pdf"), pC1, width = 5, height = 5, device = cairo_pdf)
ggsave(file.path(out_dir, "AP1_Violin.pdf"), pD1, width = 5, height = 5, device = cairo_pdf)
ggsave(file.path(out_dir, "Tfh_Violin.png"), pC1, width = 5, height = 5, dpi = 300)
ggsave(file.path(out_dir, "AP1_Violin.png"), pD1, width = 5, height = 5, dpi = 300)

# ---------------------------
# 4)Donor-level correlation
# ---------------------------
df_scores <- FetchData(obj, vars = c("donor_new", "TfhScore1", "AP1Score1"))
colnames(df_scores) <- c("donor_new", "TfhScore1", "AP1Score1")

donor_summary <- df_scores %>%
  group_by(donor_new) %>%
  summarise(
    mean_Tfh = mean(TfhScore1, na.rm = TRUE),
    mean_AP1 = mean(AP1Score1, na.rm = TRUE),
    .groups = "drop"
  )
if(nrow(donor_summary) > 1) {
  cor_val <- cor(donor_summary$mean_AP1, donor_summary$mean_Tfh)
  subtitle_text <- paste0("Pearson r = ", round(cor_val, 3))
} else {
  subtitle_text <- ""
}

pE <- ggplot(donor_summary, aes(x = mean_AP1, y = mean_Tfh, label = donor_new)) +
  geom_smooth(method = "lm", se = FALSE, color = "grey", linetype = "dashed") +
  geom_point(aes(color = donor_new), size = 6) +
  scale_color_manual(values = my_cols) +
  geom_text(nudge_y = 0.005) +
  theme_classic() +
  labs(title = "Donor Correlation", x = "Mean AP-1 Score", y = "Mean Tfh Score")

ggsave(file.path(out_dir, "Correlation.pdf"), pE, width = 6, height = 5, device = cairo_pdf)

