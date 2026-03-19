
set.seed(12345)

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  this_file <- rstudioapi::getActiveDocumentContext()$path
  base_dir <- dirname(dirname(this_file)) 
} else {
  base_dir <- getwd() 
}

rds_path <- file.path(base_dir, "data_demo", "b_seurat_object.rds")
out_dir  <- file.path(base_dir, "output_example", "pseudotime_trajectory")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
message("Loading data from: ", rds_path)
message("Saving outputs to: ", out_dir)



need_cran <- c(
  "Seurat","ggplot2","ggrepel","dplyr","scales",
  "viridis","Cairo","grid"
)
for (p in need_cran) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}
invisible(lapply(need_cran, library, character.only = TRUE))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("slingshot", quietly = TRUE)) {
  BiocManager::install("slingshot", update = FALSE, ask = FALSE)
}
if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  BiocManager::install("SingleCellExperiment", update = FALSE, ask = FALSE)
}
library(slingshot)
library(SingleCellExperiment)


seu <- readRDS(rds_path)
DefaultAssay(seu) <- "RNA"

cand_ct <- c("final_type_updated","final_type","cell_type","predicted.id")
ct_col  <- cand_ct[cand_ct %in% names(seu@meta.data)][1]
if (length(ct_col) == 0) {
  stop("NO Cell：final_type_updated/final_type/cell_type/predicted.id")
}
Idents(seu) <- seu[[ct_col, drop = TRUE]]

pick_labels <- function(pats, labs){
  keep <- rep(FALSE, length(labs))
  for (p in pats) keep <- keep | grepl(p, labs, ignore.case = TRUE)
  labs[keep]
}

labs <- levels(Idents(seu))
naive_b_labs  <- pick_labels(c("^naive\\s*[_ ]?b$", "naive.?b"), labs)
memory_b_labs <- pick_labels(c("^memory\\s*[_ ]?b$", "memory.?b"), labs)
plasma_labs   <- pick_labels(c("plasmablast", "^plasma\\s*[_ ]?b$", "^plasma$"), labs)

target_labs <- unique(c(naive_b_labs, memory_b_labs, plasma_labs))
if (!length(target_labs)) {
  stop("Idents；check table(Idents(seu)).")
}

b <- subset(seu, idents = target_labs)
b$stage <- as.character(Idents(b))
b$stage[grepl(paste(naive_b_labs,  collapse="|"), b$stage, ignore.case = TRUE)] <- "Naive B"
b$stage[grepl(paste(memory_b_labs, collapse="|"), b$stage, ignore.case = TRUE)] <- "Memory B"
b$stage[grepl(paste(plasma_labs,  collapse="|"), b$stage, ignore.case = TRUE)]  <- "Plasmablast"

Idents(b) <- b$stage
b <- subset(b, idents = c("Naive B","Memory B","Plasmablast"))


md <- b@meta.data

get_ht <- function(x){
  x <- tolower(as.character(x))
  dplyr::case_when(
    grepl("(^|[^a-z0-9])(ht[-_ ]?7|hashtag[-_ ]?7|hto[-_ ]?7)([^a-z0-9]|$)", x) ~ "HT7",
    grepl("(^|[^a-z0-9])(ht[-_ ]?8|hashtag[-_ ]?8|hto[-_ ]?8)([^a-z0-9]|$)", x) ~ "HT8",
    TRUE ~ NA_character_
  )
}

ht <- rep(NA_character_, nrow(md))
for (cn in c("sample","hash.ID","HTO_maxID")) {
  if (cn %in% names(md)) {
    tmp <- get_ht(md[[cn]])
    ht[is.na(ht)] <- tmp[is.na(ht)]
  }
}
if (any(is.na(ht)) && "orig.ident" %in% names(md)) {
  x <- tolower(as.character(md$orig.ident))
  ht[is.na(ht) & grepl("ht7|hashtag.?7", x)] <- "HT7"
  ht[is.na(ht) & grepl("ht8|hashtag.?8", x)] <- "HT8"
}
b$HT <- ht
b$condition <- ifelse(b$HT=="HT7","Vaccine", ifelse(b$HT=="HT8","Control", NA))


if (!"percent.mt" %in% colnames(b@meta.data)) {
  b[["percent.mt"]] <- PercentageFeatureSet(b, pattern = "^MT-")
}
b <- subset(b, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 & percent.mt < 15)

b <- NormalizeData(b, verbose = FALSE)
b <- FindVariableFeatures(b, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
b <- ScaleData(b, features = VariableFeatures(b), verbose = FALSE)
b <- RunPCA(b, features = VariableFeatures(b), npcs = 30, verbose = FALSE)
b <- RunUMAP(
  b,
  dims = 1:30,
  reduction.name = "refUMAP",
  reduction.key  = "refUMAP_",
  verbose = FALSE
)

umap_mat <- Embeddings(b, "refUMAP")
labels   <- factor(Idents(b), levels = c("Naive B","Memory B","Plasmablast"))

sce <- SingleCellExperiment(
  assays = list(counts = matrix(0, nrow = 1, ncol = ncol(b)))
)
colnames(sce) <- colnames(b)

reducedDims(sce)$refUMAP <- umap_mat
colData(sce)$stage <- labels

sce <- slingshot(
  sce,
  clusterLabels = "stage",
  reducedDim    = "refUMAP",
  start.clus    = "Naive B",
  end.clus      = "Plasmablast"
)

pt_mat <- slingPseudotime(sce)
pt1 <- pt_mat[,1]
pt_sc <- scales::rescale(pt1, to = c(0,30), from = range(pt1, na.rm = TRUE))


curves <- slingCurves(sce)

if (length(curves) < 1) {
  stop("slingshot")
}


curve_raw <- as.data.frame(curves[[1]]$s)
colnames(curve_raw) <- c("refUMAP_1", "refUMAP_2")


curve_raw <- curve_raw[complete.cases(curve_raw), , drop = FALSE]
curve_raw <- unique(curve_raw)


step_dist <- c(
  0,
  sqrt(diff(curve_raw$refUMAP_1)^2 + diff(curve_raw$refUMAP_2)^2)
)
curve_raw$arc <- cumsum(step_dist)
curve_raw <- curve_raw[order(curve_raw$arc), ]


spx <- smooth.spline(curve_raw$arc, curve_raw$refUMAP_1, spar = 0.7)
spy <- smooth.spline(curve_raw$arc, curve_raw$refUMAP_2, spar = 0.7)

arc_new <- seq(min(curve_raw$arc), max(curve_raw$arc), length.out = 300)
curve_df <- data.frame(
  refUMAP_1 = predict(spx, arc_new)$y,
  refUMAP_2 = predict(spy, arc_new)$y
)


arrow_df <- curve_df[(nrow(curve_df)-1):nrow(curve_df), ]


plot_df <- data.frame(
  refUMAP_1  = umap_mat[,1],
  refUMAP_2  = umap_mat[,2],
  pseudotime = pt_sc,
  stage      = as.character(labels),
  condition  = b$condition,
  row.names  = colnames(b)
)

label_df <- plot_df |>
  dplyr::group_by(stage) |>
  dplyr::summarise(
    refUMAP_1 = median(refUMAP_1),
    refUMAP_2 = median(refUMAP_2),
    .groups = "drop"
  )


p_curve <- ggplot(plot_df, aes(refUMAP_1, refUMAP_2)) +
  geom_point(
    aes(color = pseudotime),
    size = 0.8, alpha = 0.75
  ) +
  scale_color_viridis(
    option = "plasma",
    direction = -1,
    limits = c(0, 30),
    name = "pseudotime"
  ) +
  geom_path(
    data = curve_df,
    aes(refUMAP_1, refUMAP_2),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1.0,
    lineend = "round",
    linejoin = "round"
  ) +
  geom_segment(
    data = arrow_df[1, , drop = FALSE],
    aes(
      x    = arrow_df$refUMAP_1[1],
      y    = arrow_df$refUMAP_2[1],
      xend = arrow_df$refUMAP_1[2],
      yend = arrow_df$refUMAP_2[2]
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1.0,
    arrow = arrow(type = "closed", length = unit(0.14, "cm"))
  ) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(refUMAP_1, refUMAP_2, label = stage),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold",
    color = c("black", "navy", "navy")[match(label_df$stage, c("Naive B","Memory B","Plasmablast"))],
    segment.alpha = 0.4,
    max.overlaps = Inf
  ) +
  labs(
    title = "B Cell Developmental Pseudotime",
    x = "refUMAP_1",
    y = "refUMAP_2"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.key.height = unit(0.6, "cm")
  )

print(p_curve)

pdf_file <- file.path(out_dir, "Bcell_pseudotime_single_curve.pdf")
svg_file <- file.path(out_dir, "Bcell_pseudotime_single_curve.svg")
png_file <- file.path(out_dir, "Bcell_pseudotime_single_curve.png")

Cairo::CairoPDF(file = pdf_file, width = 7, height = 6, onefile = TRUE)
print(p_curve)
dev.off()

ggsave(svg_file, plot = p_curve, width = 7, height = 6, units = "in", device = "svg")
ggsave(png_file, plot = p_curve, width = 7, height = 6, units = "in", dpi = 600)

message("done：\n", pdf_file, "\n", svg_file, "\n", png_file)