## ======================= Basic configuration =======================
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  project_root <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  project_root <- getwd() 
}
RDS_DIR  <- file.path(project_root, "data_demo")
OUTDIR   <- file.path(project_root, "output_example", "GSEA")

if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

files_by_ct <- c(
  "DC"     = file.path(RDS_DIR, "DCs_annotated.rds"),
  "T cell" = file.path(RDS_DIR, "Tcells_annotated.rds"),
  "B cell" = file.path(RDS_DIR, "Bcells_annotated.rds")
)
SPECIES  <- "Homo sapiens"

## ======================= Required packages =======================
pkgs <- c("Seurat","dplyr","tibble","stringr","purrr","ggplot2",
          "data.table","msigdbr","fgsea","future")
need <- pkgs[!suppressWarnings(sapply(pkgs, requireNamespace, quietly = TRUE))]
if (length(need) > 0) install.packages(need, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

## Avoid future errors caused by object size limits
options(future.globals.maxSize = 8 * 1024^3)  # 8 GB

## ======================= Hashtag recognition rules =======================
HIGH_TAGS <- c("5","6","HT5","HT6","HTO5","HTO6","hashtag_5","hashtag_6",
               "hash5","hash6","hash_5","hash_6","Tag5","Tag6","5_6","56",
               "HTO_5","HTO_6")
LOW_TAGS  <- c("7","8","HT7","HT8","HTO7","HTO8","hashtag_7","hashtag_8",
               "hash7","hash8","hash_7","hash_8","Tag7","Tag8","7_8","78",
               "HTO_7","HTO_8")

# If you already know the column name, set it directly, e.g. "HTO_maxID";
# otherwise leave it as NULL for automatic detection
HASHTAG_COL <- NULL

## ======================= Utility functions =======================
guess_hashtag_col <- function(md, prefer = HASHTAG_COL) {
  if (!is.null(prefer) && prefer %in% colnames(md)) return(prefer)
  cand <- grep("(?:^|_)(HTO|hash|hashtag|tag)(?:$|_)", colnames(md),
               value = TRUE, ignore.case = TRUE)
  cand <- setdiff(cand, grep("nCount|nFeature|^percent", colnames(md),
                             value = TRUE, ignore.case = TRUE))
  cand <- cand[sapply(md[cand], function(x) is.character(x) || is.factor(x))]
  if (length(cand) == 0) stop("No suitable hashtag label column was found in meta.data; please set HASHTAG_COL manually.")
  message("HASHTAG_COL was not explicitly provided; automatically using column: ", cand[1])
  cand[1]
}

make_group_from_hashtag <- function(v, high_tags = HIGH_TAGS, low_tags = LOW_TAGS) {
  v <- as.character(v)
  grp <- ifelse(v %in% high_tags, "High",
                ifelse(v %in% low_tags, "Low", NA_character_))
  need_fill <- is.na(grp)
  if (any(need_fill)) {
    vv <- v[need_fill]
    grp2 <- dplyr::case_when(
      grepl("(5|6)", vv) & !grepl("(7|8)", vv) ~ "High",
      grepl("(7|8)", vv) & !grepl("(5|6)", vv) ~ "Low",
      TRUE ~ NA_character_
    )
    grp[need_fill] <- grp2
  }
  factor(grp, levels = c("Low", "High"))
}

pick_assay <- function(seu) {
  if ("RNA" %in% names(seu@assays)) return("RNA")
  if ("SCT" %in% names(seu@assays)) return("SCT")
  names(seu@assays)[1]
}

run_gsea_on_obj <- function(seu, celltype_title,
                            top_each_side = 10,
                            outdir = OUTDIR,
                            species = SPECIES) {
  stopifnot(inherits(seu, "Seurat"))
  md <- seu@meta.data

  # Force sequential mode to avoid exporting oversized objects with future
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan("sequential")

  hashtag_col <- guess_hashtag_col(md, HASHTAG_COL)
  hv <- md[[hashtag_col]]
  if (is.null(hv)) stop("Hashtag column does not exist: ", hashtag_col)

  seu$GSEA_group <- make_group_from_hashtag(hv)
  seu <- subset(seu, !is.na(GSEA_group))
  if (ncol(seu) < 50) stop("Too few valid cells (<50); please check whether the hashtag grouping is correct.")

  Idents(seu) <- "GSEA_group"
  DefaultAssay(seu) <- pick_assay(seu)

  markers <- Seurat::FindMarkers(
    seu,
    ident.1 = "High", ident.2 = "Low",
    logfc.threshold = 0, min.pct = 0, min.diff.pct = -Inf,
    test.use = "wilcox", verbose = FALSE
  ) %>% tibble::rownames_to_column("gene") %>% dplyr::distinct(gene, .keep_all = TRUE)

  ranks <- markers$avg_log2FC
  names(ranks) <- markers$gene
  ranks <- sort(ranks[!is.na(ranks)], decreasing = TRUE)
  if (length(ranks) < 50) stop("Too few ranked genes available for GSEA.")

  hallmark <- msigdbr::msigdbr(species = species, category = "H") %>%
    dplyr::select(gs_name, gene_symbol) %>% dplyr::distinct()
  pathways <- split(hallmark$gene_symbol, hallmark$gs_name)

  set.seed(1)
  fg <- fgsea::fgseaMultilevel(pathways, ranks, minSize = 10, maxSize = 500)
  fg <- tibble::as_tibble(fg) %>% dplyr::arrange(padj)

  fg_out <- fg %>%
    dplyr::mutate(leadingEdge = vapply(leadingEdge, function(x) paste(x, collapse = ";"), character(1)))

  saveRDS(fg, file.path(outdir, paste0(gsub(" ", "_", celltype_title), "_Hallmark_fgsea.rds")))
  write.csv(
    fg_out,
    file.path(outdir, paste0(gsub(" ", "_", celltype_title), "_Hallmark_fgsea.csv")),
    row.names = FALSE
  )

  top_pos <- fg %>% dplyr::filter(NES > 0) %>% dplyr::slice_min(order_by = padj, n = top_each_side)
  top_neg <- fg %>% dplyr::filter(NES < 0) %>% dplyr::slice_min(order_by = padj, n = top_each_side)

  top_df <- dplyr::bind_rows(top_neg, top_pos) %>%
    dplyr::mutate(
      gs_name   = gsub("^HALLMARK_", "", .data$pathway),
      gs_name   = gsub("_", " ", gs_name),  # Replace underscores with spaces for better readability
      Direction = factor(ifelse(NES > 0, "High", "Low"),
                         levels = c("High", "Low")),
      # Use -log10 significance so that bubble size is more intuitive
      neg_log10_padj = -log10(pmax(padj, 1e-10))
    ) %>% dplyr::arrange(NES)

  top_df <- top_df %>%
    mutate(
      # pmax 
      logP = -log10(pmax(padj, 1e-300)) 
    )
op_df <- top_df %>% 
mutate(logP = ifelse(logP < -log10(0.05), -log10(0.05), logP))
actual_breaks_count <- sum(c(-log10(0.05), -log10(0.02), -log10(0.01), -log10(0.001)) <= max(top_df$logP))
my_override_size <- switch(as.character(actual_breaks_count),
                           "4" = c(2, 3, 5, 7),  
                           "3" = c(2, 3, 5))
p <- ggplot(top_df, aes(x = NES, y = reorder(gs_name, NES))) +
    geom_point(aes(size = logP, fill = Direction), shape = 21, color = "black", stroke = 0.3) +
    scale_size_continuous(
      name = "FDR",
      # range 
      range = c(2, 8), 
      breaks = c(-log10(0.05), -log10(0.02), -log10(0.01),-log10(0.001)), 
      labels = c("0.05","0.02", "0.01", "0.001")
    ) +
   guides(size = guide_legend(override.aes = list(size = my_override_size))) +
    scale_fill_manual(values = c("High" = "#D62728", "Low" = "#1F77B4")) +
    labs(
      title = paste0("GSEA: ", celltype_title),
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      subtitle = "Size: Significance | Color: Enrichment Direction"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black"),
      legend.key.height = unit(1, "lines")
    )

  fn_png <- file.path(outdir, paste0("GSEA_", gsub(" ", "_", celltype_title), "_Hallmark_dot.png"))
  ggsave(fn_png, plot = p, width = 5.2, height = 4.2, dpi = 300)

  message("Finished: ", celltype_title, " | Plot saved to: ", fn_png)
  list(fgsea = fg, plot = p, top_table = top_df)
}

## ======================= Main workflow =======================
results_list <- list()
for (ct in names(files_by_ct)) {
  f <- files_by_ct[[ct]]
  if (!file.exists(f)) stop("File not found: ", f, "\nPlease check the path: ", RDS_DIR)
  seu <- readRDS(f)
  results_list[[ct]] <- run_gsea_on_obj(seu, celltype_title = ct)
}

pdf(file.path(OUTDIR, "GSEA_DC_Tcell_Bcell_Hallmark.pdf"), width = 7.5, height = 10)
print(results_list[["DC"]]$plot)
print(results_list[["T cell"]]$plot)
print(results_list[["B cell"]]$plot)
dev.off()