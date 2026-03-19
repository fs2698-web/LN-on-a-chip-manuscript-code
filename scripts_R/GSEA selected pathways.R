## ===================== 0. Environment Setup =====================
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  script_dir <- getwd()
}
project_root <- dirname(script_dir) 

RDS_DIR     <- file.path(project_root, "data_demo")
OBJ_PATH    <- file.path(RDS_DIR, "immune_annotated.rds")
OUTPUT_ROOT <- file.path(project_root, "output_example")
OUTDIR      <- file.path(OUTPUT_ROOT, "GSEA IL2")
FIGDIR      <- file.path(OUTDIR, "figures")

if(!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)
if(!dir.exists(FIGDIR)) dir.create(FIGDIR, recursive = TRUE)

## ===================== 1. Parameters and Packages =====================
options(stringsAsFactors = FALSE)
SPECIES  <- "Homo sapiens"
TARGET_HASHTAGS <- c("ht5","ht6","ht7","ht8")  
NEW_NAMES <- c("D3", "D7", "D9", "D10") 
FGSEA_MODE <- "multilevel" 

pkgs <- c("Seurat","presto","RcppParallel","fgsea","msigdbr",
          "data.table","dplyr","stringr","ggplot2","tidyr")

for(p in pkgs) {
  if(!require(p, character.only = TRUE)) {
    message("Installing missing package: ", p)
    if(p == "presto") {
      if(!require("devtools")) install.packages("devtools")
      devtools::install_github("immunogenomics/presto")
    } else {
      install.packages(p)
    }
    library(p, character.only = TRUE)
  }
}

ts <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

## ===================== 2. Functions =====================
detect_hashtag <- function(meta){
  cand <- grep("hash|^ht$|hashtag", names(meta), ignore.case = TRUE, value = TRUE)
  if(length(cand) == 0) stop("Error: No column containing 'hashtag' or 'hash' found!")
  cand[1]
}

map_tags <- function(avail, targets){
  sapply(targets, function(x){
    cand <- avail[stringr::str_detect(avail, stringr::regex(x, ignore_case=TRUE))]
    if(!length(cand)){
      num  <- stringr::str_extract(x, "\\d+")
      cand <- avail[stringr::str_detect(avail, stringr::regex(num, ignore_case=TRUE))]
    }
    if(!length(cand)) {
      message(paste0("Notice: Hashtag '", x, "' not found. Available: ", paste(avail, collapse=", ")))
      return(NA)
    }
    cand[1]
  }, USE.NAMES = FALSE)
}

rank_with_presto <- function(seu, group_col, tags_use){
  RcppParallel::setThreadOptions(numThreads = 1)
  DefaultAssay(seu) <- "RNA"
  seu <- NormalizeData(seu, verbose = FALSE)
  actual_tags <- tags_use[tags_use %in% unique(seu[[group_col, drop=TRUE]])]
  pre <- presto::wilcoxauc(seu, group_by = group_col, groups_use = actual_tags)
  tiny <- 1e-300
  res <- split(pre, pre$group)
  rank_list <- lapply(res, function(mk){
    mk$feature <- toupper(mk$feature)
    mk$score   <- sign(mk$logFC) * (-log10(pmax(mk$padj, tiny)))
    rnk <- mk$score; names(rnk) <- mk$feature
    sort(rnk[is.finite(rnk)], decreasing = TRUE)
  })
  list(rank = rank_list, markers = res)
}

calc_running_ES <- function(ranks, gene.set){
  ranks <- sort(ranks, decreasing = TRUE)
  N <- length(ranks); hits <- which(names(ranks) %in% gene.set)
  if (!length(hits)) return(NULL)
  Nh <- length(hits); Nm <- N - Nh
  ra <- abs(ranks); Ph <- sum(ra[hits])
  hit <- rep(FALSE, N); hit[hits] <- TRUE
  es <- cumsum(ifelse(hit, ra/Ph, -1/Nm))
  data.frame(idx = seq_len(N), es = es)
}

plot_overlay_by_hashtag <- function(curve_df, nes_tab, title, gene_set, rank_all, pad_frac = 0.08) {
  df <- curve_df
  df$tag_lab <- factor(as.character(df$tag), levels = NEW_NAMES)
  gs_up <- toupper(gene_set)
  rug_df <- lapply(unique(as.character(df$tag)), function(tr){
    rnk <- rank_all[[tr]]
    xs <- which(toupper(names(rnk)) %in% gs_up)
    if (!length(xs)) return(NULL)
    data.frame(tag_lab = tr, x = xs)
  })
  rug_df <- do.call(rbind, rug_df)
  if (!is.null(rug_df)) rug_df$tag_lab <- factor(rug_df$tag_lab, levels = NEW_NAMES)
  y_rng <- range(df$es, na.rm = TRUE); span <- diff(y_rng); if (span == 0) span <- 1e-6
  y_lim <- c(y_rng[1] - pad_frac * span, y_rng[2] + pad_frac * span)
  cols <- setNames(c("#E69F00","#56B4E9","#009E73","#CC79A7"), NEW_NAMES)
  
  ggplot(df, aes(idx, es, group = tag_lab, color = tag_lab)) +
    geom_hline(yintercept = 0, linewidth = 0.5, color = "grey35") +
    geom_line(linewidth = 1.2, lineend = "round") +
    scale_color_manual(values = cols) +
    coord_cartesian(ylim = y_lim, clip = "off") +
    labs(title = title, x = "Rank in gene list", y = "Enrichment score", color = "Sample") +
    theme_classic(base_size = 14) +
    theme(legend.position = "right", plot.title = element_text(face = "bold", size = 12),
          axis.text = element_text(color = "black"), plot.margin = margin(18, 26, 30, 30)) -> p
  
  if (!is.null(rug_df)) {
    p <- p + geom_segment(data = rug_df, aes(x = x, xend = x, y = y_lim[1], 
                                               yend = y_lim[1] + 0.04 * diff(y_lim), color = tag_lab), 
                          linewidth = 0.4, alpha = 0.6)
  }
  return(p) 
} 

do_one_pathway <- function(path_name, genes_upper, title_prefix, rank_all){
  ts(paste0("Running GSEA for: ", path_name))
  nes_tab <- data.frame(tag = names(rank_all), NES = NA_real_, FDR = NA_real_)
  for (tg in names(rank_all)) {
    stats <- rank_all[[tg]]
    res <- if (FGSEA_MODE == "fast") {
      fgsea(pathways = list(path_name = genes_upper), stats = stats, nperm = 2000, minSize = 3)
    } else {
      fgseaMultilevel(pathways = list(path_name = genes_upper), stats = stats, minSize = 3)
    }
    if (nrow(res)) nes_tab[nes_tab$tag==tg, c("NES","FDR")] <- c(round(res$NES,2), signif(res$padj,2))
  }
  
  csv_name <- file.path(OUTDIR, paste0("NES_FDR_", gsub("[^A-Za-z0-9]","_", path_name), ".csv"))
  tryCatch({ data.table::fwrite(nes_tab, csv_name) }, error = function(e){ message("Could not write CSV.") })
  
  curves <- lapply(names(rank_all), function(tg){
    cur <- calc_running_ES(rank_all[[tg]], genes_upper)
    if (is.null(cur)) return(NULL); cur$tag <- tg; cur
  }) %>% Filter(f = Negate(is.null)) %>% dplyr::bind_rows()
  
  p <- NULL
  if (nrow(curves) > 0) {
    p <- plot_overlay_by_hashtag(curves, nes_tab, title_prefix, genes_upper, rank_all)
    try({ ggsave(file.path(FIGDIR, paste0("Overlay_", gsub("[^A-Za-z0-9]","_", path_name), ".png")), p, width = 7.5, height = 5.5) })
  }
  list(nes = nes_tab, plot = p)
}

## ===================== 3. Running =====================
ts("Loading Seurat object...")
obj  <- readRDS(OBJ_PATH)

HASH_COL <- detect_hashtag(obj@meta.data)
avail_tags <- unique(as.character(obj@meta.data[[HASH_COL]]))
OLD_TAGS <- map_tags(avail_tags, TARGET_HASHTAGS)
valid_idx <- !is.na(OLD_TAGS)
OLD_TAGS <- OLD_TAGS[valid_idx]; CURRENT_NEW_NAMES <- NEW_NAMES[valid_idx]

obj@meta.data[[HASH_COL]] <- as.character(obj@meta.data[[HASH_COL]])
for(i in seq_along(OLD_TAGS)){
  obj@meta.data[[HASH_COL]][obj@meta.data[[HASH_COL]] == OLD_TAGS[i]] <- CURRENT_NEW_NAMES[i]
}
obj@meta.data[[HASH_COL]] <- factor(obj@meta.data[[HASH_COL]], levels = CURRENT_NEW_NAMES)
Idents(obj) <- obj@meta.data[[HASH_COL]]

ts("Calculating ranks...")
rank_res <- rank_with_presto(obj, HASH_COL, CURRENT_NEW_NAMES)
rank_all <- rank_res$rank

## ===================== 4. Defining IL2  =====================
ts("Defining IL2 gene set...")
hallmark <- msigdbr(species = SPECIES, collection = "H")
il2_genes <- hallmark %>% 
  filter(gs_name == "HALLMARK_IL2_STAT5_SIGNALING") %>% 
  pull(gene_symbol) %>% unique() %>% toupper()

pathways_sel <- list(IL2_STAT5_Signaling = il2_genes)

results_list <- list()
for (pn in names(pathways_sel)) {
  results_list[[pn]] <- do_one_pathway(pn, pathways_sel[[pn]], pn, rank_all)
}

## ===================== 5. Save Report =====================
while (!is.null(dev.list())) dev.off()

timestamp <- format(Sys.time(), "%H%M")
report_name <- paste0("IL2_Report_", timestamp, ".pdf")
report_full_path <- file.path(FIGDIR, report_name)

pdf(report_full_path, width = 10, height = 8)
for (pn in names(results_list)) {
  if (!is.null(results_list[[pn]]$plot)) print(results_list[[pn]]$plot)
}
dev.off()

ts(paste0("Done! IL2 analysis saved in: ", report_full_path))