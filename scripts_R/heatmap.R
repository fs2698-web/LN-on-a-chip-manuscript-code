## =========================================================
## Immune functional heatmap across target hashtags
## Fixed version for nested pb_tag structure
## Hashtag -> Donor mapping:
## HT5 = D3, HT6 = D9, HT7 = D10, HT8 = D7
## =========================================================

## ===================== 0. Packages =====================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(msigdbr)
  library(edgeR)
  library(pheatmap)
})

## ===================== 1. Auto path =====================
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  this_file <- rstudioapi::getActiveDocumentContext()$path
  base_dir <- dirname(dirname(this_file))
} else {
  base_dir <- getwd()
}

## ===================== 2. Helper functions =====================

drop_sex_genes <- function(genes) {
  sex_genes <- c(
    "XIST", "TSIX", "SRY", "DDX3Y", "USP9Y", "UTY",
    "EIF1AY", "KDM5D", "RPS4Y1", "ZFY"
  )
  setdiff(genes, sex_genes)
}

cap_z <- function(x, limit = 3) {
  x[x > limit] <- limit
  x[x < -limit] <- -limit
  x
}

safe_file <- function(...) {
  normalizePath(file.path(...), winslash = "/", mustWork = FALSE)
}

## =========================================================
## Recursively find count-like matrices inside nested pb_tag
## =========================================================
find_count_mats <- function(x, prefix = NULL) {
  mats <- list()

  if (is.matrix(x) || is.data.frame(x)) {
    m <- as.matrix(x)
    if (!is.null(rownames(m)) && !is.null(colnames(m))) {
      nm <- ifelse(is.null(prefix), "matrix", prefix)
      mats[[nm]] <- m
    }
    return(mats)
  }

  if (is.list(x)) {
    nms <- names(x)
    if (is.null(nms)) nms <- seq_along(x)

    for (i in seq_along(x)) {
      nm <- as.character(nms[i])
      sub_prefix <- if (is.null(prefix)) nm else paste(prefix, nm, sep = "/")
      mats <- c(mats, find_count_mats(x[[i]], sub_prefix))
    }
  }

  mats
}

## =========================================================
## Build global counts: gene x target_hashtags
## Works for nested pb_tag:
## pb_tag[[B/T/DC/Mac]][[subcluster]][[...]]
## =========================================================
build_global_counts <- function(pb, tags) {

  all_mats <- find_count_mats(pb)

  if (length(all_mats) == 0) {
    stop("No count-like matrices found inside pb_tag.")
  }

  cat("Found ", length(all_mats), " matrix/data.frame objects inside pb_tag\n", sep = "")
  cat("Example matrix paths:\n")
  print(utils::head(names(all_mats), 10))

  usable_mats <- lapply(all_mats, function(m) {
    valid_cols <- intersect(tags, colnames(m))
    if (length(valid_cols) == 0) return(NULL)
    m[, valid_cols, drop = FALSE]
  })
  usable_mats <- usable_mats[!vapply(usable_mats, is.null, logical(1))]

  if (length(usable_mats) == 0) {
    stop(
      "No matrices inside pb_tag contain target hashtags.\n",
      "target_hashtags = ", paste(tags, collapse = ", "), "\n",
      "Please inspect internal column names."
    )
  }

  cat("Usable matrices containing target hashtags: ", length(usable_mats), "\n", sep = "")

  all_genes <- unique(unlist(lapply(usable_mats, rownames)))
  all_genes <- all_genes[!is.na(all_genes) & nzchar(all_genes)]

  counts_global <- matrix(
    0,
    nrow = length(all_genes),
    ncol = length(tags),
    dimnames = list(all_genes, tags)
  )

  for (m in usable_mats) {
    genes <- rownames(m)
    cols  <- colnames(m)

    genes <- genes[!is.na(genes) & nzchar(genes)]
    m <- m[genes, , drop = FALSE]

    counts_global[genes, cols] <- counts_global[genes, cols] + m
  }

  keep <- rowSums(counts_global, na.rm = TRUE) > 0
  counts_global <- counts_global[keep, , drop = FALSE]

  cat("counts_global dim: ", nrow(counts_global), " x ", ncol(counts_global), "\n", sep = "")
  cat("counts_global colnames:\n")
  print(colnames(counts_global))
  cat("counts_global rownames head:\n")
  print(utils::head(rownames(counts_global)))

  return(counts_global)
}

## =========================================================
## Optional inspection helper
## =========================================================
inspect_pb_cols <- function(x, prefix = NULL) {
  if (is.matrix(x) || is.data.frame(x)) {
    cat("\n---", ifelse(is.null(prefix), "matrix", prefix), "---\n")
    cat("dim:", paste(dim(x), collapse = " x "), "\n")
    cat("colnames:\n")
    print(colnames(x))
    return(invisible(NULL))
  }

  if (is.list(x)) {
    nms <- names(x)
    if (is.null(nms)) nms <- seq_along(x)

    for (i in seq_along(x)) {
      nm <- as.character(nms[i])
      sub_prefix <- if (is.null(prefix)) nm else paste(prefix, nm, sep = "/")
      inspect_pb_cols(x[[i]], sub_prefix)
    }
  }
}

## =========================================================
## Main plotting function
## =========================================================
immune_function_heatmap <- function(
  pb_tag,
  deg_pairs,
  target_hashtags,
  out_pdf,
  donor_map = NULL,
  top_per_bucket = 8,
  drop_mt_ribo = TRUE
) {
  ## ---------- Hallmark sets ----------
  h <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H") |> as.data.frame()
  hset <- function(name) unique(h$gene_symbol[h$gs_name == name])

  bucket_list <- list(
    IFN_alpha_gamma = unique(c(
      hset("HALLMARK_INTERFERON_ALPHA_RESPONSE"),
      hset("HALLMARK_INTERFERON_GAMMA_RESPONSE")
    )),
    TNFA_NFKB       = hset("HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
    IL6_JAK_STAT3   = hset("HALLMARK_IL6_JAK_STAT3_SIGNALING"),
    Inflammatory    = hset("HALLMARK_INFLAMMATORY_RESPONSE"),
    Antigen_Present = unique(c(
      "HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","B2M","TAP1","TAP2","CIITA"
    )),
    DC_Maturation   = c(
      "CCR7","LAMP3","CD83","IL12B","CCL19","CCL21","HLA-DRA","HLA-DRB1"
    ),
    Tfh_Program     = c(
      "CXCR5","PDCD1","ICOS","BCL6","IL21","SH2D1A"
    ),
    GC_PB_Program   = c(
      "AICDA","PRDM1","XBP1","MZB1","JCHAIN","IRF4"
    ),
    Exhaustion_Inhib = c(
      "PDCD1","CTLA4","LAG3","HAVCR2","TIGIT","ENTPD1","TOX"
    )
  )

  immune_universe <- unique(unlist(bucket_list))

  ## ---------- DEG subset ----------
  if (!"contrast" %in% names(deg_pairs)) {
    stop("deg_pairs must contain a column named 'contrast'.")
  }
  if (!"gene" %in% names(deg_pairs)) {
    stop("deg_pairs must contain a column named 'gene'.")
  }
  if (!"logFC" %in% names(deg_pairs)) {
    stop("deg_pairs must contain a column named 'logFC'.")
  }

  deg_sub <- deg_pairs[
    sapply(contrast, function(x) {
      parts <- strsplit(x, "_vs_")[[1]]
      length(parts) == 2 && all(parts %in% target_hashtags)
    })
  ]

  if (!nrow(deg_sub)) {
    stop(
      "No target-only contrasts for immune buckets.\n",
      "Check contrast names in DEG file and target_hashtags."
    )
  }

  ## ---------- Gene ranking ----------
  fdr_candidates <- intersect(c("FDR","adj.P.Val","FDR.x","FDR.y"), names(deg_sub))

  if (length(fdr_candidates) > 0) {
    fdr_col <- fdr_candidates[1]
    rank_tbl <- deg_sub[, .(
      minFDR = suppressWarnings(min(get(fdr_col), na.rm = TRUE)),
      maxAbsLogFC = max(abs(logFC), na.rm = TRUE)
    ), by = gene]
  } else {
    rank_tbl <- deg_sub[, .(
      minFDR = 1,
      maxAbsLogFC = mean(abs(logFC), na.rm = TRUE)
    ), by = gene]
  }

  rank_tbl <- rank_tbl[is.finite(maxAbsLogFC)]
  rank_tbl <- rank_tbl[!is.na(gene) & nzchar(gene)]

  cand <- intersect(immune_universe, rank_tbl$gene)
  cand <- drop_sex_genes(cand)

  if (drop_mt_ribo) {
    cand <- cand[!grepl("^MT-|^RPL|^RPS", cand, ignore.case = TRUE)]
  }

  rank_tbl <- rank_tbl[gene %in% cand]

  ## ---------- Per-bucket top genes ----------
  pick <- character(0)

  for (bk in names(bucket_list)) {
    genes_bk <- intersect(bucket_list[[bk]], rank_tbl$gene)
    if (!length(genes_bk)) next

    sub_rank <- rank_tbl[gene %in% genes_bk][order(minFDR, -maxAbsLogFC)]
    pick <- c(pick, head(sub_rank$gene, top_per_bucket))
  }

  target_total <- 9 * top_per_bucket
  if (length(pick) < min(target_total, nrow(rank_tbl))) {
    rest <- setdiff(rank_tbl[order(minFDR, -maxAbsLogFC)]$gene, pick)
    pick <- unique(c(pick, head(rest, target_total - length(pick))))
  }

  pick <- unique(pick)

  if (length(pick) < 2) {
    stop("Too few immune genes selected before global count matching.")
  }

  ## ---------- Build global counts ----------
  counts_global <- build_global_counts(pb_tag, target_hashtags)

  cat("Selected genes in DEG bucket: ", length(pick), "\n", sep = "")
  genes_use <- intersect(pick, rownames(counts_global))
  cat("Genes matched in counts_global: ", length(genes_use), "\n", sep = "")

  if (length(genes_use) < 2) {
    missing_show <- head(setdiff(pick, rownames(counts_global)), 20)
    stop(
      "Too few genes remain in global counts.\n",
      "Matched genes: ", length(genes_use), "\n",
      "Example missing genes: ", paste(missing_show, collapse = ", "), "\n",
      "This suggests pb_tag rownames are not gene symbols or the internal matrices use a different naming scheme."
    )
  }

  ## ---------- CPM + z-score ----------
  counts_use <- counts_global[genes_use, target_hashtags, drop = FALSE]
  storage.mode(counts_use) <- "numeric"

  logcpm <- edgeR::cpm(counts_use, log = TRUE, prior.count = 1)

  zmat <- t(scale(t(logcpm)))
  zmat[!is.finite(zmat)] <- 0
  zmat <- cap_z(zmat, limit = 3)

  ## ---------- Rename columns to donor names ----------
  original_hashtags <- colnames(zmat)

  if (!is.null(donor_map)) {
    donor_labels <- donor_map[original_hashtags]
    donor_labels[is.na(donor_labels)] <- original_hashtags
    colnames(zmat) <- donor_labels
  }

  ## ---------- Annotation ----------
  ann_col <- data.frame(
    Hashtag = original_hashtags,
    Donor = colnames(zmat),
    row.names = colnames(zmat)
  )

  ## ---------- Plot ----------
  if (!dir.exists(dirname(out_pdf))) {
    dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  }

  pheatmap::pheatmap(
    zmat,
    filename = out_pdf,
    width = 8.8,
    height = max(6, min(30, nrow(zmat) * 0.12)),
    color = colorRampPalette(c("#2c7bb6", "#ffffbf", "#d7191c"))(100),
    scale = "none",
    clustering_distance_rows = "correlation",
    clustering_method = "complete",
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    annotation_col = ann_col,
    fontsize = 9,
    fontsize_row = 6.5,
    fontsize_col = 10,
    main = paste0(
      "Immune functional heatmap across donors (",
      paste(colnames(zmat), collapse = ", "),
      "; per-bucket top=", top_per_bucket, ")"
    )
  )

  ## ---------- Save gene list ----------
  gene_csv <- sub("\\.pdf$", "_genes.csv", out_pdf)
  fwrite(data.table(Gene = rownames(zmat)), gene_csv)

  message("✓ Immune functional heatmap: ", out_pdf)
  message("✓ Gene list saved: ", gene_csv)

  invisible(list(
    counts_global = counts_global,
    genes_use = genes_use,
    zmat = zmat
  ))
}

## ===================== 3. Load data =====================
pb_rds  <- safe_file(base_dir, "data_demo", "heatmap", "pseudobulk_by_hashtag_counts.rds")
deg_csv <- safe_file(base_dir, "data_demo", "heatmap", "DEG_edgeR_by_hashtag_pairwise.csv")

if (!file.exists(pb_rds)) stop("Missing file: ", pb_rds)
if (!file.exists(deg_csv)) stop("Missing file: ", deg_csv)

pb_tag    <- readRDS(pb_rds)
deg_pairs <- fread(deg_csv)

cat("===== pb_tag structure =====\n")
str(pb_tag, max.level = 2)

## ===================== 4. User settings =====================
target_hashtags <- c(
  "TotalSeqA-HT5",  # D3
  "TotalSeqA-HT8",  # D7
  "TotalSeqA-HT6",  # D9
  "TotalSeqA-HT7"   # D10
)

donor_map <- c(
  "TotalSeqA-HT5" = "D3",
  "TotalSeqA-HT6" = "D9",
  "TotalSeqA-HT7" = "D10",
  "TotalSeqA-HT8" = "D7"
)

out_pdf <- safe_file(
  base_dir,
  "output_example",
  "heatmap",
  "heatmap_IMMUNE_TOP_HT5678_D3_D7_D9_D10.pdf"
)

## ===================== 5. Run =====================
res <- immune_function_heatmap(
  pb_tag = pb_tag,
  deg_pairs = deg_pairs,
  target_hashtags = target_hashtags,
  out_pdf = out_pdf,
  donor_map = donor_map,
  top_per_bucket = 8,
  drop_mt_ribo = TRUE
)