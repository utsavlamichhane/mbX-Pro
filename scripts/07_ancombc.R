#!/usr/bin/env Rscript

# 07_ancombc.R
# Differential abundance using ANCOM-BC (ANCOMBC::ancombc2)
#
# Inputs:
# - microbiome table CSV (default: level-7.csv) exported from QIIME2 taxa collapse
# - metadata TSV
# - group_var (metadata variable to test)
#
# Outputs:
# - ancombc_out/ with result tables + sessionInfo.txt

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat("
Usage:
  Rscript scripts/07_ancombc.R --microbiome level-7.csv --metadata metadata.txt --group_var BMIClass --out_dir ancombc_out

Required:
  --metadata    Metadata file
  --group_var   Metadata variable name to test

Optional:
  --microbiome  Microbiome table CSV (default: level-7.csv)
  --out_dir     Output directory (default: ancombc_out)

  # ANCOM-BC cutoffs (defaults are common; edit as needed)
  --prv_cut     Prevalence cutoff (default: 0 = no filter)
  --lib_cut     Library size cutoff (default: 0 = no filter)
  --alpha       Significance level (default: 0.05)
  --p_adj       p-adjust method (default: BH)
  --global      TRUE/FALSE (default: FALSE)
  --pairwise    TRUE/FALSE (default: FALSE)
")
}

get_arg <- function(flag, default = NA) {
  ix <- which(args == flag)
  if (length(ix) == 0) return(default)
  if (ix == length(args)) stop(paste0("Missing value for ", flag))
  args[ix + 1]
}

microbiome <- get_arg("--microbiome", "level-7.csv")
metadata   <- get_arg("--metadata", NA)
group_var  <- get_arg("--group_var", NA)
out_dir    <- get_arg("--out_dir", "ancombc_out")

prv_cut <- as.numeric(get_arg("--prv_cut", "0"))
lib_cut <- as.numeric(get_arg("--lib_cut", "0"))
alpha   <- as.numeric(get_arg("--alpha", "0.05"))
p_adj   <- get_arg("--p_adj", "BH")

global  <- tolower(get_arg("--global", "FALSE")) == "true"
pairwise<- tolower(get_arg("--pairwise", "FALSE")) == "true"

if (is.na(metadata) || is.na(group_var)) {
  usage()
  stop("ERROR: --metadata and --group_var are required.")
}
if (!file.exists(microbiome)) stop(paste0("ERROR: microbiome file not found: ", microbiome))
if (!file.exists(metadata)) stop(paste0("ERROR: metadata file not found: ", metadata))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Install deps (Bioconductor) if missing ----
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc_install_if_missing <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      BiocManager::install(p, ask = FALSE, update = FALSE)
    }
  }
}

bioc_install_if_missing(c("ANCOMBC", "TreeSummarizedExperiment", "S4Vectors"))

suppressWarnings(suppressMessages({
  library(ANCOMBC)
  library(TreeSummarizedExperiment)
  library(S4Vectors)
}))

# ---- Read microbiome table (taxa x samples) ----
df <- read.csv(microbiome, check.names = FALSE, stringsAsFactors = FALSE)

tax_col <- names(df)[1]
taxa <- df[[1]]
df <- df[, -1, drop = FALSE]

# Make numeric matrix
to_num <- function(x) suppressWarnings(as.numeric(x))
mat <- as.matrix(data.frame(lapply(df, to_num), check.names = FALSE))
rownames(mat) <- taxa

# ---- Read metadata (samples x vars) ----
md <- read.delim(metadata, check.names = FALSE, stringsAsFactors = FALSE)

# Identify sample-id column (QIIME2 often uses "#SampleID" or "sample-id")
sample_col <- NULL
for (cand in c("#SampleID", "sample-id", "SampleID", "sampleid", "id", "ID")) {
  if (cand %in% names(md)) { sample_col <- cand; break }
}
if (is.null(sample_col)) {
  # fallback: first column
  sample_col <- names(md)[1]
}

rownames(md) <- md[[sample_col]]

if (!(group_var %in% names(md))) {
  stop(paste0("ERROR: group_var not found in metadata columns: ", group_var))
}

# ---- Align samples ----
common <- intersect(colnames(mat), rownames(md))
if (length(common) < 2) {
  stop("ERROR: Not enough overlapping sample IDs between microbiome table columns and metadata sample IDs.")
}

mat <- mat[, common, drop = FALSE]
md  <- md[common, , drop = FALSE]

# ---- Build TreeSummarizedExperiment ----
# Parse taxonomy string into ranks if possible
parse_tax <- function(x) {
  parts <- unlist(strsplit(x, ";\\s*"))
  parts <- trimws(parts)
  # strip prefixes like k__/p__/g__/s__
  strip <- function(z) sub("^[a-z]__","", z)
  parts2 <- vapply(parts, strip, character(1))
  # pad to 7 ranks
  if (length(parts2) < 7) parts2 <- c(parts2, rep(NA_character_, 7 - length(parts2)))
  parts2[1:7]
}

tax_mat <- t(vapply(rownames(mat), parse_tax, character(7)))
colnames(tax_mat) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
row_df <- S4Vectors::DataFrame(tax_mat, row.names = rownames(mat))

col_df <- S4Vectors::DataFrame(md, row.names = rownames(md))

tse <- TreeSummarizedExperiment(
  assays = SimpleList(counts = mat),
  rowData = row_df,
  colData = col_df
)

fix_formula <- as.formula(paste0("~ ", group_var))

# ---- Run ANCOM-BC2 ----
message("[INFO] Running ANCOM-BC2 ...")
out <- ancombc2(
  data = tse,
  assay_name = "counts",
  fix_formula = fix_formula,
  p_adj_method = p_adj,
  alpha = alpha,
  prv_cut = prv_cut,
  lib_cut = lib_cut,
  group = group_var,
  global = global,
  pairwise = pairwise,
  verbose = TRUE
)

# ---- Save outputs ----
saveRDS(out, file.path(out_dir, "ancombc2_output.rds"))

# Try to write common result tables
# ancombc2 returns a list with $res in the official tutorial.
if (!is.null(out$res)) {
  res <- out$res
  if (is.data.frame(res)) {
    write.csv(res, file.path(out_dir, "results_res.csv"), row.names = FALSE)
  } else if (is.list(res)) {
    # write each element
    for (nm in names(res)) {
      obj <- res[[nm]]
      fn <- file.path(out_dir, paste0("res_", nm, ".csv"))
      if (is.data.frame(obj)) {
        write.csv(obj, fn, row.names = FALSE)
      } else if (is.matrix(obj)) {
        write.csv(as.data.frame(obj), fn, row.names = TRUE)
      }
    }
  }
}

# Session info for reproducibility
sink(file.path(out_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()

message(paste0("[DONE] ANCOM-BC results in: ", normalizePath(out_dir, winslash = "/")))
