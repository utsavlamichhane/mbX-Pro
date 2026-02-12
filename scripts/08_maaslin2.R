#!/usr/bin/env Rscript

# 08_maaslin2.R
# Differential abundance/association using MaAsLin2
#
# Inputs:
# - microbiome table CSV (taxa x samples) from QIIME2 taxa collapse export
# - metadata TSV
# - group_var (metadata variable to test; can be comma-separated for multiple fixed effects)
#
# Outputs:
# - maaslin2_out/ with MaAsLin2 results + sessionInfo.txt

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat("
Usage:
  Rscript scripts/08_maaslin2.R --microbiome level-7.csv --metadata metadata.txt --group_var BMIClass --out_dir maaslin2_out

Required:
  --metadata     Metadata file
  --group_var    Fixed effect(s). Comma-separated allowed, e.g. BMIClass,Sex

Optional:
  --microbiome   Microbiome table CSV (default: level-7.csv)
  --out_dir      Output directory (default: maaslin2_out)

  # MaAsLin2 processing choices (defaults are conservative / minimal transforms)
  --normalization  (default: NONE)
  --transform      (default: NONE)
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
out_dir    <- get_arg("--out_dir", "maaslin2_out")
norm       <- get_arg("--normalization", "NONE")
trans      <- get_arg("--transform", "NONE")

if (is.na(metadata) || is.na(group_var)) {
  usage()
  stop("ERROR: --metadata and --group_var are required.")
}
if (!file.exists(microbiome)) stop(paste0("ERROR: microbiome file not found: ", microbiome))
if (!file.exists(metadata)) stop(paste0("ERROR: metadata file not found: ", metadata))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Install/load MaAsLin2 if missing
if (!requireNamespace("Maaslin2", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install("Maaslin2", ask = FALSE, update = FALSE)
}
suppressWarnings(suppressMessages(library(Maaslin2)))

# Read microbiome (taxa x samples), then transpose to samples x features for MaAsLin2
df <- read.csv(microbiome, check.names = FALSE, stringsAsFactors = FALSE)
taxa <- df[[1]]
df <- df[, -1, drop = FALSE]

to_num <- function(x) suppressWarnings(as.numeric(x))
mat <- as.matrix(data.frame(lapply(df, to_num), check.names = FALSE))
rownames(mat) <- taxa

# Transpose: samples rows, features columns
mat_t <- t(mat)
input_data <- as.data.frame(mat_t, check.names = FALSE)

# Read metadata; set rownames to sample IDs
md <- read.delim(metadata, check.names = FALSE, stringsAsFactors = FALSE)

sample_col <- NULL
for (cand in c("#SampleID", "sample-id", "SampleID", "sampleid", "id", "ID")) {
  if (cand %in% names(md)) { sample_col <- cand; break }
}
if (is.null(sample_col)) sample_col <- names(md)[1]
rownames(md) <- md[[sample_col]]

# Align
common <- intersect(rownames(input_data), rownames(md))
if (length(common) < 2) {
  stop("ERROR: Not enough overlapping sample IDs between microbiome table and metadata.")
}
input_data <- input_data[common, , drop = FALSE]
md <- md[common, , drop = FALSE]

fixed_effects <- trimws(unlist(strsplit(group_var, ",")))
missing_fx <- fixed_effects[!(fixed_effects %in% names(md))]
if (length(missing_fx) > 0) {
  stop(paste0("ERROR: These fixed effects are not in metadata: ", paste(missing_fx, collapse = ", ")))
}

message("[INFO] Running MaAsLin2 ...")
fit <- Maaslin2(
  input_data = input_data,
  input_metadata = md,
  output = out_dir,
  fixed_effects = fixed_effects,
  normalization = norm,
  transform = trans
)

# Session info for reproducibility
sink(file.path(out_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()

message(paste0("[DONE] MaAsLin2 results in: ", normalizePath(out_dir, winslash = "/")))
