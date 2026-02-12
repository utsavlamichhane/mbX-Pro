#!/usr/bin/env Rscript

# 05_run_mbx.R
# Runs mbX: ezclean, ezviz, ezstat
# Always use level-7.csv by default, but you can point to a different file via --microbiome.
#
# The script auto-detects the mbX `level` argument from the taxonomy strings if you don't pass --level.

suppressWarnings(suppressMessages({
  library(utils)
}))

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat("
Usage:
  Rscript scripts/05_run_mbx.R --microbiome level-7.csv --metadata metadata.txt --group_var BMIClass

Required:
  --metadata     Metadata file (tab-delimited; QIIME2-style works)
  --group_var    Metadata column name to use for ezviz/ezstat

Optional:
  --microbiome   Microbiome table CSV (default: level-7.csv)
  --level        mbX level code (k/p/c/o/f/g/s). If omitted, auto-detected from taxonomy strings.
  --top_taxa     For ezviz (default: 10)
  --flip         For ezviz (default: False)
  --out_dir      Output directory (default: mbx_out)
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
level_arg  <- get_arg("--level", NA)
top_taxa   <- as.integer(get_arg("--top_taxa", "10"))
flip       <- get_arg("--flip", "False")
out_dir    <- get_arg("--out_dir", "mbx_out")

if (is.na(metadata) || is.na(group_var)) {
  usage()
  stop("ERROR: --metadata and --group_var are required.")
}
if (!file.exists(microbiome)) stop(paste0("ERROR: microbiome file not found: ", microbiome))
if (!file.exists(metadata)) stop(paste0("ERROR: metadata file not found: ", metadata))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Install/load mbX from r-universe if missing ---
if (!requireNamespace("mbX", quietly = TRUE)) {
  message("[INFO] Installing mbX from r-universe...")
  options(repos = c(
    "https://cloud.r-project.org",
    "https://utsavlamichhane.r-universe.dev"
  ))
  install.packages("mbX")
}
suppressWarnings(suppressMessages(library(mbX)))

# --- Infer mbX level code from taxonomy strings if not provided ---
infer_level <- function(tax_strings) {
  # Accepts strings like: "d__Bacteria; p__Firmicutes; ...; g__Blautia; s__Blautia_wexlerae"
  # Returns one of: k,p,c,o,f,g,s (default to g if none found)
  tax <- tax_strings[!is.na(tax_strings) & nzchar(tax_strings)]
  if (length(tax) == 0) return("g")
  tax <- tax[1:min(500, length(tax))]
  txt <- paste(tax, collapse = " ")
  if (grepl("s__", txt, fixed = TRUE)) return("s")
  if (grepl("g__", txt, fixed = TRUE)) return("g")
  if (grepl("f__", txt, fixed = TRUE)) return("f")
  if (grepl("o__", txt, fixed = TRUE)) return("o")
  if (grepl("c__", txt, fixed = TRUE)) return("c")
  if (grepl("p__", txt, fixed = TRUE)) return("p")
  if (grepl("k__", txt, fixed = TRUE) || grepl("d__", txt, fixed = TRUE)) return("k")
  return("g")
}

level_code <- level_arg
if (is.na(level_code)) {
  df0 <- read.csv(microbiome, check.names = FALSE, stringsAsFactors = FALSE)
  first_col <- df0[[1]]
  level_code <- infer_level(first_col)
}

writeLines(paste0("mbX level used: ", level_code), file.path(out_dir, "mbx_level_used.txt"))

# Run inside output directory so plots/results land here
old_wd <- getwd()
setwd(out_dir)
on.exit(setwd(old_wd), add = TRUE)

# Copy inputs for provenance
file.copy(microbiome, "level-7.csv", overwrite = TRUE)
file.copy(metadata, "metadata.txt", overwrite = TRUE)

message("[INFO] Running ezclean...")
ezclean(microbiome_data = "level-7.csv", metadata = "metadata.txt", level = level)

message("[INFO] Running ezviz...")
ezviz("level-7.csv", "metadata.txt", level, group_var, top_taxa = top_taxa, flip = flip)

message("[INFO] Running ezstat...")
ezstat(microbiome_data = "level-7.csv", metadata = "metadata.txt", level = level, selected_metadata = group_var)

# Session info for reproducibility
sink("sessionInfo.txt")
print(sessionInfo())
sink()

message(paste0("[DONE] mbX outputs are in: ", normalizePath(out_dir, winslash = "/")))
