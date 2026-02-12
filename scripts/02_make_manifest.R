#!/usr/bin/env Rscript

# 02_make_manifest.R
# Create a QIIME 2 manifest file (PairedEndFastqManifestPhred33V2) from a FASTQ directory.
#
# Based on the overall logic of your Manifest_creator.R (enumerate files → infer sample + direction → write manifest),
# but with robust filename parsing for common patterns:
#   - Casava 1.8: <SAMPLE>_S1_L001_R1_001.fastq.gz
#   - SRA Toolkit: <SAMPLE>_1.fastq.gz / <SAMPLE>_2.fastq.gz
#   - Common: <SAMPLE>_R1.fastq.gz / <SAMPLE>_R2.fastq.gz

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat("
Usage:
  Rscript scripts/02_make_manifest.R --fastq_dir data/fastq --read_type paired --out manifests/manifest_paired.txt

Required:
  --fastq_dir   Directory containing FASTQ(.gz) files
  --read_type   paired | single
  --out         Output manifest file path

Optional:
  --pattern     Regex pattern to match FASTQ files (default: \\\\.(fastq|fq)(\\\\.gz)?$)
")
}

get_arg <- function(flag, default = NA) {
  ix <- which(args == flag)
  if (length(ix) == 0) return(default)
  if (ix == length(args)) stop(paste0("Missing value for ", flag))
  args[ix + 1]
}

fastq_dir <- get_arg("--fastq_dir", NA)
read_type <- tolower(get_arg("--read_type", NA))
out_file  <- get_arg("--out", NA)
pattern   <- get_arg("--pattern", "\\.(fastq|fq)(\\.gz)?$")

if (is.na(fastq_dir) || is.na(read_type) || is.na(out_file)) {
  usage()
  stop("ERROR: --fastq_dir, --read_type, and --out are required.")
}
if (!dir.exists(fastq_dir)) stop(paste0("ERROR: fastq_dir does not exist: ", fastq_dir))
if (!(read_type %in% c("paired","single"))) stop("ERROR: --read_type must be 'paired' or 'single'")

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

files <- list.files(fastq_dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
if (length(files) == 0) stop(paste0("ERROR: No FASTQ files found in ", fastq_dir))

abs_files <- normalizePath(files, winslash = "/", mustWork = TRUE)
base_names <- basename(abs_files)

infer_read <- function(fn) {
  f <- tolower(fn)
  if (grepl("_r1_001\\.(fastq|fq)(\\.gz)?$", f) ||
      grepl("([_\\.])r1\\.(fastq|fq)(\\.gz)?$", f) ||
      grepl("_1\\.(fastq|fq)(\\.gz)?$", f)) return("forward")
  if (grepl("_r2_001\\.(fastq|fq)(\\.gz)?$", f) ||
      grepl("([_\\.])r2\\.(fastq|fq)(\\.gz)?$", f) ||
      grepl("_2\\.(fastq|fq)(\\.gz)?$", f)) return("reverse")
  return(NA_character_)
}

infer_sample <- function(fn) {
  # Casava 1.8: SAMPLE_S1_L001_R1_001.fastq.gz
  if (grepl("^(.*)_S[0-9]+_L[0-9]+_R[12]_001\\.(fastq|fq)(\\.gz)?$", fn, ignore.case = TRUE, perl = TRUE)) {
    return(sub("^(.*)_S[0-9]+_L[0-9]+_R[12]_001\\.(fastq|fq)(\\.gz)?$", "\\1", fn, ignore.case = TRUE, perl = TRUE))
  }

  # SAMPLE_R1.fastq.gz or SAMPLE.R1.fastq.gz (optionally with extra tokens after R1)
  if (grepl("^(.*?)([_\\.])R[12].*\\.(fastq|fq)(\\.gz)?$", fn, ignore.case = TRUE, perl = TRUE)) {
    return(sub("^(.*?)([_\\.])R[12].*\\.(fastq|fq)(\\.gz)?$", "\\1", fn, ignore.case = TRUE, perl = TRUE))
  }

  # SAMPLE_1.fastq.gz / SAMPLE_2.fastq.gz
  if (grepl("^(.*)_([12])\\.(fastq|fq)(\\.gz)?$", fn, ignore.case = TRUE, perl = TRUE)) {
    return(sub("^(.*)_([12])\\.(fastq|fq)(\\.gz)?$", "\\1", fn, ignore.case = TRUE, perl = TRUE))
  }

  # Fallback: capture "S1", "Sample-12", etc.
  pat <- "(?i)(s(?:ample)?[-_]?[0-9]+)"
  m <- regexpr(pat, fn, perl = TRUE)
  if (m[1] != -1) return(toupper(regmatches(fn, m)[[1]]))

  stop(paste0("ERROR: Could not infer sample-id from filename: ", fn))
}

direction <- vapply(base_names, infer_read, character(1))
if (any(is.na(direction))) {
  bad <- base_names[is.na(direction)]
  stop(paste0("ERROR: Could not infer forward/reverse from these files:\n  - ", paste(bad, collapse = "\n  - "),
              "\nExpected patterns include: *_R1_001.fastq.gz / *_R2_001.fastq.gz, *_1.fastq.gz / *_2.fastq.gz, or *_R1.fastq.gz / *_R2.fastq.gz"))
}

sample_id <- vapply(base_names, infer_sample, character(1))

df <- data.frame(
  sample_id = sample_id,
  direction = direction,
  abs_path = abs_files,
  stringsAsFactors = FALSE
)

# If single-end: keep only forward reads
if (read_type == "single") {
  df2 <- df[df$direction == "forward", ]
  df2 <- df2[!duplicated(paste(df2$sample_id, df2$abs_path)), ]
  df2 <- df2[order(df2$sample_id), c("sample_id", "abs_path")]
  colnames(df2) <- c("sample-id", "absolute-filepath")
  write.table(df2, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste0("[DONE] Wrote single-end manifest: ", out_file, "\n"))
  quit(save = "no", status = 0)
}

# Paired-end
samples <- sort(unique(df$sample_id))
out <- data.frame(
  `sample-id` = character(0),
  `forward-absolute-filepath` = character(0),
  `reverse-absolute-filepath` = character(0),
  stringsAsFactors = FALSE
)

for (s in samples) {
  subdf <- df[df$sample_id == s, , drop = FALSE]
  fwd <- subdf$abs_path[subdf$direction == "forward"]
  rev <- subdf$abs_path[subdf$direction == "reverse"]

  if (length(fwd) != 1 || length(rev) != 1) {
    stop(paste0(
      "ERROR: Sample '", s, "' does not have exactly one forward and one reverse FASTQ.\n",
      "Forward files: ", paste(fwd, collapse = ", "), "\n",
      "Reverse files: ", paste(rev, collapse = ", ")
    ))
  }

  out <- rbind(out, data.frame(
    `sample-id` = s,
    `forward-absolute-filepath` = fwd,
    `reverse-absolute-filepath` = rev,
    stringsAsFactors = FALSE
  ))
}

write.table(out, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste0("[DONE] Wrote paired-end manifest: ", out_file, "\n"))
