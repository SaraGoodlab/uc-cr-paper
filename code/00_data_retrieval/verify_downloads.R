#!/usr/bin/env Rscript
# verify_downloads.R
# Verify downloaded ENA sequencing data.
#
# Usage:
#   Rscript verify_downloads.R

suppressPackageStartupMessages({
  library(dplyr)
  library(fs)
  library(readr)
  library(stringr)
})

script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script <- sub(file_arg, "", args[grepl(file_arg, args)])
  if (length(script) == 0) {
    return(normalizePath(getwd()))
  }
  normalizePath(script)
}

script_dir <- path_dir(script_path())
repo_root <- path_norm(path(script_dir, "..", ".."))

fmt_gb <- function(bytes) sprintf("%.2f GB", as.numeric(bytes) / 1e9)

list_fastqs <- function(directory) {
  if (!dir_exists(directory)) {
    return(character())
  }
  dir_ls(directory, regexp = "\\.(fastq|fq)(\\.gz)?$", recurse = FALSE)
}

summarize_files <- function(files) {
  tibble(
    files = length(files),
    size_gb = fmt_gb(sum(file_size(files)))
  )
}

warn_zero_sized <- function(files, label) {
  zero <- files[file_size(files) == 0]
  if (length(zero) > 0) {
    warning(sprintf("[%s] %d zero-size files detected", label, length(zero)), call. = FALSE)
    print(zero)
  }
}

base_raw <- path(repo_root, "data", "raw")
rnaseq_dir <- path(base_raw, "rnaseq")
micro_dir <- path(base_raw, "microbiome")
metadata_dir <- path(repo_root, "data", "metadata")

cat("Verifying downloaded sequencing data...\n\n")

# RNA-seq files ------------------------------------------------------------
cat("=== RNA-seq Data ===\n")
if (dir_exists(rnaseq_dir)) {
  colon_files <- list_fastqs(path(rnaseq_dir, "colon"))
  spleen_files <- list_fastqs(path(rnaseq_dir, "spleen"))

  colon_summary <- summarize_files(colon_files)
  spleen_summary <- summarize_files(spleen_files)

  cat(sprintf("Colon files: %d (%s)\n", colon_summary$files, colon_summary$size_gb))
  cat(sprintf("Spleen files: %d (%s)\n", spleen_summary$files, spleen_summary$size_gb))

  total_rnaseq <- length(colon_files) + length(spleen_files)
  if (total_rnaseq %% 2 != 0) {
    warning("Total RNA-seq FASTQs is odd; expected paired-end files.", call. = FALSE)
  }

  warn_zero_sized(c(colon_files, spleen_files), "RNA-seq")
} else {
  cat("RNA-seq directory not found. Run download_ena_rnaseq.sh first.\n")
  colon_files <- spleen_files <- character()
}

cat("\n")

# Microbiome files ---------------------------------------------------------
cat("=== Microbiome Data ===\n")
if (dir_exists(micro_dir)) {
  micro_files <- list_fastqs(micro_dir)
  micro_summary <- summarize_files(micro_files)
  cat(sprintf("FASTQ files: %d (%s)\n", micro_summary$files, micro_summary$size_gb))

  metadata_path <- path(metadata_dir, "microbiome.tsv")
  if (file_exists(metadata_path)) {
    metadata <- read_tsv(metadata_path, show_col_types = FALSE)
    expected_ids <- metadata$sample_id %>% str_to_lower()
    observed_ids <- path_file(micro_files) %>%
      str_remove("\\.(fastq|fq)(\\.gz)?$") %>%
      str_to_lower()

    missing_ids <- setdiff(expected_ids, observed_ids)
    if (length(missing_ids) > 0) {
      warning(sprintf("%d samples from metadata have no corresponding FASTQ", length(missing_ids)), call. = FALSE)
      print(missing_ids)
    }
  } else {
    cat("Microbiome metadata not found; skipping sample completeness check.\n")
  }

  warn_zero_sized(micro_files, "Microbiome")
} else {
  cat("Microbiome directory not found. Run download_ena_microbiome.sh first.\n")
  micro_files <- character()
}

cat("\n=== Verification Complete ===\n")
cat(sprintf("\nSummary:\n  RNA-seq files: %d\n  Microbiome files: %d\n",
            length(colon_files) + length(spleen_files),
            length(micro_files)))
