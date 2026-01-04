#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tximport)
  library(readr)
  library(dplyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  if (!flag %in% args) stop(flag, " not provided", call. = FALSE)
  args[match(flag, args) + 1]
}

rsem_dir <- get_arg("--rsem_dir")
metadata_path <- get_arg("--metadata")
output_path <- get_arg("--output")

metadata <- read_tsv(metadata_path, show_col_types = FALSE) %>%
  mutate(sample_id = trimws(sample_id))

if (!"sample_id" %in% names(metadata)) {
  stop("metadata is missing a sample_id column", call. = FALSE)
}

sample_order <- metadata$sample_id[metadata$sample_id != ""]
if (length(sample_order) == 0) {
  stop("no sample_id entries found in metadata", call. = FALSE)
}

rsem_files <- list.files(
  rsem_dir,
  pattern = "\\.genes\\.results$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(rsem_files) == 0) {
  stop("no RSEM .genes.results files found under ", rsem_dir, call. = FALSE)
}

sample_names <- sub("\\.genes\\.results$", "", basename(rsem_files))
if (anyDuplicated(sample_names)) {
  dupes <- unique(sample_names[duplicated(sample_names)])
  stop("duplicate RSEM sample names found: ", paste(dupes, collapse = ", "), call. = FALSE)
}

file_lookup <- setNames(rsem_files, sample_names)
missing_samples <- setdiff(sample_order, names(file_lookup))
if (length(missing_samples) > 0) {
  stop("no RSEM files found for samples: ", paste(missing_samples, collapse = ", "), call. = FALSE)
}

files <- file_lookup[sample_order]
txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE, countsFromAbundance = "no")

counts <- txi$counts
counts_tbl <- as_tibble(counts, rownames = "gene_id")

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write_tsv(counts_tbl, output_path)
