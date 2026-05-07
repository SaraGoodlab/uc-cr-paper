#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(tibble)
})

# Parse the CLI arguments
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
  if (!flag %in% args) stop(flag, " not provided", call. = FALSE)
  args[match(flag, args) + 1]
}

counts_path <- get_arg("--counts")
metadata_path <- get_arg("--metadata")
qc_output <- get_arg("--output")
size_factor_path <- get_arg("--size_factors")
excluded_path <- get_arg("--excluded")
threshold <- as.numeric(get_arg("--threshold"))

# Load the count matrix and confirm that genes are stored as rows.
counts_tbl <- read_tsv(counts_path, show_col_types = FALSE)
if (!"gene_id" %in% names(counts_tbl)) {
  stop("Counts table must contain a gene_id column", call. = FALSE)
}

genes <- counts_tbl$gene_id
count_mat <- counts_tbl |> select(-gene_id)

# Ensure metadata can be matched one-to-one to the count matrix columns.
sample_ids <- colnames(count_mat)
metadata <- read_tsv(metadata_path, show_col_types = FALSE)
required_metadata_cols <- c("sample_id", "tissue", "treatment")
missing_metadata_cols <- setdiff(required_metadata_cols, names(metadata))
if (length(missing_metadata_cols) > 0) {
  stop(
    "Metadata is missing required columns: ",
    paste(missing_metadata_cols, collapse = ", "),
    call. = FALSE
  )
}

metadata <- metadata |>
  mutate(
    sample_id = trimws(sample_id),
    tissue = tolower(trimws(tissue)),
    treatment = trimws(treatment)
  )
if (anyDuplicated(metadata$sample_id[metadata$sample_id != ""])) {
  dupes <- unique(metadata$sample_id[duplicated(metadata$sample_id)])
  stop("Duplicate sample_id values found in metadata: ",
       paste(dupes, collapse = ", "), call. = FALSE)
}

missing_metadata <- setdiff(sample_ids, metadata$sample_id)
if (length(missing_metadata) > 0) {
  stop(
    "Counts table columns missing from metadata: ",
    paste(missing_metadata, collapse = ", "),
    call. = FALSE
  )
}
missing_counts <- setdiff(metadata$sample_id[metadata$sample_id != ""],
                          sample_ids)
if (length(missing_counts) > 0) {
  stop(
    "Metadata sample_id values missing from counts table: ",
    paste(missing_counts, collapse = ", "),
    call. = FALSE
  )
}

metadata <- metadata |> filter(sample_id %in% sample_ids)
if (nrow(metadata) == 0) {
  stop("No metadata rows match the counts table columns", call. = FALSE)
}

# Estimate DESeq2 size factors in a design-free QC context.
count_mat <- as.matrix(count_mat[, metadata$sample_id, drop = FALSE])
storage.mode(count_mat) <- "numeric"
count_mat <- round(count_mat)
rownames(count_mat) <- genes

col_data <- metadata |> column_to_rownames("sample_id")

dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = col_data,
                              design = ~1)
dds <- estimateSizeFactors(dds)
size_factors <- sizeFactors(dds)
library_sizes <- colSums(count_mat)
sample_order <- rownames(col_data)

# Record both the full QC summary and the derived exclusion list
qc_tbl <- tibble(
  sample_id = sample_order,
  tissue = col_data$tissue,
  treatment = col_data$treatment,
  library_size = library_sizes[sample_order],
  size_factor = size_factors[sample_order],
  flagged = size_factors[sample_order] < threshold
)

qc_tbl <- qc_tbl |> arrange(tissue, treatment, sample_id)

excluded <- qc_tbl |>
  filter(flagged) |>
  arrange(sample_id) |>
  pull(sample_id)
write_lines(excluded, excluded_path)
write_csv(qc_tbl, qc_output)

size_factor_tbl <- qc_tbl |> select(sample_id, size_factor)
write_csv(size_factor_tbl, size_factor_path)
