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
excluded_path <- get_arg("--excluded")
tissue_filter <- tolower(get_arg("--tissue"))
comparison_name <- get_arg("--comparison")
case_label <- get_arg("--case")
control_label <- get_arg("--control")
padj_threshold <- as.numeric(get_arg("--padj"))
log2fc_threshold <- as.numeric(get_arg("--log2fc"))
results_path <- get_arg("--results")
deg_path <- get_arg("--degs")
minimum_replicates <- 2

# Load the merged count matrix and confirm the expected gene-row layout
counts_tbl <- read_tsv(counts_path, show_col_types = FALSE)
if (!"gene_id" %in% names(counts_tbl)) {
  stop("Counts table must contain a gene_id column", call. = FALSE)
}

# Validate and normalize metadata fields
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

if (file.exists(excluded_path)) {
  to_exclude <- sort(unique(trimws(read_lines(excluded_path))))
} else {
  to_exclude <- character()
}
to_exclude <- to_exclude[to_exclude != ""]
if (length(to_exclude) > 0) {
  metadata <- metadata |> filter(!(sample_id %in% to_exclude))
}

# Keep only the requested tissue and the two treatment groups for this contrast.
metadata <- metadata |>
  filter(tissue == tissue_filter, treatment %in% c(control_label, case_label))
if (nrow(metadata) == 0) {
  stop(paste("No samples remain for DESeq2 in", comparison_name,
             "for tissue", tissue_filter), call. = FALSE)
}

count_tbl_only <- counts_tbl |> select(-gene_id)
sample_ids <- colnames(count_tbl_only)
missing_counts <- setdiff(metadata$sample_id, sample_ids)
if (length(missing_counts) > 0) {
  stop(
    "Metadata sample_id values missing from counts table for ",
    comparison_name,
    ": ",
    paste(missing_counts, collapse = ", "),
    call. = FALSE
  )
}

count_mat <- as.matrix(count_tbl_only[, metadata$sample_id, drop = FALSE])
rownames(count_mat) <- counts_tbl$gene_id
count_mat <- round(count_mat)

# Set the control group as the reference level
col_data <- metadata |> column_to_rownames("sample_id")
col_data <- col_data |>
  mutate(treatment = factor(treatment, levels = c(control_label, case_label)))

# Fit the model and export both the full table and DEG hits.
dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = col_data,
                              design = ~ treatment)
dds <- DESeq(dds)
res <- results(dds, contrast = c("treatment", case_label, control_label))
res_tbl <- as_tibble(res, rownames = "gene_id") |>
  mutate(comparison = comparison_name, tissue = tissue_filter) |>
  arrange(gene_id)

sig_tbl <- res_tbl |>
  filter(!is.na(padj)) |>
  filter(padj <= padj_threshold, abs(log2FoldChange) >= log2fc_threshold) |>
  arrange(padj, desc(abs(log2FoldChange)), gene_id)

write_csv(res_tbl, results_path)
write_lines(sig_tbl$gene_id, deg_path)
