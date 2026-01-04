#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(tibble)
})

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

counts_tbl <- read_tsv(counts_path, show_col_types = FALSE)
if (!"gene_id" %in% names(counts_tbl)) {
  stop("Counts table must contain a gene_id column", call. = FALSE)
}

genes <- counts_tbl$gene_id
count_mat <- counts_tbl %>% select(-gene_id)

sample_ids <- colnames(count_mat)
metadata <- read_tsv(metadata_path, show_col_types = FALSE) %>%
  mutate(
    sample_id = trimws(sample_id),
    tissue = tolower(trimws(tissue)),
    treatment = trimws(treatment)
  )

metadata <- metadata %>% filter(sample_id %in% sample_ids)
if (nrow(metadata) == 0) {
  stop("No metadata rows match the counts table columns", call. = FALSE)
}

count_mat <- as.matrix(count_mat[, metadata$sample_id, drop = FALSE])
storage.mode(count_mat) <- "numeric"
count_mat <- round(count_mat)
rownames(count_mat) <- genes

col_data <- metadata %>% column_to_rownames("sample_id")

dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = col_data, design = ~1)
dds <- estimateSizeFactors(dds)
size_factors <- sizeFactors(dds)
library_sizes <- colSums(count_mat)
sample_order <- rownames(col_data)

qc_tbl <- tibble(
  sample_id = sample_order,
  tissue = col_data$tissue,
  treatment = col_data$treatment,
  library_size = library_sizes[sample_order],
  size_factor = size_factors[sample_order],
  flagged = size_factors[sample_order] < threshold
)

qc_tbl <- qc_tbl %>% arrange(flagged, tissue, treatment)

excluded <- qc_tbl %>% filter(flagged) %>% pull(sample_id)
write_lines(excluded, excluded_path)
write_csv(qc_tbl, qc_output)

size_factor_tbl <- qc_tbl %>% select(sample_id, size_factor)
write_csv(size_factor_tbl, size_factor_path)
