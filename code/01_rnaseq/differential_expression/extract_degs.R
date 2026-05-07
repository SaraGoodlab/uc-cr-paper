#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  if (!flag %in% args) stop(flag, " not provided", call. = FALSE)
  args[match(flag, args) + 1]
}

input_dir <- get_arg("--input_dir")
padj_threshold <- as.numeric(get_arg("--padj"))
log2fc_threshold <- as.numeric(get_arg("--log2fc"))
summary_path <- get_arg("--summary")
table_path <- get_arg("--table")

result_files <- list.files(input_dir, 
                           pattern = "_results\\.csv$",
                           full.names = TRUE) |>
  sort()

if (length(result_files) == 0) {
  stop("No result files found in ", input_dir)
}

combined <- purrr::map_dfr(result_files, read_csv, show_col_types = FALSE)
required_result_cols <- c(
  "tissue", "comparison", "gene_id", "baseMean",
  "log2FoldChange", "stat", "pvalue", "padj"
)
missing_result_cols <- setdiff(required_result_cols, names(combined))
if (length(missing_result_cols) > 0) {
  stop(
    "Combined results are missing required columns: ",
    paste(missing_result_cols, collapse = ", "),
    call. = FALSE
  )
}

summary_lines <- c("RNA-seq DEG summary")

for (file in result_files) {
  df <- read_csv(file, show_col_types = FALSE)
  missing_file_cols <- setdiff(required_result_cols, names(df))
  if (length(missing_file_cols) > 0) {
    stop(
      "Result file is missing required columns: ",
      basename(file),
      " -> ",
      paste(missing_file_cols, collapse = ", "),
      call. = FALSE
    )
  }
  tissue <- unique(df$tissue)
  comparison <- unique(df$comparison)
  if (length(tissue) != 1 || length(comparison) != 1) {
    stop("Each result file must contain one tissue and one comparison: ",
         basename(file), call. = FALSE)
  }
  degs <- df |>
    filter(!is.na(padj),
           padj <= padj_threshold,
           abs(log2FoldChange) >= log2fc_threshold)
  up <- sum(degs$log2FoldChange > 0)
  down <- sum(degs$log2FoldChange < 0)
  summary_lines <- c(summary_lines,
                     sprintf("%s | %s: %d DEGs (%d up, %d down)",
                             tissue, comparison, nrow(degs), up, down))
}

write_lines(summary_lines, summary_path)

combined |>
  select(tissue, comparison, gene_id, baseMean,
         log2FoldChange, stat, pvalue, padj) |>
  arrange(tissue, comparison, padj) |>
  write_csv(table_path)
