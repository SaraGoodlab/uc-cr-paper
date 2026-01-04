#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  if (!flag %in% args) stop(flag, " not provided", call. = FALSE)
  args[match(flag, args) + 1]
}

input_dir <- get_arg("--input_dir")
summary_path <- get_arg("--summary")
table_path <- get_arg("--table")

result_files <- list.files(input_dir, pattern = "_results\\.csv$", full.names = TRUE)
if (length(result_files) == 0) {
  stop("No result files found in ", input_dir)
}

combined <- purrr::map_dfr(result_files, read_csv, show_col_types = FALSE)

summary_lines <- c("RNA-seq DEG summary")

for (file in result_files) {
  df <- read_csv(file, show_col_types = FALSE)
  tissue <- unique(df$tissue)
  comparison <- unique(df$comparison)
  degs <- df %>% filter(!is.na(padj), padj <= 0.05, abs(log2FoldChange) >= 1)
  up <- sum(degs$log2FoldChange > 0)
  down <- sum(degs$log2FoldChange < 0)
  summary_lines <- c(summary_lines, sprintf("%s | %s: %d DEGs (%d up, %d down)", tissue, comparison, nrow(degs), up, down))
}

write_lines(summary_lines, summary_path)

combined %>%
  select(tissue, comparison, gene_id, baseMean, log2FoldChange, stat, pvalue, padj) %>%
  arrange(tissue, comparison, padj) %>%
  write_csv(table_path)
