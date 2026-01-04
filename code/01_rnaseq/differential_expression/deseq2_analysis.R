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
excluded_path <- get_arg("--excluded")
tissue_filter <- tolower(get_arg("--tissue"))
comparison_name <- get_arg("--comparison")
case_label <- get_arg("--case")
control_label <- get_arg("--control")
padj_threshold <- as.numeric(get_arg("--padj"))
log2fc_threshold <- as.numeric(get_arg("--log2fc"))
results_path <- get_arg("--results")
deg_path <- get_arg("--degs")

counts_tbl <- read_tsv(counts_path, show_col_types = FALSE)
metadata <- read_tsv(metadata_path, show_col_types = FALSE) %>%
  mutate(
    sample_id = trimws(sample_id),
    tissue = tolower(trimws(tissue)),
    treatment = trimws(treatment)
  )

to_exclude <- if (file.exists(excluded_path)) read_lines(excluded_path) else character()
if (length(to_exclude) > 0) {
  metadata <- metadata %>% filter(!(sample_id %in% to_exclude))
}

metadata <- metadata %>% filter(tissue == tissue_filter, treatment %in% c(control_label, case_label))
if (nrow(metadata) < 4) {
  stop(paste("Not enough samples for DESeq2 in", comparison_name, "- need at least 4"), call. = FALSE)
}

count_tbl_only <- counts_tbl %>% select(-gene_id)
count_mat <- as.matrix(count_tbl_only[, metadata$sample_id, drop = FALSE])
rownames(count_mat) <- counts_tbl$gene_id
count_mat <- round(count_mat)

col_data <- metadata %>% column_to_rownames("sample_id")
col_data <- col_data %>% mutate(treatment = factor(treatment, levels = c(control_label, case_label)))

dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = col_data, design = ~ treatment)
dds <- DESeq(dds)
res <- results(dds, contrast = c("treatment", case_label, control_label))
res_tbl <- as_tibble(res, rownames = "gene_id") %>% mutate(comparison = comparison_name, tissue = tissue_filter)

sig_tbl <- res_tbl %>%
  filter(!is.na(padj)) %>%
  filter(padj <= padj_threshold, abs(log2FoldChange) >= log2fc_threshold)

write_csv(res_tbl, results_path)
write_lines(sig_tbl$gene_id, deg_path)
