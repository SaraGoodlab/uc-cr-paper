#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  if (!flag %in% args) stop(flag, " not provided", call. = FALSE)
  args[match(flag, args) + 1]
}

counts_path <- get_arg("--counts")
metadata_path <- get_arg("--metadata")
excluded_path <- get_arg("--excluded")
plot_path <- get_arg("--plot")
output_path <- get_arg("--output")

counts_tbl <- read_tsv(counts_path, show_col_types = FALSE)
metadata <- read_tsv(metadata_path, show_col_types = FALSE) %>%
  mutate(
    sample_id = trimws(sample_id),
    tissue = tolower(trimws(tissue)),
    treatment = trimws(treatment)
  )

to_exclude <- if (file.exists(excluded_path)) read_lines(excluded_path) else character()
gene_ids <- counts_tbl$gene_id
count_tbl_only <- counts_tbl %>% select(-gene_id)
sample_ids <- colnames(count_tbl_only)
metadata <- metadata %>% filter(sample_id %in% sample_ids, !(sample_id %in% to_exclude))
if (nrow(metadata) < 3) {
  stop("Need at least three samples for PCA", call. = FALSE)
}

count_mat <- as.matrix(count_tbl_only[, metadata$sample_id, drop = FALSE])
rownames(count_mat) <- gene_ids
count_mat <- round(count_mat)
col_data <- metadata %>% column_to_rownames("sample_id")

dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = col_data, design = ~1)
dds <- estimateSizeFactors(dds)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
vsd_mat <- assay(vsd)
pca <- prcomp(t(vsd_mat), center = TRUE, scale. = FALSE)

pca_df <- as.data.frame(pca$x[, 1:2, drop = FALSE]) %>%
  mutate(
    sample_id = rownames(pca$x),
    tissue = col_data[sample_id, "tissue"],
    treatment = col_data[sample_id, "treatment"]
  )

write_csv(pca_df, output_path)

percent_var <- pca$sdev^2 / sum(pca$sdev^2)
percent_label <- function(x) sprintf("%.1f%%", x * 100)

plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = treatment, shape = tissue)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "RNA-seq PCA",
    x = paste0("PC1 (", percent_label(percent_var[1]), ")"),
    y = paste0("PC2 (", percent_label(percent_var[2]), ")")
  ) +
  theme_minimal()

ggsave(plot_path, plot, width = 6, height = 5, dpi = 800)
