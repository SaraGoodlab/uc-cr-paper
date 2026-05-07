#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
})

#Parse CLI arguments
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

# Load the merged count matrix and enforce the expected gene/sample layout.
counts_tbl <- read_tsv(counts_path, show_col_types = FALSE)
if (!"gene_id" %in% names(counts_tbl)) {
  stop("Counts table must contain a gene_id column", call. = FALSE)
}

# Validate that metadata has the columns needed to label PCA
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

# Apply the QC exclusion list before subsetting the count matrix for PCA.
if (file.exists(excluded_path)) {
  to_exclude <- sort(unique(trimws(read_lines(excluded_path))))
} else {
  to_exclude <- character()
}

to_exclude <- to_exclude[to_exclude != ""]
gene_ids <- counts_tbl$gene_id
count_tbl_only <- counts_tbl |> select(-gene_id)
sample_ids <- colnames(count_tbl_only)
analysis_sample_ids <- setdiff(sample_ids, to_exclude)
missing_metadata <- setdiff(analysis_sample_ids, metadata$sample_id)
if (length(missing_metadata) > 0) {
  stop(
    "Counts table columns missing from metadata after exclusions: ",
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

metadata <- metadata |> filter(sample_id %in% analysis_sample_ids)
if (nrow(metadata) < 3) {
  stop("Need at least three samples for PCA", call. = FALSE)
}

# Build a DESeq2 object for normalization/VST
count_mat <- as.matrix(count_tbl_only[, metadata$sample_id, drop = FALSE])
rownames(count_mat) <- gene_ids
count_mat <- round(count_mat)
col_data <- metadata |> column_to_rownames("sample_id")

dds <- DESeqDataSetFromMatrix(countData = count_mat, 
                              colData = col_data, design = ~1)
dds <- estimateSizeFactors(dds)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
vsd_mat <- assay(vsd)
pca <- prcomp(t(vsd_mat), center = TRUE, scale. = FALSE)

# Export both the plotted coordinates and a figure.
pca_df <- as.data.frame(pca$x[, 1:2, drop = FALSE]) |>
  mutate(
    sample_id = rownames(pca$x),
    tissue = col_data[sample_id, "tissue"],
    treatment = col_data[sample_id, "treatment"]
  ) |>
  arrange(tissue, treatment, sample_id)

write_csv(pca_df, output_path)

percent_var <- pca$sdev^2 / sum(pca$sdev^2)
percent_label <- function(x) sprintf("%.1f%%", x * 100)

plot <- ggplot(pca_df,
               aes(x = PC1, y = PC2, color = treatment, shape = tissue)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "RNA-seq PCA",
    x = paste0("PC1 (", percent_label(percent_var[1]), ")"),
    y = paste0("PC2 (", percent_label(percent_var[2]), ")")
  ) +
  theme_minimal()

ggsave(plot_path, plot, width = 6, height = 5, dpi = 800)
