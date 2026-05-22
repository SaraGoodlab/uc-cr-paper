"""
Post-processing of nf-core/rnaseq output and DESeq2 differential expression
"""

# Configuration
RNASEQ_CONFIG = config["rnaseq"]
RNASEQ_DESEQ2 = RNASEQ_CONFIG["deseq2"]
RNASEQ_CONTRASTS = {contrast["id"]: contrast for contrast in RNASEQ_CONFIG["contrasts"]}
RNASEQ_CONTRAST_IDS = list(RNASEQ_CONTRASTS)
RNASEQ_SAMPLES = RNASEQ_CONFIG["samples"]

# RNA-seq input/output locations
RNASEQ_METADATA = METADATA_DIR / "rnaseq.tsv"
RNASEQ_PROCESSED = PROCESSED_DATA_DIR / "rnaseq"
NFCORE_OUTPUT_DIR = RNASEQ_PROCESSED / "nfcore"
# nf-core/rnaseq places STAR-RSEM outputs under star_rsem/ by default.
NFCORE_RSEM_DIR = NFCORE_OUTPUT_DIR / "star_rsem"
RNASEQ_COUNTS_DIR = RNASEQ_PROCESSED / "rsem_counts"
RNASEQ_RESULTS_DIR = RESULTS_DIR / "rnaseq"
RNASEQ_QC_DIR = RNASEQ_RESULTS_DIR / "quality_metrics"
RNASEQ_DESEQ_DIR = RNASEQ_RESULTS_DIR / "deseq2_results"
RNASEQ_DEG_DIR = RNASEQ_RESULTS_DIR / "deg_lists"
RNASEQ_PCA_DIR = FIGURES_DIR / "rnaseq" / "pca_plots"
RNASEQ_TABLE_DIR = TABLES_DIR / "rnaseq"
RNASEQ_SIZE_FACTORS = RNASEQ_QC_DIR / "library_size_factors.csv"
RNASEQ_EXCLUDED_SAMPLES = RNASEQ_QC_DIR / "excluded_samples.txt"

RNASEQ_RSEM_FILES = expand(
    NFCORE_RSEM_DIR / "{sample}.genes.results",
    sample=RNASEQ_SAMPLES,
)

RNASEQ_RESULT_FILES = expand(
    RNASEQ_DESEQ_DIR / "{contrast}_results.csv",
    contrast=RNASEQ_CONTRAST_IDS,
)


rule collect_rsem_counts:
    """
    Merge per-sample nf-core RSEM gene-level files into one count matrix via tximport.
    """
    input:
        rsem_files = RNASEQ_RSEM_FILES,
        metadata = RNASEQ_METADATA
    output:
        counts = RNASEQ_COUNTS_DIR / "gene_counts.tsv"
    params:
        rsem_dir = NFCORE_RSEM_DIR,
        script = CODE_DIR / "01_rnaseq" / "differential_expression" / "tximport_rsem_counts.R"
    log:
        LOG_DIR / "01_rnaseq" / "collect_rsem_counts.log"
    container:
        R_ANALYSIS_CONTAINER
    shell:
        """
        Rscript {params.script} \
            --rsem_dir {params.rsem_dir} \
            --metadata {input.metadata} \
            --output {output.counts} \
            > {log} 2>&1
        """


rule quality_control:
    """
    Estimate size factors and record any samples flagged for exclusion.
    """
    input:
        rsem_counts = RNASEQ_COUNTS_DIR / "gene_counts.tsv",
        metadata = RNASEQ_METADATA
    output:
        qc_summary = RNASEQ_QC_DIR / "qc_summary.csv",
        size_factors = RNASEQ_SIZE_FACTORS,
        excluded_samples = RNASEQ_EXCLUDED_SAMPLES
    params:
        script = CODE_DIR / "01_rnaseq" / "qc" / "quality_control.R",
        size_factor_threshold = RNASEQ_DESEQ2["size_factor_threshold"]
    log:
        LOG_DIR / "01_rnaseq" / "quality_control.log"
    container:
        R_ANALYSIS_CONTAINER
    shell:
        """
        Rscript {params.script} \
            --counts {input.rsem_counts} \
            --metadata {input.metadata} \
            --output {output.qc_summary} \
            --size_factors {output.size_factors} \
            --excluded {output.excluded_samples} \
            --threshold {params.size_factor_threshold} \
            > {log} 2>&1
        """


rule pca_analysis:
    """
    Produce a global PCA view after applying the QC-driven exclusion list.
    """
    input:
        rsem_counts = RNASEQ_COUNTS_DIR / "gene_counts.tsv",
        metadata = RNASEQ_METADATA,
        excluded = RNASEQ_EXCLUDED_SAMPLES
    output:
        pca_plot = RNASEQ_PCA_DIR / "pca_all_samples.png",
        pca_data = RNASEQ_QC_DIR / "pca_results.csv"
    params:
        script = CODE_DIR / "01_rnaseq" / "qc" / "pca_analysis.R"
    log:
        LOG_DIR / "01_rnaseq" / "pca_analysis.log"
    container:
        R_ANALYSIS_CONTAINER
    shell:
        """
        Rscript {params.script} \
            --counts {input.rsem_counts} \
            --metadata {input.metadata} \
            --excluded {input.excluded} \
            --plot {output.pca_plot} \
            --output {output.pca_data} \
            > {log} 2>&1
        """


rule deseq2_analysis:
    """
    Run one DESeq2 contrast per tissue/comparison pair defined in the config.
    """
    input:
        counts = RNASEQ_COUNTS_DIR / "gene_counts.tsv",
        metadata = RNASEQ_METADATA,
        excluded = RNASEQ_EXCLUDED_SAMPLES
    output:
        results = RNASEQ_DESEQ_DIR / "{contrast}_results.csv",
        degs = RNASEQ_DEG_DIR / "{contrast}_degs.txt"
    params:
        script = CODE_DIR / "01_rnaseq" / "differential_expression" / "deseq2_analysis.R",
        tissue = lambda wildcards: RNASEQ_CONTRASTS[wildcards.contrast]["tissue"],
        comparison = lambda wildcards: RNASEQ_CONTRASTS[wildcards.contrast]["comparison"],
        case = lambda wildcards: RNASEQ_CONTRASTS[wildcards.contrast]["case"],
        control = lambda wildcards: RNASEQ_CONTRASTS[wildcards.contrast]["control"],
        padj_threshold = RNASEQ_DESEQ2["padj_threshold"],
        log2fc_threshold = RNASEQ_DESEQ2["log2fc_threshold"]
    log:
        LOG_DIR / "01_rnaseq" / "deseq2"/ "{contrast}.log"
    container:
        R_ANALYSIS_CONTAINER
    shell:
        """
        Rscript {params.script} \
            --counts {input.counts} \
            --metadata {input.metadata} \
            --excluded {input.excluded} \
            --tissue {params.tissue} \
            --comparison {params.comparison} \
            --case {params.case} \
            --control {params.control} \
            --padj {params.padj_threshold} \
            --log2fc {params.log2fc_threshold} \
            --results {output.results} \
            --degs {output.degs} \
            > {log} 2>&1
        """


rule summarize_rnaseq_results:
    """
    Summarize all configured DESeq2 contrasts with the same thresholds used upstream.
    """
    input:
        RNASEQ_RESULT_FILES
    output:
        summary = RNASEQ_DESEQ_DIR / "summary_degs.txt",
        table = RNASEQ_TABLE_DIR / "differential_expression.rnaseq.csv"
    params:
        script = CODE_DIR / "01_rnaseq" / "differential_expression" / "extract_degs.R",
        padj_threshold = RNASEQ_DESEQ2["padj_threshold"],
        log2fc_threshold = RNASEQ_DESEQ2["log2fc_threshold"]
    log:
        LOG_DIR / "01_rnaseq" / "summarize_rnaseq_results.log"
    container:
        R_ANALYSIS_CONTAINER
    shell:
        """
        Rscript {params.script} \
            --input_dir {RNASEQ_DESEQ_DIR} \
            --padj {params.padj_threshold} \
            --log2fc {params.log2fc_threshold} \
            --summary {output.summary} \
            --table {output.table} \
            > {log} 2>&1
        """
