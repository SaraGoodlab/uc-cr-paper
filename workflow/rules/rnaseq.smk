"""
RNA-seq analysis rules: pot-prcoessing of nfcore-rnaseq output and DESeq2 differential expression
"""

import re

# Get sample metadata (if available)
TISSUES = config["rnaseq"]["tissues"]
COMPARISONS = config["rnaseq"]["comparisons"]
RNASEQ_METADATA = DATA_DIR / "metadata" / "rnaseq.tsv"
NFCORE_OUTPUT_DIR = DATA_DIR / "processed" / "rnaseq" / "nfcore"
# nf-core/rnaseq places STAR-RSEM outputs under star_rsem/ by default
NFCORE_RSEM_DIR = NFCORE_OUTPUT_DIR / "star_rsem"
COMPARISON_LOOKUP = {c["name"]: c for c in COMPARISONS}
TISSUE_REGEX = "|".join(sorted(re.escape(t) for t in TISSUES))
COMPARISON_REGEX = "|".join(sorted(re.escape(c["name"]) for c in COMPARISONS))


rule collect_rsem_counts:
    """
    Merge nf-core RSEM gene-level results into a single counts matrix via tximport
    """
    input:
        rsem_dir = directory(NFCORE_RSEM_DIR),
        metadata = RNASEQ_METADATA
    output:
        counts = DATA_DIR / "processed" / "rnaseq" / "rsem_counts" / "gene_counts.tsv"
    params:
        rsem_dir = NFCORE_RSEM_DIR,
        script = CODE_DIR / "01_rnaseq" / "differential_expression" / "tximport_rsem_counts.R"
    log:
        "logs/rnaseq/collect_rsem_counts.log"
    conda:
        "envs/deseq2_env.yaml"
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
    Perform quality control checks and calculate size factors
    """
    input:
        rsem_counts = DATA_DIR / "processed" / "rnaseq" / "rsem_counts" / "gene_counts.tsv",
        metadata = RNASEQ_METADATA
    output:
        qc_summary = RESULTS_DIR / "rnaseq" / "quality_metrics" / "qc_summary.csv",
        size_factors = DATA_DIR / "metadata" / "library_size_factors.csv",
        excluded_samples = DATA_DIR / "metadata" / "excluded_samples.txt"
    params:
        script = CODE_DIR / "01_rnaseq" / "qc" / "quality_control.R",
        size_factor_threshold = config["rnaseq"]["deseq2"]["size_factor_threshold"]
    conda:
        "envs/r_analysis_env.yaml"
    shell:
        """
        Rscript {params.script} \
            --counts {input.rsem_counts} \
            --metadata {input.metadata} \
            --output {output.qc_summary} \
            --size_factors {output.size_factors} \
            --excluded {output.excluded_samples} \
            --threshold {params.size_factor_threshold}
        """

rule pca_analysis:
    """
    Perform PCA analysis on variance-stabilized counts
    """
    input:
        rsem_counts = DATA_DIR / "processed" / "rnaseq" / "rsem_counts" / "gene_counts.tsv",
        metadata = RNASEQ_METADATA,
        excluded = DATA_DIR / "metadata" / "excluded_samples.txt"
    output:
        pca_plot = FIGURES_DIR / "rnaseq" / "pca_plots" / "pca_all_samples.png",
        pca_data = RESULTS_DIR / "rnaseq" / "quality_metrics" / "pca_results.csv"
    params:
        script = CODE_DIR / "01_rnaseq" / "qc" / "pca_analysis.R"
    conda:
        "envs/r_analysis_env.yaml"
    shell:
        """
        Rscript {params.script} \
            --counts {input.rsem_counts} \
            --metadata {input.metadata} \
            --excluded {input.excluded} \
            --plot {output.pca_plot} \
            --output {output.pca_data}
        """

def get_comparison(name):
    if name not in COMPARISON_LOOKUP:
        raise ValueError(f"Comparison {name} not defined in config")
    return COMPARISON_LOOKUP[name]

rule deseq2_analysis:
    """
    Run DESeq2 differential expression analysis for each tissue and comparison
    """
    input:
        counts = DATA_DIR / "processed" / "rnaseq" / "rsem_counts" / "gene_counts.tsv",
        metadata = RNASEQ_METADATA,
        excluded = DATA_DIR / "metadata" / "excluded_samples.txt"
    output:
        results = RESULTS_DIR / "rnaseq" / "deseq2_results" / "{tissue}_{comparison}_results.csv",
        degs = RESULTS_DIR / "rnaseq" / "deg_lists" / "{tissue}_{comparison}_degs.txt"
    params:
        script = CODE_DIR / "01_rnaseq" / "differential_expression" / "deseq2_analysis.R",
        tissue = "{tissue}",
        comparison = "{comparison}",
        case = lambda wildcards: get_comparison(wildcards.comparison)["case"],
        control = lambda wildcards: get_comparison(wildcards.comparison)["control"],
        padj_threshold = config["rnaseq"]["deseq2"]["padj_threshold"],
        log2fc_threshold = config["rnaseq"]["deseq2"]["log2fc_threshold"]
    conda:
        "envs/deseq2_env.yaml"
    resources:
        mem_mb = config["resources"]["deseq2_mem"]
    wildcard_constraints:
        tissue = TISSUE_REGEX,
        comparison = COMPARISON_REGEX
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
            --results {outilya@ilya-Nitro-AN515-44:~/Projects/test/workflow/rules$ cp ~/Projects/ibd_mouse/uc-cr-paper/workflow/rules/rnaseq.smk .
put.results} \
            --degs {output.degs}
        """

rule summarize_rnaseq_results:
    """
    Create summary of all RNA-seq DEG results
    """
    input:
        expand(
            RESULTS_DIR / "rnaseq" / "deseq2_results" / "{tissue}_{comparison}_results.csv",
            tissue=TISSUES,
            comparison=[c["name"] for c in COMPARISONS],
        )
    output:
        summary = RESULTS_DIR / "rnaseq" / "deseq2_results" / "summary_degs.txt",
        table = TABLES_DIR / "differential_expression.rnaseq.csv"
    params:
        script = CODE_DIR / "01_rnaseq" / "differential_expression" / "extract_degs.R"
    conda:
        "envs/r_analysis_env.yaml"
    shell:
        """
        Rscript {params.script} \
            --input_dir {RESULTS_DIR}/rnaseq/deseq2_results \
            --summary {output.summary} \
            --table {output.supp_table}
        """
