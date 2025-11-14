"""
Data retrieval rules for downloading sequencing data from ENA
"""

rule download_ena_rnaseq:
    """
    Download RNA-seq FASTQ files from ENA.
    """
    output:
        colon_dir = directory(DATA_DIR / "raw" / "rnaseq" / "colon"),
        spleen_dir = directory(DATA_DIR / "raw" / "rnaseq" / "spleen")
    params:
        script = CODE_DIR / "00_data_retrieval" / "download_ena_rnaseq.sh"
    shell:
        """
        bash {params.script}
        """

rule download_ena_microbiome:
    """
    Download 16S rRNA microbiome FASTQ files from ENA.
    """
    output:
        output_dir = directory(DATA_DIR / "raw" / "microbiome")
    params:
        script = CODE_DIR / "00_data_retrieval" / "download_ena_microbiome.sh"
    shell:
        """
        bash {params.script}
        """

