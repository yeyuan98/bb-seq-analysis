# Barebone Snakemake Pipeline
#   for paired-end Illumina sequencing
# Snakemake summarize MultiQC rule
# Ye Yuan (yeyu@umich.edu)
# August, 2022

from os import path

localrules: multiqc_sample

rule multiqc_sample:
    # This only includes FASTQC and the logs folder for each sample.
    # You should configure this rule further if you want to gather information from other locations.
    #   Note - because we are specifying folders but not files, we cannot use input directive.
    params:
        out_dir = path.join("processed","{sample}",config["trim_type"]),
        fastqc = path.join("processed", "{sample}", config["trim_type"], "fastqc"),
        logs = path.join("processed", "{sample}", config["trim_type"], "logs"),
        alignment = path.join("processed", "{sample}", config["trim_type"], "alignment"),
        peak_calling = path.join("processed", "{sample}", config["trim_type"], "peak_calling")
    output:
        path.join("processed","{sample}",config["trim_type"],"multiqc_report.html")
    conda:
        "../env.yaml"
    shell:
        "multiqc --no-data-dir -o {params.out_dir} {params.fastqc} {params.logs} {params.alignment} {params.peak_calling}"
