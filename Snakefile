# Barebone Snakemake Pipeline
#   for paired-end Illumina sequencing
# Snakemake main rule
# Ye Yuan (yeyu@umich.edu)
# August, 2022


from workflow import bb_core
from os import path
configfile: "config.yaml"
ends = bb_core.get_paired_ends(verbose=True)
ext = bb_core.get_ext(verbose=True)
samples = bb_core.list_samples()

localrules: bb_all, multiqc, bb_notrim

include: "workflow/multiqc.smk"
include: "workflow/custom.smk"


# ------ INPUT RULES ------
rule bb_all:
    # all rule should include all TERMINAL steps (iter up till fastq)
    #   bb workflow provides only fastqc by default
    # add in your custom pipeline in workflow/custom.smk and add the terminal step here
    input:
        expand(
            path.join("processed", "{sample}", config["trim_type"], "fastqc", "{end}_fastqc.html"),
            sample=samples, end=ends
        )


rule multiqc:
    # this is a separate input rule, using multiqc to gather logs for all samples
    # REQUIRES modification to the multiqc_sample rule workflow/multiqc.smk to work properly
    input:
        expand(
            path.join("processed", "{sample}", config["trim_type"], "multiqc_report.html"),
            sample=samples
        )



# ------ MINIMUM WORKFLOW RULES ------
rule bb_fastqc:
    # Minimum workflow: fastqc
    input:
        path.join("processed", "{sample}", config["trim_type"], "fastq", "{end}" + ext)
    output:
        path.join("processed","{sample}",config["trim_type"],"fastqc","{end}_fastqc.html")
    conda:
        "env.yaml"
    threads:
        config["resources"]["threads"]["fastqc"]
    resources:
        mem_mb=16000,
        time="1:00:00"
    params:
        out_dir = lambda wildcards, output: path.dirname(output[0])
    shell:
        "fastqc -t {threads} -o {params.out_dir} {input}"


rule bb_trim:
    # Minimum workflow: paired-end trimming
    # IMPORTANT: You may need to modify cutadapt parameters here for your needs.
    #   Here we show settings with Nextera library
    input:
        unpack(lambda wildcards: bb_core.sample_raw_data(wildcards.sample))
    output:
        R1=path.join("processed", "{sample}", "trimmed", "fastq", ends[0] + ext),
        R2=path.join("processed", "{sample}", "trimmed", "fastq", ends[1] + ext)
    conda:
        "env.yaml"
    threads:
        config["resources"]["threads"]["trim"]
    resources:
        mem_mb=16000,
        time="3:00:00"
    params:
        out_dir = lambda wildcards, output: path.dirname(output[0]),
        R1_adapter = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC",
        R2_adapter = "CTGTCTCTTATACACATCTGACGCTGCCGACGA",
        min_length = "20",
        min_quality = "30"
    log:
        path.join("processed", "{sample}", "trimmed", "logs", "cutadapt.log")
    shell:
        "cutadapt -j {threads} -a {params.R1_adapter} -A {params.R2_adapter} " +
        "-m {params.min_length} -q {params.min_quality} " +
        "-o {output.R1} -p {output.R2} {input.R1} {input.R2} > {log}"

rule bb_notrim:
    # Minimum workflow: no trimming, use original FASTQ
    input:
        path.join("raw_data", "{sample}{end}" + ext)
    output:
        path.join("processed","{sample}", "original", "fastq", "{end}" + ext)
    params:
        abs_in = lambda wildcards, input: path.abspath(input[0]),
        abs_out = lambda wildcards, output: path.abspath(output[0])
    shell:
        "ln -s {params.abs_in} {params.abs_out}"
