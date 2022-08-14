# Barebone Snakemake Pipeline
#   for ATAC-seq with paired-end sequencing
# Bowtie2 alignment
# Ye Yuan (yeyu@umich.edu)
# August, 2022

from os import path
import bb_core

ends = bb_core.get_paired_ends()
ext = bb_core.get_ext()

rule bowtie2_align:
    input:
        fastq_R1=path.join("processed", "{sample}", config["trim_type"], "fastq", ends[0]+ext),
        fastq_R2=path.join("processed","{sample}",config["trim_type"],"fastq",ends[1]+ext),
        # index_bait is used to trigger index existence check by snakemake. Not used in command directly.
        index_bait=path.join("processed", "index", config["genome_name"], config["genome_name"] + ".1.bt2")
    params:
        index_base=path.join("processed", "index", config["genome_name"], config["genome_name"])
    output:
        path.join("processed", "{sample}", config["trim_type"], "alignment", "out.bam")
    conda:
        "../env.yaml"
    threads:
        config["resources"]["threads"]["align"]
    resources:
        mem_mb=config["resources"]["mem_mb"]["align"],
        time="3:00:00"
    log:
        path.join("processed", "{sample}", config["trim_type"], "logs", "align.log")
    shell:
        "bowtie2 -x {params.index_base} --very-sensitive --dovetail -X 1000 -p {threads} "+
        "-1 {input.fastq_R1} -2 {input.fastq_R2} -b {output} >> {log}"
