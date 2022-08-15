# Barebone Snakemake Pipeline
#   for ATAC-seq with paired-end sequencing
# Peak calling with MACS2
# Ye Yuan (yeyu@umich.edu)
# August, 2022

from os import path

# Note:
#   Both `bedtools bamtobed` and `macs2 callpeak` do not have multithreading.

rule BAMtoBED:
    input:
        path.join("processed", "{sample}", config["trim_type"], "alignment", "filtered.bam")
    output:
        path.join("processed", "{sample}", config["trim_type"], "peak_calling", "filtered.bed")
    conda:
        "../env.yaml"
    threads:
        config["resources"]["threads"]["peak_calling"]
    resources:
        mem_mb=config["resources"]["mem_mb"]["peak_calling"],
        time="3:00:00"
    log:
        path.join("processed","{sample}",config["trim_type"],"logs","peak_calling.log")
    shell:
        "bedtools bamtobed -i {input} 1> {output} 2>> {log}"


rule MACS2:
    input:
        path.join("processed","{sample}",config["trim_type"], "peak_calling", "filtered.bed")
    output:
        path.join("processed","{sample}",config["trim_type"], "peak_calling", "NA_peaks.xls")
    conda:
        "../env.yaml"
    threads:
        config["resources"]["threads"]["peak_calling"]
    resources:
        mem_mb=config["resources"]["mem_mb"]["peak_calling"],
        time="3:00:00"
    log:
        path.join("processed","{sample}",config["trim_type"],"logs","peak_calling.log")
    params:
        eff_genome=config["peak_calling"]["effective_genome_size"],
        shift=config["peak_calling"]["shift"],
        ext=config["peak_calling"]["extension_size"],
        out_dir=lambda wildcards, output: path.dirname(output[0])
    shell:
        """
            macs2 callpeak -t {input} -g {params.eff_genome} --nomodel --shift {params.shift} --extsize {params.ext} \
                -f BED --call-summits --keep-dup all --bdg --outdir {params.out_dir} >> {log} 2>&1
        """