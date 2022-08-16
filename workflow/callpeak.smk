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
        path.join("processed","{sample}",config["trim_type"], "peak_calling", "NA_peaks.xls"),
        path.join("processed","{sample}",config["trim_type"],"peak_calling","NA_treat_pileup.bdg")
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


rule bdg2bw:
    # Converts bdg to bw for visualization. clips off out-of-bound reads
    # Uses Kent tools.  TODO: PLATFORM SPECIFIC.
    # Needs chromsome sizes information. TODO: [DEV] GET CHROM.SIZES DURING INDEX BUILDING.
    input:
        path.join("processed","{sample}",config["trim_type"],"peak_calling","NA_treat_pileup.bdg")
    output:
        path.join("processed","{sample}",config["trim_type"],"peak_calling","NA_treat_pileup.bw")
    threads:
        1
    resources:
        mem_mb=7000,
        time="1:00:00"
    params:
        bedclip_url="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedClip",
        bdg2bw_url="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig",
        sizes_url="http://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.chrom.sizes",
        sizes_file="./dm6.chrom.sizes"
    shell:
        """
            wget -q -nc {params.bedclip_url}
            wget -q -nc {params.bdg2bw_url}
            wget -q -nc {params.sizes_url}
            chmod +x ./bedClip
            chmod +x ./bedGraphToBigWig
            ./bedClip {input} {params.sizes_file} {wildcards.sample}.bdg
            ./bedGraphToBigWig {wildcards.sample}.bdg {params.sizes_file} {output}
            rm {wildcards.sample}.bdg
        """
