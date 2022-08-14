# Barebone Snakemake Pipeline
#   for ATAC-seq with paired-end sequencing
# Post-alignment operations
#   sort/index of alignment, idxstats,
#   ... & filtering of alignment (samtools)
#   ... & PCR duplicate removal
# Ye Yuan (yeyu@umich.edu)
# August, 2022

from os import path

rule summarizeThenFilter_postalign:
    # Note - Picard MarkDuplicates is not multi-threaded; therefore, this filter step is separated in another rule.
    input:
        path.join("processed", "{sample}", config["trim_type"], "alignment", "out.bam")
    output:
        sorted=path.join("processed", "{sample}", config["trim_type"], "alignment", "sorted.bam"),
        sorted_idx=path.join("processed", "{sample}", config["trim_type"], "alignment", "sorted.bam.bai"),
        idxstat=path.join("processed", "{sample}", config["trim_type"], "alignment", "idxstat.txt"),
        filtered=path.join("processed", "{sample}", config["trim_type"], "alignment", "filtered.postpcr.bam")
    conda:
        "../env.yaml"
    threads:
        config["resources"]["threads"]["summarize_postalign"]
    resources:
        mem_mb=config["resources"]["mem_mb"]["summarize_postalign"],
        time="3:00:00"
    log:
        path.join("processed", "{sample}", config["trim_type"], "logs", "summarize_postalign.log")
    params:
        flag=config["post_align_filters"]["flag"],
        min_quality=config["post_align_filters"]["min_quality"],
        mito=config["post_align_filters"]["mito"]
    shell:
        '''
            now=$(date +"%r")
            echo "$now Starting sort..." >> {log}
            samtools sort -@ {threads} {input} 1> {output.sorted} 2>> {log}
            now=$(date +"%r")
            echo "$now Starting index..." >> {log}
            samtools index -@ {threads} {output.sorted} {output.sorted_idx} >>{log} 2>&1
            now=$(date +"%r")
            echo "$now Getting idxstats of alignment..." >> {log}
            samtools idxstats -@ {threads} {output.sorted} 1> {output.idxstat} 2>> {log}
            now=$(date +"%r")
            echo "$now Filtering by flag and quality..." >> {log}
            samtools view -@ {threads} -b -f {params.flag} -q {params.min_quality} {output.sorted} 1> {wildcards.sample}.intermediate.bam 2>> {log}
            now=$(date +"%r")
            echo "$now Filtering mito reads..." >> {log}
            samtools view {wildcards.sample}.intermediate.bam | grep -Pv "\t{params.mito}\t" - | samtools view -b -o {output.filtered} - >>{log} 2>&1
            rm {wildcards.sample}.intermediate.bam
        '''

rule pcrDuplicateFilter_postalign:
    input:
        path.join("processed", "{sample}", config["trim_type"], "alignment", "filtered.postpcr.bam")
    output:
        idxstat=path.join("processed", "{sample}", config["trim_type"], "alignment", "idxstat.filtered.txt"),
        metrics=path.join("processed", "{sample}", config["trim_type"], "alignment", "picard.metrics.txt"),
        filtered=path.join("processed", "{sample}", config["trim_type"], "alignment", "filtered.bam")
    conda:
        "../env.yaml"
    threads:
        3
    resources:
        mem_mb=20000,
        time="5:00:00"
    log:
        path.join("processed", "{sample}", config["trim_type"], "logs", "summarize_postalign.log")
    shell:
        '''
            now=$(date +"%r")
            echo "$now Removing PCR duplicates..." >> {log}
            picard MarkDuplicates --INPUT {input} --METRICS_FILE {output.metrics} --OUTPUT {output.filtered} --REMOVE_DUPLICATES true
            now=$(date +"%r")
            echo "$now Getting idxstats of filtered alignment..." >> {log}
            samtools idxstats -@ {threads} {output.filtered} 1> {output.idxstat} 2>> {log}
            now=$(date +"%r")
            echo "$now DONE!" >> {log}
        '''
