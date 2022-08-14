# Barebone Snakemake Pipeline
#   for ATAC-seq with paired-end sequencing
#   Genome download and index building
# Ye Yuan (yeyu@umich.edu)
# August, 2022

from os import path


localrules: genome_download_filter

seqkit_url = "https://github.com/shenwei356/seqkit/releases/download/v2.3.0/seqkit_linux_amd64.tar.gz"
seqkit_bin = "seqkit"
genome_url = "http://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz"
genome_contigs_removal = ["chrUn", "random"]
genome_name = "dm6"

def get_file_from_url(url):
    return url[url.rfind("/")+1:]

rule genome_download_filter:
    # Downloads genome for index building
    # ... and then filter by contig
    # TODO: This rule downloads external FASTQ processing kit to perform filtering.
    output:
        temp(get_file_from_url(genome_url))
    params:
        filter_str="-p " + " -p ".join(genome_contigs_removal),
        download_url=genome_url,
        out_path=lambda wildcards, output: output[0],
        seqkit_url=seqkit_url,
        seqkit_zip=get_file_from_url(seqkit_url),
        seqkit_bin=seqkit_bin
    log:
        path.join("processed", "genome_download.log")
    shell:
        """
            echo "---------- GETTING SEQKIT ----------" >> {log}
            wget -a {log} {params.seqkit_url}
            tar xzvf {params.seqkit_zip}
            chmod +x {params.seqkit_bin}
            echo "---------- WGET GENOME ----------" >> {log}
            wget -a {log} {params.download_url}
            echo "---------- ORIGINAL GENOME FASTQ ENTRIES ----------" >> {log}
            zcat {output} | grep -P '^>' - >> {log}
            echo "---------- FILTERING GENOME FASTQ ----------" >> {log}
            zcat {output} | {params.seqkit_bin} grep -r -v {params.filter_str} | gzip -c - > {params.out_path}
            echo "---------- FILTERED GENOME FASTQ ENTRIES ----------" >> {log}
            zcat {output} | grep -P '^>' - >> {log}
            echo "---------- CLEANING UP ----------" >> {log}
            rm {params.seqkit_bin}
        """


rule bowtie2_index:
    # Creates Bowtie2 index based on genome information provided
    input:
        get_file_from_url(genome_url)
    output:
        directory(path.join("processed", "index", genome_name))
    conda:
        "../env.yaml"
    shell:
        """
            
        """