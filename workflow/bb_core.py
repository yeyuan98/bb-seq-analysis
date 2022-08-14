# Barebone Snakemake Pipeline
#   for paired-end Illumina sequencing
# Ye Yuan (yeyu@umich.edu)
# August, 2022

"""
    The following file hierarchy shall be observed by all pipelines
    root
        raw_data
            {sample1}.{read1}.fastq.gz
            {sample1}.{read2}.fastq.gz
            ...... (more samples)
        processed
            {sample1}
                {fastq_type=original}
                    # If original fastq is the intended input, structure here shall be the same as trimmed below
                {fastq_type=trimmed}
                    {step=fastq}
                        # This folder shall store trimmed fastq files
                    {step=fastqc}
                        # Fastqc of the trimmed fastq files
                    logs
                        # Logs from different processing steps shall go here
                    ...... (more processing steps)
                    multiqc.html
                        # multiqc report
            ...... (more samples)
        workflow
            trim.smk
            fastqc.smk
            multiqc.smk
            ...... (more processing steps)
        Snakefile
        bb_core.py
"""


# ------------------ RAW_DATA LEVEL FUNCTIONS ---------------------

def list_raw_data():
    """
        Lists raw_data directory content
    :return: list of paths for files in raw_data folder
    """
    from os import listdir
    return [f for f in listdir("raw_data") if not f.startswith('.')]  # ignore hidden files


def get_paired_ends(verbose=False):
    """
        Returns identified paired-end read suffix
        Method: search for R1[^.]* and R2[^.]* and removes extension if exists
    :param verbose: optional, whether to print message
    :return: list of length 2 [{read1}, {read2}]. Raise error if any of R1/R2 returns non-unique search.
    """
    import re
    # Compile search patterns
    p1 = re.compile("R1[^.]*")
    p2 = re.compile("R2[^.]*")
    data_files = list_raw_data()
    # Perform searches
    r1 = [p1.search(f) for f in data_files]
    r1 = [r[0] for r in r1 if r is not None]
    r2 = [p2.search(f) for f in data_files]
    r2 = [r[0] for r in r2 if r is not None]
    # Verify consistency
    if len(set(r1)) != 1 or len(set(r2)) != 1:
        raise ValueError("BB_CORE: Inconsistent read suffix. Please verify raw_data files.")
    results = [r1[0], r2[0]]
    if verbose:
        print("BB_CORE: Paired-end read suffix are " + str(results))
    return results


def get_ext(verbose=False):
    """
        Return extension of the raw_data files.
        Supported special cases:
            .fastq.gz, .fq.gz
    :param verbose: optional, whether to print message
    :return: string of the raw_data file extension. Raise error if different extensions are found.
    """
    from os.path import splitext

    def splitext_(path):
        for ext in ['.fastq.gz', '.fq.gz']:
            if path.endswith(ext):
                return ext
        return splitext(path)[1]
    data_files = list_raw_data()
    exts = [splitext_(f) for f in data_files]
    # Verify consistency
    if len(set(exts)) != 1:
        raise ValueError("BB_CORE: Inconsistent extensions. Please verify raw_data files.")
    result = exts[0]
    if verbose:
        print("BB_CORE: Paired-end read extension is " + str(result))
    return result


def list_samples():
    """
        Lists all sample strings in the raw_data folder
        Method: identifying paired ends with get_paired_ends(),
            ... and then trim off file names from start of paired ends.
            ... finally get unique names.
    :return: list of sample name strings
    """
    import re
    # Get raw file names
    results = list_raw_data()
    print("BB_CORE: File count = " + str(len(results)))
    # Get paired ends and find the index position for each file
    ends = get_paired_ends()
    p = re.compile("(?:" + ends[0] + "|" + ends[1] + ")")
    idx = [p.search(t).start() for t in results]
    # Subset file name
    results = [r[:i] for r, i in zip(results,idx)]
    # Get unique values
    results = list(set(results))
    return results


# ------------------ SAMPLE LEVEL FUNCTIONS ---------------------

def sample_raw_data(sample):
    """
        Given sample get paths to the paired-end read raw data.
    :param sample: sample name.
    :return: dict of keys R1 / R2, paths to the paired-end raw_data files.
    """
    from os.path import join
    ends = get_paired_ends()
    ext = get_ext()
    return {'R1': join("raw_data", sample + ends[0] + ext), 'R2': join("raw_data", sample + ends[1] + ext)}
