"""
BAM Preprocessing Rules
========================
Input : data/bam/{sample}.bam   (raw aligned reads)
Output: results/bam/final/{sample}.filtered.bam(.bai)

Pipeline
--------
  sort_bam          – coordinate-sort the raw input BAM
  mark_duplicates   – Picard MarkDuplicates, REMOVE_DUPLICATES=true
  filter_bam        – samtools: MAPQ filter, flag filter, chrM removal
  remove_blacklist  – bedtools intersect -v to strip blacklist regions
                      + re-sort + index the final cleaned BAM

FLAG arithmetic (-F 1804):
  4    read unmapped
  8    mate unmapped
  256  not primary alignment
  512  read fails platform/vendor quality checks
  1024 read is PCR or optical duplicate
  Sum  = 1804
"""

_RESULTS   = config["results_dir"]
_BAM_DIR   = config.get("bam_dir", "data/bam")
_BLACKLIST = config.get("blacklist", "resources/hg38-blacklist.v2.bed")
_MIN_MAPQ  = config.get("min_mapq", 30)
_THREADS   = config.get("threads", 8)
_MITO      = config.get("mito_chr", "chrM")
_F_EX      = config.get("filter_flags_exclude", 1804)
_F_IN      = config.get("filter_flags_include", 2)
_PE        = config.get("paired_end", True)


# ── 1. Sort ────────────────────────────────────────────────────────────────────
rule sort_bam:
    """Coordinate-sort the raw input BAM (skip if already sorted)."""
    input:
        bam = _BAM_DIR + "/{sample}.bam",
    output:
        bam = temp(_RESULTS + "/bam/sorted/{sample}.sorted.bam"),
    log:
        _RESULTS + "/logs/sort_bam/{sample}.log",
    threads: _THREADS
    resources:
        mem_mb = config.get("memory_mb", 8000),
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        "samtools sort -@ {threads} -o {output.bam} {input.bam} 2>{log}"


# ── 2. Remove PCR duplicates ───────────────────────────────────────────────────
rule mark_duplicates:
    """Remove PCR / optical duplicates with Picard MarkDuplicates."""
    input:
        bam = _RESULTS + "/bam/sorted/{sample}.sorted.bam",
    output:
        bam     = temp(_RESULTS + "/bam/dedup/{sample}.dedup.bam"),
        metrics = _RESULTS + "/qc/picard_metrics/{sample}.dup_metrics.txt",
    log:
        _RESULTS + "/logs/mark_duplicates/{sample}.log",
    threads: _THREADS
    resources:
        mem_mb = config.get("memory_mb", 8000),
    conda:
        "../../envs/preprocessing.yaml"
    params:
        java_opts = "-Xmx8g",
    shell:
        """
        picard {params.java_opts} MarkDuplicates \
            INPUT={input.bam} \
            OUTPUT={output.bam} \
            METRICS_FILE={output.metrics} \
            REMOVE_DUPLICATES=true \
            ASSUME_SORTED=true \
            VALIDATION_STRINGENCY=LENIENT \
            2>{log}
        samtools index {output.bam} 2>>{log}
        """


# ── 3. Quality & chrM filter ───────────────────────────────────────────────────
rule filter_bam:
    """
    Apply samtools quality filters then remove mitochondrial reads.

    For paired-end data:
      -F {flag_exclude}  exclude unmapped/non-primary/dup/supplementary etc.
      -f {flag_include}  require proper pair
      -q {min_mapq}      minimum MAPQ

    chrM removal uses 'samtools idxstats' to enumerate all chr names,
    greps out the mitochondrial one, and extracts only those chromosomes.
    """
    input:
        bam = _RESULTS + "/bam/dedup/{sample}.dedup.bam",
    output:
        bam = temp(_RESULTS + "/bam/filtered/{sample}.nomito.bam"),
    log:
        _RESULTS + "/logs/filter_bam/{sample}.log",
    threads: _THREADS
    resources:
        mem_mb = config.get("memory_mb", 8000),
    conda:
        "../../envs/preprocessing.yaml"
    params:
        flag_exclude = _F_EX,
        flag_include = ("-f " + str(_F_IN)) if _PE else "",
        mapq         = _MIN_MAPQ,
        mito         = _MITO,
    shell:
        """
        # Step 1: MAPQ + flag filter → temp BAM
        samtools view -@ {threads} -b \
            -F {params.flag_exclude} \
            {params.flag_include} \
            -q {params.mapq} \
            {input.bam} \
            -o {_RESULTS}/bam/filtered/{wildcards.sample}.qcpassed.bam \
            2>{log}
        samtools index {_RESULTS}/bam/filtered/{wildcards.sample}.qcpassed.bam 2>>{log}

        # Step 2: remove mitochondrial chromosome reads
        samtools idxstats \
            {_RESULTS}/bam/filtered/{wildcards.sample}.qcpassed.bam \
            | awk '$3>0 {{print $1}}' \
            | grep -v '{params.mito}' \
            | xargs samtools view -@ {threads} -b \
                {_RESULTS}/bam/filtered/{wildcards.sample}.qcpassed.bam \
                -o {output.bam} 2>>{log}

        # Clean up temp
        rm -f {_RESULTS}/bam/filtered/{wildcards.sample}.qcpassed.bam \
              {_RESULTS}/bam/filtered/{wildcards.sample}.qcpassed.bam.bai
        """


# ── 4. Remove ENCODE blacklist regions ────────────────────────────────────────
rule remove_blacklist:
    """
    Exclude reads overlapping ENCODE blacklist regions via bedtools intersect.
    Re-sort and index the final cleaned BAM.

    Blacklist file: config['blacklist']
    Download from: https://github.com/Boyle-Lab/Blacklist/tree/master/lists
    """
    input:
        bam       = _RESULTS + "/bam/filtered/{sample}.nomito.bam",
        blacklist = _BLACKLIST,
    output:
        bam = _RESULTS + "/bam/final/{sample}.filtered.bam",
        bai = _RESULTS + "/bam/final/{sample}.filtered.bam.bai",
    log:
        _RESULTS + "/logs/remove_blacklist/{sample}.log",
    threads: _THREADS
    resources:
        mem_mb = config.get("memory_mb", 8000),
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        """
        bedtools intersect -v -abam {input.bam} -b {input.blacklist} \
            2>{log} \
            | samtools sort -@ {threads} -o {output.bam} 2>>{log}
        samtools index {output.bam} 2>>{log}
        """
