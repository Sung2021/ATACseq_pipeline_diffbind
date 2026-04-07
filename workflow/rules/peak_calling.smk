"""
Peak Calling Rules  (MACS3)
============================
Input : results/bam/final/{sample}.filtered.bam
Output:
  results/peaks/{sample}/{sample}_peaks.narrowPeak   – per-sample narrow peaks
  results/peaks/{sample}/{sample}_peaks.broadPeak    – per-sample broad peaks (optional)
  results/peaks/consensus/{condition}_consensus_peaks.bed – per-condition merged peaks
  results/peaks/consensus/all_samples_consensus_peaks.bed – all-sample merged peaks

Notes
-----
* --format BAMPE  uses actual fragment size from PE reads (recommended for ATAC-seq)
* --keep-dup all  duplicates have already been removed in preprocessing
* Summits are used as reference points for DiffBind re-centering
"""

_RESULTS  = config["results_dir"]
_THREADS  = config.get("threads", 8)
_BLACKLIST = config.get("blacklist", "resources/hg38-blacklist.v2.bed")

_MACS3_FMT  = config["macs3"].get("format",   "BAMPE")
_MACS3_Q    = config["macs3"].get("qvalue",   0.05)
_MACS3_DUP  = config["macs3"].get("keep_dup", "all")
_CALL_SUM   = config["macs3"].get("call_summits", True)
_GSIZE      = config.get("genome_size", "hs")
_MERGE_GAP  = config.get("merge_gap",   100)
_MIN_SAMP   = config.get("consensus_min_samples", 2)


# ── 1. MACS3 narrow-peak calling ──────────────────────────────────────────────
rule macs3_callpeak:
    """
    Call open-chromatin peaks per sample with MACS3.

    Key options for ATAC-seq:
      -f BAMPE   – use paired-end fragment size (no need for --nomodel)
      --keep-dup all – duplication already handled in preprocessing
      --call-summits – identify local summit within each peak (useful for
                       motif analysis and DiffBind re-centering)
    """
    input:
        bam = _RESULTS + "/bam/final/{sample}.filtered.bam",
    output:
        narrowpeak = _RESULTS + "/peaks/{sample}/{sample}_peaks.narrowPeak",
        summits    = _RESULTS + "/peaks/{sample}/{sample}_summits.bed",
        xls        = _RESULTS + "/peaks/{sample}/{sample}_peaks.xls",
    log:
        _RESULTS + "/logs/macs3/{sample}.log",
    threads: 1   # MACS3 is single-threaded
    conda:
        "../../envs/macs3.yaml"
    params:
        outdir      = _RESULTS + "/peaks/{sample}",
        fmt         = _MACS3_FMT,
        qval        = _MACS3_Q,
        keep_dup    = _MACS3_DUP,
        gsize       = _GSIZE,
        call_summits = "--call-summits" if _CALL_SUM else "",
    shell:
        """
        macs3 callpeak \
            -t {input.bam} \
            -f {params.fmt} \
            -n {wildcards.sample} \
            --outdir {params.outdir} \
            -g {params.gsize} \
            --keep-dup {params.keep_dup} \
            -q {params.qval} \
            {params.call_summits} \
            2>{log}
        """


# ── 2. MACS3 broad-peak calling (optional / for annotation) ───────────────────
rule macs3_callpeak_broad:
    """
    Broad-peak mode to capture larger accessible domains.
    Useful for annotation but narrow peaks are preferred for DiffBind.
    This rule is NOT in the default 'all' target – call explicitly if needed:
      snakemake results/peaks/{sample}/{sample}_peaks.broadPeak
    """
    input:
        bam = _RESULTS + "/bam/final/{sample}.filtered.bam",
    output:
        broadpeak = _RESULTS + "/peaks/{sample}/{sample}_peaks.broadPeak",
        gappedpeak = _RESULTS + "/peaks/{sample}/{sample}_peaks.gappedPeak",
    log:
        _RESULTS + "/logs/macs3_broad/{sample}.log",
    conda:
        "../../envs/macs3.yaml"
    params:
        outdir   = _RESULTS + "/peaks/{sample}",
        fmt      = _MACS3_FMT,
        qval     = _MACS3_Q,
        gsize    = _GSIZE,
    shell:
        """
        macs3 callpeak \
            -t {input.bam} \
            -f {params.fmt} \
            -n {wildcards.sample} \
            --outdir {params.outdir} \
            -g {params.gsize} \
            --keep-dup all \
            -q {params.qval} \
            --broad \
            --broad-cutoff {params.qval} \
            2>{log}
        """


# ── 3. Remove blacklist from peaks ────────────────────────────────────────────
rule filter_peaks_blacklist:
    """Strip ENCODE blacklist regions from per-sample narrow peaks."""
    input:
        peaks     = _RESULTS + "/peaks/{sample}/{sample}_peaks.narrowPeak",
        blacklist = _BLACKLIST,
    output:
        _RESULTS + "/peaks/{sample}/{sample}_peaks.filtered.narrowPeak",
    log:
        _RESULTS + "/logs/filter_peaks/{sample}.log",
    conda:
        "../../envs/macs3.yaml"
    shell:
        """
        bedtools intersect -v \
            -a {input.peaks} \
            -b {input.blacklist} \
            > {output} 2>{log}
        """


# ── 4. Per-condition consensus peak set ───────────────────────────────────────
rule consensus_peaks_per_condition:
    """
    Merge peaks from all samples in a condition into a non-redundant
    consensus peak set.

    Strategy:
      1. Concatenate all filtered narrowPeak BEDs for the condition
      2. Sort by coordinate
      3. bedtools merge -d {merge_gap}  (merge peaks within merge_gap bp)
      4. Keep only standard chromosomes (grep -v _alt, random, etc.)

    The resulting BED is used by DiffBind as the analysis peak universe.
    """
    input:
        peaks = lambda wc: expand(
            _RESULTS + "/peaks/{sample}/{sample}_peaks.filtered.narrowPeak",
            sample=CONDITION_SAMPLES[wc.condition],
        ),
    output:
        _RESULTS + "/peaks/consensus/{condition}_consensus_peaks.bed",
    log:
        _RESULTS + "/logs/consensus_peaks/{condition}.log",
    conda:
        "../../envs/macs3.yaml"
    params:
        gap = _MERGE_GAP,
    shell:
        """
        cat {input.peaks} \
            | awk 'OFS="\\t" {{print $1,$2,$3}}' \
            | grep -vE '_(random|alt|Un)|chrM|chrEBV' \
            | sort -k1,1 -k2,2n \
            | bedtools merge -d {params.gap} \
            > {output} 2>{log}
        """


# ── 5. All-samples consensus peak set ─────────────────────────────────────────
rule consensus_peaks_all_samples:
    """
    Merge peaks from every sample across all conditions.
    Used as the global reference peak universe for DiffBind.
    """
    input:
        peaks = expand(
            _RESULTS + "/peaks/{sample}/{sample}_peaks.filtered.narrowPeak",
            sample=SAMPLES,
        ),
    output:
        _RESULTS + "/peaks/consensus/all_samples_consensus_peaks.bed",
    log:
        _RESULTS + "/logs/consensus_peaks/all_samples.log",
    conda:
        "../../envs/macs3.yaml"
    params:
        gap = _MERGE_GAP,
    shell:
        """
        cat {input.peaks} \
            | awk 'OFS="\\t" {{print $1,$2,$3}}' \
            | grep -vE '_(random|alt|Un)|chrM|chrEBV' \
            | sort -k1,1 -k2,2n \
            | bedtools merge -d {params.gap} \
            > {output} 2>{log}
        """


# ── 6. Peak count summary ─────────────────────────────────────────────────────
rule peak_count_summary:
    """Print a peak count table across all samples (for QC inspection)."""
    input:
        peaks = expand(
            _RESULTS + "/peaks/{sample}/{sample}_peaks.filtered.narrowPeak",
            sample=SAMPLES,
        ),
    output:
        _RESULTS + "/qc/peak_counts/peak_count_summary.txt",
    run:
        with open(output[0], "w") as fout:
            fout.write("sample\tpeak_count\n")
            for peak_file, sample in zip(input.peaks, SAMPLES):
                with open(peak_file) as fin:
                    count = sum(1 for _ in fin)
                fout.write(f"{sample}\t{count}\n")
