"""
QC Rules
=========
Metrics collected per sample
-----------------------------
  flagstat            – samtools flagstat (basic read counts)
  picard_metrics      – duplicate rate (produced by preprocessing.smk)
  library_complexity  – NRF / PBC1 / PBC2  (scripts/compute_library_complexity.py)
  fragment_size       – paired-end fragment-size distribution (deeptools)
  bigwig              – RPKM-normalised bigWig (deeptools bamCoverage)
  tss_enrichment      – TSS enrichment heatmap (deeptools computeMatrix + plotHeatmap)
  frip_score          – Fraction of Reads in Peaks (scripts/frip_score.sh)
  multiqc             – aggregate report

All outputs under {results}/qc/
"""

_RESULTS  = config["results_dir"]
_THREADS  = config.get("threads", 8)
_TSS_BED  = config.get("tss_bed",     "resources/hg38_tss_2kb.bed")
_BLACKLIST = config.get("blacklist",   "resources/hg38-blacklist.v2.bed")


# ── 1. samtools flagstat ───────────────────────────────────────────────────────
rule flagstat:
    """Basic alignment statistics for MultiQC."""
    input:
        bam = _RESULTS + "/bam/final/{sample}.filtered.bam",
    output:
        _RESULTS + "/qc/flagstat/{sample}.flagstat",
    log:
        _RESULTS + "/logs/flagstat/{sample}.log",
    threads: 2
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        "samtools flagstat -@ {threads} {input.bam} > {output} 2>{log}"


# ── 2. samtools stats (for MultiQC insert-size section) ───────────────────────
rule samtools_stats:
    """Detailed samtools stats including insert-size histogram."""
    input:
        bam = _RESULTS + "/bam/final/{sample}.filtered.bam",
    output:
        _RESULTS + "/qc/samtools_stats/{sample}.stats",
    log:
        _RESULTS + "/logs/samtools_stats/{sample}.log",
    threads: 2
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        "samtools stats -@ {threads} {input.bam} > {output} 2>{log}"


# ── 3. Library complexity (NRF / PBC1 / PBC2) ────────────────────────────────
rule library_complexity:
    """
    Compute ENCODE ATAC-seq library complexity metrics using
    scripts/compute_library_complexity.py.

    Thresholds (ENCODE):
      NRF  ≥ 0.9   (distinct / total)
      PBC1 ≥ 0.7   (M1 / distinct)
      PBC2 ≥ 1.0   (M1 / M2)
    """
    input:
        bam = _RESULTS + "/bam/final/{sample}.filtered.bam",
    output:
        _RESULTS + "/qc/library_complexity/{sample}.complexity.txt",
    log:
        _RESULTS + "/logs/library_complexity/{sample}.log",
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        """
        python scripts/compute_library_complexity.py \
            {input.bam} \
            {wildcards.sample} \
            {output} \
            2>{log}
        """


# ── 4. Fragment-size distribution ─────────────────────────────────────────────
rule fragment_size:
    """
    Paired-end fragment-size distribution using deeptools bamPEFragmentSize.
    The characteristic nucleosomal ladder (sub-200 / 200 / 400 bp) confirms
    good ATAC-seq library quality.
    """
    input:
        bam = _RESULTS + "/bam/final/{sample}.filtered.bam",
        bai = _RESULTS + "/bam/final/{sample}.filtered.bam.bai",
    output:
        txt  = _RESULTS + "/qc/fragment_size/{sample}_fragment_size.txt",
        hist = _RESULTS + "/qc/fragment_size/{sample}_fragment_size.png",
    log:
        _RESULTS + "/logs/fragment_size/{sample}.log",
    threads: _THREADS
    conda:
        "../../envs/deeptools.yaml"
    shell:
        """
        bamPEFragmentSize \
            --bamfiles {input.bam} \
            --samplesLabel {wildcards.sample} \
            --outRawFragmentLengths {output.txt} \
            --histogram {output.hist} \
            --numberOfProcessors {threads} \
            2>{log}
        """


# ── 5. BigWig (RPKM-normalised) ───────────────────────────────────────────────
rule bam_to_bigwig:
    """
    Generate RPKM-normalised bigWig for genome browser and TSS enrichment.
    Uses deeptools bamCoverage with --centerReads for ATAC-seq cut sites.
    """
    input:
        bam = _RESULTS + "/bam/final/{sample}.filtered.bam",
        bai = _RESULTS + "/bam/final/{sample}.filtered.bam.bai",
    output:
        bw = _RESULTS + "/bigwig/{sample}.rpkm.bw",
    log:
        _RESULTS + "/logs/bam_to_bigwig/{sample}.log",
    threads: _THREADS
    conda:
        "../../envs/deeptools.yaml"
    params:
        blacklist = _BLACKLIST,
        binsize   = 10,
    shell:
        """
        bamCoverage \
            --bam {input.bam} \
            --outFileName {output.bw} \
            --outFileFormat bigwig \
            --normalizeUsing RPKM \
            --binSize {params.binsize} \
            --centerReads \
            --blackListFileName {params.blacklist} \
            --numberOfProcessors {threads} \
            2>{log}
        """


# ── 6. TSS enrichment ─────────────────────────────────────────────────────────
rule tss_enrichment:
    """
    Compute TSS enrichment score using deeptools computeMatrix + plotHeatmap.
    Reference TSS BED file: config['tss_bed']
    Download from: see resources/README.md
    """
    input:
        bw  = _RESULTS + "/bigwig/{sample}.rpkm.bw",
        tss = _TSS_BED,
    output:
        matrix  = temp(_RESULTS + "/qc/tss_enrichment/{sample}_matrix.gz"),
        heatmap = _RESULTS + "/qc/tss_enrichment/{sample}_TSS_heatmap.pdf",
        profile = _RESULTS + "/qc/tss_enrichment/{sample}_TSS_profile.pdf",
        scores  = _RESULTS + "/qc/tss_enrichment/{sample}_TSS_scores.tab",
    log:
        _RESULTS + "/logs/tss_enrichment/{sample}.log",
    threads: _THREADS
    conda:
        "../../envs/deeptools.yaml"
    shell:
        """
        computeMatrix reference-point \
            --referencePoint TSS \
            --scoreFileName {input.bw} \
            --regionsFileName {input.tss} \
            --outFileName {output.matrix} \
            --outFileNameMatrix {output.scores} \
            --beforeRegionStartLength 2000 \
            --afterRegionStartLength  2000 \
            --binSize 10 \
            --skipZeros \
            --numberOfProcessors {threads} \
            2>{log}

        plotHeatmap \
            --matrixFile {output.matrix} \
            --outFileName {output.heatmap} \
            --samplesLabel "{wildcards.sample}" \
            --plotTitle "TSS Enrichment – {wildcards.sample}" \
            --colorMap RdBu_r \
            2>>{log}

        plotProfile \
            --matrixFile {output.matrix} \
            --outFileName {output.profile} \
            --samplesLabel "{wildcards.sample}" \
            --plotTitle "TSS Profile – {wildcards.sample}" \
            2>>{log}
        """


# ── 7. FRiP score ─────────────────────────────────────────────────────────────
rule frip_score:
    """
    Fraction of Reads in Peaks (FRiP) – key ATAC-seq quality metric.
    ENCODE threshold: FRiP ≥ 0.2 for ATAC-seq.
    Uses bedtools intersect + samtools view -c.
    """
    input:
        bam   = _RESULTS + "/bam/final/{sample}.filtered.bam",
        bai   = _RESULTS + "/bam/final/{sample}.filtered.bam.bai",
        peaks = _RESULTS + "/peaks/{sample}/{sample}_peaks.narrowPeak",
    output:
        _RESULTS + "/qc/frip/{sample}.frip.txt",
    log:
        _RESULTS + "/logs/frip/{sample}.log",
    conda:
        "../../envs/preprocessing.yaml"
    shell:
        """
        bash scripts/frip_score.sh \
            {input.bam} \
            {input.peaks} \
            {wildcards.sample} \
            {output} \
            2>{log}
        """


# ── 8. MultiQC aggregate report ───────────────────────────────────────────────
rule multiqc:
    """
    Aggregate all QC outputs into a single MultiQC HTML report.
    Collects: flagstat, samtools_stats, picard duplicate metrics,
              fragment sizes, TSS enrichment scores, FRiP.
    """
    input:
        # Wait for all per-sample QC outputs first
        flagstat   = expand(
            _RESULTS + "/qc/flagstat/{sample}.flagstat",
            sample=SAMPLES,
        ),
        stats      = expand(
            _RESULTS + "/qc/samtools_stats/{sample}.stats",
            sample=SAMPLES,
        ),
        complexity = expand(
            _RESULTS + "/qc/library_complexity/{sample}.complexity.txt",
            sample=SAMPLES,
        ),
        frip       = expand(
            _RESULTS + "/qc/frip/{sample}.frip.txt",
            sample=SAMPLES,
        ),
        picard     = expand(
            _RESULTS + "/qc/picard_metrics/{sample}.dup_metrics.txt",
            sample=SAMPLES,
        ),
    output:
        html = _RESULTS + "/qc/multiqc/multiqc_report.html",
    log:
        _RESULTS + "/logs/multiqc/multiqc.log",
    conda:
        "../../envs/preprocessing.yaml"
    params:
        outdir = _RESULTS + "/qc/multiqc",
        title  = config.get("multiqc_title", "Bulk ATAC-seq QC Report"),
        indir  = _RESULTS + "/qc",
    shell:
        """
        multiqc \
            {params.indir} \
            --outdir {params.outdir} \
            --title "{params.title}" \
            --force \
            2>{log}
        """
