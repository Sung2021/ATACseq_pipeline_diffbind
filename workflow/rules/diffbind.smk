"""
DiffBind Rules
==============
Input : filtered BAMs + filtered narrowPeak files (from preprocessing + peak_calling)
Output:
  results/diffbind/samplesheet.csv              – DiffBind-format samplesheet
  results/diffbind/diffbind_analysis.RData      – full DiffBind object + results
  results/diffbind/plots/*.pdf                  – QC and result plots
  results/diffbind/results/*_peaks.csv          – per-contrast peak tables
  results/diffbind/differential_peaks_all_contrasts.csv
"""

_RESULTS  = config["results_dir"]
_THREADS  = config.get("threads", 8)

_DB_FDR   = config["diffbind"].get("fdr_threshold",  0.05)
_DB_FOLD  = config["diffbind"].get("fold_threshold", 1.0)
_DB_METH  = config["diffbind"].get("method",         "DESeq2")
_DB_SUM   = config["diffbind"].get("summits",        250)
_DB_MINOV = config["diffbind"].get("min_overlap",    2)


# ── 1. Generate DiffBind samplesheet ──────────────────────────────────────────
rule make_diffbind_samplesheet:
    """
    Create the CSV samplesheet consumed by DiffBind::dba().
    Columns: SampleID, Condition, Replicate, bamReads, Peaks, PeakCaller
    """
    input:
        samples    = config["samples"],
        bam_files  = expand(
            _RESULTS + "/bam/final/{sample}.filtered.bam",
            sample=SAMPLES,
        ),
        peak_files = expand(
            _RESULTS + "/peaks/{sample}/{sample}_peaks.filtered.narrowPeak",
            sample=SAMPLES,
        ),
    output:
        csv = _RESULTS + "/diffbind/samplesheet.csv",
    log:
        _RESULTS + "/logs/diffbind/make_samplesheet.log",
    params:
        bam_dir   = _RESULTS + "/bam/final",
        peaks_dir = _RESULTS + "/peaks",
    shell:
        """
        python scripts/make_diffbind_samplesheet.py \
            --samples   {input.samples} \
            --bam-dir   {params.bam_dir} \
            --peaks-dir {params.peaks_dir} \
            --output    {output.csv} \
            2>{log}
        """


# ── 2. Run full DiffBind analysis ─────────────────────────────────────────────
rule run_diffbind:
    """
    Execute the DiffBind R pipeline:
      dba()           – load samplesheet
      dba.count()     – count reads in peaks (SummarizeOverlaps)
      dba.normalize() – DESeq2 / TMM normalisation
      dba.contrast()  – all pairwise contrasts between conditions
      dba.analyze()   – DESeq2 differential analysis
      dba.report()    – export results tables

    Plots generated (all PDF):
      correlation_heatmap.pdf  – sample correlation heatmap (Pearson)
      pca_plot.pdf             – PCA coloured by condition
      {contrast}_volcano.pdf   – per-contrast volcano plot
      {contrast}_MA.pdf        – per-contrast MA plot
      {contrast}_annotation_bar.pdf – peak annotation bar (ChIPseeker)
      {contrast}_tss_distance.pdf   – distance to TSS histogram

    Output tables:
      results/{contrast}_all_peaks.csv         – all tested peaks
      results/{contrast}_significant_peaks.csv – FDR < threshold, |FC| ≥ threshold
      results/{contrast}_significant_annotated.csv – ChIPseeker annotation
    """
    input:
        samplesheet  = _RESULTS + "/diffbind/samplesheet.csv",
        # Ensure all BAMs and peaks exist before starting R
        bam_indices  = expand(
            _RESULTS + "/bam/final/{sample}.filtered.bam.bai",
            sample=SAMPLES,
        ),
        consensus    = _RESULTS + "/peaks/consensus/all_samples_consensus_peaks.bed",
    output:
        rdata      = _RESULTS + "/diffbind/diffbind_analysis.RData",
        pca        = _RESULTS + "/diffbind/plots/pca_plot.pdf",
        heatmap    = _RESULTS + "/diffbind/plots/correlation_heatmap.pdf",
        all_csv    = _RESULTS + "/diffbind/differential_peaks_all_contrasts.csv",
    log:
        _RESULTS + "/logs/diffbind/run_diffbind.log",
    threads: _THREADS
    resources:
        mem_mb = config.get("memory_mb", 32000),  # DiffBind / DESeq2 is memory-heavy
    conda:
        "../../envs/diffbind.yaml"
    params:
        outdir     = _RESULTS + "/diffbind",
        fdr        = _DB_FDR,
        fold       = _DB_FOLD,
        method     = _DB_METH,
        summits    = _DB_SUM,
        min_overlap= _DB_MINOV,
        threads    = _THREADS,
    shell:
        """
        Rscript scripts/diffbind_analysis.R \
            --samplesheet {input.samplesheet} \
            --outdir      {params.outdir} \
            --fdr         {params.fdr} \
            --fold        {params.fold} \
            --method      {params.method} \
            --summits     {params.summits} \
            --min-overlap {params.min_overlap} \
            --threads     {params.threads} \
            2>{log}
        """
