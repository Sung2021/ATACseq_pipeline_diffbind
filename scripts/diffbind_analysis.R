#!/usr/bin/env Rscript
# =============================================================================
# diffbind_analysis.R
# Full DiffBind differential ATAC-seq analysis pipeline
#
# Steps:
#   1. Load samplesheet → dba()
#   2. Count reads     → dba.count()
#   3. Normalise       → dba.normalize()
#   4. QC plots        → correlation heatmap, PCA, peak overlap Venn
#   5. Set contrasts   → dba.contrast()  (all pairwise)
#   6. Analyse         → dba.analyze()   (DESeq2)
#   7. Report & plots  → per-contrast volcano, MA, annotation
#   8. Peak annotation → ChIPseeker
#   9. Save RData
#
# Usage (called by Snakemake rule run_diffbind):
#   Rscript diffbind_analysis.R \
#       --samplesheet results/diffbind/samplesheet.csv \
#       --outdir      results/diffbind \
#       --fdr         0.05 \
#       --fold        1.0 \
#       --method      DESeq2 \
#       --summits     250 \
#       --min-overlap 2 \
#       --threads     8
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(DiffBind)
  library(BiocParallel)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(GenomicRanges)
  library(IRanges)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
})

# ─── CLI argument parsing ─────────────────────────────────────────────────────
option_list <- list(
  make_option("--samplesheet", type = "character", help = "DiffBind CSV samplesheet"),
  make_option("--outdir",      type = "character", help = "Output directory"),
  make_option("--fdr",         type = "double",    default = 0.05),
  make_option("--fold",        type = "double",    default = 1.0,
              help = "Minimum |log2FC|"),
  make_option("--method",      type = "character", default = "DESeq2"),
  make_option("--summits",     type = "integer",   default = 250,
              help = "Re-centre peaks ±N bp around summit (0 = disable)"),
  make_option("--min-overlap", type = "integer",   default = 2,
              dest = "min_overlap"),
  make_option("--threads",     type = "integer",   default = 1)
)
opt <- parse_args(OptionParser(option_list = option_list))

# ─── Setup ────────────────────────────────────────────────────────────────────
dir.create(file.path(opt$outdir, "plots"),   showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$outdir, "results"), showWarnings = FALSE, recursive = TRUE)

# Parallel workers (conservative: BiocParallel SnowParam for multi-core)
if (opt$threads > 1) {
  register(MulticoreParam(workers = opt$threads))
} else {
  register(SerialParam())
}

cat("=== DiffBind ATAC-seq Analysis ===\n")
cat(sprintf("Samplesheet : %s\n", opt$samplesheet))
cat(sprintf("Output dir  : %s\n", opt$outdir))
cat(sprintf("FDR cut-off : %.3g\n", opt$fdr))
cat(sprintf("Fold cut-off: %.2f\n", opt$fold))
cat(sprintf("Method      : %s\n",   opt$method))
cat(sprintf("Threads     : %d\n",   opt$threads))

txdb  <- TxDb.Hsapiens.UCSC.hg38.knownGene
annodb <- "org.Hs.eg.db"

# ─── 1. Load samplesheet ──────────────────────────────────────────────────────
cat("\n[1/9] Loading samplesheet ...\n")
dba_obj <- dba(sampleSheet = opt$samplesheet)
print(dba_obj)

# ─── 2. Count reads in peaks ─────────────────────────────────────────────────
cat("\n[2/9] Counting reads in peaks (SummarizeOverlaps) ...\n")
count_args <- list(
  DBA            = dba_obj,
  bUseSummarizeOverlaps = TRUE
)
if (opt$summits > 0) {
  count_args$summits <- opt$summits
  cat(sprintf("       Re-centering peaks around summits ± %d bp\n", opt$summits))
}
dba_obj <- do.call(dba.count, count_args)

# ─── 3. Normalise ────────────────────────────────────────────────────────────
cat("\n[3/9] Normalising ...\n")
if (opt$method == "DESeq2") {
  dba_obj <- dba.normalize(dba_obj, method = DBA_DESEQ2)
} else {
  dba_obj <- dba.normalize(dba_obj, method = DBA_EDGER)
}

# ─── 4. QC plots ─────────────────────────────────────────────────────────────
cat("\n[4/9] Generating QC plots ...\n")

## 4a. Correlation heatmap
pdf(file.path(opt$outdir, "plots", "correlation_heatmap.pdf"), width = 11, height = 10)
dba.plotHeatmap(dba_obj, main = "Sample Pearson Correlation (read counts in peaks)")
dev.off()

## 4b. PCA – all samples
pdf(file.path(opt$outdir, "plots", "pca_plot.pdf"), width = 8, height = 7)
dba.plotPCA(dba_obj, label = DBA_CONDITION,
            main = "PCA — All Samples (normalised counts)")
dev.off()

## 4c. Peak overlap Venn / UpSet (per condition)
conditions <- unique(dba_obj$samples$Condition)
if (length(conditions) <= 4) {
  pdf(file.path(opt$outdir, "plots", "peak_overlap_venn.pdf"), width = 8, height = 7)
  tryCatch(
    dba.plotVenn(dba_obj, dba_obj$masks$All,
                 main = "Peak Overlap Across All Samples"),
    error = function(e) message("Venn skipped: ", conditionMessage(e))
  )
  dev.off()
}

## 4d. Box plot of normalised counts
pdf(file.path(opt$outdir, "plots", "count_boxplot.pdf"), width = 10, height = 6)
dba.plotBox(dba_obj, main = "Normalised Read Counts in Peaks")
dev.off()

# ─── 5. Set contrasts ────────────────────────────────────────────────────────
cat("\n[5/9] Setting up pairwise contrasts ...\n")
if (opt$method == "DESeq2") {
  dba_obj <- dba.contrast(dba_obj, minMembers = opt$min_overlap,
                           method = DBA_DESEQ2)
} else {
  dba_obj <- dba.contrast(dba_obj, minMembers = opt$min_overlap,
                           method = DBA_EDGER)
}
contrast_tbl <- dba.show(dba_obj, bContrasts = TRUE)
cat("Contrasts defined:\n")
print(contrast_tbl)

# ─── 6. Differential analysis ────────────────────────────────────────────────
cat("\n[6/9] Running differential analysis ...\n")
if (opt$method == "DESeq2") {
  dba_obj <- dba.analyze(dba_obj, method = DBA_DESEQ2, bParallel = FALSE)
} else {
  dba_obj <- dba.analyze(dba_obj, method = DBA_EDGER,  bParallel = FALSE)
}

# ─── 7. Report and per-contrast plots ────────────────────────────────────────
cat("\n[7/9] Exporting results and generating contrast plots ...\n")
contrast_tbl <- dba.show(dba_obj, bContrasts = TRUE)
all_results  <- list()

method_flag <- if (opt$method == "DESeq2") DBA_DESEQ2 else DBA_EDGER

for (i in seq_len(nrow(contrast_tbl))) {
  grp1 <- as.character(contrast_tbl[i, "Group"])
  grp2 <- as.character(contrast_tbl[i, "Group2"])
  cname <- paste0(grp1, "_vs_", grp2)
  cat(sprintf("  Contrast %d: %s\n", i, cname))

  ## Full results (all peaks – for volcano)
  res_all <- tryCatch(
    dba.report(dba_obj, contrast = i, method = method_flag,
               th = 1, bNormalized = TRUE, bCalled = FALSE),
    error = function(e) { message("  ERROR dba.report: ", e$message); NULL }
  )
  if (is.null(res_all)) next

  ## Significant results
  res_sig <- dba.report(dba_obj, contrast = i, method = method_flag,
                        th = opt$fdr, fold = opt$fold,
                        bNormalized = TRUE)

  ## Convert to data frames
  res_df     <- as.data.frame(res_all)
  res_sig_df <- as.data.frame(res_sig)

  res_df$contrast    <- cname
  res_sig_df$contrast <- cname
  res_df$significant  <- res_df$FDR   <= opt$fdr & abs(res_df$Fold) >= opt$fold

  all_results[[cname]] <- res_df

  ## Write CSVs
  write.csv(res_df,
            file = file.path(opt$outdir, "results",
                             paste0(cname, "_all_peaks.csv")),
            row.names = FALSE)
  write.csv(res_sig_df,
            file = file.path(opt$outdir, "results",
                             paste0(cname, "_significant_peaks.csv")),
            row.names = FALSE)

  cat(sprintf("     Significant peaks (FDR<%.2f, |FC|≥%.1f): %d\n",
              opt$fdr, opt$fold, nrow(res_sig_df)))

  ## ── Volcano plot ──────────────────────────────────────────────────────────
  thresh_y <- -log10(opt$fdr)

  p_volcano <- ggplot(res_df, aes(x = Fold, y = -log10(FDR + 1e-300))) +
    geom_point(aes(colour = significant), alpha = 0.45, size = 0.9) +
    scale_colour_manual(
      values = c("TRUE" = "#c0392b", "FALSE" = "grey60"),
      labels = c("Not significant",
                 sprintf("FDR < %.2g  |log2FC| ≥ %.1f", opt$fdr, opt$fold)),
      name   = NULL
    ) +
    geom_vline(xintercept = c(-opt$fold, opt$fold),
               linetype = "dashed", colour = "navy", linewidth = 0.4) +
    geom_hline(yintercept = thresh_y,
               linetype = "dashed", colour = "navy", linewidth = 0.4) +
    labs(
      title = sprintf("Volcano — %s", cname),
      x = expression(log[2]~Fold~Change),
      y = expression(-log[10](FDR))
    ) +
    theme_bw(base_size = 12)

  ggsave(file.path(opt$outdir, "plots", paste0(cname, "_volcano.pdf")),
         p_volcano, width = 7, height = 6)

  ## ── MA plot (built-in DiffBind) ──────────────────────────────────────────
  pdf(file.path(opt$outdir, "plots", paste0(cname, "_MA.pdf")), width = 7, height = 5)
  dba.plotMA(dba_obj, contrast = i, method = method_flag, th = opt$fdr)
  title(main = sprintf("MA — %s", cname))
  dev.off()
}

## Combine all contrasts into one table
all_results_df <- do.call(rbind, all_results)
write.csv(all_results_df,
          file     = file.path(opt$outdir, "differential_peaks_all_contrasts.csv"),
          row.names = FALSE)

# ─── 8. Peak annotation (ChIPseeker) ─────────────────────────────────────────
cat("\n[8/9] Annotating significant peaks with ChIPseeker ...\n")

for (cname in names(all_results)) {
  res_df  <- all_results[[cname]]
  sig_df  <- res_df[res_df$significant, , drop = FALSE]

  if (nrow(sig_df) < 5) {
    cat(sprintf("  %s: only %d significant peaks — skipping annotation\n",
                cname, nrow(sig_df)))
    next
  }

  sig_gr <- GRanges(
    seqnames = sig_df$seqnames,
    ranges   = IRanges(sig_df$start, sig_df$end),
    strand   = "*"
  )

  anno <- tryCatch(
    annotatePeak(sig_gr, tssRegion = c(-3000, 3000),
                 TxDb = txdb, annoDb = annodb, verbose = FALSE),
    error = function(e) { message("  Annotation error: ", e$message); NULL }
  )
  if (is.null(anno)) next

  anno_df <- as.data.frame(anno)
  write.csv(anno_df,
            file = file.path(opt$outdir, "results",
                             paste0(cname, "_significant_annotated.csv")),
            row.names = FALSE)

  ## Annotation bar chart
  pdf(file.path(opt$outdir, "plots", paste0(cname, "_annotation_bar.pdf")),
      width = 8, height = 5)
  plotAnnoBar(anno, title = sprintf("Peak Annotation — %s", cname))
  dev.off()

  ## Distance to TSS histogram
  pdf(file.path(opt$outdir, "plots", paste0(cname, "_tss_distance.pdf")),
      width = 8, height = 5)
  plotDistToTSS(anno, title = sprintf("Distance to TSS — %s", cname))
  dev.off()
}

# ─── 9. Save RData ────────────────────────────────────────────────────────────
cat("\n[9/9] Saving RData ...\n")
save(dba_obj, all_results, file = file.path(opt$outdir, "diffbind_analysis.RData"))

# Print session info for reproducibility
cat("\n=== Session Info ===\n")
sessionInfo()

cat("\nDiffBind analysis complete.\n")
