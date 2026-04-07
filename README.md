# Bulk ATAC-seq Analysis Pipeline

Snakemake workflow for bulk ATAC-seq data analysis:
**3 conditions × 10 replicates = 30 samples**.

---

## Pipeline Overview

```
data/bam/{sample}.bam
        │
        ▼
┌─────────────────────────────────────────────────────────┐
│  1. PREPROCESSING (workflow/rules/preprocessing.smk)     │
│                                                          │
│  sort_bam → mark_duplicates (Picard) → filter_bam       │
│  (MAPQ≥30, -F 1804, -f 2, chrM removal)                │
│  → remove_blacklist (ENCODE hg38 blacklist)              │
│  → results/bam/final/{sample}.filtered.bam(.bai)         │
└─────────────────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────┐
│  2. QC (workflow/rules/qc.smk)                          │
│                                                          │
│  • samtools flagstat / stats                             │
│  • Picard duplicate metrics                              │
│  • Library complexity: NRF / PBC1 / PBC2                 │
│  • Fragment-size distribution (deeptools)                │
│  • RPKM bigWig (deeptools bamCoverage)                   │
│  • TSS enrichment heatmap (deeptools computeMatrix)      │
│  • FRiP score (bedtools + samtools)                      │
│  • MultiQC aggregate report                              │
└─────────────────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────┐
│  3. PEAK CALLING (workflow/rules/peak_calling.smk)       │
│                                                          │
│  MACS3 (BAMPE, -q 0.05, --call-summits)                  │
│  → per-sample narrowPeak                                 │
│  → blacklist-filtered narrowPeak                         │
│  → per-condition consensus BED (bedtools merge)          │
│  → all-samples consensus BED                             │
└─────────────────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────┐
│  4. DIFFERENTIAL ACCESSIBILITY (workflow/rules/           │
│     diffbind.smk + scripts/diffbind_analysis.R)          │
│                                                          │
│  DiffBind → dba.count (SummarizeOverlaps)                │
│  → dba.normalize (DESeq2)                                │
│  → QC: correlation heatmap, PCA, peak overlap Venn       │
│  → dba.contrast (all pairwise)                           │
│  → dba.analyze (DESeq2)                                  │
│  → per-contrast: volcano, MA plot, CSV tables            │
│  → ChIPseeker annotation for significant peaks           │
└─────────────────────────────────────────────────────────┘
```

---

## Repository Structure

```
atac-seq-pipeline/
├── Snakefile                        Main workflow entry point
├── config/
│   ├── config.yaml                  Pipeline parameters (genome, thresholds…)
│   └── samples.csv                  Sample table (sample, condition, replicate)
├── workflow/
│   └── rules/
│       ├── preprocessing.smk        BAM filtering rules
│       ├── qc.smk                   QC metric rules
│       ├── peak_calling.smk         MACS3 + consensus peak rules
│       └── diffbind.smk             DiffBind rules
├── scripts/
│   ├── compute_library_complexity.py  NRF / PBC1 / PBC2
│   ├── frip_score.sh                  FRiP calculation
│   ├── make_diffbind_samplesheet.py   Build DiffBind CSV
│   └── diffbind_analysis.R            Full DiffBind R pipeline
├── envs/
│   ├── preprocessing.yaml           samtools, picard, bedtools, pysam, multiqc
│   ├── macs3.yaml                   macs3, bedtools
│   ├── deeptools.yaml               deeptools
│   └── diffbind.yaml                R + DiffBind + ChIPseeker
├── resources/
│   └── README.md                   Instructions to download reference files
└── data/
    └── bam/                        Place input BAM files here: {sample}.bam
```

---

## Quick Start

### 1. Clone and prepare data

```bash
git clone <this-repo> atac-seq-pipeline
cd atac-seq-pipeline

# Place your input BAMs
mkdir -p data/bam
cp /path/to/your/*.bam data/bam/
# Expected names: condA_rep01.bam … condC_rep10.bam
```

### 2. Download reference files

See [resources/README.md](resources/README.md) for download commands.

```bash
# Quick example (hg38)
cd resources
wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
gunzip hg38-blacklist.v2.bed.gz
cd ..
```

### 3. Edit configuration (if needed)

```bash
# Adjust genome, thresholds, thread count in:
nano config/config.yaml

# Rename samples to match your BAM filenames in:
nano config/samples.csv
```

### 4. Install Snakemake

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake>=7.32 mamba
conda activate snakemake
```

### 5. Dry run (check DAG)

```bash
snakemake --use-conda --cores 1 -n -p
# Visualise DAG (optional)
snakemake --use-conda --cores 1 --dag | dot -Tsvg > dag.svg
```

### 6. Run the full pipeline

```bash
# Local multi-core
snakemake --use-conda --cores 32 -p

# SLURM cluster (create cluster profile first)
snakemake --use-conda --profile slurm --jobs 50
```

### 7. Run only a specific step

```bash
# Only preprocessing + QC
snakemake --use-conda --cores 32 -p \
    results/qc/multiqc/multiqc_report.html

# Only DiffBind
snakemake --use-conda --cores 8 -p \
    results/diffbind/differential_peaks_all_contrasts.csv
```

---

## Expected Outputs

| Path | Description |
|------|-------------|
| `results/bam/final/*.filtered.bam` | Deduplicated, filtered BAMs |
| `results/qc/multiqc/multiqc_report.html` | Aggregated QC report |
| `results/qc/library_complexity/*.txt` | NRF / PBC1 / PBC2 |
| `results/qc/frip/*.frip.txt` | FRiP scores |
| `results/qc/tss_enrichment/*_TSS_heatmap.pdf` | TSS enrichment plots |
| `results/peaks/{sample}/*_peaks.narrowPeak` | Per-sample MACS3 peaks |
| `results/peaks/consensus/*.bed` | Condition / all-sample consensus peaks |
| `results/diffbind/plots/pca_plot.pdf` | PCA of all samples |
| `results/diffbind/plots/correlation_heatmap.pdf` | Sample correlation |
| `results/diffbind/plots/*_volcano.pdf` | Per-contrast volcano plots |
| `results/diffbind/results/*_significant_peaks.csv` | Differential peak tables |
| `results/diffbind/results/*_significant_annotated.csv` | ChIPseeker annotation |
| `results/diffbind/differential_peaks_all_contrasts.csv` | All contrasts combined |

---

## Quality Control Thresholds (ENCODE ATAC-seq)

| Metric | Concerning | Acceptable | Ideal |
|--------|-----------|------------|-------|
| NRF | < 0.7 | 0.7 – 0.9 | ≥ 0.9 |
| PBC1 | < 0.7 | 0.7 – 0.9 | ≥ 0.9 |
| PBC2 | < 1.0 | 1.0 – 3.0 | ≥ 3.0 |
| FRiP | < 0.1 | 0.1 – 0.2 | ≥ 0.3 |

---

## Configuration Reference (`config/config.yaml`)

| Key | Default | Description |
|-----|---------|-------------|
| `genome` | `hg38` | Reference genome name |
| `genome_size` | `hs` | MACS3 genome size flag (`hs`/`mm`/`dm`/`ce`) |
| `blacklist` | `resources/hg38-blacklist.v2.bed` | ENCODE blacklist BED |
| `paired_end` | `true` | `false` for single-end data |
| `min_mapq` | `30` | Minimum alignment quality |
| `macs3.qvalue` | `0.05` | MACS3 q-value cut-off |
| `macs3.call_summits` | `true` | Enable summit detection |
| `diffbind.fdr_threshold` | `0.05` | FDR cut-off for significant peaks |
| `diffbind.fold_threshold` | `1.0` | Minimum \|log2FC\| |
| `diffbind.summits` | `250` | Re-center peaks ±N bp around summit |
| `threads` | `8` | Threads per rule |

---

## Adapting to Other Organisms

Edit `config/config.yaml`:
```yaml
genome:      mm10
genome_size: mm
blacklist:   resources/mm10-blacklist.v2.bed
tss_bed:     resources/mm10_tss_2kb.bed
mito_chr:    chrM
```

And update `scripts/diffbind_analysis.R`:
```r
# Replace TxDb and orgDb:
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
txdb   <- TxDb.Mmusculus.UCSC.mm10.knownGene
annodb <- "org.Mm.eg.db"
```

---

## Dependencies

All tools are managed by conda environments in `envs/`. Installed automatically
when running Snakemake with `--use-conda`.

| Tool | Version | Purpose |
|------|---------|---------|
| samtools | 1.19 | BAM manipulation |
| picard | 3.1.1 | Duplicate removal |
| bedtools | 2.31.1 | Genomic interval operations |
| MACS3 | 3.0.1 | Peak calling |
| deeptools | 3.5.4 | QC metrics, bigWig, TSS enrichment |
| MultiQC | 1.19 | QC aggregation |
| DiffBind | 3.12.0 | Differential accessibility |
| ChIPseeker | 1.38.0 | Peak annotation |
| DESeq2 | (via DiffBind) | Statistical testing |

---

## Citation

If you use this pipeline, please cite:

- **MACS3**: Zhang et al. (2008) *Genome Biology* — peak calling
- **DiffBind**: Ross-Innes et al. (2012) *Nature*; Stark & Brown (2011)
- **DESeq2**: Love et al. (2014) *Genome Biology*
- **ChIPseeker**: Yu et al. (2015) *Bioinformatics*
- **deeptools**: Ramírez et al. (2016) *Nucleic Acids Research*
- **Picard**: Broad Institute — http://broadinstitute.github.io/picard/
- **ENCODE blacklist**: Amemiya et al. (2019) *Scientific Reports*
