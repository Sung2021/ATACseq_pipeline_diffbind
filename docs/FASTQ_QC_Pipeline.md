# FASTQ QC and Preprocessing Pipeline

Standalone **Snakemake workflow** for FASTQ quality control, trimming, and preprocessing.

---

## 🎯 Pipeline Overview

```
data/raw_fastq/
       │
       ├─→ [1] FastQC ────────────────┐
       │                              │
       │                              ├─→ [2] MultiQC (Raw report)
       │                              │
       └─→ [3] Parse QC metrics ──────┤
                                       │
              [4] Determine fastp params
                                       │
       ┌──────────────────────────────┘
       │
       ├─→ [5] fastp (trim + filter)
       │
       ├─→ [6] FastQC (on trimmed reads)
       │
       └─→ [7] MultiQC (Trimmed report)

Output:
    results/qc_preprocessing/
    ├── qc/fastqc/               → HTML reports
    ├── qc/multiqc/              → Aggregated reports
    ├── trimmed/                 → fastp-processed FASTQ
    └── qc_metrics/              → Parsed metrics & recommendations
```

### Pipeline Steps

| Step | Tool | Output |
|------|------|--------|
| 1 | **FastQC** | Per-sample quality assessment (HTML + ZIP) |
| 2 | **MultiQC** | Aggregated QC report (HTML) |
| 3 | **parse_fastqc.py** | CSV with QC metrics |
| 4 | **determine_fastp_params.py** | Recommended fastp parameters (TXT + JSON) |
| 5 | **fastp** | Trimmed FASTQ + quality report |
| 6 | **FastQC** | QC on trimmed reads |
| 7 | **MultiQC** | Aggregated report (before/after comparison) |

---

## 📁 Input Data

Place your FASTQ files in `data/raw_fastq/` with standard naming:

### Paired-end sequencing:
```
data/raw_fastq/
  ├── SampleA_R1.fastq.gz
  ├── SampleA_R2.fastq.gz
  ├── SampleB_R1.fastq.gz
  ├── SampleB_R2.fastq.gz
  └── ...
```

### Single-end sequencing:
```
data/raw_fastq/
  ├── SampleA_R1.fastq.gz
  ├── SampleB_R1.fastq.gz
  └── ...
```

Then update `config/preprocessing_qc.yaml`:
```yaml
paired_end: false  # For single-end
```

---

## ⚙️ Configuration

Edit `config/preprocessing_qc.yaml`:

```yaml
raw_fastq_dir: data/raw_fastq        # Input directory
results_dir: results/qc_preprocessing # Output directory

paired_end: true                      # true for PE, false for SE

min_quality: 20                       # Quality threshold for fastp
min_length: 20                        # Minimum read length (bp)

run_qc_after_trimming: true          # Re-run FastQC on trimmed reads
```

---

## 🚀 Quick Start

### 1. Prepare input data
```bash
# Create input directory
mkdir -p data/raw_fastq

# Place your FASTQ files there
cp /path/to/your/fastq/*.fastq.gz data/raw_fastq/
```

### 2. Run pipeline (dry-run first)
```bash
# Test: see what will run without executing
snakemake -s Snakefile_qualitycontrol \
    --use-conda --cores 16 --dry-run -p
```

### 3. Run pipeline (full execution)
```bash
# Local execution with 16 cores
snakemake -s Snakefile_qualitycontrol \
    --use-conda --cores 16 -p

# With HPC/SLURM
snakemake -s Snakefile_qualitycontrol \
    --use-conda --cores 32 --profile slurm

# Faster: parallel FastQC per-sample
snakemake -s Snakefile_qualitycontrol \
    --use-conda --cores 32 \
    --jobs 16 -p
```

### 4. View reports
```bash
# Open in browser
open results/qc_preprocessing/qc/multiqc/multiqc_report_raw.html
open results/qc_preprocessing/qc/multiqc/multiqc_report_trimmed.html

# View fastp recommendations
cat results/qc_preprocessing/qc_metrics/fastp_params_recommendation.txt
```

---

## 📊 Output Files

```
results/qc_preprocessing/
│
├── qc/
│   ├── fastqc/
│   │   ├── SampleA_R1_fastqc.html
│   │   ├── SampleA_R1_fastqc.zip
│   │   ├── SampleA_R2_fastqc.html
│   │   ├── SampleA_R2_fastqc.zip
│   │   ├── SampleA_R1.trimmed_fastqc.html    ← After trimming
│   │   ├── SampleA_R2.trimmed_fastqc.html    ← After trimming
│   │   └── ...
│   │
│   └── multiqc/
│       ├── multiqc_report_raw.html           ← Before trimming
│       ├── multiqc_report_trimmed.html       ← After trimming
│       └── multiqc_data/
│           └── multiqc_general_stats.txt
│
├── trimmed/
│   ├── SampleA_R1.trimmed.fastq.gz          ← Main output!
│   ├── SampleA_R2.trimmed.fastq.gz          ← Main output!
│   ├── SampleA.fastp.json                   ← Trimming report
│   ├── SampleA.fastp.html                   ← Trimming report
│   ├── SampleB_R1.trimmed.fastq.gz
│   ├── SampleB_R2.trimmed.fastq.gz
│   └── ...
│
└── qc_metrics/
    ├── qc_metrics_summary.csv                ← QC metrics table
    ├── fastp_params_recommendation.txt       ← Recommendations
    └── fastp_params.json                     ← Parameters (JSON)
```

---

## 📈 Quality Metrics

### FastQC produces:
- Per-base quality score
- Per-sequence quality scores
- GC content distribution
- Adapter content
- Sequence length distribution
- Duplicate sequences

### fastp reports:
- Trimming statistics
- Filtering statistics
- Adapter detection results

### MultiQC aggregates:
- All FastQC results in one HTML
- Comparison across all samples
- Interactive plots

---

## 🔧 Advanced Usage

### Run only specific rules
```bash
# Only FastQC on raw reads
snakemake -s Snakefile_qualitycontrol \
    --use-conda --cores 16 \
    -R fastqc_raw

# Only fastp trimming (assumes QC already done)
snakemake -s Snakefile_qualitycontrol \
    --use-conda --cores 16 \
    -R fastp_trim
```

### Skip trimming (only QC)
Edit `Snakefile_qualitycontrol` rule `all` to remove fastp and trimmed FastQC rules:
```python
rule all:
    input:
        expand(f"{FASTQC_DIR}/{{sample}}_R1_fastqc.html", sample=SAMPLES),
        expand(f"{FASTQC_DIR}/{{sample}}_R2_fastqc.html", sample=SAMPLES),
        f"{MULTIQC_DIR}/multiqc_report_raw.html",
```

### Custom fastp parameters
Edit `Snakefile_qualitycontrol` `rule fastp_trim` shell command:
```bash
--window_size 5 \
--window_mean_quality 15 \
--poly_x_min_len 10 \
```

### Generate conda environment manually
```bash
# If not automatically installing
conda env create -f envs/qc.yaml -n atac-qc
conda activate atac-qc

# Then run without --use-conda
snakemake -s Snakefile_qualitycontrol --cores 16 -p
```

---

## ⚠️ Troubleshooting

### Error: "No FASTQ files found"
- Check `config/preprocessing_qc.yaml` → `raw_fastq_dir` path
- Verify file naming: `{sample}_R1.fastq.gz` / `{sample}_R2.fastq.gz`
- Ensure `.fastq.gz` extension (must be gzipped)

### FastQC module failures
- Some samples might have quality issues
- Check individual FastQC HTML reports
- This is informational; pipeline continues

### fastp: "Connection reset"
- Usually memory issue
- Increase `memory_mb` in config
- Reduce `--threads` in fastp rule

### conda: "Solving environment is taking too long"
- Use `mamba` instead: `--conda-frontend mamba`
- Or pre-install: `conda env create -f envs/qc.yaml`

---

## 📚 Next Steps

After preprocessing, use trimmed FASTQ files for:
1. **Alignment**: bwa, bowtie2 (for ATAC-seq BAM)
2. **Assembly**: spades, megahit (if needed)
3. **Mapping**: minimap2 (long-read)

For ATAC-seq BAM generation, see main pipeline:
```bash
snakemake --use-conda --cores 32 -p
```

---

## 📖 References

- **FastQC**: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- **MultiQC**: https://multiqc.info/
- **fastp**: https://github.com/OpenGene/fastp
- **Snakemake**: https://snakemake.readthedocs.io/

---

## 📝 Notes

- **Thread allocation**: FastQC (2 threads), fastp (4 threads) — adjust in config
- **Memory**: ~4-8 GB per sample (adjust in config)
- **Typical runtime**: 5-10 min per sample (depends on read count)
- **ATAC-seq specific**: Generally uses **more lenient** trimming (Q20, len>20bp) than RNA-seq

