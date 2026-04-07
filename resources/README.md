# Reference Resources

Place the following reference files in this directory before running the pipeline.
All files should use **hg38** (GRCh38) coordinates.

---

## 1. ENCODE Blacklist (`hg38-blacklist.v2.bed`)

Regions with anomalously high signal in sequencing data that should be excluded.

```bash
wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
gunzip hg38-blacklist.v2.bed.gz
```

---

## 2. TSS BED file (`hg38_tss_2kb.bed`)

Used for TSS enrichment score calculation (deeptools computeMatrix).
Each row is a ±1 kb window around a unique TSS.

Generate from UCSC refFlat or GENCODE:

```bash
# Option A – GENCODE v44 (recommended)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
gunzip gencode.v44.annotation.gtf.gz

awk '$3=="transcript"' gencode.v44.annotation.gtf \
  | awk 'OFS="\t" {
      match($0, /transcript_id "([^"]+)"/, tid)
      if ($7=="+") { tss=$4-1 } else { tss=$5-1 }
      start = (tss-1000 < 0) ? 0 : tss-1000
      print $1, start, tss+1000, tid[1], ".", $7
    }' \
  | grep -v "^#" \
  | sort -k1,1 -k2,2n \
  | uniq \
  > hg38_tss_2kb.bed

# Option B – UCSC knownGene (quick download)
mysql --user=genome --host=genome-mysql.ucsc.edu -A -e \
  "SELECT chrom, txStart, txEnd, name, 0, strand FROM knownGene;" hg38 \
  | awk 'NR>1 OFS="\t" {
      tss = ($6=="+") ? $2 : $3
      start = (tss-1000 < 0) ? 0 : tss-1000
      print $1, start, tss+1000, $4, ".", $6
    }' \
  | sort -k1,1 -k2,2n \
  > hg38_tss_2kb.bed
```

---

## 3. Chromosome sizes (`hg38.chrom.sizes`)

```bash
# Using fetchChromSizes (UCSC tools)
fetchChromSizes hg38 > hg38.chrom.sizes

# Or with samtools (if you have a hg38 reference FASTA indexed)
samtools faidx hg38.fa
cut -f1,2 hg38.fa.fai > hg38.chrom.sizes
```

---

## Directory layout after download

```
resources/
├── hg38-blacklist.v2.bed
├── hg38_tss_2kb.bed
├── hg38.chrom.sizes
└── README.md            ← this file
```

---

## Mouse (mm10/mm39) users

Update `config/config.yaml`:
```yaml
genome:      mm10
genome_size: mm
blacklist:   resources/mm10-blacklist.v2.bed
tss_bed:     resources/mm10_tss_2kb.bed
mito_chr:    chrM
```

Download mm10 blacklist:
```bash
wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz
gunzip mm10-blacklist.v2.bed.gz
```
