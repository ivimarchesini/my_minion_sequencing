# Nanopore Sequencing Pipeline Quality Control

This document summarizes current quality control (QC) coverage in the pipeline (`nanopore_basecalling.sh` + `nanopore_variants.r` logic) and provides recommended extensions to improve early error detection and overall sequencing data quality assessment.

---
## 1. Existing QC Elements

### Bash / Orchestration Stage
- Input directory existence check.
- Sample sheet presence (implicit via glob: `input*.csv`).
- File permission normalization (group-write, ownership).
- Quality thresholds passed to executables: `--min-qscore` (dorado), phred (`-q`) for `bamToFreq`.

### Basecalling (Dorado)
- Minimum Qscore filter (`--min-qscore 9`).
- Kit/model specification (`--kit-name SQK-NBD114-24`).
- Reference supplied for alignment enabling downstream per-base stats.

### Alignment & Frequency Extraction
- Phred score filtering in `bamToFreq` (`-q 30`).
- Generation of nucleotide, codon frequency tables (implicit basis for coverage & variation QC).

### R Variant Analysis (`nanopore_variants.r`)
- Genome-wide coverage computation and polygon plot (`coverage.pdf`).
- Per-gene mapped read counts (`mapping-tab-date-lab.csv`).
- Coverage masking in consensus: positions with depth < `min_coverage` replaced by `N`.
- Consensus sequences at multiple frequency thresholds (1–30%) show ambiguity progression.
- Codon → AA frequency conversion; merged amino acid frequencies highlight functional variability.

### Output Structuring
- Systematic renaming of files for traceability (barcode prefixing).

---
## 2. QC Gaps & Risks

| Area | Gap | Potential Risk |
|------|-----|----------------|
| Sample sheet validation | Column presence & values not validated | Mis-assigned references / barcodes |
| Read-level yield metrics | No N50, length, Qscore distributions | Flowcell/library failure hidden |
| Barcode balance | No per-barcode read/base counts summary | Uneven demultiplexing unnoticed |
| Alignment statistics | Missing `samtools flagstat/stats` summaries | Reference mismatch or contamination undetected |
| Read quality visualization | Only thresholding; no plots | Systematic low-Qscore data not flagged |
| Coverage uniformity | Single global plot only | Gene dropouts masked |
| Variant spectrum | No allele frequency distribution / Tstv ratio | Cannot distinguish signal vs noise |
| Consensus ambiguity summary | Not quantified | High ambiguity (chemistry/model issue) overlooked |
| Contamination check | None | Mixed/incorrect sample passes pipeline |
| Performance metrics | No throughput/GPU utilization | Inefficient resource use |
| Early stop criteria | None | Wasted computation on failed runs |

---
## 3. Recommended QC Additions

### 3.1 Pre-Run Validation
- Verify sample sheet columns: `barcode`, `virus`, `id`, `coverage`.
- Cross-check barcode folders vs sample sheet entries.
- Bar chart: expected coverage (from sheet) vs initial observed (after first ~5k reads if streaming optional).

### 3.2 Basecalling QC
- Collect per-barcode: read count, total bases, mean/median Qscore, read length N50.
- Plots:
  - Read length histogram (log-scaled x).
  - Qscore distribution (violin/box per barcode).
  - Cumulative yield vs time (line).
- Tools: parse FASTQ or use `sequencing_summary.txt` with NanoPlot / pycoQC if available.

### 3.3 Demultiplexing / Barcode Balance
- Table: `barcode_counts.csv` with columns (barcode, read_count, bases, mean_length, mean_qscore).
- Bar plot: reads per barcode with expected uniformity line.
- Flag underrepresented barcodes (< threshold fraction).

### 3.4 Alignment QC
- Commands per BAM:
  - `samtools flagstat` → mapping percentages.
  - `samtools stats` → mismatch/indel rates.
  - `samtools depth` (subset) → depth distribution.
- Plots:
  - MAPQ histogram.
  - Coverage heatmap per gene (position bins).
  - Coverage coefficient of variation per gene (bar chart).

### 3.5 Nucleotide / Variant Frequency QC
- Alt allele frequency spectrum histogram (positions with alt > 1%).
- Transition vs transversion counts (bar plot).
- Shannon entropy per position (line) & per gene (boxplot).
- Ambiguity vs threshold (line: count of non-ACGT bases vs frequency threshold).

### 3.6 Consensus Assessment
- Ambiguous base count and %Ns per gene (bar charts).
- SNP count vs reference per gene (scatter or bar).
- Synonymous vs nonsynonymous vs stop changes (stacked bar).

### 3.7 Amino Acid Variability
- Heatmap: AA variability (entropy) across codon positions for each gene.
- Table: counts of positions with >X% non-WT AA.

### 3.8 Contamination Screening (Optional)
- Classify subset of unmapped reads with Kraken2/Centrifuge.
- Taxonomic composition bar/pie chart.

### 3.9 Performance Metrics
- Basecalling throughput (reads/sec, bases/sec) over time.
- GPU utilization snapshot/time course (`nvidia-smi` query if accessible).

### 3.10 Aggregated Dashboard
- Single HTML/PDF QC report combining summary tables + plots.
- Highlight thresholds with conditional formatting (red/yellow/green).

---
## 4. Example R Snippets for Added QC

### 4.1 Per-Gene Coverage Distribution
```r
genes <- read.csv(regions)
gene_cov <- lapply(1:nrow(genes), function(i) coverage[genes[i, "start"]:genes[i, "stop"]])
pdf(paste0(folder, "coverage_per_gene.pdf"), width=14, height=6)
par(mfrow=c(1,2))
boxplot(gene_cov, names=genes$region, log="y", ylab="Coverage (log)", main="Per-gene coverage")
barplot(sapply(gene_cov, mean), names=genes$region, las=2, log="y", ylab="Mean coverage", main="Mean coverage per gene")
dev.off()
```

### 4.2 Consensus Ambiguity vs Threshold
```r
thresholds <- c(1,2,5,10,15,20,30)
ambigCounts <- sapply(names(seqs), function(nm) sum(grepl("[^ACGT]", seqs[[nm]])))
dfAmbig <- data.frame(threshold=thresholds, ambiguous=ambigCounts)
pdf(paste0(folder, "ambiguity_vs_threshold.pdf"), width=8, height=5)
plot(dfAmbig$threshold, dfAmbig$ambiguous, type="b", pch=19,
     xlab="Frequency threshold (%)", ylab="# Ambiguous positions",
     main="Consensus ambiguity vs threshold")
dev.off()
```

### 4.3 Alt Allele Frequency Spectrum
```r
altFrac <- apply(freqsNorm, 1, function(row) sort(row, decreasing=TRUE)[2])
pdf(paste0(folder, "alt_frequency_spectrum.pdf"), width=8, height=5)
hist(altFrac[altFrac > 0.01], breaks=50, main="Alt allele frequency spectrum (>=1%)",
     xlab="Alt frequency", col="#336699")
dev.off()
```

### 4.4 Synonymous vs Nonsynonymous Summary (Illustrative)
```r
consAA <- aaseqs[[grep("_10%", names(aaseqs), value=TRUE)[1]]]
refDNA <- unlist(strsplit(ref, ""))
impact <- data.frame()
for (j in 1:nrow(genes)) {
  geneDNA <- refDNA[genes[j, "start"]:genes[j, "stop"]]
  refAA <- unlist(strsplit(as.character(Biostrings::translate(
    Biostrings::DNAStringSet(paste(geneDNA, collapse="")))), ""))
  consAAgene <- consAA[1:length(refAA)]
  syn <- sum(consAAgene == refAA)
  nonsyn <- sum(consAAgene != refAA & consAAgene != "*")
  stopc <- sum(consAAgene == "*")
  impact <- rbind(impact, data.frame(region=genes[j, "region"], synonymous=syn,
                                     nonsynonymous=nonsyn, stop=stopc))
}
write.csv(impact, file=paste0(folder, "AA_change_impact.csv"), row.names=FALSE)
```

---
## 5. Suggested Thresholds (Adjustable)

| Metric | Example Alert Threshold | Action |
|--------|-------------------------|--------|
| Mean Qscore | < 10 | Investigate sample/library prep |
| % Mapped Reads | < 70% | Check reference / contamination |
| Mean Coverage (gene) | < `min_coverage` | Consider resequencing / exclude gene |
| % Ns in 10% Consensus | > 5% | Coverage or chemistry issue |
| Barcode Read Fraction | < (expected / 3) | Demultiplexing imbalance |
| Nonsyn/Syn Ratio | Extreme (>3 or <0.1) | Potential error inflation |

---
## 6. Implementation Roadmap
1. Create `QC/` subfolder under each barcode output.
2. Add sample sheet validation (bash or R one-liner) early.
3. Capture dorado per-barcode metrics (parse FASTQ or summary file); store CSV + plots.
4. Generate alignment stats (flagstat / stats) and depth distribution.
5. Extend R script with added coverage, ambiguity, variant spectrum, AA change plots.
6. Aggregate metrics into `QC/report.html` via Rmarkdown.
7. Add threshold-based flags (simple CSV with status column).

---
## 7. Optional Enhancements
- MultiQC integration to aggregate log/stat files.
- Contamination screening with Kraken2 (unmapped reads subset).
- mosdepth for efficient coverage profiling.
- GPU utilization logging if performance audits required.

---
## 8. Summary
The pipeline already enforces base and read quality thresholds, computes coverage, and produces per-gene mapping counts and multi-threshold consensus sequences. By layering standardized read-level, barcode-level, alignment, variant spectrum, ambiguity, and functional impact QC—supported by a concise dashboard—you gain rapid visibility into failures (low yield, poor mapping, uneven coverage, contamination) and can intervene early, saving computational resources and ensuring reliable downstream analyses.

---
**Next Step:** Integrate section 4 snippets into the R workflow and add alignment statistics collection after basecalling.
