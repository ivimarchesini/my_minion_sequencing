# MyMinion - NGS Pipeline

## Overview
This repository contains an automated analysis pipeline for processing viral whole-genome sequencing data generated using the Oxford Nanopore Technologies (ONT) MinION Mk1C platform.  
The workflow covers all steps from raw signal extraction to consensus generation, variant analysis, and reporting.  
It is optimized for **viral genomes**, where high mutation rates and quasispecies complexity require maximal accuracy.

---

# 1. Laboratory Workflow

## 1.1 Sample Preparation
- Viral target DNA amplified via **nested PCR**.
- Amplicon quality and fragment size verified using **QIAxcel**.
- Library preparation using **ONT SQK-NBD114-24 (Kit14)**:
  - End-repair & dA-tailing  
  - Ligation of barcoded adapters  
  - Sample pooling  
  - Final adapter ligation  
  - Cleanup with AMPure XP beads  
- Prepared library loaded onto the MinION flow cell.

---

# 2. Sequencing Configuration

## 2.1 Barcode Classification
During sequencing, Dorado performs barcode classification in **inline mode**:

- Reads are assigned to barcodes during basecalling.
- Barcode stored in:
  - `BC` tag of each read
  - `RG` (Read Group) in the BAM  

This produces one BAM containing all reads pre-separated by barcode.

---

# 3. Computing Environment

Raw MinION run data are transferred from the Mk1C to the university HPC cluster equipped with **GPU nodes**.  
Basecalling, alignment, QC, and consensus processing occur entirely on this infrastructure.

After computation, all results are synchronized to the diagnostics network storage for downstream clinical interpretation.

---

# 4. Basecalling

## 4.1 Dorado Basecaller
- Uses latest **SUP (Super Accuracy)** model.
- SUP outperforms HAC and Fast models.
- Required due to high viral mutation rates and quasispecies complexity.
- Dorado is actively developed by ONT and maintained on GitHub.

## 4.2 GPU Acceleration
SUP is a **Transformer model** requiring heavy floating-point computation.  
→ GPU mode is used on HPC compute nodes.

---

# 5. Read Filtering

Reads are filtered during basecalling using:

```bash
--min-qscore 9
```
Q = -10 * log10(P)

Where `P` is the basecalling error probability.

Examples:
- **Q9** → approx. 12.5% error  
- **Q30** → later used for strict filtering by bamToFreq  

Reads <Q9 are discarded.

---

# 6. Trimming

With Kit14 (SQK-NBD114-24), Dorado performs trimming automatically unless `--no-trim` is specified.

Trimmed elements:
- ONT adapters  
- Nested PCR primers  
- Native Barcoding (NBD) barcodes  

Read orientation is determined via primer detection (`TS:A:+` / `TS:A:-`).

---

# 7. Alignment

Basecalling and alignment are performed **simultaneously**:

- GPU → basecalling  
- CPU → alignment  

The correct reference (HIV-1, HBV, SARS-CoV-2, etc.) is selected based on the sample sheet.  
A **sorted BAM** is generated for downstream processing.

---

# 8. Quality Control

## 8.1 FASTQ and BAM QC
Reads are extracted from BAM using:

```bash
bedtools bamtofastq
```
Two NanoPlot reports are generated:
1. **FASTQ-based QC**  
2. **BAM-based QC**

NanoPlot outputs include:
- Read length distribution  
- N50  
- Q-score distribution  
- Yield curves  
- Alignment identity  

These metrics determine the suitability for resistance analysis.

---

# 9. Frequency Files

Using **bamToFreq**, two CSV files per sample are generated:

## 9.1 `<sample>_freqs.csv`
Contains:
- Per-base frequencies (A, C, G, T, N, deletions)
- Coverage per position  
→ Supports SNV and minority variant detection.

## 9.2 `<sample>_codonFreqs.csv`
Contains:
- Codon frequencies across the genome  
→ Used for amino-acid–level resistance interpretation.

Strict filtering with:

```bash
-q 30
```
removes low-quality bases.

---

# 10. Alignment Viewer Output

`samttools tview` generates a text-based alignment snapshot (`ngs.txt`) including:
- Reference sequence  
- Stacked reads  
- Mismatches  
- Insertions/deletions  
- Base composition  

This serves as a GUI-independent qualitative review tool.

---

# 11. Mapping Statistics

An R script:
1. Reads genomic region definitions from CSV
2. Counts total reads via:

```bash
bamtools count -in file.bam
```
3. Computes read counts per region  
4. Outputs a table `mapping-tab-date-lab.csv`

This table summarizes genomic coverage per region.

---

# 12. Coverage Plot

The R script:
- Aggregates coverage from `freqs.csv`
- Produces a **PDF coverage plot**

Uses:
- Identifying low- or high-coverage regions  
- Evaluating uniformity of sequencing  
- Documenting data quality for reporting or publications  

---

# 13. Consensus Generation

Consensus sequences are built for thresholds: `1%, 2%, 5%, 10%, 15%, 20%, 30%`


The **10% threshold** is additionally exported as a separate FASTA.

### IUPAC Ambiguity Codes  
Used to encode minority variants.

### Amino Acid Translation
Performed using the R package **Biostrings**:
- Codon frequencies → amino-acid frequencies  
- Redundant codons merged  
- Gene-specific AA tables exported  

Output files include:
- `freqs.csv` (nucleotide frequencies)  
- `freqsC.csv` (codon frequencies)  
- `freqsAA.csv` (amino acid frequencies)  

---

# 14. Pipeline Automation

A Python CLI tool automates the entire workflow.

### Features
- Text-based UI using **textual**
- SSH transfer via **paramiko**
- Automated steps:
  - Run selection (run_YYYYMMDD)
  - Sample sheet selection
  - SCP transfer from Mk1C to HPC
  - GPU basecalling
  - Postprocessing + QC
  - Export to diagnostics directory

### Benefits
- Eliminates manual errors  
- Ensures reproducible processing  
- Saves time for high-throughput scenarios  
- User-friendly interface  

---

# 15. Reporting

All results are exported to a structured output directory:

- QC HTML reports  
- FASTQ files  
- BAM files  
- Coverage plots  
- Frequency tables  
- Consensus FASTA  
- Codon & AA frequencies  
- Alignment snapshots (ngs.txt)

These outputs form the basis for downstream **resistance mutation analysis**.

---

# License
(Add your licensing information here.)

---

# Contact
(Add maintainer or laboratory contact information here.)















































# my_minion_sequencing



1. Connect to Ramses and enter interactive session.
2. cd /projects/virology/nanopore/output/run_ivi_test/BKV_42
3. ls para ver si estan los bam files: BKV_42_calls_2025-11-10_T16-06-08.bam
4. module load bio/BEDTools/2.31.0-GCC-12.3.0
5. bedtools bamtofastq -i NA18152.bam -fq NA18152.fq # replace with correct names
6. ls -lh BKV_42_calls_2025-11-10_T16-06-08.fq to confirm that the file was created
7. cargar el modulo de FastQC module load bio/FastQC/0.11.9-Java-11
8. Convertir el archivo fq en el grafico html de fastqc: BKV_42_calls_2025-11-10_T16-06-08_fastqc.html and BKV_42_calls_2025-11-10_T16-06-08_fastqc.zip
9. en la terminal de la computadora local hacer el scp para bajar el archivo html
10. abrir el report
        
```bash
cd /projects/virology/nanopore/output/run_ivi_test/BKV_42
ls
BKV_42_calls_2025-11-10_T16-06-08.bam
module load bio/BEDTools/2.31.0-GCC-12.3.0
bedtools bamtofastq \
  -i /projects/virology/nanopore/output/run_ivi_test/BKV_42/BKV_42_calls_2025-11-10_T16-06-08.bam \
  -fq /projects/virology/nanopore/output/run_ivi_test/BKV_42/BKV_42_calls_2025-11-10_T16-06-08.fq
ls -lh
module load bio/FastQC/0.11.9-Java-11
fastqc /projects/virology/nanopore/output/run_ivi_test/BKV_42/BKV_42_calls_2025-11-10_T16-06-08.fq \
  -o /projects/virology/nanopore/output/run_ivi_test/BKV_42/
```
scp imarches@ramses4.itcc.uni-koeln.de:/projects/virology/nanopore/output/run_ivi_test/BKV_42/BKV_42_calls_2025-11-10_T16-06-08_fastqc.html "C:\Users\Viro_User\Documents\my_minion_sequencing\Output\BKV_42_calls_2025-11-10_T16-06-08_fastqc.html"
https://bedtools.readthedocs.io/en/latest/content/tools/bamtofastq.html

https://nanoplot.bioinf.be/

https://www.bioconductor.org/packages/release/bioc/html/IONiseR.html

https://www.mdpi.com/2813-9038/2/3/14


## To check available modules


```shell
module avail
```
up and down with J K
quit by pressing q


## List runs with pod5_pass
```shell
ssh minit@10.212.1.44 'find /data/run_* -type d -name "pod5_pass" | sort -u'
```

### TUI helper (select a run)

I added a small Python textual TUI `select_run.py` to connect and list the
available `pod5_pass` folders and let you pick one with the arrow keys.

Install dependencies and run:

```pwsh
pip install textual paramiko
python select_run.py
```

The app will ask for a password (leave blank to try SSH key authentication).
When you press Enter on a selection the chosen path will be printed to stdout
and the TUI will exit.

## Copy results from minion to ramses

```shell
scp -r /data/run_20250803/no_sample_id/20250803_1227_MC-115199_AYL450_752e4334/pod5_pass imarches@ramses4.itcc.uni-koeln.de:/projects/virology/nanopore/input/

```


## Copy results from ramses to agkaiser


```shell
 scp -r imarches@ramses4.itcc.uni-koeln.de:/projects/virology/nanopore/output/run_ivi_test "\\10.212.1.222\agkaiser\reports\Diagnostik\NANOPORE"
```
