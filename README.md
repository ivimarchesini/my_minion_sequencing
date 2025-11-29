# NGS Pipeline – README

## Overview
This repository contains an automated analysis pipeline for processing viral whole-genome sequencing data generated using the Oxford Nanopore Technologies (ONT) MinION Mk1C platform.  
The workflow covers all steps from raw signal extraction to consensus generation, variant analysis, and reporting.  
It is optimized for **viral genomes**, where high mutation rates and quasispecies complexity require maximal accuracy.

The pipeline complies with common diagnostic guidelines (WHO, ECDC, CDC, national societies) regarding sequencing quality, genome completeness, minimum coverage, and base quality metrics.

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
