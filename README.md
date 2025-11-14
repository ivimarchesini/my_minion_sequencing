# my_minion_sequencing



1. Connect to Ramses and enter interactive session.
2. cd /projects/virology/nanopore/output/run_ivi_test/BKV_42
3. ls para ver si estan los bam files: BKV_42_calls_2025-11-10_T16-06-08.bam
4. module load bio/BEDTools/2.31.0-GCC-12.3.0
5. bedtools bamtofastq -i NA18152.bam -fq NA18152.fq # replace with correct names

```bash
cd /projects/virology/nanopore/output/run_ivi_test/BKV_42
ls
BKV_42_calls_2025-11-10_T16-06-08.bam
module load bio/BEDTools/2.31.0-GCC-12.3.0
bedtools bamtofastq \
  -i /projects/virology/nanopore/output/run_ivi_test/BKV_42/BKV_42_calls_2025-11-10_T16-06-08.bam \
  -fq BKV_42_calls_2025-11-10_T16-06-08.fq
```

https://bedtools.readthedocs.io/en/latest/content/tools/bamtofastq.html

https://nanoplot.bioinf.be/


## To check available modules


```shell
module avail
```
up and down with J K
quit by pressing q

