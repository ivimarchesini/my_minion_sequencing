# my_minion_sequencing



1. Connect to Ramses and enter interactive session.
2. cd /projects/virology/nanopore/output/ivi_tests_run
3. ls para ver si estan los bam files
4. module load bio/BEDTools/2.31.0-GCC-12.3.0
5. bedtools bamtofastq -i NA18152.bam -fq NA18152.fq # replace with correct names

```bash
cd /projects/virology/nanopore/output/ivi_tests_run
ls
module load bio/BEDTools/2.31.0-GCC-12.3.0
bedtools bamtofastq -i NA18152.bam -fq NA18152.fq
```

https://bedtools.readthedocs.io/en/latest/content/tools/bamtofastq.html
