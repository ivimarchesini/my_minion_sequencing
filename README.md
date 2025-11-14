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


## To check available modules


```shell
module avail
```
up and down with J K
quit by pressing q

