# my_minion_sequencing



1. Connect to Ramses and enter interactive session.
2. cd /projects/virology/nanopore/output/run_ivi_test/BKV_42
3. ls para ver si estan los bam files: BKV_42_calls_2025-11-10_T16-06-08.bam
4. module load bio/BEDTools/2.31.0-GCC-12.3.0
5. bedtools bamtofastq -i NA18152.bam -fq NA18152.fq # replace with correct names
6. ls -lh BKV_42_calls_2025-11-10_T16-06-08.fq to confirm that the file was created
7. lo guarde accidentalmente en otra carpeta :  scp imarches@ramses4.itcc.uni-koeln.de:/projects/virology/BKV_42_calls_2025-11-10_T16-06-08.fq "C:\Users\Viro_User\Documents\my_minion_sequencing\Output\"
Enter passphrase for key 'C:\Users\Viro_User/.ssh/id_ed25519':
BKV_42_calls_2025-11-10_T16-06-08.fq                                                                   100%   67MB  52.9MB/s   00:01
PS C:\Users\Viro_User>

```bash
cd /projects/virology/nanopore/output/run_ivi_test/BKV_42
ls
BKV_42_calls_2025-11-10_T16-06-08.bam
module load bio/BEDTools/2.31.0-GCC-12.3.0
bedtools bamtofastq \
  -i /projects/virology/nanopore/output/run_ivi_test/BKV_42/BKV_42_calls_2025-11-10_T16-06-08.bam \
  -fq BKV_42_calls_2025-11-10_T16-06-08.fq
 ls -lh BKV_42_calls_2025-11-10_T16-06-08.fq
```

https://bedtools.readthedocs.io/en/latest/content/tools/bamtofastq.html

https://nanoplot.bioinf.be/


## To check available modules


```shell
module avail
```
up and down with J K
quit by pressing q

