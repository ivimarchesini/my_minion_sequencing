#!/bin/bash

# Run this script with ./nanopore_basecalling.sh in the terminal after login to ramses.

# ramses: ssh mpirkl1@ramses4.itcc.uni-koeln.de
# minit: ssh mpirkl1@10.212.1.44

# check job status of current jobs in ramses:
sacct --format="JobID,ReqCPUS,UserCPU,Elapsed,ReqMem,MaxRSS,JobName%20,Stat"

#
# Instructions
# input arguments expected by the script here is the individual user input:

# usage:
# arguments:
# input folder / raw pod5 files
# kit name
# reference file
# 0 or 1, if the patient files are large (default: 1)
# 0 or 1, for cpu computation (default: 0)
# csv of gene regions (start and end in two columns)
# minimal qscore
# minimal phredscore
# example once the script is ready:
# nanopore/nanopore/nanopore_basecalling.sh nanopore/input/BKV_202504/ SQK-NBD114-24 ref_BKV.fasta 1 0 genes_BKV.csv 9 30

# run interactively on up to 4 gpus:
# srun -A virology -p interactive --gpus=4 --time=05:00:00 --mem=24gb -J nanopore.$project.${pat///}.$(date +%F) --pty bash -i

# notes

# check for reference and genes file and check for samplesheet
# make executeable

# Use interactive mode by default
interactive=1

# change into /projects/virology folder  (ramses storage)
cd /projects/virology/

# This is a sample of a variable to set before running the script. This set's the input folder with the pod5 files from minion.
# input=nanopore/input/run_20250716/ # do not edit, but set input=... manually in the terminal

# Here the scripts asks for user input of the input path above. User should write run_20250716 and confirm with RETURN
read -p "Enter project directory: " input

# If inpute was set the input path prepends "nanopore/input/<userinput>"
if [ ! -d "$input" ]; then
  input=nanopore/input/$input
fi

# Check if the directory exists (to check if user input makes sense)
if [ ! -d "$input" ]; then
  echo "Directory does not exists on this path."
  exit 1
fi

# Setting some variables for later use with the $variablename syntax:
kit=SQK-NBD114-24
large=1 # > 10gb?
cpu=0
qscore=9
phredscore=30
min_coverage=500

# Interactive is probably set to one so large is set to 1
if (( interactive == 1 )); then
  large=1
fi

# These are positional arguments to be used like this when running the script ./nanopore_basecalling.sh input kit ref large cpu regions qscore phredscore interactive
# this is used once a script is ready
# input=$1
# kit=$2
# ref=$3
# large=$4
# cpu=$5
# regions=$6
# qscore=$7
# phredscore=$8
# interactive=$9

# the rest is static

# again switching to the virology folder
cd /projects/virology/
project=${input/nanopore\/input\//}
project=${project/\//}

# changing user and group permission of the files in the project to current user and virology group
chown -R $USER:virology nanopore/input/$project
# make all files group read and writeable
chmod -R 775 nanopore/input/$project

# set models path
models=nanopore/models/
# set output path
output=nanopore/output/$project
# again switching to the virology folder

cd /projects/virology/
mkdir nanopore/output/$project
# again switching to the virology folder

cd /projects/virology/
# creating scratch folders (to be deleted after run)
mkdir /scratch/virology
mkdir /scratch/virology/logs

# Load required modules to run dorado (specific to Ramses and SLURM)
module load bio/BEDTools/2.31.0-GCC-12.3.0
module load bio/SAMtools/1.18-GCC-12.3.0
module load bio/BCFtools/1.18-GCC-12.3.0
module load bio/BamTools/2.5.2-GCC-12.3.0
module load bio/VCFtools/0.1.16-GCC-12.3.0
module load bio/minimap2/2.26-GCCcore-12.3.0
# module load bio/BWA/0.7.18-GCCcore-12.3.0
# module load bio/picard/3.0.0-Java-17
# module load bio/MAFFT/7.505-GCC-12.2.0-with-extensions
PATH=$PATH:/home/mpirkl1/micromamba/bin/

# start computation

# samplesheet is expected to be in the input folder with a .csv file ending
samplesheet=$(echo $input*.csv)

# $cpu is set to zero so it should take jump into the else block
if (( $cpu==1 )); then
  large=1
  device="cpu"
  partition="-c 10"
  mem=1024gb
else
  # set devices to any CUDA (GPU) type.
  device="cuda:all"
  # Use the GPU partition with one GPU,
  partition="-p gpu --gpus=1 "
  # Use 16gb of memory for the job (on the node so system memory not GPU)
  mem=16gb # mem=4gb
fi

# set the maximum amout of time before canceling the job.
time=04:00:00 # time=00:30:00

# dorado basecalling and alignment

# large is probably set to 1 (see above)
if (( $large == 1 )); then

  pat0=nanopore/input/$project/barcode13/

  iterate over all subfolders int the project input directory (/projects/virology/nanopore/input/run_20250713/<barcode1-13>)
  for pat0 in ${input}*/; do
    pat=${pat0/$input/}
    # print current folder
    echo $pat
    # Read barcode from samplesheet
    ref=$(Rscript -e "df <- read.csv('$samplesheet');virus<-df[which(df[,'barcode']=='${pat///}'),'virus'];paste0('ref_',virus,'.fasta')")
    ref=${ref//\"}
    ref=${ref//\[}
    ref=${ref//\]}
    ref=${ref/1/}
    ref=${ref/ /}
    # Read id from samplesheet
    id=$(Rscript -e "df <- read.csv('$samplesheet');id<-df[which(df[,'barcode']=='${pat///}'),'id'];id")
    id=${id//\"}
    id=${id//\[}
    id=${id//\]}
    id=${id/1/}
    id=${id/ /}
    # create output folder
    mkdir $output/$id
    # perpare dorado command to execute
    cmd="dorado/bin/dorado basecaller sup "${pat0}" -r --device "$device" -o "${output}/$id" --models-directory "${models}" --min-qscore ${qscore} --kit-name "${kit}" --reference "${ref}
    # print ID from samplesheet
    echo $id

    # interactive probably set to 1 so we expect do be inside an interactive ramses/SLURM session already (Nvidia A30 GPU)
    if (( $interactive == 1 )); then
      $cmd
    else
      sbatch -A virology $partition --mem=$mem -t ${time} -e /scratch/virology/logs/nanopore.$project.${pat///}.$(date +%F).error.txt -o /scratch/virology/logs/nanopore.$project.${pat///}.$(date +%F).output.txt -J nanopore.$project.${pat///}.$(date +%F) --wrap="$cmd"
    fi
  done

else

  bashfile=nanopore.$project.$(date +%F).sh
  rm /scratch/virology/$bashfile
  touch /scratch/virology/$bashfile
  echo '#!/bin/sh' >> /scratch/virology/
  $bashfile
  echo "\$input=$input;\$output=$output;\$samplesheet=$samplesheet;" >> /scratch/virology/$bashfile
  echo "for pat0 in \${input}*/; do pat=\${pat0/\$input/}; mkdir \$output/\$pat; echo \$pat;
  ref=\$(Rscript -e \"df <- read.csv('\$samplesheet');virus<-df[which(df[,'id']=='\${pat///}'),'virus'];paste0('ref_',virus,'.fasta')\");
    ref=\${ref//\\\"};
    ref=\${ref//\[};
    ref=\${ref//\]};
    ref=\${ref/1/};
    ref=\${ref/ /};
    dorado/bin/dorado basecaller sup \${pat0} -r --device cuda:all -o "${output}"/\$pat --models-directory "${models}" --min-qscore '${qscore}' --kit-name "${kit}" --reference "${ref}";
    done" >> /scratch/virology/$bashfile
  chmod +x /scratch/virology/$bashfile
  sbatch -A virology $partition--mem=$mem -t ${time} -e /scratch/virology/logs/nanopore.$project.$(date +%F).error.txt -o /scratch/virology/logs/nanopore.$project.$(date +%F).output.txt -J nanopore.$project.$(date +%F) /scratch/virology/$bashfile

fi

# from here on we do not need any gpus and can exit the interactive session
# if you exit the interactive session, you have to reset all variables and the nskip the dorado part

# make msa/freq files

mem=4gb
time=00:30:00
cores=1

pat0=nanopore/output/$project/barcode03/

for pat0 in ${output}/*/; do
  pat=${pat0/$output/}
  echo $pat
  cmd='bamToFreq/bin/bamToFreq -q '$phredscore' '$pat0*.bam
  if (( $interactive == 1 )); then
    $cmd
  else
    sbatch -A virology -c $cores --mem=$mem -t ${time} -e /scratch/virology/logs/nanopore.$project.${pat///}.$(date +%F).error.txt -o /scratch/virology/logs/nanopore.$project.${pat///}.$(date +%F).output.txt -J nanopore.$project.${pat///}.$(date +%F) --wrap="$cmd"
  fi
done

# consensus files freq files

mem=4gb
time=00:10:00
cores=1

pat0=nanopore/output/$project/barcode01/

for pat0 in ${output}/*/; do
  pat=${pat0/$output/}
  echo $pat
  ref=$(Rscript -e "df <- read.csv('$samplesheet');virus<-df[which(df[,'id']=='${pat///}'),'virus'];paste0('ref_',virus,'.fasta')")
  ref=${ref//\"}
  ref=${ref//\[}
  ref=${ref//\]}
  ref=${ref/1/}
  ref=${ref/ /}
  regions=$(Rscript -e "df <- read.csv('$samplesheet');virus<-df[which(df[,'id']=='${pat///}'),'virus'];paste0('genes_',virus,'.csv')")
  regions=${regions//\"}
  regions=${regions//\[}
  regions=${regions//\]}
  regions=${regions/1/}
  regions=${regions/ /}
  coverage=$(Rscript -e "df <- read.csv('$samplesheet');coverage<-df[which(df[,'id']=='${pat///}'),'coverage'];coverage")
  coverage=${coverage//\"}
  coverage=${coverage//\[}
  coverage=${coverage//\]}
  coverage=${coverage/1/}
  coverage=${coverage/ /}
  cmd='Rscript nanopore/nanopore/nanopore_variants.r '$pat0' '$ref' '$cores' '$regions' '$coverage
  if (( $interactive == 1 )); then
    $cmd
  else
    sbatch -A virology -c $cores --mem=$mem -t ${time} -e /scratch/virology/logs/nanopore.$project.${pat///}.$(date +%F).error.txt -o /scratch/virology/logs/nanopore.$project.${pat///}.$(date +%F).output.txt -J nanopore.$project.${pat///}.$(date +%F) --wrap="$cmd"
  fi
done

# rename files

Rscript -e 'folders <- list.files(paste0("nanopore/output/'$project'/")); for (barcode in folders) {; files <- list.files(paste0("nanopore/output/'$project'/",barcode,"/")); for (file in files) {; system(paste0("mv nanopore/output/'$project'/",barcode,"/",file," nanopore/output/'$project'/",barcode,"/",barcode,"_",file)); }; }'

# rename the fasta headers

Rscript -e 'library(seqinr); folders <- list.files(paste0("nanopore/output/'$project'/")); for (barcode in folders) {; files <- list.files(paste0("nanopore/output/'$project'/",barcode,"/")); for (file in files) {; if (length(grep("\\.fasta",file))>0) {; fasta <- read.fasta(paste0("nanopore/output/'$project'/",barcode,"/",file)); names(fasta) <- paste0(barcode,"_",names(fasta)); write.fasta(fasta,names=names(fasta),file=paste0("nanopore/output/'$project'/",barcode,"/",file)); } } }'

# Setting owner, rights and permission of the output files to virology user group
chown -R $USER:virology nanopore/output/$project # give write access to virology group
chmod -R 775 nanopore/output/$project

# this exit command exits this script with a non-zero exit code indicating that an error orccured although everything looks fine.
exit 1 # THIS IS THE END

# some temp code helper code (DO NOT EXECUTE FROM HERE ON!):

for pod5 in *pod5; do
  echo $pod5
  mkdir ${pod5/.pod5/}
  mv $pod5 ${pod5/.pod5/}/$pod5
done

for pod5 in *; do
  cd $pod5
  mv *.pod5 ../../
  cd ..
  rm -rf $pod5
done

for pod5 in *pod5; do
  echo $pod5
  mv $pod5 run_20250624/$pod5
done

pod5 subset *.pod5 --output subsets --table sequencing_summary_AXO223_12cba927_53d86cc0.txt --columns barcode_arrangement --missing-ok

for pod5 in *pod5; do
  echo $pod5
  mkdir ../${pod5/.pod5/}
  mv $pod5 ../${pod5/.pod5/}/$pod5
done

for dir in *; do
  echo $dir
  mv $dir ${dir/*-/}
done

for dir in *; do
  echo $dir
  mv $dir ../$dir
done

# we still need fasta files (no!)
#
# mem=4gb
# time=00:30:00
# cores=1
#
# for pat0 in ${output}/*/; do
#   pat=${pat0/$output/}
#   done
#   echo $pat
#   pat0=$(echo $pat0*.bam)
#   cmd='samtools fasta -t '$pat0' > '${pat0/bam/fasta}
#   sbatch -A virology -c $cores --mem=$mem -t ${time} -e /scratch/virology/logs/nanopore.$project.${pat///}.$(date +%F).error.txt -o /scratch/virology/logs/nanopore.$project.${pat///}.$(date +%F).output.txt -J nanopore.$project.${pat///}.$(date +%F) --wrap="$cmd"
# done

# fastx pipeline
#
# for pat0 in ${output}*; do
#   pat=${pat0/$output/}
#   # cmd='bedtools bamtofastq -i '${output}/${pat}*.bam' -fq '${output}/$pat'calls.fastq'
#   folder=${output}${pat}/
#   for bam in ${folder}*.bam; do
#     samtools fasta --reference $ref $bam > ${output}${pat}/calls_test.fasta
#     echo $bam
#   done
#   echo $pat
# done

# mem=1gb
# time=00:10:00
#
# for pat0 in ${output}*/; do
#   pat=${pat0/$output/}
#   echo $pat
#   # cmd='bedtools bamtofastq -i '${output}/${pat}*.bam' -fq '${output}/$pat'calls.fastq'
#   cmd='samtools fasta calls_2025-05-16_T07-56-25.bam > calls.fasta'
#   sbatch -A virology -c 1 --mem=$mem -t ${time} -e /scratch/virology/logs/nanopore.$project.${pat///}.$(date +%F).error.txt -o /scratch/virology/logs/nanopore.$project.${pat///}.$(date +%F).output.txt -J nanopore.$project.${pat///}.$(date +%F) --wrap="$cmd"
# done

# old and test stuff here:

for pat0 in ${output}/*/; do
  pat=${pat0/$output/}
  pat=${pat///}
  for file0 in ${output}/$pat/*; do
    file=${file0/$pat0/}
    mv ${file0} ${pat0}${pat}_${file}
  done
done

emacs /scratch/virology/logs/nanopore.run_20250623.barcode10.2025-07-10.error.txt
emacs /scratch/virology/logs/nanopore.$project.${pat///}.$(date +%F).output.txt

sacct --format="JobID,ReqCPUS,UserCPU,Elapsed,ReqMem,MaxRSS,JobName%30,Stat"
sacct --format="JobID,ReqMem,MaxRSS,JobName%50,Stat"

sstat -j 1142239 --format=AveCPU,AveRSS,MaxRSS,MaxRSSTask --allsteps
sacct -j 950724 --format="JobID,ReqCPUS,UserCPU,Elapsed,ReqMem,MaxRSS,JobName%20,State%20"

sstat -j 950838 --format=AveCPU,AveRSS,MaxRSS,MaxRSSTask --allsteps
sacct -j 950838 --format="JobID,ReqCPUS,UserCPU,Elapsed,ReqMem,MaxRSS,JobName%20,State%20"

for i in {1178597..1178620}
  do
  scancel $i
done

# copy stuff

for pat0 in ${output}/*/; do
  pat=${pat0/$output/}
  echo mkdir $pat
  echo scp mpirkl1@ramses4.itcc.uni-koeln.de:$output$pat*.csv .$pat
  echo scp mpirkl1@ramses4.itcc.uni-koeln.de:$output$pat*.fasta .$pat
  echo scp mpirkl1@ramses4.itcc.uni-koeln.de:$output$pat*.pdf .$pat
done

for pat0 in ./*; do
  echo $pat0
  pat=${pat0/.\//}
  echo $pat
  scp $pat/Consensus.fasta Consensus_$pat.fasta
  scp $pat/Consensus_10pct.fasta Consensus_10pct_$pat.fasta
done

# try medaka:

# medaka_variant -i nanopore/output/BKV_202504/barcode01/calls.fasta -r ref_BKV.fasta

#

# bcftools mpileup -d 1000 -Ou -f ref_BKV.fasta nanopore/output/BKV_202504/barcode01/calls_2025-05-16_T07-56-25.bam | bcftools call -vmO z -o nanopore/output/BKV_202504/barcode01/study.vcf.gz
#
# vcftools --gzvcf nanopore/output/BKV_202504/barcode01/study.vcf.gz --freq --out freqfile
#
# bcftools mpileup -d 1000 -Ov -o nanopore/output/BKV_202504/barcode01/study.vcf -f ref_BKV.fasta nanopore/output/BKV_202504/barcode01/calls_2025-05-16_T07-56-25.bam

# samtools: this would work for consensus, but not frequency

# cuts=( 0.99 0.98 0.95 ) # 0.15 0.2 0.3)
#
# echo $(date +%T)
# for pat0 in ${output}/*; do
#   pat=${pat0/$output/}
#   file=${pat0}/*.bam
#   if [ ! -f ${file}.bai ]; then
#     dorado aligner $ref $file -o ${pat0}/
#   fi
#   for cut in ${cuts[@]}; do
#     samtools consensus -A -o ${pat0}/consensus_${cut}.fasta -m "simple" -c $cut $file
#   done
# done
# echo $(date +%T)

# emacs /scratch/virology/logs/nanopore.BKV.2025-05-15.error.txt
# emacs /scratch/virology/logs/nanopore.BKV.align.2025-05-14.output.txt
#
# file=nanopore/output/$project/*.bam
# cmd="bedtools bamtofastq -i "$file" -fq nanopore/output/"$project"/"$(date +%F)".fastq"
#
# sbatch -A virology -c 1 --mem=$mem -t ${time} -e /scratch/virology/logs/nanopore.$project.fastq.$(date +%F).error.txt -o /scratch/virology/logs/nanopore.$project.fastq.$(date +%F).output.txt -J nanopore.$project.fastq.$(date +%F) --wrap="$cmd"
#
# bcftools mpileup -A --threads 2 -f ref_BKV.fasta -o nanopore/output/$project/pileup.vcf -Ov $file
#
# samtools sort -OBAM $file

# install fastq2codfreq

# before you do this try better filtering quality and alignment

# curl -sL https://raw.githubusercontent.com/hivdb/codfreq/main/bin-wrapper/align-all-docker -o fastq2codfreq/
# chmod +x fastq2codfreq
#
# dorado aligner ref_BKV.fasta nanopore/output/BKV/barcode01/calls_2025-05-16_T07-56-25.bam -o nanopore/output/BKV/barcode01/
#
# bcftools mpileup -Ou -f ref_BKV.fasta nanopore/output/BKV/barcode01/calls_2025-05-16_T07-56-25.bam | bcftools call -vmO z -o nanopore/output/BKV/barcode01/study.vcf.gz
#
# bcftools mpileup -Ov -o nanopore/output/BKV/barcode01/study.vcf -f ref_BKV.fasta nanopore/output/BKV/barcode01/calls_2025-05-16_T07-56-25.bam
#
#
# tabix -p vcf nanopore/output/BKV/barcode01/study.vcf.gz
#
# bcftools consensus -f ref_BKV.fasta nanopore/output/BKV/barcode01/study.vcf.gz > nanopore/output/BKV/barcode01/out.fasta
#
#
# bcftools consensus nanopore/output/BKV/barcode01/study.vcf.gz > nanopore/output/BKV/barcode01/out.fasta
#
# vcf <- fread("nanopore/output/BKV/barcode01/study.vcf.gz")
#
# tabix -p vcf nanopore/output/BKV/barcode01/study.vcf.gz
#
# bcftools stats -F <ref.fa> -s - <study.vcf.gz> > <study.vcf.gz.stats>
# mkdir plots
# plot-vcfstats -p plots/ <study.vcf.gz.stats>
