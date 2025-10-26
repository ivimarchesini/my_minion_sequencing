#!/bin/bash

#####
## Run this script with ./nanopore_basecalling.sh in the terminal 
## after login to ramses and starting an interactive session.
#####

# ramses: ssh mpirkl1@ramses4.itcc.uni-koeln.de
# minit: ssh mpirkl1@10.212.1.44

# check job status of current jobs in ramses:
sacct --format="JobID,ReqCPUS,UserCPU,Elapsed,ReqMem,MaxRSS,JobName%20,Stat"

## Use interactive mode by default
interactive=1

## change into /projects/virology folder  (ramses storage)
cd /projects/virology/

## Here the scripts asks for user input of the input path above. User should write run_20250716 and confirm with RETURN
read -p "Enter project directory: " input

## If input was set the input path prepends "nanopore/input/<userinput>"
if [ ! -d "$input" ]; then
  input=nanopore/input/$input
fi

## Check if the directory exists (to check if user input makes sense)
if [ ! -d "$input" ]; then
  echo "Directory does not exists on this path."
  exit 1
fi

## Setting some variables for later use with the $variablename syntax:
kit=SQK-NBD114-24
large=1 # > 10gb?
cpu=0
qscore=9
phredscore=30
min_coverage=500

## make sure we are still in the /project/virology folder
cd /projects/virology/
## get project name from the current folder and variable (Clemens: not sure how this works)
project=${input/nanopore\/input\//}
project=${project/\//}

## changing user and group permission of the files in the project to current user and virology group
chown -R $USER:virology nanopore/input/$project
## make all files group read and writeable
chmod -R 775 nanopore/input/$project

## set models path
models=nanopore/models/
## set output path
output=nanopore/output/$project

## make sure we are still in the /project/virology folder
cd /projects/virology/
mkdir nanopore/output/$project

## make sure we are still in the /project/virology folder
cd /projects/virology/
## creating scratch folders (to be deleted after run) for temporary output
mkdir /scratch/virology
mkdir /scratch/virology/logs

## Load required modules to run dorado (specific to Ramses and SLURM)
module load bio/BEDTools/2.31.0-GCC-12.3.0
module load bio/SAMtools/1.18-GCC-12.3.0
module load bio/BCFtools/1.18-GCC-12.3.0
module load bio/BamTools/2.5.2-GCC-12.3.0
module load bio/VCFtools/0.1.16-GCC-12.3.0
module load bio/minimap2/2.26-GCCcore-12.3.0
# module load bio/BWA/0.7.18-GCCcore-12.3.0
# module load bio/picard/3.0.0-Java-17
# module load bio/MAFFT/7.505-GCC-12.2.0-with-extensions

## Expand the PATH to tell the terminal where to find some executables.
PATH=$PATH:/home/mpirkl1/micromamba/bin/


## samplesheet is expected to be in the input folder with a .csv file ending
samplesheet=$(echo $input*.csv)

## set devices to any CUDA (GPU) type.
device="cuda:all"
# Use the GPU partition with one GPU,
partition="-p gpu --gpus=1 "
## Use 16gb of memory for the job (on the node so system memory not GPU)
mem=16gb # mem=4gb
## set the maximum amout of time before canceling the job.
time=04:00:00 # time=00:30:00

# dorado basecalling and alignment

## iterate over all subfolders int the project input directory (/projects/virology/nanopore/input/run_20250713/<barcode1-13>)
for pat0 in ${input}*/; do
  pat=${pat0/$input/}
  ## print current folder
  echo $pat
  ## Read barcode from samplesheet then some lines of cleaning up the output of the Rscript
  ref=$(Rscript -e "df <- read.csv('$samplesheet');virus<-df[which(df[,'barcode']=='${pat///}'),'virus'];paste0('ref_',virus,'.fasta')")
  ref=${ref//\"}
  ref=${ref//\[}
  ref=${ref//\]}
  ref=${ref/1/}
  ref=${ref/ /}
  ## Read ID from samplesheet then some lines of cleaning up the output of the Rscript
  id=$(Rscript -e "df <- read.csv('$samplesheet');id<-df[which(df[,'barcode']=='${pat///}'),'id'];id")
  id=${id//\"}
  id=${id//\[}
  id=${id//\]}
  id=${id/1/}
  id=${id/ /}
  ## create output folder
  mkdir $output/$id
  ## perpare dorado command to execute
  cmd="dorado/bin/dorado basecaller sup "${pat0}" -r --device "$device" -o "${output}/$id" --models-directory "${models}" --min-qscore ${qscore} --kit-name "${kit}" --reference "${ref}
  ## print ID from samplesheet
  echo $id

  ## interactive probably set to 1 so we expect do be inside an interactive ramses/SLURM session already (Nvidia A30 GPU)
  if (( $interactive == 1 )); then
    $cmd
  else
    sbatch -A virology $partition --mem=$mem -t ${time} -e /scratch/virology/logs/nanopore.$project.${pat///}.$(date +%F).error.txt -o /scratch/virology/logs/nanopore.$project.${pat///}.$(date +%F).output.txt -J nanopore.$project.${pat///}.$(date +%F) --wrap="$cmd"
  fi
done

## Here the script should wait for the the sbatch job to complete otherwise the second part would not work.


## 
# from here on we do not need any gpus and can exit the interactive session
# if you exit the interactive session, you have to reset all variables and the nskip the dorado part

# make msa/freq files

mem=4gb
time=00:30:00
cores=1

## sample how a value of pat0 might look like. 
## Overwritten on `for pat0 in ${output}/*/; do` 3 lines below
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
  ## Read virus from samplesheet then some lines of cleaning up the output of the Rscript
  ref=$(Rscript -e "df <- read.csv('$samplesheet');virus<-df[which(df[,'id']=='${pat///}'),'virus'];paste0('ref_',virus,'.fasta')")
  ref=${ref//\"}
  ref=${ref//\[}
  ref=${ref//\]}
  ref=${ref/1/}
  ref=${ref/ /}
  ## Read regions from samplesheet then some lines of cleaning up the output of the Rscript
  regions=$(Rscript -e "df <- read.csv('$samplesheet');virus<-df[which(df[,'id']=='${pat///}'),'virus'];paste0('genes_',virus,'.csv')")
  regions=${regions//\"}
  regions=${regions//\[}
  regions=${regions//\]}
  regions=${regions/1/}
  regions=${regions/ /}
  ## Read coverage from samplesheet then some lines of cleaning up the output of the Rscript
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

## rename files in output folders
Rscript -e 'folders <- list.files(paste0("nanopore/output/'$project'/")); for (barcode in folders) {; files <- list.files(paste0("nanopore/output/'$project'/",barcode,"/")); for (file in files) {; system(paste0("mv nanopore/output/'$project'/",barcode,"/",file," nanopore/output/'$project'/",barcode,"/",barcode,"_",file)); }; }'

# rename the fasta headers
Rscript -e 'library(seqinr); folders <- list.files(paste0("nanopore/output/'$project'/")); for (barcode in folders) {; files <- list.files(paste0("nanopore/output/'$project'/",barcode,"/")); for (file in files) {; if (length(grep("\\.fasta",file))>0) {; fasta <- read.fasta(paste0("nanopore/output/'$project'/",barcode,"/",file)); names(fasta) <- paste0(barcode,"_",names(fasta)); write.fasta(fasta,names=names(fasta),file=paste0("nanopore/output/'$project'/",barcode,"/",file)); } } }'

## Setting owner, rights and permission of the output files to virology user group
chown -R $USER:virology nanopore/output/$project # give write access to virology group
chmod -R 775 nanopore/output/$project

