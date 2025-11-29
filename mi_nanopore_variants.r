# Dependencies and argument parsing: parses command line arguments: input/output folders, reference file, number of CPU cores, regions file, minimum coverage threshold.
# folder <- "nanopore/output/BKV_1/barcode01/"; reffile <- "ref_BKV.fasta"; cores <- 1; regions <- "genes_BKV.csv"; min_coverage <- 500

# For each user, it's necessary to have the following package versions:
# > packageVersion("MASS")
# [1] ‘7.3.60’
# > packageVersion("seqinr")
# [1] ‘4.2.36’
# > packageVersion("Biostrings")
# [1] ‘2.70.3’

# library load
#library(MASS)
library(seqinr)
library(Biostrings)

# Command line arguments: We are calling the commandArgs() function to retrieve the command-line arguments that are provided when the R script is executed from a terminal or shell,
# and assigning them to a character vector variable named args
args <- commandArgs(TRUE)
folder <- args[1]
reffile <- args[2]
cores <- as.numeric(args[3])
regions <- args[4]
min_coverage <- as.numeric(args[5])

print(folder)
print(reffile)
print(cores)
print(regions)
print(min_coverage)

# Parallel Setup: If more than one core is specified, it creates a fork cluster for potential parallel processing
if (cores>1) {
  try(parallel::stopCluster(cl))
  cl <- parallel::makeForkCluster(cores)
}

# Reference and data input
# Reference Sequence:
# Read the FASTA file into a list of sequences. (seqinr package)
# Extract the name of the first sequence for reference.
# Convert the first sequence into a single lowercase string (nucleotide letters) for downstream analysis.
ref <- read.fasta(reffile)
refname <- names(ref)[1]
ref <- tolower(paste(ref[[1]],collapse=""))

# List all files in a folder.
# Search for files for names starting with "calls" and ending with "freqs.csv".Freqs is a data frame of single-nucleotide frequencies.
# Search for files ending with "codonFreqs.csv". FreqsC will hold codon frequencies instead of single-nucleotide frequencies.
# This subsets the freqs data frame to only include the columns "pos","A","C","G","T","N".
# It overwrites the original freqs CSV with the filtered columns only.
files <- list.files(folder)
freqs <- read.csv(paste0(folder,files[grep("^calls.*freqs.csv",files)]))
freqsC <- read.csv(paste0(folder,files[grep("^calls.*codonFreqs.csv",files)]))
freqs <- freqs[,c("pos","A","C","G","T","N")]
write.csv(freqs,file=paste0(folder,files[grep("^calls.*freqs.csv",files)]),row.names=FALSE)

# Finding all BAM files in the folder and storing their names in the variable "bamfile"
bamfile <- files[grep("bam$",files)]

# NGS FASTA generation
# system() runs a shell command from R.
# samtools tview –> Text-based viewer for BAM alignments. Output: ngs.txt which contains a text representation of the aligned reads
system(paste0("samtools tview -d T -w ",nchar(ref)*2," ",folder,bamfile," ",reffile," > ",folder,"ngs.txt"))

# Reads the ngs.txt file into R as a data frame.
# read.delim() assumes tab-delimited text.
ngs <- read.delim(paste0(folder,"ngs.txt"))

# Renames the single column to "X1"
colnames(ngs) <- "X1"

# ngs[[1]] –> Extracts the first (and only) column of ngs as a vector.
# as.vector() –> Ensures it’s a plain vector of characters.
# as.list() –> Converts the vector into a list, which write.fasta() expects.
# names=1:nrow(ngs) –> Gives each sequence a unique name (1, 2, 3, …).
# file.out=... –> Specifies the output FASTA file (ngs.fasta).
# as.string=TRUE –> Treats each element of the list as a single sequence string (rather than splitting into individual letters).
# Result: The ngs.txt content is now saved as a FASTA file, with each line becoming a separate sequence entry.
write.fasta(as.list(as.vector(ngs[[1]])),names=1:nrow(ngs),file.out=paste0(folder,"ngs.fasta"),
            as.string=TRUE)

# Gen regions and mapping
# Read gene regions and create a data frame
genes <- read.csv(regions)

# Counting total reads in a BAM file using bamtools. "readstotal" contains the total number of reads in the BAM file.
readstotal <- as.numeric(system(paste0("bamtools count -in ",folder,bamfile),
                                intern=TRUE))

# Create a empty data frame for results."mapping" is a results table where each row corresponds to a gene/region.Columns:
#-LAB and PNAME are empty strings (lab/sample).
#-TARGET is filled with the region names from genes$region.
#-MAPPED will store the number of reads mapped to that region (initialized to 0).
#-READS-TOTAL is the same for all rows, showing the total reads in the BAM.
#FASTA is filled with "Y" (possibly to indicate that a reference sequence exists).
mapping <- data.frame(LAB=rep("",nrow(genes)),PNAME=rep("",nrow(genes)),
                      TARGET=genes[,"region"],MAPPED=rep(0,nrow(genes)),
                      'READS-TOTAL'=rep(readstotal,nrow(genes)),
                      FASTA=rep("Y",nrow(genes)))

# Counting reads per region: Loop over each gene/region in the genes data frame and builds a command to count reads only in the specified region.
# Handles cases where no reads are found
# Stores the number of reads mapped to that specific region.
for (i in 1:nrow(genes)) {
  mapped <- as.numeric(system(paste0("bamtools count -in ",folder,bamfile," -region ",refname,":",genes[i,"start"],"-",genes[i,"stop"]),intern=TRUE))
  if (length(mapped)==0) mapped <- 0
  mapping$MAPPED[i] <- mapped
}

# Saving the results: Writes the mapping data frame to a "mapping-tab-date-lab.csv" file.
write.csv(mapping,file=paste0(folder,"mapping-tab-date-lab.csv"),row.names=FALSE)

# Coverage calculation and visualisation 
# Calculating coverage: remove the first column and sums across rows, giving total reads at each genomic position.
# Coverage is a vector where each element is the total number of reads mapped to that position.
coverage <- apply(freqs[,-1],1,sum)

# Open a PDF file to save the plot.
pdf(paste0(folder,"coverage.pdf"),width=40,height=10)

# Creating an empty plot
plot(1:length(coverage),rep(max(log10(coverage+1)),length(coverage)),
     xlab="Genome position",yaxt="n",ylab="Coverage (reads per position)",
     col=rgb(0,0,0,0),ylim=c(0,max(log10(coverage+1))))

# Drawing the coverage as a polygon. It draws a filled shape on the plot.
polygon(c(1:length(coverage),length(coverage):1),c(log10(coverage+1),rep(0,length(coverage))),
        col=rgb(0.5,0.5,0,0.25),border=rgb(0.5,0.5,0,0.75))
axis(2,0:10,c(0,10^(1:10)))

# Finalize and save the plot to coverage.pdf.
dev.off()

# Confirmation: Print a message to the console to indicate the plot is done.
print("Finished coverage plot")

# Consensus Sequence Construction 
# Defining IUPAC ambiguity codes and mapping nucleotide combinations to IUPAC codes:
# Single bases → themselves (a, c, g, t)
# Multiple possibilities → standard ambiguity codes (r=A/G, y=C/T, etc.)
# acgt → n (any base)
# Used later to summarize variable positions in a consensus sequence
iupac <- function(x) {
  y <- switch(x,
              a="a",c="c",g="g",t="t",
              ag="r",ct="y",gt="k",ac="m",cg="s",at="w",
              cgt="b",agt="d",act="h",acg="v",
              acgt="n")
  return(y)
}

# Defining aminoacids
# Stores all amino acid symbols including * for stop codons.
AAs <- c("A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V","*")

# Loop over each gene
for (j in 1:nrow(genes)) {
  # keep consensus based on coverage 500+
  gene <- genes[j,"start"]:genes[j,"stop"]
  # Normalizing nucleotide frequencies
  freqsNorm <- t(apply(freqs[gene,-1],1,function(x) return(x/sum(x))))
  freqsNorm[is.na(freqsNorm)] <- 0
  #Constructing consensus sequences at multiple thresholds
  seqs <- aaseqs <- list()
  for (i in c(1,2,5,10,15,20,30)) {
    seq <- apply(freqsNorm,1,function(x) {
      tmp <- paste(tolower(colnames(freqsNorm)[x>(i/100)]),collapse="")
      tmp <- iupac(tmp)
      if (is.null(tmp)) tmp <- "x"
      return(tmp)
    })
    
  # Mask low coverage positions: Replaces positions with coverage below min_coverage with n (unknown nucleotide).
    seq[coverage[gene]<min_coverage] <- "n"
    
  # Stores consensus DNA sequence in a list, named by gene and threshold.
    seqs[[paste0(genes[j,1],"_",i,"%")]] <- toupper(seq)
    
  # Translate to amino acids
    aaseq <- unlist(strsplit(as.character(Biostrings::translate(Biostrings::DNAStringSet(paste(toupper(gsub("x","n",seq)),collapse="")),if.fuzzy.codon="solve",no.init.codon=TRUE)),""))
    aaseqs[[paste0(genes[j,1],"_",i,"%")]] <- aaseq

  # Write FASTA files: Saves 10% threshold consensus sequences individually and all thresholds together in a single FASTA file.
    if (i==10) {
      write.fasta(toupper(seq),file=paste0(folder,genes[j,"region"],"_consensus_10pct.fasta"),names=names(seqs),nbchar=80)
    }
  }
  write.fasta(seqs,file=paste0(folder,genes[j,"region"],"_consensus.fasta"),names=names(seqs),nbchar=80)
}

# Genome-wide consensus: After gene-wise loops, it repeats similar steps for the whole genome.
gene <- 1:nrow(freqs)
freqsNorm <- t(apply(freqs[gene,-1],1,function(x) return(x/sum(x))))
freqsNorm[is.na(freqsNorm)] <- 0
seqs <- aaseqs <- list()
for (i in c(1,2,5,10,15,20,30)) {
  seq <- apply(freqsNorm,1,function(x) {
    tmp <- paste(tolower(colnames(freqsNorm)[x>(i/100)]),collapse="")
    tmp <- iupac(tmp)
    if (is.null(tmp)) tmp <- "x"
    return(tmp)
  })
  seq[coverage[gene]<min_coverage] <- "n"
  seqs[[paste0("Consensus_",i,"%")]] <- toupper(seq)
  aaseq <- unlist(strsplit(as.character(Biostrings::translate(Biostrings::DNAStringSet(paste(toupper(gsub("x","n",seq)),collapse="")),if.fuzzy.codon="solve",no.init.codon=TRUE)),""))
  aaseqs[[paste0(genes[j,1],"_",i,"%")]] <- aaseq
  if (i==10) {
    write.fasta(toupper(seq),file=paste0(folder,"Consensus_10pct.fasta"),names="Consensus_10%",nbchar=80)
  }
}
write.fasta(seqs,file=paste0(folder,"Consensus.fasta"),names=names(seqs),nbchar=80)

# Confirmation: Print a message to the console to indicate the consensus is done.
print("Finished consensus")

# Aminoacid frequency files
# Preparing AA frequency table.This converts codons → amino acids using Biostrings::translate:
#-"x" is replaced with "n" for unknown codons.
#-then concatenates all codons, translates to protein, and splits into individual amino acids.
#Result: Column names are now amino acids corresponding to codons.
freqsAA <- freqsC
colnames(freqsAA)[3:(ncol(freqsAA)-1)] <- unlist(strsplit(as.character(Biostrings::translate(Biostrings::DNAStringSet(paste(toupper(gsub("x","n",colnames(freqsAA)[3:(ncol(freqsAA)-1)])),collapse="")),if.fuzzy.codon="solve",no.init.codon=TRUE)),""))

# Convert all values to numeric
for (i in 1:ncol(freqsAA)) {
freqsAA[,i] <- as.numeric(freqsAA[,i])
}
                     
# Merge columns with the same AA: Multiple codons may encode the same amino acid.Sums frequencies of codons encoding the same AA and keeps only one column per AA.
for (i in unique(colnames(freqsAA))) {
  idx <- which(colnames(freqsAA)==i)
  freqsAA[,idx[1]] <- apply(freqsAA[,idx,drop=FALSE],1,sum)
}
freqsAA <- freqsAA[,!duplicated(colnames(freqsAA))]

# Loop through genes: Defines nucleotide positions of the gene and converts to codon positions using /3.
for (j in 1:nrow(genes)) {
  gene <- genes[j,"start"]:genes[j,"stop"]
  startpos <- ceiling(genes[j,"start"]/3)
  endpos <- ceiling(genes[j,"stop"]/3)

# Subset AA frequencies for gene: finds codon positions in freqsAA that match the gene start and end.Skips the gene if no positions are found.
  freqsAAtmp <- freqsAA[!is.na(freqsAA$pos),]
  startposfound <- intersect(which.min(abs(freqsAAtmp$pos-startpos)),
                             which(freqsAAtmp$pos>=startpos))
  endposfound <- intersect(which.min(abs(freqsAAtmp$pos-endpos)),
                             which(freqsAAtmp$pos<=endpos))
  if (length(startposfound)==0 | length(endposfound)==0) next()

# Determine wild-type sequence: translates the reference DNA sequence (ref) for the gene to amino acids and adds "*" at the end to indicate stop codon.
  wt <- c(unlist(strsplit(as.character(Biostrings::translate(Biostrings::DNAStringSet(paste(unlist(strsplit(ref,""))[gene],collapse="")),if.fuzzy.codon="solve")),"")),"*")
 
# Adjust positions: aligns wild-type AA sequence with the positions in the AA frequency table.Renames column "X." to "*" (stop codon).
  wtstart <- freqsAAtmp$pos[startposfound]-startpos+1
  wtend <- endpos-freqsAAtmp$pos[endposfound]+length(wt)
  startpos <- freqsAAtmp$pos[startposfound]
  endpos <- freqsAAtmp$pos[endposfound]
  freqsAAtmp <- freqsAAtmp[which(freqsAAtmp$pos==startpos):which(freqsAAtmp$pos==endpos),]
  freqsAAtmp <- data.frame(pos=1:nrow(freqsAAtmp),WT=wt[wtstart:(min(wtend,nrow(freqsAAtmp))+wtstart-1)],freqsAAtmp)
  colnames(freqsAAtmp)[colnames(freqsAAtmp)=="X."] <- "*"
  
# Save AA frequency CSV: Writes the gene-specific AA frequency table to a CSV
  write.csv(freqsAAtmp,row.names=FALSE,
            file=paste0(folder,paste0(genes[j,"region"],"_",gsub("codon","AA",files[grep("^calls.*codonFreqs.csv",files)]))))
  
# Save codon frequency CSV:Subsets codon frequency table for the gene.
  freqsCtmp <- freqsC[which(freqsC$pos==startpos):which(freqsC$pos==endpos),]
  write.csv(freqsCtmp,row.names=FALSE,
            file=paste0(folder,paste0(genes[j,"region"],"_",files[grep("^calls.*codonFreqs.csv",files)])))
  
# Save nucleotide frequency CSV:Subsets nucleotide frequency table for the gene.
  startpos <- intersect(which.min(abs(freqs$pos-genes[j,"start"])),
                             which(freqs$pos>=genes[j,"start"]))
  endpos <- intersect(which.min(abs(freqs$pos-genes[j,"stop"])),
                           which(freqs$pos<=genes[j,"stop"]))
  freqstmp <- freqs[startpos:endpos,]
  write.csv(freqstmp,row.names=FALSE,
            file=paste0(folder,paste0(genes[j,"region"],"_",files[grep("^calls.*freqs.csv",files)])))
}

# End of the pipeline
# Quits R without saving the workspace.
# Returns status=0 (success)
q(save="no",status=0)
