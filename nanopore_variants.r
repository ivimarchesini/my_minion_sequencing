# get arguements
# folder <- "nanopore/output/BKV_1/barcode01/"; reffile <- "ref_BKV.fasta"; cores <- 1; regions <- "genes_BKV.csv"; min_coverage <- 500

library(seqinr)
library(Biostrings)

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

if (cores>1) {
  try(parallel::stopCluster(cl))
  cl <- parallel::makeForkCluster(cores)
}

ref <- read.fasta(reffile)
refname <- names(ref)[1]
ref <- tolower(paste(ref[[1]],collapse=""))

files <- list.files(folder)
freqs <- read.csv(paste0(folder,files[grep("^calls.*freqs.csv",files)]))
freqsC <- read.csv(paste0(folder,files[grep("^calls.*codonFreqs.csv",files)]))
freqs <- freqs[,c("pos","A","C","G","T","N")]
write.csv(freqs,file=paste0(folder,files[grep("^calls.*freqs.csv",files)]),row.names=FALSE)

bamfile <- files[grep("bam$",files)]

# make ngs fasta

system(paste0("samtools tview -d T -w ",nchar(ref)*2," ",folder,bamfile," ",reffile," > ",folder,"ngs.txt"))
ngs <- read.delim(paste0(folder,"ngs.txt"))
colnames(ngs) <- "X1"
write.fasta(as.list(as.vector(ngs[[1]])),names=1:nrow(ngs),file.out=paste0(folder,"ngs.fasta"),
            as.string=TRUE)

# rest works well

genes <- read.csv(regions)

readstotal <- as.numeric(system(paste0("bamtools count -in ",folder,bamfile),
                                intern=TRUE))
mapping <- data.frame(LAB=rep("",nrow(genes)),PNAME=rep("",nrow(genes)),
                      TARGET=genes[,"region"],MAPPED=rep(0,nrow(genes)),
                      'READS-TOTAL'=rep(readstotal,nrow(genes)),
                      FASTA=rep("Y",nrow(genes)))
for (i in 1:nrow(genes)) {
  mapped <- as.numeric(system(paste0("bamtools count -in ",folder,bamfile," -region ",refname,":",genes[i,"start"],"-",genes[i,"stop"]),intern=TRUE))
  if (length(mapped)==0) mapped <- 0
  mapping$MAPPED[i] <- mapped
}
write.csv(mapping,file=paste0(folder,"mapping-tab-date-lab.csv"),row.names=FALSE)

coverage <- apply(freqs[,-1],1,sum)
pdf(paste0(folder,"coverage.pdf"),width=40,height=10)
plot(1:length(coverage),rep(max(log10(coverage+1)),length(coverage)),
     xlab="Genome position",yaxt="n",ylab="Coverage (reads per position)",
     col=rgb(0,0,0,0),ylim=c(0,max(log10(coverage+1))))
polygon(c(1:length(coverage),length(coverage):1),c(log10(coverage+1),rep(0,length(coverage))),
        col=rgb(0.5,0.5,0,0.25),border=rgb(0.5,0.5,0,0.75))
axis(2,0:10,c(0,10^(1:10)))
dev.off()
print("Finished coverage plot")

# do mapped reads instead of per position coverage for gene regions

iupac <- function(x) {
  y <- switch(x,
              a="a",c="c",g="g",t="t",
              ag="r",ct="y",gt="k",ac="m",cg="s",at="w",
              cgt="b",agt="d",act="h",acg="v",
              acgt="n")
  return(y)
}

AAs <- c("A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V","*")

for (j in 1:nrow(genes)) {
  # keep consensus based on coverage 500+
  gene <- genes[j,"start"]:genes[j,"stop"]
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
    seqs[[paste0(genes[j,1],"_",i,"%")]] <- toupper(seq)
    aaseq <- unlist(strsplit(as.character(Biostrings::translate(Biostrings::DNAStringSet(paste(toupper(gsub("x","n",seq)),collapse="")),if.fuzzy.codon="solve",no.init.codon=TRUE)),""))
    aaseqs[[paste0(genes[j,1],"_",i,"%")]] <- aaseq
    if (i==10) {
      write.fasta(toupper(seq),file=paste0(folder,genes[j,"region"],"_consensus_10pct.fasta"),names=names(seqs),nbchar=80)
    }
  }
  write.fasta(seqs,file=paste0(folder,genes[j,"region"],"_consensus.fasta"),names=names(seqs),nbchar=80)
}

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

print("Finished consensus")

# try to produce aa freq file

freqsAA <- freqsC
colnames(freqsAA)[3:(ncol(freqsAA)-1)] <- unlist(strsplit(as.character(Biostrings::translate(Biostrings::DNAStringSet(paste(toupper(gsub("x","n",colnames(freqsAA)[3:(ncol(freqsAA)-1)])),collapse="")),if.fuzzy.codon="solve",no.init.codon=TRUE)),""))
for (i in 1:ncol(freqsAA)) {
  freqsAA[,i] <- as.numeric(freqsAA[,i])
}
for (i in unique(colnames(freqsAA))) {
  idx <- which(colnames(freqsAA)==i)
  freqsAA[,idx[1]] <- apply(freqsAA[,idx,drop=FALSE],1,sum)
}
freqsAA <- freqsAA[,!duplicated(colnames(freqsAA))]

for (j in 1:nrow(genes)) {
  gene <- genes[j,"start"]:genes[j,"stop"]
  startpos <- ceiling(genes[j,"start"]/3)
  endpos <- ceiling(genes[j,"stop"]/3)
  freqsAAtmp <- freqsAA[!is.na(freqsAA$pos),]
  startposfound <- intersect(which.min(abs(freqsAAtmp$pos-startpos)),
                             which(freqsAAtmp$pos>=startpos))
  endposfound <- intersect(which.min(abs(freqsAAtmp$pos-endpos)),
                             which(freqsAAtmp$pos<=endpos))
  if (length(startposfound)==0 | length(endposfound)==0) next()
  wt <- c(unlist(strsplit(as.character(Biostrings::translate(Biostrings::DNAStringSet(paste(unlist(strsplit(ref,""))[gene],collapse="")),if.fuzzy.codon="solve")),"")),"*")
  wtstart <- freqsAAtmp$pos[startposfound]-startpos+1
  wtend <- endpos-freqsAAtmp$pos[endposfound]+length(wt)
  startpos <- freqsAAtmp$pos[startposfound]
  endpos <- freqsAAtmp$pos[endposfound]
  freqsAAtmp <- freqsAAtmp[which(freqsAAtmp$pos==startpos):which(freqsAAtmp$pos==endpos),]
  freqsAAtmp <- data.frame(pos=1:nrow(freqsAAtmp),WT=wt[wtstart:(min(wtend,nrow(freqsAAtmp))+wtstart-1)],freqsAAtmp)
  colnames(freqsAAtmp)[colnames(freqsAAtmp)=="X."] <- "*"
  write.csv(freqsAAtmp,row.names=FALSE,
            file=paste0(folder,paste0(genes[j,"region"],"_",gsub("codon","AA",files[grep("^calls.*codonFreqs.csv",files)]))))
  freqsCtmp <- freqsC[which(freqsC$pos==startpos):which(freqsC$pos==endpos),]
  write.csv(freqsCtmp,row.names=FALSE,
            file=paste0(folder,paste0(genes[j,"region"],"_",files[grep("^calls.*codonFreqs.csv",files)])))
  startpos <- intersect(which.min(abs(freqs$pos-genes[j,"start"])),
                             which(freqs$pos>=genes[j,"start"]))
  endpos <- intersect(which.min(abs(freqs$pos-genes[j,"stop"])),
                           which(freqs$pos<=genes[j,"stop"]))
  freqstmp <- freqs[startpos:endpos,]
  write.csv(freqstmp,row.names=FALSE,
            file=paste0(folder,paste0(genes[j,"region"],"_",files[grep("^calls.*freqs.csv",files)])))
}

q(save="no",status=0)

# old alternate pipeline

files <- files[grep("fasta$",files)]

seqs <- read.fasta(paste0(folder,files[1]))

seqlens <- unlist(lapply(seqs,length))

# times <- numeric()
# for (i in sample(1:length(seqs),100)) {
#   start <- Sys.time()
#   pwtmp <- pairwiseAlignment(tolower(paste(seqs[[i]],collapse="")),ref)
#   end <- Sys.time()
#   times <- c(times,end-start)
# }
# summary(times)

folds <- list()
seqidx <- 1:length(seqs)
for (i in 1:cores) {
  foldsize <- min(ceiling(length(seqs)/cores),length(seqidx))
  folds[[i]] <- sample(seqidx,foldsize)
  seqidx <- seqidx[!seqidx %in% folds[[i]]]
}

fname <- paste0(folder,"pwalign.rds")
if (file.exists(fname)) {
  seqal <- readRDS(fname)
} else {
  seqal <- parallel::clusterApplyLB(cl,1:length(folds),function(i,pairwiseAlignment,seqs,ref,folds) {
    seqal <- matrix("-",length(folds[[i]]),nchar(ref))
    count <- 0
    for (j in folds[[i]]) {
      count <- count + 1
      pwtmp <- pairwiseAlignment(tolower(paste(seqs[[j]],collapse="")),ref,
                                 type="local")
      seqal[count,] <- as.matrix(pwtmp)
    }
    return(seqal)
  },pairwiseAlignment,seqs,ref,folds)
  saveRDS(seqal,fname)
}
print("Finished pairwise alignment.")

# parallel::stopCluster(cl)

aligned <- do.call(rbind,seqal)

rm(seqal)
gc()

coverage <- apply(aligned,2,function(x) {
  y <- sum(x!="-")
  return(y)
})
summary(coverage)
gc()

pdf(paste0(folder,"coverage.pdf"),width=40,height=10)
plot(1:length(coverage),rep(max(log10(coverage+1)),length(coverage)),
     xlab="Genome position",yaxt="n",ylab="Coverage (reads per position)",
     col=rgb(0,0,0,0),ylim=c(0,max(log10(coverage+1))))
polygon(c(1:length(coverage),length(coverage):1),c(log10(coverage+1),rep(0,length(coverage))),
        col=rgb(0.5,0.5,0,0.25),border=rgb(0.5,0.5,0,0.75))
axis(2,0:10,c(0,10^(1:10)))
dev.off()
print("Finished coverage plot")

# find peak?

df <- diff(coverage)
which.max(df)
right <- which.min(df)
left <- rev(order(df))[2]

targets <- which(aligned[,left+1]!="-" & aligned[,right]!="-")

# do frequency

alfreq <- lapply(1:ncol(aligned),function(i) {
  x <- aligned[,i]
  y <- table(x)
  y <- y[names(y)!="-"]
  names(y) <- tolower(names(y))
  tmp <- numeric(4)
  names(tmp) <- c("a","c","g","t")
  tmp[names(y)] <- y
  if (length(y)==0) {
    y <- numeric(1)
    names(y) <- "N"
  }
  tmp <- c(tmp,tmp/sum(tmp),as.numeric(tmp==max(tmp)),names(y)[which.max(y)],
           paste(as.numeric(tmp==max(tmp)),collapse=""))
  return(tmp)
})

alfreq <- data.frame(do.call(rbind,alfreq))
for (i in 1:12) alfreq[,i] <- as.numeric(alfreq[,i])
colnames(alfreq) <- c(rep(c("a","c","g","t"),3),"Fasta","Verketten")
openxlsx::write.xlsx(alfreq,file=paste0(folder,"nucleotide_frequency.xlsx"))

iupac <- function(x) {
  y <- switch(x,
              a="a",c="c",g="g",t="t",
              ag="r",ct="y",gt="k",ac="m",cg="s",at="w",
              cgt="b",agt="d",act="h",acg="v",
              acgt="n")
  return(y)
}

cons <- list()
for (i in c(1,2,5,10,15,20,30)) {
  cons[[paste0(i,"%")]] <- lapply(1:nrow(alfreq),function(j) {
    x <- alfreq[j,]
    y <- colnames(alfreq)[5:8][which(x[5:8]>(i/100))]
    y <- iupac(paste(sort(y),collapse=""))
    if (is.null(y)) y <- "x"
    return(y)
  })
  cons[[paste0(i,"%")]][cons[[paste0(i,"%")]]==""] <- "x"
  cons[[paste0(i,"%")]] <- toupper(cons[[paste0(i,"%")]])
  gc()
}

write.fasta(cons,file=paste0(folder,"consensus.fasta"),names=names(cons),nbchar=80)
print("Finished consensus sequences (full genome).")

q(save="no",status=0)

# some copy stuff

system("mkdir zipped")
folders <- list.files(".")
for (folder in folders) {
  if (folder=="zipped") next()
  files <- list.files(paste0(folder,"/."))
  system(paste0("mkdir zipped/",folder))
  for (file in files) {
    if (file %in% c("ngs.txt","ngs.fasta") | length(grep("bam",file))>0) next()
    system(paste0("scp ",folder,"/",file," zipped/",folder,"/",file))
  }
}


# einmal based on codons einmal based on aa
genes <- read.csv(regions)
AAs <- c("A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V","*")

for (i in 1:nrow(genes)) {
  gene <- genes[i,1]:genes[i,2]
  genefreq <- as.matrix(alfreq[gene,5:8])
  aaposmat <- diag(floor(length(gene)/3))
  aaposmat <- aaposmat[,rep(1:ncol(aaposmat),each=3)]
  for (j in c(1,2,5,10,15,20,30)) {
    idx <- which(genefreq>(j/100))
    genefreq <- genefreq*0
    genefreq[idx] <- 1
    tmp <- aaposmat%*%genefreq
    codons <- apply(tmp,1,function(x) {
      colnames(tmp)[]
    })
    for (k in 1:floor(length(gene)/3)) {

    }
  }
