#Dada2 pipeline Feb-2024 for 18S data
#https://benjjneb.github.io/dada2/tutorial.html

#Isabelle van der Zanden (based on Luis Merlotti's script)

#Installing cutadapt
#https://cutadapt.readthedocs.io/en/stable/installation.html
#conda create -n cutadaptenv cutadapt

#Cutadapt local in the local terminal (macbook)
#conda activate cutadapt
#which cutadapt

###### Dada2 pipeline for protists

#Activating packages

library(Rcpp)
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
library(dplyr) #data manipulation
library(phyloseq)
library(ggvenn)
library(tidyverse)

#Setup the work folder
setwd("/home/nioo/isabellez/sequencing/18S/18S_4R")

path <- "/home/nioo/isabellez/sequencing/18S/18S_4R/rawdata"  ## CHANGE ME to the directory containing the fastq files.
list.files(path)

#Forward and Reverse in objects

fnFs <- sort(list.files(path, pattern = "_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq.gz", full.names = TRUE))

#Indentify Primers
#In this data, primers and barcodes from sequencing were removed already

FWD <- "GGCAAGTCTGGTGCCAG"  ## 3NDf
REV <- "TCCGTCAATTYCTTTAAGT"  ## 1132rmod

# In theory if you understand your amplicon sequencing setup, 
# this is sufficient to continue. However, to ensure we have the right primers, 
# and the correct orientation of the primers on the reads, 
# we will verify the presence and orientation of these primers in the data.

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# PREFILTERING
# The presence of ambiguous bases (Ns) in the sequencing reads makes accurate 
# mapping of short primer sequences difficult. 
# Next we are going to “pre-filter” the sequences just to remove those with Ns, 
# but perform no other filtering.

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = 40, matchIDs = TRUE)

# We are now ready to count the number of times the primers appear in the 
# forward and reverse read, while considering all possible primer orientations

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[2]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[2]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[2]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[2]]))


# Forward Complement Reverse RevComp
# FWD.ForwardReads  118470          0       0       0
# FWD.ReverseReads       0          0       0   44856
# REV.ForwardReads       0          0       0   37936
# REV.ReverseReads  111103          0       0       0

#Sequences have primers still

#Install cutadapt on the Bash language:
#http://cutadapt.readthedocs.io/en/stable/index.html
#Type:

# virtualenv my_virtual_env 
# source my_virtual_env/bin/activate
# pip install cutadapt

#Local of installation in my macbook
#/opt/anaconda3/bin/cutadapt

#Remove Primers
#Warning: A lot of output will be written to the screen by cutadapt!

#chmod ugo+rwx /home/nioo/luism/.conda/envs/cutadaptenv

cutadapt <- "/home/nioo/isabellez/.local/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R
#4.6

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt

#I changed the "-n" to 6. Then, the command will ru 6 times at the same sequencing
#I have primers in the middle of the sequences

for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# As a sanity check, we will count the 
# presence of primers in the first cutadapt-ed sample:

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Forward Complement Reverse RevComp
# FWD.ForwardReads       0          0       0       0
# FWD.ReverseReads       0          0       0       0
# REV.ForwardReads       0          0       0       0
# REV.ReverseReads       0          0       0       0

#Success! Primers are no longer detected in the cutadapted reads.

#Following with the dada2 analysis
# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 8)
# sample name extraction from adapted pipeline Luis
# get.sample.name <- function(fname) strsplit(basename(fname), "_16")[[1]][1]
# sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

#Inspect read quality profiles

# plotQualityProfile(cutFs[1:3])
# plotQualityProfile(cutRs[1:3]) #this one doesnt work

#The comand will not work properly with the NOVOSeq (plattaform) data
#The Phreed score is different
#Check:
#https://github.com/benjjneb/dada2/issues/791

#Filter and trim
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, maxN = 0, maxEE = 2, 
                     truncQ = 2, rm.phix = TRUE, compress = TRUE, matchIDs = TRUE, truncLen=(230), multithread=40)  # on windows, set multithread = FALSE
head(out)

#Solving the problem with Learn Error argument (due the NovoSeq plattaform)
#Check the following post issue
#https://github.com/ErnakovichLab/dada2_ernakovichlab#learn-the-error-rates
#https://github.com/benjjneb/dada2/issues/1307

library(magrittr)
library(dplyr)

# Error function1 alter loess arguments (weights and span and enforce monotonicity)

loessErrfun_mod1 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

# check what this looks like
errF_1 <- learnErrors(
  filtFs,
  multithread = 40,
  nbases = 1e8,
  errorEstimationFunction = loessErrfun_mod1,
  verbose = TRUE,
  randomize = TRUE,
  MAX_CONSIST = 20
)

#Visualize the estimated error rates as a sanity check.

plotErrors(errF_1, nominalQ = TRUE)
# plotErrors(errR_1, nominalQ = TRUE)

#Dereplicate identical reads
#~/Desktop/fungi//cutadapt/filtered/
derepFs <- derepFastq(filtFs, verbose = TRUE)
# derepRs <- derepFastq(filtRs, verbose = TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
# names(derepRs) <- sample.names

# Sample Inference
# At this step, the core sample inference algorithm is applied to the dereplicated data.
#I added the errF to the Forward and Reverse sequences
dadaFs <- dada(derepFs, err = errF_1, multithread = 40)
# dadaRs <- dada(derepRs, err = errR_1, multithread = 40)

# Merge paired reads
# mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#Construct Sequence Table
#We can now construct an amplicon sequence variant table (ASV) table, 
#a higher-resolution version of the OTU table produced by traditional methods.
# seqtab <- makeSequenceTable(mergers)
seqtab <- makeSequenceTable(dadaFs) #if only using forward
dim(seqtab)

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=40, verbose=TRUE)

#Inspect distribution of sequence lengths:
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) #closer to 1, better! 0.9733185
table(nchar(getSequences(seqtab.nochim)))

# Track reads through the pipeline
# We now inspect the the number of reads that made 
#it through each step in the pipeline to verify everything worked as expected.

# getN <- function(x) sum(getUniques(x))
# track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, 
#                                                                        getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
# colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
#                      "nonchim")
# rownames(track) <- sample.names
# head(track)

# if only using forward:
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)

#Saving the out info
write.table(track, "track.csv", sep="\t", quote=F)
saveRDS(seqtab.nochim, "seqtab.nochim.rds")

#Assign taxonomy
pr2_ref <- "/home/nioo/isabellez/sequencing/18S/PR2_database/pr2_version_5.0.0_SSU_dada2.fasta.gz"
taxa <- assignTaxonomy(seqtab.nochim, pr2_ref, multithread = 40, tryRC = TRUE)

# Removing sequence rownames for display only
taxa.print <- taxa 
rownames(taxa.print) <- NULL
head(taxa.print)

#Saving the taxa file
write.table(taxa, "taxa.csv", sep="\t", quote=F)

#Adjusting the names
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

# count table
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.csv", sep="\t", quote=F, col.names=NA)

# tax table
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy.csv", sep="\t", quote=F, col.names=NA)

# Concatenate Counts and Taxonomy tables
counts <- read.table("ASVs_counts.csv", header = TRUE, sep = "\t", dec = ".", row.names = 1)
taxonomy <- read.table("ASVs_taxonomy.csv", header = TRUE, sep = "\t", dec = ".", row.names = 1)
counts_taxonomy <- as.data.frame(c(taxonomy, counts))
write.table(counts_taxonomy, "ASVs_tax_counts.csv", sep="\t", quote=F, col.names=NA)

# (the above command did the same) Final ASV table:
# 
# #The c(1:16) is regarding the object counts
# 
# ASV_table <- summarise_at(group_by(counts_taxonomy, Kingdom, Phylum, Class, Order, Family, Genus), c(1:16), tibble::lst(sum))
# 
# write.csv(ASV_table, "ASV_table.csv") #The ASV table will be saved with all taxonomy

samples.out <- rownames(seqtab.nochim)
vb_code <- read.csv2("dna_code_18S.csv")
samdf <- data.frame(Sample=samples.out)
samdf <- samdf %>% separate(Sample, c("i5","dna_code"))
samdf$dna_code <- as.integer(samdf$dna_code)
samdf <- left_join(samdf,vb_code,by="dna_code")

rownames(samdf) <- samples.out
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(samdf),
               tax_table(taxa))
saveRDS(ps, "ps_18S.rds")
rank_names(ps)
colnames(tax_table(ps)) <- c("Domain", "Clade_1", "Clade_2", "Kingdom", "Phylum", "Subphylum", "Class", "Genus", "Species")