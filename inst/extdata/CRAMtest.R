library(Rhtslib)
library(Rsamtools) 
library(rtracklayer) 
library(GenomicAlignments)

regions <- c(ABL1="9:130713016-130887675:+", # ABL1
             NUP214="9:131125586-131234663:+", # NUP214
             XKR3="22:16783412-16821699:-", # XKR3
             BCR="22:23179704-23318037:+", # BCR
             MT="MT:1-16569:*")
sieve <- sort(as(regions, "GRanges"))
genome(sieve) <- "GRCh38"

# SCRAP == RNA + ATAC + MT reads from the same cell (in house protocol)
# now that I think about it, we could match those up for the testing 

# check equivalence of reads loaded from BED, BAM, and CRAM
refgenome <- "GRCh38.fa"
bedfile <- "SCRAP207.chr9.chr22.chrM.GRCh38.bed.gz"
bedgr <- import(bedfile, which=sieve)
length(bedgr) # additional entries due to splitting 

bamfile <- "SCRAP207.chr9.chr22.chrM.GRCh38.bam"
bamgr <- granges(readGAlignments(bamfile, param=ScanBamParam(which=sieve)))
length(bamgr) 

cramfile <- "SCRAP207.chr9.chr22.chrM.GRCh38.cram"
cramgr <- granges(readCRAM(cramfile, param=ScanCramParam(which=sieve)))
length(cramgr)
