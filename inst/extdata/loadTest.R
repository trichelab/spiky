# column 1 (UMI) is contigUMI5_UMI5 where each UMI5 is from one read in the pair
# seq is fragment sequence after adapter trimming and UMI decoding 
# rname is the contig fragment of the UMI variable
# pos is position along the unplaced/MT contig in hg38
# 6547 is one sample
# 6548 is another 
#
testfile <- "2020_Human0.01_test.csv.gz"
test <- read.table(testfile, sep=" ", head=T) # bit mysterious, see above
# really only need UMI, pos, seq, read_count for each subject/UMI combo
testGR <- parse_spike_UMI(test$UMI, pos=test$pos, seqs=test$seq) 
save(testGR, file="../../data/testGR.rda") 
