# column 1 (UMI) is contigUMI5_UMI5 where each UMI5 is from one read in the pair
# seq is fragment sequence after adapter trimming and UMI decoding 
# rname is the contig fragment of the UMI variable
# pos is position along the unplaced/MT contig in hg38
# 6547 is one sample
# 6548 is another 
#
testfile <- system.file("extdata","2020_Human0.01_test.csv.gz",package="spiky")
test <- read.table(testfile, sep=" ", head=T) # bit mysterious, see above
# really only need UMI, pos, seq, read_count for each subject/UMI combo
testGR <- parse_spike_UMI(test$UMI, pos=test$pos, seqs=test$seq) 
# how can we make it so that spiky:::.addFragInfo() works on this? 
# at present, .addFragInfo needs a `frag_grp` column, and test.csv lacks one
message("Ideally, we would like to have frag_grp information in testGR")
save(testGR, file="../../data/testGR.rda") 
