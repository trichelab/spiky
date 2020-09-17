# load
dedup <- read.table("2020_dedup_spikein_data.csv.gz", sep=" ")[, c(1, 5, 6)]
# frag_grp has all information necessary to reconstruct the rest, so
save(dedup, file="../../data/dedup.rda") 

# melt
library(reshape2) 
spike_read_counts <- melt(dedup, id.vars="frag_grp") 
colnames(spike_read_counts) <- sub("^variable$", "id", 
                                   colnames(spike_read_counts))
colnames(spike_read_counts) <- sub("^value$", "read_count", 
                                   colnames(spike_read_counts))
spike_read_counts$id <- as.factor(sub("read_count_", "run", 
                                      spike_read_counts$id))
save(spike_read_counts, file="../../data/spike_read_counts.rda") 

                              
