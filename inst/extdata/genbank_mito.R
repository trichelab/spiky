data("genbank_mito", package="spiky") 
show(genbank_mito) 

library(genbankr)
get_gbas <- function(y) lapply(lapply(y, GBAccession), readGenBank)
GBAs <- get_gbas(genbank_mito$accession)
names(GBAs) <- vapply(GBAs, function(gba) seqnames(seqinfo(gba)), character(1))

get_seqs <- function(x) DNAStringSet(lapply(lapply(x, getSeq), `[[`, 1))
genbank_mito$sequence <- get_seqs(GBAs)
mtSpikes <- process_spikes(genbank_mito) 

mtTxDbs <- lapply(GBAs, makeTxDbFromGenBank)
mtGenes <- lapply(mtTxDbs, genes)

get_mt_tx_seqs <- function(gba) {
  gxs <- genes(gba)
  names(gxs) <- gxs$gene
  mtseq <- getSeq(gba)
  mtgxs <- granges(gxs)
  mtgxs$gene <- gxs$gene
  mtgxs$sequence <- getSeq(mtseq, gxs)
  mtgxs
}

mtGeneSeqs <- lapply(GBAs, get_mt_tx_seqs)
mtTxSpikes <- lapply(mtGeneSeqs, process_spikes)
