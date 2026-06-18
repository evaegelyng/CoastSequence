# Before running this script, check seq lengths_
# conda activate bioinfo
# seqmagick info base/results/COSQ_16S.fa 
# name                alignment    min_len   max_len   avg_len  num_seqs
# results/COSQ_16S.fa FALSE            160       448    254.82   1532113
# As the sequences are above 250 nt long, we use a bootstrap threshold of 80 below (see https://benjjneb.github.io/dada2/assign.html)

# Then move to the both_silva folder and run:
# conda activate metabar_2021
library(dada2); packageVersion("dada2")

set.seed(100) # Initialize random number generator for reproducibility

# Load sequences, OTU table and reference database
seqs <- getSequences("/HPC/bioinformatics/16s/base/results/DADA2_nochim.otus")
seqtab.nochim <- read.table("/HPC/bioinformatics/16s/base/results/DADA2_nochim.table", sep="\t", header=T, row.names=1,check.names=F)
ref_fasta <- "/HPC/bioinformatics/16s/both_silva/data/silva_nr99_v138.1_train_set.fa.gz"

# Add sequences to OTU table
seqs2 <- as.data.frame(seqs)
seqtab2 <- seqtab.nochim
row.names(seqtab2) <- seqs2$seqs
seqtab <- as.matrix(t(seqtab2))

# assign taxonomy in batches of 50,000
to_split <- seq(1, ncol(seqtab), by = 50000)
to_split2 <- c(to_split[2:length(to_split)]-1, ncol(seqtab))

taxtab = NULL

for(i in 1:length(to_split)){
  seqtab2 <- seqtab[, to_split[i]:to_split2[i]]
  taxtab2 <- assignTaxonomy(seqtab2, refFasta = ref_fasta, multithread = TRUE, minBoot=80)
  #if(!is.null(ref_fasta_spp)){taxtab2 <- addSpecies(taxtab2, refFasta = ref_fasta_spp, verbose = TRUE)}
  if(!is.null(taxtab)){taxtab <- rbind(taxtab, taxtab2)}
  if(is.null(taxtab)){taxtab <- taxtab2}
}
# save files in case phylogeny does not run
saveRDS(taxtab, "taxtab.rds")

#seqtab.nochim.withtaxa1 = merge(taxa,t(seqtab.nochim),by=0)
#write.table(seqtab.nochim.withtaxa1,"results/DADA2_classified.tsv",sep='\t')