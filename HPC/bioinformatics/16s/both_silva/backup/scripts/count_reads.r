# Counting reads at each step of 18S pipeline
# Run in R from HPC/bioinformatics/16s/both_silva
# Start interactive job
# srun --mem=128g -c 1 --account "edna" --pty /bin/bash
# Activate conda environment with R
# conda activate hmsc

# Count reads before chimera removal
raw<-read.table("../base/results/DADA2_raw.table", sep="\t",row.names = 1,header=TRUE)
sum(colSums(raw))   
# 998756902 = 998 M reads

# Count reads after chimera removal
nochim<-read.table("../base/results/DADA2_nochim.table", sep="\t",row.names = 1,header=TRUE)
sum(colSums(nochim))
# 966617118 = 967 M reads
nrow(nochim)
# 1532113 ASVs

# Count reads after cleaning based on controls
clean<-read.table("results/cleaned_otu_table_ASV_wise_EES.txt", sep="\t", header=T, row.names=1,check.names=F)
sum(colSums(clean))
# 471966847 = 472 M reads
nrow(clean)
# 582,538 ASVs

# Count reads after removal per sample of sequences found in a single PCR replicate
no.sing<-read.table("results/no_sing_cleaned_otu_table_ASV_wise_EES.txt", sep="\t", header=T, row.names=1,check.names=F)
sum(colSums(no.sing))
# 422705690 = 423 M reads
nrow(no.sing)
# 67,361 ASVs

# Count reads after taxonomic filtering
tax<-read.table("results/cleaned_otu_silva_EES.txt", sep="\t", header=T, row.names=1,check.names=F)
sum(colSums(tax))
# 422705690 = 423 M reads
nrow(tax)
# 67,361 ASVs
