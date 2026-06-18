# activate hmsc conda environment
# conda activate hmsc
library("phyloseq")

#Load tables
setwd("/HPC/bioinformatics/16s/both_silva/results")

###Make phyloseq object from raw data
otu_mat<-as.matrix(read.table("no_sing_cleaned_otu_table_ASV_wise.txt", sep="\t", header=T, row.names=1,check.names=F))

tax_mat_silva<-data.frame(read.table("tax_silva_dada2.txt", sep='\t', header=T, comment=""))
tax_mat_silva_d<-as.data.frame(unclass(tax_mat_silva), stringsAsFactors=T)

cat("\n")
cat("raw_silva_tax_summary")
nrow(tax_mat_silva_d)
cat("\n")
levels(factor(tax_mat_silva_d[,1]))
cat("\n")
head(tax_mat_silva_d)
cat("\n")
summary(tax_mat_silva_d)
cat("\n")

#Select sequences according to cleaned ASV table
n_tax<- tax_mat_silva[rownames(tax_mat_silva) %in% rownames(otu_mat),]
n_tax_d<-as.data.frame(unclass(n_tax), stringsAsFactors=T)

cat("\n")
cat("silva_match_ASV_table_tax_summary")
nrow(n_tax_d)
cat("\n")
levels(factor(n_tax_d[,1]))
cat("\n")
head(n_tax_d)
cat("\n")
summary(n_tax_d)
cat("\n")

#1st round
d_tax<-data.frame(n_tax)
nrow(d_tax) # 67361

d_tax2<-d_tax[!is.na(d_tax$Phylum), ]
nrow(d_tax2) # 67361
unique(d_tax2$Phylum)

d_tax3<-d_tax2[!is.na(d_tax2$Class), ]
sort(unique(d_tax3$Phylum))
sort(unique(d_tax3$Class))
nrow(d_tax3) # 67361

d_tax4<-subset(d_tax3, !(Class=="uncultured"))
nrow(d_tax4) # 67361
sort(unique(d_tax4$Order))

d_tax5<-subset(d_tax4, !(Order=="Chloroplast"))
nrow(d_tax5) # 54435
sort(unique(d_tax5$Order))

#Filter otu_table according to taxonomy
n_otu_mat<-otu_mat[rownames(otu_mat) %in% rownames(d_tax5),]
#####

#Load metadata
setwd("/HPC/bioinformatics/16s/both_silva/results/metadata")
metadata<-read.table("no_control_no_sing_samples_cleaned_metadata_ASV_wise.txt", sep="\t", header=T)
sampledata = sample_data(data.frame(metadata, row.names=metadata$sample_ID, stringsAsFactors=FALSE))

OTU = otu_table(n_otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(d_tax5))

p_DADAwang = phyloseq(OTU, TAX)

p_DADAwang
DADAwang1 = merge_phyloseq(p_DADAwang, sampledata)

DADAwang1

cat("\n")
cat("before_empty_removal_totals")
DADAwang1
cat("\n")
sum(sample_sums(DADAwang1))
cat("\n")

tudao = filter_taxa(DADAwang1, function(x) sum(x) > 0, TRUE)
tudao = prune_samples(sample_sums(tudao)>0,tudao)

cat("\n")
cat("after_empty_removal_totals")
tudao
cat("\n")
sum(sample_sums(tudao))
cat("\n")

tax_mat_kon<-data.frame(tax_table(tudao))
new_otu_mat<-data.frame(otu_table(tudao),check.names=F)

tax_mat_b<-as.matrix(tax_mat_kon[,c("Kingdom","Phylum","Class","Order","Family","Genus")])

setwd("/HPC/bioinformatics/16s/both_silva/results/metadata")
write.table(data.frame(sample_data(tudao), check.names=F), "cleaned_silva_metadata_3.txt", sep="\t", quote=FALSE, row.names=TRUE)

setwd("/HPC/bioinformatics/16s/both_silva/results")
write.table(new_otu_mat, "cleaned_otu_silva_3.txt", sep="\t", quote=FALSE, row.names=TRUE)
write.table(tax_mat_b, "cleaned_tax_silva_3.txt", sep="\t", quote=FALSE, row.names=TRUE)