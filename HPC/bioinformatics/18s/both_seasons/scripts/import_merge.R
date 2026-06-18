#library(igraph)
library("phyloseq")
library("ggplot2")
library("vegan")
library("reshape2")
library("plyr")
#library("scales")
#library("stringr")
#library("RColorBrewer")
#library("corrplot")
#library("gllvm")
#library("gclus")
#library(ALDEx2)
#library(propr)
#library(zCompositions)
#library("ade4")
#library("psych")
#library("ggrepel")
#library("ggVennDiagram")

#Load metadata
setwd("/HPC/bioinformatics/18s/both_seasons/results/metadata")
p_metadata<-read.table("metadata_both.txt", sep="\t", header=T)
head(p_metadata)

#Load tables
#sed
setwd("/environmental_data/data")
CN_sed_both<-read.table("CNdelta_iso_sed_both_new.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))
CNT_sed_both<-read.table("CNT_sed_both_new.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))
TP_sed_autunm<-read.table("TP_org_inorg_dens_sed_autumn.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))
TP_sed_spring<-read.table("TP_org_inorg_dens_sed_spring.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))

#wat
NUT_wat_both<-read.table("NUT_ctd_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))
dOdH_wat_both<-read.table("dOdH_wat_both.txt", sep="\t", header=T, check.names=F, na.strings=c(""," ","NA"))

#Adding salinity missing values that can be estimated from dOdH isotope

merged_wat<-merge(NUT_wat_both, dOdH_wat_both, by = "Sample_ID", all = TRUE)
merged_wat$po<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), tail, 1)
merged_wat$cl<-as.integer(gsub('\\D','', merged_wat$po))
merged_wat$pn<-gsub('\\d','_', merged_wat$po)
merged_wat$pn2<-gsub(".*__(.+).*", "\\1", merged_wat$pn)
merged_wat$pn3<-gsub(".*_(.+).*", "\\1", merged_wat$pn2)
merged_wat$hb<-ifelse(merged_wat$pn3=="EW", "eelgrass", ifelse(merged_wat$pn3=="RW", "rocks", "sand"))
merged_wat$poi<- sapply(strsplit(as.character(merged_wat$Sample_ID), "2C"), head, 1)
merged_wat$sn<-ifelse(merged_wat$poi=="", "autumn", "spring")
merged_wat$snch<-paste(merged_wat$sn,merged_wat$cl,merged_wat$hb,sep="_")
merged_wat<-subset(merged_wat,!cl==2) 

cdata_sal <- ddply(merged_wat, c("sn","cl"), dplyr::summarize, mean_sal = mean(Salinity), sd_sal   = sd(Salinity), mean_d12H = mean(d12H), sd_d12H   = sd(d12H))


#####BBBBB
########
#Testing importing salinity values from d12h
########

y <- merged_wat$Salinity
x <- merged_wat$d12H
pred_sal<-predict(lm(y ~ x))
merged_wat_pred<-merged_wat
merged_wat_pred$Salinity[84:87]<-pred_sal[84:87]

cdata_sal_pred <- ddply(merged_wat_pred, c("sn","cl"), dplyr::mutate, mean_sal = mean(Salinity), sd_sal   = sd(Salinity), mean_d12H = mean(d12H), sd_d12H   = sd(d12H))

##Keep going

head(CN_sed_both)
head(CNT_sed_both)
head(TP_sed_autunm)
head(TP_sed_spring)

nrow(CN_sed_both)
nrow(CNT_sed_both)
nrow(TP_sed_autunm)
nrow(TP_sed_spring)
nrow(TP_sed_spring)+nrow(TP_sed_autunm)

TP_sed_both<-rbind(TP_sed_autunm, TP_sed_spring)
merged_sed<-Reduce(function(x,y) merge(x = x, y = y, by = "Sample_ID", all = TRUE), list(CN_sed_both, CNT_sed_both, TP_sed_both))

head(merged_sed)
head(merged_wat)
nrow(merged_wat)
nrow(merged_sed)

#Now incorporate the salinity groups back to metadata

######Here using a different way of grouping for eelgrass
#merged_wat_pred$sal_group<-ifelse(merged_wat_pred$hb=="eelgrass", both_cs_comp_prede$sgroup[match(merged_wat_pred$snch, both_cs_comp_prede$snch)], ifelse(merged_wat_pred$Salinity<15,1,ifelse(merged_wat_pred$Salinity<25,2,3)))

######Here using the same way of grouping for all habitats -> might generate small sample size for eelgrass
merged_wat_pred$sal_group<-ifelse(merged_wat_pred$Salinity<15,1,ifelse(merged_wat_pred$Salinity<25,2,3))

merged_sed$NP_ratio<-merged_sed$N_per_kg_dry_sample/merged_sed$TP

#Clean tables (preserve only meaningful variables)
merged_sed2<-merged_sed[,c("Sample_ID","Norm_d14N_15N","Norm_d12C_13C",
"Remarks.x","QC.x","MISS.x","N_per_kg_dry_sample","C_per_kg_dry_sample",
"CN_ratio","Remarks.y","QC.y","MISS.y","Grain_size","Density_wet",
"Dry_matter_per_WW","Watercontent","Organic_content","Inorganic_content","TP","NP_ratio")]

#jump
###########
#Check if is possible to fill variables using linear regressions (Answer is NO)
setwd("/environmental_data/results")

ggplot(merged_sed2, aes(log(TP), log(C_per_kg_dry_sample))) + 
geom_point() + 
theme_bw() + 
theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 90, hjust = 1, size=4.5, vjust=0.5), strip.text = element_text(size=6), legend.title=element_text(size=8), legend.text=element_text(size=6), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"))
ggsave("log_TPxTC.pdf")

ggplot(merged_sed2, aes(log(TP), log(N_per_kg_dry_sample))) + 
geom_point() +  
theme_bw() + 
theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 90, hjust = 1, size=4.5, vjust=0.5), strip.text = element_text(size=6), legend.title=element_text(size=8), legend.text=element_text(size=6), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"))
ggsave("log_TPxTN.pdf")

ggplot(merged_sed2, aes(log(N_per_kg_dry_sample), log(C_per_kg_dry_sample))) + 
geom_point() + 
theme_bw() + 
theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 90, hjust = 1, size=4.5, vjust=0.5), strip.text = element_text(size=6), legend.title=element_text(size=8), legend.text=element_text(size=6), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"))
ggsave("log_TNxTC.pdf")


merged_sed3<-merged_sed[,c("Sample_ID","Norm_d14N_15N","Norm_d12C_13C",
"N_per_kg_dry_sample","C_per_kg_dry_sample",
"CN_ratio","Density_wet",
"Dry_matter_per_WW","Grain_size","Watercontent","Organic_content","Inorganic_content","TP","NP_ratio")]

###FIlter ultra outliers
#    CN_ratio  84 [318]
#    Norm_d12C_13C   106 [234]
#    TP      5.3 [383]
#    N 179,180
#Must remove another d12C_13C outlier (something lower than -105..)

merged_sed3$TP[383]<-NA
merged_sed3$Norm_d12C_13C[234]<-NA
merged_sed3$CN_ratio[318]<-NA
merged_sed3$N_per_kg_dry_sample[179:180]<-NA
merged_sed3$NP_ratio[179:180]<-NA


colnames(merged_sed3)
colnames(merged_sed3)<-c("Sample_ID","d14N_15N","d12C_13C","N","C","CN_ratio",           
"Density","Dry_matter","Grain_size","Water_content","Organic_content","Inorganic_content","TP","NP_ratio")


#   It might be OK to import C values from Dry_matter or water content, however it wouldn't solve the other variables with NAs [N, CN]

#   Therefore, whenever using Sed data, be careful with C, N - some NAs

########
#back to water

merged_wat2<-merged_wat_pred[,c("Sample_ID","NO2NO3","NO2","NH3","PO4",
"Si","NO3","Latitude","Longitude","Chlorophyll","Salinity",
"Temperature","CTD_depth","d18O","SD_d18O","d12H","SD_d12H","sal_group","snch","Time")]

#Merge with metadata
p_md_wat<- subset(na.omit(p_metadata), substrate_type=="water")
md_wat<-unique(p_md_wat[c("season", "cluster", "habitat")])
md_wat$Sample_ID<-ifelse(md_wat$season=="spring",paste("C",md_wat$cluster,md_wat$habitat,sep=""),
paste("2C",md_wat$cluster,md_wat$habitat,sep=""))
nrow(md_wat)
nrow(merged_wat2)

fwat<-merge(md_wat, merged_wat2, by = "Sample_ID", all = TRUE)
fwat<-subset(fwat,!cluster==2) 


p_md_sed<- subset(na.omit(p_metadata), substrate_type=="sediment")
md_sed<-unique(p_md_sed[c("sample_root", "season", "cluster", "habitat", "field_replicate")])
colnames(md_sed)<-c("Sample_ID", "season", "cluster", "habitat", "field_replicate")
nrow(md_sed)
nrow(merged_sed3)
fsed<-merge(md_sed, merged_sed3, by = "Sample_ID", all = TRUE)

#Remove cluster 2
fsed<-subset(fsed,!cluster==2) 

head(fwat)
head(fsed)
nrow(fwat)
nrow(fsed)

######
##############
#################
###If needed to retrieve data for a subset of samples (example):

aalborg<-c("C1RB1","C1SB1","C2RB1","C3RB1","C3SB1","C4RB1","C4SB1","C5RB1",
"C5SB1","C6RB1","C6SB1","C7RB1","C7SB1","C8EB1","C8RB1","C8SB1",
"C9EB1","C9RB1","C9SB1","C10EB1","C10RB1","C10SB1","C11RB1","C11SB1",
"C12RB1","C12SB1","C13RB1","C13SB1","C14RB1","C14SB1","C15RB1",
"C15SB1","C16EB1","C16RB1","C16SB1","C17EB1","C17RB1","C17SB1","C18RB1",
"C18SB1","C19EB1","C19RB1","C19SB1","C20EB1","C20RB1","C20SB1",
"C21EB1","C21RB1","C21SB1","C23RB1","C23SB1","C24EB1","C24RB1","C24SB1",
"C25EB1","C25RB1","C25SB1","C26EB1","C26RB1","C26SB1","C27EB1",
"C27RB1","C27SB1","C28EB1","C28RB1","C28SB1","C29EB1","C29RB1","C29SB1",
"C30EB1","C30RB1","C30SB1","C31EB1","C31RB1","C31SB1","C32RB1",
"C32SB1","C33RB1","C33SB1","2C1RB1","2C1SB1","2C3RB1","2C3SB1","2C4RB1",
"2C4SB1","2C5RB1","2C5SB1","2C6RB1","2C6SB1","2C7RB1",
"2C7SB1","2C8EB1","2C8RB1","2C8SB1","2C9EB1","2C9RB1","2C9SB1","2C10EB1",
"2C10RB1","2C10SB1","2C11RB1","2C11SB1","2C12RB1","2C12SB1",
"2C13RB1","2C13SB1","2C14RB1","2C14SB1","2C15RB1","2C15SB1",
"2C16EB1","2C16RB1","2C16SB1","2C17EB1","2C17RB1","2C17SB1","2C18RB1",
"2C18SB1","2C19EB1","2C19RB1","2C19SB1","2C20EB1",
"2C20RB1","2C20SB1","2C21EB1","2C21RB1","2C21SB1","2C23RB1","2C23SB1",
"2C24EB1","2C24RB1","2C24SB1","2C25EB1","2C25RB1",
"2C25SB1","2C26EB1","2C26RB1","2C26SB1","2C27EB1","2C27RB1",
"2C27SB1","2C28EB1","2C28RB1","2C28SB1","2C29EB1","2C29RB1","2C29SB1",
"2C30EB1","2C30RB1","2C30SB1","2C31EB1","2C31RB1",
"2C31SB1","2C32RB1","2C32SB1","2C33RB1","2C33SB1")

aal_md<-merged_sed[merged_sed$Sample_ID %in% aalborg, ]

nrow(aal_md)==length(aalborg)

aal_md2<-aal_md[,c("Sample_ID","Grain_size","Density_wet",
"Dry_matter_per_WW","Watercontent","Organic_content","Inorganic_content","TP")]

p_mda_sed<- subset(na.omit(p_metadata), substrate_type=="sediment")
mda_sed<-unique(p_mda_sed[c("sample_root", "season", "cluster", "habitat", "extraction_refs")])
colnames(mda_sed)<-c("Sample_ID", "season", "cluster", "habitat", "extraction_refs")
nrow(mda_sed)
nrow(aal_md2)
fseda<-merge(mda_sed, aal_md2, by = "Sample_ID", all.y = TRUE)

head(fseda)
nrow(fseda)

fseda$habitat<-as.character(fseda$habitat)
fseda$hab_uni<-ifelse(fseda$habitat=="EB", "Eelgrass", ifelse(fseda$habitat=="RB", "Rocks", ifelse(fseda$habitat=="SB", "Sand", fseda$habitat)))
fseda$shc<-paste(fseda$season, fseda$hab_uni, fseda$cluster, sep="_")

fseda$Latitude<-fwat$Latitude[match(fseda$shc, fwat$shc)]
fseda$Longitude<-fwat$Longitude[match(fseda$shc, fwat$shc)]

fseda2<-fseda[,c("Sample_ID","season","cluster","habitat","extraction_refs",
"Grain_size","Density_wet","Dry_matter_per_WW","Watercontent","Organic_content",
"Inorganic_content","TP","Latitude","Longitude")]

setwd("/environmental_data/results")
write.table(fseda2, "aalborg_samples_metadata.txt", sep="\t", quote=FALSE, row.names=FALSE)

#Make distribution plots (histogram, category-associated boxplots, ~QC, ~MISS, etc)
#Water
setwd("/environmental_data/results/wat")

#Convert negative and 0 nutrient values to 0.01 or so
fwat[,c(5:10)] <- matrix(pmax(unlist(fwat[,c(5:10)]),0), nrow=nrow(fwat))

dl005<-c("NH3", "NO3", "NO2NO3", "Si")
for (g in 1:length(dl005))
{
fwat[,dl005[g]]<-ifelse(fwat[,dl005[g]]<0.05,0.05,fwat[,dl005[g]])
}


dl0005<-c("PO4","NO2")
for (g in 1:length(dl0005))
{
fwat[,dl0005[g]]<-ifelse(fwat[,dl0005[g]]<0.005,0.005,fwat[,dl0005[g]])
}


#log-transform nutrient variables

fwat<-fwat[,c("Sample_ID","season","cluster","habitat","NO2NO3",
"NO2","NH3","PO4","Si","NO3","Chlorophyll",
"Latitude","Longitude","Salinity","Temperature",
"CTD_depth","d18O","SD_d18O","d12H","SD_d12H",
"sal_group","snch","Time")]

for (g in 5:11)
{
fscn<-ncol(fwat)
fwat[,fscn+1]<-log(fwat[,g])
colnames(fwat)[fscn+1]<-paste("log", colnames(fwat)[g], sep="_")
}

p_mui <- melt(fwat, id=c("season", "cluster", "habitat","sal_group"), measure=c("log_NO2NO3", "log_NO2", "log_NH3", "log_PO4", "log_Si", "log_NO3", "log_Chlorophyll", "Salinity", "Temperature", "d18O", "d12H"))

#Remove NAs
mui<-na.omit(p_mui)

#Calculate mean
mui2 <- ddply(mui, .(variable, season), summarise, grp.mean=mean(na.omit(value)))
mui3 <- ddply(mui, .(variable, sal_group), summarise, grp.mean=mean(na.omit(value)))

mui$season <- factor(mui$season, levels = c("spring","autumn"))
mui$habitat <- factor(mui$habitat, levels = c("RW","SW","EW"))
mui2$season <- factor(mui2$season, levels = c("spring","autumn"))

ggplot(mui, aes(x=value, fill=season)) +
geom_histogram(aes(y=..density..), position="identity", alpha=0.5, bins=50, colour="black", size=0.03) + geom_vline(data=mui2, aes(xintercept=grp.mean, color=season), linetype="dashed") + theme_classic() + scale_color_manual(values=c("#00FF33","#FF6600")) + scale_fill_manual(values=c("#00FF33","#FF6600")) + scale_x_continuous(labels = comma) + scale_y_continuous(labels = comma) + facet_wrap(~variable, ncol=3, scales="free") + theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 7), axis.text.y = element_text(size=7), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=7),legend.title=element_text(size=7),legend.text=element_text(size=6)) + theme(legend.key.size = unit(0.3, "cm")) + labs(title="Water env. histogram plot", x ="Value", y = "Count", fill = "season")
ggsave("reads_hist_raw_wat_season_log.pdf")

ggplot(mui, aes(x=value, fill=as.factor(sal_group))) +
geom_histogram(aes(y=..density..), position="identity", alpha=0.5, bins=50, colour="black", size=0.03) + geom_vline(data=mui3, aes(xintercept=grp.mean, color=as.factor(sal_group)), linetype="dashed") + theme_classic() + scale_fill_manual(values=c("#00CCFF","#339966","#999900")) + scale_color_manual(values=c("#00CCFF","#339966","#999900")) + scale_x_continuous(labels = comma) + scale_y_continuous(labels = comma) + facet_wrap(~variable, ncol=3, scales="free") + theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 7), axis.text.y = element_text(size=7), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=7),legend.title=element_text(size=7),legend.text=element_text(size=6)) + theme(legend.key.size = unit(0.3, "cm")) + labs(title="Water env. histogram plot", x ="Value", y = "Count", fill = "sal_group")
ggsave("reads_hist_raw_wat_salg3_log.pdf")

##Habitat-wise
#Calculate mean
hmui2 <- ddply(mui, .(variable, season, habitat), summarise, grp.mean=mean(na.omit(value)))
hmui3 <- ddply(mui, .(variable, sal_group, habitat), summarise, grp.mean=mean(na.omit(value)))

hmui2$season <- factor(hmui2$season, levels = c("spring","autumn"))
hmui2$habitat <- factor(hmui2$habitat, levels = c("RW","SW","EW"))
hmui3$habitat <- factor(hmui3$habitat, levels = c("RW","SW","EW"))

ggplot(mui, aes(x=value, fill=habitat)) +
geom_histogram(aes(y=..density..), position="identity", alpha=0.4, bins=25) + geom_vline(data=hmui2, aes(xintercept=grp.mean, color=habitat), linetype="dashed") + geom_density(aes(colour=habitat), fill=NA) +theme_classic() + scale_x_continuous(labels = comma) + scale_y_continuous(labels = comma) + facet_wrap(variable~season, ncol=4, scales="free") + theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 7), axis.text.y = element_text(size=7), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=7),legend.title=element_text(size=7),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm")) + labs(title="Water env. histogram plot", x ="Value", y = "Count", fill = "habitat") + scale_color_brewer(palette="Set1") +
scale_fill_brewer(palette="Set1")
ggsave("reads_hist_raw_wat_habitat_season_log.pdf")

ggplot(mui, aes(x=value, fill=habitat)) +
geom_histogram(aes(y=..density..), position="identity", alpha=0.4, bins=25) + geom_vline(data=hmui3, aes(xintercept=grp.mean, color=habitat), linetype="dashed") + geom_density(aes(colour=habitat), fill=NA) +theme_classic() + scale_x_continuous(labels = comma) + scale_y_continuous(labels = comma) +
scale_fill_brewer(palette="Set1") + scale_color_brewer(palette="Set1") + facet_wrap(variable~sal_group, ncol=6, scales="free") + theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 5), axis.text.y = element_text(size=5), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=5),legend.title=element_text(size=5),legend.text=element_text(size=5)) +  theme(legend.key.size = unit(0.3, "cm")) + labs(title="Water env. histogram plot", x ="Value", y = "Count", fill = "habitat")
ggsave("reads_hist_raw_wat_habitat_salg3_log.pdf")

for (h in 1:length(unique(mui$variable)))
{
er<-subset(mui, variable==unique(mui$variable)[h])
ermui2 <-subset(hmui2, variable==unique(mui$variable)[h])
ggplot(er, aes(x=value, fill=season)) +
geom_histogram(aes(y=..density..), position="identity", alpha=0.5, bins=40, colour="black", size=0.03) + geom_vline(data=ermui2, aes(xintercept=grp.mean, color=season), linetype="dashed") + theme_classic() + scale_color_manual(values=c("#00FF33","#FF6600")) + scale_fill_manual(values=c("#00FF33","#FF6600")) + scale_x_continuous(labels = comma) + scale_y_continuous(labels = comma) + facet_wrap(~habitat, ncol=1) + theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 7), axis.text.y = element_text(size=7), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=7),legend.title=element_text(size=7),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm")) + labs(title=paste(unique(mui$variable)[h],"env. histogram plot"), x ="Value", y = "Count", fill = "season")
ggsave(paste(unique(mui$variable)[h],"_raw_wat_habitat_season_log.pdf",sep=""))
}

for (h in 1:length(unique(mui$variable)))
{
er<-subset(mui, variable==unique(mui$variable)[h])
ermui2 <-subset(hmui3, variable==unique(mui$variable)[h])
ggplot(er, aes(x=value, fill=factor(sal_group))) +
geom_histogram(aes(y=..density..), position="identity", alpha=0.5, bins=40, colour="black", size=0.03) + geom_vline(data=ermui2, aes(xintercept=grp.mean, color=factor(sal_group)), linetype="dashed") + theme_classic() + scale_x_continuous(labels = comma) + scale_y_continuous(labels = comma) + scale_fill_manual(values=c("#00CCFF","#339966","#999900")) + scale_color_manual(values=c("#00CCFF","#339966","#999900")) + facet_wrap(~habitat, ncol=1) + theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 7), axis.text.y = element_text(size=7), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=7),legend.title=element_text(size=7),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm")) + labs(title=paste(unique(mui$variable)[h],"env. histogram plot"), x ="Value", y = "Count", fill = "sal_group")
ggsave(paste(unique(mui$variable)[h],"_raw_wat_habitat_salg3_log.pdf",sep=""))
}


#Sediment
setwd("/environmental_data/results/sed")

#incorporate salinity grouping
fsed$habitat<-as.character(fsed$habitat)
fsed$hb<-ifelse(fsed$habitat=="EB", "eelgrass", ifelse(fsed$habitat=="RB", "rocks", ifelse(fsed$habitat=="SB", "sand", fsed$habitat)))
fsed$snch<-paste(fsed$season, fsed$cluster, fsed$hb, sep="_")
fsed$Salinity<-fwat$Salinity[match(fsed$snch, fwat$snch)]

fsed$sal_group<-ifelse(fsed$Salinity<15,1,ifelse(fsed$Salinity<25,2,3))

#Test log-transform IC, OC, CN_ratio, NP_ratio, TP, N and C
vsed<-c("N", "C", "Organic_content", "Inorganic_content", "CN_ratio", "NP_ratio", "TP")

for (g in 1:length(vsed))
{
fscn<-ncol(fsed)
fsed[,fscn+1]<-log(fsed[,vsed[g]])
colnames(fsed)[fscn+1]<-paste("log", vsed[g], sep="_")
}

#Test transformed d14N_15N min value summed
fsed$d14N_15N_summedmin<-fsed$d14N_15N+26.94427829
fsed$cube_d14N_15N<-fsed$d14N_15N_summedmin^3

p_smui <- melt(fsed, id=c("season", "cluster", "habitat", "field_replicate", "sal_group"), measure=c("d14N_15N", "d12C_13C", "CN_ratio", "Density", "Dry_matter", "Grain_size", "Water_content", "TP","NP_ratio", "log_N", "log_C", "log_Organic_content", "log_Inorganic_content", "log_CN_ratio", "log_NP_ratio", "log_TP","cube_d14N_15N"))

#Remove NAs
smui<-na.omit(p_smui)

#Calculate mean
smui2 <- ddply(smui, .(variable, season), summarise, grp.mean=mean(na.omit(value)))
smui3 <- ddply(smui, .(variable, sal_group), summarise, grp.mean=mean(na.omit(value)))

smui$season <- factor(smui$season, levels = c("spring","autumn"))
smui$habitat <- factor(smui$habitat, levels = c("RB","SB","EB"))
smui2$season <- factor(smui2$season, levels = c("spring","autumn"))

ggplot(smui, aes(x=value, fill=season)) +
geom_histogram(aes(y=..density..), position="identity", alpha=0.5, bins=50, colour="black", size=0.03) + geom_vline(data=smui2, aes(xintercept=grp.mean, color=season), linetype="dashed") + theme_classic() + scale_color_manual(values=c("#00FF33","#FF6600")) + scale_fill_manual(values=c("#00FF33","#FF6600")) + scale_x_continuous(labels = comma) + scale_y_continuous(labels = comma) + facet_wrap(~variable, ncol=3, scales="free") + theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 7), axis.text.y = element_text(size=7), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=7),legend.title=element_text(size=7),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm")) + labs(title="Sediment env. histogram plot", x ="Value", y = "Count", fill = "season")
ggsave("reads_hist_raw_sed_season_log.pdf")

ggplot(smui, aes(x=value, fill=as.factor(sal_group))) +
geom_histogram(aes(y=..density..), position="identity", alpha=0.5, bins=50, colour="black", size=0.03) + geom_vline(data=smui3, aes(xintercept=grp.mean, color=as.factor(sal_group)), linetype="dashed") + theme_classic() + scale_x_continuous(labels = comma) + scale_y_continuous(labels = comma) + scale_fill_manual(values=c("#00CCFF","#339966","#999900")) + scale_color_manual(values=c("#00CCFF","#339966","#999900")) + facet_wrap(~variable, ncol=3, scales="free") + theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 7), axis.text.y = element_text(size=7), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=7),legend.title=element_text(size=7),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm")) + labs(title="Sediment env. histogram plot", x ="Value", y = "Count", fill = "sal_group")
ggsave("reads_hist_raw_sed_salg3_log.pdf")

##Habitat-wise
#Calculate mean
shmui2 <- ddply(smui, .(variable, season, habitat), summarise, grp.mean=mean(na.omit(value)))
shmui3 <- ddply(smui, .(variable, sal_group, habitat), summarise, grp.mean=mean(na.omit(value)))

shmui2$season <- factor(shmui2$season, levels = c("spring","autumn"))
shmui2$habitat <- factor(shmui2$habitat, levels = c("RB","SB","EB"))
shmui3$habitat <- factor(shmui3$habitat, levels = c("RB","SB","EB"))

ggplot(smui, aes(x=value, fill=habitat)) +
geom_histogram(aes(y=..density..), position="identity", alpha=0.4, bins=25) + geom_vline(data=shmui2, aes(xintercept=grp.mean, color=habitat), linetype="dashed") + geom_density(aes(colour=habitat), fill=NA) +theme_classic() + scale_x_continuous(labels = comma) + scale_y_continuous(labels = comma) + scale_color_brewer(palette="Set1") + scale_fill_brewer(palette="Set1") + facet_wrap(variable~season, ncol=4, scales="free") + theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 5), axis.text.y = element_text(size=5), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=7),legend.title=element_text(size=7),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm")) + labs(title="Sediment env. histogram plot", x ="Value", y = "Count", fill = "habitat")
ggsave("reads_hist_raw_sed_habitat_season_log.pdf")

ggplot(smui, aes(x=value, fill=habitat)) +
geom_histogram(aes(y=..density..), position="identity", alpha=0.4, bins=25) + geom_vline(data=shmui3, aes(xintercept=grp.mean, color=habitat), linetype="dashed") + geom_density(aes(colour=habitat), fill=NA) +theme_classic() + scale_x_continuous(labels = comma) + scale_y_continuous(labels = comma) + scale_color_brewer(palette="Set1") + scale_fill_brewer(palette="Set1") + facet_wrap(variable~sal_group, ncol=6, scales="free") + theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 4), axis.text.y = element_text(size=4), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=4),legend.title=element_text(size=5),legend.text=element_text(size=4)) +  theme(legend.key.size = unit(0.3, "cm")) + labs(title="Sediment env. histogram plot", x ="Value", y = "Count", fill = "habitat")
ggsave("reads_hist_raw_sed_habitat_salg3_log.pdf")



for (h in 1:length(unique(smui$variable)))
{
ser<-subset(smui, variable==unique(smui$variable)[h])
sermui2 <-subset(shmui2, variable==unique(smui$variable)[h])
ggplot(ser, aes(x=value, fill=season)) +
geom_histogram(aes(y=..density..), position="identity", alpha=0.5, bins=40, colour="black", size=0.03) + geom_vline(data=sermui2, aes(xintercept=grp.mean, color=season), linetype="dashed") + theme_classic() + scale_color_manual(values=c("#00FF33","#FF6600")) + scale_fill_manual(values=c("#00FF33","#FF6600")) + scale_x_continuous(labels = comma) + scale_y_continuous(labels = comma) + facet_wrap(~habitat, ncol=1) + theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 7), axis.text.y = element_text(size=7), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=7),legend.title=element_text(size=7),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm")) + labs(title=paste(unique(smui$variable)[h],"env. histogram plot"), x ="Value", y = "Count", fill = "season")
ggsave(paste(unique(smui$variable)[h],"_raw_sed_habitat_season_log.pdf",sep=""))
}

for (h in 1:length(unique(smui$variable)))
{
ser<-subset(smui, variable==unique(smui$variable)[h])
sermui2 <-subset(shmui3, variable==unique(smui$variable)[h])
ggplot(ser, aes(x=value, fill=factor(sal_group))) +
geom_histogram(aes(y=..density..), position="identity", alpha=0.5, bins=40, colour="black", size=0.03) + geom_vline(data=sermui2, aes(xintercept=grp.mean, color=factor(sal_group)), linetype="dashed") + theme_classic() + scale_x_continuous(labels = comma) + scale_y_continuous(labels = comma) + scale_fill_manual(values=c("#00CCFF","#339966","#999900")) + scale_color_manual(values=c("#00CCFF","#339966","#999900")) + facet_wrap(~habitat, ncol=1) + theme(axis.text.x = element_text(hjust = 1, vjust=0, size = 7), axis.text.y = element_text(size=7), strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")), strip.text = element_text(size=7),legend.title=element_text(size=7),legend.text=element_text(size=6)) +  theme(legend.key.size = unit(0.3, "cm")) + labs(title=paste(unique(smui$variable)[h],"env. histogram plot"), x ="Value", y = "Count", fill = "sal_group")
ggsave(paste(unique(smui$variable)[h],"_raw_sed_habitat_salg3_log.pdf",sep=""))
}


###Plots
#Water
setwd("/environmental_data/results/wat")

#Inspect and, if necessary, standardize data
head(fwat)
colnames(fwat)
ncol(fwat)
str(fwat)
fwat.c<-na.exclude(fwat)
fwat.c$log_DN<-log(fwat.c$NO2NO3+fwat.c$NH3)
fwat.c<-fwat.c[,c("Sample_ID","season","cluster","habitat","NO2NO3",
"NO2","NH3","PO4","Si","NO3",
"Latitude","Longitude","log_Chlorophyll","Salinity","Temperature",
"CTD_depth","d18O","SD_d18O","d12H","SD_d12H","log_NO2NO3","log_NO2",
"log_NH3","log_PO4","log_Si","log_NO3","log_DN","sal_group","snch","Time")]
fwat.c[,5:20]<-decostand(fwat.c[,5:20],"standardize")
head(fwat.c)
rownames(fwat.c)<-fwat.c$Sample_ID

compact_fwat<-fwat.c[,c("log_PO4","log_Chlorophyll","Salinity","Temperature","log_NO2",
"log_NH3","log_NO3","log_Si")]

##Remove outlier SAMPLES  2C4RW C15RW
compact_fwat<-compact_fwat[-63,]
compact_fwat<-compact_fwat[-83,]

compact_fwat2<-fwat.c[,c("Latitude","Longitude","log_Chlorophyll","Salinity","Temperature",
"CTD_depth","d18O","d12H","log_NO2NO3","log_NO2",
"log_NH3","log_PO4","log_Si","log_NO3","log_DN")]

##Remove outlier SAMPLES  2C4RW C15RW
compact_fwat2<-compact_fwat2[-63,]
compact_fwat2<-compact_fwat2[-83,]

#RDA biplot (PCA)
vareamb2<-rda(compact_fwat)
summary(vareamb2)

pdf("pca_wat_both_v.pdf")
biplot(vareamb2, scaling = 2, type = c("text", "points"))
dev.off()

#ggplot v.

smry <- summary(vareamb2)
pca.s  <- data.frame(smry$sites[,1:2])
pca.v <- data.frame(smry$species[,1:2])
pca.s$season<-fwat.c$season[match(rownames(pca.s), rownames(fwat.c))]
pca.s$habitat<-fwat.c$habitat[match(rownames(pca.s), rownames(fwat.c))]
pca.s$cluster<-fwat.c$cluster[match(rownames(pca.s), rownames(fwat.c))]
pca.s$sal_group<-fwat.c$sal_group[match(rownames(pca.s), rownames(fwat.c))]

pca.s$sal_group<-as.character(pca.s$sal_group)
pca.s$sal_group<-factor(pca.s$sal_group, levels = c("1","2","3"))
pca.s$season <- factor(pca.s$season, levels = c("spring","autumn"))
pca.s$habitat <- factor(pca.s$habitat, levels = c("RW","SW","EW"))
pca.s$fjord<-ifelse(pca.s$cluster==8|pca.s$cluster==9|pca.s$cluster==10|pca.s$cluster==12|pca.s$cluster==16|
pca.s$cluster==29, "fjord","open")

xlim = c(-3, 1.5) 
ylim = c(-1.5, 2.5)

rda.plot <- ggplot(data=pca.s, aes(x=PC1, y=PC2)) +
geom_point(data=pca.s, aes(x=PC1, y=PC2, color=habitat, shape=season), size=2.5, alpha=0.7) +
geom_text(data=pca.s, aes(label=cluster), size=1.8) +
scale_color_brewer(palette="Set1") +
scale_shape_manual(values=c(18, 17)) +
geom_hline(yintercept=0, linetype="dotted") +
geom_vline(xintercept=0, linetype="dotted") +
coord_fixed(ratio = diff(xlim)/diff(ylim), xlim=xlim, ylim=ylim, expand=F)

rda.plot +
geom_segment(data=pca.v, aes(x=0, xend=PC1, y=0, yend=PC2), size=0.3, color="black", arrow=arrow(length=unit(0.01,"npc"))) +
geom_text(data=pca.v, aes(x=PC1,y=PC2, label=rownames(pca.v), hjust=0.5*(1-sign(PC1)), vjust=0.5*(1-sign(PC2))), color="black", size=2.8) + theme_classic()

ggsave("pca_wat_both.pdf")

rda.plot <- ggplot(data=pca.s, aes(x=PC1, y=PC2)) +
geom_point(data=pca.s, aes(x=PC1, y=PC2, color=season, shape=as.factor(fjord)), size=2.5, alpha=0.7) +
geom_text(data=pca.s, aes(label=cluster), size=1.8) +
scale_color_manual(values=c("#00FF33","#FF6600")) +
scale_shape_manual(values=c(1, 4)) +
geom_hline(yintercept=0, linetype="dotted") +
geom_vline(xintercept=0, linetype="dotted") +
coord_fixed(ratio = diff(xlim)/diff(ylim), xlim=xlim, ylim=ylim, expand=F)

rda.plot +
geom_segment(data=pca.v, aes(x=0, xend=PC1, y=0, yend=PC2), size=0.3, color="black", arrow=arrow(length=unit(0.01,"npc"))) +
geom_text(data=pca.v, aes(x=PC1,y=PC2, label=rownames(pca.v), hjust=0.5*(1-sign(PC1)), vjust=0.5*(1-sign(PC2))), color="black", size=2.8) + theme_classic()

ggsave("pca_wat_both_season.pdf")


rda.plot <- ggplot(data=pca.s, aes(x=PC1, y=PC2)) +
geom_point(data=pca.s, aes(x=PC1, y=PC2, color=as.factor(sal_group), shape=season), size=2.5, alpha=0.7) +
geom_text(data=pca.s, aes(label=cluster), size=1.8) +
scale_color_manual(values=c("#00CCFF","#339966","#999900")) +
scale_shape_manual(values=c(18, 17)) +
geom_hline(yintercept=0, linetype="dotted") +
geom_vline(xintercept=0, linetype="dotted") +
coord_fixed(ratio = diff(xlim)/diff(ylim), xlim=xlim, ylim=ylim, expand=F)

rda.plot +
geom_segment(data=pca.v, aes(x=0, xend=PC1, y=0, yend=PC2), size=0.3, color="black", arrow=arrow(length=unit(0.01,"npc"))) +
geom_text(data=pca.v, aes(x=PC1,y=PC2, label=rownames(pca.v), hjust=0.5*(1-sign(PC1)), vjust=0.5*(1-sign(PC2))), color="black", size=2.8) + theme_classic()

ggsave("pca_wat_both_salg3.pdf")

rda.plot <- ggplot(data=pca.s, aes(x=PC1, y=PC2)) +
geom_point(data=pca.s, aes(x=PC1, y=PC2, color=as.factor(sal_group), shape=fjord), size=2.5, alpha=0.7) +
geom_text(data=pca.s, aes(label=cluster), size=1.8) +
scale_color_manual(values=c("#00CCFF","#339966","#999900")) +
scale_shape_manual(values=c(1, 4)) +
geom_hline(yintercept=0, linetype="dotted") +
geom_vline(xintercept=0, linetype="dotted") +
coord_fixed(ratio = diff(xlim)/diff(ylim), xlim=xlim, ylim=ylim, expand=F)

rda.plot +
geom_segment(data=pca.v, aes(x=0, xend=PC1, y=0, yend=PC2), size=0.3, color="black", arrow=arrow(length=unit(0.01,"npc"))) +
geom_text(data=pca.v, aes(x=PC1,y=PC2, label=rownames(pca.v), hjust=0.5*(1-sign(PC1)), vjust=0.5*(1-sign(PC2))), color="black", size=2.8) + theme_classic()

ggsave("pca_wat_both_salg3_fjord.pdf")


#Corplot
p<-cor(compact_fwat2, method = c("spearman"))
testRes<-cor.mtest(compact_fwat2, conf.level = 0.95)

pdf("corrplot_water.pdf")
corrplot(p, p.mat = testRes$p, method = 'circle', type = 'lower', insig='blank',
addCoef.col ='black', number.cex = 0.8, order = 'AOE', diag=FALSE)
dev.off()

#Sediment
setwd("/environmental_data/results/sed")
#Inspect and, if necessary, standardize data
head(fsed)
colnames(fsed)
fsed.c<-fsed[,c("Sample_ID","season","cluster","habitat","cube_d14N_15N", "d12C_13C", "log_CN_ratio", "Density", "Dry_matter", "Grain_size", "Water_content", "log_TP","log_NP_ratio", "log_N", "log_C", "log_Organic_content", "log_Inorganic_content","sal_group")] 
ncol(fsed.c)
str(fsed.c)
 
fsed.c<-na.exclude(fsed.c)
fsed.c[,5:18]<-decostand(fsed.c[,5:18],"standardize")
head(fsed.c)
rownames(fsed.c)<-fsed.c$Sample_ID

compact_fsed<-fsed.c[,c("cube_d14N_15N", "log_Organic_content","log_NP_ratio", "log_CN_ratio", "log_TP", "log_N","log_C")]

compact_fsed2<-fsed.c[,c("cube_d14N_15N", "d12C_13C", "log_CN_ratio", "Density", "Dry_matter", "Grain_size", "Water_content", "log_TP","log_NP_ratio", "log_N", "log_C", "log_Organic_content", "log_Inorganic_content")]


##Remove outlier samples - check PCA result below
#compact_fsed<-compact_fsed[-224,]
#compact_fsed<-compact_fsed[-308,]
#compact_fsed<-compact_fsed[-372,]


#RDA biplot
vareamb2<-rda(compact_fsed)
summary(vareamb2)

pdf("pca_sed_both_v.pdf")
biplot(vareamb2, scaling = 2, type = c("text", "points"))
dev.off()

#ggplot v.

smry <- summary(vareamb2)
pca.s  <- data.frame(smry$sites[,1:2])
pca.v <- data.frame(smry$species[,1:2])
pca.s$season<-fsed.c$season[match(rownames(pca.s), rownames(fsed.c))]
pca.s$habitat<-fsed.c$habitat[match(rownames(pca.s), rownames(fsed.c))]
pca.s$cluster<-fsed.c$cluster[match(rownames(pca.s), rownames(fsed.c))]
pca.s$hb<-ifelse(pca.s$habitat=="EB", "eelgrass", ifelse(pca.s$habitat=="RB", "rocks", ifelse(pca.s$habitat=="SB", "sand", pca.s$habitat)))
pca.s$snch<-paste(pca.s$season, pca.s$cluster, pca.s$hb, sep="_")
pca.s$Salinity<-fwat$Salinity[match(pca.s$snch, fwat$snch)]
pca.s$sal_group<-ifelse(pca.s$Salinity<15,1,ifelse(pca.s$Salinity<25,2,3))
pca.s$sal_group<-as.character(pca.s$sal_group)
pca.s$sal_group<-factor(pca.s$sal_group, levels = c("1","2","3"))
pca.s$season <- factor(pca.s$season, levels = c("spring","autumn"))
pca.s$habitat <- factor(pca.s$habitat, levels = c("RB","SB","EB"))
pca.s$fjord<-ifelse(pca.s$cluster==8|pca.s$cluster==9|pca.s$cluster==10|pca.s$cluster==12|pca.s$cluster==16|
pca.s$cluster==29, "fjord","open")


cent2 <- aggregate(cbind(PC1,PC2)~snch,pca.s,mean)
colnames(cent2)[2:3]<-c("cPC1","cPC2")
segs <- merge(pca.s, setNames(cent2, c("snch","cPC1","cPC2")), by = "snch", sort = FALSE)
segs$sh<-paste(segs$season, segs$habitat)

xlim = c(-1.5, 3) 
ylim = c(-1.5, 2.5)

rda.plot <- ggplot(data=segs, aes(x=PC1, y=PC2)) +
geom_point(data=segs, aes(x=PC1, y=PC2, color=habitat, shape=season), size=0.5, alpha=0.7) +
geom_segment(aes(xend = cPC1, yend = cPC2, color=habitat), alpha=0.9, size=0.2) + geom_point(aes(cPC1, cPC2, color=habitat, shape=season), size = 2, alpha=0.7) +
scale_color_brewer(palette="Set1") +
geom_text(data=segs, aes(cPC1, cPC2, label=cluster), size=1.5) +
scale_shape_manual(values=c(18, 17)) +
geom_hline(yintercept=0, linetype="dotted") +
geom_vline(xintercept=0, linetype="dotted") +
coord_fixed(ratio = diff(xlim)/diff(ylim), xlim=xlim, ylim=ylim, expand=F)

rda.plot +
geom_segment(data=pca.v, aes(x=0, xend=PC1, y=0, yend=PC2), size=0.3, color="black", arrow=arrow(length=unit(0.01,"npc"))) +
geom_text(data=pca.v, aes(x=PC1,y=PC2, label=rownames(pca.v), hjust=0.5*(1-sign(PC1)), vjust=0.5*(1-sign(PC2))), color="black", size=2.8) + theme_classic()

ggsave("pca_sed_both.pdf")

#Corplot
p<-cor(compact_fsed2, method = c("spearman"))
testRes<-cor.mtest(compact_fsed2, conf.level = 0.95)

pdf("corrplot_sed.pdf")
corrplot(p, p.mat = testRes$p, method = 'circle', type = 'lower', insig='blank',
addCoef.col ='black', number.cex = 0.8, order = 'AOE', diag=FALSE)
dev.off()

###Now, inspect water-sediment correlations
#Filter variables with many NAs in sed table
#Summarize sed samples per cluster+season+hab
head(fsed)
head(fwat)

#add new label to sed
fsed$hab_uni<-ifelse(fsed$habitat=="EB", "eelgrass", ifelse(fsed$habitat=="RB", "rocks", "sand"))
p_mui_s <- melt(fsed, id=c("season", "cluster", "hab_uni"), measure=c("cube_d14N_15N", "d12C_13C", "log_N", "log_C", "log_CN_ratio","log_NP_ratio", "Dry_matter", "Grain_size", "log_Inorganic_content","log_Organic_content", "log_TP"))
p_mui_s$shc<-paste(p_mui_s$season, p_mui_s$cluster, p_mui_s$hab_uni, sep="_")

s_no_rep <- dcast(p_mui_s, shc ~ variable, mean, margins="value")
rownames(s_no_rep)<-s_no_rep$shc

#add new label to wat
fwat$habitat<-as.character(fwat$habitat)
fwat$hab_uni<-ifelse(fwat$habitat=="EW", "eelgrass", ifelse(fwat$habitat=="RW", "rocks", ifelse(fwat$habitat=="SW", "sand", fwat$habitat)))
fwat$shc<-paste(fwat$season, fwat$cluster, fwat$hab_uni, sep="_")
rownames(fwat)<-fwat$shc
head(fwat)
fwat$log_DN<-log(fwat$NO2NO3+fwat$NH3)

pc_s<-s_no_rep[,c("cube_d14N_15N", "d12C_13C", "log_N", "log_C", "log_CN_ratio", "Dry_matter", "Grain_size", "log_Inorganic_content","log_Organic_content", "log_TP","log_NP_ratio")]
pc_w<-fwat[,c("log_PO4","log_Si","log_Chlorophyll","Salinity","Temperature","log_NO2","log_NO3","log_NH3","log_DN","Time")]

#match sites
nrow(pc_s)
nrow(pc_w)

pc_w2<-pc_w[rownames(pc_w) %in% rownames(pc_s), ]
pc_s2<-pc_s[rownames(pc_s) %in% rownames(pc_w), ]

nrow(pc_s2)
nrow(pc_w2)

#Merging wat and sed data
pc_b<-merge(pc_s, pc_w, by="row.names", all=F)
rownames(pc_b)<-pc_b$Row.names
pc_b<-pc_b[,2:21]
colnames(pc_b)

#Corplot
p<-cor(pc_b, method = c("spearman"), use="complete")
testRes<-cor.mtest(pc_b, conf.level = 0.95)

setwd("/environmental_data/results/both")
pdf("corrplot_all_sub.pdf")
corrplot(p, p.mat = testRes$p, method = 'circle', type = 'lower', insig='blank',
addCoef.col ='black', number.cex = 0.6, tl.col = 'black', tl.cex=0.7, order = 'hclust', diag=FALSE)
dev.off()


####

pc_bs<-pc_b
pc_bs$season<-fwat$season[match(rownames(pc_bs), rownames(fwat))]
pc_bs$hab_uni<-fwat$hab_uni[match(rownames(pc_bs), rownames(fwat))]
pc_bs$sal_group<-ifelse(pc_bs$Salinity<15,1,ifelse(pc_bs$Salinity<25,2,3))

colnames(pc_bs)
oqbd<-pc_bs[,c("cube_d14N_15N", "d12C_13C", "log_N", "log_C", "log_CN_ratio","log_NP_ratio", "Grain_size", "log_Inorganic_content", "log_TP","log_PO4","log_Si",
"log_Chlorophyll","Salinity","Temperature",
"log_NO2","log_NO3","log_NH3","log_DN","season")]
colnames(oqbd)<-c("cube_d14N_15N", "d12C_13C", "log_N", "log_C", "log_CN_ratio","log_NP_ratio", "Grain_size", "log_Inorganic_content", "log_TP","log_PO4","log_Si",
"Chla","Sal","Temp",
"log_NO2","log_NO3","log_NH3","log_DN","season")

cols <- character(nrow(oqbd))
cols[] <- "black"
cols[oqbd$season == "autumn"] <- "red"
cols[oqbd$season == "spring"] <- "blue"

setwd("/environmental_data/results/both")
pdf("all_vars_pre_final_data_season.pdf")
pairs(oqbd[,1:18], gap=0.1, cex=0.5, col=cols)
dev.off()

colnames(pc_bs)
oqbd<-pc_bs[,c("cube_d14N_15N", "d12C_13C", "log_N", "log_C", "log_CN_ratio","log_NP_ratio", "Grain_size", "log_Inorganic_content", "log_TP","log_PO4","log_Si",
"log_Chlorophyll","Salinity","Temperature",
"log_NO2","log_NO3","log_NH3","log_DN","sal_group")]
colnames(oqbd)<-c("cube_d14N_15N", "d12C_13C", "log_N", "log_C", "log_CN_ratio","log_NP_ratio", "Grain_size", "log_Inorganic_content", "log_TP","log_PO4","log_Si",
"Chla","Sal","Temp",
"log_NO2","log_NO3","log_NH3","log_DN","sal_group")

cols <- character(nrow(oqbd))
cols[] <- "black"
cols[oqbd$sal_group ==1] <- "red"
cols[oqbd$sal_group ==2] <- "green"
cols[oqbd$sal_group ==3] <- "blue"

setwd("/environmental_data/results/both")
pdf("all_vars_pre_final_data_salg3.pdf")
pairs(oqbd[,1:18], gap=0.1, cex=0.5, col=cols)
dev.off()

oqbd<-pc_bs[,c("cube_d14N_15N", "d12C_13C", "log_N", "log_C", "log_CN_ratio","log_NP_ratio", "Grain_size", "log_Inorganic_content", "log_TP","log_PO4","log_Si",
"log_Chlorophyll","Salinity","Temperature",
"log_NO2","log_NO3","log_NH3","log_DN","hab_uni")]
colnames(oqbd)<-c("cube_d14N_15N", "d12C_13C", "log_N", "log_C", "log_CN_ratio","log_NP_ratio", "Grain_size", "log_Inorganic_content", "log_TP","log_PO4","log_Si",
"Chla","Sal","Temp",
"log_NO2","log_NO3","log_NH3","log_DN","hab_uni")


cols <- character(nrow(oqbd))
cols[] <- "black"
cols[oqbd$hab_uni == "rocks"] <- "red"
cols[oqbd$hab_uni == "sand"] <- "blue"
cols[oqbd$hab_uni == "eelgrass"] <- "green"

setwd("/environmental_data/results/both")
pdf("all_vars_pre_final_data_hab.pdf")
pairs(oqbd[,1:18], gap=0.1, cex=0.5, col=cols)
dev.off()

#Alternative ggplot with with transformed variables
pc_bs$shc<-rownames(pc_bs)
head(pc_bs)
pc_bs[,1:18]<-decostand(pc_bs[,1:18],"standardize")

pc_bs<-pc_bs[,c("cube_d14N_15N", "d12C_13C", "log_N", "log_C", "log_CN_ratio","log_NP_ratio", "Grain_size", "log_Inorganic_content", "log_TP","log_PO4","log_Si",
"log_Chlorophyll","Salinity","Temperature",
"log_DN","season", "hab_uni", "sal_group", "shc")]
colnames(pc_bs)<-c("cube_d14N_15N", "d12C_13C", "log_N", "log_C", "log_CN_ratio","log_NP_ratio", "Grain_size", "log_Inorganic_content", "log_TP","log_PO4","log_Si",
"Chla","Sal","Temp",
"log_DN","season", "hab_uni", "sal_group", "shc")

for (i in 1:15)
{
gtrets<-melt(pc_bs, id=c("season", "hab_uni", "sal_group", "shc", colnames(pc_bs)[i]), measure=c("cube_d14N_15N", "d12C_13C", "log_N", "log_C", "log_CN_ratio","log_NP_ratio", "Grain_size", "log_Inorganic_content", "log_TP","log_PO4","log_Si",
"Chla","Sal","Temp",
"log_DN"))
colnames(gtrets)[5]<-"wat_val"
gtrets$wat_var<-colnames(pc_bs)[i]
if(i==1) {bft<-gtrets}
if(i!=1) {bft<-rbind(bft, gtrets)}
}

str(bft)

setwd("/environmental_data/results/both")
ggplot(bft, aes(x=wat_val, y=value, color=hab_uni)) + geom_point(size=0.3, alpha=0.4) + facet_grid(variable~wat_var, scales="free") + theme_bw() +
  theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(size=5, vjust=0.5), strip.text = element_text(size=5), legend.title=element_text(size=4), legend.text=element_text(size=4), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
ggsave("corr_wat_sed_point_habitat.pdf")

ggplot(bft, aes(x=wat_val, y=value, color=season)) + geom_point(size=0.3, alpha=0.4) + facet_grid(variable~wat_var, scales="free") + theme_bw() +
  theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(size=5, vjust=0.5), strip.text = element_text(size=5), legend.title=element_text(size=4), legend.text=element_text(size=4), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
ggsave("corr_wat_sed_point_season.pdf")

ggplot(bft, aes(x=wat_val, y=value, color=as.factor(sal_group))) + geom_point(size=0.3, alpha=0.4) + facet_grid(variable~wat_var, scales="free") + theme_bw() +
  theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(size=5, vjust=0.5), strip.text = element_text(size=5), legend.title=element_text(size=4), legend.text=element_text(size=4), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
ggsave("corr_wat_sed_point_salg3.pdf")


#PCA
#RDA biplot


pc_bs<-pc_b
pc_bs$season<-fwat$season[match(rownames(pc_bs), rownames(fwat))]
pc_bs$hab_uni<-fwat$hab_uni[match(rownames(pc_bs), rownames(fwat))]
pc_bs$sal_group<-ifelse(pc_bs$Salinity<15,1,ifelse(pc_bs$Salinity<25,2,3))
pc_bs$shc<-rownames(pc_bs)
pc_bs$log_CN_sum<-pc_bs$log_C+pc_bs$log_N



oqbd<-pc_bs[,c("d12C_13C","log_CN_ratio","log_CN_sum",
"log_TP","log_NP_ratio",
"log_Chlorophyll","Salinity","log_DN","log_Inorganic_content")]

oqbd<-na.exclude(oqbd)
head(oqbd)
oqbd<-decostand(oqbd,"standardize")

colnames(oqbd)
colnames(oqbd)<-c("d12C_13C","log_CN_ratio","log_CN_sum",
"log_TP","log_NP_ratio",
"log_Chla","Sal","log_DN","log_Inorganic_content")

vareamb2<-rda(oqbd)
summary(vareamb2)

pdf("pca_both_v.pdf")
biplot(vareamb2, scaling = 2, type = c("text", "points"))
dev.off()

#ggplot v.

smry <- summary(vareamb2)
pca.s  <- data.frame(smry$sites[,1:2])
pca.v <- data.frame(smry$species[,1:2])
pca.s$season<-merged_wat$sn[match(rownames(pca.s), merged_wat$snch)]
pca.s$habitat<-merged_wat$hb[match(rownames(pca.s), merged_wat$snch)]
pca.s$cluster<-merged_wat$cl[match(rownames(pca.s), merged_wat$snch)]
pca.s$hb<-ifelse(pca.s$habitat=="EB", "eelgrass", ifelse(pca.s$habitat=="RB", "rocks", ifelse(pca.s$habitat=="SB", "sand", pca.s$habitat)))
pca.s$snch<-paste(pca.s$season, pca.s$cluster, pca.s$hb, sep="_")
pca.s$Salinity<-fwat$Salinity[match(pca.s$snch, fwat$snch)]
pca.s$sal_group<-ifelse(pca.s$Salinity<15,1,ifelse(pca.s$Salinity<25,2,3))
pca.s$sh<-paste(pca.s$season, pca.s$habitat)
pca.s$sal_group<-as.character(pca.s$sal_group)
pca.s$sal_group<-factor(pca.s$sal_group, levels = c("1","2","3"))
pca.s$season <- factor(pca.s$season, levels = c("spring","autumn"))
pca.s$hb <- factor(pca.s$hb, levels = c("rocks","sand","eelgrass"))
pca.s$fjord<-ifelse(pca.s$cluster==8|pca.s$cluster==9|pca.s$cluster==10|pca.s$cluster==12|pca.s$cluster==16|
pca.s$cluster==29, "fjord","open")


xlim = c(-1.8, 2.2) 
ylim = c(-1.8, 1.8)

rda.plot <- ggplot(data=pca.s, aes(x=PC1, y=PC2)) +
geom_point(data=pca.s, aes(x=PC1, y=PC2, color=hb, shape=season), size=2.5, alpha=0.7) +
geom_text(data=pca.s, aes(label=cluster), size=1.8) +
geom_hline(yintercept=0, linetype="dotted") +
geom_vline(xintercept=0, linetype="dotted") +
scale_color_brewer(palette="Set1") +
scale_shape_manual(values=c(18, 17)) +
coord_fixed(ratio = diff(xlim)/diff(ylim), xlim=xlim, ylim=ylim, expand=F)

rda.plot +
geom_segment(data=pca.v, aes(x=0, xend=PC1, y=0, yend=PC2), size=0.3, color="black", arrow=arrow(length=unit(0.01,"npc"))) +
geom_text(data=pca.v, aes(x=PC1,y=PC2, label=rownames(pca.v), hjust=0.5*(1-sign(PC1)), vjust=0.5*(1-sign(PC2))), color="black", size=2.8) + theme_classic()

ggsave("pca_both.pdf")

rda.plot <- ggplot(data=pca.s, aes(x=PC1, y=PC2)) +
geom_point(data=pca.s, aes(x=PC1, y=PC2, color=season, shape=fjord), size=2.5, alpha=0.7) +
geom_text(data=pca.s, aes(label=cluster), size=1.8) +
geom_hline(yintercept=0, linetype="dotted") +
geom_vline(xintercept=0, linetype="dotted") +
scale_color_manual(values=c("#00FF33","#FF6600")) +
scale_shape_manual(values=c(1, 4)) +
coord_fixed(ratio = diff(xlim)/diff(ylim), xlim=xlim, ylim=ylim, expand=F)

rda.plot +
geom_segment(data=pca.v, aes(x=0, xend=PC1, y=0, yend=PC2), size=0.3, color="black", arrow=arrow(length=unit(0.01,"npc"))) +
geom_text(data=pca.v, aes(x=PC1,y=PC2, label=rownames(pca.v), hjust=0.5*(1-sign(PC1)), vjust=0.5*(1-sign(PC2))), color="black", size=2.8) + theme_classic()

ggsave("pca_both_season.pdf")

rda.plot <- ggplot(data=pca.s, aes(x=PC1, y=PC2)) +
geom_point(data=pca.s, aes(x=PC1, y=PC2, color=as.factor(sal_group), shape=season), size=2.5, alpha=0.7) +
geom_text(data=pca.s, aes(label=cluster), size=1.8) +
geom_hline(yintercept=0, linetype="dotted") +
geom_vline(xintercept=0, linetype="dotted") +
scale_color_manual(values=c("#00CCFF","#339966","#999900")) +
scale_shape_manual(values=c(18, 17)) +
coord_fixed(ratio = diff(xlim)/diff(ylim), xlim=xlim, ylim=ylim, expand=F)
rda.plot +
geom_segment(data=pca.v, aes(x=0, xend=PC1, y=0, yend=PC2), size=0.3, color="black", arrow=arrow(length=unit(0.01,"npc"))) +
geom_text(data=pca.v, aes(x=PC1,y=PC2, label=rownames(pca.v), hjust=0.5*(1-sign(PC1)), vjust=0.5*(1-sign(PC2))), color="black", size=2.8) + theme_classic()

ggsave("pca_both_salg3.pdf")

rda.plot <- ggplot(data=pca.s, aes(x=PC1, y=PC2)) +
geom_point(data=pca.s, aes(x=PC1, y=PC2, color=as.factor(sal_group), shape=fjord), size=2.5, alpha=0.7) +
geom_text(data=pca.s, aes(label=cluster), size=1.8) +
geom_hline(yintercept=0, linetype="dotted") +
geom_vline(xintercept=0, linetype="dotted") +
scale_color_manual(values=c("#00CCFF","#339966","#999900")) +
scale_shape_manual(values=c(1, 4)) +
coord_fixed(ratio = diff(xlim)/diff(ylim), xlim=xlim, ylim=ylim, expand=F)
rda.plot +
geom_segment(data=pca.v, aes(x=0, xend=PC1, y=0, yend=PC2), size=0.3, color="black", arrow=arrow(length=unit(0.01,"npc"))) +
geom_text(data=pca.v, aes(x=PC1,y=PC2, label=rownames(pca.v), hjust=0.5*(1-sign(PC1)), vjust=0.5*(1-sign(PC2))), color="black", size=2.8) + theme_classic()

ggsave("pca_both_salg3_fjord.pdf")


##Saving file
#add new label to sed
fsed$hab_uni<-ifelse(fsed$habitat=="EB", "eelgrass", ifelse(fsed$habitat=="RB", "rocks", "sand"))
p_mui_s <- melt(fsed, id=c("season", "cluster", "hab_uni"), measure=c("cube_d14N_15N", "d12C_13C", "log_N", "log_C", "log_CN_ratio","log_NP_ratio", "Dry_matter", "Grain_size", "Water_content", "log_Organic_content", "log_Inorganic_content", "log_TP","Organic_content","Inorganic_content","C","N"))
p_mui_s$shc<-paste(p_mui_s$season, p_mui_s$cluster, p_mui_s$hab_uni, sep="_")

s_no_rep <- dcast(p_mui_s, shc ~ variable, mean, margins="value")
rownames(s_no_rep)<-s_no_rep$shc

#add new label to wat
fwat$habitat<-as.character(fwat$habitat)
fwat$hab_uni<-ifelse(fwat$habitat=="EW", "eelgrass", ifelse(fwat$habitat=="RW", "rocks", ifelse(fwat$habitat=="SW", "sand", fwat$habitat)))
fwat$shc<-paste(fwat$season, fwat$cluster, fwat$hab_uni, sep="_")
rownames(fwat)<-fwat$shc
head(fwat)
fwat$log_DN<-log(fwat$NO2NO3+fwat$NH3)
fwat$DN<-fwat$NO2NO3+fwat$NH3

pc_s<-s_no_rep[,c("cube_d14N_15N", "d12C_13C", "log_N", "log_C", "log_CN_ratio","log_NP_ratio", "Dry_matter", "Grain_size", "Water_content", "log_Organic_content", "log_Inorganic_content", "log_TP","Organic_content","Inorganic_content","C","N")]
pc_w<-fwat[,c("log_PO4","log_Si","log_Chlorophyll","Salinity","Temperature",
"log_NO2","log_NO3","log_NH3","log_DN","Time","NO2","NO3","NH3","DN","PO4","Si")]

#match sites
nrow(pc_s)
nrow(pc_w)

pc_w2<-pc_w[rownames(pc_w) %in% rownames(pc_s), ]
pc_s2<-pc_s[rownames(pc_s) %in% rownames(pc_w), ]

nrow(pc_s2)
nrow(pc_w2)

#Merging wat and sed data
pc_b<-merge(pc_s, pc_w, by="row.names", all=F)
rownames(pc_b)<-pc_b$Row.names
pc_b<-pc_b[,-1]

pc_s<-s_no_rep[,c("cube_d14N_15N", "d12C_13C", "log_N", "log_C", "log_CN_ratio","log_NP_ratio", 
"Dry_matter", "Grain_size", "Water_content", "log_Organic_content","log_Inorganic_content", "log_TP","Organic_content","Inorganic_content","C","N")]
pc_w<-fwat[,c("log_PO4","log_Si","log_Chlorophyll","Salinity","Temperature",
"log_NO2","log_NO3","log_NH3","log_DN","Time","Latitude","Longitude",
"NO2","NO3","NH3","DN","PO4","Si")]

#match sites
nrow(pc_s)
nrow(pc_w)

pc_w2<-pc_w[rownames(pc_w) %in% rownames(pc_s), ]
pc_s2<-pc_s[rownames(pc_s) %in% rownames(pc_w), ]

nrow(pc_s2)
nrow(pc_w2)

#Merging wat and sed data
pc_b<-merge(pc_s, pc_w, by="row.names", all=F)
rownames(pc_b)<-pc_b$Row.names
pc_b<-pc_b[,-1]
pc_bs<-pc_b
pc_bs$season<-fwat$season[match(rownames(pc_bs), rownames(fwat))]
pc_bs$habitat<-fwat$hab_uni[match(rownames(pc_bs), rownames(fwat))]
pc_bs$sal_group<-ifelse(pc_bs$Salinity<15,1,ifelse(pc_bs$Salinity<25,2,3))

colnames(pc_bs)
pc_bs$cluster<-as.integer(gsub('\\D','', rownames(pc_bs)))

setwd("/environmental_data/results/both")
write.table(data.frame(pc_bs), "merged_metadata.txt", sep="\t", quote=FALSE, row.names=TRUE)

setwd("/environmental_data/results/sed")
write.table(data.frame(fsed), "sed_metadata.txt", sep="\t", quote=FALSE, row.names=TRUE)

setwd("/environmental_data/results/wat")
write.table(data.frame(fwat), "wat_metadata.txt", sep="\t", quote=FALSE, row.names=TRUE)

#reading files
setwd("/environmental_data/results/both")
pc_bs<-read.table("merged_metadata.txt", sep="\t", header=T)
setwd("/environmental_data/results/sed")
fsed<-read.table("sed_metadata.txt", sep="\t", header=T)
setwd("/environmental_data/results/wat")
fwat<-read.table("wat_metadata.txt", sep="\t", header=T)
fwat<-subset(fwat,!cluster==2) 

pc_bs$fjord<-ifelse(pc_bs$cluster==8|pc_bs$cluster==9|pc_bs$cluster==10|pc_bs$cluster==12|pc_bs$cluster==16|
pc_bs$cluster==29, "fjord","open")


#Plot env. vs. salinity, discreetly, separately for water and sediment without salinity grouping

htdj_s <- melt(fsed, id=c("season","hab_uni","cluster"), measure=c("cube_d14N_15N", "d12C_13C", "log_N", "log_C", "log_CN_ratio","log_NP_ratio", "Water_content", "log_Organic_content", "log_Inorganic_content", "log_TP"))
htdj_s$fjord<-ifelse(htdj_s$cluster==8|htdj_s$cluster==9|htdj_s$cluster==10|htdj_s$cluster==12|htdj_s$cluster==16|htdj_s$cluster==29, "fjord","open")

htdj_s$id<-paste(htdj_s$season, htdj_s$hab_uni, htdj_s$cluster, sep="_")
htdj_s2 <- ddply(htdj_s, .(variable, season, hab_uni), dplyr::summarize, grp.mean=mean(na.omit(value)),grp.sd=sd(na.omit(value)),
n_clusters=length(unique(cluster)),n_sites=length(na.omit(value)), position_y=min(na.omit(value))-0.5*(max(na.omit(value))))

htdj_s3<- ddply(htdj_s2, .(variable), dplyr::mutate, position_y2=min(position_y))

htdj_s3$label_l<-paste(htdj_s3$n_sites,htdj_s3$n_clusters,sep=", ")
htdj_s3$hab_uni <- factor(htdj_s3$hab_uni, levels = c("rocks","eelgrass","sand"))
htdj_s3$season <- factor(htdj_s3$season, levels = c("spring","autumn"))
htdj_s$season <- factor(htdj_s$season, levels = c("spring","autumn"))
htdj_s$hab_uni <- factor(htdj_s$hab_uni, levels = c("rocks","eelgrass","sand"))

env_sed_plot <- ggplot(htdj_s, aes(hab_uni, value, fill=hab_uni)) +
geom_boxplot(position = position_dodge(width = 0.7), alpha=0.7, outlier.shape=NA, size=0.1, width = 0.5) +
scale_fill_manual(values = c("yellowgreen", "cornflowerblue","thistle3"),
                     name = "", breaks = c("eelgrass","rocks", "sand"),
                     labels = c("eelgrass","rocks", "sand")) +
geom_point(data =htdj_s, aes(hab_uni, value, shape=fjord), alpha=0.7, size=0.6, stroke=0.1, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7)) +
geom_text(data = htdj_s3, aes(hab_uni, position_y2, label = label_l, group=hab_uni), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7), size = 1.5) +
facet_grid(variable~season, scales="free") +
theme_bw() +
scale_shape_manual(values=c(1, 4)) +
theme(strip.text = element_text(size=4),
legend.title=element_blank(),legend.text=element_text(size=6), axis.title=element_blank(), plot.title = element_text(size=5), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), panel.grid.major = element_line(colour = "gray80", size=0.1, linetype=2)) + guides(fill=guide_legend(title="Habitat"))
ggsave(env_sed_plot,filename="/summary_env_sed.pdf")


htdj_w <- melt(fwat, id=c("season","hab_uni","cluster"), measure=c("log_PO4","log_Si","log_Chlorophyll","Salinity","Temperature","log_NO2","log_NO3",
"log_NH3","log_DN","Latitude","Longitude"))
htdj_w$fjord<-ifelse(htdj_w$cluster==8|htdj_w$cluster==9|htdj_w$cluster==10|htdj_w$cluster==12|htdj_w$cluster==16|htdj_w$cluster==29, "fjord","open")

htdj_w$id<-paste(htdj_w$season, htdj_w$hab_uni, htdj_w$cluster, sep="_")
htdj_w2 <- ddply(htdj_w, .(variable, season, hab_uni), dplyr::summarize, grp.mean=mean(na.omit(value)),grp.sd=sd(na.omit(value)),
n_clusters=length(unique(cluster)), position_y=min(na.omit(value))-0.2*(max(na.omit(value))))

htdj_w3<- ddply(htdj_w2, .(variable), dplyr::mutate, position_y2=min(position_y))

htdj_w3$hab_uni <- factor(htdj_w3$hab_uni, levels = c("rocks","eelgrass","sand"))
htdj_w3$season <- factor(htdj_w3$season, levels = c("spring","autumn"))

htdj_w$season <- factor(htdj_w$season, levels = c("spring","autumn"))
htdj_w$hab_uni <- factor(htdj_w$hab_uni, levels = c("rocks","eelgrass","sand"))

env_wat_plot <- ggplot(htdj_w, aes(hab_uni, value, fill=hab_uni)) +
geom_boxplot(position = position_dodge(width = 0.7), alpha=0.7, outlier.shape=NA, size=0.1, width = 0.5) +
scale_fill_manual(values = c("yellowgreen", "cornflowerblue","thistle3"),
                     name = "", breaks = c("eelgrass","rocks", "sand"),
                     labels = c("eelgrass","rocks", "sand")) +
geom_point(data =htdj_w, aes(hab_uni, value, shape=fjord), alpha=0.7, size=0.6, stroke=0.1, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7)) +
geom_text(data = htdj_w3, aes(hab_uni, position_y2, label = n_clusters, group=hab_uni), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7), size = 1.5) +
facet_grid(variable~season, scales="free") +
theme_bw() +
scale_shape_manual(values=c(1, 4)) +
theme(strip.text = element_text(size=4),
legend.title=element_blank(),legend.text=element_text(size=6), axis.title=element_blank(), plot.title = element_text(size=5), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), panel.grid.major = element_line(colour = "gray80", size=0.1, linetype=2)) + guides(fill=guide_legend(title="Habitat"))
ggsave(env_wat_plot,filename="/summary_env_wat.pdf")



#Plot env. vs. salinity, fjords separately for water and sediment without salinity grouping

htdj_s <- melt(fsed, id=c("season","hab_uni","cluster"), measure=c("cube_d14N_15N", "d12C_13C", "log_N", "log_C", "log_CN_ratio","log_NP_ratio","Water_content", "log_Organic_content", "log_Inorganic_content", "log_TP"))
htdj_s$fjord<-ifelse(htdj_s$cluster==8|htdj_s$cluster==9|htdj_s$cluster==10|htdj_s$cluster==12|htdj_s$cluster==16|htdj_s$cluster==29, "fjord","open")

htdj_s$id<-paste(htdj_s$season, htdj_s$hab_uni, htdj_s$cluster, sep="_")
htdj_s2 <- ddply(htdj_s, .(variable, season, hab_uni, fjord), dplyr::summarize, grp.mean=mean(na.omit(value)),grp.sd=sd(na.omit(value)),
n_clusters=length(unique(cluster)),n_sites=length(na.omit(value)), position_y=min(na.omit(value))-0.5*(max(na.omit(value))))

htdj_s3<- ddply(htdj_s2, .(variable), dplyr::mutate, position_y2=min(position_y))

htdj_s3$label_l<-paste(htdj_s3$n_sites,htdj_s3$n_clusters,sep=", ")
htdj_s3$hab_uni <- factor(htdj_s3$hab_uni, levels = c("rocks","eelgrass","sand"))
htdj_s3$season <- factor(htdj_s3$season, levels = c("spring","autumn"))

htdj_s$season <- factor(htdj_s$season, levels = c("spring","autumn"))
htdj_s$hab_uni <- factor(htdj_s$hab_uni, levels = c("rocks","eelgrass","sand"))

setwd("/environmental_data/results/sed")

ggplot(htdj_s, aes(hab_uni, value, fill=fjord)) +
geom_boxplot(position = position_dodge(width = 0.7), alpha=0.7, outlier.shape=NA, size=0.1, width = 0.5) +
geom_point(data =htdj_s, aes(hab_uni, value, fill=fjord, shape=fjord), alpha=0.7, size=0.6, stroke=0.1, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7)) +
geom_text(data = htdj_s3, aes(hab_uni, position_y2, label = label_l, group=fjord), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7), size = 1.5) +
facet_grid(variable~season, scales="free") +
theme_bw() +
scale_fill_manual(values=c("#993300","#336699")) +
scale_shape_manual(values=c(1, 4)) +
theme(strip.text = element_text(size=4),
legend.title=element_text(size=5),legend.text=element_text(size=4), axis.title=element_text(size=5), plot.title = element_text(size=5), axis.text.x = element_text(size=4), axis.text.y = element_text(size=4), panel.grid.major = element_line(colour = "gray80", size=0.1, linetype=2)) + guides(fill=guide_legend(title="Habitat"))
ggsave("summary_env_sed_fjords.pdf")


htdj_w <- melt(fwat, id=c("season","hab_uni","cluster"), measure=c("log_PO4","log_Si","log_Chlorophyll","Salinity","Temperature","log_NO2","log_NO3",
"log_NH3","log_DN","Latitude","Longitude"))
htdj_w$fjord<-ifelse(htdj_w$cluster==8|htdj_w$cluster==9|htdj_w$cluster==10|htdj_w$cluster==12|htdj_w$cluster==16|htdj_w$cluster==29, "fjord","open")

htdj_w$id<-paste(htdj_w$season, htdj_w$hab_uni, htdj_w$cluster, sep="_")
htdj_w2 <- ddply(htdj_w, .(variable, season, hab_uni, fjord), dplyr::summarize, grp.mean=mean(na.omit(value)),grp.sd=sd(na.omit(value)),
n_clusters=length(unique(cluster)), position_y=min(na.omit(value))-0.2*(max(na.omit(value))))

htdj_w3<- ddply(htdj_w2, .(variable), dplyr::mutate, position_y2=min(position_y))

htdj_w3$hab_uni <- factor(htdj_w3$hab_uni, levels = c("rocks","eelgrass","sand"))
htdj_w3$season <- factor(htdj_w3$season, levels = c("spring","autumn"))

htdj_w$season <- factor(htdj_w$season, levels = c("spring","autumn"))
htdj_w$hab_uni <- factor(htdj_w$hab_uni, levels = c("rocks","eelgrass","sand"))

setwd("/environmental_data/results/wat")

ggplot(htdj_w, aes(hab_uni, value, fill=fjord)) +
geom_boxplot(position = position_dodge(width = 0.7), alpha=0.7, outlier.shape=NA, size=0.1, width = 0.5) +
geom_point(data =htdj_w, aes(hab_uni, value, fill=fjord, shape=fjord), alpha=0.7, size=0.6, stroke=0.1, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7)) +
geom_text(data = htdj_w3, aes(hab_uni, position_y2, label = n_clusters, group=fjord), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7), size = 1.5) +
facet_grid(variable~season, scales="free") +
theme_bw() +
scale_fill_manual(values=c("#993300","#336699")) +
scale_shape_manual(values=c(1, 4)) +
theme(strip.text = element_text(size=4),
legend.title=element_text(size=5),legend.text=element_text(size=4), axis.title=element_text(size=5), plot.title = element_text(size=5), axis.text.x = element_text(size=4), axis.text.y = element_text(size=4), panel.grid.major = element_line(colour = "gray80", size=0.1, linetype=2)) + guides(fill=guide_legend(title="Habitat"))
ggsave("summary_env_wat_fjords.pdf")

#Plot env. vs. salinity, discreetly, separately for water and sediment

htdj_s <- melt(fsed, id=c("season","hab_uni","sal_group","cluster"), measure=c("cube_d14N_15N", "d12C_13C", "log_N", "log_C", "log_CN_ratio","log_NP_ratio", "Water_content", "log_Organic_content", "log_Inorganic_content", "log_TP"))
htdj_s$fjord<-ifelse(htdj_s$cluster==8|htdj_s$cluster==9|htdj_s$cluster==10|htdj_s$cluster==12|htdj_s$cluster==16|htdj_s$cluster==29, "fjord","open")

htdj_s$id<-paste(htdj_s$season, htdj_s$hab_uni, htdj_s$cluster, sep="_")
htdj_s2 <- ddply(htdj_s, .(variable, sal_group, season, hab_uni), dplyr::summarize, grp.mean=mean(na.omit(value)),grp.sd=sd(na.omit(value)),
n_clusters=length(unique(cluster)),n_sites=length(na.omit(value)), position_y=min(na.omit(value))-0.5*(max(na.omit(value))))

htdj_s3<- ddply(htdj_s2, .(variable), dplyr::mutate, position_y2=min(position_y))

htdj_s3$label_l<-paste(htdj_s3$n_sites,htdj_s3$n_clusters,sep=", ")
htdj_s3$hab_uni <- factor(htdj_s3$hab_uni, levels = c("rocks","eelgrass","sand"))
htdj_s3$season <- factor(htdj_s3$season, levels = c("spring","autumn"))

htdj_s$sal_group<-as.character(htdj_s$sal_group)
htdj_s$season <- factor(htdj_s$season, levels = c("spring","autumn"))
htdj_s$hab_uni <- factor(htdj_s$hab_uni, levels = c("rocks","eelgrass","sand"))

setwd("/environmental_data/results/sed")

ggplot(htdj_s, aes(sal_group, value, fill=hab_uni)) +
geom_boxplot(position = position_dodge(width = 0.7), alpha=0.7, outlier.shape=NA, size=0.1, width = 0.5) +
geom_point(data =htdj_s, aes(sal_group, value, fill=hab_uni, shape=fjord), alpha=0.7, size=0.6, stroke=0.1, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7)) +
geom_text(data = htdj_s3, aes(sal_group, position_y2, label = label_l, group=hab_uni), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7), size = 1.5) +
facet_grid(variable~season, scales="free") +
theme_bw() +
scale_shape_manual(values=c(1, 4)) +
theme(strip.text = element_text(size=4),
legend.title=element_text(size=5),legend.text=element_text(size=4), axis.title=element_text(size=5), plot.title = element_text(size=5), axis.text.x = element_text(size=4), axis.text.y = element_text(size=4), panel.grid.major = element_line(colour = "gray80", size=0.1, linetype=2)) + guides(fill=guide_legend(title="Habitat"))
ggsave("summary_env_sed_sal_g.pdf")


htdj_w <- melt(fwat, id=c("season","hab_uni","sal_group","cluster"), measure=c("log_PO4","log_Si","log_Chlorophyll","Salinity","Temperature","log_NO2",
"log_NO3","log_NH3","log_DN","Latitude","Longitude"))
htdj_w$fjord<-ifelse(htdj_w$cluster==8|htdj_w$cluster==9|htdj_w$cluster==10|htdj_w$cluster==12|htdj_w$cluster==16|htdj_w$cluster==29, "fjord","open")

htdj_w$id<-paste(htdj_w$season, htdj_w$hab_uni, htdj_w$cluster, sep="_")
htdj_w2 <- ddply(htdj_w, .(variable, sal_group, season, hab_uni), dplyr::summarize, grp.mean=mean(na.omit(value)),grp.sd=sd(na.omit(value)),
n_clusters=length(unique(cluster)), position_y=min(na.omit(value))-0.2*(max(na.omit(value))))

htdj_w3<- ddply(htdj_w2, .(variable), dplyr::mutate, position_y2=min(position_y))

htdj_w3$hab_uni <- factor(htdj_w3$hab_uni, levels = c("rocks","eelgrass","sand"))
htdj_w3$season <- factor(htdj_w3$season, levels = c("spring","autumn"))

htdj_w$sal_group<-as.character(htdj_w$sal_group)
htdj_w$season <- factor(htdj_w$season, levels = c("spring","autumn"))
htdj_w$hab_uni <- factor(htdj_w$hab_uni, levels = c("rocks","eelgrass","sand"))

setwd("/environmental_data/results/wat")

ggplot(htdj_w, aes(sal_group, value, fill=hab_uni)) +
geom_boxplot(position = position_dodge(width = 0.7), alpha=0.7, outlier.shape=NA, size=0.1, width = 0.5) +
geom_point(data =htdj_w, aes(sal_group, value, fill=hab_uni, shape=fjord), alpha=0.7, size=0.6, stroke=0.1, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7)) +
geom_text(data = htdj_w3, aes(sal_group, position_y2, label = n_clusters, group=hab_uni), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7), size = 1.5) +
facet_grid(variable~season, scales="free") +
theme_bw() +
scale_shape_manual(values=c(1, 4)) +
theme(strip.text = element_text(size=4),
legend.title=element_text(size=5),legend.text=element_text(size=4), axis.title=element_text(size=5), plot.title = element_text(size=5), axis.text.x = element_text(size=4), axis.text.y = element_text(size=4), panel.grid.major = element_line(colour = "gray80", size=0.1, linetype=2)) + guides(fill=guide_legend(title="Habitat"))
ggsave("summary_env_wat_sal_g.pdf")


#Plot correlations among other variables in continuous scale
setwd("/environmental_data/results/both")
#Salinity
pc_bs$log_CN_sum<-log(pc_bs$C+pc_bs$N)

variables_ite<-c("Temperature","log_DN","log_Si", "log_Chlorophyll","log_NP_ratio","log_Inorganic_content",
"cube_d14N_15N","log_TP")

htdj <- melt(pc_bs, id=c("season","habitat","sal_group","cluster","fjord","Salinity"), measure=variables_ite)

htdj$id<-paste(htdj$season, htdj$habitat, sep="_")
htdj$sal_group<-as.character(htdj$sal_group)
htdj$season <- factor(htdj$season, levels = c("spring","autumn"))
htdj$habitat <- factor(htdj$habitat, levels = c("rocks","sand","eelgrass"))

setwd("/environmental_data/results/both")
cont_plot<-ggplot(htdj, aes(Salinity, value, color=habitat)) +
geom_point(data=htdj, aes(Salinity, value, shape=fjord, color=habitat), size=0.9, alpha=0.7) +
geom_smooth(data=htdj, aes(group=id, color=habitat), method="lm", size=0.2, se=F) +
facet_grid(variable~season, scale="free_y") +
scale_color_manual(values = c("yellowgreen", "cornflowerblue","thistle3"),
                     name = "", breaks = c("eelgrass","rocks", "sand"),
                     labels = c("eelgrass","rocks", "sand")) +
theme_bw() +
scale_shape_manual(values=c(1, 4)) +
theme(strip.text = element_text(size=3.5),
legend.title=element_blank(),legend.text=element_text(size=6), axis.title.x=element_text(size=6), axis.title.y=element_blank(),plot.title = element_text(size=8), axis.text.x = element_text(size=6), axis.text.y = element_text(size=6), panel.grid.major = element_line(colour = "gray80", size=0.3, linetype=2))
ggsave(cont_plot,filename="/Subset_cont_Salinity_hb.pdf")

htdj$id<-paste(htdj$season, htdj$fjord, sep="_")
htdj$sal_group<-as.character(htdj$sal_group)
htdj$season <- factor(htdj$season, levels = c("spring","autumn"))
htdj$habitat <- factor(htdj$habitat, levels = c("rocks","sand","eelgrass"))

setwd("/environmental_data/results/both")
ggplot(htdj, aes(Salinity, value, color=fjord)) +
geom_point(data=htdj, aes(Salinity, value, shape=fjord, color=fjord), size=0.9, alpha=0.7) +
geom_smooth(data=htdj, aes(group=id, color=fjord), method="lm", size=0.2, se=F) +
facet_grid(variable~season, scale="free_y") +
scale_color_manual(values=c("#993300","#336699")) +
theme_bw() +
scale_shape_manual(values=c(1, 4)) +
theme(strip.text = element_text(size=3.5),
legend.title=element_text(size=4),legend.text=element_text(size=4), axis.title=element_text(size=4), plot.title = element_text(size=8), axis.text.x = element_text(size=4), axis.text.y = element_text(size=4), panel.grid.major = element_line(colour = "gray80", size=0.3, linetype=2))
ggsave("Subset_cont_Salinity_fj.pdf")


setwd("/environmental_data/results/both")
#log_PO4
variables_ite<-c("log_Si", "log_Chlorophyll", "Temperature", "log_NO2", "log_DN")

htdj <- melt(pc_bs, id=c("season","habitat","sal_group","cluster","fjord","log_PO4"), measure=variables_ite)

htdj$id<-paste(htdj$season, htdj$fjord, sep="_")
htdj$sal_group<-as.character(htdj$sal_group)
htdj$season <- factor(htdj$season, levels = c("spring","autumn"))
htdj$habitat <- factor(htdj$habitat, levels = c("rocks","sand","eelgrass"))


setwd("/environmental_data/results/both")
ggplot(htdj, aes(log_PO4, value, color=fjord)) +
geom_point(data=htdj, aes(log_PO4, value, shape=fjord, color=fjord), size=0.9, alpha=0.7) +
geom_smooth(data=htdj, aes(group=id, color=fjord), method="lm", size=0.2, se=F) +
facet_grid(variable~season, scale="free_y") +
scale_color_manual(values=c("#993300","#336699")) +
theme_bw() +
scale_shape_manual(values=c(1, 4)) +
theme(strip.text = element_text(size=3.5),
legend.title=element_text(size=4),legend.text=element_text(size=4), axis.title=element_text(size=4), plot.title = element_text(size=8), axis.text.x = element_text(size=4), axis.text.y = element_text(size=4), panel.grid.major = element_line(colour = "gray80", size=0.3, linetype=2))
ggsave("Subset_cont_log_PO4_fj.pdf")

htdj$id<-paste(htdj$season, htdj$habitat, sep="_")
htdj$sal_group<-as.character(htdj$sal_group)
htdj$season <- factor(htdj$season, levels = c("spring","autumn"))
htdj$habitat <- factor(htdj$habitat, levels = c("rocks","sand","eelgrass"))

setwd("/environmental_data/results/both")
ggplot(htdj, aes(log_PO4, value, color=habitat)) +
geom_point(data=htdj, aes(log_PO4, value, shape=fjord, color=habitat), size=0.9, alpha=0.7) +
geom_smooth(data=htdj, aes(group=id, color=habitat), method="lm", size=0.2, se=F) +
facet_grid(variable~season, scale="free_y") +
scale_color_brewer(palette="Set1") +
theme_bw() +
scale_shape_manual(values=c(1, 4)) +
theme(strip.text = element_text(size=3.5),
legend.title=element_text(size=4),legend.text=element_text(size=4), axis.title=element_text(size=4), plot.title = element_text(size=8), axis.text.x = element_text(size=4), axis.text.y = element_text(size=4), panel.grid.major = element_line(colour = "gray80", size=0.3, linetype=2))
ggsave("Subset_cont_log_PO4_hb.pdf")



setwd("/environmental_data/results/both")
#CN_sum

pc_bs$log_CN_sum<-log(pc_bs$C+pc_bs$N)

variables_ite<-c("d12C_13C", "log_TP", "log_NP_ratio", "Water_content")

htdj <- melt(pc_bs, id=c("season","habitat","sal_group","cluster","fjord","log_CN_sum"), measure=variables_ite)

htdj$id<-paste(htdj$season, htdj$habitat, sep="_")
htdj$sal_group<-as.character(htdj$sal_group)
htdj$season <- factor(htdj$season, levels = c("spring","autumn"))
htdj$habitat <- factor(htdj$habitat, levels = c("rocks","sand","eelgrass"))

setwd("/environmental_data/results/both")
ggplot(htdj, aes(log_CN_sum, value, color=habitat)) +
geom_point(data=htdj, aes(log_CN_sum, value, shape=fjord, color=habitat), size=0.9, alpha=0.7) +
geom_smooth(data=htdj, aes(group=id, color=habitat), method="lm", size=0.2, se=F) +
facet_grid(variable~season, scale="free_y") +
scale_color_brewer(palette="Set1") +
theme_bw() +
scale_shape_manual(values=c(1, 4)) +
theme(strip.text = element_text(size=3.5),
legend.title=element_text(size=4),legend.text=element_text(size=4), axis.title=element_text(size=4), plot.title = element_text(size=8), axis.text.x = element_text(size=4), axis.text.y = element_text(size=4), panel.grid.major = element_line(colour = "gray80", size=0.3, linetype=2))
ggsave("Subset_cont_log_CN_sum_hb.pdf")

htdj$id<-paste(htdj$season, htdj$fjord, sep="_")
htdj$sal_group<-as.character(htdj$sal_group)
htdj$season <- factor(htdj$season, levels = c("spring","autumn"))
htdj$habitat <- factor(htdj$habitat, levels = c("rocks","sand","eelgrass"))

setwd("/environmental_data/results/both")
ggplot(htdj, aes(log_CN_sum, value, color=fjord)) +
geom_point(data=htdj, aes(log_CN_sum, value, shape=fjord, color=fjord), size=0.9, alpha=0.7) +
geom_smooth(data=htdj, aes(group=id, color=fjord), method="lm", size=0.2, se=F) +
facet_grid(variable~season, scale="free_y") +
scale_color_manual(values=c("#993300","#336699")) +
theme_bw() +
scale_shape_manual(values=c(1, 4)) +
theme(strip.text = element_text(size=3.5),
legend.title=element_text(size=4),legend.text=element_text(size=4), axis.title=element_text(size=4), plot.title = element_text(size=8), axis.text.x = element_text(size=4), axis.text.y = element_text(size=4), panel.grid.major = element_line(colour = "gray80", size=0.3, linetype=2))
ggsave("Subset_cont_log_CN_sum_fj.pdf")



setwd("/environmental_data/results/both")
#CN_ratio
variables_ite<-c("cube_d14N_15N", "log_N", "log_NP_ratio")

htdj <- melt(pc_bs, id=c("season","habitat","sal_group","cluster","fjord","log_CN_ratio"), measure=variables_ite)
htdj$id<-paste(htdj$season, htdj$habitat, sep="_")
htdj$sal_group<-as.character(htdj$sal_group)
htdj$season <- factor(htdj$season, levels = c("spring","autumn"))
htdj$habitat <- factor(htdj$habitat, levels = c("rocks","sand","eelgrass"))

setwd("/environmental_data/results/both")
ggplot(htdj, aes(log_CN_ratio, value, color=habitat)) +
geom_point(data=htdj, aes(log_CN_ratio, value, shape=fjord, color=habitat), size=0.9, alpha=0.7) +
geom_smooth(data=htdj, aes(group=id, color=habitat), method="lm", size=0.2, se=F) +
facet_grid(variable~season, scale="free_y") +
scale_color_brewer(palette="Set1") +
theme_bw() +
scale_shape_manual(values=c(1, 4)) +
theme(strip.text = element_text(size=3.5),
legend.title=element_text(size=4),legend.text=element_text(size=4), axis.title=element_text(size=4), plot.title = element_text(size=8), axis.text.x = element_text(size=4), axis.text.y = element_text(size=4), panel.grid.major = element_line(colour = "gray80", size=0.3, linetype=2))
ggsave("Subset_cont_log_log_CN_ratio_hb.pdf")

htdj$id<-paste(htdj$season, htdj$fjord, sep="_")
htdj$sal_group<-as.character(htdj$sal_group)
htdj$season <- factor(htdj$season, levels = c("spring","autumn"))
htdj$habitat <- factor(htdj$habitat, levels = c("rocks","sand","eelgrass"))

setwd("/environmental_data/results/both")
ggplot(htdj, aes(log_CN_ratio, value, color=fjord)) +
geom_point(data=htdj, aes(log_CN_ratio, value, shape=fjord, color=fjord), size=0.9, alpha=0.7) +
geom_smooth(data=htdj, aes(group=id, color=fjord), method="lm", size=0.2, se=F) +
facet_grid(variable~season, scale="free_y") +
scale_color_manual(values=c("#993300","#336699")) +
theme_bw() +
scale_shape_manual(values=c(1, 4)) +
theme(strip.text = element_text(size=3.5),
legend.title=element_text(size=4),legend.text=element_text(size=4), axis.title=element_text(size=4), plot.title = element_text(size=8), axis.text.x = element_text(size=4), axis.text.y = element_text(size=4), panel.grid.major = element_line(colour = "gray80", size=0.3, linetype=2))
ggsave("Subset_cont_log_CN_ratio_fj.pdf")


setwd("/environmental_data/results/both")
#Time
variables_ite<-c("Temperature", "Latitude", "Longitude", "Salinity")

htdj <- melt(pc_bs, id=c("season","habitat","sal_group","cluster","fjord","Time"), measure=variables_ite)

htdj$sal_group<-as.character(htdj$sal_group)
htdj$season <- factor(htdj$season, levels = c("spring","autumn"))
htdj$habitat <- factor(htdj$habitat, levels = c("rocks","sand","eelgrass"))
htdj$Time<-as.Date(htdj$Time, format= "%d-%m-%y")

setwd("/environmental_data/results/both")
ggplot(htdj, aes(Time, value)) +
geom_point(data=htdj, aes(Time, value, shape=fjord), size=0.9, alpha=0.7) +
geom_smooth(data=htdj, aes(group=season), method="lm", size=0.2, se=F, color="red") +
facet_grid(variable~season, scale="free") +
theme_bw() +
scale_shape_manual(values=c(1, 4)) +
scale_x_date(date_breaks = "15 days", date_labels = "%d \n %m") +
theme(strip.text = element_text(size=3.5),
legend.title=element_text(size=4),legend.text=element_text(size=4), axis.title=element_text(size=4), plot.title = element_text(size=8), axis.text.x = element_text(size=4), axis.text.y = element_text(size=4), panel.grid.major = element_line(colour = "gray80", size=0.3, linetype=2))
ggsave("Subset_cont_Time.pdf")









#Plot env. vs. salinity, continuously

htdj <- melt(pc_bs, id=c("season","habitat","sal_group","cluster","fjord","Salinity"), measure=c("log_PO4","log_Si","log_Chlorophyll","Temperature","log_NO2","log_NO3","log_NH3",
"log_DN","Latitude","Longitude","cube_d14N_15N", "d12C_13C", "log_N", "log_C", "log_CN_ratio","log_NP_ratio", "Dry_matter","Water_content", "log_Organic_content", "log_Inorganic_content", "log_TP"))

htdj$id<-paste(htdj$season, htdj$habitat, sep="_")
htdj2 <- ddply(htdj, .(variable, sal_group), dplyr::mutate,, grp.mean=mean(na.omit(value)),grp.sd=sd(na.omit(value)),n_clusters=length(unique(cluster)),n_sites=length(unique(id)))

htdj$sal_group<-as.character(htdj$sal_group)
htdj$season <- factor(htdj$season, levels = c("spring","autumn"))
htdj$habitat <- factor(htdj$habitat, levels = c("eelgrass","sand","rocks"))
htdj$gp <-paste(as.character(htdj$season),as.character(htdj$fjord),sep="_")

setwd("/environmental_data/results/both")
ggplot(htdj, aes(Salinity, value, color=season)) +
geom_point(data=htdj, aes(Salinity, value, shape=fjord, color=season), size=0.8) +
facet_wrap(variable~season, scales="free", ncol=6) +
scale_color_brewer(palette="Dark2") +
theme_bw() +
scale_shape_manual(values=c(1, 4)) +
theme(strip.text = element_text(size=3.5),
legend.title=element_text(size=4),legend.text=element_text(size=4), axis.title=element_text(size=4), plot.title = element_text(size=8), axis.text.x = element_text(size=4), axis.text.y = element_text(size=4), panel.grid.major = element_line(colour = "gray80", size=0.3, linetype=2))
ggsave("summary_env_sal_c.pdf")

#Plot env. vs. Longitude, continuously

htdj <- melt(pc_bs, id=c("season","habitat","sal_group","cluster","fjord","Longitude"), measure=c("log_PO4","log_Si","log_Chlorophyll","Temperature","log_NO2","log_NO3","log_NH3",
"log_DN","Latitude","Salinity","cube_d14N_15N", "d12C_13C", "log_N", "log_C", "log_CN_ratio","log_NP_ratio", "Dry_matter","Water_content", "log_Organic_content", "log_Inorganic_content", "log_TP"))

htdj$id<-paste(htdj$season, htdj$habitat, sep="_")
htdj2 <- ddply(htdj, .(variable, sal_group), dplyr::mutate,, grp.mean=mean(na.omit(value)),grp.sd=sd(na.omit(value)),n_clusters=length(unique(cluster)),n_sites=length(unique(id)))

htdj$sal_group<-as.character(htdj$sal_group)
htdj$season <- factor(htdj$season, levels = c("spring","autumn"))
htdj$habitat <- factor(htdj$habitat, levels = c("eelgrass","sand","rocks"))
htdj$gp <-paste(as.character(htdj$season),as.character(htdj$fjord),sep="_")

setwd("/environmental_data/results/both")
ggplot(htdj, aes(Longitude, value, color=season)) +
geom_point(data=htdj, aes(Longitude, value, shape=fjord, color=season), size=0.8) +
facet_wrap(variable~season, scales="free", ncol=6) +
scale_color_brewer(palette="Dark2") +
theme_bw() +
scale_shape_manual(values=c(1, 4)) +
theme(strip.text = element_text(size=3.5),
legend.title=element_text(size=4),legend.text=element_text(size=4), axis.title=element_text(size=4), plot.title = element_text(size=8), axis.text.x = element_text(size=4), axis.text.y = element_text(size=4), panel.grid.major = element_line(colour = "gray80", size=0.3, linetype=2))
ggsave("summary_env_long_c.pdf")

#Plot env. vs. Latitude, continuously

htdj <- melt(pc_bs, id=c("season","habitat","sal_group","cluster","fjord","Latitude"), measure=c("log_PO4","log_Si","log_Chlorophyll","Temperature","log_NO2","log_NO3",
"log_NH3","log_DN","Salinity","Longitude","cube_d14N_15N", "d12C_13C", "log_N", "log_C", "log_CN_ratio","log_NP_ratio", "Dry_matter","Water_content", "log_Organic_content", "log_Inorganic_content", "log_TP"))

htdj$id<-paste(htdj$season, htdj$habitat, sep="_")
htdj2 <- ddply(htdj, .(variable, sal_group), dplyr::mutate,, grp.mean=mean(na.omit(value)),grp.sd=sd(na.omit(value)),n_clusters=length(unique(cluster)),n_sites=length(unique(id)))

htdj$sal_group<-as.character(htdj$sal_group)
htdj$season <- factor(htdj$season, levels = c("spring","autumn"))
htdj$habitat <- factor(htdj$habitat, levels = c("eelgrass","sand","rocks"))
htdj$gp <-paste(as.character(htdj$season),as.character(htdj$fjord),sep="_")

setwd("/environmental_data/results/both")
ggplot(htdj, aes(Latitude, value, color=season)) +
geom_point(data=htdj, aes(Latitude, value, shape=fjord, color=season), size=0.8) +
facet_wrap(variable~season, scales="free", ncol=6) +
scale_color_brewer(palette="Dark2") +
theme_bw() +
scale_shape_manual(values=c(1, 4)) +
theme(strip.text = element_text(size=3.5),
legend.title=element_text(size=4),legend.text=element_text(size=4), axis.title=element_text(size=4), plot.title = element_text(size=8), axis.text.x = element_text(size=4), axis.text.y = element_text(size=4), panel.grid.major = element_line(colour = "gray80", size=0.3, linetype=2))
ggsave("summary_env_lat_c.pdf")

#Plot env. vs. PO4 (log10), continuously

htdj <- melt(pc_bs, id=c("season","habitat","sal_group","cluster","fjord","log_PO4"), measure=c("log_Si","log_Chlorophyll","Temperature","log_NO2","log_NO3",
"log_NH3","log_DN","Salinity","Longitude","cube_d14N_15N", "d12C_13C", "log_N", "log_C", "log_CN_ratio","log_NP_ratio", "Dry_matter","Water_content", "log_Organic_content", "log_Inorganic_content", "log_TP"))

htdj$id<-paste(htdj$season, htdj$habitat, sep="_")
htdj2 <- ddply(htdj, .(variable, sal_group), dplyr::mutate,, grp.mean=mean(na.omit(value)),grp.sd=sd(na.omit(value)),n_clusters=length(unique(cluster)),n_sites=length(unique(id)))

htdj$sal_group<-as.character(htdj$sal_group)
htdj$season <- factor(htdj$season, levels = c("spring","autumn"))
htdj$habitat <- factor(htdj$habitat, levels = c("eelgrass","sand","rocks"))
htdj$gp <-paste(as.character(htdj$season),as.character(htdj$fjord),sep="_")

setwd("/environmental_data/results/both")
ggplot(htdj, aes(log_PO4, value, color=season)) +
geom_point(data=htdj, aes(log_PO4, value, shape=fjord, color=season), size=0.8) +
facet_wrap(variable~season, scales="free", ncol=6) +
scale_color_brewer(palette="Dark2") +
theme_bw() +
scale_shape_manual(values=c(1, 4)) +
theme(strip.text = element_text(size=3.5),
legend.title=element_text(size=4),legend.text=element_text(size=4), axis.title=element_text(size=4), plot.title = element_text(size=8), axis.text.x = element_text(size=4), axis.text.y = element_text(size=4), panel.grid.major = element_line(colour = "gray80", size=0.3, linetype=2))
ggsave("summary_env_PO4_c.pdf")

#Plot env. vs. NO3 (log10), continuously

htdj <- melt(pc_bs, id=c("season","habitat","sal_group","cluster","fjord","log_NO3"), measure=c("log_Si","log_Chlorophyll","Temperature","log_NO2","log_PO4","log_NH3",
"log_DN","Salinity","Longitude","cube_d14N_15N", "d12C_13C", "log_N", "log_C", "log_CN_ratio","log_NP_ratio", "Dry_matter","Water_content", "log_Organic_content", "log_Inorganic_content", "log_TP"))

htdj$id<-paste(htdj$season, htdj$habitat, sep="_")
htdj2 <- ddply(htdj, .(variable, sal_group), dplyr::mutate,, grp.mean=mean(na.omit(value)),grp.sd=sd(na.omit(value)),n_clusters=length(unique(cluster)),n_sites=length(unique(id)))

htdj$sal_group<-as.character(htdj$sal_group)
htdj$season <- factor(htdj$season, levels = c("spring","autumn"))
htdj$habitat <- factor(htdj$habitat, levels = c("eelgrass","sand","rocks"))
htdj$gp <-paste(as.character(htdj$season),as.character(htdj$fjord),sep="_")

setwd("/environmental_data/results/both")
ggplot(htdj, aes(log_NO3, value, color=season)) +
geom_point(data=htdj, aes(log_NO3, value, shape=fjord, color=season), size=0.8) +
facet_wrap(variable~season, scales="free", ncol=6) +
scale_color_brewer(palette="Dark2") +
theme_bw() +
scale_shape_manual(values=c(1, 4)) +
theme(strip.text = element_text(size=3.5),
legend.title=element_text(size=4),legend.text=element_text(size=4), axis.title=element_text(size=4), plot.title = element_text(size=8), axis.text.x = element_text(size=4), axis.text.y = element_text(size=4), panel.grid.major = element_line(colour = "gray80", size=0.3, linetype=2))
ggsave("summary_env_NO3_c.pdf")

#Plot env. vs. NH3 (log10), continuously

htdj <- melt(pc_bs, id=c("season","habitat","sal_group","cluster","fjord","log_NH3"), measure=c("log_Si","log_Chlorophyll","Temperature","log_NO2","log_NO3","log_PO4","log_DN",
"Salinity","Longitude","cube_d14N_15N", "d12C_13C", "log_N", "log_C", "log_CN_ratio","log_NP_ratio", "Dry_matter","Water_content", "log_Organic_content", "log_Inorganic_content", "log_TP"))

htdj$id<-paste(htdj$season, htdj$habitat, sep="_")
htdj2 <- ddply(htdj, .(variable, sal_group), dplyr::mutate,, grp.mean=mean(na.omit(value)),grp.sd=sd(na.omit(value)),n_clusters=length(unique(cluster)),n_sites=length(unique(id)))

htdj$sal_group<-as.character(htdj$sal_group)
htdj$season <- factor(htdj$season, levels = c("spring","autumn"))
htdj$habitat <- factor(htdj$habitat, levels = c("eelgrass","sand","rocks"))
htdj$gp <-paste(as.character(htdj$season),as.character(htdj$fjord),sep="_")

setwd("/environmental_data/results/both")
ggplot(htdj, aes(log_NH3, value, color=season)) +
geom_point(data=htdj, aes(log_NH3, value, shape=fjord, color=season), size=0.8) +
facet_wrap(variable~season, scales="free", ncol=6) +
scale_color_brewer(palette="Dark2") +
theme_bw() +
scale_shape_manual(values=c(1, 4)) +
theme(strip.text = element_text(size=3.5),
legend.title=element_text(size=4),legend.text=element_text(size=4), axis.title=element_text(size=4), plot.title = element_text(size=8), axis.text.x = element_text(size=4), axis.text.y = element_text(size=4), panel.grid.major = element_line(colour = "gray80", size=0.3, linetype=2))
ggsave("summary_env_NH3_c.pdf")

#Plot env. vs. NO2 (log10), continuously

htdj <- melt(pc_bs, id=c("season","habitat","sal_group","cluster","fjord","log_NO2"), measure=c("log_Si","log_Chlorophyll","Temperature","log_PO4","log_NO3","log_NH3","log_DN",
"Salinity","Longitude","cube_d14N_15N", "d12C_13C", "log_N", "log_C", "log_CN_ratio","log_NP_ratio", "Dry_matter","Water_content", "log_Organic_content", "log_Inorganic_content", "log_TP"))

htdj$id<-paste(htdj$season, htdj$habitat, sep="_")
htdj2 <- ddply(htdj, .(variable, sal_group), dplyr::mutate,, grp.mean=mean(na.omit(value)),grp.sd=sd(na.omit(value)),
n_clusters=length(unique(cluster)),n_sites=length(unique(id)))

htdj$sal_group<-as.character(htdj$sal_group)
htdj$season <- factor(htdj$season, levels = c("spring","autumn"))
htdj$habitat <- factor(htdj$habitat, levels = c("eelgrass","sand","rocks"))
htdj$gp <-paste(as.character(htdj$season),as.character(htdj$fjord),sep="_")

setwd("/environmental_data/results/both")
ggplot(htdj, aes(log_NO2, value, color=season)) +
geom_point(data=htdj, aes(log_NO2, value, shape=fjord, color=season), size=0.8) +
facet_wrap(variable~season, scales="free", ncol=6) +
scale_color_brewer(palette="Dark2") +
theme_bw() +
scale_shape_manual(values=c(1, 4)) +
theme(strip.text = element_text(size=3.5),
legend.title=element_text(size=4),legend.text=element_text(size=4), axis.title=element_text(size=4), plot.title = element_text(size=8), axis.text.x = element_text(size=4), axis.text.y = element_text(size=4), panel.grid.major = element_line(colour = "gray80", size=0.3, linetype=2))
ggsave("summary_env_NO2_c.pdf")

#Plot env. vs. Time, continuously

htdj <- melt(pc_bs, id=c("season","habitat","sal_group","cluster","fjord","Time"), measure=c("log_PO4","log_Si","log_Chlorophyll","Salinity","Temperature",
"log_NO2","log_NO3","log_NH3","log_DN","Latitude","Longitude","cube_d14N_15N", "d12C_13C", "log_N", "log_C", "log_CN_ratio","log_NP_ratio", "Dry_matter","Water_content", "log_Organic_content", "log_Inorganic_content", "log_TP"))

htdj$id<-paste(htdj$season, htdj$habitat, sep="_")
htdj2 <- ddply(htdj, .(variable, sal_group), dplyr::mutate,, grp.mean=mean(na.omit(value)),grp.sd=sd(na.omit(value)),n_clusters=length(unique(cluster)),n_sites=length(unique(id)))

htdj$sal_group<-as.character(htdj$sal_group)
htdj$season <- factor(htdj$season, levels = c("spring","autumn"))
htdj$habitat <- factor(htdj$habitat, levels = c("eelgrass","sand","rocks"))
htdj$gp <-paste(as.character(htdj$season),as.character(htdj$fjord),sep="_")
htdj$Time<-as.Date(htdj$Time, format= "%d-%m-%y")

setwd("/environmental_data/results/both")
ggplot(htdj, aes(Time, value, color=season)) +
geom_point(data=htdj, aes(Time, value, shape=fjord, color=season), size=0.8) +
facet_wrap(variable~season, scales="free", ncol=6) +
scale_color_brewer(palette="Dark2") +
scale_x_date(date_breaks = "15 days", date_labels = "%d \n %m") +
theme_bw() +
scale_shape_manual(values=c(1, 4)) +
theme(strip.text = element_text(size=3.5),
legend.title=element_text(size=4),legend.text=element_text(size=4), axis.title=element_text(size=4), plot.title = element_text(size=8), axis.text.x = element_text(size=4), axis.text.y = element_text(size=4), panel.grid.major = element_line(colour = "gray80", size=0.3, linetype=2))
ggsave("summary_env_time_c.pdf")



#Plot env. vs. salinity, discreetly

htdj <- melt(pc_bs, id=c("season","habitat","sal_group","cluster","fjord"), measure=c("log_PO4","log_Si","log_Chlorophyll","Salinity","Temperature","log_NO2","log_NO3",
"log_NH3","log_DN","Latitude","Longitude","cube_d14N_15N", "d12C_13C", "log_N", "log_C", "log_CN_ratio","log_NP_ratio", "Dry_matter", "Grain_size", "Water_content", "log_Organic_content", "log_Inorganic_content", "log_TP"))

htdj$id<-paste(htdj$season, htdj$habitat, sep="_")
htdj2 <- ddply(htdj, .(variable, sal_group), dplyr::mutate,, grp.mean=mean(na.omit(value)),grp.sd=sd(na.omit(value)),n_clusters=length(unique(cluster)),n_sites=length(unique(id)))

htdj$sal_group<-as.character(htdj$sal_group)
htdj$season <- factor(htdj$season, levels = c("spring","autumn"))
htdj$habitat <- factor(htdj$habitat, levels = c("eelgrass","sand","rocks"))
htdj$gp <-paste(as.character(htdj$season),as.character(htdj$fjord),sep="_")




















