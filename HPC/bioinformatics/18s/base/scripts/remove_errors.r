# Retrieve command‑line arguments (expects the main path as the first argument)
args = commandArgs(trailingOnly=TRUE)

library(dada2)

# Base directory containing DADA2 input/output folders
main_path <- args[1]

start.time <- Sys.time()

############################################################
# PROCESSING SENSE‑STRAND (SS) READS — forward primer in R1
############################################################

# Path to filtered and matched (R1 to R2) reads
filt_path <- file.path(main_path, "DADA2_SS/filtered/matched")

# List all files in the directory
fns <- list.files(filt_path)

# Keep only FASTQ files
fastqs <- fns[grepl(".fastq.gz", fns)]
fastqs <- sort(fastqs)

# Separate forward and reverse reads based on filename patterns
fnFs <- fastqs[grepl("_F_", fastqs)]
fnRs <- fastqs[grepl("_R_", fastqs)]

# Build full paths to the FASTQ files
filtFs <- file.path(filt_path, fnFs)
filtRs <- file.path(filt_path, fnRs)

# Remove empty or corrupted files (size > 1 byte)
filtFs <- filtFs[sapply(filtFs, file.size) > 1]
filtRs <- filtRs[sapply(filtRs, file.size) > 1]

# Extract sample names from filenames (everything before "_F_")
SSsample.names <- sapply(strsplit(basename(filtFs), "_F_"), `[`, 1)

# Learn error models for forward and reverse reads
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Dereplicate reads (collapse identical sequences)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Assign sample names to dereplicated objects
names(derepFs) <- SSsample.names
names(derepRs) <- SSsample.names

# Run the DADA algorithm to infer true biological sequences
SSdadaFs <- dada(derepFs, err=errF, multithread = TRUE)
SSdadaRs <- dada(derepRs, err=errR, multithread = TRUE)

# Merge paired reads into full amplicons
SSmergers <- mergePairs(SSdadaFs, derepFs, SSdadaRs, derepRs,
                        verbose=TRUE, minOverlap = 5)

# Construct the sequence table (samples × ASVs)
seqtab_SS <- makeSequenceTable(SSmergers[names(SSmergers)])

# Remove chimeric sequences
seqtab.nochim_SS <- removeBimeraDenovo(seqtab_SS, verbose=TRUE)

# Save results
stSS <- file.path(args[1], "seqtab_SS_RDS")
stnsSS <- file.path(args[1], "seqtab.nochim_SS_RDS")
saveRDS(seqtab_SS, stSS)
saveRDS(seqtab.nochim_SS, stnsSS)


############################################################
# PROCESSING ANTISENSE‑STRAND (AS) READS — reverse orientation
############################################################

# Path to filtered antisense reads
filt_path <- file.path(main_path, "DADA2_AS/filtered/matched")

# List FASTQ files
fns <- list.files(filt_path)
fastqs <- fns[grepl(".fastq.gz", fns)]
fastqs <- sort(fastqs)

# Separate forward and reverse reads
fnFs <- fastqs[grepl("_F_", fastqs)]
fnRs <- fastqs[grepl("_R_", fastqs)]

# Build full paths
filtFs <- file.path(filt_path, fnFs)
filtRs <- file.path(filt_path, fnRs)

# Remove empty files
filtFs <- filtFs[sapply(filtFs, file.size) > 1]
filtRs <- filtRs[sapply(filtRs, file.size) > 1]

# Extract sample names
ASsample.names <- sapply(strsplit(basename(filtFs), "_F_"), `[`, 1)

# Learn error models
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Dereplicate reads
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Assign sample names
names(derepFs) <- ASsample.names
names(derepRs) <- ASsample.names

# Infer ASVs
ASdadaFs <- dada(derepFs, err=errF, multithread = TRUE)
ASdadaRs <- dada(derepRs, err=errR, multithread = TRUE)

# Merge paired reads
ASmergers <- mergePairs(ASdadaFs, derepFs, ASdadaRs, derepRs,
                        verbose=TRUE, minOverlap = 5)

# Build sequence table
seqtab_AS <- makeSequenceTable(ASmergers[names(ASmergers)])

# Remove chimeras
seqtab.nochim_AS <- removeBimeraDenovo(seqtab_AS, verbose=TRUE)

# Save results
stAS <- file.path(args[1], "seqtab_AS_RDS")
stnsAS <- file.path(args[1], "seqtab.nochim_AS_RDS")
saveRDS(seqtab_AS, stAS)
saveRDS(seqtab.nochim_AS, stnsAS)
