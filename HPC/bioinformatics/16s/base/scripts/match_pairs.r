# Retrieve command‑line arguments
args = commandArgs(trailingOnly=TRUE)

# Parse comma‑separated input and output file lists from command‑line arguments
inputFiles = unlist(strsplit(args[1], ","))
outputFiles = unlist(strsplit(args[2], ","))

# Print the parsed file lists for logging/debugging
print("inputs")
print(inputFiles)
print("outputs")
print(outputFiles)

library(dada2)

start.time <- Sys.time()

############################################################
# PROCESSING SENSE READS (inputFiles[3] and inputFiles[4])
############################################################

# Check whether either sense FASTQ file is empty (size == 0)
if ( file.info(inputFiles[3])$size == 0 || file.info(inputFiles[4])$size == 0 ) {

  # Report empty files and skip processing
  print(paste(inputFiles[3], " has size ", file.info(inputFiles[3])$size))
  print(paste(inputFiles[4], " has size ", file.info(inputFiles[4])$size))

} else {

  # Log which files are being processed
  print(paste("processing ", inputFiles[3], "and ", inputFiles[4]))

  # Run DADA2 paired‑end filtering:
  # - minLen=160: discard reads shorter than 160 bp
  # - maxN=0: no ambiguous bases allowed
  # - maxEE=2: expected error threshold
  # - truncQ=2: truncate reads at first quality score <= 2
  # - matchIDs=TRUE: ensure R1/R2 read IDs match
  fastqPairedFilter(
    c(inputFiles[3], inputFiles[4]),
    c(outputFiles[3], outputFiles[4]),
    minLen = 160,
    maxN = 0,
    maxEE = 2,
    truncQ = 2,
    matchIDs = TRUE
  )

  print("done")
}


############################################################
# PROCESSING ANTISENSE READS (inputFiles[1] and inputFiles[2])
############################################################

# Check whether either antisense FASTQ file is empty
if ( file.info(inputFiles[1])$size == 0 || file.info(inputFiles[2])$size == 0 ) {

  # Report empty files and skip processing
  print(paste(inputFiles[1], " has size ", file.info(inputFiles[1])$size))
  print(paste(inputFiles[2], " has size ", file.info(inputFiles[2])$size))

} else {

  # Log which files are being processed
  print(paste("processing ", inputFiles[1], "and ", inputFiles[2]))

  # Apply the same paired‑end filtering parameters as above
  fastqPairedFilter(
    c(inputFiles[1], inputFiles[2]),
    c(outputFiles[1], outputFiles[2]),
    minLen = 160,
    maxN = 0,
    maxEE = 2,
    truncQ = 2,
    matchIDs = TRUE
  )

  print("done")
}
