# Input directory containing raw FASTQ files
INDIR=$1

# Output directory where trimmed, demultiplexed FASTQs will be written
OUTDIR=$2

# Sample tag name (used to name output files)
TAG_NAME=$3

# Forward tag sequence
TAG_SEQ=$4

# Reverse tag sequence
RTAG_SEQ=$5

# Batch file with primer sequences and the minimum read length required after trimming of primers and tags
BATCHFILE=$6


# Create output folders for sense (SS) and antisense (AS) reads
mkdir -p $OUTDIR/DADA2_AS
mkdir -p $OUTDIR/DADA2_SS


# Read each line of the batch file:
while read INPUT_R1 INPUT_R2 PRIMER_F PRIMER_R MIN_LENGTH ; do

    # Define cutadapt commands:
    # CUTADAPT: strict trimming, requires full primer match, discards untrimmed reads
    CUTADAPT="$(which cutadapt) --discard-untrimmed --minimum-length ${MIN_LENGTH} -e 0"

    # CUTADAPT2: second pass, trims remaining primer sequence without discarding reads
    CUTADAPT2="$(which cutadapt) -e 0"

    # vsearch binary (not used - from earlier version of the script?)
    VSEARCH=$(which vsearch)

    # Temporary file for intermediate FASTQ storage
    TMP_FASTQ=$(mktemp)

    # Minimum primer lengths
    MIN_F=$((${#PRIMER_F}))
    MIN_R=$((${#PRIMER_R}))

    # Compute reverse‑complement primers using bash + tr
    REV_PRIMER_F="$(echo $PRIMER_F | rev | tr ATUGCYRKMBDHVN TAACGRYMKVHDBN)"
    REV_PRIMER_R="$(echo $PRIMER_R | rev | tr ATUGCYRKMBDHVN TAACGRYMKVHDBN)"

    # (Unused variable — likely leftover from earlier version)
    rev="$(echo $primer | rev | tr ATUGCYRKMBDHVN TAACGRYMKVHDBN)"

    # Output file for sense‑strand R1 reads
    FINAL_FASTQ="$OUTDIR/DADA2_SS/${TAG_NAME}_R1.fastq"

    # Build full primer+tag sequences
    FTFP="$TAG_SEQ$PRIMER_F"   # Forward tag + forward primer
    RTRP="$RTAG_SEQ$PRIMER_R"  # Reverse tag + reverse primer


    ############################################################
    # 1. PROCESS SENSE‑STRAND READS (R1)
    # Trim tag+primer from 5' end, then trim reverse‑complement
    # of the reverse primer from the 3' end.
    ############################################################
    cat "$INDIR/${INPUT_R1}" | \
        ${CUTADAPT} -g "${FTFP}" -e 0 -O "${#FTFP}" - | \
        ${CUTADAPT2} -a "${REV_PRIMER_R}" - \
        > "${FINAL_FASTQ}"


    ############################################################
    # 2. PROCESS ANTISENSE‑STRAND READS (R2)
    # Same trimming logic but applied to R2.
    ############################################################
    FINAL_FASTQ="$OUTDIR/DADA2_AS/${TAG_NAME}_R2.fastq"

    cat "$INDIR/${INPUT_R2}" | \
        ${CUTADAPT} -g "${FTFP}" -e 0 -O "${#FTFP}" - | \
        ${CUTADAPT2} -a "${REV_PRIMER_R}" - \
        > "${FINAL_FASTQ}"


    ############################################################
    # 3. PROCESS SENSE‑STRAND R2 USING REVERSE TAG+PRIMER
    # This trims the reverse tag+primer from the 5' end of R2,
    # then trims the reverse‑complement of the forward primer.
    ############################################################
    FINAL_FASTQ="$OUTDIR/DADA2_SS/${TAG_NAME}_R2.fastq"

    cat "$INDIR/${INPUT_R2}" | \
        ${CUTADAPT} -g "${RTRP}" -e 0 -O "${#RTRP}" - | \
        ${CUTADAPT2} -a "${REV_PRIMER_F}" - \
        > "${FINAL_FASTQ}"


    ############################################################
    # 4. PROCESS ANTISENSE‑STRAND R1 USING REVERSE TAG+PRIMER
    ############################################################
    FINAL_FASTQ="$OUTDIR/DADA2_AS/${TAG_NAME}_R1.fastq"

    cat "$INDIR/${INPUT_R1}" | \
        ${CUTADAPT} -g "${RTRP}" -e 0 -O "${#RTRP}" - | \
        ${CUTADAPT2} -a "${REV_PRIMER_F}" - \
        > "${FINAL_FASTQ}"

done < $INDIR/$BATCHFILE
