#!/usr/bin/env bash
set -euo pipefail

################# GLOBAL VARIABLES ####################################

OUTDIR=$1
BAMDIR="${2}/align"
COUNTSDIR="${OUTDIR}/counts"
FEATURECOUNTS=featureCounts
FRAC_COUNTS=$3  #'TRUE' or 'FALSE'
STRANDEDNESS=$4   #can be 0, 1, or 2 (unstranded, forward, reverse)
GTFFILE=$5 #"ref/ncbi_dataset/data/GCF_000001405.40/genomic.gtf"
THREADS=$6
######################################################################


# Create directories and ensure files exist
mkdir -p "${COUNTSDIR}"
echo "Created Counts Directory: ${COUNTSDIR}"

echo "BAM FILES:"
ls ${BAMDIR}/SRR*.Aligned.sortedByCoord.out.bam


[[ -f "${GTFFILE}" ]] || { echo "Missing GTFFILE: ${GTFFILE}"; exit 1; }
[[ -d "${BAMDIR}"  ]] || { echo "Error BAM file Directory: ${BAMDIR} not found"; exit 1; }


## RUN FEATURE COUNTS
echo "[FEATURECOUNTS] Starting feature counts..."

## THis feature count is for paired end reads and unstranded DNA, only counts reads that have both paired ends aligned
## It aslo counts multimapping reads (fractionally) which is beneficial for the hexaploid wheat genome. (will be rounded for DESEQ2 analysis)
if [ ${FRAC_COUNTS} == "TRUE" ]; then
    echo "Counting multi-mapping reads as fractional amounts"
    ${FEATURECOUNTS} \
        -T ${THREADS} \
        -a ${GTFFILE} \
        -o "${COUNTSDIR}/read_counts.txt" \
        -p \
        --countReadPairs \
        -B \
        -s $STRANDEDNESS \
        -M \
        --fraction \
        --verbose \
        ${BAMDIR}/SRR*.Aligned.sortedByCoord.out.bam
elif [ ${FRAC_COUNTS} == "FALSE" ]; then
    echo "Discarding multimapping reads"
    ${FEATURECOUNTS} \
        -T ${THREADS} \
        -a ${GTFFILE} \
        -o "${COUNTSDIR}/read_counts.txt" \
        -p \
        --countReadPairs \
        -B \
        -s $STRANDEDNESS \
        --verbose \
        ${BAMDIR}/SRR*.Aligned.sortedByCoord.out.bam
else
    echo "FRAC_COUNTS must be TRUE or FALSE"
    exit 1
fi

echo "[FEATURECOUNTS] Done"

