#!/bin/bash
################################################################################
# METAGENOMICS WORKFLOW - SAMPLE DISCOVERY AND SETUP MODULE (01_sample_discovery.sh)
################################################################################
#
# Purpose: Identifies all samples in the working directory based on file
#          extension pattern, extracts sample names, and sets up the job array.
#
# Usage: source ./01_sample_discovery.sh (called by main orchestrator)
#
# Outputs: SAMPLE_LIST array containing unique sample names, N = number of samples
#
# Dependencies: Requires DIRECTORY and EXTENSION from 00_config.sh
#
# Supported filename conventions (separator = '.' or '_'):
#   <sample>.R1.fastq.gz        <sample>.R2.fastq.gz
#   <sample>_R1.fastq.gz        <sample>_R2.fastq.gz
#   <sample>.1.fastq.gz         <sample>.2.fastq.gz
#   <sample>_1.fastq.gz         <sample>_2.fastq.gz
#   <sample>_R1.fq.gz           <sample>_R2.fq.gz
#   <sample>-filtered.R1.fastq.gz  (the '-filtered' is part of the sample name)
#   All of the above without the trailing .gz
#
################################################################################

# Verify required variables are set
if [ -z "$DIRECTORY" ] || [ -z "$EXTENSION" ]; then
    echo "Error: DIRECTORY and/or EXTENSION not set. Run 00_config.sh first."
    exit 1
fi

################################################################################
# SAMPLE DISCOVERY
################################################################################

SEARCH_DIR="$DIRECTORY"
echo "$(date): Searching for files in $SEARCH_DIR with extension $EXTENSION..."

# Strip the read-pair suffix to recover the sample name.
#
# Pattern breakdown:
#   [._]       — separator before the read identifier (dot or underscore)
#   (R[12]|[12]) — read identifier: R1, R2, 1, or 2
#   (\.filtered)? — optional literal '.filtered' tag
#   (\.fastq|\.fq) — mandatory read-file extension
#   (\.gz)?    — optional gzip suffix
#   $          — anchored to end of string
#
# Because the pattern is anchored at the end, internal dots/underscores that
# are part of sample names (e.g. '2017.12.29_11.37.26') are left untouched.

SAMPLE_LIST=($(find "$SEARCH_DIR" -maxdepth 1 -type f -name "*$EXTENSION" -print0 | \
    xargs -0 -I {} basename {} | \
    sed -E 's/[._](R[12]|[12])(\.filtered)?(\.fastq|\.fq)(\.gz)?$//' | \
    grep -v '\.sh$' | \
    sort -u))

# Verify samples were found
if [ ${#SAMPLE_LIST[@]} -eq 0 ]; then
    echo "$(date): No samples found matching extension '$EXTENSION' in $SEARCH_DIR"
    exit 1
fi

N=${#SAMPLE_LIST[@]}

# Validate sample count
if [ -z "$N" ] || [ "$N" -lt 1 ]; then
    echo "$(date): Error: Number of samples (N) is invalid: $N"
    exit 1
fi

echo "$(date): Found $N samples: ${SAMPLE_LIST[@]}"
echo "$(date): Queueing $N jobs and setting $CONCURRENT concurrent jobs"
echo "$(date): Working directory: $DIRECTORY"

# Export variables for use in subsequent scripts
export SAMPLE_LIST
export N
export SEARCH_DIR

# Write sample names to flat file for SLURM array indexing
printf '%s\n' "${SAMPLE_LIST[@]}" > "${DIRECTORY}/samples.txt"
echo "$(date): Sample list written to ${DIRECTORY}/samples.txt"
echo "$(date): Sample discovery complete."