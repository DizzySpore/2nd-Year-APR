#!/bin/bash

################################################################################
# METAGENOMICS WORKFLOW - MAIN ORCHESTRATOR SCRIPT
################################################################################
#
# Purpose: Central entry point for the entire metagenomics workflow.
#          Orchestrates the execution of all workflow modules in proper sequence:
#          1. Gathers configuration from user
#          2. Discovers samples
#          3. Scans checkpoints and determines resume points
#          4. Generates SLURM template
#          5. Submits job array to cluster
#          6. Monitors job execution
#
# Usage: ./main_workflow.sh
#
# Structure:
#   00_config.sh              - Interactive configuration and setup
#   01_sample_discovery.sh    - Find and validate samples
#   03_checkpoint.sh          - Scan outputs and determine resume points
#   02_slurm_template.sh      - Generate SLURM batch script
#   [SLURM execution]         - Batch processing on cluster
#   [Log monitoring]          - Track job progress
#
# Binning pipeline (replaces metaWRAP):
#   stage_07a_bowtie2_mapping.sh  - Read mapping to assembly
#   stage_07b_binning.sh          - Multi-binner binning
#   stage_07c_checkm.sh           - Bin quality assessment
#   stage_07d_bin_refinement.sh   - Hybrid bin creation
#   stage_07e_reassembly.sh       - SPAdes reassembly
#   stage_07f_bin_selection.sh    - Best version selection
#
################################################################################

set -e

################################################################################
# BANNER AND INITIALIZATION
################################################################################

clear
echo "+------------------------------------------------------------------------+"
echo "|         METAGENOMICS WORKFLOW - MAIN ORCHESTRATOR (v1.0.0)            |"
echo "|                                                                        |"
echo "|  This script orchestrates a complete metagenomic analysis workflow:   |"
echo "|  - Quality Control (FastQC, fastp, KneadData)                         |"
echo "|  - Assembly (MEGAHIT, QUAST)                                          |"
echo "|  - Binning (Bowtie2, CONCOCT, MetaBAT2, COMEbin, SemiBin2,           |"
echo "|             CheckM, CheckM2, bin refinement, SPAdes reassembly)       |"
echo "|  - Taxonomy (GTDB-Tk, Kraken2, MetaPhlAn, Kaiju, EukFinder)           |"
echo "|  - Annotation (Bakta, antiSMASH)                                      |"
echo "|  - Functional Profiling (HUMAnN3)                                     |"
echo "|  - Statistical Analysis (MaAsLin3, melonnpan)                         |"
echo "+------------------------------------------------------------------------+"
echo ""

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export SCRIPT_DIR

# Verify all required module scripts exist
echo "$(date): Verifying required module scripts..."
REQUIRED_SCRIPTS=(
    "00_config.sh"
    "01_sample_discovery.sh"
    "03_checkpoint.sh"
    "02_slurm_template.sh"
)

for script in "${REQUIRED_SCRIPTS[@]}"; do
    if [ ! -f "$SCRIPT_DIR/$script" ]; then
        echo "ERROR: Required script not found: $SCRIPT_DIR/$script"
        exit 1
    fi
done

echo "$(date): All required modules verified."
echo ""

################################################################################
# STAGE 1: CONFIGURATION
################################################################################

echo "+------------------------------------------------------------------------+"
echo "| STAGE 1: INTERACTIVE CONFIGURATION                                    |"
echo "+------------------------------------------------------------------------+"
echo ""

source "$SCRIPT_DIR/00_config.sh"

if [ -z "$DIRECTORY" ] || [ -z "$EMAIL" ]; then
    echo "ERROR: Configuration failed. Exiting."
    exit 1
fi

echo ""
echo "? Configuration completed successfully"
echo ""

################################################################################
# STAGE 2: SAMPLE DISCOVERY
################################################################################

echo "+------------------------------------------------------------------------+"
echo "| STAGE 2: SAMPLE DISCOVERY AND VALIDATION                              |"
echo "+------------------------------------------------------------------------+"
echo ""

source "$SCRIPT_DIR/01_sample_discovery.sh"

if [ -z "$N" ] || [ ${#SAMPLE_LIST[@]} -eq 0 ]; then
    echo "ERROR: Sample discovery failed. No samples found. Exiting."
    exit 1
fi

echo ""
echo "? Sample discovery completed successfully"
echo "? Found $N samples to process"
echo ""

################################################################################
# STAGE 3: CHECKPOINT SCAN
################################################################################

echo "+------------------------------------------------------------------------+"
echo "| STAGE 3: CHECKPOINT SCAN AND RESUME POINT DETECTION                  |"
echo "+------------------------------------------------------------------------+"
echo ""

source "$SCRIPT_DIR/03_checkpoint.sh"

if [ "$N" -eq 0 ]; then
    echo ""
    echo "? All samples are already complete. Nothing to submit."
    echo ""
    exit 0
fi

echo ""
echo "? Checkpoint scan completed successfully"
echo "? $N samples pending submission"
echo ""

################################################################################
# STAGE 4: SLURM JOB GENERATION
################################################################################

echo "+------------------------------------------------------------------------+"
echo "| STAGE 4: SLURM JOB TEMPLATE GENERATION                               |"
echo "+------------------------------------------------------------------------+"
echo ""

source "$SCRIPT_DIR/02_slurm_template.sh"

if [ ! -f "$DIRECTORY/metagenomics_workflow.sh" ]; then
    echo "ERROR: SLURM template generation failed. Exiting."
    exit 1
fi

echo ""
echo "? SLURM job template generated successfully"
echo ""

################################################################################
# STAGE 5: JOB SUBMISSION
################################################################################

echo "+------------------------------------------------------------------------+"
echo "| STAGE 5: SLURM JOB SUBMISSION                                        |"
echo "+------------------------------------------------------------------------+"
echo ""

cd "$DIRECTORY"

echo "$(date): Submitting job array to SLURM cluster..."
echo "  - Array size: $N samples"
echo "  - Concurrent jobs: $CONCURRENT"
echo "  - Job name: metagenomics_workflow"
echo ""

SBATCH_OUTPUT=$(sbatch metagenomics_workflow.sh 2>&1)
JOB_ID=$(echo "$SBATCH_OUTPUT" | grep -oP 'Submitted batch job \K[0-9]+' || true)

if [ -z "$JOB_ID" ]; then
    echo "ERROR: Failed to submit SLURM job"
    echo "sbatch output: $SBATCH_OUTPUT"
    exit 1
fi

echo "? Job submitted successfully"
echo "  Job ID: $JOB_ID"
echo ""

################################################################################
# STAGE 6: JOB MONITORING
################################################################################

echo "+------------------------------------------------------------------------+"
echo "| STAGE 6: JOB MONITORING AND LOG TRACKING                             |"
echo "+------------------------------------------------------------------------+"
echo ""

LOG_FILE="${DIRECTORY}/workflow-execution.log"
echo "$(date): Job ID: $JOB_ID" | tee -a "$LOG_FILE"
echo ""

# Wait for the job to start
echo "$(date): Waiting for job to start..."
JOB_STARTED=false

for i in {1..60}; do
    JOB_STATE=$(squeue --job ${JOB_ID}_1 --noheader --format=%T 2>/dev/null || true)

    if [[ "$JOB_STATE" == "RUNNING" ]]; then
        echo "? Job ${JOB_ID}_1 started at $(date)"
        JOB_STARTED=true
        break
    elif [[ -z "$JOB_STATE" ]]; then
        echo "$(date): Job $JOB_ID not yet queued..."
    else
        echo "$(date): Job $JOB_ID is in state: $JOB_STATE"
    fi

    sleep 5
done

if [ "$JOB_STARTED" = false ]; then
    echo "WARNING: Job did not start within 5 minutes. Check SLURM status manually:"
    echo "  squeue --job $JOB_ID"
    echo "  sinfo -N"
fi

echo ""
echo "+------------------------------------------------------------------------+"
echo "| WORKFLOW SUBMITTED AND MONITORING IN PROGRESS                        |"
echo "+------------------------------------------------------------------------+"
echo ""
echo "  squeue -j $JOB_ID              - Check job status"
echo "  scancel $JOB_ID                - Cancel the job"
echo "  tail -f slurm-${JOB_ID}_*.out  - View live output"
echo "  tail -f $LOG_FILE              - View workflow log"
echo ""

# Monitor job completion
echo "$(date): Monitoring job array progress..."
LAST_STATUS=""

while true; do
    ACTIVE_JOBS=$(squeue -j $JOB_ID --noheader 2>/dev/null | wc -l)

    if [ "$ACTIVE_JOBS" -eq 0 ]; then
        echo ""
        echo "? Job array $JOB_ID has completed."
        break
    else
        CURRENT_TIME=$(date "+%Y-%m-%d %H:%M:%S")
        STATUS="[$CURRENT_TIME] Active jobs: $ACTIVE_JOBS"

        if [ "$STATUS" != "$LAST_STATUS" ]; then
            echo "$STATUS"
            LAST_STATUS="$STATUS"
        fi
    fi

    sleep 30
done

################################################################################
# COMPLETION SUMMARY
################################################################################

echo ""
echo "+------------------------------------------------------------------------+"
echo "| WORKFLOW EXECUTION COMPLETE                                           |"
echo "+------------------------------------------------------------------------+"
echo ""
echo "  End time:      $(date)"
echo "  Job ID:        $JOB_ID"
echo "  Total samples: $N"
echo "  Results:       $DIRECTORY"
echo ""

if [ -d "${DIRECTORY}/status" ]; then
    COMPLETED=$(ls "${DIRECTORY}/status"/*.done 2>/dev/null | wc -l)
    FAILED=$(ls "${DIRECTORY}/status"/*.failed 2>/dev/null | wc -l)

    echo "  Completed: $COMPLETED/$N"
    echo "  Failed:    $FAILED/$N"
    echo ""

    if [ "$FAILED" -gt 0 ]; then
        echo "Failed samples:"
        ls -1 "${DIRECTORY}/status"/*.failed 2>/dev/null | \
            xargs -I {} basename {} .failed || true
        echo ""
    fi
fi

echo "Thank you for using the Metagenomics Workflow!"
echo ""