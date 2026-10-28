#!/bin/bash

################################################################################
# METAGENOMICS WORKFLOW - CHECKPOINT SCRIPT (03_checkpoint.sh)
################################################################################
#
# Purpose: Pre-submission checkpoint scanner. Examines output files for each
#          sample to determine the last completed stage, then writes a .resume
#          file so the SLURM job knows where to start. Samples that are already
#          fully complete are excluded from the submission array.
#
# Usage: source ./03_checkpoint.sh (called by main orchestrator)
#
# Dependencies: Requires DIRECTORY, SAMPLE_LIST, N from previous modules
#
# Outputs:
#   - ${DIRECTORY}/status/${SAMPLE}.resume  - Stage to resume from
#   - ${DIRECTORY}/samples.txt              - Filtered sample list for submission
#   - Updates N to reflect only pending samples
#
################################################################################

# Verify required variables
if [ -z "$DIRECTORY" ] || [ -z "$N" ] || [ ${#SAMPLE_LIST[@]} -eq 0 ]; then
    echo "Error: Required variables not set. Run previous modules first."
    return 1
fi

STATUS_DIR="${DIRECTORY}/status"
mkdir -p "$STATUS_DIR"

echo "$(date): Beginning checkpoint scan for ${N} samples..."

################################################################################
# CHECKPOINT DEFINITIONS
################################################################################
# Format: "STAGE_XX:relative/path/to/expected/output"
# SAMPLENAME is substituted with the actual sample name inside the loop.
# Wildcard paths are supported for stages with variable output filenames.

declare -a CHECKPOINTS=(
    "STAGE_01:qc_reads/fastqc_raw/SAMPLENAME/SAMPLENAME_R1_001_fastqc.html"
    "STAGE_02:qc_reads/SAMPLENAME_fastp.html"
    "STAGE_03:dehuman_reads/kneaddata_output_SAMPLENAME/out.SAMPLENAME.R1_kneaddata_paired_1.fastq"
    "STAGE_04:qc_reads/fastqc_clean/SAMPLENAME/out.SAMPLENAME.R1_kneaddata_paired_1_fastqc.html"
    "STAGE_05:megahit_assembly/SAMPLENAMEassembly/final.contigs.fa"
    "STAGE_06:quast_assembly/SAMPLENAME/SAMPLENAME_quast_results/report.html"
    "STAGE_07A:mapping/SAMPLENAME/SAMPLENAME.sorted.bam"
    "STAGE_07B:binning/SAMPLENAME/metabat2_bins/*.fa"
    "STAGE_07C:bin_refinement/SAMPLENAME/metawrap_70_10_bins/*.fa"
    "STAGE_07D:bin_reassembly/SAMPLENAME/reassembled_bins/*.orig.fa"
    "STAGE_07E:final_bins/SAMPLENAME/reassembled_bins.stats"
    "STAGE_08:gtdbtk-out/gtdbtk-SAMPLENAME-out/gtdbtk.bac120.summary.tsv"
    "STAGE_09:quast_bins/SAMPLENAME/SAMPLENAME_quast_results/report.tsv"
    "STAGE_10:kraken2_bracken_out/SAMPLENAME_prok-bracken_abundance.txt"
    "STAGE_11:status/medi_complete.txt"
    "STAGE_12:metaphlan_out/SAMPLENAME_metaphlan_profile.txt"
    "STAGE_13:bakta_out/SAMPLENAME/loop_completed.txt"
    "STAGE_13B:bakta_quast_out/SAMPLENAME/quast_loop_completed.txt"
    "STAGE_13C:dbcan_out/SAMPLENAME/loop_completed.txt"
    "STAGE_13D:eggnog_out/SAMPLENAME/loop_completed.txt"
    "STAGE_13E:abricate_out/SAMPLENAME/loop_completed.txt"
    "STAGE_14:good_bins/SAMPLENAME/*.fa"
    "STAGE_15:eukfinder_out/SAMPLENAME/gunzipped_reads/Eukfinder_results/scf_SAMPLENAME_eukfinder_short_seq_.Unk.fasta"
    "STAGE_16:kraken2_bracken_out/SAMPLENAME_euk-fin_bracken_abundance.txt"
    "STAGE_17:status/panaroo_complete.txt"
    "STAGE_18:status/snippy_complete.txt"
    "STAGE_18B:status/raxml_complete.txt"
    "STAGE_19:antismash_out/SAMPLENAME/antismash_loop_completed.txt"
    "STAGE_20:humann_reads_output/SAMPLENAME/out.SAMPLENAME_combined_pathabundance.tsv"
    "STAGE_21:humann_contigs_output/SAMPLENAME/humann_loop_completed.txt"
    "STAGE_21B:humann_contigs_output/SAMPLENAME/stage21b_regroup_completed.txt"
    "STAGE_22:status/humann_combine_complete.txt"
    "STAGE_23:status/abundance_visualization_complete.txt"
    "STAGE_24:maaslin3_output_bracken/significant_results.csv"
)

################################################################################
# NEXT STAGE LOOKUP
################################################################################

get_next_stage() {
    completed="$1"
    case "$completed" in
        "STAGE_01")  echo "STAGE_02"  ;;
        "STAGE_02")  echo "STAGE_03"  ;;
        "STAGE_03")  echo "STAGE_04"  ;;
        "STAGE_04")  echo "STAGE_05"  ;;
        "STAGE_05")  echo "STAGE_06"  ;;
        "STAGE_06")  echo "STAGE_07A" ;;
        "STAGE_07A") echo "STAGE_07B" ;;
        "STAGE_07B") echo "STAGE_07C" ;;
        "STAGE_07C") echo "STAGE_07D" ;;
        "STAGE_07D") echo "STAGE_07E" ;;
        "STAGE_07E") echo "STAGE_08"  ;;
        "STAGE_08")  echo "STAGE_09"  ;;
        "STAGE_09")  echo "STAGE_10"  ;;
        "STAGE_10")  echo "STAGE_11"  ;;
        "STAGE_11")  echo "STAGE_12"  ;;
        "STAGE_12")  echo "STAGE_13"  ;;
        "STAGE_13")  echo "STAGE_13B" ;;
        "STAGE_13B") echo "STAGE_13C" ;;
        "STAGE_13C") echo "STAGE_13D" ;;
        "STAGE_13D") echo "STAGE_13E" ;;
        "STAGE_13E") echo "STAGE_14"  ;;
        "STAGE_14")  echo "STAGE_15"  ;;
        "STAGE_15")  echo "STAGE_16"  ;;
        "STAGE_16")  echo "STAGE_17"  ;;
        "STAGE_17")  echo "STAGE_18"  ;;
        "STAGE_18")  echo "STAGE_18B" ;;
        "STAGE_18B") echo "STAGE_19"  ;;
        "STAGE_19")  echo "STAGE_20"  ;;
        "STAGE_20")  echo "STAGE_21"  ;;
        "STAGE_21")  echo "STAGE_21B" ;;
        "STAGE_21B") echo "STAGE_22"  ;;
        "STAGE_22")  echo "STAGE_23"  ;;
        "STAGE_23")  echo "STAGE_24"  ;;
        "STAGE_24")  echo "COMPLETE"  ;;
        *)           echo "STAGE_01"  ;;
    esac
}

################################################################################
# READ SELECTED STAGES
################################################################################

SELECTED_STAGES_LIST=$(cat "${DIRECTORY}/selected_stages.txt")
FIRST_SELECTED=$(head -n1 "${DIRECTORY}/selected_stages.txt")
STAGE_COUNT=$(wc -l < "${DIRECTORY}/selected_stages.txt")
TOTAL_STAGES=${#CHECKPOINTS[@]}
FULL_WORKFLOW=$([ "$STAGE_COUNT" -eq "$TOTAL_STAGES" ] && echo "true" || echo "false")

echo "$(date): Workflow mode: $( [ "$FULL_WORKFLOW" = "true" ] && echo "full" || echo "custom (${STAGE_COUNT} stages)" )"

################################################################################
# PER-SAMPLE CHECKPOINT SCAN
################################################################################

PENDING_SAMPLES=()

for SAMPLE in "${SAMPLE_LIST[@]}"; do
    echo "$(date): Scanning checkpoints for ${SAMPLE}..."

    if [ "$FULL_WORKFLOW" = "true" ]; then
        # -- Full workflow: scan disk to find last completed stage -------------
        LAST_FOUND=""

        for checkpoint in "${CHECKPOINTS[@]}"; do
            step="${checkpoint%%:*}"
            rel_path="${checkpoint#*:}"
            rel_path="${rel_path//SAMPLENAME/${SAMPLE}}"
            full_path="${DIRECTORY}/${rel_path}"

            if [[ "$full_path" == *"*"* ]]; then
                if ls ${full_path} 1>/dev/null 2>&1; then
                    LAST_FOUND="$step"
                fi
            else
                if [ -f "$full_path" ]; then
                    LAST_FOUND="$step"
                fi
            fi
        done

        NEXT_STAGE=$(get_next_stage "$LAST_FOUND")

        if [ "$NEXT_STAGE" = "COMPLETE" ]; then
            echo "$(date): ${SAMPLE} - all stages complete, skipping."
            touch "${STATUS_DIR}/${SAMPLE}.done"
            continue
        fi

        echo "$(date): ${SAMPLE} - will resume from ${NEXT_STAGE}"
        echo "$NEXT_STAGE" > "${STATUS_DIR}/${SAMPLE}.resume"
        PENDING_SAMPLES+=("$SAMPLE")

    else
        # -- Custom workflow: always rerun from first selected stage -----------
        # Ignore existing outputs — user explicitly requested these stages
        echo "$(date): ${SAMPLE} - custom workflow, will run from ${FIRST_SELECTED}"
        echo "$FIRST_SELECTED" > "${STATUS_DIR}/${SAMPLE}.resume"
        PENDING_SAMPLES+=("$SAMPLE")
    fi

done

################################################################################
# WRITE PENDING SAMPLE LIST AND UPDATE N
################################################################################

if [ ${#PENDING_SAMPLES[@]} -eq 0 ]; then
    echo "$(date): All samples are already complete. Nothing to submit."
    return 0
fi

# Write pending samples to flat file for SLURM array indexing
printf '%s\n' "${PENDING_SAMPLES[@]}" > "${DIRECTORY}/samples.txt"

# Update N to reflect only pending samples
N=${#PENDING_SAMPLES[@]}
export N

echo "$(date): ${N} samples pending submission."
echo "$(date): Sample list written to ${DIRECTORY}/samples.txt"
echo "$(date): Checkpoint scan complete."