#!/bin/bash

################################################################################
# METAGENOMICS WORKFLOW - CONFIGURATION MODULE (00_config.sh)
################################################################################
# 
# Purpose: Centralized configuration and input handling for the metagenomics 
#          workflow. Gathers user inputs, validates settings, and caches 
#          frequently-used values like email addresses.
#
# Usage: source ./00_config.sh (called by main orchestrator)
#
# Outputs: Sets global variables for directories, SLURM resources, email, etc.
#
################################################################################

# Function to get valid numeric input with bounds checking
# Args: $1 = prompt message, $2 = min value, $3 = max value
# Returns: validated numeric input
get_numeric_input() {
    local prompt="$1"
    local min="$2"
    local max="$3"
    local value
    while true; do
        read -p "$prompt" value
        if [[ "$value" =~ ^[0-9]+$ ]] && [ "$value" -ge "$min" ] && [ "$value" -le "$max" ]; then
            echo "$value"
            return 0
        else
            echo "Please enter a valid number between $min and $max"
        fi
    done
}

# Function to get and validate directory input
# Creates directory if needed, uses current working directory as default
# Args: $1 = prompt message
# Returns: validated directory path
get_directory_input() {
    local prompt="$1"
    local dir
    local return_dir
    while true; do
        read -p "$prompt" dir
        if [ -z "$dir" ]; then
            dir="$(pwd)"
            echo "Using current working directory: $dir" >&2
            return_dir="$dir"
            break
        elif [ -d "$dir" ]; then
            return_dir="$dir"
            break
        else
            read -p "Directory does not exist. Create it? (y/n): " create
            if [[ "$create" =~ ^[Yy]$ ]]; then
                mkdir -p "$dir"
                return_dir="$dir"
                break
            fi
        fi
    done
    echo "$return_dir"
    return 0
}

################################################################################
# INTERACTIVE CONFIGURATION SECTION
################################################################################

echo "=== Metagenomics Workflow Configuration ==="

# ============================================================================
# SECTION 1: WORKING DIRECTORY AND FILE EXTENSION SETUP
# ============================================================================
DIRECTORY=$(get_directory_input "Enter working directory path (or press Enter for current directory): ")

read -p "Enter file extension (e.g., .fastq.gz): " EXTENSION
if [ -z "$EXTENSION" ]; then
    EXTENSION=".fastq.gz"
    echo "Using default extension: $EXTENSION"
fi

# ============================================================================
# SECTION 2: SLURM RESOURCE CONFIGURATION
# ============================================================================
echo -e "\n=== SLURM Resource Configuration ==="

WALL_TIME=$(get_numeric_input "Enter wall time in hours (1-240): " 1 240)
WALL_TIME_PADDED=$(printf "%03d" "$WALL_TIME")

MEM_SIZE=$(get_numeric_input "Enter memory in GB (4-450): " 4 450)
NUM_CORES=$(get_numeric_input "Enter number of CPU cores (1-50): " 1 50)
CONCURRENT=$(get_numeric_input "Enter maximum concurrent jobs (1-50): " 1 50)

# ============================================================================
# SECTION 3: EMAIL CONFIGURATION WITH CACHING
# ============================================================================
DEFAULT_EMAIL=$(git config --global user.email 2>/dev/null) || true

CACHE_FILE="${HOME}/.slurm_email_cache"

get_cached_email() {
    local user=$1
    if [ -f "$CACHE_FILE" ]; then
        grep "^${user}:" "$CACHE_FILE" | cut -d':' -f2 || true
    fi
}

save_email_to_cache() {
    local user=$1
    local email=$2
    if [ -f "$CACHE_FILE" ]; then
        grep -v "^${user}:" "$CACHE_FILE" > "${CACHE_FILE}.tmp" || true
        mv "${CACHE_FILE}.tmp" "$CACHE_FILE"
    fi
    echo "${user}:${email}" >> "$CACHE_FILE"
}

CACHED_EMAIL=$(get_cached_email "$USER")
if [ -n "$CACHED_EMAIL" ]; then
    DEFAULT_EMAIL="$CACHED_EMAIL"
else
    DEFAULT_EMAIL=$(git config --global user.email 2>/dev/null) || true

    if [ -z "$DEFAULT_EMAIL" ]; then
        if [[ "$USER" =~ ^[a-z]{1,3}[0-9]{3}$ ]]; then
            DEFAULT_EMAIL="${USER}@student.bham.ac.uk"
        else
            echo "Warning: Non-student user ID detected. An email cannot be auto-generated."
            read -p "Enter email for SLURM notifications: " EMAIL
            save_email_to_cache "$USER" "$EMAIL"
        fi
    fi
fi

if [ -n "$DEFAULT_EMAIL" ]; then
    read -p "Enter email for SLURM notifications (default: $DEFAULT_EMAIL): " EMAIL
    EMAIL="${EMAIL:-$DEFAULT_EMAIL}"
    save_email_to_cache "$USER" "$EMAIL"
fi

# ============================================================================
# SECTION 4: CONFIGURATION CONFIRMATION
# ============================================================================
echo -e "\n=== Configuration Summary ==="
echo "Working Directory: $DIRECTORY"
echo "File Extension: $EXTENSION"
echo "Wall Time: ${WALL_TIME_PADDED}:00:00"
echo "Memory: ${MEM_SIZE}G"
echo "CPU Cores: $NUM_CORES"
echo "Max Concurrent Jobs: $CONCURRENT"
echo "Email: $EMAIL"

read -p "Proceed with these settings? (y/n): " confirm
if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
    echo "Configuration cancelled. Exiting..."
    return 1
fi

# ============================================================================
# SECTION 4B: WORKFLOW SELECTION
# ============================================================================
echo -e "\n=== Workflow Configuration ==="

# All available stages including segmented binning pipeline
declare -A STAGE_NAMES=(
    ["STAGE_01"]="FastQC - Quality check on raw reads"
    ["STAGE_02"]="fastp - Read trimming and quality control"
    ["STAGE_03"]="KneadData - Human read decontamination"
    ["STAGE_04"]="KneadData inspection - FastQC on decontaminated reads"
    ["STAGE_05"]="MEGAHIT - Metagenomic assembly"
    ["STAGE_06"]="QUAST assembly - Assembly quality assessment"
    ["STAGE_07A"]="Bowtie2 - Read mapping to assembly"
    ["STAGE_07B"]="Binning - CONCOCT, MetaBAT2, COMEbin, SemiBin2"
    ["STAGE_07C"]="Bin refinement - Hybrid creation, quality assessment and selection"
    ["STAGE_07D"]="Reassembly - SPAdes strict and permissive"
    ["STAGE_07E"]="Bin selection - Best version selection and labelling"
    ["STAGE_08"]="GTDB-Tk - Bin taxonomic classification"
    ["STAGE_09"]="QUAST bins - Bin quality assessment"
    ["STAGE_10"]="Kraken2/Bracken - Read-level taxonomic profiling"
    ["STAGE_11"]="MEDI - Dietary/food metagenomic quantification (Nextflow, batch)"
    ["STAGE_12"]="MetaPhlAn - Marker gene taxonomic profiling"
    ["STAGE_13"]="Bakta - Prokaryotic gene annotation (unfiltered bins from STAGE_07E)"
    ["STAGE_13B"]="QUAST Bakta - Annotation quality assessment"
    ["STAGE_13C"]="dbCAN - CAZyme annotation (Bakta proteins)"
    ["STAGE_13D"]="eggNOG-mapper - Functional/orthology annotation (Bakta proteins)"
    ["STAGE_13E"]="abricate - AMR/virulence gene screening (card, vfdb, resfinder, plasmidfinder)"
    ["STAGE_14"]="Bin filtering - Quality filter high-confidence bins"
    ["STAGE_15"]="EukFinder - Eukaryotic sequence detection"
    ["STAGE_16"]="Kraken2 Euk - Eukaryotic taxonomic profiling"
    ["STAGE_17"]="Panaroo - Targeted pangenome analysis (batch, selected taxa)"
    ["STAGE_18"]="Snippy - Within-group SNP calling on Panaroo taxon groups (batch)"
    ["STAGE_18B"]="RAxML-NG - Phylogenetic tree per taxon group from Snippy core-SNP alignment (batch)"
    ["STAGE_19"]="antiSMASH - Secondary metabolite detection"
    ["STAGE_20"]="HUMAnN3 reads - Functional profiling on reads"
    ["STAGE_21"]="HUMAnN3 contigs - Functional profiling on contigs"
    ["STAGE_21B"]="HUMAnN3 regroup - CPM renorm + regroup/rename to rxn/EC/KO (reads + per-bin)"
    ["STAGE_22"]="Combining - Merge HUMAnN3 outputs [Currently out of order- Will be modified]"
    ["STAGE_23"]="Visualization - Kraken2/Bracken visualizations [Currently out of order- Will be modified to support metadata]"
    ["STAGE_24"]="R analyses - MaAsLin3/melonnpan statistical analysis [Currently out of order- Still being implemented]"
)

# Ordered list of all stages
ALL_STAGES=(
    STAGE_01 STAGE_02 STAGE_03 STAGE_04 STAGE_05 STAGE_06
    STAGE_07A STAGE_07B STAGE_07C STAGE_07D STAGE_07E
    STAGE_08 STAGE_09 STAGE_10 STAGE_11 STAGE_12
    STAGE_13 STAGE_13B STAGE_13C STAGE_13D STAGE_13E STAGE_14 STAGE_15 STAGE_16
    STAGE_17
    STAGE_18
    STAGE_18B
    STAGE_19 STAGE_20 STAGE_21 STAGE_21B STAGE_22 STAGE_23 STAGE_24
)

read -p "Run the full default workflow? (y/n): " full_workflow

if [[ "$full_workflow" =~ ^[Yy]$ ]]; then
    SELECTED_STAGES=("${ALL_STAGES[@]}")
    echo "$(date): Full workflow selected — all stages will run."
else
    echo ""
    echo "Available stages:"
    echo "-----------------------------------------------------------------"
    for stage in "${ALL_STAGES[@]}"; do
        printf "  %-12s %s\n" "${stage}:" "${STAGE_NAMES[$stage]}"
    done
    echo "-----------------------------------------------------------------"
    echo ""
    echo "Enter the stages you want to run as a space-separated list."
    echo "Example: STAGE_01 STAGE_02 STAGE_03 STAGE_10 STAGE_12"
    echo "Binning stages: STAGE_07A STAGE_07B STAGE_07C STAGE_07D STAGE_07E"
    echo "Full Binning stages: STAGE_01 STAGE_02 STAGE_03 STAGE_04 STAGE_05 STAGE_06 STAGE_07A STAGE_07B STAGE_07C STAGE_07D STAGE_07E STAGE_08 STAGE_09"
    echo ""

    SELECTED_STAGES=()
    while true; do
        read -p "Stages to run: " -a user_input

        INVALID=()
        VALID=()
        for entry in "${user_input[@]}"; do
            # Normalise to uppercase
            entry="${entry^^}"
            # Accept numeric shorthand for non-binning stages (e.g. 1 -> STAGE_01)
            if [[ "$entry" =~ ^[0-9]+$ ]]; then
                entry=$(printf "STAGE_%02d" "$entry")
            fi

            if [[ -n "${STAGE_NAMES[$entry]}" ]]; then
                VALID+=("$entry")
            else
                INVALID+=("$entry")
            fi
        done

        if [ ${#INVALID[@]} -gt 0 ]; then
            echo "Unrecognised stages: ${INVALID[*]}"
            echo "Please enter valid stage names."
            continue
        fi

        if [ ${#VALID[@]} -eq 0 ]; then
            echo "No stages selected. Please enter at least one stage."
            continue
        fi

        # Sort selected stages into pipeline order
        SELECTED_STAGES=()
        for stage in "${ALL_STAGES[@]}"; do
            for selected in "${VALID[@]}"; do
                if [ "$stage" = "$selected" ]; then
                    SELECTED_STAGES+=("$stage")
                    break
                fi
            done
        done

        echo ""
        echo "Selected stages (in order):"
        for stage in "${SELECTED_STAGES[@]}"; do
            printf "  %-12s %s\n" "${stage}:" "${STAGE_NAMES[$stage]}"
        done
        echo ""

        # Dependency warnings
        declare -A DEPENDS=(
            ["STAGE_03"]="STAGE_02"
            ["STAGE_04"]="STAGE_03"
            ["STAGE_05"]="STAGE_03"
            ["STAGE_06"]="STAGE_05"
            ["STAGE_07A"]="STAGE_05"
            ["STAGE_07B"]="STAGE_07A"
            ["STAGE_07C"]="STAGE_07B"
            ["STAGE_07D"]="STAGE_07C"
            ["STAGE_07E"]="STAGE_07D"
            ["STAGE_08"]="STAGE_07E"
            ["STAGE_09"]="STAGE_08"
            ["STAGE_10"]="STAGE_03"
            ["STAGE_11"]="STAGE_03"
            ["STAGE_12"]="STAGE_03"
            ["STAGE_13"]="STAGE_07E"
            ["STAGE_13B"]="STAGE_13"
            ["STAGE_13C"]="STAGE_13"
            ["STAGE_13D"]="STAGE_13"
            ["STAGE_13E"]="STAGE_13"
            ["STAGE_14"]="STAGE_13"
            ["STAGE_15"]="STAGE_03"
            ["STAGE_16"]="STAGE_15"
            ["STAGE_17"]="STAGE_14"
            ["STAGE_18"]="STAGE_17"
            ["STAGE_18B"]="STAGE_18"
            ["STAGE_19"]="STAGE_13"
            ["STAGE_20"]="STAGE_03"
            ["STAGE_21"]="STAGE_05"
            ["STAGE_21B"]="STAGE_21"
            ["STAGE_22"]="STAGE_21B"
            ["STAGE_23"]="STAGE_10"
            ["STAGE_24"]="STAGE_20"
        )

        WARNED=false
        for stage in "${SELECTED_STAGES[@]}"; do
            dep="${DEPENDS[$stage]}"
            if [ -n "$dep" ]; then
                found=false
                for s in "${SELECTED_STAGES[@]}"; do
                    [ "$s" = "$dep" ] && found=true && break
                done
                if [ "$found" = false ]; then
                    echo "WARNING: ${stage} depends on ${dep} which is not selected."
                    echo "         Outputs from ${dep} must already exist in ${DIRECTORY}."
                    WARNED=true
                fi
            fi
        done

        [ "$WARNED" = true ] && echo ""

        read -p "Proceed with these stages? (y/n): " confirm_stages
        if [[ "$confirm_stages" =~ ^[Yy]$ ]]; then
            break
        fi
    done
fi

echo "$(date): ${#SELECTED_STAGES[@]} stages selected."

printf '%s\n' "${SELECTED_STAGES[@]}" > "${DIRECTORY}/selected_stages.txt"
echo "$(date): Selected stages written to ${DIRECTORY}/selected_stages.txt"
export SELECTED_STAGES

# ============================================================================
# SECTION 4C: MEDI RESOURCE WARNING
# ============================================================================
for stage in "${SELECTED_STAGES[@]}"; do
    if [ "$stage" = "STAGE_11" ]; then
        echo ""
        echo "*** WARNING: STAGE_11 (MEDI) selected ***"
        echo "MEDI's reference database is ~500GB on disk, and quantification against"
        echo "it is memory- and CPU-intensive. Recommended: as much memory as your"
        echo "allocation allows (ideally 450G+, the current configurable maximum)."
        echo "Your current settings: ${MEM_SIZE}G memory, ${NUM_CORES} cores."
        if [ "$MEM_SIZE" -lt 400 ]; then
            echo ""
            echo "WARNING: ${MEM_SIZE}G may be insufficient for MEDI and could cause the job to fail or be OOM-killed."
            read -p "Continue anyway? (y/n): " medi_confirm
            if [[ ! "$medi_confirm" =~ ^[Yy]$ ]]; then
                echo "Cancelled. Re-run and increase memory allocation before selecting STAGE_11."
                return 1
            fi
        fi
        break
    fi
done

# ============================================================================
# SECTION 5: DATABASE PATH SETUP
# ============================================================================
DATABASES="/rds/projects/h/hallly-microbiome/software/databases"

# ============================================================================
# EXPORT ALL CONFIGURATION VARIABLES
# ============================================================================
export DIRECTORY
export EXTENSION
export WALL_TIME
export WALL_TIME_PADDED
export MEM_SIZE
export NUM_CORES
export CONCURRENT
export EMAIL
export DATABASES

echo "$(date): Configuration complete and variables exported."