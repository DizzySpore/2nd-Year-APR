#!/bin/bash

################################################################################
# METAGENOMICS WORKFLOW - SLURM TEMPLATE GENERATOR (02_slurm_template.sh)
################################################################################
#
# Purpose: Generates the SLURM batch script for submission to the HPC cluster.
#          Variables known at generation time (config) are expanded directly.
#          Variables that must be evaluated at runtime on the compute node
#          are escaped so they are written literally into the generated script.
#
# Usage: source ./02_slurm_template.sh (called by main orchestrator)
#
# Dependencies: Requires all variables from 00_config.sh, 01_sample_discovery.sh
#               and 03_checkpoint.sh
#
################################################################################

# Verify required variables
if [ -z "$DIRECTORY" ] || [ -z "$N" ] || [ -z "$EMAIL" ]; then
    echo "Error: Required variables not set. Run previous modules first."
    return 1
fi

echo "$(date): Generating SLURM job template..."

cat > "${DIRECTORY}/metagenomics_workflow.sh" << EOF
#!/bin/bash

################################################################################
# METAGENOMICS WORKFLOW - SLURM EXECUTION SCRIPT
################################################################################
#
# Pipeline Stages:
# STAGE_01   FastQC               - Quality check on raw reads
# STAGE_02   fastp                - Read trimming and quality control
# STAGE_03   KneadData            - Human read decontamination
# STAGE_04   KneadData inspection - FastQC on decontaminated reads
# STAGE_05   MEGAHIT              - Metagenomic assembly
# STAGE_06   QUAST assembly       - Assembly quality assessment
# STAGE_07A  Bowtie2              - Read mapping to assembly
# STAGE_07B  Binning              - CONCOCT, MetaBAT2, COMEbin, SemiBin2
# STAGE_07C  Bin refinement       - Bin quality assessment and Hybrid bin creation and selection
# STAGE_07D  Reassembly           - SPAdes strict and permissive
# STAGE_07E  Bin selection        - Best version selection and labelling
# STAGE_08   GTDB-Tk              - Bin taxonomic classification
# STAGE_09   QUAST bins           - Bin quality assessment
# STAGE_10   Kraken2/Bracken      - Read-level taxonomic profiling
# STAGE_11   MEDI                 - Dietary/food metagenomic quantification (Nextflow, batch)
# STAGE_12   MetaPhlAn            - Marker gene taxonomic profiling
# STAGE_13   Bakta                - Prokaryotic gene annotation
# STAGE_13B  QUAST Bakta          - Annotation quality assessment
# STAGE_13C  dbCAN                - CAZyme annotation (Bakta proteins)
# STAGE_13D  eggNOG-mapper        - Functional/orthology annotation (Bakta proteins)
# STAGE_13E  abricate             - AMR/virulence gene screening (card, vfdb, resfinder, plasmidfinder)
# STAGE_14   Bin filtering        - Quality filter high-confidence bins
# STAGE_15   EukFinder            - Eukaryotic sequence detection
# STAGE_16   Kraken2 Euk          - Eukaryotic taxonomic profiling
# STAGE_17   Panaroo              - Targeted pangenome analysis (Nextflow-free, batch)
# STAGE_18   Snippy               - Within-group SNP calling on Panaroo taxon groups (batch)
# STAGE_18B  RAxML-NG             - Phylogenetic tree per taxon group (batch)
# STAGE_19   antiSMASH            - Secondary metabolite detection
# STAGE_20   HUMAnN3 reads        - Functional profiling on reads
# STAGE_21   HUMAnN3 contigs      - Functional profiling on contigs
# STAGE_21B  HUMAnN3 regroup      - CPM renorm + regroup/rename to rxn/EC/KO
# STAGE_22   Combining            - Merge HUMAnN3 outputs
# STAGE_23   Visualization        - Kraken2/Bracken visualizations
# STAGE_24   R analyses           - MaAsLin3/melonnpan statistical analysis
#
################################################################################

#SBATCH --time=${WALL_TIME_PADDED}:00:00
#SBATCH --nodes=1
#SBATCH --mem=${MEM_SIZE}G
#SBATCH --ntasks=${NUM_CORES}
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${EMAIL}
#SBATCH --array=1-${N}%${CONCURRENT}
#SBATCH -J metagenomics_workflow
#SBATCH --error=slurm-%A_%a.err
#SBATCH --output=slurm-%A_%a.out
#SBATCH --chdir=${DIRECTORY}

# -- Environment ---------------------------------------------------------------
CONDA_BASE=\$(conda info --base 2>/dev/null)
source \${CONDA_BASE}/etc/profile.d/conda.sh

# -- Selected stages -----------------------------------------------------------
SELECTED_STAGES_STR=\$(cat "${DIRECTORY}/selected_stages.txt" | tr '\n' ' ')
echo "\$(date): Selected stages: \${SELECTED_STAGES_STR}"

# -- Sample resolution ---------------------------------------------------------
SAMPLE=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "${DIRECTORY}/samples.txt")

if [ -z "\$SAMPLE" ]; then
    echo "Error: Could not resolve sample for task \${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

echo "\$(date): Task \${SLURM_ARRAY_TASK_ID} processing sample: \${SAMPLE}"

# -- Resource variables --------------------------------------------------------
NUM_THREADS=\${SLURM_NTASKS}
MEMORY=\$(( SLURM_MEM_PER_NODE / 1024 ))

echo "\$(date): Using \${NUM_THREADS} threads and \${MEMORY}G memory"

# -- Directory structure -------------------------------------------------------
DIRECTORY="${DIRECTORY}"
DATABASES="${DATABASES}"
EXTENSION="${EXTENSION}"
STATUS_DIR="${DIRECTORY}/status"

QC_READS="${DIRECTORY}/qc_reads"
DEHUMAN_READS="${DIRECTORY}/dehuman_reads"
KNEADDATA_OUT="${DIRECTORY}/dehuman_reads/kneaddata_output_\${SAMPLE}"
ASSEMBLY_DIR="${DIRECTORY}/megahit_assembly/\${SAMPLE}assembly"
ASSEMBLY_FILE="${DIRECTORY}/megahit_assembly/\${SAMPLE}assembly/final.contigs.fa"
QUAST_ASSEMBLY_OUT="${DIRECTORY}/quast_assembly/\${SAMPLE}/\${SAMPLE}_quast_results"
FASTQC_RAWS_DIR="${DIRECTORY}/qc_reads/fastqc_raw/\${SAMPLE}"
FASTQC_CLEAN_DIR="${DIRECTORY}/qc_reads/fastqc_clean/\${SAMPLE}"
MULTIQC_RAWS_DIR="${DIRECTORY}/qc_reads/multiqc_raw"
MULTIQC_CLEAN_DIR="${DIRECTORY}/qc_reads/multiqc_clean"

# -- Binning pipeline directories ----------------------------------------------
MAPPING_DIR="${DIRECTORY}/mapping/\${SAMPLE}"
BINNING_DIR="${DIRECTORY}/binning/\${SAMPLE}"
REFINEMENT_DIR="${DIRECTORY}/bin_refinement/\${SAMPLE}"
REASSEMBLY_DIR="${DIRECTORY}/bin_reassembly/\${SAMPLE}"
FINAL_BINS_DIR="${DIRECTORY}/final_bins/\${SAMPLE}"
UTILS="${SCRIPT_DIR}"

# -- Downstream analysis directories ------------------------------------------
QUAST_BINS_OUT="${DIRECTORY}/quast_bins/\${SAMPLE}/\${SAMPLE}_quast_results"
METAQUAST_REFS_FILE="${DIRECTORY}/final_bins/\${SAMPLE}/metaquast_references.txt"
KRAKEN2_OUT="${DIRECTORY}/kraken2_bracken_out"
METAPHLAN_OUT="${DIRECTORY}/metaphlan_out"
GTDBTK_OUT="${DIRECTORY}/gtdbtk-out/gtdbtk-\${SAMPLE}-out"
CROSSCHECK_OUT="${DIRECTORY}/crosscheck_out"
GOOD_BINS_DIR="${DIRECTORY}/good_bins/\${SAMPLE}"
EUKFINDER_OUT="${DIRECTORY}/eukfinder_out/\${SAMPLE}"
BAKTA_OUT="${DIRECTORY}/bakta_out/\${SAMPLE}"
BAKTA_QUAST_OUT="${DIRECTORY}/bakta_quast_out/\${SAMPLE}"
DBCAN_OUT="${DIRECTORY}/dbcan_out/\${SAMPLE}"
EGGNOG_OUT="${DIRECTORY}/eggnog_out/\${SAMPLE}"
ABRICATE_OUT="${DIRECTORY}/abricate_out/\${SAMPLE}"
ANTISMASH_OUT="${DIRECTORY}/antismash_out/\${SAMPLE}"
HUMANN_READS_DIR="${DIRECTORY}/humann_reads_output/\${SAMPLE}"
HUMANN_CONTIGS_DIR="${DIRECTORY}/humann_contigs_output/\${SAMPLE}"
VISUALIZATION_OUT="${DIRECTORY}/visualizations"
PANAROO_OUT="${DIRECTORY}/panaroo_out"
SNIPPY_OUT="${DIRECTORY}/snippy_out"
RAXML_OUT="${DIRECTORY}/raxml_out"
COMPLETENESS=50
CONTAMINATION=10
STAGE_SCRIPTS="${SCRIPT_DIR}/stages"

# -- Create directories --------------------------------------------------------
mkdir -p \${QC_READS} \${DEHUMAN_READS} \${KNEADDATA_OUT} \${ASSEMBLY_DIR} \
         \${QUAST_ASSEMBLY_OUT} \${FASTQC_RAWS_DIR} \${FASTQC_CLEAN_DIR} \
         \${MULTIQC_RAWS_DIR} \${MULTIQC_CLEAN_DIR} \
         \${MAPPING_DIR} \${BINNING_DIR} \
         \${REFINEMENT_DIR} \${REASSEMBLY_DIR} \${FINAL_BINS_DIR} \
         \${QUAST_BINS_OUT} \${KRAKEN2_OUT} \${METAPHLAN_OUT} \
         \${GTDBTK_OUT} \${GOOD_BINS_DIR} \${EUKFINDER_OUT} \
         \${BAKTA_OUT} \${BAKTA_QUAST_OUT} \${DBCAN_OUT} \${EGGNOG_OUT} \${ABRICATE_OUT} \${ANTISMASH_OUT} \
         \${HUMANN_READS_DIR} \${HUMANN_CONTIGS_DIR} \
         \${VISUALIZATION_OUT} \${PANAROO_OUT} \${SNIPPY_OUT} \${RAXML_OUT} \${STATUS_DIR} \${DIRECTORY}/medi_run_output

# -- Input file identification -------------------------------------------------
echo "\$(date): Searching for paired files for sample: \${SAMPLE}"

shopt -s nullglob

READ1=""
READ2=""

# === IMPROVED DISCOVERY ===
echo "DEBUG: Sample name  = [\${SAMPLE}]"
echo "DEBUG: Extension   = [\${EXTENSION}]"
echo "DEBUG: Directory   = [\${DIRECTORY}]"

# Note the escaped variables below so they evaluate at cluster RUNTIME
mapfile -t FILES < <(find "\${DIRECTORY}" -maxdepth 1 -type f -name "\${SAMPLE}*\${EXTENSION}" -print)

echo "DEBUG: Files found by find = \${#FILES[@]}"

for file in "\${FILES[@]}"; do
    echo "DEBUG: Checking file -> \$file"


# Separator before the read ID can be '.' or '_'
    if [[ "\$file" =~ [._]R1(\.filtered)?(\.fastq|\.fq)?(\.gz)?$ ]] || \
       [[ "\$file" =~ [._]1(\.filtered)?(\.fastq|\.fq)?(\.gz)?$ ]]; then
        READ1="\$file"
        echo "DEBUG: Matched as R1"
    elif [[ "\$file" =~ [._]R2(\.filtered)?(\.fastq|\.fq)?(\.gz)?$ ]] || \
         [[ "\$file" =~ [._]2(\.filtered)?(\.fastq|\.fq)?(\.gz)?$ ]]; then
        READ2="\$file"
        echo "DEBUG: Matched as R2"
    fi
done

# Validation
if [ -z "\$READ1" ] || [ -z "\$READ2" ]; then
    echo "Error: Could not find paired files for \${SAMPLE}" >&2
    echo "Files found by find:" >&2
    printf '%s\n' "\${FILES[@]}" >&2
    exit 1
fi

echo "\$(date): Input files identified: \$READ1 and \$READ2"

# -- Error handling ------------------------------------------------------------
trap 'echo "Error at line \$LINENO for \${SAMPLE}"; touch "\${STATUS_DIR}/\${SAMPLE}.failed"; exit 1' ERR

# -- Resume point --------------------------------------------------------------
RESUME_FILE="\${STATUS_DIR}/\${SAMPLE}.resume"
if [ -f "\$RESUME_FILE" ]; then
    RESUME_FROM=\$(cat "\$RESUME_FILE")
else
    RESUME_FROM="STAGE_01"
fi

echo "\$(date): Resuming from stage: \${RESUME_FROM}"

# -- Stage execution -----------------------------------------------------------
run_stage() {
    local stage="\$1"
    local script="\$2"
    if [[ " \${SELECTED_STAGES_STR} " =~ " \${stage} " ]]; then
        source "\${script}"
    else
        echo "\$(date): Skipping \${stage} — not in selected workflow."
    fi
}

case "\$RESUME_FROM" in
    "STAGE_01")  run_stage "STAGE_01"  "\${STAGE_SCRIPTS}/stage_01_fastqc.sh"               ;&
    "STAGE_02")  run_stage "STAGE_02"  "\${STAGE_SCRIPTS}/stage_02_fastp.sh"                ;&
    "STAGE_03")  run_stage "STAGE_03"  "\${STAGE_SCRIPTS}/stage_03_kneaddata.sh"            ;&
    "STAGE_04")  run_stage "STAGE_04"  "\${STAGE_SCRIPTS}/stage_04_kneaddata_inspection.sh" ;&
    "STAGE_05")  run_stage "STAGE_05"  "\${STAGE_SCRIPTS}/stage_05_megahit.sh"              ;&
    "STAGE_06")  run_stage "STAGE_06"  "\${STAGE_SCRIPTS}/stage_06_quast_assembly.sh"       ;&
    "STAGE_07A") run_stage "STAGE_07A" "\${STAGE_SCRIPTS}/stage_07a_bowtie2_mapping.sh"     ;&
    "STAGE_07B") run_stage "STAGE_07B" "\${STAGE_SCRIPTS}/stage_07b_binning.sh"             ;&
    "STAGE_07C") run_stage "STAGE_07C" "\${STAGE_SCRIPTS}/stage_07c_refinement.sh"  	    ;&
    "STAGE_07D") run_stage "STAGE_07D" "\${STAGE_SCRIPTS}/stage_07d_reassembly.sh"  	    ;&
    "STAGE_07E") run_stage "STAGE_07E" "\${STAGE_SCRIPTS}/stage_07e_bin_selection.sh"       ;&
    "STAGE_08")  run_stage "STAGE_08"  "\${STAGE_SCRIPTS}/stage_08_gtdbtk.sh"               ;&
    "STAGE_09")  run_stage "STAGE_09"  "\${STAGE_SCRIPTS}/stage_09_quast_bins.sh"           ;&
    "STAGE_10")  run_stage "STAGE_10"  "\${STAGE_SCRIPTS}/stage_10_kraken2.sh"              ;&
    "STAGE_11")  run_stage "STAGE_11"  "\${STAGE_SCRIPTS}/stage_11_medi.sh"                 ;&
    "STAGE_12")  run_stage "STAGE_12"  "\${STAGE_SCRIPTS}/stage_12_metaphlan.sh"            ;&
    "STAGE_13")  run_stage "STAGE_13"  "\${STAGE_SCRIPTS}/stage_13_bakta.sh"                ;&
    "STAGE_13B") run_stage "STAGE_13B" "\${STAGE_SCRIPTS}/stage_13b_quast_bakta.sh"         ;&
    "STAGE_13C") run_stage "STAGE_13C" "\${STAGE_SCRIPTS}/stage_13c_dbcan.sh"               ;&
    "STAGE_13D") run_stage "STAGE_13D" "\${STAGE_SCRIPTS}/stage_13d_eggnog.sh"              ;&
    "STAGE_13E") run_stage "STAGE_13E" "\${STAGE_SCRIPTS}/stage_13e_abricate.sh"            ;&
    "STAGE_14")  run_stage "STAGE_14"  "\${STAGE_SCRIPTS}/stage_14_bin_filter.sh"           ;&
    "STAGE_15")  run_stage "STAGE_15"  "\${STAGE_SCRIPTS}/stage_15_eukfinder.sh"            ;&
    "STAGE_16")  run_stage "STAGE_16"  "\${STAGE_SCRIPTS}/stage_16_kraken2_euk.sh"          ;&
    "STAGE_17")  run_stage "STAGE_17"  "\${STAGE_SCRIPTS}/stage_17_panaroo.sh"              ;&
    "STAGE_18")  run_stage "STAGE_18"  "\${STAGE_SCRIPTS}/stage_18_snippy.sh"               ;&
    "STAGE_18B") run_stage "STAGE_18B" "\${STAGE_SCRIPTS}/stage_18b_raxml.sh"               ;&
    "STAGE_19")  run_stage "STAGE_19"  "\${STAGE_SCRIPTS}/stage_19_antismash.sh"            ;&
    "STAGE_20")  run_stage "STAGE_20"  "\${STAGE_SCRIPTS}/stage_20_humann_reads.sh"         ;&
    "STAGE_21")  run_stage "STAGE_21"  "\${STAGE_SCRIPTS}/stage_21_humann_contigs.sh"       ;&
    "STAGE_21B") run_stage "STAGE_21B" "\${STAGE_SCRIPTS}/stage_21b_humann_regroup.sh"      ;&
    "STAGE_22")  run_stage "STAGE_22"  "\${STAGE_SCRIPTS}/stage_22_combining.sh"            ;&
    "STAGE_23")  run_stage "STAGE_23"  "\${STAGE_SCRIPTS}/stage_23_visualization.sh"        ;&
    "STAGE_24")  run_stage "STAGE_24"  "\${STAGE_SCRIPTS}/stage_24_r_analyses.sh"           ;&
    *)
        echo "\$(date): All stages complete for \${SAMPLE}"
        touch "\${STATUS_DIR}/\${SAMPLE}.done"
        ;;
esac

EOF

echo "$(date): SLURM job template written to ${DIRECTORY}/metagenomics_workflow.sh"
echo "$(date): Ready for submission."