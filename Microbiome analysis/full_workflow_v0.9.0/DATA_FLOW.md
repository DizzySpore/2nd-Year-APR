# Metagenomics Workflow v0.9.0 — Data Flow

**Last updated:** 04/06/2026

---

## Module execution order

```
main_workflow.sh (local machine)
│
├── source 00_config.sh
│       Prompts user → exports config vars
│       Writes: selected_stages.txt
│
├── source 01_sample_discovery.sh
│       Scans DIRECTORY for *EXTENSION files
│       Strips R1/R2 suffixes → unique sample names
│       Writes: samples.txt
│
├── source 03_checkpoint.sh
│       Scans output dirs for completed samples
│       Removes done samples from samples.txt
│       Updates N
│
├── source 02_slurm_template.sh
│       Substitutes all config/sample vars into template
│       Writes: metagenomics_workflow.sh
│
└── sbatch metagenomics_workflow.sh
        │
        └── [SLURM array job — 1 task per sample]
                SAMPLE  = samples.txt[SLURM_ARRAY_TASK_ID]
                RESUME  = status/${SAMPLE}.resume
                STAGES  = selected_stages.txt
                │
                └── Fallthrough case → source stage scripts
```

---

## Per-sample pipeline

```
R1.fastq.gz + R2.fastq.gz
        │
        ▼
   STAGE_01  FastQC (raw QC)
        │
        ▼
   STAGE_02  fastp (trimming)
        │
        ▼
   STAGE_03  KneadData (decontamination)
        │
   ┌────┴────────────────────────────────┐
   ▼                                     ▼
STAGE_04                            STAGE_05
FastQC (clean QC)                   MEGAHIT (assembly)
                                         │
                                         ▼
                                    STAGE_06  QUAST (assembly QC)
                                         │
                                         ▼
                                    STAGE_07A  Bowtie2 (read mapping)
                                         │
                                         ▼
                                    STAGE_07B  Binning
                                    CONCOCT / MetaBAT2 / COMEbin / SemiBin2
                                         │
                                         ▼
                                    STAGE_07C  Bin refinement
                                    CheckM2 QC → hybrid bin creation → selection
                                         │
                                         ▼
                                    STAGE_07D  Reassembly
                                    SPAdes strict + permissive
                                         │
                                         ▼
                                    STAGE_07E  Bin selection
                                    choose_best_bin.py → labelled final bins
                                         │
                              ┌──────────┼───────────┐
                              ▼          ▼           ▼
                         STAGE_08   STAGE_09   STAGE_14
                         GTDB-Tk    QUAST      Bin filter
                         taxonomy   bins QC    (≥70%/≤10%)
                                                    │
                              ┌─────────────────────┤
                              ▼                     ▼
                         STAGE_15              STAGE_17
                         EukFinder             Bakta annotation
                              │                     │
                              ▼                ┌────┴────┐
                         STAGE_16         STAGE_18   STAGE_19
                         Kraken2 Euk      QUAST      antiSMASH
                                          Bakta

STAGE_03 clean reads also feed:
        ├── STAGE_10  Kraken2/Bracken (standard DB)
        ├── STAGE_11  Kraken2 cFMD (custom DB)
        ├── STAGE_12  MetaPhlAn
        ├── STAGE_13  Kaiju
        └── STAGE_20  HUMAnN3 reads
                          └── STAGE_22  Combining
                          └── STAGE_24  MaAsLin3 / melonnpan

STAGE_05 contigs also feed:
        └── STAGE_21  HUMAnN3 contigs
                          └── STAGE_22  Combining

STAGE_10 output feeds:
        └── STAGE_23  Visualisation (Krona / ggplot2)
```

---

## Variable substitution: generation-time vs runtime

`02_slurm_template.sh` writes `metagenomics_workflow.sh` using a heredoc. Variables fall into two categories:

| Type | Example | How written | When evaluated |
|------|---------|-------------|----------------|
| Config (known now) | `${DIRECTORY}`, `${EMAIL}`, `${N}` | Unescaped `${VAR}` | At generation time on login node |
| Runtime (SLURM/per-sample) | `${SAMPLE}`, `${SLURM_ARRAY_TASK_ID}` | Escaped `\${VAR}` | At execution time on compute node |

---

## Checkpoint system

```
Stage completes
      │
      └── Output file written to known path
                │
                └── 03_checkpoint.sh scans all paths before next submission
                          │
                          ├── Sample complete?  → removed from samples.txt
                          └── Sample partial?   → resume point written to
                                                   status/${SAMPLE}.resume

On compute node:
      RESUME_FROM=$(cat status/${SAMPLE}.resume)   # e.g. STAGE_05
      case "$RESUME_FROM" in
          "STAGE_01") run_stage ... ;&   # fallthrough
          "STAGE_02") run_stage ... ;&
          ...                            # executes from resume point onward
      esac
```

---

## Directory structure (runtime)

```
DIRECTORY/
├── samples.txt                  Sample names (one per line, SLURM array index)
├── selected_stages.txt          Stages to run (written by 00_config.sh)
├── metagenomics_workflow.sh     Generated SLURM script
├── slurm-JOBID_TASKID.out/.err  Per-task logs
├── status/
│   ├── sample1.done             Written on successful completion
│   ├── sample1.resume           Last completed stage (for resume)
│   └── sample2.failed           Written on error trap
├── qc_reads/
├── dehuman_reads/
├── megahit_assembly/
├── mapping/
├── binning/
├── bin_refinement/
├── bin_reassembly/
├── final_bins/
├── good_bins/
├── gtdbtk-out/
├── kraken2_bracken_out/
├── metaphlan_out/
├── kaiju_out/
├── eukfinder_out/
├── bakta_out/
├── antismash_out/
├── humann_reads_output/
├── humann_contigs_output/
└── visualizations/
```
