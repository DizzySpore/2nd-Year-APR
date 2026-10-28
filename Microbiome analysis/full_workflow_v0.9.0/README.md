# Metagenomics Workflow v0.9.0

**Last updated:** 04/06/2026

A modular SLURM-based pipeline for end-to-end metagenomic analysis of paired-end sequencing data.

---

## Before you begin
- Paired-end reads (R1/R2) in a single directory
- SLURM cluster with Conda available
- Required databases at `<DATABASES_PATH>`

---

## Run the workflow

```bash
chmod +x *.sh
./main_workflow.sh
```

You will be prompted for:

| Prompt | Example |
|--------|---------|
| Working directory | `<WORKING_DIRECTORY>` or Enter for current |
| File extension | `.fastq.gz` |
| Wall time (hours) | `48` |
| Memory (GB) | `128` |
| CPU cores | `16` |
| Concurrent jobs | `8` |
| Email | `<EMAIL_ADDRESS>` |
| Run full workflow? | `y` / `n` (then select stages) |

Review the summary, type `y` to submit.

---

## Stage selection

When asked *"Run the full default workflow?"*, typing `n` brings up an interactive stage selector. You can enter any combination of stages as a space-separated list:

```
Available stages:
-----------------------------------------------------------------
  STAGE_01     FastQC - Quality check on raw reads
  STAGE_02     fastp - Read trimming and quality control
  STAGE_03     KneadData - Human read decontamination
  STAGE_04     KneadData inspection - FastQC on decontaminated reads
  STAGE_05     MEGAHIT - Metagenomic assembly
  STAGE_06     QUAST assembly - Assembly quality assessment
  STAGE_07A    Bowtie2 - Read mapping to assembly
  STAGE_07B    Binning - CONCOCT, MetaBAT2, COMEbin, SemiBin2
  STAGE_07C    Bin refinement - Hybrid creation, quality assessment and selection
  STAGE_07D    Reassembly - SPAdes strict and permissive
  STAGE_07E    Bin selection - Best version selection and labelling
  STAGE_08     GTDB-Tk - Bin taxonomic classification
  STAGE_09     QUAST bins - Bin quality assessment
  STAGE_10     Kraken2/Bracken - Read-level taxonomic profiling
  STAGE_11     Kraken2 cFMD - Food/custom database profiling
  STAGE_12     MetaPhlAn - Marker gene taxonomic profiling
  STAGE_13     Kaiju - Protein-level taxonomic profiling
  STAGE_14     Bin filtering - Quality filter high-confidence bins
  STAGE_15     EukFinder - Eukaryotic sequence detection
  STAGE_16     Kraken2 Euk - Eukaryotic taxonomic profiling
  STAGE_17     Bakta - Prokaryotic gene annotation
  STAGE_18     QUAST Bakta - Annotation quality assessment
  STAGE_19     antiSMASH - Secondary metabolite detection
  STAGE_20     HUMAnN3 reads - Functional profiling on reads
  STAGE_21     HUMAnN3 contigs - Functional profiling on contigs
  STAGE_22     Combining - Merge HUMAnN3 outputs
  STAGE_23     Visualization - Kraken2/Bracken visualizations
  STAGE_24     R analyses - MaAsLin3/melonnpan statistical analysis
-----------------------------------------------------------------
Stages to run: STAGE_01 STAGE_02 STAGE_03 STAGE_10 STAGE_12
```

You can also enter stage numbers as shorthand (e.g. `1 2 3` instead of `STAGE_01 STAGE_02 STAGE_03`). Stages are always executed in pipeline order regardless of input order.

If a selected stage depends on one that wasn't selected, you'll see a warning — this is allowed if the required outputs already exist in your working directory from a previous run.

---

## Monitor your jobs

```bash
squeue -j $JOB_ID                  # Check status
tail -f slurm-${JOB_ID}_1.out      # Live output
ls status/*.done                   # Completed samples
ls status/*.failed                 # Failed samples
```

---

## Outputs

| Directory | Contents |
|-----------|----------|
| `qc_reads/` | FastQC and MultiQC reports |
| `megahit_assembly/` | Assembled contigs |
| `final_bins/` | All binning outputs |
| `good_bins/` | Filtered high-quality bins (≥70% complete, ≤10% contamination) |
| `gtdbtk-out/` | GTDB-Tk taxonomy |
| `kraken2_bracken_out/` | Kraken2/Bracken profiles |
| `metaphlan_out/` | MetaPhlAn profiles |
| `bakta_out/` | Gene annotations |
| `antismash_out/` | Secondary metabolite clusters |
| `humann_reads_output/` | Functional profiles from reads |
| `humann_contigs_output/` | Functional profiles from contigs |
| `visualizations/` | Abundance plots |
| `status/` | `*.done` and `*.failed` per sample |

---

## Troubleshooting

**No samples found** — check your extension matches your files: `ls *.fastq.gz`

**Job fails mid-run** — the workflow auto-resumes from the last completed stage on resubmission. Check `slurm-${JOB_ID}_1.err` for the cause, fix it, then rerun `./main_workflow.sh`.

**Out of memory** — reduce concurrent jobs or increase memory allocation.

---

## Further documentation

| File | Contents |
|------|----------|
| `WORKFLOW_SUMMARY.md` | Module architecture, all stages, dependency map, checkpoint system |
| `DATA_FLOW.md` | Execution flow diagrams, variable substitution, directory structure |
