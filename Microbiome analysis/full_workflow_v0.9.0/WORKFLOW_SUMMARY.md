# Metagenomics Workflow v0.9.0 — Architecture Summary

**Last updated:** 04/06/2026

---

## Module overview

| Script | Role |
|--------|------|
| `main_workflow.sh` | Orchestrator — runs all modules in sequence, submits to SLURM, monitors jobs |
| `00_config.sh` | Interactive configuration: directories, SLURM resources, email, stage selection |
| `01_sample_discovery.sh` | Finds paired read files, extracts sample names, writes `samples.txt` |
| `03_checkpoint.sh` | Scans output directories to find the last completed stage per sample |
| `02_slurm_template.sh` | Generates `metagenomics_workflow.sh` with all config values substituted |
| `stages/stage_XX_*.sh` | Individual stage scripts, sourced at runtime on the compute node |

---

## Execution flow

```
main_workflow.sh
    │
    ├─ source 00_config.sh          → exports all config vars + selected_stages.txt
    ├─ source 01_sample_discovery.sh → exports SAMPLE_LIST, N, writes samples.txt
    ├─ source 03_checkpoint.sh       → updates samples.txt to pending-only, updates N
    ├─ source 02_slurm_template.sh   → writes metagenomics_workflow.sh
    └─ sbatch metagenomics_workflow.sh
           │
           └─ [SLURM array: 1 job per sample]
                  ├─ Resolve SAMPLE from samples.txt[TASK_ID]
                  ├─ Read RESUME_FROM from status/${SAMPLE}.resume
                  ├─ Read SELECTED_STAGES from selected_stages.txt
                  └─ Fallthrough case statement → source stage scripts in order
```

---

## Pipeline stages

### Quality control
| Stage | Tool | Purpose |
|-------|------|---------|
| STAGE_01 | FastQC | Raw read quality reports |
| STAGE_02 | fastp | Adapter trimming and quality filtering |
| STAGE_03 | KneadData | Human read decontamination |
| STAGE_04 | FastQC | QC on decontaminated reads |

### Assembly
| Stage | Tool | Purpose |
|-------|------|---------|
| STAGE_05 | MEGAHIT | Metagenomic assembly |
| STAGE_06 | QUAST | Assembly quality assessment |

### Binning pipeline
| Stage | Tool | Purpose |
|-------|------|---------|
| STAGE_07A | Bowtie2 | Map reads back to assembly |
| STAGE_07B | CONCOCT, MetaBAT2, COMEbin, SemiBin2 | Multi-tool binning |
| STAGE_07C | CheckM2 + custom refinement | Bin quality assessment, hybrid bin creation and selection |
| STAGE_07D | SPAdes (strict + permissive) | Bin reassembly |
| STAGE_07E | choose_best_bin.py | Best version selection and labelling |

### Taxonomy
| Stage | Tool | Purpose |
|-------|------|---------|
| STAGE_08 | GTDB-Tk | Bin taxonomic classification |
| STAGE_09 | QUAST | Bin quality assessment |
| STAGE_10 | Kraken2 + Bracken | Read-level taxonomic profiling (standard DB) |
| STAGE_11 | Kraken2 | Custom food/cFMD database profiling |
| STAGE_12 | MetaPhlAn | Marker gene taxonomic profiling |
| STAGE_13 | Kaiju | Protein-level taxonomic profiling |
| STAGE_14 | Custom filter | Quality-filter bins (≥70% complete, ≤10% contamination) |
| STAGE_15 | EukFinder | Eukaryotic sequence detection |
| STAGE_16 | Kraken2 | Eukaryotic taxonomic profiling |

### Annotation
| Stage | Tool | Purpose |
|-------|------|---------|
| STAGE_17 | Bakta | Prokaryotic gene annotation |
| STAGE_18 | QUAST | Annotation quality assessment |
| STAGE_19 | antiSMASH | Secondary metabolite gene cluster detection |

### Functional profiling and analysis
| Stage | Tool | Purpose |
|-------|------|---------|
| STAGE_20 | HUMAnN3 | Functional profiling on reads |
| STAGE_21 | HUMAnN3 | Functional profiling on contigs |
| STAGE_22 | humann_join_tables | Merge HUMAnN3 outputs |
| STAGE_23 | R (Krona, ggplot2) | Kraken2/Bracken visualisations |
| STAGE_24 | MaAsLin3, melonnpan | Differential abundance and metabolic prediction |

---

## Stage selection

At configuration time (`00_config.sh`) the user chooses to run the full workflow or a custom subset of stages. The selection is written to `selected_stages.txt` and read by each compute job at runtime.

**Dependency map (direct upstream requirement):**

```
STAGE_03 → STAGE_02        STAGE_07A → STAGE_05
STAGE_04 → STAGE_03        STAGE_07B → STAGE_07A
STAGE_05 → STAGE_03        STAGE_07C → STAGE_07B
STAGE_06 → STAGE_05        STAGE_07D → STAGE_07C
STAGE_08 → STAGE_07E       STAGE_07E → STAGE_07D
STAGE_09 → STAGE_08        STAGE_10-13 → STAGE_03
STAGE_14 → STAGE_07E       STAGE_15 → STAGE_14
STAGE_16 → STAGE_15        STAGE_17-19 → STAGE_14
STAGE_20 → STAGE_03        STAGE_21 → STAGE_05
STAGE_22 → STAGE_20        STAGE_23 → STAGE_10
STAGE_24 → STAGE_20
```

---

## Checkpoint and resume system

Each stage writes its outputs to a known path. `03_checkpoint.sh` scans these paths before submission and removes already-complete samples from `samples.txt`, reducing the array size. On the compute node, each task reads the resume stage and continues from there.

```
Run 1:  STAGE_01 ✓  STAGE_02 ✓  STAGE_03 ✓  STAGE_05 ✗  (interrupted)
Run 2:  checkpoint detected → RESUME_FROM=STAGE_05 → continues from there
```

---

## Variable flow

```
00_config.sh  ──exports──►  DIRECTORY, EXTENSION, WALL_TIME_PADDED, MEM_SIZE,
                              NUM_CORES, CONCURRENT, EMAIL, DATABASES,
                              SELECTED_STAGES → selected_stages.txt

01_sample_discovery.sh  ──exports──►  SAMPLE_LIST, N → samples.txt

02_slurm_template.sh  ──substitutes all of the above into──►  metagenomics_workflow.sh
                         (runtime vars like $SAMPLE, $SLURM_ARRAY_TASK_ID are escaped)
```

---

## Utility scripts

| Script | Purpose |
|--------|---------|
| `binning_refiner.py` | Reconcile overlapping bin sets |
| `choose_best_bin.py` | Select highest-quality version of each bin post-reassembly |
| `consolidate_two_sets_of_bins.py` | Merge two binner outputs |
| `tournament_consolidate_bins.py` | Round-robin quality tournament across all binners |
| `deduplicate_bins.py` | Remove redundant bins by ANI |
| `select_best_bins.py` | Apply completeness/contamination thresholds |
| `smart_dereplicate_contigs.py` | Contig-level dereplication |
| `separate_strain_bins.py` | Split strain-level bin clusters |
| `identify_strain_pairs.py` | Detect near-identical strain bin pairs |
