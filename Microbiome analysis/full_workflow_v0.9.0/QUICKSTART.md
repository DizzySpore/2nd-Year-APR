# Metagenomics Workflow — Quick Start

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

## Monitor your jobs

```bash
squeue -j $JOB_ID                  # Check status
tail -f slurm-${JOB_ID}_1.out      # Live output
ls status/*.done                   # Completed samples
ls status/*.failed                 # Failed samples
```

---

## Troubleshooting

**No samples found** — check your extension matches your files: `ls *.fastq.gz`

**Job fails mid-run** — the workflow auto-resumes from the last completed stage on resubmission. Check `slurm-${JOB_ID}_1.err` for the cause, fix it, then rerun `./main_workflow.sh`.

**Out of memory** — reduce concurrent jobs or increase memory allocation.
