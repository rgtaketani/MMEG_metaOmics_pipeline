# Per-sample pipelines

[← Back to the main README](../README.md)

| Script | Data type | Steps |
|---|---|---|
| [`mtx_profile_sample.sh`](mtx_profile_sample.sh) | Metatranscriptome | fastp → RiboDetector → DIAMOND blastx → `daa-meganizer` → `daa2info` |
| [`mgx_profile_sample.sh`](mgx_profile_sample.sh) | Metagenome | fastp with deduplication → DIAMOND blastx → `daa-meganizer` → `daa2info` |

```bash
bash mtx_profile_sample.sh <reads_R1.fastq.gz> <reads_R2.fastq.gz>
```

Each script processes **one sample**. The two read files are the only
arguments. Reference databases and run-time options are set with environment
variables, which are documented in each script's header:

| Variable | Default | Meaning |
|---|---|---|
| `DIAMOND_DB` | *(required)* | DIAMOND index of NCBI-nr |
| `MEGAN_MAPDB` | *(required)* | MEGAN Ultimate mapping database (`-ue`) |
| `MEGAN_SYN2TAXA` | *(empty)* | Optional SILVA synonyms file |
| `OUTDIR` | current directory | Parent folder; the sample's outputs go to `<OUTDIR>/<sample>_MTX/` or `_MGX/` |
| `THREADS` | all cores | Threads for every step |
| `SAMPLE` | from R1 file name | Overrides the sample name |
| `DIAMOND_BLOCK` / `DIAMOND_CHUNKS` | 8 / 1 | DIAMOND memory–speed trade-off |
| `MAX_JAVA_JOBS` | 2 | Concurrent `daa2info` exports |
| `KEEP_INTERMEDIATE` | 0 | Metatranscriptomes only: keep trimmed reads |

All tools must be on `PATH` ([conda environment](../envs/pipeline.yml),
plus MEGAN Ultimate's `tools/` folder). Completed steps leave `.done`
markers, so a re-run resumes where it stopped.

The `SHARED PARAMETERS` block must stay identical in both scripts. Its values
and the reasons are given in the main README under
[Design notes](../README.md#design-notes). The output files are listed under
[Pipeline outputs](../README.md#pipeline-outputs-per-sample).
