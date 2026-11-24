# Tests and simulated data

[← Back to the main README](../README.md)

This folder contains two tests that run without sequencing data, reference
databases or a MEGAN licence.

| File | Purpose |
|---|---|
| [`run_pipeline_smoke_test.sh`](run_pipeline_smoke_test.sh) | Runs both [pipeline scripts](../workflow/) on a five-read-pair input. The bioinformatics tools are replaced by the stubs in [`stubs/`](stubs/). The test checks that every output file is created, that both arms report the same shared-parameter fingerprint, and that MEGAN receives the shared parameters. It takes a few seconds. |
| [`run_analysis_test.sh`](run_analysis_test.sh) | Generates the simulated data set and renders both [R Markdown analyses](../analysis/) on it with quick settings. It takes a few minutes. |
| [`simulate_test_data.R`](simulate_test_data.R) | Writes synthetic per-sample files in the exact formats produced by the pipeline (`daa2info` tables and SeqKit read counts), plus `metadata.csv`. Requires base R only. |
| [`stubs/`](stubs/) | Stand-ins for `fastp`, `ribodetector_cpu`, `diamond`, `daa-meganizer` and `daa2info`, used only by the smoke test. |

## Running the tests

From the repository root:

```bash
# Pipeline control flow (needs seqkit and pigz; see ../envs/pipeline.yml)
bash test/run_pipeline_smoke_test.sh

# R analyses on simulated data (needs the R packages; see ../envs/install_R_packages.R)
bash test/run_analysis_test.sh
```

The analysis test writes `test/data/`, `test/results_MTX/`,
`test/results_MGX/` and two HTML reports. All of these are ignored by git
(see [`.gitignore`](../.gitignore)).

## Simulated design

The simulated data set has 72 samples: 2 countries × 2 locations per country
× 2 compartments × 3 treatments × 3 replicates. The same samples are used for
both the metatranscriptome and the metagenome arm.

| Feature | Simulated effect |
|---|---|
| Treatment | shifts about 30% of taxa and KOs |
| Compartment | shifts all taxa and KOs |
| Location | shifts taxa (random location-to-location variation) |
| Country | **none** beyond location variation |
| Sample `S05` | sequenced about 100× more shallowly (triggers the low-depth flag) |

Because country has no effect of its own, the location-level country test in
the analysis should not be significant. The free-permutation PERMANOVA, which
treats samples within a location as independent, typically reports a small
*p*-value for Country. The contrast illustrates why the analysis tests
country with whole locations as replicates (see Section 8.1 of either Rmd).
