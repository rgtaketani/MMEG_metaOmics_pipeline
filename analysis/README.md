# Statistical analyses

[← Back to the main README](../README.md)

| Report | Input | Code as PDF |
|---|---|---|
| [`analysis_MTX.Rmd`](analysis_MTX.Rmd) | Outputs of [`workflow/mtx_profile_sample.sh`](../workflow/mtx_profile_sample.sh) (`<sample>_MTX/` folders) | [`analysis_MTX_code.pdf`](../docs/analysis_MTX_code.pdf) |
| [`analysis_MGX.Rmd`](analysis_MGX.Rmd) | Outputs of [`workflow/mgx_profile_sample.sh`](../workflow/mgx_profile_sample.sh) (`<sample>_MGX/` folders) | [`analysis_MGX_code.pdf`](../docs/analysis_MGX_code.pdf) |

The two reports run the same analysis. They differ only in the
quality-control section (rRNA fraction for RNA, duplicate removal for DNA)
and in the interpretation notes. **A change to the shared code must be made
in both files.**

## Steps

1. Import of the `daa2info` tables and the [sample metadata](../config/README.md), with checks for mismatched or incomplete samples.
2. Quality control: reads retained per stage, assigned fraction, and a flag for low-depth samples.
3. Feature tables: taxa at a chosen rank (full lineages, so homonyms stay separate) and KOs mapped to KEGG categories.
4. Trimming: a feature is kept when it is detected in ≥ 20% of samples **and** has ≥ 20 reads in total.
5. Zero replacement (`zCompositions::cmultRepl`, GBM) and CLR transformation.
6. Aitchison PCA, PERMANOVA (global; treatment and compartment within locations; country with whole locations permuted; pairwise between treatments), ANOSIM and PERMDISP, for taxa and for KOs.
7. Stacked composition plots of taxa and KEGG categories.
8. DivNet Shannon and Simpson diversity with `betta_random` (location as a random effect), Wald tests and pairwise treatment comparisons within compartments.

## Rendering

Run from the repository root. Paths given in `params` are resolved from `analysis/`.

```r
rmarkdown::render("analysis/analysis_MTX.Rmd",
                  params = list(input_dir = "../megan_output",
                                metadata  = "../metadata.csv"))
```

| Parameter | Default | Meaning |
|---|---|---|
| `input_dir` | `megan_output` | Folder holding the per-sample output folders |
| `metadata` | `metadata.csv` | Sample sheet ([template](../config/metadata_template.csv)) |
| `out_dir` | `results_MTX` / `results_MGX` | Tables, figures and cache |
| `col_*` | `SampleID`, `Treatment`, … | Metadata column names |
| `treatment_levels` | *(empty)* | Treatment order, e.g. `"Control,T1,T2"`; the first is the reference |
| `min_prevalence`, `min_reads` | 0.20, 20 | Trimming thresholds |
| `beta_rank`, `plot_rank`, `divnet_rank` | Genus, Phylum, Genus | Taxonomic ranks for beta diversity, plots and DivNet |
| `kegg_plot_level` | 2 | KEGG hierarchy level shown in the function plot |
| `permutations` | 999 | Permutations for all permutation tests |
| `permanova_formula` | *(empty)* | Overrides the automatic global PERMANOVA model |
| `divnet_tuning` | `fast` | Use `careful` for final results |
| `divnet_functions` | `true` | Also estimate KO-level alpha diversity |
| `cores` | 0 | Workers; 0 = all cores but one |
| `use_cache` | `true` | Reuse expensive results while their inputs are unchanged |
| `run_analysis` | `true` | `false` renders the annotated code only (no data needed) |
| `make_plots` | `true` | `false` runs all statistics without drawing figures |

## Outputs (in `out_dir`)

* `tables/`: CSV files with trimming summaries, trimmed counts, CLR values, PCA scores, all PERMANOVA, ANOSIM and PERMDISP results, relative abundances, DivNet estimates, and `betta` term tests, pairwise differences and marginal means.
* `figures/`: Aitchison PCA (taxa, KOs), taxa and KEGG bar plots, DivNet alpha diversity (PDF and 300 dpi PNG, transparent background).
* `cache/`: cached CLR, permutation and DivNet results (`.rds`).
