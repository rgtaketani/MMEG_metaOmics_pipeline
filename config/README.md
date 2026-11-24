# Configuration

[← Back to the main README](../README.md)

## Sample metadata

[`metadata_template.csv`](metadata_template.csv) shows the sample sheet
expected by the [R analyses](../analysis/). The file has one row per sample
and can be CSV or TSV.

| Column | Content |
|---|---|
| `SampleID` | Sample name as produced by the pipeline: the R1 file name without the read tag and extension (`Plot01_S1_R1_001.fastq.gz` → `Plot01_S1`). |
| `Treatment` | Treatment (three levels in the reference design). The first level is the reference unless `treatment_levels` is set. |
| `Compartment` | Soil compartment (e.g. bulk soil, rhizosphere). |
| `Country` | Country. |
| `Location` | Sampling location (site) within a country. If two countries share a location name, the analysis renames it to `Country_Location`. |

Column names can differ from these; set them with the `col_*` parameters in
the YAML header of each Rmd. Extra columns are kept but not used. Missing
values in any design column stop the analysis with a list of the affected
samples.

The example rows are placeholders. The same `SampleID`s can be used for the
metatranscriptome and the metagenome of a sample, because the two arms are
analysed separately.
