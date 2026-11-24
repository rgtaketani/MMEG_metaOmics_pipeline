# Changelog

[← Back to the main README](README.md)

## 1.0.0 — 2026-11-23

First public release.

* Per-sample pipelines for metatranscriptomes and metagenomes
  ([`workflow/`](workflow/)): fastp, RiboDetector (RNA only), DIAMOND blastx
  and MEGAN Ultimate `daa-meganizer`, with identical shared parameters and a
  logged parameter fingerprint.
* R Markdown analyses ([`analysis/`](analysis/)): trimming, CLR, Aitchison
  PCA, PERMANOVA, ANOSIM, PERMDISP, taxa and KEGG composition plots, DivNet
  alpha diversity with `betta_random`.
* Methods text and workflow figure ([`docs/`](docs/)), simulated data and
  tests ([`test/`](test/)).
