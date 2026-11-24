#!/usr/bin/env Rscript
# =============================================================================
#  simulate_test_data.R
#
#  Generates a small, fully synthetic data set in exactly the formats written
#  by the pipeline scripts (../workflow/*_profile_sample.sh), so that the R
#  analyses in ../analysis/ can be tested without sequencing data or MEGAN.
#
#  Design (72 samples, shared by both arms):
#    2 countries x 2 locations per country x 2 compartments x 3 treatments x
#    3 replicates. Built-in effects: treatment (on ~30 % of taxa / KOs),
#    compartment (all features) and location (taxa only). There is NO
#    country effect beyond location-to-location variation, which makes the
#    data useful for checking that the location-level country test does not
#    report pseudo-replicated significance. Sample S05 is sequenced ~100x
#    more shallowly than the others to exercise the low-depth QC flag.
#
#  Files per sample (in <out>/megan_output/<sample>_MTX/ and ..._MGX/):
#    <sample>.Taxonomy.paths.counts.tsv  "[D] Bacteria; [P] ...;<TAB>count"
#    <sample>.KEGG.counts.tsv            "K00001 name<TAB>count"
#    <sample>.KEGG.paths.counts.tsv      "KEGG;L1;L2;L3;K00001 name;<TAB>count"
#    <sample>.read_counts.tsv            SeqKit "stats -T" table
#  plus <out>/metadata.csv.
#
#  Usage (from the repository root):
#    Rscript test/simulate_test_data.R test/data
#
#  Only base R is required. Samples are written in parallel (forked workers)
#  on Linux/macOS and serially on Windows.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
out  <- if (length(args)) args[1] else file.path("test", "data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Reproducible random numbers, also across parallel workers (L'Ecuyer streams).
RNGkind("L'Ecuyer-CMRG")
set.seed(1)
n_cores <- max(1L, parallel::detectCores() - 1L)

# ---- Sample design ------------------------------------------------------------
meta <- expand.grid(Rep = 1:3, Treatment = c("Control", "T1", "T2"),
                    Compartment = c("Bulk", "Rhizosphere"), Loc = 1:2,
                    Country = c("UK", "Brazil"), stringsAsFactors = FALSE)
# Location names contain "-" on purpose: compact letter displays split on "-",
# and the analysis must cope with it.
meta$Location <- paste0(meta$Country, "-L", meta$Loc)
meta$SampleID <- sprintf("S%02d", seq_len(nrow(meta)))
write.csv(meta[, c("SampleID", "Treatment", "Compartment", "Country", "Location")],
          file.path(out, "metadata.csv"), row.names = FALSE)

# ---- Taxonomy: 10 phyla x 2 classes x 3 genera, plus higher-rank nodes --------
# Each element is a MEGAN "--paths --majorRanksOnly" lineage. Reads can stop
# at phylum level, at domain level, or at the root (empty lineage). One class
# is "unknown" (rank missing in NCBI), as MEGAN writes it.
tax_paths <- character(0)
for (p in c(paste0("Bacteria|Phylum", 1:8), paste0("Archaea|Phylum", 9:10))) {
  dom <- sub("\\|.*", "", p); ph <- sub(".*\\|", "", p)
  tax_paths <- c(tax_paths, sprintf("[D] %s; [P] %s;", dom, ph))
  for (cc in 1:2) {
    cl  <- if (ph == "Phylum3" && cc == 2) "unknown" else paste0(ph, "_C", cc)
    fam <- paste0(ph, "_C", cc, "_F")
    tax_paths <- c(tax_paths, sprintf(
      "[D] %s; [P] %s; [C] %s; [O] %s_O; [F] %s; [G] %s_G%d;", dom, ph, cl, fam, fam, fam, 1:3))
  }
}
tax_paths <- c(tax_paths, "[D] Bacteria;", "")
nt <- length(tax_paths)

# Baseline abundances and effects (log scale).
base_t  <- rlnorm(nt, 3, 1.3)
eff_trt <- matrix(rnorm(nt * 3, 0, 0.5), nt, 3) * (runif(nt) < 0.3)
eff_cmp <- matrix(rnorm(nt * 2, 0, 0.6), nt, 2)
eff_loc <- matrix(rnorm(nt * 4, 0, 0.4), nt, 4)

# ---- KEGG: 400 KOs in a small three-level hierarchy ---------------------------
hier <- data.frame(
  L1 = c(rep("Metabolism", 6), rep("Genetic Information Processing", 2),
         "Environmental Information Processing", "Cellular Processes"),
  L2 = c("Carbohydrate metabolism", "Carbohydrate metabolism", "Energy metabolism",
         "Energy metabolism", "Amino acid metabolism", "Lipid metabolism",
         "Translation", "Replication and repair", "Membrane transport", "Cell motility"),
  L3 = c("Glycolysis / Gluconeogenesis [PATH:ko00010]", "Citrate cycle (TCA cycle) [PATH:ko00020]",
         "Nitrogen metabolism [PATH:ko00910]", "Methane metabolism [PATH:ko00680]",
         "Arginine biosynthesis [PATH:ko00220]", "Fatty acid degradation [PATH:ko00071]",
         "Ribosome [PATH:ko03010]", "DNA replication [PATH:ko03030]",
         "ABC transporters [PATH:ko02010]", "Flagellar assembly [PATH:ko02040]"))
nk     <- 400
kos    <- sprintf("K%05d", sort(sample(1:25000, nk)))
labels <- sprintf("%s gene%d; hypothetical enzyme %d [EC:1.1.1.%d]", kos, 1:nk, 1:nk, 1:nk)
plab   <- gsub(";", "_", labels)          # MEGAN replaces ";" inside path elements

# KOs sit in 1-3 pathways; ~10 % only in "Not Included in Pathway or Brite";
# ~20 % additionally in a BRITE hierarchy (exercises the category weighting).
kpath_lines <- unlist(lapply(seq_len(nk), function(i) {
  r <- runif(1)
  if (r < 0.1) {
    return(sprintf("KEGG;Not Included in Pathway or Brite;Unclassified: metabolism;Poorly characterized;%s;", plab[i]))
  }
  idx <- sample(nrow(hier), sample(1:3, 1, prob = c(0.6, 0.3, 0.1)))
  p <- sprintf("KEGG;%s;%s;%s;%s;", hier$L1[idx], hier$L2[idx], hier$L3[idx], plab[i])
  if (r > 0.8) p <- c(p, sprintf("KEGG;Brite Hierarchies;Protein families: metabolism;Enzymes [BR:ko01000];%s;", plab[i]))
  p
}))
kpath_ko <- regmatches(kpath_lines, regexpr("K[0-9]{5}", kpath_lines))

base_k   <- rlnorm(nk, 3, 1.2)
keff_trt <- matrix(rnorm(nk * 3, 0, 0.4), nk, 3) * (runif(nk) < 0.3)
keff_cmp <- matrix(rnorm(nk * 2, 0, 0.4), nk, 2)

# Sequencing depth per sample (identical for both arms); S05 is shallow.
depth <- rlnorm(nrow(meta), log(40000), 0.3)
depth[5] <- 400

# ---- Writer for one sample of one arm -----------------------------------------
fmt <- function(lab, cnt) paste0(lab, "\t", sprintf("%.1f", cnt))

write_sample <- function(job) {
  arm <- job$arm; i <- job$i
  s <- meta$SampleID[i]
  d <- file.path(out, "megan_output", paste0(s, "_", arm))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  ti <- match(meta$Treatment[i], c("Control", "T1", "T2"))
  ci <- match(meta$Compartment[i], c("Bulk", "Rhizosphere"))
  li <- match(meta$Location[i], unique(meta$Location))

  # Taxonomy counts: negative binomial around the expected composition.
  mu  <- base_t * exp(eff_trt[, ti] + eff_cmp[, ci] + eff_loc[, li])
  cnt <- rnbinom(nt, mu = depth[i] * mu / sum(mu), size = 2)
  writeLines(fmt(tax_paths[cnt > 0], cnt[cnt > 0]),
             file.path(d, paste0(s, ".Taxonomy.paths.counts.tsv")))

  # KO counts; the paths file repeats a KO's count on each of its paths.
  muk <- base_k * exp(keff_trt[, ti] + keff_cmp[, ci])
  kc  <- rnbinom(nk, mu = 0.6 * depth[i] * muk / sum(muk), size = 3)
  writeLines(fmt(labels[kc > 0], kc[kc > 0]), file.path(d, paste0(s, ".KEGG.counts.tsv")))
  kpc <- kc[match(kpath_ko, kos)]
  writeLines(fmt(kpath_lines[kpc > 0], kpc[kpc > 0]),
             file.path(d, paste0(s, ".KEGG.paths.counts.tsv")))

  # SeqKit read counts: raw, trimmed, and (RNA only) mRNA / rRNA fractions.
  rawp  <- round(depth[i] * 3); qcp <- round(rawp * 0.93)
  files <- c(sprintf("/data/%s_R1_001.fastq.gz", s), sprintf("/data/%s_R2_001.fastq.gz", s),
             file.path(d, paste0(s, c(".qc_R1.fq.gz", ".qc_R2.fq.gz"))))
  n <- c(rawp, rawp, qcp, qcp)
  if (arm == "MTX") {
    mr <- round(qcp * 0.25)
    files <- c(files, file.path(d, paste0(s, c(".mrna_R1.fq.gz", ".mrna_R2.fq.gz",
                                               ".rrna_R1.fq.gz", ".rrna_R2.fq.gz"))))
    n <- c(n, mr, mr, qcp - mr, qcp - mr)
  }
  rc <- data.frame(file = files, format = "FASTQ", type = "DNA", num_seqs = n,
                   sum_len = n * 150, min_len = 75, avg_len = 148, max_len = 150)
  write.table(rc, file.path(d, paste0(s, ".read_counts.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(s)
}

# ---- Write all samples (parallel where forking is available) ------------------
jobs <- c(lapply(seq_len(nrow(meta)), function(i) list(arm = "MTX", i = i)),
          lapply(seq_len(nrow(meta)), function(i) list(arm = "MGX", i = i)))
if (.Platform$OS.type == "unix" && n_cores > 1) {
  res <- parallel::mclapply(jobs, write_sample, mc.cores = n_cores, mc.set.seed = TRUE)
  if (any(vapply(res, inherits, logical(1), what = "try-error"))) stop("A worker failed.")
} else {
  invisible(lapply(jobs, write_sample))
}
cat("Simulated", nrow(meta), "samples x 2 arms in", normalizePath(out), "\n")
