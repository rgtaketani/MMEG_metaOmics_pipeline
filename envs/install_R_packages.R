#!/usr/bin/env Rscript
# =============================================================================
#  install_R_packages.R
#
#  Installs the R packages used by the analyses in ../analysis/ (R >= 4.2).
#  Works inside an renv project as well as in a user library.
#
#  Usage (from the repository root):
#    Rscript envs/install_R_packages.R
#
#  Sources
#    CRAN:          data.table, dplyr, tidyr, tibble, ggplot2, cowplot, scales,
#                   rlang, vegan, permute, zCompositions, multcompView,
#                   breakaway, knitr, rmarkdown, remotes, BiocManager
#    Bioconductor:  phyloseq
#    GitHub:        DivNet (adw96/DivNet; not on CRAN)
#  PDF rendering additionally needs a LaTeX installation with xelatex and the
#  fvextra package (e.g. tinytex::install_tinytex()).
# =============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install packages in parallel where supported (Ncpus = cores minus one).
ncpus <- max(1L, parallel::detectCores() - 1L)

cran <- c("data.table", "dplyr", "tidyr", "tibble", "ggplot2", "cowplot",
          "scales", "rlang", "vegan", "permute", "zCompositions",
          "multcompView", "breakaway", "knitr", "rmarkdown", "remotes",
          "BiocManager")
missing_cran <- setdiff(cran, rownames(installed.packages()))
if (length(missing_cran)) install.packages(missing_cran, Ncpus = ncpus)

if (!requireNamespace("phyloseq", quietly = TRUE)) {
  BiocManager::install("phyloseq", update = FALSE, ask = FALSE, Ncpus = ncpus)
}

if (!requireNamespace("DivNet", quietly = TRUE)) {
  remotes::install_github("adw96/DivNet", upgrade = "never", dependencies = TRUE)
}

# Report the installed versions (useful for the methods section).
pk <- c(cran, "phyloseq", "DivNet")
print(data.frame(package = pk,
                 version = vapply(pk, function(p) as.character(packageVersion(p)), "")),
      row.names = FALSE)
