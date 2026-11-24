#!/usr/bin/env bash
# =============================================================================
#  run_analysis_test.sh
#
#  End-to-end test of both R Markdown analyses on simulated data:
#    1. generates the synthetic data set (test/simulate_test_data.R);
#    2. renders analysis/analysis_MTX.Rmd and analysis/analysis_MGX.Rmd to
#       HTML with quick settings (199 permutations, DivNet tuning "test").
#
#  Requires R with the packages listed in ../envs/install_R_packages.R.
#  Takes a few minutes on a laptop; uses all cores but one.
#
#  Usage (from the repository root):
#    bash test/run_analysis_test.sh
#  Outputs: test/data/ (simulated inputs) and test/results_{MTX,MGX}/ plus
#  the rendered HTML reports in test/. All are ignored by git.
# =============================================================================
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/.." && pwd)

# ---- 1. Simulated inputs ------------------------------------------------------
Rscript "$HERE/simulate_test_data.R" "$HERE/data"

# ---- 2. Render both analyses ----------------------------------------------------
# rmarkdown resolves relative paths from the Rmd's folder (analysis/), so all
# paths are passed as absolute paths. Reports are written to test/.
for arm in MTX MGX; do
    echo "Rendering analysis_${arm}.Rmd ..."
    Rscript -e "rmarkdown::render(
        '$REPO/analysis/analysis_${arm}.Rmd',
        output_format = 'html_document',
        output_dir    = '$HERE',
        params = list(input_dir     = '$HERE/data/megan_output',
                      metadata      = '$HERE/data/metadata.csv',
                      out_dir       = '$HERE/results_${arm}',
                      permutations  = 199,
                      divnet_tuning = 'test'),
        quiet = TRUE)"
    echo "PASS: analysis_${arm}.Rmd rendered -> test/analysis_${arm}.html"
done
