#!/usr/bin/env bash
# =============================================================================
#  run_pipeline_smoke_test.sh
#
#  Checks the control flow of both per-sample pipeline scripts in seconds,
#  without reference databases or real bioinformatics tools. The heavy tools
#  (fastp, RiboDetector, DIAMOND, daa-meganizer, daa2info) are replaced by the
#  stubs in ./stubs/, which only copy or echo their inputs. seqkit and pigz
#  must be installed (both are in ../envs/pipeline.yml).
#
#  What is verified
#    * both scripts run to completion on a tiny synthetic read pair;
#    * every expected output file is created;
#    * the metatranscriptome (MTX) and metagenome (MGX) runs report the SAME
#      shared-parameter fingerprint (i.e. the harmonised settings agree);
#    * the MEGAN step receives the shared parameters.
#
#  Usage (from the repository root):
#    bash test/run_pipeline_smoke_test.sh
#  Results are written to a temporary directory that is removed on exit.
# =============================================================================
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# Stubs first in PATH; real seqkit and pigz are taken from the environment.
export PATH="$HERE/stubs:$PATH"
for t in seqkit pigz; do
    command -v "$t" >/dev/null || { echo "FAIL: '$t' is required for the smoke test." >&2; exit 1; }
done

# ---- Tiny synthetic input: 5 read pairs of 150 nt ----------------------------
cd "$WORK"
for mate in 1 2; do
    for i in 1 2 3 4 5; do
        printf '@SIM:1:FC:%d %d:N:0:ACGT\n%s\n+\n%s\n' "$i" "$mate" \
            "$(printf 'ACGT%.0s' {1..37})AC" "$(printf 'I%.0s' {1..150})"
    done | gzip > "Test_S1_R${mate}_001.fastq.gz"
done

# Placeholder reference files (contents are never read by the stubs).
echo stub > nr.dmnd; echo stub > megan.mdb
export DIAMOND_DB="$WORK/nr.dmnd" MEGAN_MAPDB="$WORK/megan.mdb" THREADS=2

# ---- Run both arms -------------------------------------------------------------
status=0
for arm in mtx mgx; do
    if bash "$REPO/workflow/${arm}_profile_sample.sh" Test_S1_R1_001.fastq.gz Test_S1_R2_001.fastq.gz \
            > "${arm}.out" 2>&1; then
        echo "PASS: ${arm}_profile_sample.sh completed"
    else
        echo "FAIL: ${arm}_profile_sample.sh exited with an error; last lines:"; tail -n 15 "${arm}.out"
        exit 1
    fi
done

# ---- Expected outputs -------------------------------------------------------------
common=(.daa .Taxonomy.counts.tsv .KEGG.counts.tsv .Taxonomy.read2class.tsv.gz
        .KEGG.read2class.tsv.gz .Taxonomy.paths.counts.tsv .KEGG.paths.counts.tsv
        .read_counts.tsv .software_versions.txt .pipeline.log)
for arm in MTX MGX; do
    extra=(); [[ $arm == MTX ]] && extra=(.rrna_R1.fq.gz .rrna_R2.fq.gz .mrna_R1.fq.gz)
    for suf in "${common[@]}" "${extra[@]}"; do
        f="Test_S1_${arm}/Test_S1${suf}"
        if [[ -e "$f" ]]; then :; else echo "FAIL: missing $f"; status=1; fi
    done
done
(( status == 0 )) && echo "PASS: all expected output files present"

# ---- Harmonisation: identical fingerprints -------------------------------------------
fp() { grep -h "Shared-parameter fingerprint" "Test_S1_$1/Test_S1.pipeline.log" | tail -n 1 | awk '{print $NF}'; }
if [[ "$(fp MTX)" == "$(fp MGX)" ]]; then
    echo "PASS: shared-parameter fingerprints match ($(fp MTX))"
else
    echo "FAIL: fingerprints differ (MTX $(fp MTX), MGX $(fp MGX))"; status=1
fi

# ---- MEGAN receives the shared parameters ---------------------------------------------
for opt in "--lcaAlgorithm weighted" "--minSupport 10" "--minScore 55"; do
    if grep -q -- "$opt" "Test_S1_MTX/Test_S1.pipeline.log"; then :; else
        echo "FAIL: daa-meganizer call lacks '$opt'"; status=1
    fi
done
(( status == 0 )) && echo "PASS: daa-meganizer called with the shared parameters"

exit $status
