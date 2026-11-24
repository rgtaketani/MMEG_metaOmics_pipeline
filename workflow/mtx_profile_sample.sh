#!/usr/bin/env bash
# =============================================================================
#  mtx_profile_sample.sh
#
#  Part of the soil-mtx-mgx-readbased repository; see ../README.md for the
#  overall workflow, ../docs/methods.md for a methods description and
#  ../analysis/ for the downstream statistical analysis in R.
#
#  Read-based taxonomic and functional (KEGG) profiling of ONE soil
#  METATRANSCRIPTOME sample sequenced as Illumina paired-end reads.
#  No assembly is performed; every read pair is classified individually.
#
#  Steps
#    1. Read QC and trimming ........................ fastp
#    2. rRNA removal (rRNA reads kept aside) ........ RiboDetector
#    3. Query construction (mates interleaved) ...... awk + pigz
#    4. Translated search vs NCBI-nr ................ DIAMOND blastx
#    5. Taxonomy (weighted LCA) + KEGG assignment ... MEGAN Ultimate (daa-meganizer)
#    6. Count tables and read-level assignments ..... MEGAN (daa2info)
#
#  Companion script: mgx_profile_sample.sh (metagenomes). The block headed
#  "SHARED PARAMETERS" MUST be identical in both scripts, otherwise RNA and
#  DNA profiles are not comparable (e.g. for RNA/DNA ratios). Each run writes
#  a parameter fingerprint to its log; RNA and DNA runs of the same project
#  must report the same fingerprint.
#
#  Usage
#    bash mtx_profile_sample.sh <reads_R1.fastq.gz> <reads_R2.fastq.gz>
#
#  The ONLY arguments are the two read files. The sample name is derived from
#  the R1 file name (e.g. Plot12_S3_R1_001.fastq.gz -> Plot12_S3).
#
#  ASSUMPTION: all tools below are installed and available in $PATH
#  (e.g. through the conda environment in ../envs/pipeline.yml, or modules):
#    fastp, ribodetector_cpu, seqkit, pigz,
#    diamond, daa-meganizer, daa2info (both ship with MEGAN;
#                                     KEGG assignment requires MEGAN Ultimate
#                                     and an "-ue" mapping database)
#
#  Reference databases are NOT arguments. Set them once below, or export the
#  corresponding environment variables before calling the script:
#    export DIAMOND_DB=/refs/nr.dmnd
#    export MEGAN_MAPDB=/refs/megan-map-Feb2022-ue.db
#    export MEGAN_SYN2TAXA=/refs/SSURef_Nr99_132_tax_silva_to_NCBI_synonyms.map.gz  # optional
#
#  Optional run-time settings (environment variables):
#    THREADS          number of CPU threads   (default: all cores, nproc)
#    OUTDIR           parent output directory (default: current directory)
#    SAMPLE           override derived sample name
#    DIAMOND_BLOCK    DIAMOND --block-size    (default 8;  larger = faster, more RAM)
#    DIAMOND_CHUNKS   DIAMOND --index-chunks  (default 1;  larger = less RAM, slower)
#    MAX_JAVA_JOBS    concurrent daa2info exports (default 2; each is a JVM)
#    KEEP_INTERMEDIATE keep intermediate FASTQ files (default 0)
#
#  Memory note: MEGAN command-line tools read their Java heap size from the
#  MEGAN.vmoptions file in the MEGAN installation (-Xmx). Increase it for
#  large samples if daa-meganizer/daa2info report OutOfMemoryError.
#
#  Re-running: each completed step leaves a ".done" marker in the sample's
#  output directory, so an interrupted run resumes at the first unfinished
#  step. Delete the marker (or the directory) to force a step to re-run.
# =============================================================================

# Abort on any error (-e), on use of an unset variable (-u), and if any
# command within a pipeline fails (pipefail). Bash >= 4.4 is assumed.
set -euo pipefail

# -----------------------------------------------------------------------------
# 0. Arguments and derived names
# -----------------------------------------------------------------------------
if [[ $# -ne 2 ]]; then
    echo "Usage: $(basename "$0") <reads_R1.fastq.gz> <reads_R2.fastq.gz>" >&2
    exit 1
fi

# Resolve absolute paths so the script can change directory safely.
R1_IN=$(readlink -f "$1")
R2_IN=$(readlink -f "$2")

for f in "$R1_IN" "$R2_IN"; do
    [[ -s "$f" ]] || { echo "ERROR: read file not found or empty: $f" >&2; exit 1; }
done
[[ "$R1_IN" != "$R2_IN" ]] || { echo "ERROR: R1 and R2 are the same file." >&2; exit 1; }

# Derive the sample name from the R1 file name: strip the FASTQ extension,
# then a trailing read tag such as _R1, _R1_001, .1 or _1.
if [[ -z "${SAMPLE:-}" ]]; then
    SAMPLE=$(basename "$R1_IN" | sed -E 's/\.(fastq|fq)(\.gz)?$//; s/[._-]R?1(_001)?$//')
fi

# -----------------------------------------------------------------------------
# 1. Reference databases (edit the defaults or export the variables)
# -----------------------------------------------------------------------------
# DIAMOND index of the NCBI-nr protein database (diamond makedb).
DIAMOND_DB=${DIAMOND_DB:-/path/to/nr.dmnd}
# MEGAN Ultimate Edition mapping database. Must be an "-ue" file (only these
# contain KEGG mappings) and should match the nr release used for DIAMOND.
MEGAN_MAPDB=${MEGAN_MAPDB:-/path/to/megan-map-Feb2022-ue.db}
# Optional SILVA synonyms-to-taxonomy file (daa-meganizer -s2t), as in the
# group's array script. It only affects hits to SILVA accessions, so it has no
# effect on searches against NCBI-nr; leave empty to omit.
MEGAN_SYN2TAXA=${MEGAN_SYN2TAXA:-}

# -----------------------------------------------------------------------------
# 2. Run-time settings
# -----------------------------------------------------------------------------
THREADS=${THREADS:-$(nproc)}                   # use every available core
FASTP_THREADS=$(( THREADS > 16 ? 16 : THREADS ))  # fastp uses at most 16 workers
DIAMOND_BLOCK=${DIAMOND_BLOCK:-8}
DIAMOND_CHUNKS=${DIAMOND_CHUNKS:-1}
MAX_JAVA_JOBS=${MAX_JAVA_JOBS:-2}
KEEP_INTERMEDIATE=${KEEP_INTERMEDIATE:-0}

OUTDIR=$(readlink -f "${OUTDIR:-$PWD}")/${SAMPLE}_MTX
TMPDIR_RUN=${OUTDIR}/tmp
mkdir -p "$OUTDIR" "$TMPDIR_RUN"
P=${OUTDIR}/${SAMPLE}                          # prefix for all output files

# =============================================================================
#  SHARED PARAMETERS — keep IDENTICAL in mtx_profile_sample.sh and
#  mgx_profile_sample.sh. These control query geometry, search sensitivity
#  and MEGAN assignment, i.e. everything that affects RNA/DNA comparability.
# =============================================================================
TRUNC_LEN=150          # truncate all reads to this length; set to the SHORTEST
                       # read length among ALL RNA and DNA runs of the project
MIN_LEN=75             # discard reads shorter than this after trimming
CUT_Q=20               # fastp sliding-window mean quality (3' end)
# --- DIAMOND blastx (as in the group's diam_megan_array.sh) ----------------
# Default sensitivity mode (no --sensitive flag), default number of target
# sequences per query (--max-target-seqs 25), and a minimum bit score of 55,
# which in DIAMOND replaces the e-value threshold.
DIAMOND_MIN_SCORE=55   # --min-score: minimum bit score to report a hit
# --- MEGAN daa-meganizer (as in the group's diam_megan_array.sh) -----------
MIN_SCORE=55           # -ms   minimum bit score
MAX_EXPECTED=0.01      # -me   maximum e-value
MIN_PCT_IDENTITY=0     # -mpi  minimum percent identity (0 = off)
TOP_PERCENT=10         # -top  hits within 10 % of the best score enter the LCA
MIN_SUPPORT=10         # -sup  a taxon needs >= 10 reads, otherwise its reads
                       #       are pushed up to the parent taxon. Setting -sup
                       #       (and not -supp) disables the percentage filter.
MIN_READ_COVER=0       # -mrc  minimum % of read covered by alignments (0 = off)
LCA_ALG=weighted       # -alg  weighted LCA
LCA_COVER_PCT=80       # -lcp  weighted LCA must cover 80 % of the hit weight
READ_ASSIGN=readCount  # -ram  each read counts as one
CLASSIFICATIONS=(Taxonomy KEGG)   # add EGGNOG, EC, INTERPRO2GO, SEED if needed
# =============================================================================

# -----------------------------------------------------------------------------
# Logging: everything printed from here on is also written to the log file.
# -----------------------------------------------------------------------------
LOG=${P}.pipeline.log
exec > >(tee -a "$LOG") 2>&1

msg() { printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*"; }

# Step checkpointing: run_step <name> <function>. A completed step leaves
# <prefix>.<name>.done; existing markers cause the step to be skipped.
run_step() {
    local name=$1; shift
    if [[ -f ${P}.${name}.done ]]; then
        msg "Skipping step '${name}' (already completed)."
        return 0
    fi
    msg "Starting step '${name}'."
    "$@"
    touch "${P}.${name}.done"
    msg "Finished step '${name}'."
}

# -----------------------------------------------------------------------------
# Pre-flight checks: tools in PATH and databases present
# -----------------------------------------------------------------------------
REQUIRED_TOOLS=(fastp ribodetector_cpu seqkit pigz diamond daa-meganizer daa2info)
for t in "${REQUIRED_TOOLS[@]}"; do
    command -v "$t" >/dev/null 2>&1 || { msg "ERROR: '$t' not found in PATH."; exit 1; }
done
[[ -s "$DIAMOND_DB"  ]] || { msg "ERROR: DIAMOND database not found: $DIAMOND_DB"; exit 1; }
[[ -s "$MEGAN_MAPDB" ]] || { msg "ERROR: MEGAN mapping DB not found: $MEGAN_MAPDB"; exit 1; }
if [[ -n "$MEGAN_SYN2TAXA" && ! -s "$MEGAN_SYN2TAXA" ]]; then
    msg "ERROR: MEGAN synonyms file not found: $MEGAN_SYN2TAXA"; exit 1
fi

# Parameter fingerprint: a short hash of every shared setting plus the
# reference file names. RNA and DNA runs must report the same value.
PARAM_STRING="TRUNC_LEN=${TRUNC_LEN};MIN_LEN=${MIN_LEN};CUT_Q=${CUT_Q};\
DIAMOND_MIN_SCORE=${DIAMOND_MIN_SCORE};MIN_SCORE=${MIN_SCORE};MAX_EXPECTED=${MAX_EXPECTED};\
MIN_PCT_IDENTITY=${MIN_PCT_IDENTITY};TOP_PERCENT=${TOP_PERCENT};MIN_SUPPORT=${MIN_SUPPORT};\
MIN_READ_COVER=${MIN_READ_COVER};LCA_ALG=${LCA_ALG};LCA_COVER_PCT=${LCA_COVER_PCT};\
READ_ASSIGN=${READ_ASSIGN};CLASSIFICATIONS=${CLASSIFICATIONS[*]};\
DIAMOND_DB=$(basename "$DIAMOND_DB");MEGAN_MAPDB=$(basename "$MEGAN_MAPDB");\
MEGAN_SYN2TAXA=$(basename "${MEGAN_SYN2TAXA:-none}")"
FINGERPRINT=$(printf '%s' "$PARAM_STRING" | md5sum | cut -c1-12)

msg "Sample: ${SAMPLE}  (metatranscriptome)"
msg "R1: ${R1_IN}"
msg "R2: ${R2_IN}"
msg "Output directory: ${OUTDIR}"
msg "Threads: ${THREADS}"
msg "Shared-parameter fingerprint: ${FINGERPRINT}"
msg "Shared parameters: ${PARAM_STRING}"

# Record software versions for the methods section. "|| true" stops a tool
# that returns a non-zero exit code on --version from aborting the run.
{
    echo "fastp:        $(fastp --version 2>&1 | head -n1 || true)"
    echo "RiboDetector: $(ribodetector_cpu --version 2>&1 | head -n1 || true)"
    echo "SeqKit:       $(seqkit version 2>&1 | head -n1 || true)"
    echo "pigz:         $(pigz --version 2>&1 | head -n1 || true)"
    echo "DIAMOND:      $(diamond version 2>&1 | head -n1 || true)"
    echo "MEGAN:        record manually (see the MEGAN installation / About dialog)"
    echo "DIAMOND_DB:   ${DIAMOND_DB}"
    echo "MEGAN_MAPDB:  ${MEGAN_MAPDB}"
    echo "MEGAN_SYN2TAXA: ${MEGAN_SYN2TAXA:-none}"
    echo "fingerprint:  ${FINGERPRINT}"
} > "${P}.software_versions.txt"

# =============================================================================
# STEP 1 — Read QC and trimming (fastp)
# =============================================================================
# * Adapter detection for paired-end data and removal.
# * Sliding-window quality trimming from the 3' end (mean Q >= CUT_Q).
# * Poly-G tail trimming (artefact of two-colour chemistry, NovaSeq/NextSeq).
# * Truncation to TRUNC_LEN so query length is identical across all runs.
# * NO deduplication: in RNA-seq, identical reads mostly reflect genuine
#   high expression rather than PCR duplication (UMIs would be needed).
step_qc() {
    fastp \
        --in1 "$R1_IN" --in2 "$R2_IN" \
        --out1 "${P}.qc_R1.fq.gz" --out2 "${P}.qc_R2.fq.gz" \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --cut_right --cut_right_mean_quality "$CUT_Q" \
        --max_len1 "$TRUNC_LEN" --max_len2 "$TRUNC_LEN" \
        --length_required "$MIN_LEN" \
        --thread "$FASTP_THREADS" \
        --json "${P}.fastp.json" --html "${P}.fastp.html" \
        --report_title "${SAMPLE} (metatranscriptome)"
}
run_step qc step_qc

# =============================================================================
# STEP 2 — rRNA removal (RiboDetector; metatranscriptome only)
# =============================================================================
# RiboDetector classifies each read pair as rRNA or non-rRNA with a deep
# learning model. "-e norrna" makes the non-rRNA output the high-confidence
# class: ambiguous pairs are sent to the rRNA output, which is conservative
# for mRNA-based analyses. The rRNA pairs are KEPT (active-community SSU
# profiling; the rRNA fraction is also the depletion QC metric).
# -l should be the mean read length; it is computed from the input reads.
step_rrna() {
    local mean_len
    mean_len=$(seqkit stats -j "$THREADS" -T "${P}.qc_R1.fq.gz" \
               | awk 'NR==2 {printf "%d", $7}')
    msg "Mean read length passed to RiboDetector: ${mean_len}"

    # RiboDetector writes uncompressed FASTQ; outputs are compressed with
    # pigz (parallel gzip) afterwards.
    ribodetector_cpu \
        -t "$THREADS" \
        -l "$mean_len" \
        -e norrna \
        --chunk_size 256 \
        -i "${P}.qc_R1.fq.gz" "${P}.qc_R2.fq.gz" \
        -o "${TMPDIR_RUN}/mrna_R1.fq" "${TMPDIR_RUN}/mrna_R2.fq" \
        -r "${TMPDIR_RUN}/rrna_R1.fq" "${TMPDIR_RUN}/rrna_R2.fq" \
        --log "${P}.ribodetector.log"

    pigz -p "$THREADS" -c "${TMPDIR_RUN}/mrna_R1.fq" > "${P}.mrna_R1.fq.gz"
    pigz -p "$THREADS" -c "${TMPDIR_RUN}/mrna_R2.fq" > "${P}.mrna_R2.fq.gz"
    pigz -p "$THREADS" -c "${TMPDIR_RUN}/rrna_R1.fq" > "${P}.rrna_R1.fq.gz"
    pigz -p "$THREADS" -c "${TMPDIR_RUN}/rrna_R2.fq" > "${P}.rrna_R2.fq.gz"
    rm -f "${TMPDIR_RUN}"/{mrna,rrna}_R{1,2}.fq
}
run_step rrna step_rrna

# =============================================================================
# STEP 3 — Query construction (awk + pigz)
# =============================================================================
# DIAMOND takes one query file, so mates are interleaved into ONE file and
# renamed <read>/1 and <read>/2, which keeps every read name unique.
# daa-meganizer has no paired-read mode: each MATE is classified and counted
# as a separate read (see the MEGAN step). The awk program also verifies that
# R1 and R2 are in sync and aborts if mate names differ.
# Each FASTQ record is first flattened onto one line (paste - - - -), the two
# files are joined side by side, and awk writes the interleaved records.
step_query() {
    paste \
        <(pigz -dc -p 4 "${P}.mrna_R1.fq.gz" | paste - - - -) \
        <(pigz -dc -p 4 "${P}.mrna_R2.fq.gz" | paste - - - -) \
    | awk -F'\t' 'BEGIN { OFS = "\n" }
        {
            # Fields: 1-4 = R1 header/seq/plus/qual, 5-8 = R2 record.
            n1 = $1; n2 = $5
            sub(/ .*$/, "", n1);  sub(/ .*$/, "", n2)        # drop comments
            sub(/\/[12]$/, "", n1); sub(/\/[12]$/, "", n2)   # drop old /1 /2
            if (n1 != n2) {
                printf("ERROR: mate names differ at pair %d: %s vs %s\n",
                       NR, n1, n2) > "/dev/stderr"
                exit 1
            }
            print n1 "/1", $2, "+", $4, n2 "/2", $6, "+", $8
        }' \
    | pigz -p "$THREADS" > "${P}.query.fq.gz"
}
run_step query step_query

# =============================================================================
# STEP 4 — Translated search against NCBI-nr (DIAMOND blastx)
# =============================================================================
# * Both strands are searched: the query file mixes R1 (antisense to the
#   transcript in dUTP-stranded libraries) and R2 (sense), so a single
#   --strand setting cannot be applied. This also keeps the search identical
#   to the metagenome arm.
# * Parameters follow the group's array script: default sensitivity mode,
#   default --max-target-seqs (25 hits per query) and --min-score 55 (a bit
#   score threshold; when it is set DIAMOND ignores the e-value cut-off).
#   Output is DAA (--outfmt 100), which daa-meganizer annotates in place.
# * --block-size / --index-chunks trade memory for speed (see header).
step_diamond() {
    diamond blastx \
        --db "$DIAMOND_DB" \
        --query "${P}.query.fq.gz" \
        --out "${P}.daa" \
        --outfmt 100 \
        --min-score "$DIAMOND_MIN_SCORE" \
        --threads "$THREADS" \
        --block-size "$DIAMOND_BLOCK" \
        --index-chunks "$DIAMOND_CHUNKS" \
        --tmpdir "$TMPDIR_RUN"
    # DIAMOND's progress messages go to stderr and are captured in the log.
}
run_step diamond step_diamond

# =============================================================================
# STEP 5 — Taxonomic (LCA) and KEGG assignment (MEGAN Ultimate, daa-meganizer)
# =============================================================================
# Parameters follow the group's array script (diam_megan_array.sh).
# daa-meganizer writes the classifications INTO the .daa file (in place), so
# the meganized DAA is the main output and can be opened directly in MEGAN.
# * Taxonomy: weighted LCA over hits passing MIN_SCORE / MAX_EXPECTED /
#   TOP_PERCENT; the LCA must cover LCA_COVER_PCT % of the hit weight.
# * KEGG: assigned from the best-scoring hit that has a KEGG mapping
#   (MEGAN's standard rule for functional classifications; not an LCA).
# * --minSupport 10: taxa with fewer than 10 reads are merged into their parent.
#   Because --minSupportPercent is not given, MEGAN switches the percentage
#   filter off and applies the absolute threshold only.
# * No paired-read handling: daa-meganizer does not offer it, so R1 and R2 of
#   a fragment are assigned and counted independently.
# * --parseTaxonNames true and the optional synonyms file (-s2t) replicate the
#   array script; neither changes assignments of NCBI-nr accessions.
step_megan() {
    # Optional synonyms file is added only when MEGAN_SYN2TAXA is set.
    local syn_opt=()
    [[ -n "$MEGAN_SYN2TAXA" ]] && syn_opt=(--syn2taxa "$MEGAN_SYN2TAXA")

    daa-meganizer \
        --in "${P}.daa" \
        --mapDB "$MEGAN_MAPDB" \
        "${syn_opt[@]}" \
        --minScore "$MIN_SCORE" \
        --maxExpected "$MAX_EXPECTED" \
        --minPercentIdentity "$MIN_PCT_IDENTITY" \
        --topPercent "$TOP_PERCENT" \
        --minSupport "$MIN_SUPPORT" \
        --minPercentReadCover "$MIN_READ_COVER" \
        --lcaAlgorithm "$LCA_ALG" \
        --lcaCoveragePercent "$LCA_COVER_PCT" \
        --readAssignmentMode "$READ_ASSIGN" \
        --parseTaxonNames true \
        --threads "$THREADS"
}
run_step megan step_megan

# =============================================================================
# STEP 6 — Export count tables and read-level assignments (daa2info)
# =============================================================================
# For each classification:
#   <prefix>.<CLS>.counts.tsv         class -> number of reads assigned to it
#   <prefix>.<CLS>.read2class.tsv.gz  read  -> class (for taxon x KO tables)
# plus two path tables used by the R analysis:
#   <prefix>.Taxonomy.paths.counts.tsv  lineage at major ranks -> reads
#   <prefix>.KEGG.paths.counts.tsv      KEGG hierarchy path -> reads. A KO that
#                                       sits in several pathways appears on
#                                       several lines with the SAME count, so
#                                       this file is for the hierarchy only and
#                                       must not be summed directly.
# Exports run concurrently (each is a separate JVM); MAX_JAVA_JOBS caps how
# many run at once to bound memory use.
export_class2count() {      # $1 = classification
    daa2info --in "${P}.daa" --class2count "$1" --names \
        > "${P}.$1.counts.tsv"
}
export_read2class() {       # $1 = classification
    daa2info --in "${P}.daa" --read2class "$1" --names \
        | pigz -p 2 > "${P}.$1.read2class.tsv.gz"
}
export_tax_paths() {
    daa2info --in "${P}.daa" --class2count Taxonomy --paths --majorRanksOnly \
        > "${P}.Taxonomy.paths.counts.tsv"
}
export_kegg_paths() {
    daa2info --in "${P}.daa" --class2count KEGG --paths \
        > "${P}.KEGG.paths.counts.tsv"
}

step_export() {
    local pids=() jobs_list=() fail=0 j
    for CLS in "${CLASSIFICATIONS[@]}"; do
        jobs_list+=("export_class2count $CLS" "export_read2class $CLS")
    done
    jobs_list+=("export_tax_paths")
    # KEGG hierarchy paths only when KEGG is among the exported classifications.
    [[ " ${CLASSIFICATIONS[*]} " == *" KEGG "* ]] && jobs_list+=("export_kegg_paths")

    for j in "${jobs_list[@]}"; do
        # Wait while the maximum number of concurrent exports is running.
        while (( $(jobs -rp | wc -l) >= MAX_JAVA_JOBS )); do sleep 2; done
        $j &                    # word-splitting is intended: "function arg"
        pids+=($!)
    done
    # Collect every exit status so that a failed export stops the pipeline.
    for j in "${pids[@]}"; do wait "$j" || fail=1; done
    (( fail == 0 )) || { msg "ERROR: one or more daa2info exports failed."; return 1; }
}
run_step export step_export

# =============================================================================
# Read-count summary across stages (SeqKit, multithreaded)
# =============================================================================
# Reports reads surviving each stage; the rRNA fraction for this sample is
# rrna / (rrna + mrna). Only files that exist are summarised.
step_stats() {
    local files=()
    for f in "$R1_IN" "$R2_IN" \
             "${P}".qc_R{1,2}.fq.gz \
             "${P}".mrna_R{1,2}.fq.gz "${P}".rrna_R{1,2}.fq.gz; do
        [[ -e "$f" ]] && files+=("$f")
    done
    seqkit stats -j "$THREADS" -T -a "${files[@]}" > "${P}.read_counts.tsv"
}
run_step stats step_stats

# =============================================================================
# Clean-up
# =============================================================================
rm -f "${P}.query.fq.gz"                      # large, fully reproducible
if [[ "$KEEP_INTERMEDIATE" -eq 0 ]]; then
    rm -f "${P}".qc_R{1,2}.fq.gz      # mRNA and rRNA fractions are kept
fi
rm -rf "$TMPDIR_RUN"

msg "Done. Main outputs:"
msg "  ${P}.daa  (meganized: alignments + classifications; opens in MEGAN)"
msg "  ${P}.{Taxonomy,KEGG}.counts.tsv and .read2class.tsv.gz"
msg "  ${P}.read_counts.tsv, ${P}.software_versions.txt"
msg "Shared-parameter fingerprint: ${FINGERPRINT}"
