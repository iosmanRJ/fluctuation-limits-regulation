#!/bin/zsh
# Continue existing R-Aevol runs from a saved population with the instrumented
# binary (instrument_signal_links.py), which appends the influences of signal
# proteins to every row of stats/stats_simd_best.csv. Each run is COPIED to
# $OUT/<tier>/<name> first; the original run directories are never modified.
#
#   decomp_launch.sh <out_root> <maxpar> <begin_generation> <end_generation> <run_dir>...
#
# Resuming is NOT bit-reproducible in this aevol build (tested 2026-09-28: a run
# resumed at 20,000 differed from the original from generation 20,001), so the
# continuation is a new stochastic path from the saved population. We continue
# from the final population (generation 30,000) and validate by comparing the
# continuation's cell means of S with the original analysis window.
set -u
setopt NULL_GLOB   # an unmatched glob must not abort the cleanup
OUT=$1; MAXPAR=$2; BEGIN=$3; END=$4; shift 4
B=/Users/ivanosman/raevol7_instr/build/bin
mkdir -p $OUT
echo "decomp launcher $(date '+%F %T') begin=$BEGIN end=$END par=$MAXPAR runs=$#"
for src in "$@"; do
  while [ $(pgrep -x aevol_run | wc -l) -ge $MAXPAR ]; do sleep 20; done
  tier=$(basename $(dirname $src)); name=$(basename $src)
  dst=$OUT/$tier/$name
  mkdir -p $OUT/$tier
  rm -rf $dst && cp -R $src $dst
  (
    cd $dst
    for w in world/world_*.ae; do [ "$w" = "world/world_$(printf %09d $BEGIN).ae" ] || rm -f "$w"; done
    rm -f DONE_*
    $B/aevol_run -b $BEGIN -e $END > run.log 2> run.err
    echo "exit=$? gen=$(tail -1 stats/stats_simd_best.csv | cut -d, -f1)" > DONE_decomp
  ) &
  echo "  started $tier/$name at $(date '+%T')"
  sleep 2
done
wait
echo "decomp launcher done $(date '+%F %T')"
