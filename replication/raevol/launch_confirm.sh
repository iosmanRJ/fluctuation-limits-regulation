#!/bin/zsh
# Launch the pre-registered confirmation (PAPER.md section 11.19) with the
# INSTRUMENTED build (instrument_signal_links.py) for both aevol_create and
# aevol_run, seed-major so that stopping early still leaves a balanced design.
#   launch_confirm.sh <root> <maxpar>
# Concurrency is counted with `pgrep -x aevol_run` (the binary only).
set -u
ROOT=$1; MAXPAR=$2; NGEN=30000
B=/Users/ivanosman/raevol7_instr/build/bin
echo "launcher confirm: $(date '+%F %T')  root=$ROOT  par=$MAXPAR  ngen=$NGEN"
for s in $(seq 1 12); do
  for d in $(ls -d $ROOT/primary/*_s$s 2>/dev/null | sort); do
    while [ $(pgrep -x aevol_run | wc -l) -ge $MAXPAR ]; do sleep 20; done
    (
      cd $d
      rm -f DONE_* 2>/dev/null
      $B/aevol_create -f param.in > create.log 2>&1 || { echo "CREATE_FAILED" > DONE_create; exit; }
      $B/aevol_run -n $NGEN > run.log 2> run.err
      r=$?
      echo "exit=$r gen=$(tail -1 stats/stats_simd_best.csv | cut -d, -f1)" > DONE_$r
    ) &
    echo "  started $(basename $d) at $(date '+%T')"
    sleep 3
  done
done
wait
echo "launcher confirm done: $(date '+%F %T')"
