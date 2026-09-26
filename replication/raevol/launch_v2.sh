#!/bin/zsh
# Launch the pre-registered R-Aevol sweep (PAPER.md section 11.4).
#
#   launch_v2.sh <root> <max_parallel>
#
# PRIMARY cells first, then SECONDARY. Within each tier the order is SEED-MAJOR
# (seed 1 of every cell, then seed 2, ...), so stopping early at any point
# still leaves a balanced design.
#
# Concurrency is counted with `pgrep -x aevol_run` -- the binary only. `pgrep -f`
# also matches the polling commands' own shells and deadlocks the launcher.
#
# Every run writes DONE_<exit> with the generation reached, so a run that died
# can never be mistaken for one that finished.
set -u
ROOT=$1
MAXPAR=$2
NGEN=30000
B=/Users/ivanosman/raevol7/build/bin

echo "launcher v2: $(date '+%F %T')  root=$ROOT  par=$MAXPAR  ngen=$NGEN"

for tier in primary secondary; do
  for s in $(seq 1 12); do
    for d in $(ls -d $ROOT/$tier/*_s$s 2>/dev/null | sort); do
      while [ $(pgrep -x aevol_run | wc -l) -ge $MAXPAR ]; do
        sleep 20
      done
      (
        cd $d
        rm -f DONE_* 2>/dev/null
        $B/aevol_create -f param.in > create.log 2>&1
        c=$?
        if [ $c -ne 0 ]; then
          echo "CREATE_FAILED_$c" > DONE_create_$c
          exit
        fi
        $B/aevol_run -n $NGEN > run.log 2> run.err
        r=$?
        echo "exit=$r gen=$(tail -1 stats/stats_simd_best.csv | cut -d, -f1)" > DONE_$r
      ) &
      echo "  started $tier/$(basename $d) at $(date '+%T')"
      sleep 3
    done
  done
done

wait
echo "launcher v2: ALL_DONE $(date '+%F %T')"
