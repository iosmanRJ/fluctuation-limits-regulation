#!/bin/zsh
# EvoSim arms for the strength-vs-count test.
#   launch.sh <root> <ticks> <nseeds> <parallel>
#
# Ordered SEED-MAJOR so stopping early still leaves a balanced design rather
# than some arms complete and others empty.
#
# Concurrency is counted with `pgrep -x julia`, which matches only the binary.
# NOT `pgrep -f`, which also matches the shell wrappers of the polling
# commands and deadlocks the launcher -- and note a live viewer process would
# also be counted, so check for one before trusting the count.
set -u
ROOT=$1; TICKS=$2; NSEEDS=$3; PAR=$4
J=/opt/homebrew/bin/julia
cd /Users/ivanosman/evosim
COMMON=(--p-edge-silence 0.2 --maint-edge 0.00005 --grid 128 --cond-inflow-frac 0.4 --cond-period 2000)

typeset -A ARMS
ARMS[fluc_info]="--cond-states 6"
ARMS[fluc_yoked]="--cond-states 6 --cond-scramble true"
ARMS[fluc_none]="--cond-states 6 --cond-strength 0.0"
ARMS[froz_info]="--cond-states 1"
ARMS[froz_yoked]="--cond-states 1 --cond-scramble true"
ARMS[froz_none]="--cond-states 1 --cond-strength 0.0"

mkdir -p $ROOT
echo "launcher $(date '+%F %T')  ticks=$TICKS seeds=$NSEEDS par=$PAR"

for s in $(seq 1 $NSEEDS); do
  for a in fluc_info fluc_yoked fluc_none froz_info froz_yoked froz_none; do
    while [ $(pgrep -x julia | wc -l) -ge $PAR ]; do sleep 5; done
    out=$ROOT/${a}_s${s}
    (
      rm -rf $out
      $J --project=. -t 2 scripts/run.jl --gens $TICKS --out $out --seed $s \
         $COMMON ${=ARMS[$a]} --no-snapshots > ${out}.log 2> ${out}.err
      echo "exit=$?" > ${out}.DONE
    ) &
    sleep 1
  done
done
wait
echo "launcher ALL_DONE $(date '+%F %T')"
