#!/bin/zsh
# P5 launcher (PAPER.md section 11.5).
#
#   p5_launch.sh pilot   <root>             grids 80 96 112, 6 seeds/arm (1001-1006)
#   p5_launch.sh confirm <root> <grid>      chosen grid, 30 seeds/arm (2001-2030)
#
# Settings are those of section 11.1 except grid. Seed-major order, 6 concurrent
# runs x 2 threads. Wrap the whole launcher in caffeinate: the Mac mini sleeps
# after one minute of user idle and running processes do not count as activity.
set -u
MODE=$1; ROOT=$2
J=/opt/homebrew/bin/julia
cd /Users/ivanosman/evosim
COMMON=(--p-edge-silence 0.2 --maint-edge 0.00005 --cond-inflow-frac 0.4 --cond-period 2000 --cond-strength 0.0)
typeset -A ARM
ARM[fluc_none]="--cond-states 6"
ARM[froz_none]="--cond-states 1"
PAR=6

if [ $MODE = pilot ]; then
  GRIDS=(80 96 112); SEEDS=($(seq 1001 1006))
else
  GRIDS=($3); SEEDS=($(seq 2001 2030))
fi
mkdir -p $ROOT
echo "p5 $MODE: $(date '+%F %T') grids=${GRIDS[*]} seeds=${SEEDS[1]}-${SEEDS[-1]}"

for s in $SEEDS; do
  for g in $GRIDS; do
    for a in fluc_none froz_none; do
      while [ $(pgrep -x julia | wc -l) -ge $PAR ]; do sleep 5; done
      if [ $MODE = pilot ]; then out=$ROOT/g${g}/${a}_s${s}; else out=$ROOT/${a}_s${s}; fi
      mkdir -p ${out:h}
      (
        rm -rf $out
        $J --project=. -t 2 scripts/run.jl --gens 150000 --out $out --seed $s \
           --grid $g $COMMON ${=ARM[$a]} --no-snapshots > ${out}.log 2> ${out}.err
        echo "exit=$?" > ${out}.DONE
      ) &
      sleep 1
    done
  done
done
wait
echo "p5 $MODE: ALL_DONE $(date '+%F %T')"
