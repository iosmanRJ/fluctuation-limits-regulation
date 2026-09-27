#!/bin/zsh
# Runs P5 (PAPER.md 11.5) and the grain experiment G (11.8) end to end on the
# Mac mini, so nothing depends on an ssh session staying up. Launch under
# caffeinate (the Mac mini idle-sleeps after one minute).
#
#   1. wait for the P5 viability pilot to finish
#   2. apply the pre-registered rule: the smallest grid of 80/96/112 at which
#      6/6 fluc_none AND 6/6 froz_none pilot runs survive 150k ticks
#   3. run the P5 confirmatory batch at that grid; analyse ONCE
#   4. generate and run the grain sweep; analyse ONCE
set -u
E=/Users/ivanosman/evosim_exp
X=/Users/ivanosman/raevol_exp
J=/opt/homebrew/bin/julia
log () { echo "chain: $* ($(date '+%F %T'))"; }

log "waiting for P5 pilot"
until grep -q ALL_DONE $E/p5_pilot.log 2>/dev/null; do sleep 60; done
log "pilot done"

CHOSEN=""
for g in 80 96 112; do
  line="grid $g:"
  ok=1
  for a in fluc_none froz_none; do
    n=0; s=0
    for f in $E/p5_pilot/g$g/${a}_s*.DONE(N); do
      n=$((n+1))
      grep -q extinct ${f%.DONE}.err 2>/dev/null || s=$((s+1))
    done
    line="$line $a $s/$n"
    [ $n -eq 6 ] && [ $s -eq 6 ] || ok=0
  done
  echo "chain: pilot $line"
  if [ $ok -eq 1 ] && [ -z "$CHOSEN" ]; then CHOSEN=$g; fi
done

if [ -z "$CHOSEN" ]; then
  log "P5: NO GRID QUALIFIES under the pre-registered rule -- P5 not run"
else
  log "P5: rule selects grid $CHOSEN"
  zsh $E/p5_launch.sh confirm $E/p5_confirm $CHOSEN > $E/p5_confirm.log 2>&1
  log "P5 confirm runs done"
  (cd /Users/ivanosman/evosim && $J --project=. $E/p5_analyse.jl $E/p5_confirm > $E/p5_results.txt 2>&1)
  log "P5 analysed once -> $E/p5_results.txt"
fi

cd $X
if [ -e grain ]; then log "grain/ exists -- refusing to overwrite"; exit 1; fi
/opt/homebrew/bin/python3 gen_params_grain.py grain/primary || { log "grain generator failed"; exit 1; }
log "grain sweep launching"
zsh launch_v2.sh $X/grain 10 > $X/grain_launcher.log 2>&1
log "grain runs done"
/opt/homebrew/bin/python3 analyse_grain.py grain/primary > grain_results.txt 2>&1
log "grain analysed once -> $X/grain_results.txt"
log "ALL_DONE"
