#!/bin/zsh
# Run the pre-registered confirmation end to end, unattended:
#   chain_confirm.sh <expected sha256 prefix of analyse_confirm.py>
# 1. keep the machine awake for exactly as long as the launcher runs
# 2. launch all runs, wait for them
# 3. verify the analysis script is the committed one, then run it ONCE
set -u
EXPECT=$1
R=/Users/ivanosman/raevol_exp
ROOT=$R/confirm
echo "chain confirm: start $(date '+%F %T')"
$R/launch_confirm.sh $ROOT 10 > $R/launch_confirm.log 2>&1 &
LPID=$!
caffeinate -i -s -w $LPID &
wait $LPID
echo "chain confirm: runs finished $(date '+%F %T')"
GOT=$(shasum -a 256 $R/analyse_confirm.py | cut -c1-16)
if [ "$GOT" != "$EXPECT" ]; then
  echo "chain confirm: CHECKSUM MISMATCH ($GOT != $EXPECT); analysis NOT run"
  exit 1
fi
echo "chain confirm: checksum verified ($GOT); running the pre-registered analysis once"
/opt/homebrew/bin/python3 $R/analyse_confirm.py $ROOT > $R/results_confirm_preregistered.txt 2>&1
echo "EXIT=$?" >> $R/results_confirm_preregistered.txt
echo "chain confirm: done $(date '+%F %T')"
