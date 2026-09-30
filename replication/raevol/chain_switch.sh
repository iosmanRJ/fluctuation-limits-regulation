#!/bin/zsh
# Wait for the switch launcher to finish, verify the analysis script is the
# committed one, then run the pre-registered analysis ONCE.
#   chain_switch.sh <expected sha256 prefix of analyse_switch.py>
set -u
EXPECT=$1
R=/Users/ivanosman/raevol_exp
echo "chain switch: waiting $(date '+%F %T')"
while pgrep -f launch_switch.sh > /dev/null; do sleep 60; done
echo "chain switch: runs finished $(date '+%F %T')"
GOT=$(shasum -a 256 $R/analyse_switch.py | cut -c1-16)
if [ "$GOT" != "$EXPECT" ]; then
  echo "chain switch: CHECKSUM MISMATCH ($GOT != $EXPECT); analysis NOT run"; exit 1
fi
echo "chain switch: checksum verified ($GOT); running the pre-registered analysis once"
/opt/homebrew/bin/python3 $R/analyse_switch.py $R/switch > $R/results_switch_preregistered.txt 2>&1
echo "EXIT=$?" >> $R/results_switch_preregistered.txt
echo "chain switch: done $(date '+%F %T')"
