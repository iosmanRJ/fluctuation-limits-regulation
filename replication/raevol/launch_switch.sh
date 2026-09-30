#!/bin/zsh
# Launch the environment-switch experiment (PAPER.md section 11.24).
#   launch_switch.sh <root> <maxpar> [<ngen>]
# For each run in <root>/manifest.csv:
#  1. copy the source run's generation-30,000 backup into <run>/src (the source
#     directory is never written to) and extract its best genome there with
#     aevol_post_extract (plain build; the backup format is the same);
#  2. strip the FASTA header, so aevol_create -C reads the sequence itself;
#  3. aevol_create -C with the new param.in, then aevol_run for <ngen> generations,
#     both with the INSTRUMENTED build. Seed-major order across cells.
set -u
setopt NULL_GLOB
ROOT=$1; MAXPAR=$2; NGEN=${3:-10000}
B=/Users/ivanosman/raevol7_instr/build/bin
X=/Users/ivanosman/raevol7/build/bin/aevol_post_extract
echo "launcher switch: $(date '+%F %T')  root=$ROOT  par=$MAXPAR  ngen=$NGEN"
for s in $(seq 1 12); do
  for cell in FF FU UU UF; do
    run=${cell}_s$s
    d=$ROOT/primary/$run
    src=$(awk -F, -v r=$run '$1==r {print $7}' $ROOT/manifest.csv)
    while [ $(pgrep -x aevol_run | wc -l) -ge $MAXPAR ]; do sleep 20; done
    (
      cd $d
      rm -f DONE_*
      rm -rf src && mkdir -p src/world src/exp_setup
      cp $src/world/world_000030000.ae src/world/
      cp $src/exp_setup/exp_setup_const.ae $src/exp_setup/exp_setup_000030000.ae src/exp_setup/
      cp -R $src/output_profile src/
      cp $src/param.in $src/binding_matrix.rae src/
      ( cd src && $X -t 30000 -S ../genome.fa > ../extract.log 2>&1 ) || { echo EXTRACT_FAILED > DONE_extract; exit; }
      sed -n 2p genome.fa | tr -d ' \r\n' > chrom.txt; echo >> chrom.txt
      [ $(wc -c < chrom.txt) -gt 100 ] || { echo EMPTY_GENOME > DONE_extract; exit; }
      rm -rf src
      $B/aevol_create -f param.in -C chrom.txt > create.log 2>&1 || { echo CREATE_FAILED > DONE_create; exit; }
      $B/aevol_run -n $NGEN > run.log 2> run.err
      r=$?
      echo "exit=$r gen=$(tail -1 stats/stats_simd_best.csv | cut -d, -f1)" > DONE_$r
    ) &
    echo "  started $run at $(date '+%T')"
    sleep 3
  done
done
wait
echo "launcher switch done: $(date '+%F %T')"
