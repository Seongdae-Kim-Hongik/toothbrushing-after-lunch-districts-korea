#!/bin/bash
S=${ST_ROOT}/01_스크립트; cd $S; T=BRSH; L=${ST_ROOT}/04_분석/$T
export ST_GRAPH=$L/graph_593.graph ST_THREADS=1 ST_VER=std3 ST_FBDZIP=${ST_FBDZIP}
ST_FH_ONLY=S7 ST_MV_STOP=1 ./rslot.sh Rscript 81fh_brsh.R $T > $L/run_fh_s7spec.log 2>&1
( ST_FH_ONLY=S0,S2,S3 ./rslot.sh Rscript 81fh_brsh.R $T > $L/run_fh_main.log 2>&1
  ST_FH_ONLY=S4,S5,S6,S8,S9,S10 ./rslot.sh Rscript 81fh_brsh.R $T > $L/run_fh_sensA.log 2>&1
  ST_FH_ONLY=S11,S12,S13 ./rslot.sh Rscript 81fh_brsh.R $T >> $L/run_fh_sensA.log 2>&1 ) &
( ST_FH_ONLY=S7 ST_MV_RANGE=1-100 ./rslot.sh Rscript 81fh_brsh.R $T > $L/run_fh_s7_1-100.log 2>&1
  ST_FH_ONLY=S7 ST_MV_RANGE=101-200 ./rslot.sh Rscript 81fh_brsh.R $T > $L/run_fh_s7_101-200.log 2>&1
  ST_FH_ONLY=S14,S16,S18,S19,S20 ./rslot.sh Rscript 81fh_brsh.R $T > $L/run_fh_sensB.log 2>&1 ) &
wait
ST_FH_ONLY=S7,S17,S44 ./rslot.sh Rscript 81fh_brsh.R $T >> $L/run_fh_main.log 2>&1
echo "chain done $(date)" >> $L/run_fh_main.log; echo ALLDONE
