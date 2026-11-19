#!/bin/bash
# 사용: ~/sk_sports/01_스크립트/rslot.sh Rscript 81_fit.R ...   (빈 슬롯이 날 때까지 기다렸다 실행)
N=${SK_RSLOTS:-7}; D=/tmp/sk_rslots; mkdir -p "$D"
while :; do
  for i in $(seq 1 "$N"); do
    if mkdir "$D/s$i" 2>/dev/null; then
      echo $$ > "$D/s$i/pid"; trap 'rm -rf "$D/s'$i'"' EXIT INT TERM
      nice -n 10 "$@"; exit $?
    fi
    # 죽은 프로세스가 잡고 있는 슬롯 회수
    p=$(cat "$D/s$i/pid" 2>/dev/null); [ -n "$p" ] && ! kill -0 "$p" 2>/dev/null && rm -rf "$D/s$i"
  done
  sleep 5
done
