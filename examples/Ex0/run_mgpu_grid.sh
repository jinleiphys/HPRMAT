#!/bin/bash
# Multi-GPU (cusolverMg, FP64) capacity/scaling grid on the 4x RTX 3090 node.
# For each (N, #GPUs) that fits in aggregate device memory, records wall time,
# accuracy vs the manufactured solution, and the peak per-GPU memory (nvidia-smi).
# Predicted-OOM cells are skipped so the Fortran wrapper never drops to the slow
# host ZGESV fallback.
set -u
cd "$(dirname "$0")/../.."
BIN=./examples/Ex0/benchmark_mgpu
NCH=${NCH:-64}
REFINE=${REFINE:-2}
OUT=${OUT:-mgpu_grid_results.txt}
MEMCAP_BYTES=$((23*1000*1000*1000))   # ~23 GB usable per 24 GB card

printf "# N\tnGPU\ttime_s\taccuracy\tpeakMB/GPU\tstatus\n" | tee "$OUT"

# Mixed-precision multi-GPU: FP32 factors -> 8 bytes/element on the GPUs, so a single
# 24 GB card reaches N ~ 51200 and the multi-GPU runs push beyond.
runs=(
  "12800:0" "12800:0,1" "12800:0,1,2,3"
  "25600:0" "25600:0,1" "25600:0,1,2,3"
  "38400:0" "38400:0,1" "38400:0,1,2,3"
  "51200:0" "51200:0,1" "51200:0,1,2,3"
  "64000:0,1" "64000:0,1,2,3"
  "76800:0,1,2,3"
)

for r in "${runs[@]}"; do
  N="${r%%:*}"; VIS="${r##*:}"
  NG=$(echo "$VIS" | tr ',' ' ' | wc -w)
  # predicted per-GPU bytes: FP32 A panel n^2/NG + workspace n*256, 8 B/element
  perGPU=$(python3 -c "print(int((${N}*${N}/${NG} + ${N}*256)*8))")
  if [ "$perGPU" -gt "$MEMCAP_BYTES" ]; then
    printf "%d\t%d\t%s\t%s\t%s\t%s\n" "$N" "$NG" "-" "-" "-" "SKIP(>mem)" | tee -a "$OUT"
    continue
  fi
  # background peak-memory sampler over the visible physical GPUs
  peakfile=$(mktemp)
  echo 0 > "$peakfile"
  ( while true; do
      m=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i "$VIS" 2>/dev/null | sort -rn | head -1)
      cur=$(cat "$peakfile")
      if [ -n "$m" ] && [ "$m" -gt "$cur" ]; then echo "$m" > "$peakfile"; fi
      sleep 0.3
    done ) &
  SAMP=$!
  # OpenBLAS threads matter: the FP64 refinement residual (zgemm_) runs on the host.
  res=$(CUDA_VISIBLE_DEVICES="$VIS" OPENBLAS_NUM_THREADS=32 OMP_NUM_THREADS=1 \
        "$BIN" "$N" "$NCH" 1 "$REFINE" 2>&1)
  kill "$SAMP" 2>/dev/null; wait "$SAMP" 2>/dev/null
  peak=$(cat "$peakfile"); rm -f "$peakfile"
  t=$(echo "$res"  | grep "solve time" | grep -oE "[0-9]+\.[0-9]+" | head -1)
  acc=$(echo "$res"| grep "max|X"      | awk '{print $NF}')
  info=$(echo "$res"| grep -E "info " | grep -oE "[-0-9]+" | tail -1)
  status=OK
  echo "$res" | grep -q "Failed to allocate" && status=OOM
  [ "${info:-0}" != "0" ] && status="FAIL(info=${info})"
  printf "%d\t%d\t%s\t%s\t%s\t%s\n" "$N" "$NG" "${t:-NA}" "${acc:-NA}" "${peak:-NA}" "$status" | tee -a "$OUT"
done
echo "=== grid complete -> $OUT ==="
