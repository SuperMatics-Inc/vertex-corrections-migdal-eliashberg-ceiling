#!/bin/bash
# Run inside a Batch task; the controller supplies BATCHRUN_OUTPUT_DIR.
set -eu
printf 'runner HOME=%s\n' "${HOME-unset}"
export HOME=/root OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
printf 'HOME=%s\n' "$HOME"
src="$PWD"
cd "$BATCHRUN_OUTPUT_DIR"
curl -LsSf https://astral.sh/uv/install.sh | sh
uv=/root/.local/bin/uv
if [ "${1:-}" = smoke ]; then
  exec "$uv" run --no-project --with numpy==2.2.6 --with scipy==1.15.3 python -B "$src/lattice.py" 66:8:2:0.2:1
fi
case "$BATCH_TASK_INDEX" in
  0) l=0.5; o=0.05;; 1) l=0.5; o=0.1;; 2) l=0.5; o=0.2;;
  3) l=1; o=0.05;; 4) l=1; o=0.1;; 5) l=1; o=0.2;;
  6) l=2; o=0.05;; 7) l=2; o=0.1;; 8) l=2; o=0.2;;
  *) exit 1;;
esac
failed=0
for grid in 250:16 510:16 510:24 1014:24 1014:32; do
  "$uv" run --no-project --with numpy==2.2.6 --with scipy==1.15.3 python -B "$src/lattice.py" "$grid:$l:$o:1" || failed=1
done
exit "$failed"
