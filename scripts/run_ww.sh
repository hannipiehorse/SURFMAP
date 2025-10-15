#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./scripts/run_ww.sh /home/ugfmo/surfmap_runs/pdbs/1A8E.pdb 1A8E_ww_gly0_full
#   (PDB defaults to /home/ugfmo/surfmap_runs/pdbs/1A8E.pdb, OUTDIR defaults to 1A8E_ww_gly0_full)

IMG="${IMG:-surfmap_ma:csv3d}"
RUNS="${RUNS:-/home/ugfmo/surfmap_runs}"
PDB="${1:-$RUNS/pdbs/1A8E.pdb}"
OUT="${2:-1A8E_ww_gly0_full}"

echo "Image:  $IMG"
echo "RUNS:   $RUNS"
echo "PDB:    $PDB"
echo "OUT:    $OUT"

# basic checks
ls -ld "$RUNS"
ls -l "$PDB"

# Run WW with 2D (pdf+png), 3D (png), and CSV
docker run --rm \
  -v "$RUNS":/runs \
  "$IMG" \
  python -m surfmap.bin.surfmap \
    -pdb "/runs${PDB#$RUNS}" \
    -tomap wimley_white \
    -d "/runs/$OUT" \
    --plot3d \
    --csv3d \
    --png

echo "== Outputs =="
ls -1 "$RUNS/$OUT/maps/"*wimley_white_map.*
ls -1 "$RUNS/$OUT/plots_3d/"*wimley_white_3d.png
ls -1 "$RUNS/$OUT/points_3d/"*wimley_white_3d_points.csv
