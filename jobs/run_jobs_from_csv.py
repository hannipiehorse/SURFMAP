#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Read jobs from a CSV and run SURFMAP (in Docker) for each row.
# Each run of this script writes all outputs into one batch folder:
#   /home/ugfmo/surfmap_runs/<BATCH>/<per-job-outdir>/
#
# Usage:
#   python3 -u run_jobs_from_csv.py --batch 241026_ph8_elution /path/to/jobs.csv
#
# CSV columns (header required, comma-separated):
#   protein,pdb,property,ph,salt,proj,cellsize,outdir
#     - pdb: relative to RUNS (e.g. "pdbs/1A8E.pdb") OR absolute
#     - property: wimley_white or electrostatics
#     - ph/salt: used only when property == electrostatics
#     - proj: flamsteed|mollweide|lambert (optional; default flamsteed)
#     - cellsize: e.g. 5 (optional; default 5)
#     - outdir: optional per-job name; if empty we auto-generate
#
import csv
import os
import shlex
import subprocess
import sys
import argparse
from pathlib import Path
from datetime import datetime

IMAGE = os.environ.get("SURFMAP_IMAGE", "surfmap_ma:csv3d")
RUNS = os.environ.get("RUNS", "/home/ugfmo/surfmap_runs")

DEFAULT_PROJ = "flamsteed"   # flamsteed|mollweide|lambert
DEFAULT_CELLSIZE = "5"
PLOT3D_POINT_SIZE = "3.0"    # nice default

def norm_prop(s: str) -> str:
    """Normalize property names; accept short aliases."""
    x = (s or "").strip().lower()
    if x in ("ww", "wimleywhite", "wimley-white"):
        return "wimley_white"
    if x in ("elec", "electrostatic", "electrostatic_potential"):
        return "electrostatics"
    return x

def abs_pdb_path(pdb_field: str) -> str:
    """Resolve PDB to an absolute *host* path. If relative, resolve against RUNS."""
    p = (pdb_field or "").strip()
    if not p:
        return ""
    if p.startswith("/"):
        return p
    return str(Path(RUNS) / p)

def pdb_host_to_container(pdb_abs: str) -> str:
    """
    Convert an absolute host path under RUNS to the container-visible path (/runs/...).
    This guarantees the container can read the PDB.
    """
    # Ensure pdb_abs is under RUNS
    try:
        rel = Path(pdb_abs).resolve().relative_to(Path(RUNS).resolve())
    except Exception:
        # If it's outside RUNS, we can't see it in the container; warn loudly.
        raise SystemExit(f"PDB must be inside RUNS ({RUNS}). Got: {pdb_abs}")
    return f"/runs/{rel.as_posix()}"

def auto_outdir(protein: str, prop: str, ph: str, salt: str) -> str:
    """Build a friendly per-job output folder name."""
    bits = [protein, prop]
    if prop == "electrostatics":
        if ph:
            bits.append(f"pH{ph}")
        if salt:
            bits.append(f"salt{salt}")
    return "_".join(bits)

def build_cmd(pdb_in_container: str,
              out_subdir_rel: str,
              prop: str,
              proj: str,
              cellsize: str,
              ph: str,
              salt: str):
    """
    Compose the 'docker run ... python -m surfmap.bin.surfmap ...' command.
    out_subdir_rel is a path relative to /runs (e.g. '<BATCH>/<job_outdir>').
    """
    outdir_in_container = f"/runs/{out_subdir_rel}"

    base = [
        "docker", "run", "--rm",
        "-v", f"{RUNS}:/runs",
        IMAGE,
        "python", "-m", "surfmap.bin.surfmap",
        "-pdb", pdb_in_container,
        "-tomap", prop,
        "-d", outdir_in_container,
        "-proj", proj,
        "-s", cellsize,
        "--plot3d",
        "--plot3d-point-size", PLOT3D_POINT_SIZE,
        "--csv3d",
        "--png",
        "-verbose", "2",
    ]
    if prop == "electrostatics":
        if ph:
            base += ["--ph", str(ph)]
        if salt:
            base += ["--salt", str(salt)]
    return base

def run_job(row: dict, idx: int, batch_name: str) -> int:
    protein = (row.get("protein") or "").strip()
    pdb_field = (row.get("pdb") or "").strip()
    prop = norm_prop(row.get("property"))
    ph = (row.get("ph") or "").strip()
    salt = (row.get("salt") or "").strip()
    proj = (row.get("proj") or DEFAULT_PROJ).strip().lower()
    cellsize = (row.get("cellsize") or DEFAULT_CELLSIZE).strip()
    outdir = (row.get("outdir") or "").strip().lstrip("/")  # force relative inside batch

    if not protein:
        print(f"[SKIP row {idx}] Missing 'protein'")
        return 0
    if not pdb_field:
        print(f"[SKIP row {idx}] Missing 'pdb' for protein {protein}")
        return 0
    if prop not in ("wimley_white", "electrostatics"):
        print(f"[SKIP row {idx}] Invalid property '{prop}' for protein {protein}")
        return 0
    if prop == "electrostatics" and not ph:
        print(f"[SKIP row {idx}] APBS needs 'ph' for protein {protein}")
        return 0

    pdb_abs = abs_pdb_path(pdb_field)
    if not Path(pdb_abs).is_file():
        print(f"[SKIP row {idx}] PDB not found: {pdb_abs}")
        return 0
    pdb_in_container = pdb_host_to_container(pdb_abs)

    per_job = outdir if outdir else auto_outdir(protein, prop, ph, salt)
    out_subdir_rel = f"{batch_name}/{per_job}"  # <BATCH>/<job-folder>

    print("\n=== JOB", idx, "===")
    print("Protein:", protein)
    print("PDB    :", pdb_abs)
    print("Prop   :", prop, ("(pH=" + ph + ")") if ph else "", ("(salt=" + salt + "M)") if salt else "")
    print("Outdir :", str(Path(RUNS) / out_subdir_rel))
    cmd = build_cmd(
        pdb_in_container=pdb_in_container,
        out_subdir_rel=out_subdir_rel,
        prop=prop,
        proj=proj,
        cellsize=cellsize,
        ph=ph,
        salt=salt,
    )
    print("Cmd    :", " ".join(shlex.quote(c) for c in cmd))

    try:
        res = subprocess.run(cmd, check=False)
        if res.returncode != 0:
            print(f"[JOB {idx}] FAILED with code {res.returncode}")
        else:
            print(f"[JOB {idx}] OK")
        return res.returncode
    except KeyboardInterrupt:
        print("\n[ABORTED]")
        return 130

def main():
    parser = argparse.ArgumentParser(
        description="Run multiple SURFMAP jobs described in a CSV, storing all outputs under one batch folder."
    )
    parser.add_argument("csv", help="Path to jobs CSV file")
    parser.add_argument(
        "--batch",
        help=("Name of the batch folder under RUNS (e.g. 241026_ph8_elution). "
              "If omitted, a timestamped folder like 'batch_YYYYMMDD_HHMMSS' is used.")
    )
    args = parser.parse_args()

    # Resolve CSV path
    csv_path = Path(args.csv).expanduser().resolve()
    if not csv_path.exists():
        sys.exit(f"CSV not found: {csv_path}")

    # Batch name / folder
    batch_name = args.batch.strip() if args.batch else f"batch_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    batch_dir = Path(RUNS) / batch_name
    batch_dir.mkdir(parents=True, exist_ok=True)

    print("Image :", IMAGE)
    print("RUNS  :", RUNS)
    print("CSV   :", str(csv_path))
    print("Batch :", batch_name)
    print("Batch outputs will go to:", batch_dir)

    # --- Read CSV robustly (handles BOM; auto-detect , ; or tab) ---
    with open(csv_path, "r", encoding="utf-8-sig", newline="") as f:
        # Skip comment lines starting with '#'
        lines = [ln for ln in f if not ln.lstrip().startswith("#")]

    if not lines:
        print("No jobs found in CSV.")
        sys.exit(0)

    # Try to detect delimiter from the header + first data line
    sample = "\n".join(lines[:2])
    delim = ","  # default
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=[",", ";", "\t", "|"])
        delim = dialect.delimiter
    except Exception:
        # Fallback: inspect header
        header = lines[0]
        for cand in ("\t", ";", "|", ","):
            if cand in header:
                delim = cand
                break

    rdr = csv.DictReader(lines, delimiter=delim)
    rows = list(rdr)

    print(f"Detected delimiter: {repr(delim)}")

    if not rows:
        print("No jobs found in CSV.")
        sys.exit(0)

    failures = 0
    for i, row in enumerate(rows, start=1):
        rc = run_job(row, i, batch_name=batch_name)
        if rc != 0:
            failures += 1

    print("\n--- SUMMARY ---")
    print(f"Batch   : {batch_name}")
    print(f"Outputs : {batch_dir}")
    print(f"Total rows: {len(rows)}  |  Failures: {failures}")
    sys.exit(1 if failures else 0)

if __name__ == "__main__":
    main()
