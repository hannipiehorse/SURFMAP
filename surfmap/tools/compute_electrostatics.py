import argparse
import os
import sys
import re
import shutil
from pathlib import Path
from shutil import copyfile
import subprocess
from typing import Union, Optional

from surfmap.bin.multival_csv_to_pdb import run as run_multival_csv2pdb
from surfmap.lib.logs import get_logger

logger = get_logger(name=__name__)


def get_args():
    parser = argparse.ArgumentParser(
        description="Script to compute electrostatics. Returns PDB files with coordinates of the shell surface particles and electrostatic."
    )

    parser.add_argument(
        "-pdb",
        required=True,
        type=str,
        help="Path to the pdbfile of the protein complex",
    )

    parser.add_argument(
        "-csv",
        required=True,
        type=str,
        help="CSV file with coordinates of the shell particles.",
    )

    parser.add_argument(
        "-pqr",
        required=False,
        type=str,
        default="",
        help="PQR file used to build input file of APBS",
    )

    parser.add_argument(
        "-out",
        type=str,
        default=Path("."),
        help="Path of the output (root) directory",
    )

    parser.add_argument(
        "-subdir",
        type=str,
        default=Path("tmp-elec"),
        help="Subdirectory name inside -out for electrostatics content",
    )

    parser.add_argument(
        "-salt",
        type=float,
        default=None,
        help="Monovalent salt concentration in mol/L to pass to APBS (ionic strength). Example: -salt 0.15",
    )

    return parser.parse_args()


def edit_inputgen(inputgen_file: Union[str, Path], pqrfile: Union[str, Path]):
    with open(inputgen_file, "r", encoding="utf-8") as _file:
        data = _file.read()

    # Replace some data
    data = data.replace(Path(pqrfile).name, str(pqrfile))
    data = data.replace("pot dx pot", f"pot dx {pqrfile}")

    # Write the file out
    with open(inputgen_file, "w", encoding="utf-8") as _file:
        _file.write(data)

def normalize_apbs_in(in_path: Union[str, Path], pqr_abs: Union[str, Path]) -> None:
    """
    Ensure the APBS .in file is valid:
      - exactly ONE absolute path on the 'mol pqr ...' line
      - a 'write pot dx <absolute PQR>' line exists (and is normalized)
    """
    in_path = str(in_path)
    pqr_abs = str(pqr_abs)
    txt = Path(in_path).read_text(encoding="utf-8")

    # 1) Normalize 'mol pqr ...' to the exact absolute PQR path (once)
    txt = re.sub(
        r'^(?P<indent>\s*mol\s+pqr\s+).*$',
        lambda m: f"{m.group('indent')}{pqr_abs}",
        txt,
        flags=re.IGNORECASE | re.MULTILINE,
    )

    # 2) Ensure/normalize 'write pot dx ...' inside the first ELEC block
    if re.search(r'^\s*write\s+pot\s+dx', txt, flags=re.IGNORECASE | re.MULTILINE):
        txt = re.sub(
            r'^\s*write\s+pot\s+dx.*$',
            f"  write pot dx {pqr_abs}",
            txt,
            flags=re.IGNORECASE | re.MULTILINE,
        )
    else:
        # Inject just before the first 'end' in the first ELEC block
        def _inject_write(match):
            block = match.group(0)
            return re.sub(
                r'^\s*end\s*$',
                f"  write pot dx {pqr_abs}\nend",
                block,
                count=1,
                flags=re.IGNORECASE | re.MULTILINE,
            )
        txt = re.sub(
            r'(?is)elec\b.*?end',   # first ELEC..end
            _inject_write,
            txt,
            count=1,
        )

    Path(in_path).write_text(txt, encoding="utf-8")


def write_minimal_apbs_in(outfile_in: Union[str, Path],
                          pqrfile: Union[str, Path],
                          salt_m: Optional[float] = None,
                          temperature_k: float = 298.15) -> None:
    """
    Minimal APBS input as a last-resort fallback.
    Note: some APBS builds still want explicit grid params; whenever possible
    we prefer using `inputgen.py` so it computes them.
    """
    pqrfile = str(pqrfile)
    base = str(Path(pqrfile))
    lines = []
    lines.append("read")
    lines.append(f"  mol pqr {pqrfile}")
    lines.append("end")
    lines.append("")
    lines.append("elec")
    lines.append("  mg-auto")
    lines.append("  lpbe")
    lines.append("  mol 1")
    lines.append("  bcfl sdh")
    lines.append("  pdie 2.0")
    lines.append("  sdie 78.54")
    lines.append("  srfm smol")
    lines.append("  chgm spl4")
    lines.append("  srad 1.4")
    lines.append("  swin 0.3")
    lines.append(f"  temp {temperature_k}")
    if salt_m is not None and float(salt_m) > 0.0:
        s = float(salt_m)
        lines.append(f"  ion charge  1 conc {s} radius 2.0")
        lines.append(f"  ion charge -1 conc {s} radius 2.0")
    lines.append("  calcenergy no")
    lines.append("  calcforce  no")
    lines.append(f"  write pot dx {base}")
    lines.append("end")
    lines.append("")
    lines.append("quit")
    Path(outfile_in).write_text("\n".join(lines) + "\n", encoding="utf-8")

def run(
    pdb_filename: Union[str, Path],
    csv_filename: Union[str, Path],
    force_field: str = "CHARMM",
    pqr_filename: Union[str, Path] = "",
    out_dir: Union[str, Path] = ".",
    out_subdir: Union[str, Path] = "tmp-elec",
    ph: Optional[float] = None,
    salt: Optional[float] = None,
) -> str:
    """Compute electrostatics potential and return the shell PDB with values in B-factor."""
    PDB2PQR_FORCE_FIELDS = ["AMBER", "CHARMM", "PARSE", "TYL06", "PEOEPB", "SWANSON"]

    # check if force field is allowed
    if force_field.upper() not in PDB2PQR_FORCE_FIELDS:
        logger.error(f"Error, the force field {force_field} is not accepted.")
        logger.error(f"Accepted force fields are: {' '.join(PDB2PQR_FORCE_FIELDS)}.")
        logger.error("Exiting now.\n")
        exit(1)

    # define access to useful APBS executables
    apbs_bin = shutil.which("apbs")
    env_apbs = os.getenv("APBS")

    if apbs_bin:
        PATH_APBS = Path(apbs_bin).resolve().parent.parent  # .../bin/apbs → prefix
    elif env_apbs:
        PATH_APBS = Path(env_apbs)
    else:
        raise RuntimeError(
            "APBS not found. Install it so 'apbs' is on PATH, or set APBS=/path/to/prefix"
        )
    apbs = str(PATH_APBS / "bin" / "apbs")
    inputgen = f"{str(PATH_APBS)}/share/apbs/tools/manip/inputgen.py"
    multivalue = f"{str(PATH_APBS)}/share/apbs/tools/bin/multivalue"

    # Create directory containing electrostatics calculation related files
    outdir_root = Path(out_dir)
    outdir = (
        outdir_root
        if (out_subdir and str(outdir_root).endswith(str(out_subdir)))
        else (outdir_root / out_subdir)
    )
    outdir.mkdir(exist_ok=True, parents=True)
    outfile_basename = str(outdir / Path(pdb_filename).stem)

    # compute pqr file if not given as input
    outfile_pqr = f"{outfile_basename}.pqr"
    if not pqr_filename:
        cmd_pdb2pqr = [
            "pdb2pqr30",
            "--ff", force_field,
            "--whitespace",
        ]
        # optional pH handling
        if ph is not None:
            cmd_pdb2pqr += ["--with-ph", str(ph)]

        # (sometimes helps on tricky PDBs; safe to leave on)
        cmd_pdb2pqr += ["--keep-chain", "--drop-water"]

        # positional args (add ONCE): input pdb, then output pqr
        cmd_pdb2pqr += [str(pdb_filename), outfile_pqr]

        logger.debug("Convert PDB to PQR format: " + " ".join(cmd_pdb2pqr))
        status = subprocess.run(cmd_pdb2pqr, capture_output=True)
        if status.returncode != 0:
            logger.error(
                f"PDB2PQR failed (status {status.returncode}). "
                f"STDERR:\n{status.stderr.decode(errors='ignore')}"
            )
            exit(1)

        # NEW: ensure the file is actually there and non-empty
        if not Path(outfile_pqr).is_file() or Path(outfile_pqr).stat().st_size == 0:
            logger.error(f"PDB2PQR reported success but '{outfile_pqr}' was not created.")
            exit(1)
    else:
        logger.debug(f"Copying user-given PQR file as {outfile_pqr}")
        copyfile(src=pqr_filename, dst=outfile_pqr)


    # --- Generate APBS input file robustly ---
    pqr_path = Path(outfile_pqr)
    pqr_dir = pqr_path.parent          # run inputgen IN the PQR directory
    pqr_base = pqr_path.name           # pass only the basename to inputgen
    outfile_in = str(pqr_path.with_suffix(".in"))

    cmd_inputgen = [sys.executable, inputgen, "--method=auto", "--potdx"]
    if salt is not None:
        try:
            s = float(salt)
            cmd_inputgen += ["--istrng", str(s)]
        except Exception:
            logger.error("Invalid -salt value; must be a number in mol/L (e.g., 0.15).")
            exit(1)
    cmd_inputgen += [pqr_base]

    logger.debug(f"Running inputgen in cwd={pqr_dir}: {' '.join(cmd_inputgen)}")
    status = subprocess.run(cmd_inputgen, cwd=str(pqr_dir), capture_output=True)

    # Track whether WE wrote the file ourselves (so we can skip edit_inputgen)
    wrote_minimal = False
    wrote_from_stdout = False

    # Case A: inputgen wrote an .in file next to the PQR
    if Path(outfile_in).is_file():
        pass  # good

    else:
        # Case B: try to capture stdout
        stdout_text = status.stdout.decode(errors="ignore")
        looks_like_apbs_in = ("read" in stdout_text.lower()) or ("elec" in stdout_text.lower())
        if looks_like_apbs_in and stdout_text.strip():
            # Sanitize leading banners (drop until first 'read' or 'elec' at line start)
            lines = stdout_text.splitlines()
            def _is_start(line: str) -> bool:
                ls = line.strip().lower()
                return ls.startswith("read") or ls.startswith("elec")
            while lines and not _is_start(lines[0]):
                lines.pop(0)
            if lines:
                Path(outfile_in).write_text("\n".join(lines) + "\n", encoding="utf-8")
                wrote_from_stdout = True
            else:
                # fallback below
                pass

        # Case C: nothing usable from stdout → last resort minimal file
        if not Path(outfile_in).is_file():
            logger.warning(f"inputgen produced no usable file or content; writing minimal APBS .in: {outfile_in}")
            write_minimal_apbs_in(outfile_in=outfile_in, pqrfile=outfile_pqr, salt_m=salt)
            wrote_minimal = True

    # If we did NOT write the file ourselves, normalize paths via edit_inputgen
    # (Only needed when inputgen used basenames in the file.)
    if not wrote_minimal and not wrote_from_stdout:
        # Only rewrite if the file still contains just the basename; if it's already absolute, skip.
        txt = Path(outfile_in).read_text(encoding="utf-8")
        if pqr_base in txt and str(pqr_path) not in txt:
            logger.debug(f"Editing generated inputgen file: {outfile_in} (making PQR path absolute)")
            edit_inputgen(inputgen_file=outfile_in, pqrfile=outfile_pqr)
        else:
            logger.debug(f"Skipping edit_inputgen: PQR path already absolute in {outfile_in}")
    # --- end robust inputgen ---



    # final normalize (fixes 'mol pqr' absolute path and guarantees 'write pot dx')
    logger.debug(f"Normalizing APBS .in file: {outfile_in}")
    normalize_apbs_in(outfile_in, outfile_pqr)


    # --- Inject ionic strength (salt) into the ELEC block (only if missing) ---
    if salt is not None:
        # Only inject if inputgen did NOT already add ions
        txt = Path(outfile_in).read_text(encoding="utf-8")
        if re.search(r'^\s*ion\s+charge', txt, flags=re.IGNORECASE | re.MULTILINE):
            logger.debug("Skipping ion injection: ions already present in ELEC block.")
        else:
            try:
                s = float(salt)
            except Exception:
                logger.error("Invalid -salt value; must be a number in mol/L (e.g., 0.10).")
                exit(1)

            lines = txt.splitlines(keepends=False)
            new_lines = []
            in_elec = False
            injected = False
            for line in lines:
                stripped = line.strip().lower()
                if stripped.startswith("elec"):
                    in_elec = True
                    new_lines.append(line); continue
                if in_elec and stripped == "end":
                    new_lines.append(f"    ion charge  1 conc {s} radius 2.0")
                    new_lines.append(f"    ion charge -1 conc {s} radius 2.0")
                    injected = True
                    in_elec = False
                    new_lines.append(line); continue
                new_lines.append(line)
            if not injected:
                new_lines.append("\n# injected by SurfMap: ionic strength")
                new_lines.append("elec")
                new_lines.append(f"    ion charge  1 conc {s} radius 2.0")
                new_lines.append(f"    ion charge -1 conc {s} radius 2.0")
                new_lines.append("end")
            Path(outfile_in).write_text("\n".join(new_lines) + "\n", encoding="utf-8")


    

    # compute APBS electrostatics calculation
    outfile_pot = f"{outfile_pqr}.dx"  # default output name of APBS
    cmd_apbs = [apbs, outfile_in]
    logger.debug(f"Running APBS command: {cmd_apbs}")
    status = subprocess.run(args=cmd_apbs, capture_output=True)
    if status.returncode != 0:
        logger.error(
            f"Error occured during the APBS command, the process will stop. Status: {status.returncode}"
        )
        logger.error(f"STDERR: {status.stderr.decode(errors='ignore')}")
        exit(1)

    # cross the MSMS generated CSV file with the potential gridfile created by APBS
    outfile_multivalue = f"{outfile_basename}.mult"
    cmd_multivalue = [multivalue, str(csv_filename), outfile_pot, outfile_multivalue]
    logger.debug(f"Running multivalue command: {cmd_multivalue}")
    status = subprocess.run(args=cmd_multivalue, capture_output=True)
    if status.returncode != 0:
        logger.error(
            f"Error occured during the multivalue command, the process will stop. Status: {status.returncode}"
        )
        logger.error(f"STDERR: {status.stderr.decode(errors='ignore')}")
        exit(1)

    # convert .csv file into a PDB file
    outfile_pdb = f"{outfile_basename}_shell.pdb"
    logger.debug("Converting CSV to PDB file")
    run_multival_csv2pdb(ifile=outfile_multivalue, ofile=outfile_pdb)

    # remove io.mc generated by APBS
    io_apbs = Path.cwd() / "io.mc"
    logger.debug(f"Removing APBS derived log file: {io_apbs}")
    io_apbs.unlink(missing_ok=True)

    return outfile_pdb


def main():
    args = get_args()
    return run(
        pdb_filename=args.pdb,
        csv_filename=args.csv,
        pqr_filename=args.pqr,
        out_dir=args.out,
        out_subdir=args.subdir,
        salt=args.salt,
    )
