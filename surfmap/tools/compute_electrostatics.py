import argparse
import os
import sys
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


    # generate APBS input file
    outfile_in = f"{outfile_basename}.in"  # autogenerated by inputgen in same path as the pqr file
    cmd_inputgen = [sys.executable, inputgen, "--method=auto", outfile_pqr]
    logger.debug(f"Running inputgen command: {' '.join(cmd_inputgen)}")
    status = subprocess.run(args=cmd_inputgen, capture_output=True)
    if status.returncode != 0:
        logger.error(
            f"Error occured during the inputgen command, the process will stop. Status: {status.returncode}"
        )
        logger.error(f"STDERR: {status.stderr.decode(errors='ignore')}")
        exit(1)

    # edit inputgen outfile to use absolute PQR path
    logger.debug(f"Editing generated inputgen file: {outfile_in}")
    edit_inputgen(inputgen_file=outfile_in, pqrfile=outfile_pqr)

    # --- Inject ionic strength (salt) into the ELEC block ---
    if salt is not None:
        try:
            s = float(salt)
        except Exception:
            logger.error(
                "Invalid -salt value; must be a number in mol/L (e.g., 0.10)."
            )
            exit(1)

        try:
            with open(outfile_in, "r", encoding="utf-8") as f:
                lines = f.readlines()
        except FileNotFoundError:
            logger.error(f"APBS input file not found: {outfile_in}")
            exit(1)

        new_lines = []
        in_elec = False
        injected = False

        for line in lines:
            stripped = line.strip().lower()

            if stripped.startswith("elec"):
                in_elec = True
                new_lines.append(line)
                continue

            if in_elec and stripped == "end":
                # Add monovalent ions just before ELEC 'end'
                new_lines.append(f"    ion charge  1 conc {s} radius 2.0\n")
                new_lines.append(f"    ion charge -1 conc {s} radius 2.0\n")
                injected = True
                in_elec = False
                new_lines.append(line)
                continue

            new_lines.append(line)

        if not injected:
            # Fallback: append a minimal ELEC block if none was detected
            new_lines.append("\n# injected by SurfMap: ionic strength\n")
            new_lines.append("elec\n")
            new_lines.append(f"    ion charge  1 conc {s} radius 2.0\n")
            new_lines.append(f"    ion charge -1 conc {s} radius 2.0\n")
            new_lines.append("end\n")

        with open(outfile_in, "w", encoding="utf-8") as f:
            f.writelines(new_lines)
    # --- end salt injection ---

    

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
