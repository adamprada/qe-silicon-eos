"""Download the SSSP efficiency PBE Si pseudopotential and verify pw.x is available."""

import hashlib
import os
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

PSEUDO_DIR = Path("pseudo")
PSEUDO_FILENAME = "Si.pbe-n-rrkjus_psl.1.0.0.UPF"
PSEUDO_URL = (
    "https://pseudopotentials.quantum-espresso.org/upf_files/"
    "Si.pbe-n-rrkjus_psl.1.0.0.UPF"
)
PSEUDO_MD5 = "fa25574f73a70a4139f2adfbefec430c"


def verify_md5(path: Path, expected: str) -> bool:
    h = hashlib.md5()
    h.update(path.read_bytes())
    return h.hexdigest() == expected


def download_pseudopotential() -> None:
    PSEUDO_DIR.mkdir(exist_ok=True)
    dest = PSEUDO_DIR / PSEUDO_FILENAME

    if dest.exists() and dest.stat().st_size > 0:
        print(f"Pseudopotential already present: {dest}")
    else:
        print(f"Downloading {PSEUDO_FILENAME} ...")
        urllib.request.urlretrieve(PSEUDO_URL, dest)
        if dest.stat().st_size == 0:
            print(f"ERROR: downloaded file is empty: {dest}")
            sys.exit(1)
        print(f"Downloaded: {dest} ({dest.stat().st_size} bytes)")

    if not verify_md5(dest, PSEUDO_MD5):
        print(f"ERROR: MD5 mismatch for {dest}. File may be corrupt.")
        sys.exit(1)

    print(f"Filename : {PSEUDO_FILENAME}")
    print(f"Path     : {dest.resolve()}")


def check_pw_x() -> None:
    pw_path = shutil.which("pw.x")
    if pw_path is None:
        print(
            "ERROR: pw.x not found on PATH.\n"
            "Activate the conda environment first:\n"
            "  conda activate qe-silicon-eos"
        )
        sys.exit(1)

    result = subprocess.run(
        ["pw.x", "--version"],
        capture_output=True,
        text=True,
    )
    version_output = (result.stdout + result.stderr).strip()
    first_line = version_output.splitlines()[0] if version_output else "(unknown)"
    print(f"pw.x     : {pw_path}")
    print(f"Version  : {first_line}")


if __name__ == "__main__":
    download_pseudopotential()
    check_pw_x()
