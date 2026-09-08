"""Run our witness through the OTHER author's verifier, unmodified.

    python crosscheck.py
    python crosscheck.py certificates/ours_m128_n384_W3.json

A verifier written by the same people who wrote the witness is worth something,
but not much: if both share a mistaken idea of what the constant is, both agree
and both are wrong. So the strongest check available here is somebody else's
code, which was written before ours existed and knows nothing about it.

`third_party/verify_bound.py` is Russell's, kept byte for byte as published
(`techno-optimist/minimum-autocorrelation-bound`, LICENSE beside it). It reads a
fixed path -- `<repo>/certs/certificate_n480.json` -- so this script builds that
layout in a temporary directory, drops our certificate in under the name it
expects, and runs the file untouched. Nothing is patched, monkey-patched, or
passed in: if his script rejects our witness, it says so.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
THEIRS = HERE / "third_party" / "verify_bound.py"
DEFAULT = HERE / "certificates" / "ours_evolved_m370_n481_W1.3.json"


def main() -> int:
    cert = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT
    if not cert.is_file():
        print(f"no such certificate: {cert}")
        return 1
    if not THEIRS.is_file():
        print(f"the third-party verifier is missing: {THEIRS}")
        return 1

    print(f"witness   {cert.name}")
    print(f"verifier  {THEIRS.name}, unmodified, from techno-optimist/minimum-autocorrelation-bound")
    print()

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / "certs").mkdir()
        (root / "tools").mkdir()
        # The name his script looks for. The file is ours; only the filename is
        # borrowed, because his path is hardcoded.
        shutil.copy(cert, root / "certs" / "certificate_n480.json")
        shutil.copy(THEIRS, root / "tools" / "verify_bound.py")

        finished = subprocess.run(
            [sys.executable, str(root / "tools" / "verify_bound.py")],
            capture_output=True, text=True,
        )

    print(finished.stdout.rstrip())
    if finished.stderr.strip():
        print(finished.stderr.rstrip())
    print()
    if finished.returncode == 0:
        print("his verifier accepts our witness (exit 0)")
    else:
        print(f"his verifier REFUSED our witness (exit {finished.returncode})")
    return finished.returncode


if __name__ == "__main__":
    raise SystemExit(main())
