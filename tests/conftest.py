"""Make the package and the shared case table importable from a bare checkout."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

for location in (ROOT / "src", Path(__file__).resolve().parent):
    if str(location) not in sys.path:
        sys.path.insert(0, str(location))
