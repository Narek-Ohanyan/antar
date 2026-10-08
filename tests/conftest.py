"""Test set-up: the web data (ui/data/manifest.json and the grid files) is generated from the committed results by ``ui/build_data.py`` and is not itself committed, so a fresh checkout has none.
Several tests (the site bundle, the publisher, the photograph credits, the decision-page consistency checks) read it, and some decide at import time whether to skip, so it is built here, before the tests are
collected, when it is missing. The build needs only numpy and pyyaml and takes seconds."""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def pytest_configure(config):
    if not (ROOT / "ui" / "data" / "manifest.json").exists():
        subprocess.run([sys.executable, str(ROOT / "ui" / "build_data.py")], check=True, cwd=ROOT, stdout=subprocess.DEVNULL)
