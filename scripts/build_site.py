"""Builds the deployable copy of the web interface and checks it.

    python3 scripts/build_site.py [--out _site] [--skip-data]

1. runs ui/build_data.py (unless --skip-data), so the data files and the asset stamps are current;
2. copies ui/ without Python files and caches;
3. adds deploy/.htaccess (the Apache/LiteSpeed rules for antar.narekohanyan.com) and robots.txt;
4. checks the result: every local file index.html loads exists, every file the scripts fetch exists, every ?v= stamp is the build stamp, the
   photographs in the manifest are present, and nothing that should stay behind (Python, caches) was copied.
Exit status 1 if anything is wrong. Nothing is uploaded here; scripts/publish_site.py puts the result into the website repository.
"""
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
UI = ROOT / "ui"
HTACCESS = ROOT / "deploy" / ".htaccess"
ROBOTS = "User-agent: *\nAllow: /\n"
MAX_BYTES = 16_000_000

# files the scripts fetch at run time (they are not named in index.html)
FETCHED = ["data/manifest.json", "data/methodology.json", "assets/relief.svg", "assets/architecture.svg", "assets/katex/katex.min.js", "assets/katex/katex.min.css",
           "assets/map/grid.json", "assets/map/checksums.json", "assets/map/borders.geojson", "assets/map/elevation.bin", "assets/map/region.bin", "assets/map/forest.bin", "assets/map/woodland.bin",
           "assets/map/water.bin", "assets/map/forest_broadleaf.bin", "assets/map/forest_oak.bin", "assets/map/forest_pine.bin", "assets/map/forest_juniper.bin",
           "assets/hero.mp4", "assets/hero-poster.jpg", "assets/author.jpg", "assets/antar_logo.jpeg", "assets/photos/credits.json", ".htaccess", "robots.txt"]


def build_data():
    subprocess.run([sys.executable, str(UI / "build_data.py")], check=True, cwd=ROOT, stdout=subprocess.DEVNULL)


def assemble(out):
    out = Path(out)
    if out.exists():
        shutil.rmtree(out)
    shutil.copytree(UI, out, ignore=shutil.ignore_patterns("__pycache__", "*.py", "*.pyc", ".DS_Store"))
    shutil.copy2(HTACCESS, out / ".htaccess")
    (out / "robots.txt").write_text(ROBOTS)
    return out


def check_bundle(out):
    """Returns a list of problems; empty means the bundle is complete and consistent."""
    out = Path(out)
    bad = []
    index = (out / "index.html").read_text() if (out / "index.html").exists() else ""
    if not index:
        return ["index.html is missing"]
    build = re.search(r'window\.ANTAR_BUILD = "([A-Za-z0-9]+)"', index)
    if not build:
        bad.append("index.html has no ANTAR_BUILD stamp")
    for ref in re.findall(r'(?:src|href)="([^"]+)"', index):
        if re.match(r"(https?:|data:|#|mailto:)", ref):
            continue
        path, _, query = ref.partition("?")
        if not (out / path).is_file():
            bad.append(f"index.html loads {path}, which is not in the bundle")
        if query.startswith("v=") and build and query[2:] != build.group(1):
            bad.append(f"{path} carries stamp {query[2:]}, not the build stamp {build.group(1)}")
    for f in FETCHED:
        if not (out / f).is_file():
            bad.append(f"{f} is missing")
    manifest = json.loads((out / "data" / "manifest.json").read_text()) if (out / "data" / "manifest.json").exists() else {}
    if (manifest.get("aegis") or {}).get("units_file") and not (out / "data" / manifest["aegis"]["units_file"]).is_file():
        bad.append(f"data/{manifest['aegis']['units_file']} is missing (the Decision page map needs it)")
    for gid, info in (manifest.get("grids") or {}).items():
        if not (out / "data" / info["file"]).is_file():
            bad.append(f"grid {gid}: data/{info['file']} is missing")
    for ph in manifest.get("photos") or []:
        for w in (800, 1600):
            if not (out / "assets" / "photos" / f"{ph['file']}-{w}.jpg").is_file():
                bad.append(f"photograph {ph['file']}-{w}.jpg is missing")
    stray = [str(p.relative_to(out)) for p in out.rglob("*") if p.suffix in {".py", ".pyc"} or p.name == "__pycache__"]
    if stray:
        bad.append("copied by mistake: " + ", ".join(stray[:5]))
    total = sum(p.stat().st_size for p in out.rglob("*") if p.is_file())
    if total > MAX_BYTES:
        bad.append(f"the bundle is {total / 1e6:.1f} MB; the limit is {MAX_BYTES / 1e6:.0f} MB")
    return bad


def main():
    out = Path(sys.argv[sys.argv.index("--out") + 1]) if "--out" in sys.argv else ROOT / "_site"
    if "--skip-data" not in sys.argv:
        build_data()
    assemble(out)
    problems = check_bundle(out)
    files = [p for p in out.rglob("*") if p.is_file()]
    print(f"{out}: {len(files)} files, {sum(p.stat().st_size for p in files) / 1e6:.1f} MB")
    for p in problems:
        print("PROBLEM:", p)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
