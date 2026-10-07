"""Puts the built web interface into the website repository's antar/ folder, which Hostinger deploys to antar.narekohanyan.com.

    python3 scripts/publish_site.py --site-repo <clone of Narek-Ohanyan/narekohanyan.com>            # show what would change
    python3 scripts/publish_site.py --site-repo <clone> --push                                       # commit and push (this goes live)

Guards, each of them tested (tests/test_publish_site.py): the clone's origin must be the website repository; the working tree must be clean and on main;
the only paths that may change are under antar/; nothing is pushed without --push; when nothing differs nothing is committed. The commit message is plain
and names the build stamp and the ANTAR source commit.
"""
import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_site  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
WEBSITE = "Narek-Ohanyan/narekohanyan.com"
FOLDER = "antar"


def git(repo, *args, check=True):
    r = subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True)
    if check and r.returncode:
        raise SystemExit(f"git {' '.join(args)} failed in {repo}: {r.stderr.strip()}")
    return r.stdout.strip()


def publish(repo, push=False, expect=WEBSITE, skip_data=False, message=None):
    repo = Path(repo).resolve()
    origin = git(repo, "remote", "get-url", "origin")
    if expect not in origin:
        raise SystemExit(f"refusing: origin of {repo} is {origin}, not {expect}")
    if git(repo, "rev-parse", "--abbrev-ref", "HEAD") != "main":
        raise SystemExit("refusing: the website clone is not on main")
    if git(repo, "status", "--porcelain"):
        raise SystemExit("refusing: the website clone has uncommitted changes")
    git(repo, "fetch", "origin", "main")
    if git(repo, "rev-parse", "HEAD") != git(repo, "rev-parse", "origin/main"):
        raise SystemExit("refusing: the website clone is not up to date with origin/main (pull first)")

    with tempfile.TemporaryDirectory() as tmp:
        if not skip_data:
            build_site.build_data()
        bundle = build_site.assemble(Path(tmp) / "site")
        problems = build_site.check_bundle(bundle)
        if problems:
            raise SystemExit("refusing: the bundle has problems:\n  " + "\n  ".join(problems))
        target = repo / FOLDER
        if target.exists():
            shutil.rmtree(target)
        shutil.copytree(bundle, target)
        stamp = [l for l in (bundle / "index.html").read_text().splitlines() if "ANTAR_BUILD" in l][0].split('"')[1]

    git(repo, "add", "-A", FOLDER)
    changed = [l[3:] for l in git(repo, "status", "--porcelain").splitlines()]
    outside = [p for p in changed if not p.startswith(FOLDER + "/")]
    if outside:
        raise SystemExit("refusing: paths outside antar/ would change: " + ", ".join(outside[:5]))
    if not changed:
        print("nothing to publish: antar/ already matches this build")
        return None
    print(f"{len(changed)} files would change in {FOLDER}/ (build {stamp}):")
    for line in git(repo, "diff", "--cached", "--stat", "--", FOLDER).splitlines()[-6:]:
        print("  " + line)
    if not push:
        git(repo, "reset", "-q")
        shutil.rmtree(repo / FOLDER, ignore_errors=True)
        git(repo, "checkout", "--", FOLDER, check=False)
        print("dry run: nothing committed (add --push to publish)")
        return None
    src = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "--", "ui", "src", "configs", "deploy"], capture_output=True, text=True).stdout.strip())
    if dirty:
        print("WARNING: the ANTAR working tree has uncommitted changes under ui/, src/, configs/ or deploy/; the published build is not exactly commit " + src)
    git(repo, "commit", "-q", "-m", message or f"Update the ANTAR web interface (build {stamp}, source {src})")
    git(repo, "push", "origin", "main")
    sha = git(repo, "rev-parse", "--short", "HEAD")
    print(f"pushed {sha}; Hostinger deploys it to antar.narekohanyan.com")
    return sha


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--site-repo", required=True)
    ap.add_argument("--push", action="store_true")
    ap.add_argument("--skip-data", action="store_true")
    a = ap.parse_args()
    publish(a.site_repo, push=a.push, skip_data=a.skip_data)


if __name__ == "__main__":
    main()
