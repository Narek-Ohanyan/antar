"""scripts/publish_site.py only ever changes antar/ in the website repository, only pushes when told to, and refuses a wrong or dirty clone."""
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import publish_site  # noqa: E402


def run(*args, cwd):
    return subprocess.run(args, cwd=cwd, check=True, capture_output=True, text=True).stdout.strip()


@pytest.fixture()
def site(tmp_path):
    origin = tmp_path / "website-origin.git"
    run("git", "init", "--bare", "-b", "main", str(origin), cwd=tmp_path)
    clone = tmp_path / "website"
    run("git", "clone", str(origin), str(clone), cwd=tmp_path)
    run("git", "config", "user.name", "Test Author", cwd=clone)
    run("git", "config", "user.email", "test@example.org", cwd=clone)
    (clone / "index.html").write_text("<h1>personal site</h1>")
    (clone / "about.html").write_text("about")
    run("git", "add", "-A", cwd=clone)
    run("git", "commit", "-m", "Start the website", cwd=clone)
    run("git", "push", "origin", "main", cwd=clone)
    return clone, origin


def publish(clone, **kw):
    return publish_site.publish(clone, expect="website-origin.git", skip_data=True, **kw)


def test_a_dry_run_changes_nothing(site):
    clone, origin = site
    head = run("git", "rev-parse", "HEAD", cwd=clone)
    assert publish(clone) is None
    assert run("git", "rev-parse", "HEAD", cwd=clone) == head and run("git", "status", "--porcelain", cwd=clone) == ""
    assert not (clone / "antar").exists() and run("git", "rev-parse", "main", cwd=origin) == head


def test_a_push_adds_only_antar_with_a_plain_message_and_leaves_the_rest_alone(site):
    clone, origin = site
    sha = publish(clone, push=True)
    assert sha and run("git", "rev-parse", "--short", "main", cwd=origin) == sha
    changed = run("git", "show", "--name-only", "--format=", "HEAD", cwd=clone).splitlines()
    assert changed and all(p.startswith("antar/") for p in changed)
    assert "antar/index.html" in changed and "antar/.htaccess" in changed and "antar/data/manifest.json" in changed
    msg = run("git", "log", "-1", "--format=%B", cwd=clone)
    assert msg.startswith("Update the ANTAR web interface (build ") and "Co-Authored" not in msg and "Generated" not in msg
    assert (clone / "index.html").read_text() == "<h1>personal site</h1>" and (clone / "about.html").read_text() == "about"


def test_publishing_twice_is_a_no_op_and_stale_files_are_removed(site):
    clone, _ = site
    publish(clone, push=True)
    head = run("git", "rev-parse", "HEAD", cwd=clone)
    assert publish(clone, push=True) is None and run("git", "rev-parse", "HEAD", cwd=clone) == head
    (clone / "antar" / "old-file.txt").write_text("left over")
    run("git", "add", "-A", cwd=clone)
    run("git", "commit", "-m", "Add a stale file", cwd=clone)
    run("git", "push", "origin", "main", cwd=clone)
    publish(clone, push=True)
    assert not (clone / "antar" / "old-file.txt").exists()
    assert (clone / "about.html").exists()


def test_it_refuses_a_dirty_clone_a_wrong_remote_and_a_stale_clone(site, tmp_path):
    clone, origin = site
    (clone / "about.html").write_text("edited")
    with pytest.raises(SystemExit, match="uncommitted"):
        publish(clone)
    run("git", "checkout", "--", "about.html", cwd=clone)
    with pytest.raises(SystemExit, match="origin"):
        publish_site.publish(clone, expect="Narek-Ohanyan/narekohanyan.com", skip_data=True)
    other = tmp_path / "other"
    run("git", "clone", str(origin), str(other), cwd=tmp_path)
    run("git", "config", "user.name", "T", cwd=other)
    run("git", "config", "user.email", "t@example.org", cwd=other)
    (other / "new.html").write_text("n")
    run("git", "add", "-A", cwd=other)
    run("git", "commit", "-m", "Elsewhere", cwd=other)
    run("git", "push", "origin", "main", cwd=other)
    with pytest.raises(SystemExit, match="up to date"):
        publish(clone)
