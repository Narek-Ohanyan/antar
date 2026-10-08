import gzip

import yaml

from antar.io import fitted


def test_plain_files_are_read_and_written_as_before(tmp_path):
    p = tmp_path / "a.yaml"
    fitted.write_yaml(p, {"x": [1, 2], "y": {"z": "é"}})
    assert p.read_text().startswith("x:") and fitted.read_yaml(p) == {"x": [1, 2], "y": {"z": "é"}}
    assert fitted.resolve(p) == p and fitted.write_target(p) == p


def test_a_file_named_as_compressed_is_stored_gzipped_and_found_by_its_plain_name(tmp_path, monkeypatch):
    monkeypatch.setattr(fitted, "COMPRESSED", {"big.yaml"})
    p = tmp_path / "big.yaml"
    written = fitted.write_yaml(p, {"members": {"a": list(range(1000))}})
    assert written == tmp_path / "big.yaml.gz" and not p.exists()
    with gzip.open(written, "rt") as f:
        assert yaml.safe_load(f)["members"]["a"][999] == 999
    assert fitted.resolve(p) == written and fitted.read_yaml(p)["members"]["a"][5] == 5
    assert written.stat().st_size < 0.5 * len(yaml.safe_dump({"members": {"a": list(range(1000))}}))      # the numbers compress to a fraction of the text


def test_writing_removes_a_stale_plain_copy_so_there_is_only_one(tmp_path, monkeypatch):
    monkeypatch.setattr(fitted, "COMPRESSED", {"big.yaml"})
    (tmp_path / "big.yaml").write_text("old: 1\n")
    fitted.write_yaml(tmp_path / "big.yaml", {"new": 2})
    assert not (tmp_path / "big.yaml").exists() and fitted.read_yaml(tmp_path / "big.yaml") == {"new": 2}


def test_a_missing_file_resolves_to_its_plain_name(tmp_path):
    assert fitted.resolve(tmp_path / "nothing.yaml") == tmp_path / "nothing.yaml"
