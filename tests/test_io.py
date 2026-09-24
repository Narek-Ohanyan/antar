from antar.io import manifest


def _entry(local_path=None, checksum=None):
    return manifest.ManifestEntry(
        variable="tas",
        source="CHELSA V2.1",
        version="2.1",
        url="https://chelsa-climate.org/downloads/",
        citation="Karger et al. 2017, Scientific Data",
        license="CC BY 4.0",
        access_date="2026-09-24",
        local_path=local_path,
        checksum_sha256=checksum,
    )


def test_write_then_load_manifest_roundtrips(tmp_path):
    path = tmp_path / "manifest.yaml"
    entries = [_entry(), _entry(local_path="tas/armenia.tif")]
    manifest.write_manifest(entries, path)
    loaded = manifest.load_manifest(path)
    assert loaded == entries


def test_load_manifest_missing_file_is_empty_list(tmp_path):
    path = tmp_path / "does_not_exist.yaml"
    (tmp_path / "placeholder").write_text("")  # keep tmp_path non-empty; irrelevant to the assertion
    manifest.write_manifest([], path)
    assert manifest.load_manifest(path) == []


def test_sha256_checksum_matches_known_value(tmp_path):
    f = tmp_path / "data.txt"
    f.write_bytes(b"antar")
    # sha256("antar") computed independently
    import hashlib

    assert manifest.sha256_checksum(f) == hashlib.sha256(b"antar").hexdigest()


def test_verify_entry_true_only_when_file_present_and_checksum_matches(tmp_path):
    (tmp_path / "data").mkdir()
    f = tmp_path / "data" / "tas.tif"
    f.write_bytes(b"raster-bytes")
    good = _entry(local_path="data/tas.tif", checksum=manifest.sha256_checksum(f))
    bad_checksum = _entry(local_path="data/tas.tif", checksum="0" * 64)
    missing_file = _entry(local_path="data/missing.tif", checksum="0" * 64)
    no_local_path = _entry()

    assert manifest.verify_entry(good, base_dir=tmp_path) is True
    assert manifest.verify_entry(bad_checksum, base_dir=tmp_path) is False
    assert manifest.verify_entry(missing_file, base_dir=tmp_path) is False
    assert manifest.verify_entry(no_local_path, base_dir=tmp_path) is False
