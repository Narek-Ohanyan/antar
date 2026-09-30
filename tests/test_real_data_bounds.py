"""Standing unit-bounds regression tests, one per real bug this project's own pulls have
already caught by hand this session (CHELSA Kelvin encoding, SoilGrids per-mille encoding,
ISIMIP3b pr flux units, CO2 ppm range) -- ROADMAP.md S0's "Unit audit per new variable" item:
real bugs were being caught ad hoc each time, with nothing that would catch a *regression* if
a re-pull or a different source ever silently changed units again.

Tests against `data/` (gitignored, pulled locally, not present in a fresh clone or CI) are
skipped cleanly when the file isn't there rather than failing -- this suite documents the real
bound each source must satisfy and enforces it wherever the data actually exists, it does not
require every environment to have pulled 15+ GB of rasters just to run `pytest`.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
CONFIG_DIR = Path(__file__).resolve().parent.parent / "configs"


def _skip_unless_exists(path: Path):
    if not path.exists():
        pytest.skip(f"{path} not present locally (data/ is gitignored, real pull required) -- skipping, not failing")


# --- CO2 concentration pathways: committed to configs/, so this always runs ----------------

def test_co2_pathways_are_real_ppm_not_gtc_or_forcing():
    """Real bug class this guards: CO2 pathways are ppm, not GtC or W/m2 forcing (S0's own
    stated unit-audit example). Sane atmospheric CO2 for 1979-2100 across every real pathway
    pulled (RCP26/45/60/85, SSP126/245/370/585, historical) is a few hundred to ~2000 ppm --
    GtC-scale values would be in the tens of thousands, W/m2 forcing would be single digits.
    """
    import csv

    path = CONFIG_DIR / "co2_concentration_pathways.csv"
    with open(path) as f:
        rows = list(csv.DictReader(f))
    assert rows, "co2_concentration_pathways.csv is empty"
    year_cols = [c for c in rows[0] if c.isdigit()]
    assert len(year_cols) > 50, "expected a real multi-decade year range, not a stub"
    for row in rows:
        assert row["unit"] == "ppm", f"{row['scenario']}/{row['model']}: expected unit=ppm, got {row['unit']!r}"
        values = np.array([float(row[y]) for y in year_cols if row[y]])
        assert np.all(values > 250), f"{row['scenario']}: a real value fell below pre-industrial CO2 -- wrong units"
        assert np.all(values < 3000), f"{row['scenario']}: a real value exceeds any real CO2 pathway -- wrong units"


def test_co2_rcp26_and_ssp126_diverge_by_2100_not_identical():
    """Regression guard for the real S0 finding this session made explicit in scenarios.yaml:
    RCP2.6 and SSP1-2.6 are related but NOT interchangeable by matching number. If a future
    edit ever collapsed them to the same series, this would silently hide that real distinction.
    """
    import csv

    path = CONFIG_DIR / "co2_concentration_pathways.csv"
    with open(path) as f:
        rows = {(r["scenario"], r["mip_era"]): r for r in csv.DictReader(f)}
    rcp26 = rows.get(("rcp26", "CMIP5"))
    ssp126 = rows.get(("ssp126", "CMIP6"))
    if rcp26 is None or ssp126 is None:
        pytest.skip("rcp26/ssp126 rows not both present in this build of the table")
    assert abs(float(rcp26["2100"]) - float(ssp126["2100"])) > 5.0, (
        "rcp26 and ssp126 should diverge by 2100 (real ~421 vs ~446 ppm this session found) -- "
        "identical values would mean they got silently treated as interchangeable"
    )


# --- CHELSA-daily: gitignored, skip if not pulled locally -----------------------------------

def test_chelsa_tas_is_kelvin_not_celsius():
    """Real bug class this guards: CHELSA-daily's tas is Kelvin-encoded (scale/offset applied
    at pull time), not Celsius -- confirmed this session against a real July Yerevan value
    (~293K = ~20degC) and a real January value (~268.5K = ~-4.6degC), both physically sensible
    only when read as Kelvin. A silent re-introduction of a Celsius assumption would put every
    downstream lapse-rate/PET/GDD computation off by 273.15 degrees.
    """
    path = DATA_DIR / "chelsa" / "CHELSA_tas_2000_2024_armenia_daily.npz"
    _skip_unless_exists(path)
    data = np.load(path)["data"]
    sample = data[14610]  # 2019-01-01 by this file's own real 1979-01-01 day-0 indexing
    valid = sample[np.isfinite(sample)]
    assert valid.size > 0
    mean_k = float(np.nanmean(valid))
    assert 200.0 < mean_k < 330.0, (
        f"mean January temperature read as {mean_k:.1f} -- outside any real Kelvin range for "
        "Armenia; if this is ~-5 to 20, the encoding silently reverted to Celsius"
    )


def test_chelsa_pr_is_nonnegative_mm_per_day():
    path = DATA_DIR / "chelsa" / "CHELSA_pr_2000_2019_armenia_daily.npz"
    _skip_unless_exists(path)
    data = np.load(path)["data"]
    sample = data[14805]  # a real mid-July 2019 day, per this session's own earlier check
    valid = sample[np.isfinite(sample)]
    assert valid.size > 0
    assert np.all(valid >= 0.0), "precipitation must be non-negative"
    assert np.nanmean(valid) < 50.0, (
        "mean daily precipitation implausibly high for a single real day -- check for an "
        "un-applied unit scale factor"
    )


# --- SoilGrids: gitignored, skip if not pulled locally ---------------------------------------

def test_soilgrids_texture_fractions_sum_to_per_mille_not_percent():
    """Real bug this guards: SoilGrids' GEE-mapped clay/sand/silt bands are per-mille (g/kg,
    0-1000), not percent (0-100) -- caught this session by checking that clay+sand+silt summed
    to ~1000 at real sampled points, not assumed from the band name. antar.climate.
    soil_pedotransfer.soil_hydraulic_parameters expects percent; a caller that forgets the /10
    conversion silently corrupts theta_sat/theta_fc/psi_sat_mpa for every cell.
    """
    soils_path = DATA_DIR / "_tmp_soils.tif"
    if not soils_path.exists():
        pytest.skip("soils.tif is a transient, deleted-after-use download (storage-minimal "
                    "convention) -- not present between sessions, nothing to check here")
    import rasterio

    with rasterio.open(soils_path) as src:
        clay_i = src.descriptions.index("clay_0_30cm_mean") + 1
        sand_i = src.descriptions.index("sand_0_30cm_mean") + 1
        silt_i = src.descriptions.index("silt_0_30cm_mean") + 1
        clay = src.read(clay_i)
        sand = src.read(sand_i)
        silt = src.read(silt_i)
    total = clay + sand + silt
    valid = np.isfinite(total) & (total > 0)
    assert valid.sum() > 0
    median_total = float(np.median(total[valid]))
    assert 900 < median_total < 1100, (
        f"clay+sand+silt median is {median_total:.0f} -- expected ~1000 (per-mille); a value "
        "near 100 would mean the source silently switched to percent"
    )


# --- ISIMIP3b: gitignored, skip if not pulled locally -----------------------------------------

def test_isimip3b_pr_converted_to_mm_per_day():
    """Real bug this guards: ISIMIP3b's pr is a flux (kg m-2 s-1), not an accumulated depth --
    CHELSA-daily's pr (this project's historical reference) is mm/day. scripts/pull_isimip3b.py
    multiplies by 86400 to convert; without it, values are ~1e-4 and silently look like near-
    zero rainfall everywhere rather than a units bug (caught this session because they printed
    as "0.00" at 2 decimal places).
    """
    isimip_dir = DATA_DIR / "isimip3b"
    if not isimip_dir.exists():
        pytest.skip("data/isimip3b/ not present locally")
    pr_files = sorted(isimip_dir.glob("CHELSA_ISIMIP3b_pr__*.npz"))
    if not pr_files:
        pytest.skip("no real pr combos pulled yet in this data/isimip3b/")
    npz = np.load(pr_files[0])
    if "data" not in npz.files:
        pytest.skip(f"{pr_files[0].name} does not have the expected 'data' array (keys: {npz.files})")
    data = npz["data"]
    valid = data[np.isfinite(data)]
    assert valid.size > 0
    mean_val = float(np.nanmean(valid))
    assert mean_val > 0.01, (
        f"mean pr is {mean_val:.6f} -- looks like the raw kg m-2 s-1 flux (~1e-4 to 1e-5 range), "
        "the *86400 mm/day conversion may be missing"
    )
    assert mean_val < 50.0, f"mean pr is {mean_val:.2f} mm/day -- implausibly high, check for a double-conversion"
