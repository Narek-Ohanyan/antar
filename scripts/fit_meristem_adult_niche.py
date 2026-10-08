"""MERISTEM's adult-niche model (concept note Sec 7/Module D, sec:D): "For each
species we fit a penalised presence-background model on physically meaningful
predictors (CWD, growing degree days, winter minimum, VPD, soil), with the
background drawn from a target-group sample to absorb collection bias, and
score it under spatial-block CV with the continuous Boyce index."

Nothing implementing this existed before this script -- antar.niche.adult had only
the Boyce index validation metric, not the model it validates. Adds
fit_presence_background to antar/niche/adult.py (penalized logistic regression,
the standard practical equivalent of a presence-background point-process/MaxEnt-style
model -- Fithian & Hastie 2013) and uses it here against real data end to end.

PER-SPECIES, not pooled: a first version of this script pooled all seven target
species into one presence class and got a weak, unstable result (mean Boyce index
0.25 across spatial-block folds, one fold strongly negative) -- diagnosed as a
real methodological error, not just noise, once the concept note's own wording
("for each species we fit a...") was checked directly: beech, oak, pine and
juniper have very different climate niches, and pooling them blurs exactly the
signal a niche model is supposed to recover. Refit per species instead.
Pinus kochiana (n=8 real GBIF presence points) is too sparse for a stable 8-
feature penalised fit and is explicitly reported as skipped, not forced.

CWD (climatic water deficit) is not yet a pulled product -- it needs TOPOHYDRO's
full water-balance run (PET ensemble, AET), not yet done against real data. Used
here as CWD_approx = PET - P (annual basis, from CHELSA-BIOCLIM+'s petmean and
bio12) -- a standard simplified definition, not the model's real full AET-based
term, and named accordingly so it is never confused with the real thing once
TOPOHYDRO's water balance exists.
"""
import io
import json
import sys
import time
from pathlib import Path

import numpy as np
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from antar.niche.adult import boyce_index, fit_presence_background
from antar.validation.splits import block_kfold

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
BIOCLIM_DIR = DATA_DIR / "chelsa_bioclim"
OUT_PATH = Path(__file__).resolve().parent.parent / "configs" / "fitted" / "meristem_adult_niche.yaml"
CACHE_PATH = DATA_DIR / "_cache_niche_features.npz"  # avoids re-downloading soils.tif on reruns

BBOX = (43.4, 38.8, 46.7, 41.4)  # lon_min, lat_min, lon_max, lat_max
GRID_SHAPE = (312, 396)
PIXEL_SIZE = 2.6 / 312  # confirmed identical on both axes against the real CHELSA transform
MIN_PRESENCE_FOR_FIT = 25  # below this, an 8-feature penalised fit is not reliable -- report, don't force

TARGET_SPECIES = [
    "Fagus orientalis", "Carpinus betulus", "Quercus macranthera", "Quercus iberica",
    "Pinus kochiana", "Juniperus polycarpos", "Juniperus excelsa",
]

SOILS_DRIVE_FILE_ID = "18LRnI4Nsnlj6ks1UClGDe4BdufqUG24n"  # configs/manifests/gee_exports.yaml
FEATURE_NAMES = ["bio06_winter_min_C", "gdd5", "vpdmean_Pa", "cwd_mm",
                  "clay_0_30cm_mean", "sand_0_30cm_mean", "silt_0_30cm_mean", "soc_0_30cm_mean"]
CWD_COLUMN = FEATURE_NAMES.index("cwd_mm")
REAL_CWD_PATH = DATA_DIR / "_real_cwd_for_meristem.npz"
PRESENCE_OMISSION = 0.10      # the score below which a tenth of a species' own records fall is the presence threshold used when the model is applied (antar.niche.application)
SPECIES_GROUP = {"Fagus orientalis": "mesic_diffuse_porous_broadleaf", "Carpinus betulus": "mesic_diffuse_porous_broadleaf", "Quercus macranthera": "ring_porous_oak",
                 "Quercus iberica": "ring_porous_oak", "Pinus kochiana": "pine", "Juniperus polycarpos": "juniper_arid_conifer", "Juniperus excelsa": "juniper_arid_conifer"}


def latlon_to_rc(lat, lon):
    row = np.floor((BBOX[3] - lat) / PIXEL_SIZE).astype(int)
    col = np.floor((lon - BBOX[0]) / PIXEL_SIZE).astype(int)
    return np.clip(row, 0, GRID_SHAPE[0] - 1), np.clip(col, 0, GRID_SHAPE[1] - 1)


def load_bioclim_historical(var):
    with np.load(BIOCLIM_DIR / f"CHELSA_BIOCLIM_{var}_armenia.npz") as f:
        data = f["data"]
    keys = json.loads((BIOCLIM_DIR / f".{var}_progress.json").read_text())["keys"]
    return data[keys.index("historical")]


def fetch_species_presence(species):
    resp = requests.get("https://api.gbif.org/v1/occurrence/search", params={
        "country": "AM", "scientificName": species, "hasCoordinate": "true", "limit": 300,
    }, timeout=30)
    resp.raise_for_status()
    pts = []
    for r in resp.json().get("results", []):
        lat, lon = r.get("decimalLatitude"), r.get("decimalLongitude")
        if lat is not None and lon is not None and BBOX[0] <= lon <= BBOX[2] and BBOX[1] <= lat <= BBOX[3]:
            pts.append((lat, lon))
    return np.array(pts)


def fetch_target_group_background(exclude_species, n_target=3000):
    """GBIF Plantae occurrences in Armenia, excluding the target species themselves --
    this IS the "target-group sample" the concept note specifies, not a uniform
    random background: collectors who record the target species also tend to
    record other plants on the same trips, so this controls for spatial collection
    bias in a way a uniform background would not.
    """
    records, offset, page = [], 0, 300
    while len(records) < n_target:
        resp = requests.get("https://api.gbif.org/v1/occurrence/search", params={
            "country": "AM", "kingdomKey": 6, "hasCoordinate": "true",
            "limit": page, "offset": offset,
        }, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        for r in data.get("results", []):
            sp = r.get("species")
            lat, lon = r.get("decimalLatitude"), r.get("decimalLongitude")
            if sp in exclude_species or lat is None or lon is None:
                continue
            if not (BBOX[0] <= lon <= BBOX[2] and BBOX[1] <= lat <= BBOX[3]):
                continue
            records.append((lat, lon))
        if data.get("endOfRecords", True):
            break
        offset += page
        time.sleep(0.2)
    return np.array(records[:n_target])


def download_soils_transiently():
    import ee
    from google.oauth2.credentials import Credentials
    import google.auth.transport.requests as gareq
    from googleapiclient.discovery import build
    from googleapiclient.http import MediaIoBaseDownload

    creds_path = Path.home() / ".config" / "earthengine" / "credentials"
    d = json.loads(creds_path.read_text())
    creds = Credentials(
        None, refresh_token=d["refresh_token"], token_uri="https://oauth2.googleapis.com/token",
        client_id=ee.oauth.CLIENT_ID, client_secret=ee.oauth.CLIENT_SECRET, scopes=d["scopes"],
    )
    creds.refresh(gareq.Request())
    drive = build("drive", "v3", credentials=creds)
    request = drive.files().get_media(fileId=SOILS_DRIVE_FILE_ID)
    local_path = DATA_DIR / "_tmp_soils.tif"
    fh = io.FileIO(local_path, "wb")
    downloader = MediaIoBaseDownload(fh, request)
    done = False
    while not done:
        _, done = downloader.next_chunk()
    fh.close()
    return local_path


def extract_features(ll, bioclim, soils_path):
    import rasterio
    from rasterio.warp import transform as warp_transform
    rows, cols = latlon_to_rc(ll[:, 0], ll[:, 1])
    feats = np.column_stack([bioclim[name][rows, cols] for name in ("bio06", "gdd5", "vpdmean", "cwd_approx")])
    with rasterio.open(soils_path) as src:
        xs, ys = warp_transform("EPSG:4326", src.crs, ll[:, 1].tolist(), ll[:, 0].tolist())
        soil_vals = np.array(list(src.sample(zip(xs, ys))))
    return np.column_stack([feats, soil_vals])


def build_or_load_features():
    """Returns dict: species -> (ll, X); plus 'background' -> (ll, X). Cached across reruns
    so a script bug downstream of feature extraction doesn't force re-downloading soils.tif."""
    if CACHE_PATH.exists():
        print(f"Loading cached features from {CACHE_PATH}", flush=True)
        cached = np.load(CACHE_PATH, allow_pickle=True)
        return {k: (cached[f"{k}_ll"], cached[f"{k}_X"]) for k in cached["keys"]}

    print("=== Fetching real GBIF presence records (per species) ===", flush=True)
    presence = {sp: fetch_species_presence(sp) for sp in TARGET_SPECIES}
    for sp, ll in presence.items():
        print(f"  {sp}: {len(ll)} real presence points", flush=True)

    print("=== Fetching target-group GBIF background (not uniform random) ===", flush=True)
    background_ll = fetch_target_group_background(set(TARGET_SPECIES))
    print(f"  {len(background_ll)} real target-group background points", flush=True)

    print("=== Loading CHELSA-BIOCLIM+ historical predictors ===", flush=True)
    bioclim = {name: load_bioclim_historical(name) for name in ("bio06", "gdd5", "vpdmean", "petmean", "bio12")}
    bioclim["cwd_approx"] = bioclim["petmean"] - bioclim["bio12"]

    print("=== Downloading SoilGrids transiently (deleted after extraction) ===", flush=True)
    soils_path = download_soils_transiently()

    result = {}
    for sp, ll in presence.items():
        result[sp] = (ll, extract_features(ll, bioclim, soils_path) if len(ll) else np.empty((0, 8)))
    result["background"] = (background_ll, extract_features(background_ll, bioclim, soils_path))
    soils_path.unlink()  # storage-minimal: local copy deleted immediately after extraction

    save_kwargs = {"keys": np.array(list(result.keys()), dtype=object)}
    for k, (ll, X) in result.items():
        save_kwargs[f"{k}_ll"] = ll
        save_kwargs[f"{k}_X"] = X
    np.savez_compressed(CACHE_PATH, **save_kwargs)
    return result


def with_real_cwd(features):
    """Replace the climatic-water-deficit column (the proxy ``petmean - bio12`` in the cache) by the real water-balance CWD of each point for the species' own functional group
    (``compute_real_cwd_for_meristem.py``). The cache holds the points in a fixed order (each species' records, then the shared background); the real-CWD file keeps that order."""
    z = np.load(REAL_CWD_PATH, allow_pickle=True)
    labels = list(z["groups"])
    out = {}
    order = list(TARGET_SPECIES) + ["background"]
    offsets, start = {}, 0
    for k in order:
        n = len(features[k][0])
        offsets[k] = (start, start + n)
        start += n
    if start != len(labels):
        raise ValueError(f"the real-CWD file has {len(labels)} points but the feature cache has {start}")
    for k in order:
        a, b = offsets[k]
        if set(labels[a:b]) != {k} and b > a:
            raise ValueError(f"the real-CWD file's points {a}:{b} are not '{k}' records")
    for sp in TARGET_SPECIES:
        ll, X = features[sp]
        a, b = offsets[sp]
        X = X.copy()
        if len(X):
            X[:, CWD_COLUMN] = z[f"real_cwd_mm__{SPECIES_GROUP[sp]}"][a:b]
        bll, bX = features["background"]
        a, b = offsets["background"]
        bX = bX.copy()
        bX[:, CWD_COLUMN] = z[f"real_cwd_mm__{SPECIES_GROUP[sp]}"][a:b]
        out[sp] = (ll, X, bll, bX)
    return out


def fit_one_species(species, presence_ll, presence_X, background_ll, background_X):
    valid_p = np.all(np.isfinite(presence_X), axis=1)
    presence_ll, presence_X = presence_ll[valid_p], presence_X[valid_p]
    valid_b = np.all(np.isfinite(background_X), axis=1)
    background_ll, background_X = background_ll[valid_b], background_X[valid_b]

    if len(presence_X) < MIN_PRESENCE_FOR_FIT:
        return {
            "status": "skipped_insufficient_data",
            "n_presence_after_nan_drop": int(len(presence_X)),
            "reason": f"fewer than {MIN_PRESENCE_FOR_FIT} valid presence points -- an "
                      f"{len(FEATURE_NAMES)}-feature penalised fit would not be reliable "
                      f"(risk of near-perfect separation / unstable coefficients), not forced.",
        }

    all_ll = np.vstack([presence_ll, background_ll])
    all_X = np.vstack([presence_X, background_X])
    all_y = np.concatenate([np.ones(len(presence_X)), np.zeros(len(background_X))])
    block_size = 0.25  # degrees, roughly 25-28 km at this latitude
    boyce_scores = []
    for fold, (tr, te) in enumerate(block_kfold(all_ll[:, 1], all_ll[:, 0], block_size, k=5, seed=0)):
        if all_y[tr].sum() < 5 or (all_y[tr] == 0).sum() < 5:
            continue
        model = fit_presence_background(all_X[tr][all_y[tr] == 1], all_X[tr][all_y[tr] == 0])
        pred_all = model.predict_proba(all_X[te])[:, 1]
        pred_presence = pred_all[all_y[te] == 1]
        if len(pred_presence) < 3:
            continue
        cbi = boyce_index(pred_presence, pred_all)
        boyce_scores.append(None if np.isnan(cbi) else round(float(cbi), 3))

    final_model = fit_presence_background(presence_X, background_X)
    valid_scores = [s for s in boyce_scores if s is not None]
    presence_scores = final_model.predict_proba(presence_X)[:, 1]
    threshold = float(np.quantile(presence_scores, PRESENCE_OMISSION))
    return {
        "status": "fitted",
        "n_presence_after_nan_drop": int(len(presence_X)),
        "n_background_after_nan_drop": int(len(background_X)),
        "spatial_block_cv": {
            "block_size_deg": block_size, "k": 5,
            "boyce_index_per_fold": boyce_scores,
            "boyce_index_mean": round(float(np.mean(valid_scores)), 3) if valid_scores else None,
        },
        "coefficients": {n: float(c) for n, c in zip(FEATURE_NAMES, final_model.coef_[0])},
        "intercept": float(final_model.intercept_[0]),
        "presence_threshold": {"omission": PRESENCE_OMISSION, "score": threshold, "n_records_below": int((presence_scores < threshold).sum())},
    }


def main():
    features = build_or_load_features()
    per_species = with_real_cwd(features)
    background_ll = features["background"][0]
    print(f"\n=== Fitting per-species models against {len(background_ll)} shared background points (real water-balance CWD of each species' group) ===", flush=True)

    results = {}
    for sp in TARGET_SPECIES:
        presence_ll, presence_X, bg_ll, bg_X = per_species[sp]
        print(f"\n--- {sp} ({len(presence_ll)} raw presence points) ---", flush=True)
        r = fit_one_species(sp, presence_ll, presence_X, bg_ll, bg_X)
        results[sp] = r
        if r["status"] == "fitted":
            print(f"  Boyce index per fold: {r['spatial_block_cv']['boyce_index_per_fold']}", flush=True)
            print(f"  Boyce index mean: {r['spatial_block_cv']['boyce_index_mean']}", flush=True)
        else:
            print(f"  {r['reason']}", flush=True)

    output = {
        "fit_date": __import__("datetime").date.today().isoformat(),
        "method": "Per-species penalised (L2) presence-background logistic regression "
                  "(antar.niche.adult.fit_presence_background), spatial-block CV (5-fold, "
                  "0.25deg blocks), scored with the continuous Boyce index.",
        "background_type": "GBIF target-group (Plantae, Armenia, excluding target species) -- "
                            "not uniform random, shared across all species' fits.",
        "n_background_total": int(len(background_ll)),
        "features": FEATURE_NAMES,
        "cwd_caveat": ("cwd_mm is the real water-balance climatic water deficit (TOPOHYDRO, PM-FAO56) for 2019 at each record, computed with the species group's own rooting depth (Canadell et al. 1996). "
                       "The other climate predictors are 1981-2010 CHELSA-BIOCLIM+ normals, so one year of water deficit sits beside climatological predictors; 2019 is the project's single reference year."),
        "pooling_note": "An earlier version of this fit pooled all 7 species into one presence "
                        "class and got a weak, unstable result (mean Boyce 0.25, one fold "
                        "negative) -- diagnosed as pooling ecologically disparate species "
                        "(beech/oak/pine/juniper have very different niches), not noise; refit "
                        "per species to match the concept note's actual specification.",
        "species": results,
    }
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    import yaml
    output["species_group"] = SPECIES_GROUP
    OUT_PATH.write_text(yaml.safe_dump(output, sort_keys=False))
    print(f"\nWrote {OUT_PATH}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
