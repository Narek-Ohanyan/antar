"""Register the completed CHELSA-BIOCLIM+ pull (scripts/pull_chelsa_bioclim.py) as
manifest entries, using its own manifest_rows.json (real checksums, real per-variable
completeness) rather than re-deriving anything.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from antar.io.manifest import ManifestEntry, write_manifest

OUT_DIR = Path(__file__).resolve().parent.parent / "data" / "chelsa_bioclim"
ACCESS_DATE = "2026-09-30"

CITATION = ("Karger, D.N., Conrad, O., Bohner, J., Kawohl, T., Kreft, H., Soria-Auza, R.W., "
            "Zimmermann, N.E., Linder, H.P., Kessler, M. (2022). Global downscaled projections "
            "for climate impacts studies. Earth System Science Data, 14, 5573-5603. "
            "https://doi.org/10.5194/essd-14-5573-2022")

DESCRIPTIONS = {
    "bio01": "Annual mean temperature (degC). Standard bioclim predictor -- MERISTEM adult-niche Boyce index.",
    "bio02": "Mean diurnal range (degC). MERISTEM niche predictor.",
    "bio03": "Isothermality (bio02/bio07 x 100). MERISTEM niche predictor.",
    "bio04": "Temperature seasonality (std dev x 100). MERISTEM niche predictor.",
    "bio05": "Max temperature of warmest month (degC). MERISTEM niche predictor.",
    "bio06": "Min temperature of coldest month (degC). MERISTEM niche predictor.",
    "bio07": "Temperature annual range (bio05-bio06, degC). MERISTEM niche predictor.",
    "bio08": "Mean temperature of wettest quarter (degC). MERISTEM niche predictor.",
    "bio09": "Mean temperature of driest quarter (degC). MERISTEM niche predictor.",
    "bio10": "Mean temperature of warmest quarter (degC). MERISTEM niche predictor.",
    "bio11": "Mean temperature of coldest quarter (degC). MERISTEM niche predictor.",
    "bio12": "Annual precipitation (mm). MERISTEM niche predictor.",
    "bio13": "Precipitation of wettest month (mm). MERISTEM niche predictor.",
    "bio14": "Precipitation of driest month (mm). MERISTEM niche predictor.",
    "bio15": "Precipitation seasonality (coefficient of variation). MERISTEM niche predictor.",
    "bio16": "Precipitation of wettest quarter (mm). MERISTEM niche predictor.",
    "bio17": "Precipitation of driest quarter (mm). MERISTEM niche predictor.",
    "bio18": "Precipitation of warmest quarter (mm). MERISTEM niche predictor.",
    "bio19": "Precipitation of coldest quarter (mm). MERISTEM niche predictor.",
    "gdd0": "Growing degree days above 0 degC. MERISTEM's GDD modifier cross-check (gdd_modifier).",
    "gdd5": "Growing degree days above 5 degC. MERISTEM's GDD modifier cross-check.",
    "gdd10": "Growing degree days above 10 degC. MERISTEM's GDD modifier cross-check.",
    "gddlgd0": "GDD to last growing degree day, 0 degC base. Growing-season timing for MERISTEM.",
    "gddlgd5": "GDD to last growing degree day, 5 degC base. Growing-season timing for MERISTEM.",
    "gddlgd10": "GDD to last growing degree day, 10 degC base. Growing-season timing for MERISTEM.",
    "gsl": "Growing season length (days), TREELIM methodology. MERISTEM's thermal_treeline_modifier cross-check.",
    "gsp": "Growing season precipitation (mm), TREELIM methodology. MERISTEM treeline cross-check.",
    "gst": "Growing season mean temperature (degC), TREELIM methodology. MERISTEM treeline cross-check.",
    "fcf": "Frost change frequency. MERISTEM's late_frost_modifier cross-check.",
    "fgd": "First growing degree day (day of year). Frost/growing-season timing for MERISTEM.",
    "lgd": "Last growing degree day (day of year). Frost/growing-season timing for MERISTEM.",
    "vpdmean": "Mean vapour pressure deficit (Pa). XYLEM water-stress terms. Historical only -- no future projections exist in this archive (confirmed by directory listing, not a pull failure).",
    "vpdmax": "Max vapour pressure deficit (Pa). XYLEM water-stress terms. Historical only -- same as vpdmean.",
    "petmean": "Mean potential evapotranspiration. TOPOHYDRO water-balance cross-check. Historical only.",
    "petmax": "Max potential evapotranspiration. TOPOHYDRO water-balance cross-check. Historical only.",
    "sfcWindmean": "Mean surface wind speed. ERA5-Land cross-check. Historical only.",
    "rsdsmean": "Mean shortwave radiation. Slope-radiation energy term cross-check. Historical only.",
}


def main():
    rows = json.loads((OUT_DIR / "chelsa_bioclim_manifest_rows.json").read_text())
    entries = []
    for r in rows:
        if r["local_path"] is None:
            print(f"WARNING: {r['variable']} has no local file (zero scenarios found), skipping")
            continue
        historical_only = r["n_missing_404"] == 45 and r["n_scenarios"] == 1
        notes = DESCRIPTIONS.get(r["variable"], "")
        notes += (f" {r['n_scenarios']}/46 scenarios present, {r['n_missing_404']} confirmed-missing (404), "
                  f"{r['n_failed_read']} failed-read.")
        if historical_only:
            notes += " Historical (1981-2010) only; all 45 future GCM x SSP combinations genuinely absent."
        entries.append(ManifestEntry(
            variable=f"chelsa_bioclim_{r['variable']}",
            source="CHELSA-BIOCLIM+",
            version="1.0",
            url="https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/",
            citation=CITATION,
            license="CC0 1.0",
            access_date=ACCESS_DATE,
            local_path=r["local_path"],
            checksum_sha256=r["sha256"],
            spatial_extent="Armenia (clipped)",
            notes=notes,
        ))
    write_manifest(entries, "configs/manifests/chelsa_bioclim.yaml")
    print(f"Wrote {len(entries)} entries to configs/manifests/chelsa_bioclim.yaml")


if __name__ == "__main__":
    main()
