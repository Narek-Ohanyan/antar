# Roadmap to a working framework

A living checklist, not a log (see `IMPLEMENTATION_LOG.md` for the append-only history of
decisions and bugs). Update items in place as they close; don't duplicate them into the log
unless something non-obvious was learned closing them.

Every entry below marked **CONFIRMED** was checked directly against the provider's own
catalog/API/documentation this session. Everything marked **UNVERIFIED** was found by search
but not independently confirmed (a tool outage cut several checks short) — verify before
depending on it, the same way the CHELSA/CORDEX mix-up and the GLO-30 coverage gap were caught
by checking the actual data rather than trusting the name of a dataset.

## 0. The mandatory cross-cutting constraint: consistency before coverage

Every new source below must clear this bar before it's wired in — adding a source that isn't
harmonised is worse than not having it, because it produces numbers that look plausible and
aren't.

- [ ] **RCP vs SSP scenario correspondence.** We now have RCP-labelled data (CORDEX on hyperion:
      RCP2.6/RCP8.5) and SSP-labelled data (CHELSA-BIOCLIM+, NASA/GDDP-CMIP6, ISIMIP3b: all
      SSP1-2.6/SSP3-7.0/SSP5-8.5). RCPs (CMIP5 forcing pathways) and SSPs (CMIP6
      socioeconomic-plus-forcing pathways) are not the same framework and are not
      interchangeable by number alone (RCP2.6 and SSP1-2.6 target similar end-of-century forcing
      but via different socioeconomic assumptions; treating them as identical would be exactly
      the kind of silent unit-like error this project's conventions exist to prevent). Before any
      SSP-sourced data is mixed with the RCP-sourced CORDEX data in the same ensemble or the same
      figure, write down the explicit correspondence being assumed (or don't mix them — run them
      as separate, labelled ensemble members).
- [ ] **Resolution reconciliation, documented per source, not assumed:** CHELSA-daily (1 km,
      daily) vs TerraClimate (~4.6 km, monthly) vs NASA/GDDP-CMIP6 (0.25°, ~28 km, daily) vs
      ISIMIP3b (0.5°, daily) vs CORDEX WAS-22 (already resampled on hyperion) vs SoilGrids
      (250 m) vs the master grid (30 m, EPSG:32638). Every one of these needs an explicit
      regridding/aggregation method to the master grid, stated in `configs/`, not left implicit
      in whatever a script happens to do.
- [ ] **Unit audit per new variable**: TerraClimate's `tmmn`/`tmmx` are °C already (confirm, not
      CHELSA's typical scaled-integer K-offset encoding); GHCN-Daily is tenths of °C / tenths of
      mm in the raw files (a classic silent-error source — the v1 audit this project's own
      `scripts/audit_v1.py` exists to catch was exactly this class of mistake); CO2 pathways are
      ppm, not GtC or W/m² forcing. Add a unit-bounds test for each new variable the same way
      existing TOPOHYDRO/XYLEM variables are tested (`ΔVPD/ΔT_max` physical-bound pattern).
- [ ] **One forest/treeline definition, one reference period** — already a stated convention
      (`reference_period: [1991, 2020]`); any new historical source (TerraClimate, GHCN-Daily)
      must be checked against that period's actual coverage before being used for anomalies, the
      same check that caught the first CHELSA-daily pull's 2000-2024 gap.

## 1. Data — real sources confirmed this session, not yet wired in

### 1.1 Future climate projections (highest-value finding)
- [ ] **NASA/GDDP-CMIP6 on Earth Engine** — **CONFIRMED**. Asset `ee.ImageCollection("NASA/GDDP-CMIP6")`
      (not "NEX-GDDP-CMIP6" — the ID changed). Daily pr/tas/tasmax/tasmin/hurs/huss/rlds/rsds/
      sfcWind, 1950-2100, 0.25°, 34 GCMs, already bias-corrected (BCSD). Free EE account, no
      institutional approval. **UNVERIFIED**: EE catalog page only documents historical/ssp245/
      ssp585 — confirm ssp126/ssp370 are actually present in this ingested collection (NASA's own
      archive has all four) before relying on it in the Code Editor. Full-SSP fallback with zero
      registration: `s3://nex-gddp-cmip6` (AWS Open Data, NetCDF, `--no-sign-request`).
- [ ] **ISIMIP3b direct access** — **CONFIRMED** real and accessible: `data.isimip.org`, 5 core
      GCMs (GFDL-ESM4, IPSL-CM6A-LR, MPI-ESM1-2-HR, MRI-ESM2-0, UKESM1-0-LL) × ssp126/ssp370/
      ssp585 + historical, bias-adjusted against W5E5. No account needed for public downloads
      (accounts are for DKRZ HPC / protocol participants only). Recommended mechanism: official
      `isimip-client` Python package, `cutout_bbox(west, east, south, north)` for a server-side
      regional cutout. **UNVERIFIED**: the portal sits behind an "Anubis" anti-bot gate that
      blocked a direct fetch this session — run a live `pip install isimip-client` +
      `cutout_bbox` test before committing the pipeline to this path.
- [ ] **Decision needed**: given both of the above are real, decide whether to (a) replace the
      CORDEX/CMIP5 hyperion ensemble with the concept note's original ISIMIP3b + NEX-GDDP-CMIP6
      design now that it's confirmed reachable, (b) run both as separate labelled ensembles (more
      work, more honest about structural uncertainty), or (c) keep CORDEX as primary and use
      NASA/GDDP-CMIP6 only as a cross-check. This is a scientific-design decision, not just a
      data-access one — make it deliberately, not by default.

### 1.2 Historical cross-checks and validation data
- [ ] **TerraClimate** — **CONFIRMED**, `IDAHO_EPSCOR/TERRACLIMATE` on Earth Engine (matches the
      concept note's own spec as an independent water-balance benchmark). 14 bands incl. pr,
      tmmn/tmmx, pet, aet, soil moisture, runoff, PDSI, VPD. 1958-2024 monthly, ~4.6 km, CC0.
      Direct-download alternative (no EE account): `climatologylab.org/terraclimate.html`
      (THREDDS/OPeNDAP/wget).
- [x] **GHCN-Daily Armenian stations — pulled and registered (2026-09-30)**: all 53 real
      stations, 1,325,996 QC-passed observations, PRCP/TMAX/TMIN/TAVG/SNOW/SNWD, units
      converted and validated against known physical values (Yerevan July ~26-30°C, January
      ~0.7-6°C) before trusting the pull at scale. Drive-hosted, `configs/manifests/
      ghcnd_armenia.yaml`. **This replaces "Armhydromet stations are institution-only" from the
      earlier manual-data list — they were not.** ECA&D also covers Armenia but far more thinly
      (2-3 stations) — not pulled, supplementary only if ever needed.

### 1.3 Terrain and accessibility
- [x] **Copernicus GLO-30 Armenia gap — resolved definitively (2026-09-30), SRTM confirmed
      correct.** The "licensing exclusion lifted ~Dec 2023" theory was checked directly against
      the authoritative source itself (Copernicus Data Space Ecosystem's own STAC API,
      `stac.dataspace.copernicus.eu`, collection `cop-dem-glo-30-dged-cog`) and against the AWS
      Open Data mirror (`copernicus-dem-30m`, tile-by-tile listing) — both return the exact same
      3 tiles as Earth Engine's copy (`N38/E043`, `N38/E044`, `N39/E043`), covering only the
      southwest corner. This is not a stale-mirror artifact: **Armenia genuinely has zero
      Copernicus GLO-30 coverage anywhere in the Copernicus system as of today.** The earlier
      "should have full coverage now" note (from the initial research pass) was wrong — corrected
      here after checking the primary source directly rather than trusting it. SRTM 30 m stays
      the right substitute; no further action needed on this item.
- [x] **Road network / accessibility — pulled (2026-09-30), via the Earth Engine path.**
      Geofabrik's binary downloads (`.shp.zip`, `.osm.pbf`) hit a persistent redirect loop from
      this environment when actually attempted — real server/network behavior, not a data
      problem (the HTML listing page itself loads fine) — not pursued further once the EE
      alternative confirmed working. Used the Malaria Atlas Project's global travel-time-to-
      cities layer instead: confirmed real, asset id `Oxford/MAP/accessibility_to_cities_2015_v1_0`
      (verified directly, single `accessibility` band, minutes-to-nearest-city), Armenia bbox
      values 0-517 minutes (physically sensible for the terrain). Exported to Drive
      (`antar_accessibility_to_cities`, task `3GOANAHUA6I5AINU42GBPF6M`, queued under restricted
      mode as of this writing — register in `configs/manifests/` once it actually completes).
      OSM via Geofabrik remains a real fallback if finer road-network detail (vs. travel-time)
      is ever needed — retry from a different network/environment, or use the Overpass API
      instead of the static mirror.

### 1.4 Bioclimate and occurrence
- [x] **CHELSA-BIOCLIM+ — access confirmed (2026-09-30), not yet pulled.** `bioclim/` is a real
      sibling directory to `daily/` on the same server (`os.unil.cloud.switch.ch/chelsa02/chelsa/
      global/bioclim/`), same `/vsicurl/` access pattern, listed directly: bio01-19, cltmax/mean/
      min/range, cmimax/mean/min/range, fcf, fgd, gdd0/5/10 and more. This is a properly scoped
      pull job on its own (many variables x historical + 5 GCMs x SSP126/370/585 x 3 future
      windows) — deliberately not rushed alongside the smaller items in this batch; next real
      data-pull task once scoped.
- [x] **GBIF occurrence data — pulled (2026-09-30).** 524 real presence records for the 7 target
      species (`configs/species_traits.csv`'s example taxa), `configs/manifests/
      gbif_occurrences.yaml`. A real data-quality bug was caught: `country=AM` alone let through
      4 `(0,0)` "null island" bad-georeferencing records; fixed with an explicit bbox filter.

### 1.5 CO2 concentration pathways
- [x] **RCMIP / Meinshausen et al. (2020) — pulled (2026-09-30).** Real, direct, no registration
      (`zenodo.org/record/4589756/files/rcmip-concentrations-annual-means-v5-1-0.csv`, 21 MB
      confirmed downloadable). Filtered to exactly this project's scenarios — both RCP26/45/60/85
      (CORDEX) and SSP126/245/370/585 (everything else) plus historical — into one small table,
      `configs/co2_concentration_pathways.csv`, `configs/manifests/co2_pathways.yaml`. Real
      numbers confirm RCP2.6 and SSP1-2.6 are close but not identical by 2100 (~421 vs. ~446 ppm)
      — the §0 correspondence risk is real, not hypothetical.

### 1.6 Confirmed genuinely not obtainable as open data
- [ ] Tree-ring chronologies for Armenia — **CONFIRMED absent from ITRDB**; the two real, recent
      papers found (Opała-Owczarek et al. 2021, *Atmosphere*; Stepanyan et al. 2026, *Ecologies*)
      both state their raw data is available "on request from the corresponding author" only.
      This is now a concrete task, not a vague one: **email the corresponding authors of both
      papers** if the raw ring-width series matter for growth-climate calibration. ("Voss
      dendrochronology Armenia," referenced earlier as a citation, could not be located — recheck
      the author name/spelling against the concept note's own bibliography.)
  - [ ] National Forest Monitoring/Inventory plots, provenance/genetic trial locations,
        insect/pathogen outbreak records, treeline field-survey transects, and prior-project
        plantation survey data (survival/height by planting year and method) — still
        institution/field-survey-only, unchanged from the earlier assessment.
  - [ ] TRY and XFT trait databases — real, but each needs your own free registration; not an
        access-method problem, just needs doing.

## 2. Data — already submitted but never verified/finished

- [x] **Six-plus Earth Engine exports, checked and registered (2026-09-30)**: with the real
      GCP project ID, `terrain`, `soils`, `era5land_forcing`, `snow`, `land_tenure`,
      `vitality_composites`, `landtrendr_segmentation`, `structure`, and
      `disturbance_ancillary` were all confirmed COMPLETED and registered in
      `configs/manifests/gee_exports.yaml` (38 entries, ~85 GB, left on Drive by design —
      cloud-first, nothing downloaded). `vegetation_state` had genuinely FAILED (`Image.select:
      Parameter 'input' is required and may not be null`) — real bug: MCD12Q2 land-surface
      phenology only covers 2001 onward, and the export loop's 2000 start hit `.first()` on an
      empty collection. Fixed (`src/antar/io/gee_export.py`: phenology bands genuinely omitted
      for years outside real coverage, not fabricated) and resubmitted — **your GCP project is
      in Earth-Engine restricted mode (noncommercial compute quota exceeded)**, so the resubmit
      is queued (`READY`), not yet run; check `ee.data.getTaskList()` again once it clears.
- [x] **GHCN-Daily Armenian stations pulled and registered (2026-09-30)** — see §1.2 above;
      `configs/manifests/ghcnd_armenia.yaml`, Drive-hosted, 5 MB.
- [ ] **CORDEX future extraction on hyperion**: only 2 of 4 variables (`pr`, `tas`) done for 1 of
      8 GCM×RCP combinations (MOHC-HadGEM2-ES/rcp26). `tasrange`/`tasskew` for that combination,
      plus all 7 remaining combinations, still pending — blocked on hyperion network access.
- [ ] Future `tasmax`/`tasmin` derivation from CORDEX `tas`+`tasrange` (needed only if CORDEX
      stays in the ensemble per the §1.1 decision) — flagged, not implemented.

## 3. Model — per engine, real fitting against real data

Nothing below has been fit against real data yet; every engine is implemented and tested against
synthetic/placeholder data only (`configs/species_traits.csv` is entirely
`ILLUSTRATIVE_PLACEHOLDER`, no fit/calibration script exists in the repo yet).

- [ ] **TOPOHYDRO**: fit the monthly lapse rate (`fit_monthly_lapse_rate`) and cold-air-pooling
      coefficient against real station data (now unblocked — GHCN-Daily, §1.2) and CHELSA-daily;
      fit the precipitation-elevation gradient the same way; run the PET ensemble and QDM against
      real ERA5-Land once §2's export is confirmed.
- [ ] **XYLEM**: replace `species_traits.csv` placeholders with real trait priors (TRY/XFT once
      registered, §1.6) for the four functional groups; re-run the Monte Carlo failure engine and
      re-fit the monotone emulator against real trait draws; re-run the emulator's held-out
      release gate with real data, not synthetic.
- [ ] **MNEME**: build the real person-period panel from the vitality/disturbance exports (§2);
      fit the cloglog GLM + monotone GBM + NNLS stack against real dieback labels; run the actual
      nested blocked-CV validation designed in `configs/cv.yaml`.
- [ ] **MERISTEM**: fit the adult-niche Boyce index against real GBIF occurrence data (§1.4,
      genus-filtered) and the attainable-height quantile model against real GEDI/canopy-height
      data (once §2's structure export is confirmed); calibrate growth modifiers against whatever
      real growth-climate data ends up available (§1.6's tree-ring email outreach, or the
      prior-project plantation survey data if obtained).
- [ ] **REFUGIUM**: run the full robust-refugium criteria against real MERISTEM/hazard outputs
      once both feed it real fits, not placeholders.
- [ ] **AEGIS**: run the CVaR portfolio optimisation and efficient frontier against a real
      scenario/model/parameter ensemble once §1.1's projection-source decision is made.
- [ ] **Pipeline wiring**: `workflow/Snakefile` is still 100% placeholder `echo` commands for
      every rule — replace each with the real `antar.*` entry point once its inputs are real,
      so the pipeline can actually run end-to-end rather than existing only as tested library code.
- [ ] **Dev environment**: no `.venv` with the full `.[dev]` extras currently exists in this
      checkout (only the narrow `.venv_chelsa_pull`) — recreate it (`pip install -e ".[dev]"`)
      and re-run the full test suite to confirm it's still green before real fitting starts.
- [ ] **CO2 physiological effects**: once §1.5's CO2 pathway table is in hand, wire it into
      whatever growth/water-use-efficiency term in MERISTEM/XYLEM is meant to use it (check the
      concept note for where CO2 fertilisation or WUE scaling is specified, if at all — flagged
      here as unconfirmed whether the concept note even calls for this mechanistically, rather
      than assuming it does).

## 4. Validation

- [ ] Design-based accuracy assessment against a stratified validation sample (concept note's own
      stated minimum: ≥300 plots) — blocked entirely on real ground-truth data (§1.6).
  Structural gap or not: the pipeline's `validation_report` rule (AOA maps, calibration,
  model card) is also still an `echo` placeholder in `workflow/Snakefile` — same §3 item.

## 5. Interface — deferred until real fits exist, scoped for when that decision is made

Not started by design (explicitly deferred pending real model fits — see
`IMPLEMENTATION_LOG.md`). Scope as you described it, for when the go-ahead is given:

- [ ] Home page.
- [ ] Interactive map with buttons/dropdowns to select what's being predicted (refugia score,
      viability, hazard, attainable height, etc.), scenario/GCM/RCP-or-SSP member, and horizon
      (2050/2080/2100).
- [ ] Time-series and graph views per selected site/pixel (climate trajectory, hazard components,
      growth curve).
- [ ] Detailed methodology page — needs to be the *actual* methodology as implemented (which
      differs from the concept note in several logged, deliberate ways: CORDEX/CMIP5 vs.
      ISIMIP3b/NEX-GDDP-CMIP6 per §1.1's decision, SRTM vs. GLO-30 per §1.3, the treeline/GDD
      modifier split, the Eq. 8.6 conjunction implementation, etc.) — not a restatement of the
      concept note.
  - [ ] A plain-language toggle/button on this page ("explain it simply" or similar) that
        switches each formula/mechanism to a non-academic explanation built on concrete,
        everyday analogies (e.g. a tree's water column under drought stress explained the way
        you'd explain suction through a straw, or the two-phase hydraulic-failure threshold
        explained the way you'd explain a dam overtopping) — no jargon, no equations, one
        analogy per mechanism. Runs alongside the technical version, not instead of it: every
        engine's page needs both a formula/equation view and a plain-language view, switchable
        without leaving the page. Content has to be written per-mechanism once the real
        methodology text exists (§3) — this is a content task as much as a UI one, not just a
        toggle to build.
- [ ] References page — every dataset in this roadmap and `IMPLEMENTATION_LOG.md`'s manifests,
      with its actual required citation (already recorded per-entry in `configs/manifests/*.yaml`
      — this page can likely be generated from those rather than hand-written).
- [ ] Acknowledgments — the placeholder already left in `README.md` for the prior project's own
      attribution, to be filled in by hand.
- [ ] Stack decision (not yet made): what serves the map/API — needs the real gridded outputs to
      exist first (§3) before this can be scoped concretely rather than guessed at.
