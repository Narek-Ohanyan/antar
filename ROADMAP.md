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
      values 0-517 minutes (physically sensible for the terrain). Original export queued
      indefinitely under the first GCP project's restricted mode; cancelled and resubmitted on
      the second project (`antar-armenia-2`, set up specifically for this), completed and
      registered in `configs/manifests/gee_exports.yaml`.
      OSM via Geofabrik remains a real fallback if finer road-network detail (vs. travel-time)
      is ever needed — retry from a different network/environment, or use the Overpass API
      instead of the static mirror.

### 1.4 Bioclimate and occurrence
- [x] **CHELSA-BIOCLIM+ — pulled (2026-09-30).** 37 of the archive's 74 real variables (confirmed
      by listing the bucket), scoped to what an ANTAR engine actually uses — bio01-19 (MERISTEM
      niche), gdd0/5/10 + gddlgd0/5/10 (MERISTEM GDD modifier), gsl/gsp/gst (MERISTEM treeline,
      TREELIM methodology), fcf/fgd/lgd (MERISTEM late-frost modifier), vpdmean/max (XYLEM),
      petmean/max (TOPOHYDRO), sfcWindmean, rsdsmean (ERA5-Land/radiation cross-checks) —
      historical (1981-2010) + 3 future windows x 5 GCMs x 3 SSPs where available. A real nodata
      bug (int32 sentinel leaking through unmasked) was caught on the first smoke test and fixed
      before the full pull ran; see IMPLEMENTATION_LOG.md. 6 of the 37 variables (vpdmean, vpdmax,
      petmean, petmax, sfcWindmean, rsdsmean) turned out historical-only in the archive itself —
      confirmed by directory listing, not a pull failure. Registered in `configs/manifests/
      chelsa_bioclim.yaml`, 223 MB, kept local (small enough that the cloud-first convention
      doesn't apply the way it does to the multi-GB GEE rasters).
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

### 1.6 Ground-truth and traits — real findings, not a blanket "institution-only"

A dedicated search (2026-09-30, not done in the original research pass) overturned part of the
earlier "institution-only" assessment. Checked live, not from search snippets alone:

- [x] **Skipped by user decision**: tree-ring chronologies. Confirmed absent from ITRDB; the two
      real papers found (Opała-Owczarek et al. 2021, *Atmosphere*; Stepanyan et al. 2026,
      *Ecologies*) both say "on request from the corresponding author" only. User chose not to
      pursue the outreach.
- [ ] **Real, pullable, previously missed: `data.opendata.am` (Data Catalog Armenia)**, a genuine
      CKAN open-data portal, no account needed, 173 forest-related datasets. Most useful for
      ANTAR: **forest degradation by forestry branch (Artsvaberd, Ijevan, Eghegnut), 2016-2020**
      (REC Caucasus, 6 resources each) — a real disturbance-label complement to Hansen GFC/MODIS
      burned area for MNEME; a general forest-cover layer already in **EPSG:32638** (matches the
      master grid exactly); and Armenian National Agrarian University soil property layers
      (clay/silt/sand/bulk density/CEC/pH) as a possible SoilGrids cross-check. Not yet pulled —
      next concrete task.
- [ ] **Real, pullable: EPPO pest/pathogen occurrence records** for Armenia — confirmed live,
      `gd.eppo.int/country/AM/organisms`, CSV/Excel export, no login. Not yet pulled.
- [x] **Confirmed genuinely absent, not a search failure**: EUFGIS provenance/genetic
      conservation units — Armenia has zero registered units (checked live at eufgis.org).
- [ ] **Treeline field-survey transects**: no direct public download found. Closest match is the
      Transcaucasian Vegetation Database (2,882 Braun-Blanquet plots, Zenodo record 10412100) but
      only descriptive metadata is there — raw plot access would need contacting the database's
      custodians directly, a real but manual next step if this matters enough to pursue.
  - [ ] National Forest Inventory plots specifically (not the degradation/cover layers above) and
        prior-project plantation survey data (survival/height by planting year and method) —
        confirmed still institution/field-survey-only: Armenia has had no full NFI since the
        1980s, a new FAO-supported one is mid-first-cycle with data not yet released, and the
        Ministry of Environment's ArcGIS Hub requires sign-in despite search results suggesting
        otherwise (checked live, not assumed from a snippet).
- [ ] **TRY and XFT trait databases — registration process confirmed, still needs you to do it.**
      TRY: `try-db.org/TryWeb/RegStart.php`, free account (institution field is optional), no
      approval wait; published trait data is CC BY/open access since v5. XFT: real standalone site
      at `xylemfunctionaltraits.org` (the Choat et al. 2012 hydraulic trait dataset), no
      registration barrier found, also cross-hosted inside TRY's archive.

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
      for years outside real coverage, not fabricated), resubmitted on the second GCP project
      (`antar-armenia-2`) after the first project's restricted mode left it queued indefinitely,
      completed, and registered.
- [x] **Second GCP project (`antar-armenia-2`) set up for fresh Earth Engine quota (2026-09-30)**:
      the first project (`pure-highlander-495708-a9`) hit its noncommercial compute quota and
      entered restricted mode, throttling exports to indefinite `READY` queueing. Rather than
      wait it out, created and registered a second project — both stuck exports (`vegetation_
      state`, `accessibility_to_cities`) moved there and completed within minutes. Use this
      project for future EE work; fall back to creating a third if this one also fills up.
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

- [x] **TOPOHYDRO monthly lapse rate + precip gradient — real fit, done (2026-09-30).** The
      first real fit anywhere in ANTAR: `scripts/fit_topohydro_lapse_rate.py`, real GHCN-Daily
      station data (41 of 53 stations with observations in the 1991-2020 reference period),
      45,665 station-days of temperature, 50,607 of precipitation. Results physically checked,
      not just "the fit ran": lapse rate -4.9 to -7.6 K/km with the correct seasonal pattern
      (steeper in summer, shallower in winter), R²=0.32-0.70. Precipitation's elevation gradient
      came back with R²≈0 every month — a real, informative finding (Armenia's precipitation
      isn't elevation-driven the way temperature is), not a bug; confirms the `exposure_factor`
      correction term will matter more than the base gradient. Saved to
      `configs/fitted/topohydro_lapse_rate.yaml`.
  - [x] **Saxton & Rawls (2006) soil pedotransfer — built, done (2026-09-30).**
        `src/antar/climate/soil_pedotransfer.py`: SoilGrids clay/sand/SOC -> theta_sat/
        theta_fc/theta_lim/psi_sat_mpa/b_clapp_hornberger, the parameters
        `topoclimate_forcing` needs and nothing previously computed. Coefficients read
        directly from a source reproducing Saxton & Rawls's own Table 1 (the primary PDF
        was down; a wrong "2006" link was caught resolving to the different 1986 paper
        before being trusted) rather than recalled from memory — a wrong pedotransfer
        constant would silently corrupt every cell's water balance. Cross-checked against
        Clapp & Hornberger's own 1978 b-value table for loam (5.23 fit vs. 5.39 published),
        not just internal self-consistency. 4 new tests, full suite 147 passed.
  - [ ] **Still needed**: the cold-air-pooling coefficient (`k_cap`) — needs a terrain concavity
        index that isn't computed yet (flagged since `export_terrain`: needs neighbourhood/flow-
        routing operations on the downloaded SRTM DEM, not a per-pixel Earth Engine operation).
        Also: run the PET ensemble and QDM against real ERA5-Land once that export (§2, already
        registered) is actually loaded and checked. A real gridded `topoclimate_forcing` run
        still needs ERA5-Land loaded, terrain slope/aspect/concavity computed, and a per-cell
        orchestration script — the pedotransfer piece above unblocks but doesn't complete this.
- [x] **XYLEM P50/slope/gmin — real XFT data, done (2026-09-30).** `configs/species_traits.csv`'s
      p50_mpa and slope_pct_per_mpa are now real values from the XFT database (congeneric
      proxies for all four groups — no exact-species match exists for any target taxon; see
      IMPLEMENTATION_LOG.md for the full per-group species/sample-size provenance, including a
      real mid-analysis correction on the juniper group). gmin25 updated from real data only
      where XFT had it (pine). psi_close_mpa, tp_c, lethal_plc, capacitance_mmol_m2_mpa remain
      placeholders — no clean, unit-compatible XFT field exists for any of them.
  - [x] **hyper_sd (trait uncertainty) — real, computed from XFT's actual spread (2026-09-30).**
        `configs/fitted/xylem_trait_hyper_sd.yaml`: real between-record P50/slope standard
        deviations per group (broadleaf 0.84 MPa, oak 1.10, pine 0.61; juniper undefined,
        n=1). Previously every call site only ever used an illustrative guessed value.
  - [ ] **Still blocked**: actually re-running the two-level Monte Carlo needs a `simulate`
        callable driven by real per-cell daily climate forcing (VPD, soil water) — that needs
        TOPOHYDRO's `topoclimate_forcing` run against real gridded CHELSA-daily data using
        today's newly-fit lapse rate, not yet done (today's TOPOHYDRO fit produced the
        monthly coefficients, not a gridded forcing run). Re-fitting the monotone emulator and
        its held-out release gate both wait on that. TRY registration (§1.6) could still
        supplement psi_close_mpa/capacitance if it has compatible data — not checked yet.
- [ ] **MNEME**: build the real person-period panel from the vitality/disturbance exports (§2);
      fit the cloglog GLM + monotone GBM + NNLS stack against real dieback labels; run the actual
      nested blocked-CV validation designed in `configs/cv.yaml`.
- [x] **MERISTEM adult-niche model — built and fit per species, done (2026-09-30).**
      `antar.niche.adult` held only the Boyce index (a validation metric) before this — the
      presence-background model it validates didn't exist. Added `fit_presence_background`
      (penalised logistic regression) and fit it per species (not pooled — a first pooled
      attempt gave a weak, unstable result and was diagnosed as a real methodological error
      against the concept note's own "for each species" wording, not noise) against real GBIF
      presences, a real GBIF target-group background (not uniform-random), real CHELSA-
      BIOCLIM+ climate predictors, and real SoilGrids soil fractions. Results are honestly
      mixed: good signal for Fagus orientalis (Boyce 0.43) and both junipers (~0.34 each,
      consistent across folds); weak for the oaks (0.14, -0.07); Pinus kochiana (n=8)
      explicitly skipped as too sparse for a reliable fit rather than forced. See
      `configs/fitted/meristem_adult_niche.yaml` and IMPLEMENTATION_LOG.md for full detail,
      including a real CRS bug caught before trusting the soil extraction (soils.tif is in
      EPSG:32638, not WGS84).
  - [ ] **Still needed**: replace `CWD_approx = petmean - bio12` with TOPOHYDRO's real AET-
        based CWD once that water balance is run against real data. The attainable-height
        quantile model (real GEDI/canopy-height data, §2's structure export is registered but
        not yet loaded/processed) and growth-modifier calibration (needs tree-ring or
        plantation-survey growth-climate data, neither obtained) are both still open.
- [ ] **REFUGIUM**: run the full robust-refugium criteria against real MERISTEM/hazard outputs
      once both feed it real fits, not placeholders.
- [ ] **AEGIS**: run the CVaR portfolio optimisation and efficient frontier against a real
      scenario/model/parameter ensemble once §1.1's projection-source decision is made.
- [ ] **Pipeline wiring**: `workflow/Snakefile` is still 100% placeholder `echo` commands for
      every rule — replace each with the real `antar.*` entry point once its inputs are real,
      so the pipeline can actually run end-to-end rather than existing only as tested library code.
- [x] **Dev environment — installed (2026-09-30).** `pip install -e ".[dev]"` run against
      system Python (plus `scipy`/`scikit-learn`/`rasterio`/`netCDF4`/`isimip-client` installed
      separately as this session needed them for real pulls/fits). Full suite re-run: **143
      passed**, 0 failures — includes the new `fit_presence_background` test, everything else
      still green after today's real-data changes.
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
