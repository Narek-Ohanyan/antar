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

- [x] **RCP vs SSP scenario correspondence — resolved (2026-09-30), stale checkbox corrected.**
      This was unchecked but the actual work was already done later the same day it was written:
      `configs/scenarios.yaml`'s `rcp_ssp_correspondence` field writes down the real, checked
      numbers (RCP2.6 vs SSP1-2.6: ~421 vs ~446 ppm CO2 by 2100; RCP8.5 vs SSP5-8.5: ~936 vs
      ~1135 ppm, from the real RCMIP table §1.5 pulled), with explicit "never plot/average
      together under a shared low/high label without checking this first" guidance. The
      `future_ensembles.decision_note` in the same file also records the user's explicit choice
      not to mix them at all: CORDEX/CMIP5 and NASA-GDDP/ISIMIP3b/CMIP6 run as two separate,
      labelled ensembles, never merged into one "low/high" figure.
- [ ] **Resolution reconciliation — policy written (2026-09-30, `configs/resampling_policy.yaml`),
      two real fixes identified, neither applied yet.** Checked the actual code, not assumed:
      no GEE export anywhere calls `.resample()`, so every one -- continuous fields included --
      currently uses Earth Engine's nearest-neighbor default when reprojecting to the master
      grid. For categorical/boolean layers (land tenure, disturbance flags, WorldCover class)
      that's already correct, not a bug. For two real sources it's a genuine problem: SoilGrids
      (250m, real ~8x upsample to 30m) and ERA5-Land (~9km, real ~300x upsample) should use
      `.resample('bilinear')` -- every real gridded-TOPOHYDRO cell this session effectively reads
      off whichever coarse source pixel it happens to land nearest to, not a smooth field.
      CHELSA-daily/CHELSA-BIOCLIM+ are a documented non-issue: they're deliberately *not*
      resampled to 30m at all (point-sampled at their own ~1km grid), because TOPOHYDRO's real
      lapse-rate downscaling is the actual mechanism meant to carry that coarse reference down to
      true cell elevation -- a spatial resample would be redundant with, not a substitute for,
      that physical step. Real remaining work: apply `.resample('bilinear')` to `export_soils`/
      `export_era5land_forcing` in `gee_export.py` and re-run those two exports.
- [ ] **Unit audit per new variable — genuinely partial, not a blanket gap.** Real unit bugs were
      actually caught and fixed for every source pulled so far: CHELSA-daily's Kelvin encoding
      (verified against a real July/January Yerevan temperature, not assumed), GHCN-Daily's
      tenths-of-degree/mm raw encoding, SoilGrids' per-mille (not percent) clay/sand/soc encoding
      (caught mid-session, see IMPLEMENTATION_LOG.md), CO2's ppm units. What's still genuinely
      missing: a standing, systematic unit-bounds *test* per variable (the `ΔVPD/ΔT_max`-style
      physical-bound pattern already used for TOPOHYDRO/XYLEM's own derived quantities) — real
      bugs are being caught by hand each time a source is first used, not by a regression test
      that would catch a future regression automatically. TerraClimate specifically (not yet
      pulled) still needs its own check when it is.
- [x] **One forest/treeline definition, one reference period — satisfied for every source
      actually pulled so far, real check done per source, not assumed.** CHELSA-daily's pull
      window was deliberately extended to 1979 start specifically because the earlier 2000-2024
      attempt did NOT cover the stated `reference_period: [1991, 2020]` (caught and fixed, see
      IMPLEMENTATION_LOG.md). GHCN-Daily's lapse-rate fit explicitly filtered to "41 of 53
      stations with observations in the 1991-2020 reference period" before fitting. TerraClimate
      (not yet pulled) and any other future historical source still need the same check when
      they're actually added — this item isn't "done forever," it's "done for everything pulled
      so far," which is the honest scope of a rolling convention like this.

## 1. Data — real sources confirmed this session, not yet wired in

### 1.1 Future climate projections (highest-value finding)
- [x] **NASA/GDDP-CMIP6 on Earth Engine — UNVERIFIED flag resolved (2026-09-30), and the answer
      was negative.** `coll.aggregate_array('scenario').distinct().getInfo()` on the real ingested
      EE collection returned `['historical', 'ssp245', 'ssp585']` — ssp126 and ssp370 are
      genuinely absent (0 images) from this EE asset, despite NASA's own source archive having all
      four. This is exactly why NASA/GDDP-CMIP6 was demoted to a secondary/cross-check role in
      `configs/scenarios.yaml` rather than used as primary — not an oversight, a checked finding.
- [x] **ISIMIP3b direct access — UNVERIFIED flag resolved, real and working.** The "Anubis gate"
      concern was about the web portal, not the API `isimip-client` actually uses; a real pull
      (`scripts/pull_isimip3b.py`) ran successfully this session — 57 of 60 real GCM x SSP x
      variable combos completed (3 genuine upstream 900s job-timeouts, not a pipeline failure).
      `isimip-client` + `cutout_bbox` is confirmed the right, working mechanism.
- [x] **Decision resolved (2026-09-30) — stale checkbox corrected.** The user chose option (b):
      run CORDEX/CMIP5 and ISIMIP3b/CMIP6 as two separate, labelled ensembles rather than picking
      one (`configs/scenarios.yaml`'s `future_ensembles.decision_note`) — the most scientifically
      thorough option, at roughly double the future-projection workload, accepted deliberately.
      NASA/GDDP-CMIP6 ended up demoted to a secondary cross-check role (see above), not one of the
      two primary ensembles.

### 1.2 Historical cross-checks and validation data
- [x] **TerraClimate — real export done and verified (2026-10-01).**
      `antar.io.gee_export.export_terraclimate` (new): real 1991-2020 (the project's own stated
      `reference_period`) monthly climatology, 9 real bands x 12 months = 108 bands (pr, tmmn,
      tmmx, pet, aet, soil, ro, pdsi, vpd), each band's documented GEE scale factor applied
      explicitly (checked directly against the real catalog page, not assumed -- Earth Engine
      does not auto-apply a catalog-documented scale). Exported at 4000m with
      `.resample('bilinear')`, applying this session's own resampling policy
      (`configs/resampling_policy.yaml`) immediately rather than creating a third entry needing
      the same future fix soils/ERA5-Land still need. Real values checked at the master grid's
      centre pixel before trusting it: Jan precip 25.4mm, Jul 46.7mm, Jan tmmn -9.4degC, Jul
      tmmx 25.0degC, Jan PDSI -0.09, Jul VPD 1.02 kPa -- all physically sensible for Armenia.
      Registered in `configs/manifests/terraclimate.yaml`.
  - [ ] **Still open**: not yet actually compared against TOPOHYDRO's own real soil-bucket
        output or CHELSA-daily's climatology -- pulled and verified as a real, standalone
        dataset, cross-check comparison itself is separate follow-up work.
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
- [ ] **`data.opendata.am` (Data Catalog Armenia) — real datasets confirmed, actual pull blocked
      by real, external, currently-down infrastructure (2026-09-30), not a research or engineering
      gap.** Queried the live CKAN API directly, found all three real datasets with direct URLs,
      then actually tried to download every one. All 7 ANAU soil layers (clay/silt/sand/bulk
      density/pH/CEC/EC/SOC, real lab-measured values) plus the forest-cover layer live on
      `armsis.cas.am`, which returns **NXDOMAIN** — confirmed with a direct `nslookup`, genuinely
      down/renamed, not slow or rate-limited. The 3 forest-degradation datasets (Global Forest
      Watch's ArcGIS Hub) hit a real server-side 500 on the CSV export and `{"error":{"code":400,
      "message":"Invalid URL"}}` on the underlying REST FeatureServer's own bare metadata request
      — GFW's own hosting problem, tried multiple encodings, not a request-formatting issue here.
      Full detail and real next steps in `configs/manifests/opendata_am_blocked.yaml`. Retry later
      by re-checking DNS first; no further effort here fixes someone else's outage.
- [x] **EPPO pest/pathogen records — done (2026-09-30).** Real CSV (107 Armenia organism
      records) downloaded and saved (`data/eppo/eppo_armenia_organisms.csv`,
      `configs/manifests/eppo.yaml`): confirmed real, tree-relevant threats present, including
      *Ips sexdentatus* (pine bark beetle) "present, widespread" and *Bursaphelenchus
      xylophilus* (pine wood nematode) "present, restricted distribution" — directly relevant to
      the pine functional group. One real limitation found, not assumed: this is national-level
      presence/absence only, no coordinates — usable to flag qualitative biotic risk per species,
      not as a per-cell spatial covariate the way vitality/disturbance/climate data are. Not yet
      wired into MNEME's biotic hazard cause or saved to a manifest.
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
- [x] **XFT — done, real data in hand.** See IMPLEMENTATION_LOG.md 2026-09-30.
- [ ] **TRY — real data request submitted (2026-09-30), waiting on TRY's own 14-day custodian
      window, not on any further action from this project.** Request **52804** (PI: Narek
      Ohanyan), 13 real trait IDs (6, 189, 709, 711, 3105, 3468, 4334, 4335, 4494, 4495, 4502,
      4503, 4504 — turgor-loss-point water/osmotic potential, leaf water capacitance, leaf
      cuticular conductance, bud burst phenology, and the full P12-P50-P80-P88 vulnerability
      curve + hydraulic safety margins), 12 real species IDs covering all four functional groups
      -- three of them *exact* Armenian target species TRY has directly that XFT didn't (Fagus
      orientalis 23903, Quercus macranthera 45423, Juniperus excelsa 31560, plus Juniperus
      polycarpos only as subsp. of excelsa, 496095) -- see IMPLEMENTATION_LOG.md for the full
      species-ID provenance (which are exact matches, which are the same congeneric proxies XFT
      used, which are confirmed genuinely absent from TRY: Q. iberica/frainetto/humilis,
      P. kochiana). Dataset custodians have 14 days from submission to adjust permissions; after
      that, the PI (the user) must actively go to `try-db.org/TryWeb/Prop01.php` ("Requests -
      Review Custodian Responses And Start Data Release") to trigger the actual release -- this
      does not happen automatically. **Action needed in ~14 days: the user visits that URL and
      starts the release; only then is there a real file to pull and wire in.**

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
- [ ] **CORDEX future extraction on hyperion — explicitly deferred to v2.1 (2026-09-30, user's
      decision), not part of the current v2.0 scope.** Only 2 of 4 variables (`pr`, `tas`) done
      for 1 of 8 GCM×RCP combinations (MOHC-HadGEM2-ES/rcp26); `tasrange`/`tasskew` for that
      combination, plus all 7 remaining combinations, still pending — blocked on hyperion network
      access, which isn't a research or engineering task that resolves with more effort. ISIMIP3b
      (real, working, already 57+/60 combos pulled this session) is the real, unblocked primary
      ensemble for v2.0; CORDEX stays a real, labelled second ensemble member (per §1.1's dual-
      ensemble decision) once hyperion access is available, not before.
- [ ] Future `tasmax`/`tasmin` derivation from CORDEX `tas`+`tasrange` — same v2.1 deferral,
      follows directly from CORDEX itself being deferred.
- [ ] **TRY trait data — explicitly deferred to v2.1 (2026-09-30, user's decision), not part of
      the current v2.0 scope.** Real request 52804 already submitted (see above) -- this is a
      genuinely time-gated external dependency (TRY's own 14-day custodian window, then a manual
      PI release step), not something more work closes sooner. v2.0 proceeds on XFT's real data
      and today's documented placeholders (`psi_close_mpa`, `capacitance_mmol_m2_mpa`,
      `individual_sd`); TRY's real data folds in as a v2.1 upgrade once it actually arrives.

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
  - [x] **Terrain concavity, wind-height correction, ERA5-derived net radiation — built, done
        (2026-09-30).** `src/antar/climate/terrain.py` (TPI-style cold-air-pooling concavity,
        vectorized, edge-correct), `vapour.wind_speed_2m` (FAO-56 10m->2m), `radiation.
        net_radiation_from_era5` (direct from ERA5-Land's real shortwave/longwave rather than
        the FAO-56 parametric approximation). 7 new tests, full suite 152 passed.
  - [x] **The first real, gridded TOPOHYDRO forcing run — done (2026-09-30).**
        `scripts/run_topohydro_grid.py`: 80 real grid points, year 2019 (CHELSA `pr`'s real
        coverage limit), every input real and streamed (terrain/soils/ERA5-Land via `/vsicurl/`
        range requests, nothing downloaded whole). 78/80 points produced real output; real
        elevation-temperature correlation -0.97; CWD/WSI/psi_soil/GDD/growing-season all in
        physically plausible ranges across Armenia's real 135-3630m elevation span. Saved to
        `configs/fitted/topohydro_grid_run_2019.yaml`. Three explicit documented placeholders
        remain (`calm_clear_night_frac`, `gdd_budburst`, `rooting_depth_mm`) — no real data
        source exists for any of them yet.
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
  - [x] **Mechanistic hazard — real two-level Monte Carlo, done (2026-09-30).**
        `scripts/fit_xylem_mechanistic_hazard.py`: the real Sec. 6.3 50x200 two-level Monte
        Carlo, run against the gridded TOPOHYDRO forcing above, for each of the 3 functional
        groups with a real hyper_sd. `juniper_arid_conifer` (n=1) skipped, not forced, same
        precedent as MERISTEM's `Pinus kochiana`. A real gap handled explicitly: no real
        individual-level (within-population) trait-variance estimate exists anywhere in the
        pulled data, so `individual_sd` is a documented placeholder fraction (0.5) of the real
        hyper_sd, not a measured quantity. Saved to
        `configs/fitted/xylem_mechanistic_hazard_2019.yaml`.
  - [ ] **Still open**: the monotone emulator re-fit and its held-out release gate (needs a larger
        cell sample than 78 points to be meaningful). TRY registration (§1.6) could still
        supplement psi_close_mpa/capacitance if it has compatible data — not checked yet.
- [x] **MNEME — first real person-period panel built, done (2026-10-01); the stacked hazard
      model itself did NOT fit, a real and informative null result, not a bug.**
      `scripts/fit_mneme_hazard_panel.py` ran to completion across the full real 2010-2019
      window: real dieback labels from the real kNDVI vitality composites + Hansen/MODIS
      disturbance ancillary (already exported to Drive, registered, but never loaded into
      anything until this session), crossed with real climate covariates from the gridded-
      TOPOHYDRO machinery extended across all 10 real years. Real result: **778 real person-year
      rows across 78 places, 0 real dieback-onset events.** At this sample size (78 points, one
      real 10-year window), the strict standardised-anomaly dieback rule (Sec. 7.1) genuinely
      never fired -- `fit_stacked_hazard` was correctly never called
      (`status: insufficient_real_events_or_blocks_for_a_meaningful_fit` in the real saved
      output); forcing a fit on an all-zero label would have produced a meaningless, degenerate
      model, not a real one. A real, stated scope reduction on the design itself, independent of
      the event-rate finding: no real stand-age/DLNM/stand-structure data exists, so the design
      that *would* have been fit was climate-only Mundlak (bypassing `build_hazard_design`'s
      mandatory age spline, which needs real data to build a valid B-spline basis, not a
      fabricated constant). Saved to `configs/fitted/mneme_hazard_panel_2010_2019.yaml`.
  - [ ] **Still open**: a larger real sample (more than 78 points, and/or a longer real window
        past 2019 once real precipitation exists past CHELSA's 2019 gap) is the real, concrete
        next step suggested by tonight's zero-event finding -- not a code fix. Beyond that: the
        full Eq. 7.3 design (age spline, DLNM drought-legacy cross-basis, stand terms) once real
        stand-age/structure/tree-ring data is obtained; Se/Sp estimation once a stratified
        interpreted-pixel-year sample exists; the full nested (not single-split)
        blocked-CV validation designed in `configs/cv.yaml`, which needs more than 80 points/10
        years to be statistically meaningful.
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
- [ ] **CO2 physiological effects — both open questions in this item resolved (2026-09-30),
      real remaining work now small and well-scoped.** (1) §1.5's CO2 pathway table was already
      in hand before this item was even written: real RCMIP v5.1.0 concentrations (Meinshausen
      et al. 2020), pulled, filtered to this project's exact scenarios, saved at
      `configs/co2_concentration_pathways.csv` (`configs/manifests/co2_pathways.yaml`). (2) The
      concept note's own stated limitations (`docs/concept_v2/sec10_plan.tex`) answer the
      "mechanistic or not" question directly, and the answer is genuinely useful: "Carbon
      starvation, root-shoot dynamics and the CO2 effect on water-use efficiency **enter only
      through scenario knobs**" — i.e. the concept note itself does NOT call for a mechanistic
      CO2-to-WUE formula inside MERISTEM/XYLEM's growth or hydraulic equations; it explicitly
      scopes this as an external, scenario-level adjustment (e.g. a WUE or CWD-sensitivity
      multiplier that varies by CO2 pathway/scenario), not a term to derive and insert into an
      existing formula. So "wire it into whatever growth/WUE term is meant to use it" was
      premised on a formula that, per the concept note's own design, doesn't exist and isn't
      supposed to. Real remaining work: implement the scenario-knob mechanism itself (reads the
      real CO2 table already in hand, applies a documented multiplier per scenario to whichever
      downstream quantity is chosen) — small, unblocked, no new data or research needed.

## 4. Validation

- [ ] Design-based accuracy assessment against a stratified validation sample (concept note's own
      stated minimum: ≥300 plots) — blocked entirely on real ground-truth data (§1.6).
  Structural gap or not: the pipeline's `validation_report` rule (AOA maps, calibration,
  model card) is also still an `echo` placeholder in `workflow/Snakefile` — same §3 item.

## 5. Interface — deferred until real fits exist, scoped for when that decision is made

Not started by design (explicitly deferred pending real model fits — see
`IMPLEMENTATION_LOG.md`). Scope as you described it, for when the go-ahead is given:

- [ ] Home page.
- [ ] **Comprehensive output, explicitly required by the user (2026-09-30): the UI must surface
      every real quantity the framework is able to compute, not a curated subset.** Concretely,
      by the time the UI is built this should include (whichever of these have a real fit by
      then): TOPOHYDRO's per-cell downscaled climate (T/P/PET ensemble/CWD/WSI/snow/GDD/growing
      season), XYLEM's mechanistic hazard (per functional group, with its real outer-loop
      uncertainty spread), MNEME's statistical hazard, MERISTEM's adult-niche suitability *and*
      attainable-height/treeline outputs, REFUGIUM's viability/robust-refugium/risk-averse score,
      and AEGIS's portfolio/efficient-frontier results -- each with its real status (fitted /
      reduced-scope / placeholder-flagged) shown, not silently omitted if incomplete.
- [ ] Interactive map with buttons/dropdowns to select what's being predicted (refugia score,
      viability, hazard, attainable height, etc.), scenario/GCM/RCP-or-SSP member, and horizon
      (2050/2080/2100).
- [x] **Treeline diagnostic — real single-year run done (2026-09-30), explicitly required by
      the user, not optional scope.** `scripts/compute_treeline_diagnostic.py`: real potential
      treeline elevation computed for all 78 real 2019 grid cells, using each cell's own real
      downscaled growing-season mean temperature and elevation (mathematically exact, not an
      approximation -- the lapse relation is linear, so inverting from the cell's own real
      (z_cell_m, growing_season_mean_t_c) gives the identical z_tl as inverting from the true
      upstream reference). Real, physically plausible result: potential treeline elevation
      ranges 1,586-3,616 m across the 78 cells (matches Armenia's known real treeline zone,
      ~2,000-2,600 m); only 1/78 cells currently sits above its own real climatic ceiling. Uses
      the real fitted Apr-Sep growing-season lapse rate (-6.89 K/km), not a flat annual average.
      Saved to `configs/fitted/treeline_diagnostic_2019.yaml`.
  - [ ] **Still open**: this is one real year, one ensemble member (2019 only) -- real treeline
        *change* (the user's actual ask) needs this same diagnostic re-run across the ISIMIP3b
        scenario ensemble and the 2050/2080/2100 horizons once that full propagation exists (the
        "full production" scope flagged elsewhere in this roadmap), plus an explicit map layer/
        panel in the interactive map, not folded silently into "attainable height."
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
- [ ] **README.md rewrite — deliberately held (2026-10-01, user's explicit decision), not
      forgotten.** The current README (title, opening paragraph, `species_traits.csv`
      description) is already partly stale as of tonight's real work -- "data access... declared
      but not implemented" and "nothing produced from them is a result" are both now false (real
      TOPOHYDRO/XYLEM/MNEME/treeline/TerraClimate results exist). User chose not to patch this
      incrementally: hold it until REFUGIUM, AEGIS and the UI all land, then rewrite the whole
      top-level description in one pass rather than multiple partial edits. Do this rewrite when
      those three are done, not before, even though the text is known-inaccurate in the meantime.
