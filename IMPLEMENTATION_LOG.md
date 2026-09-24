# Implementation log

Chronological record of assumptions and decisions not fully specified by the concept note.
Each entry states the assumption and the rationale; nothing here is a result.

## 2026-09-24 — package rename `foracca2` -> `antar`

Package directory, `pyproject.toml` name/description, all imports (`src/`, `tests/`, `scripts/`,
`workflow/Snakefile`), and `environment.yml` renamed from `foracca2` to `antar`. Full test suite
(77 tests) verified green under the new name before any further change.

Engine naming: the five scientific engines use the names TOPOHYDRO (climate/), XYLEM (hydraulics/),
MNEME (hazard/), MERISTEM (niche/), REFUGIUM (viability/), AEGIS (decision/) in README, docstrings,
and other public-facing text, per the naming brief. Internal package/folder names are unchanged
for code clarity and to avoid a larger diff than the interfaces require. `validation/`,
`uncertainty/` and `io/` are cross-cutting support packages, not named engines, and keep plain names.

**Prior-project name.** The concept note's own LaTeX source under `docs/concept_v2/` (title,
running headers, `main.tex` etc.) still names the prior project throughout, since that document
*is* the prior source material being read for its content. It has been left unmodified rather than
edited or removed: it is archival background, not documentation this project authors going forward.
Every other file in the repository (code, docstrings, comments, config, `README.md`, package
metadata) has had the prior project's name removed or reworded without inventing new attribution
text — e.g. "the v1 FORACCA training table" -> "the prior-project v1 training table". A single
placeholder was left in `README.md`'s Acknowledgments section for manual completion. This scoping
of "documentation" (excludes the archived prior document, includes everything authored here) is an
assumption, not stated explicitly in the naming brief — flagging it for confirmation.

**No git repository exists yet** in this directory (`foracca2/`) or its parent. Author identity in
the global git config is already set correctly (not a placeholder). Repository initialisation was
not performed and is left for explicit instruction, since it is a one-way decision about where
history starts.

## 2026-09-24 — TOPOHYDRO implementation (concept note Sec. 5)

New, tested code: `climate/downscale.py` (monthly lapse-rate and precipitation-gradient
regression, cold-air-pooling index, temperature/precipitation downscaling), two additions
to `climate/indices.py` (`growing_season_length`, `growing_season_mean_temperature`, Eq. 8.3's
own definition -- days with T_mean >= T0, not a contiguous-run rule), two additions to
`climate/pet.py` (`stomatal_vpd_response`, `canopy_conductance`, Eq. 5.5's `g_c` term),
`climate/stand_coupling.py` + `configs/stand_defaults.yaml` (placeholder LAI), `io/manifest.py`
(dataset manifest schema/reader/writer/checksum, no download logic), and `climate/forcing.py`
(`topoclimate_forcing`, the per-cell daily orchestration). Full test suite green throughout.

Assumptions made where the concept note states the physics' *role* but not a closed form,
each chosen for the least additional invention and flagged rather than silently resolved:

* **Cold-air pooling.** No functional form is given beyond "delta_CAP <= 0, scales a terrain
  concavity index by the share of clear, calm nights." Implemented as
  `delta_CAP = -k_cap * max(concavity, 0) * calm_clear_frac`, the simplest form with that
  qualitative behaviour (zero on convex/ridge terrain, zero with no calm/clear nights,
  monotonic in both factors), `k_cap` a fittable coefficient with a neutral default of 1.0.
  Confirmed from the text: delta_CAP corrects T_min only (nocturnal drainage), not T_mean/T_max
  -- `forcing.topoclimate_forcing` applies it that way.
* **Stomatal VPD response, `f(VPD)` in Eq. 5.5.** No functional form or coefficient is given.
  Implemented as the standard Jarvis-type exponential `f(VPD) = exp(-k_vpd * VPD)`; `k_vpd`
  (kPa-1) has **no default** and must be supplied, since no citation-backed number for it was
  available to put in confidently.
* **Canopy aerodynamic conductance `g_a` for `canopy_pm`.** The concept note never specifies a
  canopy-scale `g_a` formula (only the reference-crop FAO-56 form `u2/208` is precedented in
  this codebase, for a short reference crop, not a forest canopy). Rather than reuse that
  formula outside its validated scope, `forcing.topoclimate_forcing` treats `g_a` as an optional
  caller-supplied input: `canopy_pm` is computed only if both `g_a` and `g_c` are given, and is
  silently absent from the PET ensemble otherwise. This is a real gap, not resolved here.
* **Slope-corrected net radiation.** `antar.climate.radiation` documents itself as geometry-only
  (a shortwave beam-direction ratio; no cast shadows, no beam/diffuse split). `forcing.py` takes
  net radiation `rn_mj_m2` as a direct upstream input (per the module's own stated boundary) and
  reports the shortwave slope ratio only as a separate diagnostic (`rs_slope_mj_m2`) -- it is
  *not* multiplied into `rn_mj_m2`, because net radiation carries a longwave term the beam-
  geometry ratio does not apply to. The full Eq. 5.4 beam/diffuse/terrain-shading/reflected
  decomposition is not implemented.
* **PET stays an ensemble through the soil bucket.** Rather than pick one PET formulation to
  drive the water balance, `forcing.topoclimate_forcing` runs the soil bucket, CWD and WSI once
  per PET formulation, keyed by formulation name. This follows directly from the concept note's
  own framing of PET as a structural uncertainty to be carried, not collapsed.
* **`gdd_budburst`** (the GDD threshold used by `late_frost_days`) is species/functional-group
  specific and is not given a numeric value anywhere in the concept note; `forcing.py` requires
  it as a mandatory argument rather than defaulting it.
* **Validation harness for Sec. 5.5** (leave-one-station-out CV of downscaled T/P against
  ERA5-Land/TerraClimate/soil moisture) needs no new code: `antar.validation.splits.leave_group_out`
  already implements leave-one-group-out generically, and station ID is a valid group. Running it
  against real station data waits on the hyperion data pull.
* Bias adjustment (QDM, `antar.climate.bias`) is treated as an upstream step, applied to the
  reference forcing before it reaches `topoclimate_forcing`, not called from inside it -- the two
  operate at different temporal granularities (distributional-by-month vs. per-day-per-cell).

## 2026-09-24 — outer directory renamed `foracca2/` -> `antar/`

The repository root (previously `foracca2/`, sibling of the Python package renamed earlier)
is now `antar/`, at the same location under `FORACCA_v2/`. Full test suite re-verified green
from the new path.

## 2026-09-24 — XYLEM implementation (concept note Sec. 6)

The mechanistic kernels (vulnerability curve, biphasic g_min(T), two-phase HFI engine,
one-level Monte Carlo, generic monotone-emulator tools) already existed and were tested.
New, tested code completes what Sec. 6 specifies beyond those kernels:

* **Two-level Monte Carlo** (`hydraulics/monte_carlo.py`: `sample_trait_hyperparameters`,
  `two_level_failure_probability`), closing Eq. 6.5's stated structure exactly: an outer
  loop of hyperparameter draws (knowledge uncertainty) wrapping the existing inner loop of
  individual draws.
* **Stress-spell summaries** (`hydraulics/spells.py`: `longest_true_run`,
  `stress_spell_summary`), the four named emulator features from the "Emulation" paragraph
  ("longest closed spell, mean VPD24 and T during it, number of closed days"). "During it"
  was read as referring to the longest spell specifically, not the whole record.
* **Emulator training + the note's own release gate** (`hydraulics/emulator.py`:
  `build_hazard_emulator`, `emulator_max_abs_error`), a generic LHS-over-named-axes -> run
  engine -> fit monotone GBM -> held-out max-abs-error wrapper. The concept note's 0.02
  absolute-hazard release gate is implemented as `emulator_max_abs_error` itself; the *tests*
  assert a looser bound because they deliberately use a fast-test-sized LHS (hundreds to a
  few thousand samples), not a production budget -- hitting 0.02 in practice is a sample-size/
  capacity tuning question for the real model-fitting phase, not something a unit test should
  claim to have proven.
* **Link-scale recalibration** (`hydraulics/calibration.py`: `recalibrate_hazard`), Sec. 6.4's
  `cloglog(h~) = a_g + b_g cloglog(h_mech)` transform. Duplicates (rather than imports) the
  small cloglog/cloglog-inverse pair already in `hazard/models.py`, matching this codebase's
  existing convention of each engine package staying import-independent of its siblings.
* **Cross-engine wiring** (`hydraulics/pipeline.py`: `mechanistic_hazard_for_cell`): XYLEM's
  actual entry point from TOPOHYDRO output, per Fig. 2. This surfaced a real gap on the
  TOPOHYDRO side -- `climate/forcing.py`'s `CellTopoclimate` did not yet expose soil water
  potential (Phase 1's closure trigger) or atmospheric pressure, both needed downstream.
  Extended `topoclimate_forcing` to compute `psi_soil_mpa` (Eq. 5.8, Clapp-Hornberger, once
  per PET formulation like the rest of the water balance) and to return `pressure_kpa`; this
  is a required-argument change (five new Clapp-Hornberger soil parameters), updated
  everywhere `topoclimate_forcing` is called. `pipeline.py` uses T_max as the leaf-temperature
  proxy and VPD24 for Phase 2, per Sec. 6.4's own documented simplification.

**Deferred, not built:**
* The ABC pattern-oriented calibration step (Sec. 6.4) needs real tree-ring/dieback data
  (hyperion) for two of its three acceptance criteria and MERISTEM's Boyce index for the
  third -- cross-engine and data blocked, not implementable in isolation. `calibration.py`
  documents this boundary; only the (data-agnostic) recalibration transform is implemented.
* Separate trait classes for seedlings/saplings/mature trees (Sec. 6.3) need no new code:
  `Traits` is already a plain dataclass, so each ontogeny stage is just a separate `Traits`
  instance with its own capacitance/rooting-depth/psi_close values.

## 2026-09-24 — MNEME implementation (concept note Sec. 7)

The learners and gate already existed and were tested: `CloglogGLM`, `monotone_hgb`,
`stack_nnls` (Sec. 7.3's three learners plus stacking), `hazard.gating` (Eq. 7.4's blend) and
`validation.aoa.AreaOfApplicability` (Sec. 7.3's DI + upper-whisker threshold, computed
against other CV folds, exactly as specified) -- these needed no new code, just confirming
they already cover Sec. 7.3 in full. `panel.mundlak_decompose` and `panel.person_period`
(within/between split, censored panel expansion) also already existed. New, tested code:

* **Misclassification correction** (`hazard/observation.py`: `correct_for_misclassification`),
  Eq. 7.1 exactly, given Se/Sp. Estimating Se/Sp themselves needs a stratified sample of
  interpreted pixel-years (real data, not available) -- only the correction formula is built.
* **Dieback event rule** (`hazard/observation.py`: `dieback_event`): the standardised-anomaly
  threshold rule from Sec. 7.1, with the >=2-growing-season persistence requirement from
  Table 3's caveat column. The note gives the thresholds (-2, -1, 10-year window, 2-season
  recovery) numerically but not the spread estimator for "standardised anomaly relative to
  its trailing median" -- implemented as trailing-window standard deviation, the direct
  reading of that phrase. Years too near either end of the record to evaluate the trailing
  or recovery window are left unflagged (undetermined), not asserted either way. The upstream
  *segmentation* that would produce a cleaned kNDVI series (LandTrendr/CCDC) is a substantial,
  data-dependent algorithm in its own right and is not implemented -- this function takes an
  already-segmented (or raw) annual-max-kNDVI series and a caller-supplied disturbance flag.
* **Rare-event weighting** (`hazard/panel.py`: `rare_event_weights`), Sec. 7.2's stated recipe
  exactly: weight 1 for kept events, 1/pi0 for non-events sampled at known rate pi0 (King &
  Zeng 2001).
* **Space-for-time bracket** (`hazard/panel.py`: `no_adapt_bracket`, `full_acclimation_bracket`,
  `space_for_time_bracket`), Eq. 7.5's two projection brackets.
* **Age/time-since-disturbance spline** (`hazard/dlnm.py`: `spline_basis`), exposing the
  existing private B-spline machinery as f(a_i,t)'s basis (Eq. 7.3), reused rather than
  duplicated from `crossbasis`'s exposure axis.

**Scoped out, not built:** a "fit the full Eq. 7.3 stacked model" orchestrator. Unlike
TOPOHYDRO's forcing pipeline or XYLEM's Monte Carlo (deterministic physics, genuinely
testable against synthetic weather), Sec. 7.3's learners are *fitted statistical models* --
exercising them meaningfully needs the real dieback-labelled panel (blocked on the hyperion
pull), not synthetic labels that would validate nothing real. All the pieces it would call
(`CloglogGLM`, `monotone_hgb`, `stack_nnls`, `dlnm.crossbasis`, `dlnm.spline_basis`,
`panel.mundlak_decompose`, `validation.splits.block_kfold`, `validation.aoa`) are already
built and tested individually; assembling Eq. 7.3's design matrix from them, once real data
exists, is direct column concatenation, not a hidden decision -- also true of the existing
"elastic-net cloglog GLM" naming in Sec. 7.3: `CloglogGLM` (pre-existing, not touched here) is
ridge-penalised (L2) only, not true elastic net (L1+L2); it keeps the stated
interpretable/linearly-extrapolating property but not exact sparsity. Flagged, not rebuilt.

## 2026-09-24 — Earth Engine access and the real vitality-composite / disturbance pull

`antar.io.gee_export` implemented for real (was a stub): kNDVI (Camps-Valls et al. 2021,
kNDVI = tanh(NDVI^2), exact under their recommended adaptive sigma, not an approximation)
from HLS v2.0 (`NASA/HLS/HLSL30`/`HLSS30` v002) merged with Landsat Collection 2 SR
(L5/7/8/9), each with its own QA-band cloud/shadow/snow masking; yearly growing-season
(Jul-Aug) median composites 2000-2024; LandTrendr (Kennedy et al. 2010) segmentation with
the published LT-GEE default parameters (not tuned); and a disturbance-ancillary layer
(Hansen GFC `lossyear` OR any MODIS MCD64A1 burn that year) feeding
`hazard.observation.dieback_event`'s `no_disturbance` input directly, computed
independently of the vitality signal so dieback is never confused with harvest/fire.

Google Earth Engine project: `pure-highlander-495708-a9` (registered by the user).
Google Cloud Storage was the first export destination tried but the project has no
billing account linked (creating one needs a payment method -- the user's decision, not
made here); switched to Google Drive export instead (`ee.batch.Export.image.toDrive`),
confirmed to have ample quota (~1.7 TB free). A small (2-year, dry-run) submission was
verified end-to-end and cancelled before committing to the real one.

**Real 2000-2024, full-Armenia, 30 m export submitted** (Drive folder `antar_gee_exports`):

| Task | Earth Engine task ID | Product | Status |
|---|---|---|---|
| Vitality composites | `YYAQVTPB5MTSMQVSJKIX7PRO` | `kndvi_{year}`, `valid_count_{year}`, 2000-2024 | running |
| LandTrendr segmentation (v1, failed) | `UI7IB32PW6L5KF7CVKNDXWRA` | -- | **FAILED**, see below |
| LandTrendr segmentation (v2, fixed) | `4ZELRUNLTE7CY25JBEMG57DN` | `{year,original,fitted,vertex}` per year, flattened | running |
| Disturbance ancillary | `XLH4TJOJYOVCT4P7SNOE672A` | `no_disturbance_{year}`, 2000-2024 | **completed** |
| Structure (GEDI + CHM) | `Y2WR5SXAA2DXPCEVU7ZVI2KL` | `gedi_rh98_median`, `gedi_valid_count`, `eth_chm_2020` | running |

These are asynchronous server-side batch jobs; submission does not mean completion. Once
they finish and the GeoTIFFs are downloaded into the local, gitignored `data/` directory,
each becomes an `antar.io.manifest.ManifestEntry` (source, version, checksum) rather than
a hardcoded path, per the manifest convention -- not yet done, pending completion.

**hyperion connectivity, ongoing**: the CHELSA remote-extraction pull (see the entry
below on the ncks approach) hit repeated network interruptions -- some transient (Google
Drive uploads recovered via retry), one not: hyperion.wsl.ch itself became unreachable
mid-transfer (confirmed via ping: general internet fine, hyperion specifically timing
out), killing the `pr` variable pull at 87.6 of ~108 MB. The incomplete file was deleted
rather than left looking valid. This is on the institution's network, not something
fixable from here; retry once hyperion is reachable again.

**Bug found and fixed: LandTrendr export failed** with `Input array has length 24 on axis
1, but 25 labels provided`. Root cause: LandTrendr's per-pixel output array drops years
where that pixel's composite was masked (no valid satellite observation in that Jul-Aug
window), so the array's length along the "years" axis varies across the raster --
confirmed directly (a single test point had length 25, but `arrayFlatten` with a
fixed 25-label list failed region-wide, meaning some pixels elsewhere in Armenia had
fewer). Fixed by padding every pixel to a constant length with `arrayPad([4, len(years)],
-9999)` before `arrayFlatten`; -9999 marks "this pixel had fewer valid years than the
record," not a real observation -- verified region-wide (not just at one point) before
resubmitting. `export_vitality_composites` in `io/gee_export.py` updated; the other three
export functions build plain band-concatenated images, not per-pixel arrays, so they were
never affected by this failure mode.

**Real structure export** (`export_structure`, GEDI L2A RH98 2019-2023 + Lang et al. 2023
ETH global canopy height, `users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1`) also
implemented and submitted. WSL S2-VHM, named in the module's original stub docstring, has
no confirmed public Earth Engine asset -- checked directly under several plausible ids, not
found.

**WSL S2-VHM found**, not on Earth Engine but on EnviDat (WSL's own data repository), user-
supplied link: Rueetschi & Ginzler (2025), "Sentinel-2 Vegetation Height Model Armenia,"
EnviDat, doi:10.16904/envidat.690, CC-BY-SA. A CNN trained on Lidar-NFI-referenced VHM in
Switzerland, spatially transferred to Armenia; annual maximum vegetation height at 10 m
from Sentinel-2 (May-Sept composites), 2017-2025, 9 GeoTIFFs (~28 GB total). This is an
Armenia-specific, purpose-built product -- more directly relevant to MERISTEM's attainable-
height target than the global GEDI/ETH-CHM pair, which stay as a complementary/cross-check
source, not replaced. Downloaded directly (not via Earth Engine); exact resource URLs obtained via EnviDat's
CKAN API (`package_show`).

**Local disk correction**: the first attempt downloaded straight into `data/structure/
s2_vhm/` for all 9 years sequentially, but local disk had only ~13-16 GB free against a
~28 GB series -- caught before it filled (the user flagged it; confirmed with `df -h`, only
~2.5 GB of the first file had landed). Killed that download, deleted the partial file, and
switched to download-one / upload-to-Drive / delete-local-copy per year, so peak local
usage stays ~1 file (~3.2 GB) instead of the full series. Each file's SHA-256 is computed
locally before upload and recorded (script:
`/private/tmp/.../scratchpad/pull_s2_vhm_to_drive.py`, run outside the repo). Manifest
entries (source, version/DOI, checksum, Drive file id in place of a local path) to be added
once the run completes. Note: this dataset's own tags include the *prior* project's name
(EnviDat's page for this dataset, "within the framework of the FORACCA project" — that is
the dataset's real, correct citation context and is not this repository's branding; recorded
here as citation only, consistent with the acknowledgments-only placement decided at the
package rename).

## 2026-09-24 — closing the two flagged MNEME gaps

**Elastic-net `CloglogGLM`.** Added an `l1` parameter alongside the existing `l2`
(default `l1=0.0`, exactly reproducing the prior ridge-only behaviour via the same
closed-form solve). `l1 > 0` switches the per-IRLS-iteration weighted-least-squares
subproblem to cyclic coordinate descent with soft-thresholding (Friedman, Hastie &
Tibshirani 2010) -- the standard IRLS-outer/coordinate-descent-inner algorithm for
penalised GLMs, not a new method invented here. Verified against the closed-form ridge
solve at `l1=0` (exact match) and against the KKT stationarity conditions directly at
`l1>0` (implementation-independent correctness check, not a comparison against another
library's possibly-differently-scaled elastic-net convention).

**The Sec. 7.3 stacked-model orchestrator** (`hazard/pipeline.py`): `build_hazard_design`
(column-concatenates the age spline, Mundlak within/between climate columns, DLNM
cross-basis and stand terms -- no new decision beyond what Eq. 7.3 itself specifies),
`fit_stacked_hazard` (nested blocked-CV: GLM + monotone GBM out-of-fold link predictions
stacked by NNLS, RF fit alongside for benchmark comparison only, never in the stack, per
Sec. 7.3's own "cannot extrapolate" reasoning), and `gate_hazard` (thin composition of the
already-existing `validation.aoa.AreaOfApplicability` + `hazard.gating` into Eq. 7.4).

One judgment call: `stack_nnls`'s own docstring calls for "a smoothed empirical link
target," not the per-observation cloglog of a raw 0/1 label (undefined at y in {0,1}
anyway). Implemented as the cloglog of the empirical event rate within quantile bins of
the two learners' averaged OOF link prediction -- the same binning idiom
`validation.metrics.murphy_decomposition` already uses, not a new smoothing technique.

All of this is tested against synthetic panels (signal recovery, KKT conditions, stack
weights summing to 1 and non-negative, AOA gate favouring h_stat inside / h_mech outside).
None of it has been run against the real dieback panel yet -- that now depends on the
Earth Engine pull below finishing and being turned into person-period rows with real
events, which is the actual remaining blocker, not the code.

## 2026-09-24 — MERISTEM implementation (concept note Sec. 8.1)

Adult niche (Boyce index), attainable-height quantile regression, Chapman-Richards growth,
the water modifier f_W, and establishment hazard already existed and were tested. Two
cross-cutting pieces MERISTEM needs also already existed elsewhere and needed no new code:
conformalised quantile regression for attainable-height intervals (`validation.conformal.
cqr_interval`, Romano et al. 2019, already generic) and the establishment window's g(tau)
smooth term / tau_max truncation (usable directly via `hazard.dlnm.spline_basis` and the
growth module's own outputs). New, tested code, all in `niche/growth.py` unless noted:

* **`gdd_modifier`** (f_T) and **`late_frost_modifier`** (f_F): the concept note names these
  as Liebig-type modifiers ("as in 3-PG" for f_T) but gives no closed form for either.
  Implemented as, respectively, the increasing-saturating (Michaelis-Menten/Monod) and
  decreasing-saturating counterpart of the module's own pre-existing `water_modifier` --
  documented choices, not quoted formulas.
* **`thermal_treeline_modifier`** (f_tl, Eq. 8.3) and **`potential_treeline_elevation`**
  (diagnostic z_tl): kept strictly separate from `gdd_modifier`, per the explicit
  instruction that these are physically distinct mechanisms and must never be merged or
  aliased -- f_tl is a hard ceiling (growth stops), f_T is a heat-sum ramp (growth slows).
  Default T_tl/LGS_min are the concept note's own cited Koerner & Paulsen values; the
  transition widths s_T/s_L are, per the note, "fitted, not assumed" and have no default.
* **`probability_reaches_height`** (`niche/growth.py`): the ensemble evaluator for
  P[H_T >= H_min] = P[sum phi_y >= a*] named explicitly in Sec. 8.1, over a caller-supplied
  ensemble of phi_y series and Chapman-Richards parameter draws.
* **`water_limited_height_fallback`** (`niche/height.py`): the AOA fallback for H* outside
  the training domain. The note names the mechanism (scale by AET/PET) but not its exact
  form; implemented as a simple ratio scaling of a reference H*, capped at 1 -- documented,
  not quoted.

**Still blocked on real data**: the attainable-height quantile model itself needs real
GEDI/canopy-height training data (`antar.io.gee_export.export_structure`, still a stub) to
fit against -- same category as XYLEM's ABC calibration and MNEME's real dieback panel:
the formulas are built and tested, the real fit waits on a data pull not yet done.

## 2026-09-24 — REFUGIUM implementation (concept note Sec. 8.2)

Most of Sec. 8.2 already existed and was tested: competing-hazard combination and cohort
viability (`viability/cohort.py`, matches Eq. 8.5/2.1 exactly, confirmed against the
existing `test_competing_risks_and_viability`), the buffer index and risk-averse refugium
score (`viability/refugia.py`), and climate-analogue matching (`antar.climate.analogs`
-- Mahalanobis distance, novelty percentile, Mahony et al. 2017 -- already covers the
"climate analogues" subsection in full, just housed under `climate/` since it is climate-
space math reused elsewhere too).

One real gap, already fixed in the same pass as the MERISTEM work above:
`refugia.robust_refugium` only implemented criterion (a) of the module's own documented
Eq. 8.6 definition (ensemble viability), silently missing (b) topographic buffering and (c)
area-of-applicability. Extended to accept both as optional arguments (defaulting to `None`
= criterion skipped, so the existing single-criterion call sites keep working unchanged)
and AND them together when supplied.

Confirmed, not rebuilt: the "product assumes independence conditional on covariates ...
tested by comparing the modelled all-cause hazard with observed all-cause plantation
mortality" validation Sec. 8.2 calls for is exactly what `validation.metrics.
calibration_slope_intercept` (or `brier_skill`) already does generically -- no
REFUGIUM-specific code needed, just that usage. Compound events (drought followed by fire)
are explicitly stated as entering "as interactions in Module C" (MNEME), not REFUGIUM's own
job.

New test: `test_viability_composes_with_meristem_height_probability` -- the actual
Module D -> E handoff from Fig. 2 (MERISTEM's `growth.probability_reaches_height` feeding
REFUGIUM's `cohort.viability`) had never been exercised together before; now it is.

REFUGIUM was the last engine with meaningful spec gaps to close without real data. AEGIS
(the decision layer) is the only engine not yet reviewed against the concept note.

## 2026-09-24 — AEGIS implementation (concept note Sec. 10.3)

`decision/optimize.py` already implemented essentially all of Sec. 10.3's core: the
CVaR-robust MILP (Rockafellar & Uryasev 2000 linear form, solved via `scipy.optimize.milp`
-- which uses HiGHS, matching "solved with HiGHS in the skeleton"), the one-option-per-unit,
budget, diversity-cap and eligibility constraints, the optional per-basin water-use cap, and
`evaluate_portfolio` for the out-of-sample CVaR check the note explicitly calls for ("a
large gap between the in-sample and out-of-sample values is a warning" of overfitting to
the scenario sample). Confirmed correct against the exact constraint forms in the text
before adding anything.

Two things the text names explicitly but that were not yet implemented:

* **`efficient_frontier`**: "sweeping lambda from 0 to 1 traces the mean-CVaR frontier, and
  the distance between its ends is the price of robustness, expressed in units of benefit."
  Solves `robust_portfolio` once per lambda in a caller-supplied grid and reports
  `price_of_robustness` as the drop in expected benefit from the first to the last lambda
  solved -- exactly the note's quantity when the grid runs 0 to 1 (the default), well-defined
  but not that literal quantity for any other endpoints (documented in the docstring).
* **`extrapolation_footprint`**: "areas with a large extrapolation footprint are not
  silently planted or excluded: they are flagged, and the plan states how much of its
  expected benefit comes from cells outside the area of applicability." This is a
  *reporting* diagnostic on an already-solved plan, distinct from the `eligible` (e_uj)
  constraint that hard-excludes legally-forbidden or out-of-region options from the
  optimisation itself -- the note describes both mechanisms (an exclusion constraint and a
  separate transparency report), and both are now present.

This closes the last engine with clear, buildable-without-real-data gaps. All five
scientific engines (TOPOHYDRO, XYLEM, MNEME, MERISTEM, REFUGIUM) plus AEGIS have now been
checked against the concept note section by section; remaining work is either (a) fitting
against real data once the hyperion/Earth Engine/EnviDat pulls finish, or (b) genuinely
data-dependent items already flagged in earlier entries (LandTrendr/CCDC segmentation
itself, XYLEM's ABC calibration, the real stacked-hazard fit).

## 2026-09-24 — full Table 4 variable inventory: what's realistically obtainable

Went through Table 4's full variable list systematically rather than continuing to pull
one dataset at a time. `antar.io.gee_export` extended with six more export functions,
each checked for real Earth Engine availability before being written (not assumed):
`export_terrain`, `export_soils`, `export_era5land_forcing`, `export_vegetation_state`,
`export_snow`, `export_land_tenure`. All six submitted as real Drive export tasks
(2000-2024 where the variable has a meaningful yearly value, single-layer where it
doesn't -- terrain, soils, land tenure).

**Two things caught during verification, before submitting real exports:**

* **DEM source substituted.** The concept note specifies Copernicus GLO-30; Earth
  Engine's copy of it (`COPERNICUS/DEM/GLO30`) has a real coverage gap over Armenia --
  checked directly by probing tile ids, only 3 of the ~14 tiles the study bbox needs
  exist there, covering just the southwest corner. Switched to `USGS/SRTMGL1_003`
  (SRTM 30 m), confirmed full bbox coverage (106,256 valid 1 km-sampled pixels, sensible
  70-4978 m range) before using it. A substitution, not a silent assumption.
* **Wind speed bug caught by sanity-checking the actual number, not just that the code
  ran.** First version averaged the u/v wind components separately over the year, then
  took the magnitude of the averages -- since wind direction varies throughout the year,
  opposite-direction days cancel toward zero in that averaging order, and the test point
  came back at 0.11 m/s (implausibly calm). Fixed to compute daily speed (hypot of that
  day's u,v) first, then average the daily speeds; the same point now returns 0.71 m/s,
  a physically reasonable number for a sheltered mountain valley. This is exactly the
  kind of error that runs without crashing and only shows up if the actual value is
  checked against physical expectation -- worth remembering for every other new pull.

**Confirmed NOT realistically obtainable** (genuine gaps, matching the concept note's own
"not every row is available at the start... documented gap, not a silent zero"):
national forest inventory plots, provenance/genetic trial locations, insect/pathogen
outbreak records, treeline field-survey transects (institution/field-survey-only data);
road network/accessibility (no public Earth Engine asset found under any plausible id for
this region -- checked, not assumed); trait databases XFT/TRY (specialist data services
requiring separate registration/API access, not Earth Engine assets).

CO2 concentration trajectories (Table 4's "atmospheric CO2 trajectory" row) are not a
spatial dataset -- they are the standard published SSP concentration pathways (a small,
well-established reference table), not yet added as a `configs/` entry; a follow-up item,
not a data-access gap.

## 2026-09-24 — SSH access to hyperion.wsl.ch

Non-interactive key-based auth check (`ssh -o BatchMode=yes ohanyann@hyperion.wsl.ch true`) failed
with "Permission denied (publickey,keyboard-interactive)": no key is authorized on the remote yet.
No data pull has been attempted. Per the security convention, `ssh-copy-id` must be run by the user
in their own terminal; this session will re-check non-interactively once told it's done.

Later the same day: key-based access confirmed working (non-interactive check exits 0). A read-only
directory listing of `/storage/zilker/chelsaV2/output/FORACCA/Armenia/` (2.8 TB) found:

* `pr`, `tas`, `tasrange`, `tasskew` populated; `tasmax` and `tasmin` present but **empty** (no
  files under either). `tasrange` (diurnal range) and `tasskew` sitting alongside empty
  `tasmax`/`tasmin` is consistent with CHELSA's own published method of reconstructing daily
  tasmax/tasmin from tas+tasrange+tasskew rather than storing them directly -- a plausible
  explanation, not confirmed.
* The populated variables are organised by **GCM x RCP x RCM** (e.g.
  `pr/MOHC-HadGEM2-ES/rcp26/..._CORDEX_WAS-22_GERICS_MOHC-HadGEM2-ES_rcp26_r1i1p1_REMO2015_v1_...`),
  i.e. **CMIP5 GCMs (HadGEM2-ES, MPI-ESM-LR/MR, NorESM1-M) bias-adjusted via CORDEX WAS-22
  regional models (REMO2015, RegCM4-7) under RCP2.6/RCP8.5** -- not the CMIP6/SSP,
  ISIMIP3b + NEX-GDDP-CMIP6 ensemble the concept note's Table 3/`configs/scenarios.yaml`
  specify (5 ISIMIP3b models, SSP1-2.6/2-4.5/3-7.0/5-8.5).
* A `version1/` subtree holds derived indices (tmean/pr seasonal and annual summaries, GDD
  days, r10mm/r20mm, CDD/CWD) per GCM x RCP -- reads as the prior project's own v1 pipeline
  output, not primary forcing.

This is a real mismatch between what the concept note's data plan assumes and what is
actually on hyperion, surfaced rather than silently reconciled -- flagged to the user; no
manifest or pull has been built against this yet pending their decision on how to proceed
(adapt the ensemble design to the CMIP5/RCP/CORDEX data that exists, treat it as legacy/v1
and look for a CMIP6 source elsewhere, or derive tasmax/tasmin locally from tas+tasrange+
tasskew).

**Decision: adapt.** `configs/scenarios.yaml` rewritten to the CMIP5/RCP/CORDEX ensemble
actually on hyperion (4 GCMs, RCP2.6/RCP8.5, GERICS-REMO2015 and ORNL-RegCM4-7 regional
models) rather than the concept note's ISIMIP3b/NEX-GDDP-CMIP6/SSP design. No CMIP6/SSP
source was substituted -- `secondary_source` is now `null`. Two things flagged in the new
config rather than resolved: (1) the CMIP6-specific ECS/hot-model screening (Zelinka et al.
2020; Hausfather et al. 2022) does not carry over to these CMIP5 models and needs its own
sourced values before an equivalent policy can be written; (2) filenames on hyperion
(`*_sim_hist_ba.nc`, `weights_pr.nc`, `*_qm_transfer_plot.png`) suggest this RCM output is
already quantile-mapping bias-adjusted upstream, which may conflict with Sec. 5.1's
assumption that raw GCM output still needs `antar.climate.bias.quantile_delta_mapping`
applied against ERA5-Land -- needs confirming before TOPOHYDRO's bias-adjustment step is
wired to this data source. `tasmax`/`tasmin` are still absent from hyperion; not resolved
here.
