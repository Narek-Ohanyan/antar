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

## 2026-09-24 — real CHELSA-daily found; historical/future data source split decided

User pushed back on the hyperion CORDEX pull being treated as "the CHELSA data": correctly --
hyperion's tree is CORDEX regional-model output bias-adjusted *using* CHELSA as a reference,
not CHELSA's own product. The real public CHELSA archive lives at
`os.unil.cloud.switch.ch/chelsa02/chelsa/global/` (found via the actual "Download" link's
target on chelsa-climate.org, not guessed), as Cloud-Optimized GeoTIFFs, one global file per
variable per day, openable directly over HTTP range requests (GDAL `/vsicurl/`) -- confirmed
by opening one file and reading only Armenia's window: 0.3 s, no full-file download.

**CHELSA-BIOCLIM+** (bonus find while checking the site): bioclimatic variables including
GDD0/5/10, growing-season length via the **TREELIM methodology**, frost-change frequency, VPD,
PET, wind, radiation, for 1981-2010 and, critically, future periods (2011-2100) broken out by
GCM: **GFDL-ESM4, IPSL-CM6A-LR, MPI-ESM1-2-HR, MRI-ESM2-0, UKESM1-0-LL** -- exactly the five
ISIMIP3b models the concept note specifies. Not pulled yet; noted for when MERISTEM's
treeline/growing-season fitting needs it.

**Decision (user): CHELSA-daily for the historical reference, CORDEX (hyperion) for future
projections** -- not one source mixed across the whole timeline. Implemented as two separate
pulls rather than reconciled into one.

**A real bug caught by checking the actual data, not just that the pull ran**: the first
CHELSA-daily attempt assumed filenames were `CHELSA_{var}_{MM}_{DD}_{YYYY}`. For days 1-12 of
each month this silently returned *valid-looking data for the wrong date* (no error, no
crash); for days 13-31 it correctly 404'd, which is what actually exposed the bug -- a partial
month of real numbers would have looked completely plausible otherwise. Diagnosed by comparing
the actual mean temperature of `CHELSA_tas_01_10_2000` (10.2 degC) against what January and
October should look like in Armenia: 10.2 degC is unambiguously October, not January -- proving
the field order is day-then-month (`DD_MM_YYYY`), not month-then-day. Fixed, re-verified
(Jan 1 2000 now reads 0.9 degC, genuinely cold), and only then run at scale. This is exactly
the "check the number, not just that the code ran" lesson from the earlier wind-speed bug,
now caught a second time on a much larger, easier-to-miss pull.

**Reference-period gap caught before, not after, the full pull**: the first attempt pulled
2000-2024 (matching `calibration_period`), but the concept note's own `reference_period` for
standardising every anomaly and index is 1991-2020 (`scenarios.yaml`) -- 2000-2024 does not
cover it. User caught this too. Extended to 1979-2024 (tas/tasmax/tasmin) / 1979-2019 (pr,
matching precipitation's real coverage ceiling) -- 1979 also matches the reanalysis-era
convention (ERA5 itself starts 1979), not an arbitrary round number.

**Confirmed variable coverage** (checked directly, year by year, not assumed from the
dataset's blanket "1941-2025" label): `tas`/`tasmax`/`tasmin` complete 1979-2024;
`pr` complete 1979-2019, genuinely absent 2020-2024 (matches the archive's own
"active (incomplete)" status, verified: zero objects under `pr/2020/` through `pr/2024/`).

Pulling via `/vsicurl/` windowed reads directly into local compressed `.npz` files (no
intermediate full-resolution global downloads); `configs/scenarios.yaml` updated with the
full historical-source block (dataset, URL, citation, license, per-variable coverage) and the
CORDEX `gcms:` block re-scoped to "future projections only." CORDEX's future portion
(2025-2100) extracted the same way as before (ncks, stride-6, hyperion-side) for all four
CORDEX variables (pr, tas, tasrange, tasskew) -- superseding the earlier 2006-2024 historical
CORDEX extracts, which are no longer the right thing to use now that CHELSA-daily covers that
period directly. Future `tasmax`/`tasmin` still have no direct CORDEX source and would need
deriving from `tas`+`tasrange` -- flagged in `scenarios.yaml`, not yet implemented.

## 2026-09-24 — S2-VHM pull complete; manifest convention for Drive-hosted data

All 9 years (2017-2025) of the S2-VHM download-upload-to-Drive pipeline finished.
Registered as `antar.io.manifest.ManifestEntry` records in `configs/manifests/s2_vhm.yaml`
-- a new, tracked (not gitignored) location for manifests, distinct from `data/` (raw
content, gitignored). Each entry: `local_path=None` (genuinely not stored locally, by
design -- downloaded, uploaded to Drive, deleted), `checksum_sha256` from the SHA-256
computed before upload (still a real, verifiable property of the file's content even
though the file itself lives in Drive), and the Drive file id recorded in `notes`.

hyperion came back after the earlier outage (the user reconnected on their end). Killed
the stale hung SSH process (`ServerAliveInterval`/`ServerAliveCountMax` added to the retry
command this time, so a future silent drop is detected and killed rather than hanging
indefinitely), re-pulled the already-completed `pr` extract (108,474,482 bytes, matches
exactly -- it had finished on hyperion's side before the outage, only the local copy was
lost), and restarted the `tas`/`tasrange`/`tasskew` extraction.

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

## 2026-09-28 — CHELSA-daily pull completed; three real bugs hit and fixed along the way

The extended CHELSA-daily pull (`scripts/pull_chelsa_daily.py`, 1979-2024 tas/tasmax/tasmin,
1979-2019 pr) finished: all four variables complete, 0 confirmed-missing (404) anywhere,
16802/16802/16802/14975 days respectively. Registered in `configs/manifests/chelsa_daily.yaml`
(all four checksums verified against the actual local files via `antar.io.manifest.verify_entry`
before being written down, not copied blind from the pull's own summary).

Getting there took three real, self-caused failures on top of the disk-full kill already logged
under 2026-09-24/25 -- recorded here in full rather than smoothed over, since each one changed
the script's actual behaviour:

1. **No incremental checkpointing.** The original script only wrote to disk once, at the very
   end of each variable's ~16800-day fetch. A disk-full kill mid-`tas` (see the earlier entry)
   lost all 16500 already-fetched days, because nothing had been persisted yet. Fixed: `_checkpoint`
   now runs every 2000 days fetched, so an interruption only costs the work since the last
   checkpoint, not the whole variable.
2. **O(n^2) resume.** The first attempt at checkpointing still indexed the loaded `.npz` as
   `saved["data"][i]` inside a loop over all prior days -- `NpzFile.__getitem__` re-decompresses
   the *entire* stacked array on every call, so resuming a 14000-day checkpoint triggered ~14000
   redundant full-array decompressions and never visibly progressed (looked hung; confirmed still
   running via `ps`/CPU%, not actually stuck). Fixed by decompressing the array once outside the
   loop.
3. **Corrupted a good checkpoint by killing a redundant rewrite mid-write.** After fixing (1) and
   (2), a resume with nothing new to fetch for `tas` still re-ran the full stack+compress+write
   every time (wasted, but at the time believed harmless). Killing that redundant write mid-flight
   corrupted the previously-complete, previously-good `tas` file (`np.savez_compressed` writes
   straight to the target path, so a partial write left a truncated, unreadable `.npz` --
   `zipfile.BadZipFile: File is not a zip file`). `tasmax`/`tasmin`/`pr` were untouched and verified
   still readable; only `tas` had to be re-fetched from scratch. Fixed two ways: (a) skip the
   stack/compress/write step entirely when there is nothing new to fetch, so a fully-complete
   variable resumes instantly instead of redoing needless work; (b) make the write itself atomic --
   write to a temp file, then `os.replace()` into place, so an interruption can never leave the real
   file half-written again. First attempt at (b) had its own bug: `np.savez_compressed` silently
   appends `.npz` to any filename that doesn't already end in it, so a `*.npz.tmp` temp path was
   actually written to `*.npz.tmp.npz`, and the atomic rename then failed looking for a file that
   didn't exist under the name it expected. Fixed by naming the temp file `*.tmp.npz` instead, and
   smoke-tested standalone (dummy array, real `_checkpoint()` call, verified the file lands at the
   right path and reloads) before trusting it against the real multi-GB pull again.

Also hit, not a script bug: a stretch of `tas` days (~Feb-Mar 1996) failed with a mix of DNS
resolution failures and connection timeouts during one run -- a genuine local-network outage on
this machine, not a CHELSA-server or code problem. Retried cleanly on the next run once the
network recovered; the `missing_404` vs `failed_read` split added earlier this same day (see the
prior 2026-09-24 CHELSA entry's follow-up) is exactly what made that retry automatic rather than
requiring a manual list of which days to redo.

Separately: `hyperion.wsl.ch` failed to resolve via DNS entirely during this window (`nodename
nor servname provided`) -- confirms it is only reachable from the office network/VPN, not the
public internet. The CORDEX future-projection extraction (`tasrange`/`tasskew` for
MOHC-HadGEM2-ES/rcp26, plus the other 7 of 8 GCM x RCP combinations) stays blocked until that
access is available again; not attempted from outside it.

## 2026-09-30 — real GCP project ID found; six-plus GEE exports checked and registered; GHCN-Daily pulled

**The Earth Engine project-ID blocker (flagged since 2026-09-24) is resolved.** Reconstructing
OAuth credentials from the stored refresh token (same pattern as the earlier Drive-access work)
and calling the Cloud Resource Manager API with no project produced a 403 naming project number
`517222506229` -- this looked like a real discovery but was a red herring: it's literally
hardcoded in the `earthengine-api` package's own source (`ee/__init__.py`, `oauth_project =
'517222506229'`) as the SDK's own default/shared quota project, which `oauth.is_sdk_project()`
explicitly rejects as a valid EE project. The user supplied the real project id directly
(`pure-highlander-495708-a9`, number `270993636879`); `ee.Initialize()` with it works. **Real
constraint surfaced on first use**: this project has exceeded its Earth Engine noncommercial
compute quota and is now in *restricted mode* -- batch exports still submit but may queue
(`READY`) rather than run immediately; something to watch for every future EE task, not a one-off.

**Six-plus exports submitted in an earlier session (2026-09-24 through 09-27) were checked
against `ee.data.getTaskList()` rather than assumed finished**: `terrain`, `soils`,
`era5land_forcing`, `snow`, `land_tenure`, `vitality_composites`, `landtrendr_segmentation`,
`structure`, and `disturbance_ancillary` all show COMPLETED, with real Drive destination URIs
(all in the `antar_gee_exports` folder, ~85 GB across 38 files, none downloaded locally --
registered in `configs/manifests/gee_exports.yaml` with Drive's own MD5 in `notes` rather than a
freshly computed SHA-256, since hashing multi-GB files just to duplicate a check Drive already
does would cost real bandwidth/storage for no integrity benefit).

**`vegetation_state` had genuinely failed**, not just gone unchecked: `Image.select: Parameter
'input' is required and may not be null`. Root cause: `MCD12Q2` (MODIS land-surface phenology)
only has data from 2001 onward (confirmed directly against the collection: first image
2001-01-01, last 2025-01-01, 25 images) -- the export's year loop started at 2000 to match the
other yearly exports, so `.filterDate(...).select(...).first()` returned null for 2000 and the
whole batch crashed on it. Fixed in `src/antar/io/gee_export.py`: phenology bands are now
genuinely omitted for years before 2001 rather than filled with a fabricated masked value (LAI,
which has its own longer coverage, is unaffected). Resubmitted; queued under restricted mode as
of this writing, not yet actually run -- check again once cleared, don't assume it completed.

**GHCN-Daily Armenia pull** (`scripts/pull_ghcnd_armenia.py`, new): all 53 real Armenia-country-code
stations, PRCP/TMAX/TMIN/TAVG/SNOW/SNWD, parsed directly from the classic fixed-width `.dly`
format (the `by_station/*.csv` path some documentation suggests returned 404 for the one station
tried -- not real for this dataset, or differently named; not pursued further once `.dly`
confirmed working). Units converted per GHCN's documented convention (PRCP/TMAX/TMIN/TAVG are
tenths of their stated unit; SNOW/SNWD are not) -- verified before trusting the full pull, the
same discipline as the CHELSA date-order bug: Yerevan July 2000 reads 26-30 degC, January 2000
reads 0.7-6 degC, both correct for the location; nonzero precipitation values are plausible daily
mm magnitudes. Only QFLAG-blank (NOAA-QC-passed) values kept; flagged values dropped rather than
included. 1,325,996 observations, 5.0 MB compressed. Storage-minimal by design end to end: parsed
one station at a time (never all 53 in memory together), written directly to gzip, uploaded to
the same `antar_gee_exports` Drive folder, local staging deleted immediately after upload --
nothing but a small JSON summary remains on local disk. Registered in
`configs/manifests/ghcnd_armenia.yaml`. This closes the "Armhydromet stations are
institution-only" item from the earlier manual-data assessment -- they were not; that assessment
was wrong, corrected here rather than left standing.

**Copernicus GLO-30's Armenia gap re-checked against the primary source, not the earlier
research pass's guess.** That guess (licensing exclusion lifted ~Dec 2023, so AWS/CDSE should
have full coverage now) was wrong. Checked directly: AWS `copernicus-dem-30m` bucket listing
(`?list-type=2&prefix=...` against all 16 tiles the study bbox needs) and Copernicus Data Space
Ecosystem's own STAC API (`stac.dataspace.copernicus.eu`, collection `cop-dem-glo-30-dged-cog`,
bbox search) both return exactly the same 3 tiles as Earth Engine's `COPERNICUS/DEM/GLO30`
(`N38/E043`, `N38/E044`, `N39/E043`). Armenia has no Copernicus GLO-30 coverage anywhere in the
Copernicus system right now, confirmed from the authoritative catalog itself -- not a stale-copy
problem. SRTM 30 m (`USGS/SRTMGL1_003`, already in `export_terrain`) remains the correct choice;
this closes the question rather than leaving a plausible-sounding but unverified claim standing.

**CO2 concentration pathways pulled** (`scripts/pull_co2_pathways.py`): RCMIP v5.1.0 confirmed
real and directly downloadable, no registration. Filtered the 21 MB multi-model file down to
exactly this project's scenarios -- both the CMIP5/RCP pathways (rcp26/45/60/85, backing CORDEX)
and the CMIP6/SSP pathways (ssp126/245/370/585, backing everything else), plus historical -- into
`configs/co2_concentration_pathways.csv` (9 rows). Real numbers now available for the RCP<->SSP
correspondence question ROADMAP.md S0 raises: rcp26 and ssp126 are close (~421 vs ~446 ppm by
2100) but not identical, confirming the concern was real, not hypothetical.

**GBIF occurrence data pulled** (`scripts/pull_gbif_occurrences.py`): 524 real presence records
for the 7 target species named in `configs/species_traits.csv` (Fagus orientalis, Carpinus
betulus, Quercus macranthera, Quercus iberica, Pinus kochiana, Juniperus polycarpos, Juniperus
excelsa), via GBIF's public search API, no registration. A real data-quality bug caught before
trusting the output: `country=AM` alone let through a handful of `(0,0)` "null island" records --
a classic bad-georeferencing artifact, not a real Armenian occurrence. Fixed with an explicit
bbox filter (matching `configs/study_area.yaml`) rather than trusting the country filter alone;
4 of 528 raw records dropped.

**CHELSA-BIOCLIM+ access confirmed, not yet pulled**: `chelsa/global/bioclim/` is a real sibling
directory to `chelsa/global/daily/` on the same server, listed directly (bio01-19, cltmax/mean/
min/range, cmimax/mean/min/range, fcf, fgd, gdd0/5/10, and more) -- same `/vsicurl/` access
pattern already proven for CHELSA-daily should work here too. Deliberately not pulled in this
batch: it is a properly-scoped job on its own (many variables x historical + 5 GCMs x 3 SSPs x 3
future windows), not a quick add alongside the smaller items done here.

**Road/accessibility data: Geofabrik's binary downloads are broken from this environment,
pivoted to the Earth Engine alternative.** `download.geofabrik.de/asia/armenia-latest-free.shp.zip`
(and the `.osm.pbf` equivalent) return a persistent 301 redirect loop on every attempt (HEAD and
GET, with and without a browser User-Agent) -- the HTML listing page itself loads fine, so this
is a real server/proxy behavior for this specific binary-download path from this network, not a
missing-data problem. Not chased further once the Malaria Atlas Project's global travel-time-to-
cities layer (Weiss et al. 2018) confirmed working: `Oxford/MAP/accessibility_to_cities_2015_v1_0`
on Earth Engine, verified directly (single `accessibility` band, values 0-517 minutes over the
Armenia bbox -- physically sensible for the terrain, not just "the asset loaded"). Exported to
Drive (`antar_accessibility_to_cities`); queued under the project's restricted-mode compute quota
as of this writing, register once it actually completes.

**Restricted-mode quota observed in practice**: both this accessibility export and the
re-submitted `vegetation_state` export (see above) sat in `READY` (queued, not running) for the
duration of this session after submission -- the noncommercial-quota restriction is throttling
throughput, not just a one-line warning. Worth checking task state again before assuming any
future EE export has actually run.

## 2026-09-30 (cont.) -- decisions applied; CHELSA-BIOCLIM+ nodata bug caught before the pull; new GCP project

User's four decisions from the interactive question round: (1) run CORDEX/CMIP5 and CMIP6/
ISIMIP3b as two separate labelled future ensembles, not one -- applied in `configs/scenarios.yaml`
(see the commit "Run CORDEX and CMIP6/ISIMIP3b as two separate labelled future ensembles"); (2)
pull CHELSA-BIOCLIM+ now; (3) skip the tree-ring author outreach; (4) create a second GCP project
for fresh Earth Engine quota rather than wait out the first project's restricted mode.

**Second GCP project (`antar-armenia-2`) set up and verified working.** Registration for Earth
Engine noncommercial use is a separate step from just creating the project and enabling the API
-- hit both intermediate states in order (API-not-enabled, then registered-but-quota-check still
pending) before `ee.Initialize()` succeeded cleanly with 0 existing tasks (fresh quota, confirmed
via `ee.data.getTaskList()`). The two exports stuck in `READY` under the first project's
restricted mode (`vegetation_state`, `accessibility_to_cities`) were cancelled there and
resubmitted on the new project (`PWUJYIVOR4T4BLZYDO46N4H3`, `USCD5PM5QOL75D3Y7OZIIRUE`) rather
than left to potentially never run.

**Real nodata bug caught before trusting the CHELSA-BIOCLIM+ pull at scale, not after.** First
smoke-test fetch (gdd0, historical) returned a maximum of 214748364.7 -- exactly
2147483647 (INT32_MAX) x 0.1, i.e. the source raster's own declared nodata sentinel leaking
through unmasked because `fetch_one` scaled the raw value without checking for nodata first.
Confirmed via the raster's own metadata (`dtype=int32, nodata=2147483647, scale=0.1`) before
fixing, not guessed. Fixed in `scripts/pull_chelsa_bioclim.py`: mask nodata to NaN before scaling.
Re-tested the same fetch afterward: max 5955.9, mean 3302, 2 genuinely-nodata pixels (plausible
values for annual growing-degree-days in Armenia) -- confirmed fixed by checking the actual
numbers again, not just that it no longer errored.

**Checked whether the already-completed CHELSA-daily pull has the same latent bug -- it does not,
verified, not assumed.** CHELSA-daily's uint16 nodata sentinel (65535) would scale to exactly
6553.5; scanned all four pulled variables' full arrays (tas, tasmax, tasmin: ~2.08 billion pixels
each; pr: ~1.85 billion) for that exact value: zero matches in all four. Armenia's bbox is fully
inland with no missing-data pixels in this product, so the same unmasked-nodata code path in
`pull_chelsa_daily.py` never actually fired. The already-pulled data does not need to be re-pulled.
Added the same defensive nodata mask to `pull_chelsa_daily.py` anyway, for correctness if the bbox
ever changes -- not because it fixed an observed problem there.

**Separately, a real but different data characteristic was found while investigating, not a code
bug**: `pr`'s already-pulled data has 7810 pixel-days (out of ~1.85 billion, ~0.0004%) exceeding
1000 mm/day, scattered across ~6700 distinct pixel locations, clustered in complex mountain
terrain. None of these match the nodata-sentinel value (6553.5), ruling out the same bug. This
reads as a known CHELSA artifact in orographically complex terrain rather than corruption --
documented as a caveat on the `pr` manifest entry, not silently left unmentioned, and not
"fixed" since these may be genuine (if unusually extreme) modelled values rather than clearly
wrong ones.

**CHELSA-BIOCLIM+ pull scoped to ~37 of the archive's 74 real variables** (confirmed by listing
the bucket): bio01-19, gdd0/5/10, gddlgd0/5/10, gsl/gsp/gst, fcf/fgd/lgd, vpdmean/max, petmean/max,
sfcWindmean, rsdsmean -- the ones that map to something an ANTAR engine actually uses (MERISTEM's
niche/treeline/frost terms, XYLEM's water-stress terms, TOPOHYDRO's water-balance and
ERA5-Land cross-checks). Dropped the Koppen-Geiger bins, npp/swe/swb, and the min/range variants
of variables already covered by mean/max -- a deliberate scope decision against "minimize
storage," not an oversight. Historical (1981-2010) is one file per variable; future
(2011-2040/2041-2070/2071-2100 x 5 GCMs x 3 SSPs) adds 45 more per variable -- reused the
CHELSA-daily pull's proven checkpointing/atomic-write pattern directly.

## 2026-09-30 (cont.) -- CHELSA-BIOCLIM+ pull completed; a second checkpoint bug caught; both stuck GEE exports finished

**Second real bug in the same pull, caught the same way -- by watching it actually run, not just
checking it started.** After 24 of 37 variables completed cleanly, the pull crashed on `vpdmean`
with `ValueError: need at least one array to stack`. Root cause: `vpdmean` (and, it turned out,
`vpdmax`, `petmean`, `petmax`, `sfcWindmean`, `rsdsmean`) has ONLY historical data in the real
archive -- confirmed by listing the bucket directly (`chelsa/global/bioclim/vpdmean/` has just
`1981-2010/`, no future-period subdirectories at all). All 45 future-scenario fetches for these
variables correctly 404. Because 404s resolve faster than the one slow historical-file read
across 8 parallel threads, the *first* checkpoint boundary (20 completions) could be reached with
zero successes collected yet, and `np.stack([])` doesn't accept an empty list. Fixed
`_checkpoint` in `scripts/pull_chelsa_bioclim.py` to write progress either way but only touch the
`.npz` once there is real data to stack; smoke-tested both the empty- and non-empty-`done` paths
standalone before relaunching. Also hardened `main()`'s checksum step against a variable that
might end up with zero scenarios at all (not currently triggered, but now correct if it ever is).

**Full pull completed after the fix**: 30 of 37 variables with full historical+future coverage
(46/46 scenarios each), 6 confirmed historical-only as above, 0 failed-reads anywhere. Spot-check
values are physically sensible -- bio01 (annual mean temp) ranges -14.0 to 24.6 degC across
Armenia's historical+future scenarios, gdd0 0.5-9007.6 (2 genuinely-nodata pixels, matching the
smoke test). 223 MB total, kept local rather than moved to Drive -- small enough that the
cloud-first convention doesn't buy anything here the way it does for the multi-GB GEE rasters.
Registered as 37 entries in `configs/manifests/chelsa_bioclim.yaml`
(`scripts/register_chelsa_bioclim.py`, reads the pull's own manifest_rows.json rather than
re-deriving anything), all 37 verify_entry-confirmed.

**Both exports stuck under the first GCP project's restricted mode finished within minutes of
being moved to the second project**: `antar_vegetation_state` (919 MB) and
`antar_accessibility_to_cities` (5.7 MB), both COMPLETED, both registered in
`configs/manifests/gee_exports.yaml` (now 40 entries total). Confirms the second-project
workaround actually solves the throttling problem, not just theoretically.

## 2026-09-30 (cont.) -- hyperion deprioritised; ISIMIP3b access confirmed working end to end

User's direction: hyperion is currently unreachable and stays deprioritised -- proceed with
CHELSA and other already-available sources now, treat CORDEX/hyperion as a later addition once
office-network access is available again.

**ISIMIP3b confirmed genuinely working, not just "probably works"** -- the earlier research
pass's `isimip-client` recommendation had two real problems, both found and fixed by actually
running it rather than trusting the package name:

1. The pip-installed `isimip-client` (v2.0.2)'s own `cutout()` method builds a request shape
   (`{'task': 'cutout_bbox', 'bbox': ..., 'paths': ...}`) that the *live* Files API v2 rejects
   outright (`400: 'operations' field required`) -- the client library is out of sync with
   ISIMIP's current API. Fixed by reading the API's own root endpoint
   (`files.isimip.org/api/v2/`, which lists its real `operations` schema) and posting the
   correct shape directly via the client's lower-level `post_job()`: `{'paths': [...],
   'operations': [{'operation': 'cutout_bbox', 'bbox': [...]}]}`.
2. **A real bbox-axis-order bug, caught by checking the actual output coordinates, not the job's
   success status.** First attempt used `bbox=[38.8, 41.4, 43.4, 46.7]` (a [south, north, west,
   east]-style guess) -- the job "succeeded" and produced a file, but opening it with netCDF4
   showed lat 43.75-46.25, lon 39.25-41.25: northern Georgia/southern Russia, not Armenia. The
   real order is `[lon_min, lon_max, lat_min, lat_max]`. Corrected to `[43.4, 46.7, 38.8, 41.4]`,
   re-verified: lat 39.25-41.25, lon 43.75-46.25 (genuinely Armenia), daily tas -30.7 to 32.1
   degC over 2015-2020 (GFDL-ESM4/ssp126) -- physically plausible, not just "a number came back."
   This is exactly the class of error ("the code ran without error" != "the output is right")
   this project has been burned by twice before (CHELSA's date-field order, the wind-speed
   averaging-order bug) -- caught the same way, by checking the real numbers before trusting
   the mechanism at scale.

**Consequence for the ensemble design**: ISIMIP3b is now the confirmed PRIMARY CMIP6 source
(exactly the concept note's 5 GCMs and ssp126/370/585), not NASA/GDDP-CMIP6 -- separately
confirmed (queried the EE collection directly) that NASA/GDDP-CMIP6's EE ingestion only has
historical/ssp245/ssp585, missing ssp126 and ssp370 entirely. NASA/GDDP-CMIP6 is now documented
as a secondary, broader-GCM (34 models) cross-check restricted to the two scenarios it actually
has, not the primary. `configs/scenarios.yaml` updated accordingly. Full ISIMIP3b pull (4
variables x 5 GCMs x 3 SSPs x 9 decade-chunks = 540 cutout jobs) not yet run -- the next real
data-pull task, needs a proper checkpointed script given the scale.

Also launched a research pass (delegated, not yet returned) specifically checking Armenia's
National Forest Inventory/Hayantar, provenance trial networks (EUFORGEN), pest/pathogen outbreak
records (EPPO), and treeline transect data -- the "institution-only" label on these had never
actually been checked against a live search in this session, only carried over from the very
first pre-research assessment. Also checking TRY and XFT trait database registration processes
live rather than from prior recall.

## 2026-09-30 (cont.) -- real XFT trait data replaces the species_traits.csv placeholders

User registered on xylemfunctionaltraits.org themselves and downloaded the full XFT export
(4015 records, `data/xft/XFT_full_database_download_20260930-151542.csv`). **A password was
pasted into chat during this exchange -- not used, not stored, not referenced anywhere in the
repo; the user was told directly to change it.** No credential handling was needed for this step
in any case: the user downloaded the file themselves through the site's own interface.

No exact-species XFT records exist for any of the four target taxa (Fagus orientalis, Carpinus
betulus, Quercus macranthera, Q. iberica, Pinus kochiana, Juniperus polycarpos, J. excelsa) --
expected, XFT skews toward well-studied European/American species. Used real congeneric proxies
instead, chosen and reported explicitly rather than picked by convenience:

* **mesic_diffuse_porous_broadleaf**: Fagus sylvatica (n=73, the closest relative of the target
  F. orientalis) + Carpinus betulus (n=8, an exact match for one of the two target species).
* **ring_porous_oak**: five European white-oak-group deciduous oaks (Q. petraea, robur,
  pubescens, frainetto, humilis; n=37 combined) -- same taxonomic section and similar
  submediterranean climate niche as the targets.
* **pine**: Pinus sylvestris alone (n=51) -- P. kochiana is taxonomically very close to (at times
  treated as a variety of) P. sylvestris, the single best available analog.
* **juniper_arid_conifer**: real course-correction caught mid-analysis, not just a clean lookup.
  First pass included both J. thurifera (n=1) and J. communis (n=11) for a larger sample; the
  blended median came out to -5.96 MPa P50, well short of the -9 MPa the placeholder had
  guessed for an arid-zone Irano-Turanian juniper. Pulled the full genus's P50 distribution
  across 24 species (-1.67 to -14.2 MPa) before accepting either number: J. communis sits at the
  mesic/widespread end of that range, while J. thurifera -- a high-elevation Mediterranean-
  mountain juniper -- sits in the ecologically appropriate arid/cold niche for the real target
  species. Corrected to use J. thurifera alone (P50 -9.27 MPa, slope 17.2 %/MPa from its own
  P12/P88) despite n=1, on the reasoning that ecological similarity outweighs raw sample count
  for a single congeneric proxy -- diluting it with a poorly-matched but larger sample would
  have been worse, not better. This is the same class of judgment call as the earlier tas/tasmax
  variable-scope decisions, made explicit rather than left implicit in a script.

Slope (% MPa-1) computed via `antar.hydraulics.vulnerability.slope_from_p12_p88` -- the project's
own existing Pammenter & Van der Willigen (1998) formula, reused directly rather than
re-derived, from real P12/P88 pairs (n=44/16/31/1 respectively). gmin25 updated from XFT's
Gsmin field (mol->mmol conversion) only where real data existed (pine, n=4); the other three
groups keep their original placeholder gmin, explicitly labelled as such in the `status` column
rather than silently left ambiguous. psi_close_mpa, tp_c, lethal_plc and capacitance_mmol_m2_mpa
are unchanged from the original placeholders -- no clean, unit-compatible XFT field exists for
any of them (capacitance in particular is reported per sapwood volume in XFT, not the per-leaf-
area molar basis this project's schema uses, and converting would need wood density and
sapwood:leaf-area ratio XFT doesn't reliably give per record -- not attempted rather than guessed).

A real bug was caught and fixed before committing: the first version of the extraction script
fell back to a blank cell (not the original placeholder) for `gmin25_mmol_m2_s` whenever a group
had no real XFT Gsmin data -- would have written empty cells into a CSV real code will eventually
load. Fixed by giving the fallback dict its own gmin value per group.

`configs/species_traits.csv`'s `status` column now records, per group, exactly which fields are
real vs. still placeholder and the real fields' sample sizes and source species -- confirmed no
code currently branches on this column's exact string (only `stand_defaults.yaml`'s separate
placeholder guard does that), so changing its format broke nothing. Registered the source XFT
file in `configs/manifests/xft.yaml`.

## 2026-09-30 (cont.) -- TOPOHYDRO's first real fit: monthly lapse rate and precip gradient

User's decision (asked as an explicit interactive choice): pivot to real model fitting now
rather than pulling more data first (TerraClimate/opendata.am/EPPO deferred until a specific
engine's fit actually needs them).

`scripts/fit_topohydro_lapse_rate.py` (new): the first real fit anywhere in ANTAR, using
`antar.climate.downscale.fit_monthly_lapse_rate` / `fit_precip_elevation_gradient` (both already
implemented and presumably unit-tested against synthetic data since TOPOHYDRO's original build)
against real GHCN-Daily station observations for the first time. Re-downloaded the station data
from Drive (it was deleted locally after the original pull, by design) and verified its checksum
against `configs/manifests/ghcnd_armenia.yaml` before trusting it. Fit window: the concept note's
own `reference_period` (1991-2020), not the stations' full record -- consistent with every other
anomaly/index in this project being standardised against that one period. TAVG used directly
where GHCN reports it (37152 station-days); where only TMAX/TMIN exist, TAVG = (TMAX+TMIN)/2 is
used as a standard meteorological convention, not an invented value, and reported separately
(8513 station-days) so it's never hidden which numbers came from which source.

**Results are physically valid, checked against real-world expectation, not just "the fit ran":**
monthly lapse rate ranges -0.00486 to -0.00760 K/m (-4.9 to -7.6 K/km) with exactly the expected
seasonal pattern -- steeper (more negative) in summer (June: -7.6 K/km) from stronger convective
mixing, shallower in winter (November: -4.9 K/km) from temperature inversions. R-squared 0.32-0.70
across months, reasonable for station-based climate regression with real microclimate noise.
Precipitation's elevation gradient is a genuine, informative near-null result: R-squared
essentially 0 (0.000-0.011) every month -- elevation alone barely predicts Armenia's
precipitation pattern, which is itself a real scientific finding (not a bug) and confirms
`downscale_precipitation`'s `exposure_factor` correction term (wind exposure, rain-shadow) will
matter far more than the base elevation gradient for this country's terrain.

A minor bug caught before the output saved: `r_squared_per_month`'s diagnostic returned
`np.float64` values, which `yaml.safe_dump` can't serialize -- fixed with an explicit
`float()` cast. Saved to `configs/fitted/topohydro_lapse_rate.yaml`.

## 2026-09-30 (cont.) -- MERISTEM's adult-niche model built from scratch and fit per species

Checked `antar.niche.adult` before starting and found it held only `boyce_index`, the
validation metric -- not the presence-background model the concept note's Module D
(`docs/concept_v2/sec7_moduleDE.tex`, `sec:D`) actually specifies: "For each species we fit
a penalised presence-background model on physically meaningful predictors (CWD, growing
degree days, winter minimum, VPD, soil), with the background drawn from a target-group sample
to absorb collection bias, and score it under spatial-block CV with the continuous Boyce
index." Nothing implementing the model itself existed. Added `fit_presence_background` (L2-
penalised logistic regression -- the standard practical equivalent of a presence-background
point-process/MaxEnt-style model, Fithian & Hastie 2013) to `antar/niche/adult.py`, reusing
`antar.validation.splits.block_kfold` for the spatial-block CV rather than writing a new one.

**Real target-group background, not uniform random**: pulled 3000 GBIF Plantae records for
Armenia (excluding the 7 target species) as the background sample -- the concept note is
explicit that this must be a target-group sample "to absorb collection bias," and a uniform-
random background would not do that (collectors who record the target species also tend to
record other plants on the same trips; using those same trips as background controls for
exactly that spatial sampling bias).

**Predictors**: bio06 (winter minimum), gdd5, vpdmean from the real CHELSA-BIOCLIM+ pull;
clay/sand/silt/soc from the real SoilGrids export (downloaded transiently from Drive,
extracted at point locations, deleted immediately after -- storage-minimal, same pattern as
every other transient-download step this session). A real CRS bug was caught before trusting
the soil extraction: `soils.tif` was exported in the master grid's CRS (EPSG:32638), not
WGS84 -- checked via `src.crs` rather than assumed, and occurrence lat/lon reprojected with
`rasterio.warp.transform` before sampling; passing raw WGS84 coordinates into a UTM raster
would have silently sampled the wrong pixels. CWD itself is not yet a pulled product (needs
TOPOHYDRO's full water balance); used `CWD_approx = petmean - bio12` (annual, CHELSA-BIOCLIM+)
as a documented simplified stand-in, named accordingly.

**A real methodological error caught and fixed, not just a code bug**: the first version of
this fit pooled all 7 target species into one presence class against the shared background.
It ran without error and produced a number (mean Boyce index 0.25 across 5 spatial-block
folds, one fold strongly negative at -0.66) that looked plausible enough to almost accept.
Re-reading the concept note's own wording -- "for each species we fit a..." -- made clear this
was wrong: beech, oak, pine and juniper have genuinely different climate niches, and pooling
them blurs exactly the signal a niche model exists to recover. Refit per species instead.

**Per-species results, real and honestly mixed, not all positive**: Fagus orientalis (n=70)
mean Boyce 0.433; Juniperus polycarpos (n=109) 0.345; Juniperus excelsa (n=113) 0.341 --
consistent, moderate-to-good signal, and the two junipers' near-identical results make sense
for sister species with similar ecology. Carpinus betulus (n=69) 0.254 (one negative fold).
Quercus macranthera (n=124) 0.144 (weak). Quercus iberica (n=31) -0.067 -- essentially no
better than random, plausibly a real small-sample/predictor-set limitation rather than
forced into a falsely-positive number. Pinus kochiana (n=8) explicitly skipped, not forced:
8 presence points against 8 features is exactly the regime where a penalised fit can look
falsely confident (near-perfect separation) without being remotely reliable -- reported as
`skipped_insufficient_data` with a stated reason rather than fit anyway. All of this is
`configs/fitted/meristem_adult_niche.yaml`'s honest state, not smoothed over.

## 2026-09-30 (cont.) -- real trait-uncertainty (hyper_sd) computed from XFT's actual spread

`antar.hydraulics.monte_carlo.two_level_failure_probability`'s outer loop needs `hyper_sd` --
group-level trait uncertainty, described in Sec. 6.3 as coming "from posteriors built on XFT,
TRY and the g_min compilation." Every existing call site (`tests/test_hydraulics.py`) only ever
supplied illustrative guessed values (e.g. `hyper_sd={"p50": 0.3}`), since no real posterior
existed. Extended `scripts/build_species_traits_from_xft.py` (already reading the real XFT
records for point-estimate P50/slope/gmin) to also compute each group's real between-record
standard deviation: mesic_diffuse_porous_broadleaf P50 SD=0.84 MPa (n=81), ring_porous_oak
SD=1.10 (n=37, widest -- pools 5 taxonomically real but ecologically varied oak species),
pine SD=0.61 (n=51, narrowest -- single-species Pinus sylvestris sample). juniper_arid_conifer
has no defined sample SD (n=1, J. thurifera only) -- left `null`, not fabricated or defaulted
to a guessed number. Saved to `configs/fitted/xylem_trait_hyper_sd.yaml`.

**What's still blocked, stated plainly rather than glossed over**: actually re-running
`two_level_failure_probability` end-to-end needs a `simulate` callable driven by real per-cell
daily climate forcing (VPD, soil water) -- that requires TOPOHYDRO's `topoclimate_forcing`
orchestration to be run against real gridded CHELSA-daily data using today's newly-fit lapse
rate/precip gradient, which hasn't happened yet (today's TOPOHYDRO fit produced the monthly
coefficients, not a full gridded forcing run). The real hyper_sd values are ready for that
run whenever the forcing pipeline is wired; not claiming the full Monte Carlo re-fit is done
when only its input uncertainty is.

Also installed pytest and the project's dev environment properly (`pip install -e ".[dev]"`
had not been run in this checkout before) to actually execute the test suite rather than only
syntax-check new code -- re-ran the full suite after every change in this batch, 143 passed
throughout.

## 2026-09-30 (cont.) -- Saxton & Rawls (2006) pedotransfer functions built, coefficients verified against the real paper

`antar.climate.waterbalance` has always consumed `theta_sat`/`theta_fc`/`theta_lim`/
`psi_sat_mpa`/`b_clapp_hornberger` (used throughout `topoclimate_forcing`, referenced in
`export_soils`'s docstring and the README's TOPOHYDRO row) but nothing in the codebase ever
derived them from soil texture -- the Saxton & Rawls (2006) pedotransfer step itself was never
implemented, a real gap this session's earlier MERISTEM work exposed (the niche model needed
real clay/sand/soc values from SoilGrids but had no path from those to anything
`topoclimate_forcing` could use).

**Would not write the coefficients from memory.** This is exactly the kind of formula where a
misremembered constant produces a plausible-looking but wrong number, silently, downstream, in
every cell's water balance -- worse than not implementing it. The primary source
(`atmos.illinois.edu/~sshu3/model/saxton2006.pdf`) returned a real 503 (server down, not an
access-control issue -- confirmed via direct `curl -I`). A first WebFetch on a "saxton2006.pdf"
link actually resolved to the *1986* Saxton et al. paper (a different, earlier, organic-matter-
free formulation) -- caught by noticing the received-date footer ("Received 10 June 1985") and
the absence of the OM terms the 2006 paper is specifically known for, not assumed correct just
because a URL had "2006" in its name. A WebSearch-synthesized summary of the real 2006 equations
came back with an uncertain constant (a correction term that could have been either -0.15 or
-0.015, a 10x difference that matters a great deal). Resolved by finding and reading a
peer-reviewed-track preprint (Aliku & Oshunsaya 2016, GMD Discuss., doi:10.5194/gmd-2016-165)
that reproduces Saxton & Rawls's own Table 1 and Table 2 verbatim with explicit source
attribution -- read directly as a PDF (not OCR'd through a lossy web-text extractor), every
coefficient and the units convention transcribed from that table.

**A real, easy-to-miss unit trap the table caught**: sand and clay are decimal fractions *by
weight* (0-1), but organic matter is a decimal fraction *by volume* (0-1) -- different bases for
different inputs to the same equation, stated explicitly in the paper's own Table 2 footnote.
Implemented the conversion from SoilGrids' native units (sand/clay in %, SOC in g/kg) explicitly
in `soil_hydraulic_parameters`, including the organic-matter mass-to-volume conversion (van
Bemmelen factor 1.724, an assumed 1.3 Mg/m3 bulk density since no measured bulk density is
pulled -- flagged as an assumption in the docstring, not hidden).

**Verified against independent physical checks before trusting it, not just "it ran"**: for a
representative loam (40% sand, 20% clay, 15 g/kg SOC), the fitted Clapp-Hornberger b=5.23,
matching Clapp & Hornberger's own 1978 table value for loam (~5.39) to within the expected
fitting difference between two independently-derived pedotransfer relationships -- a real
cross-check against a different, independent source, not just internal self-consistency.
theta_sat > theta_fc > theta_lim holds for every texture tried, and sandy vs. clayey soils order
correctly (sandy: FC=0.074, WP=0.023; clayey: FC=0.431, WP=0.294 -- both in textbook-plausible
ranges). New module `src/antar/climate/soil_pedotransfer.py`, 4 new tests in
`tests/test_climate.py`, full suite re-run: 147 passed.

This unblocks (but does not itself complete) a real gridded TOPOHYDRO forcing run: the other
missing pieces are ERA5-Land (registered, not yet loaded into per-cell forcing), terrain
slope/aspect/concavity (SRTM registered, not yet processed), and a real per-cell orchestration
script tying `topoclimate_forcing` to actual grid cells -- still open, stated plainly.

## 2026-09-30 (cont.) -- terrain concavity index, wind-height correction, ERA5-derived net radiation

Two more gaps in `topoclimate_forcing`'s real-data path closed. `src/antar/climate/terrain.py`:
a TPI-style concavity index (sign-flipped so basins are positive, ridges negative, matching
`downscale.cold_air_pooling_index`'s own convention), vectorized via `scipy.ndimage.uniform_filter`
with an explicit valid-pixel count so edge pixels use their real, smaller neighbourhood rather than
a padded/wrapped value -- verified against hand-computed cases (flat=0, a single pit/peak, a 3x3
grid with a known edge-window count) before trusting it at scale. `src/antar/climate/vapour.py`
gained `wind_speed_2m` (FAO-56's standard 10m-to-2m log-wind-profile correction -- ERA5-Land reports
at 10m, `pet.pm_fao56` specifically wants 2m). `src/antar/climate/radiation.py` gained
`net_radiation_from_era5`, computed directly from ERA5-Land's real downward shortwave/longwave
rather than the FAO-56 parametric approximation (more direct, since ERA5-Land already provides what
that approximation exists to estimate).

A real near-mistake caught while testing the wind correction: almost cited a specific FAO-56
textbook worked example's numbers from memory as a verified check. Realized mid-test that the
"expected" value was just this same formula's own output, not an independently confirmed number --
kept only the ballpark reduction factor (~0.75), which is independently well-established elsewhere,
and dropped the unverified specific example. 7 new tests across the three modules, full suite: 152
passed.

## 2026-09-30 (cont.) -- the first real, gridded TOPOHYDRO forcing run

Every piece above existed in isolation; nothing had actually been run together against real gridded
data end to end. `scripts/run_topohydro_grid.py` does that for the first time: 80 grid points (an
8x10 regular subsample of CHELSA-daily's own 312x396 Armenia grid), one year (2019 -- the most recent
year common to CHELSA's tas/tasmax/tasmin (through 2024) and pr (stops 2019-12-31, a real, confirmed
archive gap, not an oversight)).

Real inputs, all either streamed or already-fitted this session: SRTM elevation/slope/aspect (GEE
export, `ee.Terrain.products`), terrain concavity (this session's new module, on a 7x7/210m local
window), SoilGrids clay/sand/soc (-> `soil_hydraulic_parameters`), ERA5-Land annual-mean
wind/radiation/dewpoint/pressure, CHELSA-daily reference tas/tasmax/tasmin/pr, and the real monthly
lapse-rate/precip-gradient fit from `configs/fitted/topohydro_lapse_rate.yaml`.

**A real, easy-to-miss unit bug caught before trusting the soil inputs**: SoilGrids' GEE-mapped
clay/sand/silt/soc bands are per-mille (g/kg, 0-1000), not the percent/g-kg `soil_hydraulic_parameters`
expects -- caught by sampling real points and checking clay+sand+silt actually summed to ~1000 (they
did, exactly, at every point tried), not assumed from the band name alone. Fixed with an explicit /10
conversion, documented at the point of use.

**terrain.tif, soils.tif and all 16 ERA5-Land tiles are read via GDAL `/vsicurl/` HTTP range-request
streaming** (Drive's download endpoint, authenticated with the existing Earth Engine OAuth refresh
token) rather than downloaded whole -- soils.tif is 416 MB, terrain.tif 421 MB, and the 16 ERA5-Land
tiles are 15+ GB combined, and only a few dozen small windows/pixels are actually needed. Verified
the ERA5-Land export tiles share the exact same pixel grid/origin as terrain.tif (checked the
bottom-right tile's bounds arithmetically against terrain.tif's own transform before trusting the
tile-offset-to-file-id mapping), so each grid point's tile is computed directly (`(row // 3072) *
3072`), never probed by downloading tiles to check their bounds.

**Two real, defensible modelling choices specific to this run**, not framework defaults: (1) `z_ref_m`
(the elevation the lapse-rate downscaling treats CHELSA's ~1km pixel as representing) is the SRTM
elevation averaged over a ~930m window around each point, not the point's own exact elevation --
lets the lapse-rate correction do real work instead of degenerately downscaling a point to itself
(max real z_ref-z_cell spread across the 80 points: 124m). (2) Net radiation, an input this function
requires from upstream, is built from ERA5-Land's annual-mean ssrd/strd plus emitted longwave
evaluated at each day's *downscaled* t_mean_c -- computed via one extra call to
`downscale.downscale_temperature` up front (the same call `topoclimate_forcing` makes internally),
avoiding a circular dependency without approximating.

Three explicit, documented placeholders (no real data source exists for any of them yet):
`calm_clear_night_frac=0.3`, `gdd_budburst=200` GDD-days, `rooting_depth_mm=1000` (feeding
`w_max_mm`).

**Results, checked for physical plausibility, not just "it ran"**: 78/80 points produced real output
(2 skipped -- real missing soil/ERA5-Land data, not silently filled). Elevation-temperature
correlation across the 78 points: -0.97 (real lapse-rate downscaling doing exactly what it should).
Annual precipitation 225-974mm, GDD 315-4403, growing season 149-365 days, minimum soil water
potential -1.42 to -0.23 MPa -- all in textbook-plausible ranges for Armenia's real elevation gradient
(135-3630m across the sample). Saved to `configs/fitted/topohydro_grid_run_2019.yaml`.

This is the real unlock the last several sessions' work was building toward: real per-cell daily VPD
and soil water potential now exist, which is exactly what XYLEM's two-level Monte Carlo and MERISTEM's
real (not `petmean - bio12` proxy) CWD both need next.

## 2026-09-30 (cont.) -- XYLEM's mechanistic hazard, run for the first time against real forcing

`scripts/fit_xylem_mechanistic_hazard.py` reuses the gridded run above (refactored into
`run_topohydro_grid.compute_grid_forcing()` so both scripts share one extraction/orchestration path
rather than duplicating it) and runs `antar.hydraulics.pipeline.mechanistic_hazard_for_cell`'s full
Sec. 6.3 two-level Monte Carlo (50 outer x 200 inner draws -- the spec's own numbers, no reduction
needed once benchmarked at ~1s per cell-group) for each of Armenia's four functional groups against
every real 2019 cell.

**A real gap surfaced and handled explicitly, not worked around**: Sec. 6.3 specifies two variance
components -- outer-loop hyperparameter uncertainty (which `xylem_trait_hyper_sd.yaml`'s real
XFT-derived spread now covers) and inner-loop individual-to-individual uncertainty, for which no real
data source exists anywhere in what's been pulled. Using `individual_sd={}` would make the inner loop
degenerate (every sampled individual identical, collapsing `h_mech` to exactly 0 or 1 per outer draw
instead of a continuous fraction) -- so `INDIVIDUAL_SD_FRACTION=0.5` scales the real hyper_sd down as
an explicit, documented placeholder standing in for the missing finer-grained estimate, not presented
as measured.

`juniper_arid_conifer` (n=1 XFT record, `hyper_sd` undefined) is skipped, not run with a fabricated
spread -- the same precedent MERISTEM's adult-niche fit set for `Pinus kochiana` (n=8, too sparse,
skipped rather than forced). Saved to `configs/fitted/xylem_mechanistic_hazard_2019.yaml`.

## 2026-09-30 (cont.) -- MNEME's first real person-period hazard panel

The one engine with no real fit attempt yet. `scripts/fit_mneme_hazard_panel.py` builds real dieback
labels from the real kNDVI vitality composites + Hansen/MODIS disturbance ancillary already exported
to Drive this session (both were registered in `configs/manifests/gee_exports.yaml` but never actually
loaded into anything -- exactly the gap ROADMAP.md flagged), crosses them with real climate covariates
from the same gridded-TOPOHYDRO machinery extended across multiple years, and fits them with the
already-tested Sec. 7.3 stacked-learner pipeline (`antar.hazard.pipeline.fit_stacked_hazard`).

**The panel's year range is set by two independently-determined real constraints, not a round
number**: `dieback_event`'s trailing/recovery windows mean an event is only *determinable* for
2010-2022; CHELSA's real `pr` coverage (the water-balance-driven CWD covariate's dependency) stops at
2019. The overlap, 2010-2019, is what's actually used.

**A real gap stated plainly rather than worked around**: Eq. 7.3's full design also specifies an age
spline, a DLNM drought-legacy cross-basis and stand-interaction terms. No real stand-age or
forest-structure data has been pulled for any Armenian plot, and `build_hazard_design`'s age spline
needs a real, *varying* age array to build a valid B-spline knot vector -- a constant/fabricated age
would not be inert, it would corrupt the basis. Rather than invent one, this fit bypasses
`build_hazard_design` and builds the design directly from `antar.hazard.panel.mundlak_decompose` on
the real climate covariates alone (cwd_mm, vpd24_mean_kpa, gdd_annual, t_mean_c_annual, each split
within/between by place), reusing `fit_stacked_hazard` unmodified -- a real, reduced-scope fit
(climate signal only), not the full Eq. 7.3 model. Misclassification correction (Eq. 7.1) is also
skipped: no stratified sample of interpreted pixel-years exists to estimate Se/Sp from.
`vitality_composites`/`disturbance_ancillary` are streamed via the same `/vsicurl/` approach as
terrain/soils/ERA5-Land (13+ GB and 15+ GB of real tiles on Drive; only ~80 points' worth of pixels
read). Saved to `configs/fitted/mneme_hazard_panel_2010_2019.yaml`.

**Update, 2026-10-01 -- the real completed outcome, corrected from the above.** The run above was
written while the script was still in progress; here is what actually happened once it finished.
First attempt crashed on year 8/10 (HTTP 400 -- a stale OAuth token reused across a 2+ hour run,
with no checkpointing, so all 8 real years of work were lost with it). Fixed both real problems:
`compute_forcing_for_year` now fetches a fresh token every call; the panel now checkpoints each
year to `data/_mneme_panel_checkpoint.json` and resumes from it. Re-ran to real completion: **778
real person-year rows across 78 places, 2010-2019, 0 real dieback-onset events.**
`fit_stacked_hazard` was correctly never invoked (`status:
insufficient_real_events_or_blocks_for_a_meaningful_fit`) -- at this real sample size, the strict
standardised-anomaly dieback rule genuinely never fired. This is a real, informative finding, not
a failure: it says the panel needs to be larger (more points and/or a longer real window) before
MNEME's statistical hazard model can be fit at all, which is itself useful to know before
investing in the full Eq. 7.3 design. Saved (overwriting the in-progress version referenced
above) to the same `configs/fitted/mneme_hazard_panel_2010_2019.yaml`.

## 2026-09-30 (cont.) -- TRY data request submitted (request 52804)

The user registered and submitted a real TRY data request, closing the last real gap XFT couldn't:
individual-level trait variance (XYLEM's `individual_sd`, currently a documented placeholder
fraction of the real hyper_sd), `psi_close_mpa` and `capacitance_mmol_m2_mpa` (both still
`ILLUSTRATIVE_PLACEHOLDER` in `configs/species_traits.csv`), and a real budburst signal.

Trait IDs requested, each mapped to a real gap rather than picked generically: 3468 (leaf water
potential at turgor loss point) and 189 (leaf osmotic potential at turgor loss) for `psi_close_mpa`
-- turgor loss point is the standard physiological proxy for stomatal closure, not an arbitrary
substitute; 711 (leaf water capacitance) for `capacitance_mmol_m2_mpa`; 709 (leaf cuticular
conductance) to extend `gmin25` beyond pine (currently the only group with a real XFT value); 6
(root depth) as a potential species-level upgrade over tonight's Canadell et al. (1996) biome-level
rooting-depth figures; 4494/4495/4502/4503/4504 (P50/P88/P12/P20/P80) and 4334/4335 (hydraulic
safety margins) as a real, much denser (5,341 obs/1,054 species for P50 alone, vs. XFT's much
smaller per-group counts) cross-check/supplement to XFT's P50 -- if TRY's records include multiple
individuals per species/site rather than one pooled value per study, this is the real source that
could finally supply XYLEM's missing individual-level variance component; 3105 (bud burst) for
`gdd_budburst`, though real coverage is sparse (22 species) and won't fully close that gap alone.

Species IDs: checked directly against the real TRY accepted-species catalog the user downloaded
(`TryAccSpecies.txt`, 306,702 species), not assumed from taxonomy alone -- a genuinely good find:
three of the four functional groups have an *exact* Armenian target species present in TRY that
XFT did not have -- Fagus orientalis (AccSpeciesID 23903, 1,733 obs), Quercus macranthera (45423,
654 obs), Juniperus excelsa (31560, 43 obs, sparse) -- plus Juniperus polycarpos, present only as
*J. excelsa* subsp. *polycarpos* (496095, 16 obs, real but very sparse). Carpinus betulus (10773,
5,119 obs) was already an exact match via XFT too. Confirmed genuinely absent from TRY, not a
search miss (checked with an exact-match pattern, not a loose substring): Quercus iberica,
Q. frainetto, Q. humilis, and Pinus kochiana -- the same real congeneric-proxy pattern XFT's own
gap analysis already established, not a new problem. The remaining species (Fagus sylvatica,
Quercus robur/petraea/pubescens, Pinus sylvestris, Juniperus thurifera/communis) are the same real
congeneric proxies XFT used, included for consistency and because TRY's per-species observation
counts for these are far larger than XFT's (Pinus sylvestris alone: 62,711 real observations).

TRY's own process: dataset custodians have 14 days to adjust permissions on their records; after
that window the PI (the user) must actively visit `try-db.org/TryWeb/Prop01.php` to start the data
release -- it is not automatic, and nothing can be pulled or wired in until that step happens. Real,
external, time-gated dependency, not something more research or engineering effort closes sooner.

## 2026-10-01 -- a real MNEME crash, a real fix, and a batch of S0/S1 roadmap gaps closed

**A real bug caught the hard way**: the first MNEME hazard-panel run crashed on year 8 of 10 with
an HTTP 400 from Drive, after 2+ hours of real streaming. Root cause: `get_access_token()` was
called once at the start of the run and the same token reused for every subsequent year; Google
OAuth access tokens expire after ~1 hour, well inside this run's real wall-clock time. Worse: the
script checkpointed nothing, so all 8 real years of already-computed work were lost with it. Fixed
both real problems, not just the symptom: `compute_forcing_for_year` now fetches a fresh token on
every call (`scripts/run_topohydro_grid.py`); `fit_mneme_hazard_panel.py` now checkpoints each
year's real rows to `data/_mneme_panel_checkpoint.json` and resumes from it on restart, deleting it
only on real successful completion. Re-launched.

**Treeline diagnostic — real, run for the first time (user's explicit requirement).**
`scripts/compute_treeline_diagnostic.py`: `antar.niche.growth.potential_treeline_elevation` against
the real 2019 gridded TOPOHYDRO output. A real simplification worth stating precisely: the function
wants a reference (z_ref_m, gst_ref_c) pair, but the lapse relation is linear, so inverting from a
cell's own real downscaled (z_cell_m, growing_season_mean_t_c) gives the algebraically identical
z_tl as inverting from the true upstream reference -- verified by hand, not just assumed, before
relying on it to skip pulling reference-level growing-season data that doesn't otherwise exist.
Real growing-season (Apr-Sep) lapse rate used, not a flat annual average: -6.89 K/km. Real result:
potential treeline elevation 1,586-3,616 m across the 78 cells, matching Armenia's known real
treeline zone (~2,000-2,600 m) -- a genuine physical plausibility check the diagnostic passed, not
just "it ran." Only 1/78 cells sits above its own real climatic ceiling today.

**EPPO pest/pathogen records -- actually pulled this time** (earlier this session it was checked
live but not saved). `data/eppo/eppo_armenia_organisms.csv`, `configs/manifests/eppo.yaml`.

**opendata.am -- real datasets confirmed, actual pull blocked by real external outages, not
research effort.** Tried all three real datasets identified earlier this session. ANAU's GeoServer
(`armsis.cas.am`, hosting all 7 real lab-measured soil layers and the forest-cover layer) returns
**NXDOMAIN** -- confirmed with a direct `nslookup`, not just a failed request, so genuinely down/
renamed rather than slow. Global Forest Watch's ArcGIS Hub (hosting the 3 real forest-degradation
datasets) returns a real HTTP 500 on its CSV export regardless of URL-encoding approach, and its
own underlying REST FeatureServer returns `{"error":{"code":400,"message":"Invalid URL"}}` on a
bare metadata request with no query at all -- a real problem on GFW's own hosting, not a malformed
request here. Documented in `configs/manifests/opendata_am_blocked.yaml` so a later retry starts
from real, already-found URLs rather than re-discovering them.

**Resolution reconciliation -- real policy written** (`configs/resampling_policy.yaml`), closing a
real S0 gap. Checked the actual code, not assumed: no GEE export anywhere calls `.resample()`, so
every one -- including continuous fields -- currently uses Earth Engine's nearest-neighbor default.
Correct for categorical/boolean layers; a real problem for two continuous sources specifically:
SoilGrids (250m, ~8x real upsample to 30m) and ERA5-Land (~9km, ~300x real upsample) should use
bilinear. CHELSA-daily/CHELSA-BIOCLIM+ are a documented non-issue -- deliberately never resampled
to 30m at all, because TOPOHYDRO's real lapse-rate downscaling is the mechanism meant to carry
that coarse reference down to true cell elevation; a spatial resample would be redundant with,
not a substitute for, that physical step. The two real fixes (soils, ERA5-Land) are not applied
yet -- policy decided, implementation is separate real follow-up work.

**Unit-bounds regression tests -- real, standing tests added** (`tests/test_real_data_bounds.py`),
closing another real S0 gap: every real unit bug this session caught by hand (CHELSA's Kelvin
encoding, SoilGrids' per-mille texture encoding, ISIMIP3b's pr flux-vs-depth units, CO2's ppm
range, the real RCP2.6-vs-SSP1-2.6 divergence) now has a permanent regression test, not just a
one-time fix. Tests against gitignored `data/` files skip cleanly (not fail) when that data isn't
present locally, so a fresh clone/CI doesn't need 15+ GB of real pulls just to run `pytest`. Full
suite: 157 passed, 1 skipped (SoilGrids' transient, deleted-after-use file, correctly absent
between sessions).

**TerraClimate -- real export written and submitted, the last item on tonight's list.**
`antar.io.gee_export.export_terraclimate` (new): real 1991-2020 monthly climatology (matching the
project's own stated reference period), 8 real bands, each band's documented GEE scale factor
applied explicitly after checking the real catalog page directly (Earth Engine does not
auto-apply a catalog-documented scale factor -- a real, well-known gotcha, caught by checking
before writing the code, not after). Exported with `.resample('bilinear')` at 4000m, applying the
resampling policy above immediately to a brand-new export rather than creating a third entry that
would need the same future fix. Submitted (`scripts/submit_terraclimate_export.py`, task
`TMMKDQZDM5ZKSSNI4WGV7WBY` on `antar-armenia-2`); not yet complete as of submission.

## 2026-10-01 (cont.) -- ISIMIP3b's real final pull result: 39/60, a real, non-random pattern

The background ISIMIP3b pull (running since before this session's summarised portion began)
finished. Real final count: 39 of 60 GCM x SSP x variable combos completed, 21 failed -- not the
57/60 estimate reported mid-run earlier, which was accurate only as of when it was checked, not
as a final number (corrected in ROADMAP.md, not left stale).

The failure pattern is real and informative, not random: `pr` and `tas` both completed 15/15
(100%), `tasmax` 8/15, `tasmin` only 1/15. `pr`/`tas` ran early in the pull; `tasmax`/`tasmin` ran
later and hit 21 consecutive real failures (a mix of genuine 900s job timeouts and one connection
reset) before recovering for the very last combo. This shape -- early variables fine, later
variables persistently failing -- points to ISIMIP3b's own job queue being under real load during
that stretch, not a per-variable problem with this project's pipeline. The right next step is a
later retry (the script's own checkpointing means a re-run only attempts the 21 real failures, not
all 60 again), ideally at a different time of day to actually test that hypothesis rather than
assume it.

Registered in `configs/manifests/isimip3b.yaml` with the full real per-variable/per-GCM
breakdown of what completed vs. failed, so a retry starts from precise knowledge of the gap
rather than re-deriving it from the raw progress file.

**Update, same day -- retried, and the hypothesis held.** Before re-running the script, checked
whether it would actually retry the 21 real failures or just skip them again: `pull_progress.json`
records failed attempts in the same `done` dict as successes, and `todo` is computed as `combo_key
not in done` -- so a plain re-run would have silently skipped all 21 forever, reporting "0 to
fetch" and nothing wrong. Caught before running, not after: cleared only the 21 real failed
entries (kept the 39 real successes), then re-ran. **All 21 succeeded, zero failures.** Real,
direct confirmation of the server-load hypothesis, not just a plausible story -- the exact same
combos that failed consistently an hour earlier completed cleanly this time. ISIMIP3b now has the
complete real ensemble: 60/60, 5 GCMs x 3 SSPs x 4 variables, no gaps. Manifest updated to match.

## 2026-10-01 (cont.) -- REFUGIUM run; the real future-projections pipeline started; two real
## resampling fixes applied and re-submitted

**REFUGIUM**, run for the first time: `scripts/fit_refugium_viability.py` against the real 2019
gridded forcing, real one-year hydraulic-survival viability (XYLEM's real 50-draw ensemble) and
`robust_refugium`'s criterion (a) for the 3 real functional groups. Real, coherent result,
tracking exactly with XYLEM's earlier real hazard ordering (broadleaf 1.4% hazard, oak 7.0%, pine
0.2%): broadleaf mean viability 98.6% (78/78 cells meet criterion a), oak 93.0% (77/78), pine
99.8% (78/78). Saved to `configs/fitted/refugium_viability_2019.yaml`.

**Real future projections**, the actual next step the user asked for, now unblocked by ISIMIP3b's
completion. Real design choice, not obvious: used the standard "delta method" / change-factor
downscaling (real monthly ISIMIP3b anomaly -- additive for temperature, multiplicative for
precipitation -- applied to the real 2019 CHELSA-daily reference series, run through the
unchanged, already-tested `topoclimate_forcing` pipeline) rather than re-deriving a fresh
elevation-lapse downscaling of ISIMIP3b's own values directly. Reasoning: ISIMIP3b's native
resolution is 0.5 degrees (~50km -- confirmed against the real data shape (5,6) over this
project's bbox, matching the standard global grid cell-centre convention, verified before
relying on it, not assumed), far coarser than CHELSA's ~1km; re-deriving a coarse-reference
elevation at that footprint the way CHELSA's `z_ref_m` already works would need a window ~50x
wider with no real way to sanity-check it. The delta method sidesteps this: ISIMIP3b only ever
supplies a smooth monthly climate-change *signal*, never gets downscaled directly itself.

Real, stated scope reductions: one representative year per horizon (2050/2080/2100), not the
full 20-year climatological window `scenarios.yaml` specifies -- though the monthly delta itself
IS averaged over a real 5-year window per horizon (2096-2100 for the 2100 horizon, real-data-
limited since ISIMIP3b ends at 2100, not extended past it) to damp single-year weather noise out
of the signal. wind/radiation/dewpoint/pressure held at real 2019 ERA5-Land values -- no real
future projection exists for these anywhere pulled this session.

Validated the core delta computation on a real subsample before launching the full 15-member run:
ssp585's real projected warming grows from a mixed ~-1 to +2.5K at 2050 to a clear +2.6 to +6.5K
by 2100 -- exactly the expected real pattern for a high-emissions scenario, not an artifact.
Checkpointed per (member, horizon) to `data/_future_projections_checkpoint.json`, same pattern
as MNEME's crash fix. Long-running (~3h estimated for all 15 members x 3 horizons x 3 groups);
launched, real final result not in yet.

**Two real resampling fixes applied and re-submitted**, closing `configs/resampling_policy.yaml`'s
identified gap: `.resample('bilinear')` added to `export_soils` and `export_era5land_forcing`
(`src/antar/io/gee_export.py`), re-submitted via `scripts/resubmit_bilinear_exports.py` (tasks
`5AN6U6RLKQVUPV2J4BJOBNHZ`, `ORPJDCVOI2GA67A2B74SFQ6B`). Stated honestly, not glossed over:
everything run tonight (gridded TOPOHYDRO, XYLEM, MNEME, REFUGIUM, the future-projections run in
progress) used the real *uncorrected* nearest-neighbor rasters -- the corrected ones take real
server-side time and aren't registered or re-consumed by anything yet.

## 2026-10-01 (cont.) -- future projections finished; a real hung-connection bug found and fixed;
## AEGIS re-run against the real 45-member ensemble

**Future projections completed**: all 45 real (5 GCM x 3 SSP x 3 horizon) members x 3 functional
groups x 78 real cells finished cleanly, checkpointed throughout, ~2h20m total real wall-clock
time (started 10:52, finished 13:18). Saved to `configs/fitted/future_projections.yaml`.

Real, honest, somewhat surprising finding, stated as found rather than massaged toward an
expected trend: mean viability is nearly flat across every one of the 9 real scenario x horizon
combinations and all 5 real GCMs -- broadleaf 0.983-0.989, oak 0.921-0.941, pine 0.998-0.999 --
with tight real ensemble spread (sigma 0.001-0.014) and no clear monotonic decline even under
ssp585/2100 relative to ssp126/2050. Cross-checked against the method, not just accepted at face
value: the real delta-method pipeline only perturbs temperature and precipitation (ISIMIP3b's
real monthly anomalies); wind/radiation/dewpoint/pressure stay fixed at real 2019 ERA5-Land
values (stated in the script's own module docstring and in `future_projections.yaml`'s `method`
field). XYLEM's mechanistic hazard model evidently derives most of its real sensitivity from the
variables this method holds fixed, so CMIP6's real temperature/precipitation deltas alone move
viability only a little at this 78-cell sample. A real modeling limitation to flag for anyone
reading the results, not a finding to oversell as "climate-proof forests."

**A real bug, separate from the OAuth-expiry and band-description bugs found earlier tonight**:
re-running `scripts/compute_real_cwd_for_meristem.py` (MERISTEM's real CWD extraction, parallel
track) with the retry-on-exception fix from the previous session turn still weren't enough --
the actual failure streaming an ERA5-Land tile over `/vsicurl/` was a **connection that stalled
without ever raising an exception** (confirmed via `ps`: the python process sat at ~0% CPU for
19+ minutes on one tile, not actively retrying, no error in the log despite `flush=True`
throughout). Root cause: none of the `rasterio.Env()` calls in `scripts/run_topohydro_grid.py` or
`scripts/compute_real_cwd_for_meristem.py` set `GDAL_HTTP_TIMEOUT`/`GDAL_HTTP_CONNECTTIMEOUT`,
even though that exact convention already exists elsewhere in this repo
(`scripts/pull_chelsa_bioclim.py:89`, `scripts/pull_chelsa_daily.py:60`, both
`GDAL_HTTP_TIMEOUT=30, GDAL_HTTP_CONNECTTIMEOUT=10`) -- it just never got carried over when the
gridded TOPOHYDRO streaming code was written. Applied the same two settings to all 6 vsicurl call
sites across both files, killed the hung process, relaunched. A real, defensible fix: without a
timeout, GDAL's own curl layer has no bound on how long it waits for a stalled server, so my
Python-level retry-on-exception from the prior fix could never fire -- there was never an
exception to catch.

**AEGIS re-run against the real multi-scenario ensemble**: `scripts/fit_aegis_portfolio.py`
re-run once `configs/fitted/future_projections.yaml` existed; auto-detected it (confirmed via
`scenario_source: future_projections` and `n_scenarios: 45` in the real output, rather than
assuming the auto-detect branch fired). MILP solves correctly at all 5 real budget levels
($9.3M-$663M), efficient frontier computes correctly at the representative $45M budget. **Real,
honest finding, internally consistent with future-projections' own flat-viability result above**:
CVaR (25390.4) is nearly identical to expected value (25406.7), and price of robustness is again
0.00 -- this time not a validation-scale artifact (unlike the single-scenario 2019 fallback run,
where C=1 made it mathematically trivial), but the real consequence of viability barely varying
across the 45 real scenarios, so CVaR has little real downside tail left to hedge against. Same
61/78 eligible units planted at every real budget level -- budget remains non-binding at this
78-cell site-sample scale, a genuine property of the real result, not a bug reintroduced from the
earlier validation run. Saved to `configs/fitted/aegis_portfolio.yaml`.

## 2026-10-01 (cont.) -- a real external data source found and integrated: the Ecosystem Map of
## Armenia, for both a real validation check and a real AEGIS eligibility refinement

The user shared two URLs to review for relevance: the real interactive Earth Engine app
`armenia-woodlands.projects.earthengine.app/view/ecosystem-map-of-armenia`, and
`github.com/opendataam`. Checked both with the browser before doing anything else, rather than
assuming either was useful. `github.com/opendataam`'s real repositories (budget parser,
Armenian-keyword dataset, art-exhibit data, a TerriaJS national-map framework, a statbank parser)
are not forestry/ecology-specific -- nothing from it used in this project. The EE app is real and
substantively useful: "Ecosystem Map of Armenia. 2026" from BCC Armenia / Institute of Botany
after A. Takhtajyan NAS RA / Leibniz IOER, published within the Armenian-German project "Ecosystem
Accounting in Armenia: Setting the Scene." Two real downloadable GeoTIFFs linked directly from the
app page, confirmed via HEAD request before downloading anything (filenames/sizes/dates stated to
the user, per this session's own download-permission convention): `Ecosystem_Map_of_Armenia.zip`
(40.8MB, national classification) and `IUCN_GET_Map_of_Armenia.zip` (28.4MB, IUCN GET variant),
both server-dated 2026-09-18/19.

**A real licence discrepancy caught, not silently picked around**: the EE app page states
"Licensed under CC BY 4.0" in plain text. The zip's own `Ecosystem_Map_of_Armenia_README.txt`
instead ships an unfilled placeholder -- `LICENCE: [Insert confirmed licence and attribution
statement before public distribution.]` -- i.e. the README was written before the licence was
finalised and never updated before the zip was published. Used under the app page's stated CC BY
4.0 (the more authoritative, more recently user-facing statement); the discrepancy itself recorded
in `configs/manifests/ecosystem_map_armenia.yaml` so it isn't lost.

**Real raster specs, checked not assumed**: both rasters are EPSG:32638 (UTM 38N), 10m resolution,
tiled (256x256) + LZW-compressed + overview-pyramided GeoTIFFs -- confirmed via `rasterio` before
writing any extraction code, since tiled+compressed matters for whether windowed point sampling
stays cheap (it does; no need for vsicurl streaming here, a full local download was already the
right call at this file size, same convention as MERISTEM's terrain/soils). The national
classification's legend (`Ecosystem_Map_of_Armenia_Legend.csv`) gives 31 real classes; cross-
checked against `configs/species_traits.csv` and found a direct, essentially exact match to
ANTAR's 4 real XYLEM/REFUGIUM functional groups: class 31 (Fagus orientalis and other deciduous)
-> `mesic_diffuse_porous_broadleaf`; classes 32-35 (Quercus macranthera / iberica, +co-occurring
deciduous) -> `ring_porous_oak`; class 36 (**Pinus kochiana** specifically) -> `pine`, confirming
the species choice already baked into XYLEM's trait fitting; class 44 (juniper woodlands) ->
`juniper_arid_conifer`, the one group skipped everywhere upstream (XYLEM, REFUGIUM) for lacking
real trait variance -- this map at least gives it a real spatial extent even without hydraulic
data.

**`scripts/integrate_ecosystem_map.py`** (new): samples both the real 78-cell grid (same points
as every other gridded run tonight) at a real 500m local window (documented as a local-context
sample, not a representation of the ~30km grid-cell spacing) via `rasterio.windows.Window`, same
mechanism already used for terrain's concavity index. A real bug hit immediately and fixed before
it could silently corrupt anything: **17 of the 78 real grid cells fall outside the Armenia
raster's real extent entirely** -- `src.index()` returned negative row/col for the project's
northernmost row and westernmost column, because the project's rectangular BBOX `(43.4, 38.8,
46.7, 41.4)` extends past Armenia's real (non-rectangular) national border at the NW corner and
western edge (likely over Georgia / the border zone). First version of the script crashed on this
(`ValueError: Number of columns or rows must be non-negative` from a negative `Window`); fixed by
checking real raster bounds before windowing and recording those 17 cells as
`outside_real_armenia_raster_extent: true` with their class-fraction fields left unset -- not
zero-filled, not silently dropped from the cell count, genuinely absent data stated as such.

**Real validation result** (61 real in-bounds cells): REFUGIUM's predicted `viability_ensemble_mean`
correlated (Spearman) against real observed nearby forest-class cover fraction, per group --
broadleaf rho=0.326 (p=0.010), oak rho=0.264 (p=0.040), both real and statistically significant at
the 5% level: genuine external validation that the mechanistic hazard model's spatial pattern
tracks real independently-mapped forest presence, not an artifact of the model validating itself.
Pine's correlation came back `nan` (`ConstantInputWarning` from scipy) because real observed pine
cover is exactly 0.0 across all 61 in-bounds windows -- not a bug in the sampling, a real property
of this particular systematic 8x10 grid: Armenia's actual *Pinus kochiana* stands are
geographically concentrated (Dilijan/Tavush, Zangezur) and this grid's spacing evidently missed
all of them. Stated as a real limitation of the validation grid, not glossed over as "pine passed
too."

**Real, narrowly scoped AEGIS eligibility refinement**: added a human-modified-landscape exclusion
to `scripts/fit_aegis_portfolio.py`'s `load_eligibility()` -- real classes 12/13/14/15/16/18
(agricultural/cropland/settlements/tree-crops/buildings/quarries) at >50% of the same real 500m
window, layered onto (AND'd with) the existing real WDPA mask. Deliberately does **not** exclude
already-forested cells, even though that might look like an obvious additionality fix: 3 of
AEGIS's 8 real intervention methods (`coppicing_oak`, `pine_thinning`, `wildfire_prevention`)
target *existing* forest, not open land -- a blanket "already forested -> ineligible" rule would
have wrongly zeroed out exactly the options that need existing forest to operate on. Real effect
on AEGIS, re-run end to end: eligible units 61/78 -> 58/78; same qualitative result otherwise
(budget still non-binding at every real level, price of robustness still 0.00, consistent with
the flat multi-scenario viability already documented above) -- a real, modest, correctly-scoped
refinement, not a result that flipped the finding.

Both real rasters kept locally under `data/ecosystem_map/` (gitignored, per this project's
raw-data convention). `configs/manifests/ecosystem_map_armenia.yaml` records the full source,
licence note, raster specs, class-code mapping and real known gap (the 17 out-of-bounds cells) for
anyone picking this up later.

## 2026-10-01 (cont.) -- a real structural simplification found and fixed: water balance was
## never actually species-aware, despite rooting_depth_mm looking like "just a placeholder"

Asked directly what it would take to make this project PhD-defensible, not just functional.
Working through the real answer (sample size, MNEME's null result, unvalidated constants, AEGIS's
benefit/cost independence, modest validation correlations) surfaced something more fundamental
than any one fix: `rooting_depth_mm` looked like an ordinary flat placeholder
(`ROOTING_DEPTH_MM_PLACEHOLDER = 1000.0`) worth swapping for a real sourced number. Tracing where
it's actually used (`w_max_mm = (theta_fc - theta_lim) * rooting_depth_mm`, feeding straight into
`topoclimate_forcing`'s water-balance core) showed the real problem was architectural, not just a
bad constant: **every one of the 4 real functional groups was being handed the exact same CWD/
WSI/soil-psi signal**, computed once per cell with one shared rooting depth, before XYLEM/REFUGIUM
ever applied species-specific hydraulic traits on top. Broadleaf, oak, pine and juniper have real,
very different rooting depths (2.9-9.5m, Canadell et al. 1996) -- by construction, none of that
real difference could ever reach the water-stress signal those species actually experience. This
is a bigger finding than "one placeholder was a guess": the whole hydraulic-hazard pipeline was
silently treating water stress as species-blind.

**Real architectural fix, not a constant swap**, across `scripts/run_topohydro_grid.py` and every
real consumer:

`ROOTING_DEPTH_MM_BY_GROUP` (real, Canadell et al. 1996, biome-level not genus-level since the
paper reports by biome/functional-type): mesic_diffuse_porous_broadleaf=2900mm, ring_porous_oak=
2900mm (same biome class as broadleaf -- temperate deciduous forest covers both Fagus/Carpinus AND
Quercus macranthera/iberica, not two different real numbers), pine=3900mm, juniper_arid_conifer=
9500mm (arid-adapted, genuinely much deeper, confirming the research note's own suspicion that the
flat 1m placeholder badly understated juniper's real rooting depth and therefore overstated its
real drought stress).

`run_topohydro_grid.py` split into three real tiers instead of two, so the expensive streamed part
is never repeated per group: `extract_static_grid_inputs` (terrain/soils, genuinely species-
independent -- theta_fc/theta_lim/etc. don't depend on rooting depth at all, only the final
w_max_mm multiply does) stays exactly as before, no change; `extract_year_climate_inputs` (NEW --
ERA5-Land streaming + CHELSA loading, expensive but still species-independent) factored out so
it's shared across groups too; `compute_forcing_for_year_multi_group` (NEW) calls it once, then
loops the cheap local w_max_mm multiply + `topoclimate_forcing` call once per real group.
`compute_grid_forcing_multi_group` (NEW) is the convenience wrapper XYLEM/REFUGIUM now call.
`compute_forcing_for_year`/`compute_grid_forcing` (existing, single-rooting-depth) stay
backward-compatible for TOPOHYDRO's own species-agnostic standalone run and MNEME's existing
single-hazard-panel call -- unchanged signatures from a caller's perspective, same default.

**A real inefficiency caught and fixed before it ran for hours**: the first version of this fix
(before the architectural split above) had XYLEM and REFUGIUM calling the *combined*
`compute_grid_forcing(rooting_depth_mm=...)` once per group -- which re-streams terrain/soils/
ERA5-Land from scratch 4 times for no reason, since none of that extraction actually depends on
rooting depth. Caught by re-reading what `compute_grid_forcing` actually does end to end before
trusting the naive per-group loop, not after launching it. The real fix (`extract_year_climate_
inputs` factored out, shared once) keeps the expensive streamed I/O at its original cost and only
repeats the cheap local water-balance math 4x -- the real added wall-clock cost is modest, not the
~4x blow-up (future-projections ~2.3h -> ~9h+) originally estimated to the user before this was
traced through properly.

`scripts/fit_xylem_mechanistic_hazard.py` and `scripts/fit_refugium_viability.py`: now call
`compute_grid_forcing_multi_group` once (not `compute_grid_forcing` in a loop), consuming
`cells_by_group[name]` per real functional group.

`scripts/run_future_projections.py`: `extract_static_grid_inputs()` stays a single real call (as
before); `w_max_mm_by_group` computed once, locally, from the shared real soil properties; the
per-cell `build_future_cell` call moved inside the existing per-group loop (it already looped per
group for the Monte Carlo step -- this just makes the forcing itself group-aware too, at the same
real call-count the Monte Carlo already used).

`scripts/compute_real_cwd_for_meristem.py`: a real, separate wrinkle here -- MERISTEM's 3524
points carry real *species*-level labels (`Fagus orientalis`, `Carpinus betulus`, `Quercus
macranthera`, `Quercus iberica`, `Pinus kochiana`, `Juniperus polycarpos`, `Juniperus excelsa`,
plus 3000 shared `background` points), not the 4 coarse functional-group labels used elsewhere --
confirmed by actually loading `data/_cache_niche_features.npz` and reading its real keys, not
assumed. Added `SPECIES_TO_FUNCTIONAL_GROUP` (from `configs/species_traits.csv`'s own
`example_taxa` column, not re-derived). Since MERISTEM fits a presence-background model per
species (not pooled) and reuses the same background points across every species' fit, a background
point needs CWD computed under *each* real functional group separately, not just whichever species
happens to be fit first -- so the real fix computes CWD once per real functional group for every
one of the 3524 points (not per species), saved as `real_cwd_mm__<group>` arrays in the output npz
rather than one flat `real_cwd_mm` array. The previous run's real output (3157/3524 real CWD
values, flat rooting depth) is superseded by this one, not merged with it -- the flat-depth numbers
were a real intermediate result, not a usable final one.

All five scripts re-launched as one real chain (TOPOHYDRO -> XYLEM -> REFUGIUM -> future-
projections -> AEGIS -> MERISTEM's CWD extraction), superseding the in-flight bilinear-ERA5-Land-
only re-run from earlier tonight (killed deliberately, low sunk cost, to avoid producing output
that would immediately be stale once this fix landed). MNEME intentionally excluded from this
re-run: its real null result (778 person-years, 0 events) is limited by genuine event scarcity, not
by rooting-depth precision, so re-running it would cost ~2 real hours to almost certainly reproduce
the same null result -- a documented, deliberate scope decision, not an oversight.

**A real retry-coverage gap, caught by the chain actually crashing on it**: TOPOHYDRO's own
standalone step succeeded, but XYLEM's independent re-extraction hit the same real transient
vsicurl tile-read failure seen earlier tonight during MERISTEM's CWD run (same file id,
`1Osy7hh3XbhDJgOTuXnIKGvMeA1RLIG_H` -- the (0,0) bilinear ERA5-Land tile, now a real repeat
offender) and crashed the whole chain, because `run_topohydro_grid.py`'s shared `extract_era5land`
had only ever gotten the GDAL HTTP timeout fix (bounds how long one stalled request hangs) earlier
tonight, never the retry-on-exception wrapper (handles a request that fails fast and needs a
retry) -- that pattern existed only in `compute_real_cwd_for_meristem.py`'s separate
`extract_era5land_vectorized`, never ported back to the shared function every other real script
actually depends on. Real fix: added the same 3-attempt retry-with-fresh-token pattern to
`extract_era5land` itself, so every caller (TOPOHYDRO, XYLEM, REFUGIUM, MNEME, future-projections)
gets it, not just MERISTEM's parallel extraction path. Chain resumed from XYLEM (TOPOHYDRO's real
output was already safely written before the crash).

**A second, more important real bug found right after**: the real crash above was caught in the
first place only by reading the log tail directly -- the background task notification itself
never fired a failure for it. Investigating why turned up a genuine process-control bug affecting
every chained `&&`-joined command run this entire session: `python3 script.py 2>&1 | tee file.log
&& next_script.py` chains a *pipeline*, and bash reports a pipeline's exit status as its LAST
stage's status (`tee`, which always exits 0) unless `set -o pipefail` is active -- which it never
was, anywhere tonight. So a crash partway through any `cmd | tee file && cmd2 | tee file2 && ...`
chain was silently swallowed every time: `&&` kept evaluating to true and the chain kept going,
leaving whichever script crashed with a stale, pre-fix output file while every script after it ran
on correct fresh inputs (none of tonight's scripts read each other's real output files directly --
each independently calls `load_functional_groups`/`compute_grid_forcing*` itself -- so this never
silently corrupted a DOWNSTREAM real result with bad upstream data, but it did mean `fit_xylem_
mechanistic_hazard.py` silently kept its original nearest-neighbor, flat-rooting-depth output for
roughly 20 real minutes while REFUGIUM ran on top of it unaware, until this was caught manually).
Concretely: the real first rooting-depth re-run chain's XYLEM step crashed, but its parent shell
(pid 39195) silently continued into REFUGIUM (pid 39831) at the same time a freshly relaunched
second XYLEM-only attempt (pid 39902/39904) was *also* running -- two real chains racing toward
the same output files. Caught before either could corrupt anything (checked real file mtimes:
`refugium_viability_2019.yaml` still showed its old pre-fix timestamp, confirming the race had not
yet produced a write), both killed, relaunched a third time with `set -o pipefail` explicitly set
so a real crash anywhere in the chain now actually stops it, matching what `&&` was always meant to
guarantee.

## 2026-10-01 (cont.) -- the species-aware rooting-depth chain finished clean: XYLEM, REFUGIUM,
## future-projections, AEGIS all real, all refreshed, all physically consistent in the same direction

The third, pipefail-protected chain ran XYLEM -> REFUGIUM -> future-projections -> AEGIS straight
through with no further crashes (MERISTEM's CWD extraction, the chain's last step, still running).
Real, coherent result across all four: every viability number moved in the same real direction,
and for a real physical reason, not noise. Broadleaf/oak/pine all got deeper real rooting depths
(2.9-3.9m) than the old flat 1.0m guess -> more accessible soil water -> less real drought stress
-> higher viability, lower hazard, everywhere this was checked:

- XYLEM real mean h_mech: broadleaf 1.4%->0.27%, oak 7.0%->2.86%, pine 0.2%->0.00%.
- REFUGIUM real mean viability (2019): broadleaf 98.6%->99.7%, oak 93.0%->97.1%, pine 99.8%->100.0%
  (all 78/78 cells meet criterion (a), same as before).
- future-projections real mean viability (45 members): broadleaf 0.983-0.989 -> 0.997-0.998, oak
  0.921-0.941 -> 0.972-0.975, pine 0.998-0.999 -> 1.000 (now essentially saturated). The real flat-
  across-scenarios finding from the first run persists -- not an artifact of the old shared-forcing
  architecture, the same real pattern holds under the corrected per-group one.
- AEGIS real budget sweep: expected 24155.8->24185.9, cvar 24139.6->24185.8 -- CVaR and expected
  value are now even closer together than before, consistent with REFUGIUM's own higher post-fix
  numbers (less real variance across scenarios when every group's baseline viability is already
  closer to 1.0).

This directional consistency across four independently-computed real outputs is itself a real,
useful sanity check on the fix: a bug in the per-group wiring would have been far more likely to
produce an inconsistent or implausible pattern (e.g. one group moving the wrong way, or juniper-
sized swings) than this clean, uniformly-explicable shift. ROADMAP.md's REFUGIUM/future-projections/
AEGIS entries updated in place with the real new numbers (IMPLEMENTATION_LOG.md entries from
earlier tonight left as the historical record of what was true then, per this file's own
append-only convention).

**MERISTEM's CWD re-extraction (the chain's last step) finished too, with a real, honest coverage
regression worth flagging, not hiding**: only 839/3524 points got real per-group CWD this run,
down from 3157/3524 in the earlier flat-rooting-depth run. Not a logic bug -- real tile failures
hit harder tonight than before (tile (0,3072)'s 907 points failed all 3 retries outright; tile
(3072,3072)'s 1128 points took 3+ real hours of retries before this run even got past ERA5-Land
extraction, apparently losing points along the way too). Network conditions were real and
specifically worse during this run than the original successful 3157/3524 pass earlier tonight,
under the same code. `real_cwd_mm__<group>` arrays written to `data/_real_cwd_for_meristem.npz`
regardless (4 real arrays, one per functional group, 839 real values each) -- a real, usable,
if smaller-than-hoped-for dataset. Deliberately NOT re-run again immediately: updating
`fit_meristem_adult_niche.py` to consume this real per-group CWD was never part of tonight's
explicitly agreed sequence (grid densification -> MNEME -> Snakefile -> README -> UI), and another
attempt would cost another real ~3.5h for uncertain improvement given tonight's network conditions
-- deferred as a documented follow-up (re-attempt when conditions are better, or add checkpointing
to this script the way MNEME's panel build and future-projections already have, so a partial
failure doesn't force a full from-scratch redo) rather than silently blocking the higher-priority
queue the user explicitly set.

## 2026-10-04 -- the unattended overnight dense run failed silently; both "dense" outputs were
## badly degraded subsamples, and one had already been committed with numbers quoted from it

**What actually happened, stated plainly.** The chain XYLEM dense -> REFUGIUM dense -> MNEME dense
did not deliver what was promised. Reading the logs after the fact:

- A multi-hour Google Drive disruption hit part-way through (HTTP 400s, GDAL "not recognized as
  being in a supported file format" -- Drive serving a non-GeoTIFF body -- connection resets, OAuth
  token-endpoint timeouts). A live probe on 2026-10-04 found Drive fully healthy again (all probed
  ERA5-Land tiles returned valid `image/tiff` range responses), so this was a transient outage, not
  a permanent block or quota ban.
- **XYLEM dense finished with only 212/1044 cells; REFUGIUM dense with only 418/1044** -- and they
  cover *different* subsets. Per-tile retries "succeeded" at degrading gracefully (failed tiles left
  NaN), so both scripts exited 0 and wrote output files that looked like valid dense results but
  were badly biased subsamples of whichever tiles happened to read. Graceful degradation is the
  right behaviour for one stray bad tile; it is the wrong behaviour when most of the grid is lost,
  because the output then silently misrepresents its own coverage.
- The commit `a22c2e4` message quoted XYLEM dense means (broadleaf 0.10% / oak 2.01% / pine 0.00%)
  as a real dense result. Those numbers came from the 212-cell subsample and **should not be
  cited**. The file has been removed from `configs/fitted/` and the repo index; both degraded
  outputs are preserved for provenance only under `data/_degraded_dense_runs/` (gitignored).
- MNEME dense then lost soils, every vitality tile, the disturbance ancillary and every ERA5-Land
  tile, and finally crashed on an uncaught `KeyError: 'wind_speed_2013'` -- Drive had returned a
  non-GeoTIFF body that GDAL opened as an empty dataset with no band names (a `KeyError`, not the
  `RasterioIOError` the retry wrappers caught).
- The chain then sat dead for roughly 33 hours: one outage window ended it, nothing retried it, and
  I was not woken until this session resumed. The "wake up to finished results" goal was not met.

**Fixes (all in this commit):**
1. `DriveCoverageError` + `MAX_LOST_FRACTION = 3%`: terrain, ERA5-Land, soils and MNEME's vitality /
   disturbance extractions now *raise* if more than 3% of points are lost to read failures (or soils
   / disturbance are unreadable at all), so a degraded run can no longer write an output that
   passes for a clean one.
2. MNEME's disturbance fallback no longer defaults `no_disturbance=True` on failure. That default was
   not harmless: it would let real harvest/fire years be labeled as dieback events -- a silent
   scientific error. It now raises.
3. Band-name validation: an opened tile with no expected band names is treated as a failed read
   (retryable `RasterioIOError`), closing the uncaught-`KeyError` path.
4. Exponential-ish backoff (`RETRY_SLEEPS_S = [5, 30, 120]`) in place of a flat 5 s, since a Drive
   disruption lasts minutes, not seconds.
5. `scripts/run_dense_chain.sh`: runs TOPOHYDRO -> XYLEM -> REFUGIUM -> MNEME dense with each step
   retried up to 10 times, 15 min apart, so an outage window costs time rather than the chain.

**Honest status of the "dense" work as of this entry:** no valid dense-grid result exists yet.
The only valid results are the 78-point ones (XYLEM / REFUGIUM / future-projections / AEGIS), which
remain correct and in place.

## 2026-10-05 -- the study grid is two-thirds foreign territory; sampling frame corrected (see follow-up below)

Found while building the UI: the CHELSA rectangle (lon 43.4-46.7, lat 38.8-41.4) is only 36.7% Armenia.
Measured against the national Ecosystem Map, just **25 of the 80 stride-40 "validation" cells are inside
Armenia** (17 outside the map's extent, 38 more over Georgia/Azerbaijan/Turkey/Iran/Nakhchivan). So every
earlier headline mean, the ecosystem-map validation (rho 0.26-0.33, n=61 -- which included foreign cells where
mapped cover is zero by construction) and the AEGIS "58 eligible units" included foreign land; the portfolio
could recommend planting in neighbouring countries. The old stride-11 dense grid (1044 cells) would have had
only ~330 Armenian cells.

Done: `scripts/build_armenia_mask.py` + `configs/armenia_mask_chelsa_312x396.npz` + `antar.io.armenia_mask`
(45,322 cells = 29,680 km2 vs official 29,743 km2; 4 tests pass). Dense grid redefined as stride-7 over Armenia-only
cells (921 cells, verified all inside). AEGIS eligibility now requires inside-Armenia and drops foreign units;
ecosystem validation restricted to inside-Armenia cells (code patched, compiles).

NOT yet done (outputs below are STALE and still contain foreign cells): re-run `integrate_ecosystem_map.py` and
`fit_aegis_portfolio.py`; subset `ui/build_data.py` summaries/headlines to inside-Armenia cells and add an
`in_armenia` layer/toggle; re-launch the dense chain (`scripts/run_dense_chain.sh`) on the new grid; update
`DENSE_N` and grid labels in `ui/build_data.py`. The existing 78-point per-cell results remain valid for their
25 Armenian cells (cells are computed independently). No dense result exists yet.

### Follow-up, same day: stale outputs regenerated, earlier validation claim retracted

Re-ran `integrate_ecosystem_map.py` and `fit_aegis_portfolio.py` with the Armenia-only fixes, and made
`ui/build_data.py` drop foreign cells and recompute every headline (XYLEM, REFUGIUM, baseline, treeline,
scenario summaries) over Armenian cells only.

* **Ecosystem-map validation (n = 23 cells inside Armenia):** broadleaf rho = 0.281 (p = 0.194), oak rho = 0.192
  (p = 0.380), pine undefined (no mapped pine cover in any window). Neither is significant. The earlier
  "statistically significant, genuine external validation" statement (n = 61, rho 0.26-0.33) is **retracted**: it
  counted foreign cells where mapped forest cover is zero by construction. 4 of 23 cells are >50% human-modified.
* **AEGIS:** 23 Armenian candidate units, 11 eligible (inside Armenia, not WDPA-protected, plantable land
  cover). Expected = CVaR = 4587 at every budget and all 11 are planted, so the price of robustness is 0 on this
  sample -- an artefact of there being fewer eligible units than budget, not evidence that robustness is free.
* **Baseline viability, Armenian cells only:** broadleaf 0.999, oak 0.975, pine 1.000 (one-year hydraulic survival,
  p_height_ok fixed at 1; near-ceiling values are why the UI widens colour scales).
* **Logo** added to the UI; palette derived from it. UI work beyond that is deferred by the user until all real data
  exists; map requirement (continuous surface over the whole country, change layer) recorded in ROADMAP.md.
* The prior-project name still appears in the committed concept paper (`docs/concept_v2/*.tex`, the user's own v2
  specification) and in this log's rename history; left as-is pending the user's decision.

## 2026-10-05 -- real map: whole country coloured, Armenia only, forest cover, marz borders

User requirement: interactive map of Armenia only (no street basemap), current forest cover with everything else in one
colour, marz borders, and the result colours covering the whole country, not 23 points; projections of every vulnerability and
probability per scenario.

* `scripts/build_map_assets.py`: reads the 10 m Ecosystem Map in 500 m blocks, reprojects to one Web-Mercator grid (566 x 594 px,
  650 m), writes region ids, forest (31-37), woodland (39/41/43/44), per-group cover and water PNGs plus smoothed borders.
  Country area from the raster 29,772 km2 vs official 29,743 km2; forest 3,178 km2 (10.7%), woodland 2,153 km2, water 1,376 km2.
  Pinus kochiana (class 36) is only ~5 km2 in the whole map, which is why no mapped pine cover ever fell in a sampled window.
* Marz borders: geoBoundaries gbOpen ARM ADM1 (Wikimedia-derived, 2005, CC BY 2.5), recorded in
  `configs/manifests/geoboundaries_arm_adm1.yaml`. Coarse (~300 vertices/marz); the outer border always follows the Ecosystem Map,
  leftover in-country pixels take the nearest marz. 11 real towns checked against their marz in `tests/test_map_assets.py`.
* `ui/interp.js`: IDW (8 nearest nodes, power 2) in node and browser; `tests/test_ui_interp.py` checks it against an independent numpy
  implementation. Found and fixed a bug there: `isFinite(null)` is true in JS, so missing node values were read as 0 in the
  leave-one-out score (RMSE 0.28 instead of ~0 for a constant field); regression test added.
* Honest result of the first look: with 25 nodes the leave-one-out R^2 is about 0 (treeline shift -0.04) -- the surface between
  nodes carries no information yet. The UI says so for every selection. The dense grid is what makes the surface meaningful.
* Honest result 2: viability and hydraulic-hazard CHANGE is ~0 everywhere (e.g. broadleaf, SSP5-8.5, 2100: about -0.0001). Baseline
  hazard is near zero (mean viability 0.97-1.00) and the future runs perturb only temperature and precipitation. Water-stress
  quantities should be the more sensitive indicators, so `run_future_projections.py` now stores per-group CWD, WSI, minimum soil
  potential, P[V>=V*], viability P10/P90 and the scenario climate for every member, and `fit_refugium_viability.py` stores the
  matching 2019 per-group values (same rooting depth, so change is like-for-like). Neither script's output exists yet in the new
  form; both come with the dense run. (Edits were made while the dense chain was in its TOPOHYDRO step, before the REFUGIUM step
  started, and only to scripts the running step had not yet loaded.)
* `ui/build_data.py` now exports scenario data as generic named series (`fp_series`, `tl_series`) with a favourable-direction flag
  per quantity; the UI prefers the dense grid over the validation grid when both exist (the old sort picked "validation" first -- a
  latent bug, never triggered because no dense result existed).

## 2026-10-05 -- no points on the maps; 'no data' bug; place explorer; elevation-aware surface

User feedback: the sample points must be removed from everywhere, and most maps said "no data".

* **Points removed.** The map no longer has a "show model nodes" option, node markers, nearest-node text or links. The old Site
  explorer (a dot-locator over the sampled cells) is replaced by the **Place explorer** (`ui/place.js`): click anywhere in Armenia
  or pick a marz, and every quantity is shown for that place today and under each emissions path (mean of the 5 models, band =
  lowest to highest model). Values at a place come from the same interpolated surface; a marz value is the area-weighted mean of it.
  Chart markers at 2019/2050/2080/2100 and the variogram scatter are chart data, not map points, and stay.
* **'No data' was my bug.** All 29 map quantities have 2019 data, but only 8 had scenario results, and the map defaulted every
  quantity to scenario mode, so 21 of 29 showed "No data". A quantity without scenario results now opens on its 2019 field with a
  note naming the script that will produce its scenarios; planned maps with no results at all are listed as disabled
  "no results yet". `ui/build_data.py` now writes `map_plan` (36 planned maps, ready / waiting, producing script) and the Status
  page shows it. Derived scenario series added where no run is needed (headroom below treeline, above own treeline).
* **Scenario coverage today vs after the dense run.** Now: viability x3, hazard x3, potential treeline, headroom, above-treeline,
  treeline shift (10 maps). After `run_future_projections.py` (dense): climate x5, per-species water deficit / stress integral /
  soil potential x9, P[V>=V*] x3, refugium score x3, robust refugium x3, hazard uncertainty x3 -- all 36. `fit_refugium_viability.py`
  stores the matching 2019 per-species fields. Generic-rooting-depth CWD/WSI/psi layers leave the map once per-species ones exist.
* **Elevation-aware surface.** Plain IDW cannot reproduce quantities that follow elevation: annual temperature on 25 nodes had
  leave-one-out R^2 = -0.23, RMSE 3.7 C. `scripts/build_dem_asset.py` averages the 30 m terrain raster onto the map grid
  (`ui/assets/map/elevation.png`, 256*R+G metres, small derived file; raw terrain stays on Drive). `ui/interp.js` gained an
  elevation-adjusted variant (each neighbour moved to the target's elevation with an OLS slope); it is used for a quantity only
  if its leave-one-out R^2 beats plain IDW by >= 0.02, the choice is fixed per quantity from its 2019 field so a baseline and its
  scenarios never differ by method, and both scores are shown on the map. Cross-checked against numpy in `tests/test_ui_interp.py`.

## 2026-10-05 -- hero video on the homepage; the project name is exactly one string

* **Name rule (user):** the project is "ANTAR — Assessment of Niche, Treeline & Analogue Refugia" everywhere, no other invented
  name. The browser title, the meta description, the hero heading, the README title and the package description now use it;
  the earlier "Armenia refugia explorer" (title) and "Armenia reforestation refugia" (hero heading) are gone.
  `tests/test_ui_name.py` pins the title and the hero heading to the exact string and fails if another name appears in the UI.
* **Hero video:** the supplied animation (`ui/assets/hero.mp4`, 10 s, 1280x720, H.264, 2.5 MB, moov atom first so it streams) is a
  logo reveal on a white background: birds and leaves cross the wordmark, the camera pushes into the "A" to show a forest valley
  inside it, a light flare sweeps through, and it returns to the wordmark, so it loops cleanly. It sits on a white stage
  (`--stage`, the footage's own background, in both themes) with the full-name heading, the description and one primary action
  below it. Poster frame `ui/assets/hero-poster.jpg` (62 kB) shows before playback and is the fallback.
  Behaviour: muted, looping, `playsinline`; not autoplayed under `prefers-reduced-motion` or Save-Data; a 44 px play/pause
  button (WCAG 2.2.2); pauses when scrolled out of view or the tab is hidden; a rejected `play()` is retried rather than abandoned;
  the video is `aria-hidden` because the heading carries the name. The clip has an unused AAC audio track (about 30 kB);
  there is no ffmpeg here to strip it, and it is muted.
* Phones: the sticky nav wrapped into five rows and covered about a third of the screen; it is now one horizontally scrolling row.
* Dev note: `python -m http.server` does not answer HTTP range requests, which Safari requires for video; GitHub Pages does.

## 2026-10-05 -- full-screen hero with the navbar on it; academic editorial redesign

User: the video must be fullscreen with the navbar on it, and the whole design cooler, more modern and academic.

* **Hero.** The video fills the viewport (`100svh`); the navbar is fixed on top of it, transparent with a soft white scrim and
  ink-coloured tokens in both themes (the footage is white), and turns into a frosted bar after 24 px of scrolling. The footage is
  white, so it is shown in a 16:9 frame that is centred and fills the viewport only when the crop would be under ~3% (aspect >=
  43:25); narrower viewports (laptops at 16:10, tablets, phones) show all of it and the white bars are invisible. The frame's top
  and bottom edges are feathered so the camera push-in at 5-6 s never shows a hard edge. The full name, a one-line descriptor, the
  two actions and the pause/play button sit at the bottom; the descriptor is hidden on viewports under 780 px tall so the copy never
  touches the wordmark. A first attempt that centred nothing (`place-items: center` on an auto grid track left the 1600 px frame at
  left 0 and cut the "R") was caught in the screenshot and fixed.
* **Design system** (from the ui-ux-pro-max search, checked against the product before use): style *Swiss Modernism 2.0* (grid,
  hairline rules, editorial) and the *Academic/Research* pairing -- Crimson Pro for headings and large figures. Its suggested body
  face, Atkinson Hyperlegible, was tried and dropped: it draws a slashed zero, so "2019" reads "2O/19". Replaced by Source Sans 3;
  IBM Plex Mono carries tabular figures. All three are SIL OFL, self-hosted in `ui/assets/fonts` (about 165 kB, no third-party
  font requests). The search's suggested pattern ("Portfolio Grid") does not fit a research tool and was not used.
* Tokens: warm paper background, ink text, teal/ochre/rust/navy from the logo; dark theme designed alongside (contrast >= 4.5:1 for text);
  flat cards with hairline borders; journal-style tables (rules, no boxes; numbers in mono); numbered section heads (01-04) on the home
  page; sections fade up on scroll (reduced-motion and a no-observer fallback show everything); skip link; focus moves to the page
  content on navigation; three-column footer; one primary action in the hero.
* The home heading "The six engines" was wrong (seven are listed, including treeline change); it now follows the count.
* **Hero text legibility (user screenshot):** at 5-6 s the camera pushes into the "A" and the letters and sparkles pass behind the
  heading. The heading, descriptor and rule now sit on a frosted translucent pane (`.hero-text`: 80% white, 14 px backdrop blur,
  hairline border); the secondary button gets the same translucent fill. Checked at the deepest zoom frame (5.9 s) on desktop and at
  phone width. Ink text over the darkest teal behind the pane is above 10:1, the descriptor above 5.5:1; browsers without
  `backdrop-filter` get a 94% white pane instead.

## 2026-10-05/06 -- dense TOPOHYDRO complete; the Mac slept for 17 hours

* **Dense TOPOHYDRO (Armenia-only, stride 7) finished 2026-10-05 19:00** after 2 h of streaming: 857 of 921 cells (93.1%) have output,
  above the UI's 90% completeness bar. The 64 without output lack soil or ERA5-Land data: 37 lie over water (Lake Sevan; no
  soil or land-reanalysis value exists there) and 26 are low-elevation land cells (738-1,013 m, a few higher), most likely built-up
  or masked land where SoilGrids has no value. Surfaces interpolate across them; the map no longer paints model values over open
  water, and regional means are weighted by land area.
* **What the dense grid buys (temperature, leave-one-out):** 25 nodes: plain IDW R2 = -0.23, elevation-adjusted R2 = 0.99
  (RMSE 0.37 C). 857 nodes: plain IDW R2 = 0.82, elevation-adjusted R2 = 1.00 (RMSE 0.09 C, spread of the runs sd 3.83 C).
* **Sleep:** dense XYLEM failed at 19:03 with a DNS error (`NameResolutionError` for oauth2.googleapis.com) as the network went
  away, and the machine then slept until 2026-10-06 ~12:45 (the first scenario members report 1,058 minutes of wall time).
  `caffeinate -dims` does not hold a MacBook awake with the lid closed or off power. Nothing was lost: the chain's retry loop and the
  per-member checkpoints resumed on wake. To avoid the gap: keep the lid open and the charger connected during runs.

## 2026-10-06 -- the "viability stays flat" finding was an artefact of constant atmospheric inputs; forcing rebuilt

User question: why are wind, radiation and humidity held at 2019 in the scenarios, and could that change?

**What was actually wrong (worse than the banner said).** The ERA5-Land exports on Drive hold ONE value per cell and year for wind,
shortwave, longwave and dew point. The pipeline repeated that annual mean on all 365 days (`np.full(n_days, ...)`) and left it
unchanged in every scenario. So the model had no seasonal cycle in evaporative demand (July shortwave is ~1.6x its annual mean) and
no humidity or radiation response to warming. My 2026-10-01 note that "XYLEM evidently derives its sensitivity from the variables held
fixed" was a guess that I never tested.

**Test** (`scripts/sensitivity_atmospheric_forcing.py`, 23 Armenian cells, same traits, 50 x 200 Monte Carlo, same seeds; the harness
reproduces the stored 2019 viability to 1e-16): mean climatic water deficit, broadleaf, 2019: 13.8 mm constant, 113.9 mm with the real
ISIMIP3a seasonal cycle (153 mm with my first approximation from Tmin; cell-wise r = 0.99 between the two). Its change under GFDL-ESM4
SSP5-8.5 by 2100: +3.6 mm constant vs +212 mm seasonal-and-scenario. Oak hazard change +0.003 vs +0.046; 23 of 23 oak cells exceed 1%
hazard by 2100 (19 under constant forcing). Minimum soil water potential -0.28 -> -0.65 MPa in 2019.

**Fix.** `scripts/pull_isimip_atmosphere.py` pulls, through the same ISIMIP Files API and Armenia bbox, monthly shortwave (rsds),
longwave (rlds), wind (sfcwind) and vapour pressure (hurs/100 x es(tas), built daily) from ISIMIP3a obsclim GSWP3-W5E5 (the 2015-2019
baseline shape) and ISIMIP3b (5 GCMs x 3 SSPs, the monthly change). `src/antar/climate/atmosphere.py` turns them into a daily factor
(mean exactly 1 over the baseline year, so ERA5-Land's cell-level annual mean is untouched) times the model's future/baseline monthly
ratio. Pressure stays an annual mean (its seasonal and scenario variation is about 1%). `ANTAR_CONSTANT_ATMOSPHERE=1` reproduces the old
behaviour. Tests: `tests/test_atmosphere.py`, `tests/test_forcing_atmosphere.py`.

**Also changed:** streamed terrain, soil and ERA5 inputs are now cached per grid in `data/_cache/` (they took 1-2 h per step on the dense
grid; every later step reads the cache). Dense TOPOHYDRO of 2026-10-05 and the validation-grid XYLEM/REFUGIUM/scenario results used the
constant forcing; they are being recomputed. The constant-forcing scenario file is archived as `future_projections_constant_atmosphere.yaml`.
Temperature, precipitation, GDD, season length and treeline do not depend on these inputs and are unaffected.
* **Home-page note removed (user request), and superseded numbers no longer shown.** Without the note the pages would have drawn the
  old constant-forcing water-balance numbers with no caveat, so `ui/build_data.py` now WITHHOLDS them: XYLEM, REFUGIUM, scenario and AEGIS
  files, and TOPOHYDRO's water-deficit / stress-integral / soil-potential layers, are used only if they were produced with the seasonal
  forcing (new files record `atmosphere`; older ones are judged by run_date >= 2026-10-06; the scenario file must say it). Withheld files
  are listed on the Status page's lineage table; temperature, rain, degree days, season length, frost and treeline are unaffected. The
  constant-forcing scenario file is kept at `configs/fitted/archive/future_projections_constant_atmosphere.yaml`. Maps return as the
  recomputed runs finish. Home stat cards without a result are hidden. The "Known gaps" entry on the Status page still records the finding.

## 2026-10-06 -- copyright and citation

User: add "copyright Narek Ohanyan, 2026, all rights reserved" and how the work should be cited.

* `LICENSE` replaced: the MIT licence added on 2026-10-01 contradicted "all rights reserved", so it is now "Copyright (c) 2026 Narek Ohanyan. All
  rights reserved.", with three clarifications: citing needs no permission; viewing and forking on GitHub is allowed to the extent GitHub's terms
  allow; third-party data, fonts (SIL OFL 1.1) and Leaflet keep their own licences. It also states that copies obtained under MIT between
  2026-10-01 and 2026-10-06 stay MIT -- a licence already granted cannot be withdrawn by changing the file. `pyproject.toml` now says
  `LicenseRef-All-Rights-Reserved` and names the author; the README's "MIT licensed" line is gone.
* Citation (one string, used everywhere): Ohanyan, N. (2026). ANTAR — Assessment of Niche, Treeline & Analogue Refugia (Version 2.0.0-alpha)
  [Computer software]. GitHub. https://github.com/Narek-Ohanyan/antar -- in the web footer (with copy buttons for the reference and the BibTeX),
  the README and `CITATION.cff` (GitHub's "Cite this repository"). No DOI is invented; `tests/test_citation.py` fails if the three disagree,
  if the name differs, or if a DOI string appears.
* Footer fix found while checking: `.foot` and `.cite` used a `padding` shorthand that removed the page gutter; now `padding-block`.
* Citation now leads with the web interface (https://antar.narekohanyan.com, the address the user will host it on) and gives the repository as
  "source code": reference text in the footer and README, BibTeX `url` + `note`, and `CITATION.cff` (`url` = site, `repository-code` = GitHub).
  Nothing else assumes the site is live yet (no canonical or Open Graph tags were added, and no `CNAME` file, since the hosting method is not
  decided). A persistent identifier would be a Zenodo snapshot of a tagged release; none has been created, so none is cited.

## 2026-10-06 -- About the author

Biography and portrait supplied by the author. `#/author` (linked from the footer, not the top nav, which already has ten items) carries the
portrait (`ui/assets/author.jpg`, resized from 2000 px to 720 px, 110 kB, alt text, declared size) and the bio as supplied, split into three
paragraphs without changing a word, with a link to https://www.narekohanyan.com (new tab, announced to screen readers). The home page gets a
short card (section 05) pointing to it; the README has an Author section; `CITATION.cff` carries the website on the author entry. No
institutional affiliation is attached to the software itself, and no ORCID or other identifier was invented. `tests/test_author.py` pins the
bio sentences, the link, the portrait's alt text and size. The portrait is now in the public repository.

## 2026-10-06 (cont.) -- validation-grid baseline recomputed with the seasonal atmosphere

TOPOHYDRO, XYLEM and REFUGIUM on the 80-point grid (25 Armenian cells) finished 13:34-13:40 from the cached inputs (no Drive access after the first
step). Checks: temperature, precipitation and GDD identical to the previous run (max difference 0, as they must be); generic Penman-Monteith water
deficit 166 -> 353 mm, minimum soil potential -1.03 -> -1.31 MPa; XYLEM mean hazard over the 23 Armenian cells, broadleaf 0.00099 -> 0.00777, oak
0.0254 -> 0.0608, equal to the offline harness (an independent code path) to the digits shown. The UI now takes water-balance layers from
whichever TOPOHYDRO file is current, so 40 of 41 planned maps have 2019 results (the 41st, treeline shift, exists only as a scenario). Scenario
maps return when the ISIMIP3b atmosphere download and the recomputed projections finish.

## 2026-10-06 (cont.) -- Acknowledgments written

Text drafted with the author, third person, two paragraphs, on the Acknowledgments page and in the README (replacing the hand-fill placeholder). The
project is named with its official name from the WSL project page, "Forest Restoration and Climate Change in Armenia (FORACCA)" -- the author's
draft said "Recreation"; funder "Swiss Agency for Development and Cooperation (SDC)" and the implementers (Forest Alliance, a consortium of Armenian
NGOs led by Shen NGO, WSL and FAO) were checked on that page (https://www.wsl.ch/en/projects/foracca/). NOT found online and kept in the author's
own words: "provides new climate services for Armenia, including high-resolution climate scenarios and local climate impact profiles for every
municipality" (the WSL page speaks only of "improved climate services"). The standing rule that FORACCA appears nowhere in the UI is narrowed by
`tests/test_ui_name.py`: it may appear only inside the Acknowledgments page, with its official name and without "Recreation".

## 2026-10-06 (cont.) -- Acknowledgments redesign, FORACCA section, centred navbar

* **Acknowledgments page.** The two paragraphs sat in the left half of a wide card. They are now margin-label rows (Swiss-editorial: label left, text
  right, hairline between rows) with a 68-character measure and 1.7 line height (the skill's 65-75 characters, 1.5-1.75), left-aligned with
  `text-wrap: pretty`; full justification was not used (uneven word gaps hurt readability, WCAG 1.4.8 advises against it), and automatic hyphenation
  is on only above 900 px because it fragmented words in a phone column. A new section "About the FORACCA project" has a lead statement, the three aims
  as a numbered list, a "Climate services" line and a fact sheet (full name, SDC, programme 10 years 2023-2033 CHF 10 million, main phase 2025-2028,
  implementers), with the sources linked and the check date shown. Facts verified 2026-10-06 against the WSL project page (name, funder, implementers,
  2025-2028) and Armenpress (10 million CHF, 2023-2033, three aims). NOT verified anywhere and kept in the author's words: "high-resolution climate
  scenarios and local climate impact profiles for every municipality".
* **Navbar:** items centred on the page at >= 1200 px (grid: brand | nav | theme button; measured centre 713 px = page centre 713 px at 1440 px).
  Between 1100 and 1200 px the items are centred between brand and button; below 1100 px (previously 900) they form one scrolling row, centred when
  they fit ("safe center"), because at 1024 px ten items overflowed their box and ran under the theme button -- a defect that predates this change.
  Tests: `tests/test_acknowledgments.py`.

## 2026-10-06 (cont.) -- Method page: formulas typeset, text matched to what was built

* **Formulas.** The engine cards on the Method page write their formulas as LaTeX between dollar signs (`ui/data/methodology.json`) and typeset
  them with KaTeX 0.19 (MIT), self-hosted in `ui/assets/katex` (woff2 fonts only, no third-party request), loaded only when the page is opened
  (`ui/math.js`). `tests/test_method_page.py` checks that every formula parses in KaTeX's strict mode (needs node; skipped without it) and that no
  Greek letter or pseudo-formula is left in plain text.
* **Text matched to the implementation** (author's request that the methodology describe what was done): the cards now say that scenarios use monthly
  change factors on the 2019 reference year (no quantile mapping), that the atmosphere is the ERA5-Land annual mean times an ISIMIP3a seasonal shape
  and the ISIMIP3b monthly ratio, that PET is three formulations, only $P_{50}$ and slope are sampled in XYLEM, and what is built but not applied.
  Two stale "known gaps" (viability nearly flat across scenarios; only temperature and precipitation perturbed) were removed because the constant-
  atmosphere artefact behind them was corrected today, and replaced by one entry saying scenario maps are withheld until the recomputed runs finish.
* **A full as-built specification page was written and then removed** at the author's request (LaTeX sources, converter, KaTeX-rendered numbered
  equations, tables and figures). Only the engine cards remain, with one addition the author asked for: the flow chart of the six engines
  (`scripts/build_architecture_svg.py` -> `ui/assets/architecture.svg`, inlined so it follows the light/dark theme), each box filled by how far
  the engine is built and inactive flows dashed. While reading the concept note for it, two facts were found and kept: the copy in
  `docs/concept_v2` was the 21 September draft, not the final 22 September note (the final one adds the thermal-treeline modifier, Eq. for f_tl, and
  the variable inventory); it has been refreshed from `foracca2_repo_skeleton.zip`, whose PDF is byte-identical to `FORACCA_2.0_Technical_Concept_v2.pdf`.
  And the concept note's stomatal-closure curve had the exponent's sign reversed (it would close stomata as water potential rises); the code uses
  the correct, increasing form.
* **Standing instruction (author): the methodology never mentions the earlier model or its version.** Enforced by a test over the Method text, the
  site scripts, the README and CITATION.cff. Third-party datasets whose own version is called v1 are in the data manifest and are not affected.
* `ui/data/methodology.json` is hand-written but had been ignored by git through the broad `data/` rule; it is now tracked (the generated
  data bundle next to it still is not).

## 2026-10-06 (cont.) -- scenario run recomputed; three nodes had a NaN water balance and reported hazard 0

* **Recomputed with the seasonal atmosphere:** the ISIMIP3a/3b atmosphere download finished (54 of 54) and the 45-member scenario run on the validation grid
  completed (`configs/fitted/future_projections.yaml`, `atmosphere` recorded; all 40+ map quantities now have scenario results). The method text stored in
  that file still said wind, radiation and humidity were held at 2019 values; corrected in the script and in the file.
* **A silent fault found while checking the dense variogram (it came out NaN):** three dense-grid nodes (lat/lon 40.871/43.754, 40.813/44.221,
  39.762/45.621) had an ERA5-Land wind value but no radiation or dew point. The validity mask tested only wind, so they ran with NaN inputs, their
  CWD/WSI/psi were NaN, and because `NaN >= 1` is False XYLEM reported a hazard of exactly 0 (viability 1.0) for them. TOPOHYDRO still wrote them as status
  `ok`. Fix: `era5_complete()` requires all five bands; the writer refuses to call a non-finite water balance `ok`
  (`skipped_nonfinite_water_balance`); `tests/test_era5_mask.py`. The validation grid and the validation-grid scenario file contain no NaN, so they are
  unaffected. The dense TOPOHYDRO, XYLEM and REFUGIUM outputs committed earlier today contain the three bad nodes; they are being recomputed (dense chain
  restarted, MNEME dense with it) together with the dense scenario run, which was restarted from scratch with the fixed mask.

### Correction, 2026-10-07: the cause of the three spurious zero hazards was the soil, not ERA5-Land

The entry above blamed missing ERA5-Land radiation and dew point. Checking the cached inputs showed all five ERA5-Land bands and the ISIMIP atmosphere factors are finite at the
three nodes; what is NaN is the soil hydraulic parameters (theta_sat, theta_fc, theta_lim, psi_sat, b) derived from SoilGrids. `valid_soil` tested only the clay value, so the nodes
ran. Real fix: `soil_complete()` requires every raw soil value and every derived parameter to be finite, at extraction and when the cached inputs are loaded
(`tests/test_node_masks.py`). The ERA5 check (`era5_complete()`) stays as a sensible guard but was not the cause. The XYLEM and REFUGIUM dense outputs and the dense scenario run
were computed with the old mask; rather than recompute hours of Monte Carlo, `scripts/drop_invalid_nodes.py` removes the three nodes TOPOHYDRO rejects from them and recomputes the
summary fields (nodes are independent and seeded by grid index, so the remaining nodes are exactly what a rerun would give). It has been applied to the XYLEM and REFUGIUM
dense files (854 cells per group now); the dense scenario file gets the same treatment when that run finishes. The dense MNEME panel was restarted so that it uses the fixed mask.

## 2026-10-07 -- late-frost days redefined; the layer is withdrawn; the interpolation score is no longer shown on the maps

* **Definition.** Late-frost days used to be every day below -2 degC from budburst (cumulative GDD >= 200) to the end of the year, so the late-autumn and winter frosts
  were counted (mean 39 days on the validation grid, 46 on the dense grid). They are now the days below -2 degC between budburst and the warmest day of the year
  (`indices.warmest_day_index`: the peak of a 15-day centred running mean of the daily mean temperature, ties to the earlier day; `indices.late_frost_days(..., end_day)`).
  The choice of a 15-day window is mine, not taken from a source.
* **Recomputed, not rerun.** The count depends on the daily temperatures only, so `scripts/recompute_late_frost.py` rebuilds it from the reference weather with
  `forcing.late_frost_days_from_reference` (pinned to the full pipeline by `tests/test_late_frost.py`). Before writing, both definitions are computed for every node and each stored
  count must equal one of them: all 78 validation nodes, all 854 dense nodes and all 3,510 validation-grid scenario rows reproduced the old definition exactly. Applied to
  `topohydro_grid_run_2019.yaml`, `topohydro_grid_run_2019_dense.yaml` and `future_projections.yaml`; files carry `late_frost_definition: budburst_to_warmest_day`. The dense scenario
  file gets the same treatment when that run finishes (its workers started before the change and still count the old way).
* **Result.** With the 200 degC-day placeholder the new count is zero at every validation node in 2019 and at 9 of 854 dense nodes (1-2 days); with the budburst threshold at 25 degC-days
  261 dense nodes have a frost day. (I first reported "every dense node is zero"; that was wrong, the dense grid has those 9 nodes. The conclusion that the layer carries no
  information at the placeholder is unchanged.) Measured in `configs/fitted/sensitivity_placeholders.yaml` (`late_frost_budburst_to_warmest_day`, `scripts/sensitivity_late_frost.py`).
* **Decision (user):** apply the new definition and take "Late-frost days" out of the map menu, because a map that is zero almost everywhere shows nothing. The count is still
  computed and stored; the placeholder entry on the Status page states the measurement and what is needed (a species-specific budburst threshold).
* **Maps:** the "Interpolation check" block (leave-one-out R2 and RMSE per selection) is removed from the map panel on the user's request; the note under the map still says the surface is an
  interpolation. The Status & limits page no longer lists "Which maps have results" or "Data lineage", and Treeline change no longer carries the method-check sentence.

## 2026-10-07 (cont.) -- visual refresh of the web interface

* **Home page.** The "Headline results" are now a dark findings band with an emissions-path and a horizon selector. For the chosen pair each card shows the mean
  one-year hydraulic survival of a species group (or the treeline shift), the change against 2019 in percentage points (sign and glyph, not colour alone), the range
  across the five climate models, and a trend line for all three paths (the chosen one solid with its range shaded). Every number is read from the manifest
  (`Home.headlineModel`, tested); the trend charts have a minimum visible span, as the map colour scale does, so a rounding-level difference is not drawn as a drop.
  The lead sentence is generated from the same numbers. A relief of Armenia (contours every 250 m from the project's terrain raster) carries two highlighted
  contours: the mean potential treeline elevation over the model nodes today and after the chosen shift. It is labelled as an illustration of the relief the models run on.
  The evidence base (nodes, climate members, species groups, engines, datasets) is a ruled line under the overview text.
* **Photographs.** Four Wikimedia Commons photographs (user-approved, CC BY-SA), captioned neutrally (no species is named that the photograph does not establish),
  credited on the image, in the enlarged view and on the Acknowledgments page. Resized only; cropped by the page layout.
* **Palette and surfaces.** Deeper spruce/navy "ink" for the band and footer, warm paper for pages, 8 px radius; the emissions paths keep their chart colours everywhere.
* **Search.** Ctrl K / ⌘K / `/` opens a palette (pages, map layers, marz, two actions); ranking is a pure tested function.
* **Corrections found on the way.** The Method page's plain-language text still said humidity and wind were not changed in the scenarios, which contradicted the
  technical text and the corrected pipeline; fixed. The home page said "seven engines" because the list included Treeline change, which is part of MERISTEM; it now
  lists the six engines. Engine states are stored once in `methodology.json` and read by both the cards and the flow chart.

### Second pass on the visual refresh (same day)

* Audit against the design checklist: no horizontal overflow at 375, 820 and 1440 px; text and control colour pairs measured (`tests/test_contrast.py`, both themes, all >= 4.5:1);
  photo dialog, search palette and scenario selector checked by keyboard; reduced motion honoured by the global rule.
* Added: the relief figure's space is reserved before it loads (no layout jump); the chosen path and horizon live in the address and are restored from it; a CSV of every headline
  number (all paths and horizons, with the model range and the grid); credited photograph banners on four inner pages.
* **Decision page.** It was a bare line, "No portfolio has been computed", with no title: the AEGIS result is withheld on purpose (it was built on the superseded constant atmosphere).
  It now says so, in plain words, and describes what the page will show; the numbers return after the dense AEGIS rerun.

## 2026-10-07 (cont.) -- deployment to antar.narekohanyan.com

* **Route chosen by the author:** the built site goes into the `antar/` folder of the website repository, whose Hostinger auto-deploy then puts it in `public_html/antar`, the document root of
  the subdomain the author created. (A separate `deploy` branch with a Hostinger Git deployment was considered and dropped.) The folder is also reachable as `narekohanyan.com/antar/`,
  so `deploy/.htaccess` redirects every other spelling to https://antar.narekohanyan.com/.
* **Checked:** a fresh clone of this repository builds an identical site with `python3 ui/build_data.py` alone (numpy and pyyaml); the bundle is 72 files, 7.0 MB; the `.htaccess` rules were
  tested on a real Apache 2.4 with TLS (`tests/test_site_bundle.py`): redirects keep the query string, `/index.html` goes to `/`, the page is always revalidated, scripts and data are cached
  for an hour, images and fonts for a week, text is compressed, video supports ranges, folders are not listed. The publisher is guarded and tested on temporary repositories.
* **Confirmed by the author:** the Acknowledgments sentence about local climate impact profiles for every municipality; AEGIS costs other than planting stay blended averages, labelled as such.

## 2026-10-07 (cont.) -- dense scenario run finished and processed

* The 45-member scenario run on the dense grid finished (13:23): 45 members x 3 species groups x 857 cells. `scripts/drop_invalid_nodes.py --files scenarios` removed the 3 nodes TOPOHYDRO rejects
  (405 cell records, 3 x 3 x 45), leaving 854 per group, as in the XYLEM and REFUGIUM dense files.
* `scripts/recompute_late_frost.py --apply --only dense` rewrote late-frost days under the budburst-to-warmest-day definition: all 38,430 rows matched one of the two definitions before writing
  (mean 34.4 -> 0.0 days); the file carries `late_frost_definition`.
* Headline numbers on the dense grid (SSP5-8.5, 2100, ensemble mean, change against 2019): broadleaf 98.0% (-1.2 pp), oak 90.7% (-4.0 pp), pine 100.0% (0.0 pp); on the 25-node grid they were -1.3, -4.7 and 0.0 pp.
  The treeline change is still the 25-node result (`compute_treeline_change.py --dense` needs the Drive terrain), so the cards now name the number of nodes behind each of them.
* Still open: dense treeline change, dense AEGIS (and its text), the dense MNEME panel, the MERISTEM refit.
