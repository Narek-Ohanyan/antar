# antar — ANTAR repository skeleton

ANTAR (Assessment of Niche, Treeline & Analogue Refugia) is a hybrid process-statistical framework for
predicting where climate-resilient reforestation will hold in Armenia. This repo fixes the interfaces
of the framework and implements, with tests, the numerical kernels the rest depends on.
It is a **skeleton**, not the finished framework: data access and the full-scale pipeline steps are declared
but not implemented, and the species traits in `configs/species_traits.csv` are **placeholders** that exercise
the code. Nothing produced from them is a result.

## Quick start

```bash
pip install -e ".[dev]"      # numpy, scipy, pandas, scikit-learn, pyyaml (+ pytest, matplotlib)
python -m pytest             # a few seconds
```

Optional extras: `.[geo]` (rasters, zarr), `.[bayes]` (PyMC), `.[ml]` (boosting, SHAP, conformal), `.[workflow]` (Snakemake, DVC).

## The five engines

| Engine | Package | What it holds |
|---|---|---|
| **TOPOHYDRO** | `antar.climate` | station-fitted lapse-rate/cold-air-pooling and precipitation-gradient downscaling, dew-point VPD, slope radiation, snow, soil bucket with mass-balance closure, a PET *ensemble* (not one formulation) carried through to CWD/WSI, growing-season/GDD/late-frost indices, SPEI, quantile delta mapping, climate analogues, and a per-cell daily orchestration (`climate.forcing.topoclimate_forcing`) composing all of the above |
| **XYLEM** | `antar.hydraulics` | vulnerability curve and inverse, biphasic `g_min(T)`, two-phase plant-hydraulic failure engine (closed form and daily HFI recursion), the full two-level Monte Carlo failure probability (outer trait-hyperparameter loop + inner individual loop), stress-spell summary features, a monotone-emulator trainer with the spec's own held-out release gate, link-scale recalibration, and `hydraulics.pipeline.mechanistic_hazard_for_cell` wiring a TOPOHYDRO cell directly into it |
| **MNEME** | `antar.hazard` | person-period panel with Mundlak within/between split, rare-event (King & Zeng) weighting, misclassification-corrected observed hazard, a persistence-checked dieback-event rule, an elastic-net cloglog GLM, monotone GBM, NNLS stacking, DLNM cross-basis (drought-legacy lags) and an age-spline term, the space-for-time adaptation bracket, `hazard.pipeline` wiring the full Sec. 7.3 design-assembly + nested-CV-stacked-learners + gate end to end, and the area-of-applicability gate (via `antar.validation.aoa`) handing over to XYLEM outside the training domain |
| **MERISTEM** | `antar.niche` | adult-niche Boyce index, attainable-height quantile model with a water-limited (AET/PET) area-of-applicability fallback, Chapman–Richards growth in physiological age with all four Eq. 8.2 modifiers (water, growing-degree-day, late-frost, and the thermal-treeline ceiling — kept as separate, independently testable functions), the diagnostic potential treeline elevation, an ensemble `P[H_T ≥ H_min]` evaluator, and establishment hazard |
| **REFUGIUM** | `antar.viability` | competing-hazard cohort viability (composes directly with MERISTEM's `P[H_T≥H_min]`), robust refugia under the full Eq. 8.6 conjunction (ensemble viability, topographic buffering, area of applicability), risk-averse refugium scoring; climate-analogue matching lives in `antar.climate.analogs` |
| **AEGIS** | `antar.decision` | CVaR-robust portfolio optimisation (MILP via HiGHS) across the scenario/model/parameter ensemble, out-of-sample evaluator, the mean–CVaR efficient frontier (price of robustness), and an extrapolation-footprint diagnostic on the solved plan |

Supporting, cross-cutting packages:

| Package | What it holds |
|---|---|
| `antar.validation` | block / forward-chaining / group splits, variogram range, dissimilarity index and AOA, proper scoring rules, calibration, conformal intervals, design-based estimators |
| `antar.uncertainty` | ANOVA variance fractions, Sobol' indices |
| `antar.io` | master grid (EPSG:32638, 30 m); data manifests (source, version, URL, citation, licence, checksum) so no dataset enters through a hardcoded path; Earth Engine exports covering most of Table 4's realistically-obtainable variables (vitality/kNDVI + LandTrendr, disturbance, canopy structure, terrain, soils, ERA5-Land forcing, vegetation state, snow, land tenure) |

`configs/` holds the study area, scenario/ensemble design, validation design, hydraulic-trait placeholders
(`species_traits.csv`) and the LAI placeholder standing in for MERISTEM until it exists (`stand_defaults.yaml`);
`workflow/Snakefile` declares the contract between engines (commands are placeholders). Raw and derived rasters are
never committed: they live in a local, gitignored `data/` directory and are tracked by manifest
(`antar.io.manifest`) rather than by path.

## Reproducing the v1 audit and the figures

```bash
python scripts/audit_v1.py --parquet <path to Armenia_ML_Training_Data.parquet> --out docs/audit_v1.json   # ~2 min
python scripts/make_figures.py --audit docs/audit_v1.json --out docs/figures
```

`audit_v1.py` re-runs a prior-project random forest and the checks behind Table 1 and Fig. 1 of the concept
note (random vs. blocked CV, coordinates-only baseline, features predicting location, flat extrapolation, unit
audit). Figures 3 and 4 use placeholder traits and synthetic data by design.

## Conventions

* Every formula lives in tested code; tests include closed-form checks (two-phase engine against `B / E_min`, Ishigami
  function for Sobol' indices, mass balance to machine precision), physical bounds (e.g. `ΔVPD/ΔTmax`) and regression
  tests that encode known prior-project errors.
* Validation is specified before fitting (`configs/cv.yaml`); every fitted quantity, including the AOA threshold and
  stacking weights, is fitted inside the training fold only.
* Data enter through manifests (source, version, checksum, licence); nothing is downloaded implicitly.

## Status

Version 2.0.0-alpha. Licence and authorship statement to be finalised before any public release.

## Acknowledgments

<!-- prior project name/attribution to be filled in by hand -->
