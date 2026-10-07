# ANTAR — Assessment of Niche, Treeline & Analogue Refugia

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

## The six engines

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

## Reproducing the figures

```bash
python scripts/make_figures.py --out docs/figures  # the hydraulic-engine and extrapolation-experiment figures
```

The figures of the hydraulic engine and of the controlled extrapolation experiment use placeholder traits and synthetic
data by design: they demonstrate mechanisms, not skill.

## The web interface: artwork and photographs

* `ui/assets/relief.svg` (the contour lines behind page titles, the footer and the home page) is drawn from the project's own terrain raster by
  `python scripts/build_relief_art.py`; it is terrain, not a model result.
* `ui/assets/photos/` holds four third-party landscape photographs from Wikimedia Commons (CC BY-SA 3.0/4.0), resized to 800 and 1600 px wide.
  `credits.json` records author, licence, date and source page for each; the build refuses a photograph without an entry, and the credit is shown on
  the image, in the enlarged view, on the Acknowledgments page and (as a general notice) in the footer. To add one, put `name-800.jpg` and
  `name-1600.jpg` next to it, add its entry to `credits.json`, and run `python ui/build_data.py`; `tests/test_photo_credits.py` checks the rest.
* The state of each engine (run, reduced scope, built but not applied) is written once, in `ui/data/methodology.json`; the home-page cards and the flow
  chart on the Method page (`scripts/build_architecture_svg.py`) both read it.
* Ctrl K or ⌘K (or `/`) opens a search box that jumps to any page, map layer or marz.
* The Treeline, Decision, Models and Method pages open with one of the same photographs behind the title, credited on the banner. The headline band's
  emissions path and horizon are kept in the address (`#/home?ssp=ssp370&hz=2080`), and a button downloads every headline number as CSV.
* `tests/test_contrast.py` reads the colour tokens from `ui/styles.css` and checks WCAG contrast for the pairs the pages use, in both themes.

## Conventions

* Every formula lives in tested code; tests include closed-form checks (two-phase engine against `B / E_min`, Ishigami
  function for Sobol' indices, mass balance to machine precision), physical bounds (e.g. `ΔVPD/ΔTmax`) and regression
  tests that encode known construction errors.
* Validation is specified before fitting (`configs/cv.yaml`); every fitted quantity, including the AOA threshold and
  stacking weights, is fitted inside the training fold only.
* Data enter through manifests (source, version, checksum, licence); nothing is downloaded implicitly.

## Status

Version 2.0.0-alpha.

## Copyright

© 2026 Narek Ohanyan. All rights reserved. See `LICENSE`. Third-party data, fonts and libraries keep their own licences (listed in
`configs/manifests/` and on the References page of the web interface).

## How to cite

Citing the work needs no permission. Please cite it as:

> Ohanyan, N. (2026). *ANTAR — Assessment of Niche, Treeline & Analogue Refugia* (Version 2.0.0-alpha) [Computer software and web interface].
> https://antar.narekohanyan.com (source code: https://github.com/Narek-Ohanyan/antar)

```bibtex
@misc{ohanyan2026antar,
  author       = {Ohanyan, Narek},
  title        = {{ANTAR} --- Assessment of Niche, Treeline \& Analogue Refugia},
  year         = {2026},
  version      = {2.0.0-alpha},
  howpublished = {Computer software and web interface},
  url          = {https://antar.narekohanyan.com},
  note         = {Source code: https://github.com/Narek-Ohanyan/antar}
}
```

GitHub's "Cite this repository" button reads `CITATION.cff`. Please also cite the datasets your use depends on (see the References page).

## Author

Narek Ohanyan — climate & environmental researcher at the AUA Acopian Center for the Environment, working on modelling forest climate
resilience. Full biography on the web interface ("About the author") and at https://www.narekohanyan.com.

## Acknowledgments

The author would like to express sincere gratitude to the Swiss Federal Institute for Forest, Snow and Landscape Research WSL, and in particular to Franziska Zilker, Tobias Kühnhanss and Dr. Michael James McCarthy of the Dynamic Macroecology group, for providing datasets and bias-corrected environmental data and for their technical feedback on the methodology. The author also thanks his supervisor, PD Dr. Marco Pütz, and the coordinator, Dr. Dominik Braunschweiger, for making his guest scientist visit at WSL possible.

The author further acknowledges Alen Amirkhanian, Director of the AUA Acopian Center for the Environment, and the wider team of the Forest Restoration and Climate Change in Armenia (FORACCA) project for their logistical support and collaborative insights during the conceptualization of ANTAR.

### About the FORACCA project

[Forest Restoration and Climate Change in Armenia (FORACCA)](https://www.wsl.ch/en/projects/foracca/) is a Swiss-funded programme that supports reforestation on community lands, climate-smart forest management and climate-resilient development in Armenia. It is funded by the Swiss Agency for Development and Cooperation (SDC) (10 years, 2023–2033, CHF 10 million; main phase 2025–2028) and implemented by the Forest Alliance, a consortium of Armenian NGOs led by Shen NGO, the Swiss Federal Research Institute WSL and the Food and Agriculture Organization of the United Nations (FAO). Its aims are to advance scientific understanding of Armenia's capacity to address climate change and sustainably manage its forests; to promote climate-smart practices in rural areas; and to ensure evidence-based policymaking for climate adaptation and efficient forest management. The project provides new climate services for Armenia, including high-resolution climate scenarios and local climate impact profiles for every municipality.

Sources: [WSL](https://www.wsl.ch/en/projects/foracca/), [Armenpress](https://armenpress.am/en/article/1126549); checked 2026-10-06.
