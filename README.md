# ANTAR — Assessment of Niche, Treeline & Analogue Refugia

ANTAR is a hybrid process-statistical framework for finding where climate-resilient forest restoration will hold in Armenia. Six engines are chained: a daily water balance
(TOPOHYDRO), a mechanistic hydraulic-failure hazard (XYLEM), a statistical dieback hazard (MNEME), species niches and the climatic treeline (MERISTEM), viability and robust refugia
under 45 climate scenarios (REFUGIUM), and a budget-constrained, scenario-robust planting portfolio (AEGIS). Results are on a dense grid of 854 Armenian nodes (921 for the treeline),
about 6.5 km apart, and are explored at **https://antar.narekohanyan.com**. Every result carries its caveats: what is a placeholder, what is not fitted, and why.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23236944.svg)](https://doi.org/10.5281/zenodo.23236944)

Version 2.0.0. The honest summary of what it can and cannot say is on the site's *Models & validation* and *Status & limits* pages and in `IMPLEMENTATION_LOG.md`; in short:

* **Run on real data:** the water balance and hydraulic hazard for 2019, 45 scenario members (5 climate models × 3 emissions paths × 3 horizons), the climatic treeline shift, species niches
  (refitted with the real water deficit), refugium criteria (a), (b) and (c), the variance partition of the ensemble, and the CVaR-robust planting portfolio with species niches and an optional cap on the share of one species group.
* **Built but not fitted:** MNEME's hazard model. The satellite vitality record finds 162 dieback onsets in 570,175 forest pixel-years, but cannot tell them from changes in the record itself (92 pixels with a persistent decline against 1,370 with a persistent rise), so no hazard is claimed. In its place a vitality-response analysis finds no detectable drought response overall and supports XYLEM's drought-stress ranking for oak only.
* **Placeholders:** several hydraulic traits, the budburst threshold behind late frost, the calm-clear-night fraction; they are listed with their measured effect on the Status page.
* **Not modelled:** growth and height (MERISTEM), stand connectivity, method-specific survival, water use. Deferred to version 2.1: TRY plant traits and CORDEX regional projections.

## Quick start

```bash
pip install -e ".[dev]"      # numpy, scipy, pandas, scikit-learn, pyyaml (+ pytest, matplotlib)
python -m pytest             # about a minute on a fresh checkout (it builds the web data first and starts a real Apache to test the deployment rules); tests that need the raw rasters are skipped
python ui/build_data.py      # builds the web data from configs/fitted/ (needs only pyyaml and numpy)
python scripts/serve_ui.py   # http://localhost:8765
```

Optional extras: `.[geo]` (rasters, zarr), `.[bayes]` (PyMC), `.[ml]` (boosting, SHAP, conformal), `.[workflow]` (Snakemake, DVC).
The fitted results the site shows are committed in `configs/fitted/`; the raw and derived rasters are not (see Data and reproduction).

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

## Data and reproduction

Raw and derived rasters live in a local, gitignored `data/` directory and are tracked by manifest (`configs/manifests/`: source, version, URL, citation, licence, checksum), never by path.
Large rasters are streamed from Google Drive (exports of Google Earth Engine; credentials are the Earth Engine OAuth file in `~/.config/earthengine/`, never stored in the repository).
Reproducing everything from scratch takes days of compute; the order on the dense grid is:

```bash
python scripts/run_topohydro_grid.py --dense                       # daily water balance of the 2019 climate at the nodes
python scripts/fit_xylem_mechanistic_hazard.py --dense             # hydraulic-failure hazard
python scripts/fit_refugium_viability.py --dense                   # viability, risk-averse score, criterion (a)
python scripts/run_future_projections.py --dense --workers 5       # the 45 scenario members (about a day)
python scripts/compute_treeline_change.py --dense                  # climatic treeline shift
python scripts/drop_invalid_nodes.py                               # removes nodes the water balance rejected
python scripts/estimate_variogram.py --dense && python scripts/integrate_ecosystem_map.py
python scripts/compute_refugium_criteria.py && python scripts/compute_uncertainty_partition.py
```

MNEME and MERISTEM (the Drive-heavy steps; `scripts/run_meristem_then_mneme_all.sh` chains them so that two jobs never read Drive at once):

```bash
python scripts/extract_mneme_forest_pixels.py --dense --pixels-per-node=100000   # forest pixels of every climate cell: kNDVI, harvest/fire layer, height
python scripts/fit_mneme_hazard_panel.py --dense --pixels-per-node=100000        # climate of the cells, event labels, mirrored-series check, hazard fit when it is allowed
python scripts/fit_mneme_vitality_response.py                                     # vitality response to drought and the check of XYLEM's ranking (local)
python scripts/compute_real_cwd_for_meristem.py                                   # real water-balance deficit at the 3,524 species points
python scripts/fit_meristem_adult_niche.py && python scripts/apply_meristem_niche.py
python scripts/fit_aegis_portfolio.py --dense                                     # the portfolio, with the niche and the species-group caps
python ui/build_data.py
```

The harvest/fire layer is `antar_disturbance_ancillary_v2` (the first export was wrong; see `IMPLEMENTATION_LOG.md`). `docs/mneme_data_request.md` says what observed dieback records would let MNEME be fitted and who may hold them.

### Figures

```bash
python scripts/make_figures.py --out docs/figures  # the hydraulic-engine and extrapolation-experiment figures
```

They use placeholder traits and synthetic data by design: they demonstrate mechanisms, not skill.

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

## Publishing the web interface

The site is served at https://antar.narekohanyan.com from the folder `antar/` of the website repository
(https://github.com/Narek-Ohanyan/narekohanyan.com), which Hostinger deploys on every push to `main` into `public_html/antar`.

```bash
python scripts/publish_site.py --site-repo <clone of narekohanyan.com>          # builds, checks, shows what would change
python scripts/publish_site.py --site-repo <clone of narekohanyan.com> --push   # commits only antar/ and pushes: this goes live
```

`scripts/build_site.py` assembles the bundle (`ui/` without Python files, plus `deploy/.htaccess` and `robots.txt`) and refuses an incomplete one.
`deploy/.htaccess` keeps one address (https, the subdomain, no `/index.html`; the same folder is also reachable as `narekohanyan.com/antar/`), sets the cache
rules and compression, and is tested against a real Apache in `tests/test_site_bundle.py`. `tests/test_publish_site.py` checks that the publisher changes only
`antar/`, never pushes without `--push`, and refuses a dirty, stale or wrong clone.

## Conventions

* Every formula lives in tested code; tests include closed-form checks (two-phase engine against `B / E_min`, Ishigami
  function for Sobol' indices, mass balance to machine precision), physical bounds (e.g. `ΔVPD/ΔTmax`) and regression
  tests that encode known construction errors.
* Validation is specified before fitting (`configs/cv.yaml`); every fitted quantity, including the AOA threshold and
  stacking weights, is fitted inside the training fold only.
* Data enter through manifests (source, version, checksum, licence); nothing is downloaded implicitly.

## Status

Version 2.0.0 (2026-10-08).

## Copyright

© 2026 Narek Ohanyan. All rights reserved. See `LICENSE`. Third-party data, fonts and libraries keep their own licences (listed in
`configs/manifests/` and on the References page of the web interface).

## How to cite

Citing the work needs no permission. Please cite it as:

> Ohanyan, N. (2026). *ANTAR — Assessment of Niche, Treeline & Analogue Refugia* (Version 2.0.0) [Computer software and web interface].
> https://antar.narekohanyan.com (source code: https://github.com/Narek-Ohanyan/antar). https://doi.org/10.5281/zenodo.23236944

```bibtex
@misc{ohanyan2026antar,
  author       = {Ohanyan, Narek},
  title        = {{ANTAR} --- Assessment of Niche, Treeline \& Analogue Refugia},
  year         = {2026},
  version      = {2.0.0},
  howpublished = {Computer software and web interface},
  url          = {https://antar.narekohanyan.com},
  doi          = {10.5281/zenodo.23236944},
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

