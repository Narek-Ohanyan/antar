"""ANTAR (Assessment of Niche, Treeline & Analogue Refugia): a hybrid
process-statistical framework for predicting climate-resilient reforestation
in Armenia.

Engines
-------
climate      TOPOHYDRO - topoclimate, PET, snow, soil-water balance, drought indices
hydraulics   XYLEM     - two-phase plant-hydraulic failure engine
hazard       MNEME     - discrete-time mortality hazard with drought-legacy lag
niche        MERISTEM  - niche, attainable height, growth, establishment
viability    REFUGIUM  - cohort viability, competing risks, robust refugia
validation   Blocked CV, metrics, area of applicability, conformal, design-based accuracy
uncertainty  Ensemble variance decomposition and Sobol' indices
decision     AEGIS     - robust (CVaR) site-species portfolio optimisation

Everything in here is written for clarity and testability first.  Heavy I/O
(GEE exports, Zarr/COG handling) lives in ``antar.io`` as documented stubs.
"""

__version__ = "2.0.0-alpha"
