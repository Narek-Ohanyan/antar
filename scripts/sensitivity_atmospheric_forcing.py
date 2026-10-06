"""How much do the atmospheric inputs the model currently holds fixed matter?

The model's wind, shortwave, longwave, dew point and pressure are ERA5-Land ANNUAL means: one number per cell,
repeated for all 365 days, and unchanged in every future scenario. This script re-runs the same hydraulic model
(same traits, same 50 x 200 Monte Carlo, same seeds) with a physically based seasonal / warming-responsive version
of the three that drive evaporative demand, and compares:

  V0  current:  Rs, Rl, ea constant through the year, and in the future
  V1  proposed: * ea(day) = es(Tmin(day) - k)   (FAO-56 Annex 6: dew point tracks the daily minimum), with k fitted
                  per cell so that the annual mean of ea equals ERA5-Land's annual dew-point value
                * Rs(day) = annual Rs x Ra(day, lat) / mean(Ra)   (seasonal shape from the extraterrestrial radiation,
                  constant cloudiness)
                * Rl(day) = c x eps_clear(ea, T) sigma T^4   (Brutsaert 1975 clear-sky emissivity), c fitted per cell so
                  the annual mean equals ERA5-Land's annual longwave
              In a future climate Tmin and T rise, so ea and Rl rise with them (Rs keeps its shape).
  wind and pressure: unchanged in both.

Both are run for 2019 and for one future member (GFDL-ESM4, SSP5-8.5, 2100). Inputs: the pickled context
`run_future_projections.py` writes (--ctx), so no network is needed.
"""
import json
import pickle
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_future_projections import compute_deltas  # noqa: E402
from run_topohydro_grid import wind_speed_2m  # noqa: E402
from fit_xylem_mechanistic_hazard import PET_FORMULATION, OUTER_DRAWS, INNER_DRAWS  # noqa: E402

from antar.climate import downscale  # noqa: E402
from antar.climate.atmosphere import AtmosphereShape  # noqa: E402
from antar.climate.forcing import topoclimate_forcing  # noqa: E402
from antar.climate.radiation import extraterrestrial_radiation, net_radiation_from_era5  # noqa: E402
from antar.climate.vapour import saturation_vapour_pressure  # noqa: E402
from antar.hydraulics.monte_carlo import two_level_failure_probability  # noqa: E402
from antar.hydraulics.twophase import simulate_two_phase  # noqa: E402
from antar.io.armenia_mask import inside_armenia  # noqa: E402

SIGMA_MJ = 4.903e-9      # Stefan-Boltzmann, MJ m-2 day-1 K-4


def fit_k(tmin_c, ea_target_kpa):
    """Offset k (K) such that mean_day es(tmin - k) = ea_target (bisection; es is monotone in k)."""
    lo, hi = -10.0, 40.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if saturation_vapour_pressure(tmin_c - mid).mean() > ea_target_kpa:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def forcing(c, i, deltas, horizon, w_max, variant):
    t_ref, tx_ref, tn_ref, p_ref = (c[k][:, i] for k in ("t_mean_ref_all", "t_max_ref_all", "t_min_ref_all", "p_ref_all"))
    doy, month = c["doy"], c["month"]
    wind10, ssrd, strd, dewpoint_k, pressure_pa = (a[i] for a in c["era5_2019"])
    ea_ann = float(saturation_vapour_pressure(dewpoint_k - 273.15))
    lat = float(c["lats"][i])
    gamma_k_per_m, precip_gradient_per_m = c["lapse"]
    gamma_of_day = gamma_k_per_m[month - 1]
    if deltas is None:
        t, tx, tn, p = t_ref, tx_ref, tn_ref, p_ref
    else:
        m = month - 1
        t, tx, tn = t_ref + deltas[(horizon, "tas")][:, i][m], tx_ref + deltas[(horizon, "tasmax")][:, i][m], tn_ref + deltas[(horizon, "tasmin")][:, i][m]
        p = p_ref * deltas[(horizon, "pr")][:, i][m]
    n_days = len(doy)
    t_for_rn = downscale.downscale_temperature(t, c["elevation"][i], c["z_ref_m"][i], gamma_of_day)
    if variant == "V0":
        ea = np.full(n_days, ea_ann)
        rs_j, rl_j = np.full(n_days, ssrd), np.full(n_days, strd)
    elif variant == "V2":                                           # real ISIMIP seasonal shape and scenario change
        atm = AtmosphereShape(Path(__file__).resolve().parent.parent / "data" / "isimip3b")
        f = atm.baseline_factors(lat, c["lons"][i], doy) if deltas is None else atm.scenario_factors(lat, c["lons"][i], doy, "gfdl-esm4", "ssp585", horizon)
        ea, rs_j, rl_j = ea_ann * f["ea"], ssrd * f["rs"], strd * f["rl"]
        wind10 = wind10 * f["wind"]
    else:
        k = fit_k(tn_ref, ea_ann)                                    # fitted on the 2019 reference, reused in the future
        ea = saturation_vapour_pressure(tn - k)
        ra = extraterrestrial_radiation(np.radians(lat), doy)
        ra_ref = extraterrestrial_radiation(np.radians(lat), doy)
        rs_j = ssrd * ra / ra_ref.mean()
        def rl_clear(tair_c, ea_kpa):
            tk = tair_c + 273.15
            return 1.24 * (ea_kpa * 10.0 / tk) ** (1.0 / 7.0) * SIGMA_MJ * tk ** 4          # MJ m-2 day-1, ea in hPa
        cfac = (strd / 1e6) / rl_clear(t_ref, saturation_vapour_pressure(tn_ref - k)).mean()   # cloud enhancement, from 2019
        rl_j = cfac * rl_clear(t, ea) * 1e6
    rn = net_radiation_from_era5(rs_j, rl_j, t_for_rn)
    u2 = wind_speed_2m(wind10, z_m=10.0) * np.ones(n_days)
    s = c["soil"]
    return topoclimate_forcing(
        doy=doy, month=month, t_mean_ref_c=t, t_max_ref_c=tx, t_min_ref_c=tn, p_ref_mm=p, ea_ref_kpa=ea,
        u2_m_s=u2, rn_mj_m2=rn, z_cell_m=c["elevation"][i], z_ref_m=c["z_ref_m"][i], lat_deg=lat,
        slope_deg=c["slope"][i], aspect_deg=c["aspect"][i], gamma_k_per_m=gamma_k_per_m,
        precip_gradient_per_m=precip_gradient_per_m, w_max_mm=w_max[i], theta_sat=s["theta_sat"][i],
        psi_sat_mpa=s["psi_sat_mpa"][i], b_clapp_hornberger=s["b_clapp_hornberger"][i], theta_fc=s["theta_fc"][i],
        theta_lim=s["theta_lim"][i], gdd_budburst=200.0, concavity_index=c["concavity"][i], calm_clear_night_frac=0.3,
        pressure_kpa=pressure_pa / 1000.0)


VARIANTS = ("V0", "V1")


def main():
    global VARIANTS
    if "--variants" in sys.argv:
        VARIANTS = tuple(sys.argv[sys.argv.index("--variants") + 1].split(","))
    if "--periods" in sys.argv:
        PERIODS = sys.argv[sys.argv.index("--periods") + 1].split(",")
    else:
        PERIODS = ["2019", "ssp585_2100"]
    ctx_path = sys.argv[sys.argv.index("--ctx") + 1]
    out_path = Path(sys.argv[sys.argv.index("--out") + 1])
    c = pickle.load(open(ctx_path, "rb"))
    cells = [i for i in range(c["n"]) if c["valid_mask"][i] and bool(inside_armenia(c["lats"][i], c["lons"][i]))]
    print(f"{len(cells)} Armenian cells", flush=True)
    deltas = compute_deltas(c["lats"], c["lons"], "gfdl-esm4", "ssp585")
    rows = []
    for i in cells:
        for period, d, hz in [x for x in (("2019", None, None), ("ssp585_2100", deltas, 2100)) if x[0] in PERIODS]:
            for variant in VARIANTS:
                row = {"cell": int(i), "lat": float(c["lats"][i]), "lon": float(c["lons"][i]), "elev": float(c["elevation"][i]), "period": period, "variant": variant, "groups": {}}
                for gname, g in c["groups"].items():
                    cell = forcing(c, i, d, hz, c["w_max_mm_by_group"][gname], variant)
                    psi = cell.psi_soil_mpa[PET_FORMULATION]

                    def simulate(traits, psi=psi, cell=cell):
                        return simulate_two_phase(psi, cell.t_max_c, cell.vpd_24h_kpa, cell.pressure_kpa, traits, 1.0)["hfi_max"]
                    h = two_level_failure_probability(simulate, g["base"], g["hyper_sd"], g["individual_sd"], outer_draws=OUTER_DRAWS, inner_draws=INNER_DRAWS, seed=i)
                    row["groups"][gname] = {"h": float(h.mean()), "cwd": float(cell.cwd_mm[PET_FORMULATION]), "psi_min": float(np.min(psi)),
                                            "vpd_max_p95": float(np.percentile(cell.vpd_max_kpa, 95)), "vpd24_summer": float(np.mean(cell.vpd_24h_kpa[(c["month"] >= 6) & (c["month"] <= 8)]))}
                rows.append(row)
        print(f"  cell {i} done", flush=True)
    out_path.write_text(json.dumps(rows))
    print("wrote", out_path)


if __name__ == "__main__":
    main()
