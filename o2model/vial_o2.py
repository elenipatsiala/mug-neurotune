"""Oxygen supply to a single organoid in a sealed cryovial: 1 g (ground) vs µg (ISS).

Question: could the LEO organoids of Marotta et al. 2024 (CRS-19) have been more
oxygen-limited than ground controls *in identical vials*, because microgravity removes
buoyancy-driven convection?

Facts used from the paper (Marotta et al. 2024, Stem Cells Transl Med 13:1186):
  * one organoid per Nunc cryovial (Thermo #374088), 1 mL medium, sealed, 37 °C,
    30 days static culture, no medium exchange; organoids ~200-500 µm diameter,
    formed from 1e5 neural precursor cells.
  * authors report "no signs of hypoxia in either group" (method not specified).
Everything else is an ASSUMPTION with a range (see PRIORS) - the model is a
hypothesis generator, not a measurement.

Model, three levels:
  L0  whole-vial O2 budget (sealed, impermeable worst case): days until headspace O2 is used up.
  L1  vial-scale transport from O2 source (gas headspace) to the organoid:
        1 g : organoid sediments to the bottom, source = top interface, distance = medium height;
              thermal convection adds mixing  D_eff = D + c_mix * U * r_vial
        µg  : no sedimentation and no buoyant convection; the gas bubble and the organoid sit at
              random positions, so distance ~ U(a, h); residual mixing from g-jitter/vibration only.
      resistance R_far = L / (D_eff * A_cross)
  L2  organoid scale: external boundary layer + internal reaction-diffusion
        external: Sherwood number Sh = 1 + (1 + Pe)^(1/3)  (creeping flow past a sphere),
                  Pe = U * 2a / D;  pure diffusion -> Sh = 2.
                  in 1 g a fraction of the surface rests on the vial wall (contact factor f).
        internal: spherical diffusion with Michaelis-Menten uptake, solved by finite differences.
Output: fraction of organoid volume below hypoxia thresholds, in each gravity condition, under
matched parameters (paired Monte Carlo) -> probability that µg is MORE hypoxic than 1 g.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

R_GAS = 8.314
T_K = 310.15
ALPHA_MM_PER_MMHG = 1.3e-3          # O2 solubility in medium at 37 °C, mM per mmHg (≈ mol/m³ per mmHg)
MMHG_PER_PCT = 7.13                  # 1 % O2 at 37 °C, humidified, 760 mmHg -> ~7.1 mmHg
C_1PCT = ALPHA_MM_PER_MMHG * MMHG_PER_PCT          # ≈ 0.0093 mol/m³ (1 % O2): strong HIF activation
C_ANOXIC = ALPHA_MM_PER_MMHG * MMHG_PER_PCT * 0.1  # 0.1 % O2: anoxia / necrosis risk

# ----------------------------------------------------------------------------- priors
# (low, high, scale) ; 'log' = log-uniform. Ranges are deliberately wide.
PRIORS = {
    "radius_um":        (100, 250, "lin"),     # paper: 200-500 µm diameter
    "n_cells":          (5e4, 2e5, "log"),     # 1e5 at formation; proliferation/death over 30 d
    "ocr_amol_s":       (5, 100, "log"),       # O2 consumption per cell, amol/cell/s (neural cells in vitro)
    "km_uM":            (0.5, 5, "log"),       # Michaelis-Menten constant for O2
    "D_med":            (2.6e-9, 3.3e-9, "lin"),  # m²/s, O2 in medium at 37 °C
    "D_tissue_frac":    (0.5, 0.9, "lin"),     # D_tissue = frac * D_med
    "C_far_frac":       (0.4, 1.0, "lin"),     # source O2 as fraction of air-saturation (headspace slowly depletes)
    "U_1g_um_s":        (0.5, 50, "log"),      # thermal-convection velocity in an incubated vial, 1 g
    "U_ug_um_s":        (0.001, 0.5, "log"),   # residual mixing on ISS (quasi-steady ~1e-6 g + vibration)
    "contact_1g":       (0.7, 0.95, "lin"),    # fraction of organoid surface exposed when resting on wall
    "vial_radius_mm":   (4.0, 5.5, "lin"),     # inner radius of a ~1.8-2 mL cryovial (assumed)
    "c_mix":            (0.05, 0.3, "lin"),    # mixing-length factor for D_eff
}
MEDIUM_ML = 1.0
C_SAT_AIR = ALPHA_MM_PER_MMHG * 0.21 * (760 - 47)   # ≈ 0.195 mol/m³ at 21 % O2, 37 °C, humidified


def sample_priors(n: int, seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    out = {}
    for k, (lo, hi, sc) in PRIORS.items():
        u = rng.random(n)
        out[k] = np.exp(np.log(lo) + u * (np.log(hi) - np.log(lo))) if sc == "log" else lo + u * (hi - lo)
    df = pd.DataFrame(out)
    df["u_pos"] = rng.random(n)   # µg: random position of organoid relative to bubble
    return df


def solve_sphere(a: float, Vmax: float, Km: float, D_t: float, C_out_fn, n: int = 160,
                 iters: int = 200, tol: float = 1e-10) -> tuple[np.ndarray, np.ndarray, float]:
    """Steady spherical reaction-diffusion with a Robin outer boundary.

    D_t (1/r²) d/dr(r² dC/dr) = Vmax C/(Km + C),  dC/dr(0)=0,
    total uptake Q = (C_far - C(a)) / R_ext  supplied through the surface.
    C_out_fn is (C_far, R_ext). Newton iteration on a cell-centred finite-volume grid.
    Returns radii, concentration profile, total uptake Q (mol/s).
    """
    C_far, R_ext = C_out_fn
    dr = a / n
    rc = (np.arange(n) + 0.5) * dr            # cell centres
    rf = np.arange(n + 1) * dr                # faces
    vol = 4 / 3 * np.pi * (rf[1:] ** 3 - rf[:-1] ** 3)
    area = 4 * np.pi * rf ** 2
    C = np.full(n, C_far * 0.5)
    # outer face conductance in series: half-cell diffusion + external resistance
    g_half = D_t * area[-1] / (dr / 2)
    g_out = 1.0 / (1.0 / g_half + R_ext)
    g_int = D_t * area[1:-1] / dr               # between cells i and i+1
    for _ in range(iters):
        upt = Vmax * C / (Km + C)
        dupt = Vmax * Km / (Km + C) ** 2
        F = -upt * vol
        F[:-1] += g_int * (C[1:] - C[:-1])
        F[1:] -= g_int * (C[1:] - C[:-1])
        F[-1] += g_out * (C_far - C[-1])
        # Jacobian (tridiagonal)
        main = -dupt * vol
        main[:-1] -= g_int
        main[1:] -= g_int
        main[-1] -= g_out
        J = np.diag(main) + np.diag(g_int, 1) + np.diag(g_int, -1)
        step = np.linalg.solve(J, -F)
        Cn = np.maximum(C + step, 0.0)
        if np.max(np.abs(Cn - C)) < tol:
            C = Cn
            break
        C = Cn
    Q = g_out * (C_far - C[-1])
    return rc, C, Q


def run_condition(p: pd.Series, gravity: str) -> dict:
    a = p.radius_um * 1e-6
    vol_org = 4 / 3 * np.pi * a ** 3
    q_cell = p.ocr_amol_s * 1e-18
    Vmax = q_cell * p.n_cells / vol_org          # mol/m³/s
    Km = p.km_uM * 1e-3                           # mol/m³
    D_t = p.D_med * p.D_tissue_frac
    r_v = p.vial_radius_mm * 1e-3
    A_cross = np.pi * r_v ** 2
    h = MEDIUM_ML * 1e-6 / A_cross                 # medium column height
    if gravity == "1g":
        U = p.U_1g_um_s * 1e-6
        L = h                                      # organoid at bottom, gas at top
        contact = p.contact_1g
    else:
        U = p.U_ug_um_s * 1e-6
        L = a + p.u_pos * (h - a)                  # bubble anywhere
        contact = 1.0
    D_eff = p.D_med + p.c_mix * U * r_v
    R_far = L / (D_eff * A_cross)
    Pe = U * 2 * a / p.D_med
    Sh = 1 + (1 + Pe) ** (1 / 3)
    k = Sh * p.D_med / (2 * a)
    R_near = 1.0 / (k * 4 * np.pi * a ** 2 * contact)
    C_far = p.C_far_frac * C_SAT_AIR
    r, C, Q = solve_sphere(a, Vmax, Km, D_t, (C_far, R_far + R_near))
    shell_vol = 4 / 3 * np.pi * (np.append(r[1:] - (r[1] - r[0]) / 2, a) ** 3
                                 - np.append(0, r[1:] - (r[1] - r[0]) / 2) ** 3)
    frac = lambda thr: float(shell_vol[C < thr].sum() / shell_vol.sum())  # noqa: E731
    return {"surface_pct_o2": C[-1] / C_1PCT, "core_pct_o2": C[0] / C_1PCT,
            "hypoxic_frac_1pct": frac(C_1PCT), "anoxic_frac": frac(C_ANOXIC),
            "uptake_mol_s": Q, "demand_mol_s": q_cell * p.n_cells,
            "R_far": R_far, "R_near": R_near, "Sh": Sh, "Pe": Pe, "L_mm": L * 1e3}


def vial_budget_days(p: pd.Series, vial_total_ml: float = 1.8) -> float:
    """L0: days until a sealed, gas-impermeable vial runs out of O2 (headspace + dissolved)."""
    head_m3 = (vial_total_ml - MEDIUM_ML) * 1e-6
    n_head = 0.21 * 101325 * head_m3 / (R_GAS * T_K)
    n_dis = C_SAT_AIR * MEDIUM_ML * 1e-6
    demand = p.ocr_amol_s * 1e-18 * p.n_cells
    return (n_head + n_dis) / demand / 86400


def monte_carlo(n: int = 3000, seed: int = 0) -> pd.DataFrame:
    P = sample_priors(n, seed)
    rows = []
    for i, p in P.iterrows():
        g1 = run_condition(p, "1g")
        ug = run_condition(p, "ug")
        row = p.to_dict()
        row.update({f"{k}_1g": v for k, v in g1.items()})
        row.update({f"{k}_ug": v for k, v in ug.items()})
        row["budget_days"] = vial_budget_days(p)
        rows.append(row)
    df = pd.DataFrame(rows)
    df["delta_hypoxic"] = df.hypoxic_frac_1pct_ug - df.hypoxic_frac_1pct_1g
    df["delta_surface_pct"] = df.surface_pct_o2_ug - df.surface_pct_o2_1g
    return df
