"""
=============================================================================
  Herd-Specific Lactation Curve Estimation
  A module of the RuFaS (Ruminant Farm Systems) model
=============================================================================

  This script demonstrates how RuFaS predicts the daily milk yield of each
  cow group (by parity) over a 305-day lactation period using Wood's curve.

  Two methods are compared:
    - Base Method   : Uses national-average parameters from published data.
    - Cal Method    : Calibrates the curve to YOUR herd's actual annual milk
                      production (AHMP), parity structure, and herd size.

  The Cal Method is the key contribution of Gong et al. (2025, J. Dairy Sci.)
  It reduces prediction error from ~40% (Base) down to ~2% (Cal).

  Reference:
    Gong Y., et al. (2025). Herd-specific lactation curve estimation for
    the Ruminant Farm Systems model. J. Dairy Sci.
    https://doi.org/10.3168/jds.2024-25809

  Authors: Yijing Gong (gong44@wisc.edu), Haowen Hu (hh598@cornell.edu)
=============================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.integrate import quad

# ─────────────────────────────────────────────────────────────────────────────
# ADJUSTMENT DICTIONARY
# ─────────────────────────────────────────────────────────────────────────────
# These values shift the three Wood's curve parameters (a, b, c) based on
# herd characteristics.  Source: Li et al. (2022), J. Dairy Sci. 105:7525.
#
# Each entry is [Δa, Δb (×10⁻²), Δc (×10⁻⁴)]
# ─────────────────────────────────────────────────────────────────────────────

adjustment_dict = {
    # Parity group
    "1":               [-4.18, -0.37,  -9.31],
    "2":               [ 2.16, -1.20,   2.66],
    "3":               [ 2.02,  1.57,   6.65],
    # Year of calving (2006–2016)
    "2006":            [-0.37,  0.72,   0.83],
    "2007":            [-0.59,  1.00,   1.23],
    "2008":            [-0.31,  0.47,   0.98],
    "2009":            [-0.24,  0.24,   0.60],
    "2010":            [-0.11, -0.14,   0.31],
    "2011":            [ 0.10, -0.58,  -0.56],
    "2012":            [ 0.33, -0.71,  -0.83],
    "2013":            [ 0.27, -0.51,  -0.73],
    "2014":            [ 0.12,  0.069, -0.37],
    "2015":            [ 0.28, -0.12,  -0.68],
    "2016":            [ 0.52, -0.44,  -0.78],
    # Month of calving
    "January":         [-0.46,  1.81,   3.13],
    "February":        [ 0.18,  0.76,   2.43],
    "March":           [ 1.05, -0.77,   1.04],
    "April":           [ 1.58, -2.03,  -0.56],
    "May":             [ 1.49, -2.47,  -1.95],
    "June":            [ 0.74, -2.01,  -2.75],
    "July":            [-0.41, -0.81,  -2.68],
    "August":          [-0.96,  0.11,  -2.06],
    "September":       [-1.08,  0.78,  -1.08],
    "October":         [-0.85,  1.20,   0.27],
    "November":        [-0.63,  1.45,   1.51],
    "December":        [-0.65,  1.98,   2.70],
    # U.S. region
    "Appalachian":     [-0.22, -0.042, -0.89],
    "Corn Belt":       [ 0.55, -0.58,  -1.12],
    "Delta":           [-2.56,  0.59,   1.47],
    "Lake":            [ 0.61, -0.40,  -0.64],
    "Mountain":        [-0.96,  3.13,   1.50],
    "Northeast":       [ 1.04, -1.99,  -1.13],
    "Northern Plains": [-0.26,  0.19,  -0.79],
    "New York":        [ 0.67, -1.21,  -0.45],
    "Pennsylvania":    [ 1.15, -0.96,   0.06],
    "Southeast":       [-2.00,  2.59,   2.60],
    "Southern Plains": [-0.51, -1.02,  -0.93],
    "West Coast":      [ 1.09,  0.53,   0.52],
    "Wisconsin":       [ 1.4,  -0.83,  -0.2 ],
    # Milking frequency
    "2x/d":            [-0.74,  0.090,  0.15],
    "3x/d":            [ 0.74, -0.090, -0.15],
}

# ─────────────────────────────────────────────────────────────────────────────
# CORE FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

def wood_curve(t, a, b, c):
    """
    Wood's lactation curve: y(t) = a * t^b * exp(-c * t)

    Parameters
    ----------
    t : array-like  Days in Milk (DIM), typically 1–305
    a : float       Scale parameter  (controls overall yield level)
    b : float       Rising-phase parameter (controls how fast yield climbs)
    c : float       Declining-phase parameter (controls how fast yield drops)

    Returns
    -------
    Daily milk yield in kg/d
    """
    return a * np.power(t, b) * np.exp(-c * t)


def integrate_wood_curve(a, b, c, dim_start=1, dim_end=305):
    """Numerically integrate Wood's curve from dim_start to dim_end (kg total)."""
    result, _ = quad(wood_curve, dim_start, dim_end, args=(a, b, c))
    return result


def calculate_parity_m305(AHMP, num_milking_cows,
                           parity_structure_pct=(0.386, 0.281, 0.333)):
    """
    Back-calculate the 305-day milk yield (M305) for each parity group.

    The herd-average M305 is split into parity-specific M305 values using
    fixed parity multipliers (1.18× for P2, 1.25× for P3+) from Li et al. (2022).

    Parameters
    ----------
    AHMP                : float  Annual Herd Milk Production (kg/yr)
    num_milking_cows    : int    Number of milking cows in the herd
    parity_structure_pct: tuple  Fraction of cows in parity 1, 2, 3+

    Returns
    -------
    (P1_M305, P2_M305, P3_M305) in kg
    """
    herd_M305 = AHMP * 305 / (365 * num_milking_cows)
    P1_pct, P2_pct, P3_pct = parity_structure_pct
    P1_M305 = herd_M305 / (P1_pct + 1.18 * P2_pct + 1.25 * P3_pct)
    P2_M305 = P1_M305 * 1.18
    P3_M305 = P1_M305 * 1.25
    return P1_M305, P2_M305, P3_M305


def get_base_parameters(parity, year=None, region=None, milking_frequency=None):
    """
    Compute BASE Wood's curve parameters for a parity group.

    Starts from the national-average baseline (a=19.9, b=0.247, c=0.003376)
    and applies adjustments for parity, year, region, and milking frequency.

    Returns
    -------
    (a, b, c, M305_predicted)
    """
    # National-average baseline parameters (Li et al. 2022, Table 2)
    a = 19.9
    b = 24.7  * 1e-2
    c = 33.76 * 1e-4

    for key in [str(parity), year, region, milking_frequency]:
        if key and key in adjustment_dict:
            da, db, dc = adjustment_dict[key]
            a += da
            b += db * 1e-2
            c += dc * 1e-4

    M305 = integrate_wood_curve(a, b, c)
    return a, b, c, M305


def calibrate_parameter_a(target_M305, a_base, b, c, search_range=(-20, 20), step=0.01):
    """
    Calibrate parameter 'a' so that the integrated curve matches target_M305.

    This is the core of the Cal Method: we search for the value of 'a' that
    makes the Wood's curve integrate to exactly the herd's known M305.

    Parameters
    ----------
    target_M305 : float  The desired 305-day yield (kg)
    a_base      : float  Starting point for the search (from Base Method)
    b, c        : float  Fixed shape parameters (kept from Base Method)
    search_range: tuple  Search window for Δa around a_base
    step        : float  Search resolution

    Returns
    -------
    (a_calibrated, M305_achieved)
    """
    best_a    = a_base
    best_M305 = 0
    min_diff  = float('inf')

    for delta_a in np.arange(search_range[0], search_range[1], step):
        a_try    = a_base + delta_a
        M305_try = integrate_wood_curve(a_try, b, c)
        diff = abs(M305_try - target_M305)
        if diff < min_diff:
            min_diff  = diff
            best_a    = a_try
            best_M305 = M305_try

    return best_a, best_M305


def get_calibrated_parameters(parity, target_M305,
                               year=None, region=None, milking_frequency=None):
    """
    Compute CALIBRATED Wood's curve parameters for a parity group.

    Uses the Base Method to get b and c, then searches for the 'a' that
    makes the curve integrate to the herd-specific target_M305.

    Returns
    -------
    (a_cal, b, c, M305_achieved)
    """
    # Get base a, b, c (we only keep b and c from the base)
    a_base, b, c, _ = get_base_parameters(parity, year, region, milking_frequency)

    # Search for the best 'a' around the base value
    best_a, best_M305 = calibrate_parameter_a(target_M305, a_base, b, c,
                                               search_range=(-20, 20), step=0.01)
    return best_a, b, c, best_M305


# ─────────────────────────────────────────────────────────────────────────────
# HIGH-LEVEL DEMO FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

def run_base_method(year, region, milking_frequency, verbose=True):
    """
    BASE METHOD: Predict lactation curves using national-average parameters.

    This is what RuFaS does by default when no herd-specific data is provided.
    It uses adjustments for year, region, and milking frequency, but does NOT
    use the actual herd's milk production records.

    Parameters
    ----------
    year              : str  e.g. '2016'
    region            : str  e.g. 'New York'
    milking_frequency : str  '2x/d' or '3x/d'
    verbose           : bool Print detailed output

    Returns
    -------
    dict with keys 'P1', 'P2', 'P3', each containing (t, y, a, b, c, M305)
    """
    if verbose:
        print()
        print("=" * 65)
        print("  BASE METHOD  —  National-average parameters")
        print("=" * 65)
        print(f"  Herd profile : {region}, {year}, milked {milking_frequency}")
        print("-" * 65)

    results = {}
    t = np.arange(1, 305)

    for parity_label, parity_key in [("Parity 1", "1"),
                                      ("Parity 2", "2"),
                                      ("Parity 3+", "3")]:
        a, b, c, M305 = get_base_parameters(parity_key, year, region, milking_frequency)
        y = wood_curve(t, a, b, c)
        peak_dim = int(b / c)
        peak_yield = wood_curve(peak_dim, a, b, c)

        if verbose:
            print(f"\n  {parity_label}")
            print(f"    Wood's parameters :  a = {a:.4f},  b = {b:.6f},  c = {c:.6f}")
            print(f"    Predicted M305    :  {M305:,.0f} kg")
            print(f"    Peak yield        :  {peak_yield:.1f} kg/d  at DIM {peak_dim}")

        results[parity_key] = {"t": t, "y": y, "a": a, "b": b, "c": c, "M305": M305}

    return results


def run_cal_method(year, region, milking_frequency,
                   AHMP, num_milking_cows, parity_structure_pct,
                   verbose=True):
    """
    CAL METHOD: Calibrate lactation curves to herd-specific production data.

    This is the key innovation: instead of relying on national averages,
    we use three pieces of real farm data to anchor the curve:
      1. Annual Herd Milk Production (AHMP, kg/yr)
      2. Number of milking cows
      3. Parity structure (% of cows in each parity group)

    Parameters
    ----------
    year                 : str    e.g. '2016'
    region               : str    e.g. 'New York'
    milking_frequency    : str    '2x/d' or '3x/d'
    AHMP                 : float  Annual Herd Milk Production (kg/yr)
    num_milking_cows     : int    Number of milking cows
    parity_structure_pct : tuple  (P1_fraction, P2_fraction, P3_fraction)
    verbose              : bool   Print detailed output

    Returns
    -------
    dict with keys 'P1', 'P2', 'P3', each containing (t, y, a, b, c, M305)
    """
    P1_M305, P2_M305, P3_M305 = calculate_parity_m305(
        AHMP, num_milking_cows, parity_structure_pct
    )
    herd_M305 = AHMP * 305 / (365 * num_milking_cows)

    if verbose:
        print()
        print("=" * 65)
        print("  CAL METHOD  —  Herd-specific calibration")
        print("=" * 65)
        print(f"  Herd profile   : {region}, {year}, milked {milking_frequency}")
        print(f"  Herd size      : {num_milking_cows:,} milking cows")
        print(f"  AHMP           : {AHMP:,.0f} kg/yr")
        print(f"  Herd-avg M305  : {herd_M305:,.0f} kg")
        print(f"  Parity split   : P1={parity_structure_pct[0]*100:.1f}%  "
              f"P2={parity_structure_pct[1]*100:.1f}%  "
              f"P3+={parity_structure_pct[2]*100:.1f}%")
        print("-" * 65)
        print(f"  Back-calculated M305 targets:")
        print(f"    Parity 1  →  {P1_M305:,.0f} kg")
        print(f"    Parity 2  →  {P2_M305:,.0f} kg")
        print(f"    Parity 3+ →  {P3_M305:,.0f} kg")
        print("-" * 65)

    results = {}
    t = np.arange(1, 305)
    targets = {"1": P1_M305, "2": P2_M305, "3": P3_M305}
    labels  = {"1": "Parity 1", "2": "Parity 2", "3": "Parity 3+"}

    for parity_key, target_M305 in targets.items():
        a, b, c, M305 = get_calibrated_parameters(
            parity_key, target_M305, year, region, milking_frequency
        )
        y = wood_curve(t, a, b, c)
        peak_dim   = int(b / c)
        peak_yield = wood_curve(peak_dim, a, b, c)

        if verbose:
            print(f"\n  {labels[parity_key]}")
            print(f"    Wood's parameters :  a = {a:.4f},  b = {b:.6f},  c = {c:.6f}")
            print(f"    Calibrated M305   :  {M305:,.0f} kg  (target: {target_M305:,.0f} kg)")
            print(f"    Peak yield        :  {peak_yield:.1f} kg/d  at DIM {peak_dim}")

        results[parity_key] = {"t": t, "y": y, "a": a, "b": b, "c": c, "M305": M305}

    return results


def plot_comparison(base_results, cal_results,
                    herd_label="Example Herd",
                    save_path=None):
    """
    Plot Base vs. Cal lactation curves side by side for all three parities.

    Parameters
    ----------
    base_results : dict  Output from run_base_method()
    cal_results  : dict  Output from run_cal_method()
    herd_label   : str   Title label for the herd
    save_path    : str   If provided, save the figure to this path
    """
    colors = {"1": "#3B82F6", "2": "#10B981", "3": "#F59E0B"}
    parity_labels = {"1": "Parity 1", "2": "Parity 2", "3": "Parity 3+"}

    fig = plt.figure(figsize=(14, 6))
    fig.patch.set_facecolor("white")
    gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.35)

    for col, (method_label, results) in enumerate(
        [("Base Method\n(National-average parameters)", base_results),
         ("Cal Method\n(Herd-specific calibration)",    cal_results)]
    ):
        ax = fig.add_subplot(gs[col])
        for parity_key in ["1", "2", "3"]:
            r = results[parity_key]
            ax.plot(r["t"], r["y"],
                    color=colors[parity_key],
                    linewidth=2.5,
                    label=f"{parity_labels[parity_key]}  (M305 = {r['M305']:,.0f} kg)")

        ax.set_xlabel("Days in Milk (DIM)", fontsize=12)
        ax.set_ylabel("Daily Milk Yield (kg/d)", fontsize=12)
        ax.set_title(method_label, fontsize=13, fontweight="bold", pad=12)
        ax.legend(fontsize=10, framealpha=0.9)
        ax.grid(linestyle=":", alpha=0.6)
        ax.set_xlim(0, 305)
        ax.set_ylim(bottom=0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.suptitle(
        f"Lactation Curve Estimation — {herd_label}",
        fontsize=15, fontweight="bold", y=1.02
    )

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"\n  Figure saved → {save_path}")

    plt.show()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# DEMO  —  Run this script directly to see a worked example
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":

    print()
    print("╔══════════════════════════════════════════════════════════════╗")
    print("║   Herd-Specific Lactation Curve Estimation  (RuFaS module)  ║")
    print("╚══════════════════════════════════════════════════════════════╝")
    print()
    print("  This demo uses a real New York Holstein dairy farm (Farm F1)")
    print("  from Gong et al. (2025, J. Dairy Sci.) to illustrate the")
    print("  difference between the Base and Cal methods.")

    # ── Farm F1 (New York) inputs ──────────────────────────────────────────
    YEAR              = "2016"
    REGION            = "New York"
    MILKING_FREQ      = "3x/d"
    AHMP              = 18_405_332   # kg/yr  (actual farm record)
    NUM_MILKING_COWS  = 1_284
    PARITY_STRUCTURE  = (0.3595, 0.2635, 0.3770)  # P1, P2, P3+

    # ── Step 1: Run Base Method ────────────────────────────────────────────
    base = run_base_method(YEAR, REGION, MILKING_FREQ)

    # ── Step 2: Run Cal Method ─────────────────────────────────────────────
    cal = run_cal_method(YEAR, REGION, MILKING_FREQ,
                         AHMP, NUM_MILKING_COWS, PARITY_STRUCTURE)

    # ── Step 3: Compare M305 predictions ──────────────────────────────────
    print()
    print("=" * 65)
    print("  COMPARISON SUMMARY  —  Predicted M305 (kg per cow per 305 d)")
    print("=" * 65)
    print(f"  {'Parity':<12}  {'Base Method':>16}  {'Cal Method':>16}")
    print(f"  {'-'*12}  {'-'*16}  {'-'*16}")
    for pk, plabel in [("1","Parity 1"), ("2","Parity 2"), ("3","Parity 3+")]:
        print(f"  {plabel:<12}  {base[pk]['M305']:>14,.0f} kg"
              f"  {cal[pk]['M305']:>14,.0f} kg")
    print()
    print("  NOTE: The Cal Method anchors each parity curve to the herd's")
    print("  actual annual production record, making it far more accurate")
    print("  for farms that deviate from the national average.")
    print()

    # ── Step 4: Plot ───────────────────────────────────────────────────────
    print("  Generating comparison plot …")
    plot_comparison(base, cal,
                    herd_label="Farm F1 — New York Holstein Dairy (2016)",
                    save_path="lactation_curve_comparison.png")
