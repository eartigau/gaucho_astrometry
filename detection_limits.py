#!/usr/bin/env python3
"""
detection_limits.py
───────────────────
Plot the minimum detectable planet mass as a function of orbital period for:
  • Ground-based radial velocity (RV) — one curve per star, using the
    per-star RV precision defined in the YAML.
  • Gaia DR4 astrometry — one curve per star, using the nominal
    single-epoch along-scan precision and number of transits from the YAML.

Usage
-----
    python detection_limits.py [stars.yaml]

If no argument is given the script looks for ``stars.yaml`` in the
current working directory.
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import yaml

# ── physical constants ────────────────────────────────────────────────────────
MJUP_TO_MSUN   = 9.54588e-4   # 1 M_Jup  in M_sun
MEARTH_TO_MSUN = 3.00348e-6   # 1 M_Earth in M_sun
MEARTH_TO_MJUP = MEARTH_TO_MSUN / MJUP_TO_MSUN   # ~1/317.8
MNEP_TO_MJUP   = 17.147 * MEARTH_TO_MJUP
MSAT_TO_MJUP   = 0.2994        # Saturn / Jupiter

# RV semi-amplitude coefficient [m/s] for a circular orbit:
#   K = RV_K0 * (Mp/MJup) * (Mstar/Msun)^{-1/2} * (P/yr)^{-1/3}
RV_K0 = 28.4329  # m/s


# ── physics ───────────────────────────────────────────────────────────────────

def rv_min_mass(period_yr, mstar_msun, sigma_rv_ms, snr):
    """Minimum detectable M_p sin i [M_Jup] via RV (circular orbit)."""
    k_min = snr * sigma_rv_ms
    return (k_min / RV_K0) * mstar_msun**0.5 * period_yr**(1.0 / 3.0)


def gaia_min_mass(period_yr, mstar_msun, dist_pc,
                  epoch_precision_muas, n_transits, snr):
    """
    Minimum detectable planet mass [M_Jup] via Gaia astrometry.

    The astrometric semi-amplitude is

        α [µas] = (Mp/M★) × (a_p [AU] / d [pc]) × 1e6

    with  a_p = (M★ [M☉] × P [yr]²)^{1/3}  from Kepler's third law.

    The effective astrometric noise on a fitted orbital solution scales as
    σ_eff ≈ epoch_precision / √(n_transits / 2), following Lattanzi et al.
    and Perryman et al. (2014) [the /2 accounts for the two astrometric
    parameters needed per baseline projection in the simplified model].
    """
    a_au = (mstar_msun * period_yr**2) ** (1.0 / 3.0)
    alpha_per_mjup = (MJUP_TO_MSUN / mstar_msun) * (a_au / dist_pc) * 1e6  # µas/MJup

    sigma_eff_muas = epoch_precision_muas / np.sqrt(n_transits / 2.0)
    alpha_min = snr * sigma_eff_muas
    return alpha_min / alpha_per_mjup  # M_Jup


# ── plot helpers ──────────────────────────────────────────────────────────────

# Reference planet masses [M_Jup] and labels
REF_PLANETS = [
    ("Earth",   MEARTH_TO_MJUP,         "#8B9EB7"),
    ("Neptune", MNEP_TO_MJUP,            "#5A8F7B"),
    ("Saturn",  MSAT_TO_MJUP,            "#C09A50"),
    ("Jupiter", 1.0,                     "#B05A2F"),
]


def _add_reference_lines(ax, period_range, mass_range):
    """Draw horizontal reference lines for solar-system planet masses."""
    for name, mass_mjup, color in REF_PLANETS:
        if mass_range[0] <= mass_mjup <= mass_range[1]:
            ax.axhline(mass_mjup, color=color, lw=0.9, ls=":", alpha=0.75, zorder=1)
            ax.text(
                period_range[1] * 0.97, mass_mjup * 1.06,
                name, color=color, fontsize=8.5,
                ha="right", va="bottom", style="italic", zorder=5,
            )


def _secondary_yaxis_labels(ax, mass_range):
    """Add a right-hand y-axis with Earth-mass labels at the same ticks."""
    ax2 = ax.twinx()
    ax2.set_yscale("log")
    ax2.set_ylim([m / MEARTH_TO_MJUP for m in mass_range])
    ax2.set_ylabel(r"Planet mass  ($M_\oplus$)", fontsize=12, labelpad=8)
    ax2.tick_params(labelsize=10)
    return ax2


# ── main ──────────────────────────────────────────────────────────────────────

def make_plot(cfg: dict) -> None:
    stars        = cfg["stars"]
    p_range      = cfg.get("period_range_yr",  [0.02, 30.0])
    m_range      = cfg.get("mass_range_mjup",  [1e-4, 30.0])
    snr          = cfg.get("snr_threshold",    3.0)
    ep_prec      = cfg.get("gaia_dr4_epoch_precision_muas", 7.0)
    n_transits   = cfg.get("gaia_dr4_n_transits", 420)
    out_file     = cfg.get("output_file", "detection_limits.pdf")

    periods = np.logspace(np.log10(p_range[0]), np.log10(p_range[1]), 1000)

    # ── figure style ──────────────────────────────────────────────────────────
    mpl.rcParams.update({
        "font.family":      "sans-serif",
        "font.sans-serif":  ["Helvetica Neue", "Arial", "DejaVu Sans"],
        "axes.linewidth":   1.2,
        "xtick.direction":  "in",
        "ytick.direction":  "in",
        "xtick.major.size": 6,
        "ytick.major.size": 6,
        "xtick.minor.size": 3,
        "ytick.minor.size": 3,
        "xtick.top":        True,
        "ytick.right":      False,  # secondary axis takes the right side
    })

    fig, ax = plt.subplots(figsize=(11, 7))

    # ── color palette (colorblind-friendly) ───────────────────────────────────
    colors = [
        "#E63946", "#457B9D", "#2A9D8F", "#E9C46A", "#F4A261",
        "#6A0572", "#1D3557", "#A8DADC",
    ]
    colors = (colors * 4)[: len(stars)]   # repeat if more stars than colors

    # ── detection curves ──────────────────────────────────────────────────────
    for star, color in zip(stars, colors):
        name   = star["name"]
        d_pc   = star["distance_pc"]
        m_star = star["mass_msun"]
        rv_acc = star["rv_accuracy_ms"]

        rv_mass   = rv_min_mass(periods, m_star, rv_acc,    snr)
        gaia_mass = gaia_min_mass(periods, m_star, d_pc,
                                  ep_prec, n_transits, snr)

        ax.plot(periods, rv_mass,   color=color, lw=2.2, ls="-",  zorder=4)
        ax.plot(periods, gaia_mass, color=color, lw=2.2, ls="--", zorder=4)

        # Label the RV curve at mid-period to avoid legend crowding
        idx = len(periods) // 3
        ax.annotate(
            name,
            xy=(periods[idx], rv_mass[idx]),
            xytext=(0, 5), textcoords="offset points",
            color=color, fontsize=8.5, fontweight="bold",
            ha="center", va="bottom", zorder=6,
        )

    # ── reference planet lines ────────────────────────────────────────────────
    _add_reference_lines(ax, p_range, m_range)

    # ── axes ──────────────────────────────────────────────────────────────────
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(p_range)
    ax.set_ylim(m_range)

    # x-tick labels in years with nice formatting
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(
        lambda x, _: f"{x:g} yr" if x >= 1 else f"{x*365.25:.0f} d"
    ))
    ax.xaxis.set_minor_locator(mpl.ticker.LogLocator(subs="all", numticks=100))
    ax.yaxis.set_minor_locator(mpl.ticker.LogLocator(subs="all", numticks=100))

    ax.set_xlabel("Orbital period", fontsize=13, labelpad=6)
    ax.set_ylabel(r"Planet mass  ($M_{\rm Jup}$)", fontsize=13, labelpad=6)
    ax.set_title(
        "Planet detection limits: ground-based RV vs. Gaia DR4 astrometry\n"
        r"(solid = RV,  dashed = Gaia DR4 astrometry,  $"
        + f"{snr:.0f}$" + r"\,\sigma$ threshold)",
        fontsize=12, pad=10,
    )

    ax.grid(True, which="major", lw=0.5, alpha=0.35, color="gray")
    ax.grid(True, which="minor", lw=0.25, alpha=0.18, color="gray")

    # ── secondary y-axis (Earth masses) ───────────────────────────────────────
    _secondary_yaxis_labels(ax, m_range)

    # ── legend ────────────────────────────────────────────────────────────────
    style_handles = [
        Line2D([0], [0], color="k", lw=2.2, ls="-",
               label=f"Ground-based RV  ({snr:.0f}σ)"),
        Line2D([0], [0], color="k", lw=2.2, ls="--",
               label=f"Gaia DR4 astrometry  ({snr:.0f}σ,  {ep_prec:.0f} µas/epoch)"),
    ]
    star_handles = [
        mpatches.Patch(facecolor=c, label=s["name"] +
                       f"  ({s['distance_pc']:.2f} pc,  {s['rv_accuracy_ms']:.1f} m/s)")
        for s, c in zip(stars, colors)
    ]

    leg1 = ax.legend(
        handles=style_handles, fontsize=9.5,
        loc="upper left", framealpha=0.9, edgecolor="0.7",
        title="Detection method", title_fontsize=9,
    )
    ax.legend(
        handles=star_handles, fontsize=9,
        loc="lower right", framealpha=0.9, edgecolor="0.7",
        title="Stars (dist., RV acc.)", title_fontsize=9,
    )
    ax.add_artist(leg1)

    # ── save ──────────────────────────────────────────────────────────────────
    fig.tight_layout()
    fig.savefig(out_file, dpi=200, bbox_inches="tight")
    print(f"Saved → {out_file}")
    plt.close(fig)


if __name__ == "__main__":
    cfg_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("stars.yaml")
    if not cfg_path.exists():
        sys.exit(f"Config file not found: {cfg_path}")
    with cfg_path.open() as fh:
        cfg = yaml.safe_load(fh)
    make_plot(cfg)
