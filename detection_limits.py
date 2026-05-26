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
MEARTH_TO_MJUP = MEARTH_TO_MSUN / MJUP_TO_MSUN   # ~1/317.8  (used for axis conversion)
MNEP_TO_MJUP   = 17.147 * MEARTH_TO_MJUP          # Neptune = 17.147 M_Earth
MSAT_TO_MJUP   = 0.2994                            # Saturn  = 0.2994 M_Jup

# Prefactor in the RV semi-amplitude formula for a circular orbit:
#   K [m/s] = RV_K0 × (Mp sin i / MJup) × (M★/M☉)^{-1/2} × (P/yr)^{-1/3}
# Derived from Kepler's 3rd law + Newton; see e.g. Butler et al. (2006).
RV_K0 = 28.4329  # m/s


# ── physics ───────────────────────────────────────────────────────────────────

def rv_min_mass(period_yr, mstar_msun, sigma_rv_ms, snr):
    """
    Minimum detectable M_p sin i [M_Jup] via radial velocity (circular orbit).

    Inverts the RV semi-amplitude formula:
        K = RV_K0 × (Mp sin i / MJup) × (M★/M☉)^{-1/2} × (P/yr)^{-1/3}
    Setting K = snr × σ_RV gives the minimum detectable projected mass.

    Parameters
    ----------
    period_yr   : float or array  – orbital period in years
    mstar_msun  : float           – stellar mass in solar masses
    sigma_rv_ms : float           – single-measurement RV precision in m/s
    snr         : float           – detection S/N threshold (e.g. 3)

    Returns
    -------
    float or array : minimum M_p sin i in Jupiter masses
    """
    # Detection requires K ≥ snr × σ_RV
    k_min = snr * sigma_rv_ms
    # Invert K formula: Mp = K/RV_K0 × M★^{1/2} × P^{1/3}
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

    Parameters
    ----------
    period_yr            : float or array – orbital period in years
    mstar_msun           : float          – stellar mass in solar masses
    dist_pc              : float          – distance to the star in parsecs
    epoch_precision_muas : float          – Gaia single-epoch AL precision [µas]
    n_transits           : int            – total number of Gaia transits (epochs)
    snr                  : float          – detection S/N threshold

    Returns
    -------
    float or array : minimum detectable planet mass in Jupiter masses
    """
    # Planet semi-major axis [AU] from Kepler's 3rd law: a³ = M★ P²  (SI-free form)
    a_au = (mstar_msun * period_yr**2) ** (1.0 / 3.0)

    # Astrometric signal per Jupiter mass [µas]:
    #   α = (Mp/M★) × (a/d)  →  convert mass ratio to M_Jup, distance to arcsec
    alpha_per_mjup = (MJUP_TO_MSUN / mstar_msun) * (a_au / dist_pc) * 1e6  # µas/MJup

    # Effective noise on the fitted astrometric orbit (averages down with √N);
    # the factor /2 reflects that two free parameters absorb half the transits
    sigma_eff_muas = epoch_precision_muas / np.sqrt(n_transits / 2.0)

    # Minimum detectable signal, then convert to planet mass
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
        # Only draw the line if it falls inside the current y-axis range
        if mass_range[0] <= mass_mjup <= mass_range[1]:
            ax.axhline(mass_mjup, color=color, lw=0.9, ls=":", alpha=0.75, zorder=1)
            # Place the label just above the line, flush with the right edge
            ax.text(
                period_range[1] * 0.97, mass_mjup * 1.06,
                name, color=color, fontsize=8.5,
                ha="right", va="bottom", style="italic", zorder=5,
            )


def _secondary_yaxis_labels(ax, mass_range):
    """
    Add a right-hand y-axis whose scale is identical to the left axis but
    labelled in Earth masses instead of Jupiter masses.
    The twin axis shares the same y-coordinate system; we just rescale the
    displayed values using the M_Jup → M_Earth conversion factor.
    """
    ax2 = ax.twinx()
    ax2.set_yscale("log")
    # Convert the M_Jup limits to M_Earth for the right-side axis
    ax2.set_ylim([m / MEARTH_TO_MJUP for m in mass_range])
    ax2.set_ylabel(r"Planet mass  ($M_\oplus$)", fontsize=12, labelpad=8)
    ax2.tick_params(labelsize=10)
    return ax2


# ── main ──────────────────────────────────────────────────────────────────────

def make_plot(cfg: dict) -> None:
    """Build and save the detection-limit plot from a parsed YAML config dict."""

    # ── unpack configuration with sensible defaults ────────────────────────────
    stars        = cfg["stars"]
    p_range      = cfg.get("period_range_yr",  [0.02, 30.0])  # years
    m_range      = cfg.get("mass_range_mjup",  [1e-4, 30.0])  # Jupiter masses
    snr          = cfg.get("snr_threshold",    3.0)            # sigma level
    ep_prec      = cfg.get("gaia_dr4_epoch_precision_muas", 7.0)   # µas/epoch
    n_transits   = cfg.get("gaia_dr4_n_transits", 420)         # total Gaia transits
    out_file     = cfg.get("output_file", "detection_limits.pdf")

    # Logarithmically spaced period grid spanning the requested range
    periods = np.logspace(np.log10(p_range[0]), np.log10(p_range[1]), 1000)

    # ── figure style ──────────────────────────────────────────────────────────
    # Apply presentation-friendly rcParams: clean sans-serif font, inward ticks
    # on all four sides (top/bottom/left; right reserved for the twin axis).
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
        "ytick.right":      False,  # right ticks handled by the secondary y-axis
    })

    fig, ax = plt.subplots(figsize=(11, 7))

    # ── color palette (colorblind-friendly) ───────────────────────────────────
    # Eight manually chosen colours with sufficient contrast between adjacent
    # hues; the list is tiled if there are more stars than colours.
    colors = [
        "#E63946", "#457B9D", "#2A9D8F", "#E9C46A", "#F4A261",
        "#6A0572", "#1D3557", "#A8DADC",
    ]
    colors = (colors * 4)[: len(stars)]   # tile to cover any number of stars

    # ── detection curves ──────────────────────────────────────────────────────
    # Each star contributes two curves sharing the same colour:
    #   solid  (-)  = minimum mass detectable by ground RV
    #   dashed (--) = minimum mass detectable by Gaia DR4 astrometry
    # Everything above a curve is detectable by that method.
    for star, color in zip(stars, colors):
        name   = star["name"]
        d_pc   = star["distance_pc"]
        m_star = star["mass_msun"]
        rv_acc = star["rv_accuracy_ms"]

        # Compute detection-limit arrays over the period grid
        rv_mass   = rv_min_mass(periods, m_star, rv_acc, snr)
        gaia_mass = gaia_min_mass(periods, m_star, d_pc,
                                  ep_prec, n_transits, snr)

        ax.plot(periods, rv_mass,   color=color, lw=2.2, ls="-",  zorder=4)
        ax.plot(periods, gaia_mass, color=color, lw=2.2, ls="--", zorder=4)

        # Place the star name along the RV curve at 1/3 of the period range
        # (avoids the crowded short-period end and the right-edge labels)
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

    # Custom x-tick formatter: show periods < 1 yr in days, ≥ 1 yr in years
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(
        lambda x, _: f"{x:g} yr" if x >= 1 else f"{x*365.25:.0f} d"
    ))
    # Dense minor ticks on both axes for easier reading on a log scale
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
    # Two separate legend boxes:
    #   • top-left:     line-style legend (detection method)
    #   • bottom-right: colour legend (one patch per star, with key parameters)
    # leg1 must be added back as an artist after the second legend call,
    # because each ax.legend() call replaces the previous one by default.
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
    ax.add_artist(leg1)  # re-attach the first legend (second call would remove it)

    # ── save ──────────────────────────────────────────────────────────────────
    fig.tight_layout()                               # remove excess whitespace
    fig.savefig(out_file, dpi=200, bbox_inches="tight")  # vector + raster at 200 dpi
    print(f"Saved → {out_file}")
    plt.close(fig)                                   # free memory


if __name__ == "__main__":
    cfg_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("stars.yaml")
    if not cfg_path.exists():
        sys.exit(f"Config file not found: {cfg_path}")
    with cfg_path.open() as fh:
        cfg = yaml.safe_load(fh)
    make_plot(cfg)
