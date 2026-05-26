# gaucho_astrometry

Produces a publication/presentation-quality PDF comparing **ground-based radial-velocity (RV)** and **Gaia DR4 astrometric** planet detection limits on a single log–log mass vs. orbital-period diagram.
All physical parameters and plot settings live in a single YAML file — no code edits required to add new stars or adjust thresholds.

---

## What the code does

For each star listed in the YAML the script draws two curves on the
*planet mass (M_Jup) vs. orbital period (yr)* plane:

| Line style | Method | Sensitivity drivers |
|---|---|---|
| **Solid** | Ground-based RV | RV single-measurement precision σ_RV; stellar mass |
| **Dashed** | Gaia DR4 astrometry | Single-epoch along-scan precision; number of Gaia transits; stellar mass; distance |

Everything **above** a curve is detectable by that method at the chosen S/N threshold.

### Physical formulas

**RV** — minimum detectable projected mass (circular orbit):

$$M_p \sin i \;[\text{M}_{\rm Jup}] \;=\; \frac{N_\sigma\,\sigma_{\rm RV}}{28.43\;\text{m/s}} \;\times\; \left(\frac{M_\star}{M_\odot}\right)^{1/2} \;\times\; \left(\frac{P}{1\;\text{yr}}\right)^{1/3}$$

**Gaia astrometry** — astrometric semi-amplitude and minimum detectable mass:

$$\alpha\;[\mu\text{as}] \;=\; \frac{M_p}{M_\star} \;\times\; \frac{a_p\;[\text{AU}]}{d\;[\text{pc}]} \;\times\; 10^6, \qquad a_p = \left(M_\star\,P^2\right)^{1/3}$$

$$M_{p,\min} \;=\; \frac{N_\sigma\,\sigma_{\rm eff}}{\alpha / M_p}, \qquad \sigma_{\rm eff} = \frac{\sigma_{\rm epoch}}{\sqrt{N_{\rm transits}/2}}$$

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/<your-username>/gaucho_astrometry.git
cd gaucho_astrometry
```

### 2. Create a virtual environment (recommended)

```bash
python -m venv .venv
source .venv/bin/activate   # macOS / Linux
# .venv\Scripts\activate    # Windows
```

### 3. Install dependencies

The script only requires standard scientific Python packages:

```bash
pip install numpy matplotlib pyyaml
```

There is no `setup.py` or package to install — the script runs directly.

---

## Usage

```bash
python detection_limits.py [config.yaml]
```

If no argument is given the script looks for `stars.yaml` in the current working directory.

The output PDF path is set by `output_file` in the YAML (default: `detection_limits.pdf`).

---

## Configuration file (`stars.yaml`)

All parameters are documented inline in `stars.yaml`.  The key sections are:

```yaml
# ── Detection thresholds ──────────────────────────────────────────────────────
snr_threshold: 3.0                    # sigma level required for a detection

# ── Gaia DR4 instrument parameters ───────────────────────────────────────────
gaia_dr4_epoch_precision_muas: 7.0   # single-epoch AL precision [µas] for G~10-12
gaia_dr4_n_transits: 420             # total transits over DR4 baseline (~10.5 yr)

# ── Plot axes ─────────────────────────────────────────────────────────────────
period_range_yr: [0.005, 30.0]       # x-axis limits [yr]
mass_range_mjup: [1.0e-4, 30.0]     # y-axis limits [M_Jup]

# ── Output ───────────────────────────────────────────────────────────────────
output_file: detection_limits.pdf

# ── Stars ────────────────────────────────────────────────────────────────────
stars:
  - name: "Barnard's Star"
    distance_pc:    1.83    # distance in parsecs
    mass_msun:      0.16    # stellar mass in solar masses
    rv_accuracy_ms: 1.5     # RV single-measurement precision [m/s]
```

Add or remove entries under `stars:` to include any target.
There is no hard limit on the number of stars; colours cycle automatically.

### Gaia DR4 precision note

The default `gaia_dr4_epoch_precision_muas: 7.0` µas corresponds to the nominal
Gaia end-of-mission single-epoch along-scan precision for bright stars (G ≈ 10–12).
For fainter targets increase this value (e.g. ~24 µas at G ≈ 14).
References: Gaia Collaboration (2016, A&A 595 A1); de Bruijne (2012, Ap&SS 341 31).

---

## Output

A single PDF (`detection_limits.pdf` by default) containing:

- Log–log mass vs. period axes with dual y-axis (M_Jup left, M_⊕ right).
- One solid + one dashed curve per star (RV and Gaia DR4 respectively).
- Horizontal reference lines for Earth, Neptune, Saturn, and Jupiter.
- Two legend boxes: detection method (line style) and star list (colour + key parameters).
- Sub-period axis labels in days for P < 1 yr, years for P ≥ 1 yr.

---

## File overview

```
gaucho_astrometry/
├── detection_limits.py   # main script
├── stars.yaml            # input config (edit this)
└── README.md             # this file
```

---

## License

MIT — see `LICENSE` if present, otherwise use freely with attribution.
