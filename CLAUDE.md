# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

**FRBAO** is a minimal simulation framework for validating **galaxy × FRB dispersion measure (DM)** clustering measurements in a controlled periodic box. The scientific goal is to test whether the BAO feature can be detected in the free-electron field using FRB DMs cross-correlated with galaxy large-scale structure.

Current scope: a validation layer to answer whether the smooth galaxy × DM cross-spectrum can be reliably recovered and its covariance estimated in a clean, reproducible mock setup.

## Minimal Model

- **Galaxy field:** δ_g = b_g δ_m
- **Electron field:** δ_e(k) = b_e exp[-(kR_e)²/2] δ_m(k)  (Gaussian-smoothed, same underlying realization as galaxies)
- **FRB observable:** ΔDM(x,y,z_s) ∝ ∫₀^z_s δ_e(x,y,z) dz  (line-of-sight projection along +z)

Fields are generated in periodic 3D boxes from wiggle and no-wiggle input power spectra. Spectra are measured with `pypower`. Covariance is estimated across multiple realizations.

## What is explicitly NOT modeled

BAO fitting, host/MW DM contributions, instrumental noise, FRB redshift evolution, RSDs, cut-sky geometry, selection effects, lightcone structure, baryonic gas physics.

## Key Dependencies

This package lives alongside other cosmodesi packages at `/Users/epaillas/code/`:

| Dependency | Role |
|---|---|
| `pypower` | FFT-based power spectrum / cross-spectrum measurement |
| `cosmoprimo` | Input power spectra (wiggle/no-wiggle via CLASS/CAMB) |
| `pycorr` | Two-point correlation functions (if needed) |

## Development Setup

```bash
# Install in editable mode (once setup.py/pyproject.toml exists)
pip install -e /Users/epaillas/code/frbao/

# DESI/cosmodesi environment at NERSC
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
```

## Running Tests

```bash
pytest /Users/epaillas/code/frbao/tests/
pytest /Users/epaillas/code/frbao/tests/test_<module>.py::test_name
```
