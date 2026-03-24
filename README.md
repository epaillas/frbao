# FRBAO: Fast Radio Bursts meet Baryon Acoustic Oscillations

Minimal simulation framework to validate **galaxy × FRB-DM** clustering measurements in a controlled periodic box.

This repository is the first step toward testing whether the **BAO feature can be detected in the free-electron field** using **FRB dispersion measures (DMs)** cross-correlated with galaxy large-scale structure. The current implementation is intentionally simple and focuses on a baseline **validation stage**.

## Installation

```bash
pip install -e /path/to/frbao/
```

Dependencies (`cosmoprimo`, `mockfactory`) are from the [cosmodesi](https://github.com/cosmodesi) ecosystem and must be installed separately. At NERSC, source the cosmodesi environment first:

```bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
```

## Current goal

A sanity-check layer designed to answer one question:

> Can we reliably recover the smooth galaxy × DM cross-spectrum and estimate its covariance in a clean, reproducible mock setup?

## What is modeled

- Periodic 3D matter fields generated from  **wiggle** and **no-wiggle** input power spectra
- Foreground galaxy catalogs drawn from those fields
- A continuous free-electron field constructed from the same realization
- A dense FRB catalog with source positions in the box
- A mock FRB observable
  \[
  \Delta \mathrm{DM}(x,y,z_s) \propto \int_0^{z_s} \delta_e(x,y,z)\,dz
  \]
  using a fixed Cartesian line of sight along `+z`
- Measurement of the galaxy auto-spectrum and galaxy × DM cross-spectrum with `pypower`
- Covariance estimation from multiple realizations

## What is not modeled yet

- BAO fitting
- Host-galaxy DM
- Milky Way DM
- Instrumental noise
- FRB redshift evolution
- Redshift-space distortions
- Cut-sky geometry or survey masks
- Selection effects
- Lightcone structure
- Realistic baryonic gas physics

## Minimal model

The current implementation uses:

- **Galaxy field**
  \[
  \delta_g = b_g\,\delta_m
  \]

- **Electron field**
  \[
  \delta_e(\mathbf{k}) = b_e \exp\!\left[-\frac{(kR_e)^2}{2}\right]\delta_m(\mathbf{k})
  \]

- **FRB observable**
  \[
  \Delta \mathrm{DM} \propto \int \delta_e\,dz
  \]

This is a controlled estimator-validation model, not a full astrophysical description.
