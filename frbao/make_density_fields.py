"""Generate paired Gaussian matter, galaxy, and electron density fields with and without BAO wiggles."""

import numpy as np
from cosmoprimo.fiducial import DESI
from cosmoprimo import PowerSpectrumBAOFilter
from mockfactory import EulerianLinearMock


def get_power_spectra(z=0.5):
    """Return (pk_wiggle, pk_nowiggle) callables at redshift z.

    Parameters
    ----------
    z : float, optional
        Redshift at which to evaluate the power spectrum.

    Returns
    -------
    pk_wiggle : PowerSpectrumInterpolator1D
        Full linear matter power spectrum (with BAO wiggles).
    pk_nowiggle : PowerSpectrumInterpolator1D
        De-wiggled (smooth) power spectrum via Wallish et al. 2018 BAO filter.
    """
    cosmo = DESI()
    pk = cosmo.get_fourier().pk_interpolator().to_1d(z=z)
    filt = PowerSpectrumBAOFilter(pk, engine='wallish2018', cosmo=cosmo)
    pk_nw = filt.smooth_pk_interpolator()
    return pk, pk_nw


def make_matter_field(power, nmesh, boxsize, seed):
    """Generate a Gaussian matter density field on a periodic mesh.

    Parameters
    ----------
    power : callable
        Power spectrum P(k), callable as power(k) where k is in h/Mpc.
    nmesh : int
        Number of mesh cells per side.
    boxsize : float
        Periodic box size in Mpc/h.
    seed : int
        Random seed for the white noise realization.

    Returns
    -------
    mesh_delta_r : RealField
        Real-space matter density field delta_m on the mesh.
    """
    mock = EulerianLinearMock(power, nmesh=nmesh, boxsize=boxsize, seed=seed)
    mock.set_real_delta_field(bias=1.0)
    return mock.mesh_delta_r


def make_galaxy_field(power, nmesh, boxsize, seed, bias_g, nbar, sample_seed):
    """Generate a galaxy density field and Poisson-sample a galaxy catalog.

    The galaxy field is delta_g = bias_g * delta_m. Galaxies are then
    Poisson-sampled from the density field with mean number density nbar.

    Parameters
    ----------
    power : callable
        Matter power spectrum P(k).
    nmesh : int
        Number of mesh cells per side.
    boxsize : float
        Periodic box size in Mpc/h.
    seed : int
        Random seed for the matter field realization (must match the matter field seed).
    bias_g : float
        Linear galaxy bias.
    nbar : float
        Mean galaxy number density in (h/Mpc)^3.
    sample_seed : int
        Random seed for the Poisson sampling step.

    Returns
    -------
    mesh_delta_r : RealField
        Real-space galaxy density field on the mesh.
    positions : array of shape (N, 3)
        Poisson-sampled galaxy positions in Mpc/h.
    """
    mock = EulerianLinearMock(power, nmesh=nmesh, boxsize=boxsize, seed=seed)
    mock.set_real_delta_field(bias=bias_g)
    mock.set_analytic_selection_function(nbar=nbar)
    mock.poisson_sample(seed=sample_seed)
    return mock.mesh_delta_r, mock.position


def make_electron_field(power, nmesh, boxsize, seed, bias_e, smoothing_radius):
    """Generate a free-electron density field with scale-dependent bias.

    The electron field is delta_e(k) = bias_e * exp[-(k * R_e)^2 / 2] * delta_m(k),
    where R_e is the Gaussian smoothing radius that suppresses small-scale power.

    Parameters
    ----------
    power : callable
        Matter power spectrum P(k).
    nmesh : int
        Number of mesh cells per side.
    boxsize : float
        Periodic box size in Mpc/h.
    seed : int
        Random seed for the matter field realization (must match the matter field seed).
    bias_e : float
        Electron bias amplitude.
    smoothing_radius : float
        Gaussian smoothing radius R_e in Mpc/h.

    Returns
    -------
    mesh_delta_r : RealField
        Real-space electron density field on the mesh.
    """
    mock = EulerianLinearMock(power, nmesh=nmesh, boxsize=boxsize, seed=seed)
    # Apply scale-dependent bias b_e * S_e(k) = b_e * exp[-(kR_e)^2/2] in Fourier space
    for kslab, slab in zip(mock.mesh_delta_k.slabs.x, mock.mesh_delta_k.slabs):
        k2 = sum(kk**2 for kk in kslab)
        slab[...] *= bias_e * np.exp(-k2 * smoothing_radius**2 / 2)
    mock.set_real_delta_field(bias=1.0)
    return mock.mesh_delta_r


if __name__ == '__main__':
    from frbao.validation_plots import plot_power_spectra

    nmesh = 256
    boxsize = 1000.   # Mpc/h
    seed = 42
    z = 0.5
    bias_g = 2.0
    bias_e = 1.0
    R_e = 5.0         # Mpc/h
    nbar = 1e-3       # (h/Mpc)^3

    pk, pk_nw = get_power_spectra(z=z)

    # Matter fields
    delta_m_w = make_matter_field(pk, nmesh, boxsize, seed)
    delta_m_nw = make_matter_field(pk_nw, nmesh, boxsize, seed)
    print(f'Matter (wiggle):    mean={np.mean(delta_m_w):.4f}, std={np.std(delta_m_w):.4f}')
    print(f'Matter (no-wiggle): mean={np.mean(delta_m_nw):.4f}, std={np.std(delta_m_nw):.4f}')

    # Galaxy fields
    delta_g_w, pos_g_w = make_galaxy_field(pk, nmesh, boxsize, seed, bias_g, nbar, sample_seed=seed + 1)
    delta_g_nw, pos_g_nw = make_galaxy_field(pk_nw, nmesh, boxsize, seed, bias_g, nbar, sample_seed=seed + 1)
    print(f'Galaxy (wiggle):    mean={np.mean(delta_g_w):.4f}, std={np.std(delta_g_w):.4f}, N={len(pos_g_w)}')
    print(f'Galaxy (no-wiggle): mean={np.mean(delta_g_nw):.4f}, std={np.std(delta_g_nw):.4f}, N={len(pos_g_nw)}')

    # Electron fields
    delta_e_w = make_electron_field(pk, nmesh, boxsize, seed, bias_e, R_e)
    delta_e_nw = make_electron_field(pk_nw, nmesh, boxsize, seed, bias_e, R_e)
    print(f'Electron (wiggle):    mean={np.mean(delta_e_w):.4f}, std={np.std(delta_e_w):.4f}')
    print(f'Electron (no-wiggle): mean={np.mean(delta_e_nw):.4f}, std={np.std(delta_e_nw):.4f}')

    plot_power_spectra(pk, pk_nw)
