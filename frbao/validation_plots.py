"""Validation plots for frbao."""

import numpy as np
import matplotlib.pyplot as plt


def plot_power_spectra(pk_wiggle, pk_nowiggle, kmin=1e-3, kmax=0.5, nk=500):
    """Plot wiggle vs no-wiggle P(k) as k * P(k) vs k (linear axes).

    Parameters
    ----------
    pk_wiggle : callable
        Wiggle power spectrum, pk(k) -> P(k).
    pk_nowiggle : callable
        No-wiggle power spectrum.
    kmin : float, optional
        Minimum wavenumber in h/Mpc.
    kmax : float, optional
        Maximum wavenumber in h/Mpc.
    nk : int, optional
        Number of k points.
    """
    k = np.geomspace(kmin, kmax, nk)
    fig, ax = plt.subplots()
    ax.plot(k, k * pk_wiggle(k), label='Wiggle')
    ax.plot(k, k * pk_nowiggle(k), label='No wiggle')
    ax.set_xlabel(r'$k$ [$h\,\mathrm{Mpc}^{-1}$]')
    ax.set_ylabel(r'$k \, P(k)$ [$h^{-2}\,\mathrm{Mpc}^2$]')
    ax.legend()
    plt.show()
