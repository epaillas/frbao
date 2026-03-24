"""Steps 0–3 of the FRB × galaxy pipeline.

Step 0: Place FRBs as uniform random background sightlines.
Step 1: Assign a dispersion measure (DM) to each FRB by integrating delta_e along +z.
Step 2: Subtract the mean DM to produce a centred delta_DM per FRB.
Step 3: Select foreground galaxies within the chosen redshift slab.
"""

import numpy as np


def place_frbs(boxsize, n_frb, z_source, seed=None):
    """Place FRBs on a fixed source plane behind the foreground galaxy slab.

    Transverse positions (x, y) are drawn uniformly across the box; all
    sources are placed at the same depth z = z_source.

    Parameters
    ----------
    boxsize : float
        Periodic box size in Mpc/h.
    n_frb : int
        Number of FRBs to place.
    z_source : float
        Fixed source-plane depth in Mpc/h. Must satisfy z_source > z_g_max + dz_buffer.
    seed : int, optional
        Random seed for reproducibility.

    Returns
    -------
    positions : ndarray of shape (n_frb, 3)
        FRB positions [x, y, z] in Mpc/h.
    """
    rng = np.random.default_rng(seed)
    x = rng.uniform(0., boxsize, size=n_frb)
    y = rng.uniform(0., boxsize, size=n_frb)
    z = np.full(n_frb, z_source)
    return np.column_stack([x, y, z])


def assign_dm(frb_positions, delta_e, boxsize, a_dm=1.0):
    """Assign a dispersion measure to each FRB by line-of-sight integration.

    For each FRB at (x_i, y_i, z_i), integrates the electron overdensity
    along the fixed +z LOS from the observer plane z = 0 to z_i:

        DeltaDM_i = a_dm * sum_{j=0}^{iz_max - 1} delta_e[ix, iy, j] * dz

    The transverse lookup uses nearest-grid-point (NGP) assignment.

    Parameters
    ----------
    frb_positions : ndarray of shape (n_frb, 3)
        FRB positions [x, y, z] in Mpc/h.
    delta_e : RealField or ndarray of shape (nmesh, nmesh, nmesh)
        Electron overdensity field. Axis order is [x, y, z].
    boxsize : float
        Periodic box size in Mpc/h.
    a_dm : float, optional
        Overall DM amplitude prefactor (default 1.0).

    Returns
    -------
    dm : ndarray of shape (n_frb,)
        Raw DM value per FRB.
    """
    delta = np.array(delta_e)   # (nmesh, nmesh, nmesh)
    nmesh = delta.shape[0]
    dz = boxsize / nmesh

    x, y, z = frb_positions[:, 0], frb_positions[:, 1], frb_positions[:, 2]

    ix = (x / boxsize * nmesh).astype(int) % nmesh
    iy = (y / boxsize * nmesh).astype(int) % nmesh
    iz_max = np.clip((z / boxsize * nmesh).astype(int), 0, nmesh)

    dm = np.array([
        a_dm * delta[ix[i], iy[i], :iz_max[i]].sum() * dz
        for i in range(len(frb_positions))
    ])
    return dm


def subtract_mean_dm(dm):
    """Subtract the ensemble mean DM to produce a centred fluctuation.

    Parameters
    ----------
    dm : ndarray of shape (n_frb,)
        Raw DM values.

    Returns
    -------
    ddm : ndarray of shape (n_frb,)
        Centred DM fluctuations delta_DM = DM - <DM>.
    """
    return dm - dm.mean()


def select_galaxies_in_slab(positions, zmin, zmax):
    """Select galaxies whose z-coordinate falls within [zmin, zmax).

    Parameters
    ----------
    positions : ndarray of shape (N, 3)
        Galaxy positions in Mpc/h.
    zmin : float
        Minimum z boundary of the foreground slab in Mpc/h.
    zmax : float
        Maximum z boundary of the foreground slab in Mpc/h.

    Returns
    -------
    slab_positions : ndarray of shape (M, 3)
        Positions of galaxies inside the slab.
    """
    z = positions[:, 2]
    mask = (z >= zmin) & (z < zmax)
    return positions[mask]


if __name__ == '__main__':
    from frbao.make_density_fields import (
        get_power_spectra, make_electron_field, make_galaxy_field,
    )

    # --- Configuration ---
    nmesh = 256
    boxsize = 1000.0    # Mpc/h
    seed = 42
    z = 0.5
    bias_g = 2.0
    bias_e = 1.0
    R_e = 5.0           # Mpc/h
    nbar = 1e-3         # (h/Mpc)^3
    n_frb = 10_000
    galaxy_slab = (200.0, 400.0)   # Mpc/h
    dz_buffer = 50.0               # Mpc/h
    z_source = galaxy_slab[1] + dz_buffer + 100.0

    # --- Generate fields ---
    pk, pk_nw = get_power_spectra(z=z)
    delta_e_w  = make_electron_field(pk,    nmesh, boxsize, seed, bias_e, R_e)
    delta_e_nw = make_electron_field(pk_nw, nmesh, boxsize, seed, bias_e, R_e)
    _, pos_g_w  = make_galaxy_field(pk,    nmesh, boxsize, seed, bias_g, nbar, sample_seed=seed + 1)
    _, pos_g_nw = make_galaxy_field(pk_nw, nmesh, boxsize, seed, bias_g, nbar, sample_seed=seed + 1)

    # --- Step 0: place FRBs ---
    frb_pos = place_frbs(boxsize, n_frb, z_source, seed=seed + 2)
    print(f'Step 0 | FRBs: N={len(frb_pos)}, z_source={z_source:.1f} Mpc/h, '
          f'x in [{frb_pos[:, 0].min():.1f}, {frb_pos[:, 0].max():.1f}]')

    # --- Step 1: assign DM ---
    dm_w  = assign_dm(frb_pos, delta_e_w,  boxsize)
    dm_nw = assign_dm(frb_pos, delta_e_nw, boxsize)
    print(f'Step 1 | DM (wiggle):    mean={dm_w.mean():.4f}, std={dm_w.std():.4f}')
    print(f'Step 1 | DM (no-wiggle): mean={dm_nw.mean():.4f}, std={dm_nw.std():.4f}')

    # --- Step 2: subtract mean ---
    ddm_w  = subtract_mean_dm(dm_w)
    ddm_nw = subtract_mean_dm(dm_nw)
    print(f'Step 2 | δDM (wiggle):    mean={ddm_w.mean():.2e}, std={ddm_w.std():.4f}')
    print(f'Step 2 | δDM (no-wiggle): mean={ddm_nw.mean():.2e}, std={ddm_nw.std():.4f}')

    # --- Step 3: foreground galaxy slab ---
    gals_w  = select_galaxies_in_slab(pos_g_w,  *galaxy_slab)
    gals_nw = select_galaxies_in_slab(pos_g_nw, *galaxy_slab)
    print(f'Step 3 | Galaxies in slab (wiggle):    N={len(gals_w)} '
          f'(z in [{galaxy_slab[0]}, {galaxy_slab[1]}])')
    print(f'Step 3 | Galaxies in slab (no-wiggle): N={len(gals_nw)}')
