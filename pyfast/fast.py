'''
Created on 12/02/2014

@author: Steven

This is supposed to run the whole idea
'''

from hmf import MassFunction
import numpy as np
from scipy.integrate import cumtrapz
from scipy.interpolate import UnivariateSpline as spline

def fast(filename, cell_size, **mf_kwargs):

    # Read in the dm particle possitions from GADGET file
    dm_pos, mpart = _get_gadget(filename)

    # Assign particles to a grid. At this stage, don't get rid of original pos
    grid, edges = _make_grid(dm_pos, cell_size)
    grid *= mpart

    cdf, icdf = prepare_mf(mpart, grid, mf_kwargs)

    halomasses = choose_halo_masses(cdf, icdf, mpart, grid)

    halopos = place_halos(halomasses, grid, dm_pos)

    return halopos

def _get_gadget(filename):
    """
    Reads a gadget binary specified by filename
    """
    positions_dtype = np.dtype([("x", 'f'),
                                ("y", 'f'),
                                ("z", 'f')
                                ])

    header_dtype = np.dtype([('npart', ('i', 6)),
                             ('massarr', ('d', 6)),
                             ('time', 'd'),
                             ('z', 'd'),
                             ('FlagSfr', 'i'),
                             ('FlagFeedback', 'i'),
                             ('n_all', ('i', 6)),
                             ('FlagCooling', 'i'),
                             ('num_files', 'i'),
                             ('boxsize', 'd'),
                             ('omegam', 'd'),
                             ('omegav', 'd'),
                             ('h', 'd'),
                             ('FlagAge', 'i'),
                             ('FlagMetals', 'i'),
                             ('n_all_HW', ('i', 6)),
                             ('flag_entr_ics', 'i')
                             ])

    with open(filename, "rb") as f:
        # Header block
        f.read(4)
        header = np.fromfile(f, dtype=header_dtype, count=1)
        # So far only 196 out of 256 so read the rest
        f.read(256 - 196)
        f.read(4)

        # Positions Block
        f.read(4)
        pos = np.fromfile(f, dtype=positions_dtype, count=header['npart'][0][1])
        f.read(4)

        # Velocities Block (don't save it anywhere)
        f.read(4)
        np.fromfile(f, dtype=positions_dtype, count=header['npart'][0][1])
        f.read(4)

        # ID's block
        f.read(4)
        ids = np.fromfile(f, dtype='i', count=header['npart'][0][1])


    ids = np.array(ids - 1)

    indices = np.argsort(ids)
    pos = pos[indices]
    ids = sorted(ids)

    header_dict = {}
    for name in header.dtype.names:
        header_dict[name] = header[name][0]

    return pos.view(np.float32).reshape(pos.shape + (-1,)) , header_dict['boxsize'], header_dict['massarr'][1]


def _make_grid(dm_pos, boxsize, cell_size):
    bins = int(boxsize / cell_size)
    H, edges = np.histogramdd(dm_pos, bins)
    # cell_assigned = np.array((dm_pos.shape[0], 3))
    return H, edges

def prepare_mf(mpart, grid, mf_kwargs):
    # The maximum mass we will care about.
    max_mass = grid.max() * mpart

    M = np.linspace(np.log10(mpart), np.log10(max_mass), 2000)

    mf = MassFunction(M=M, **mf_kwargs).dndm

    cumfunc = cumtrapz(M, 10 ** M * mf)

    cdf = spline(M, np.log10(cumfunc))
    icdf = spline(np.log10(cumfunc), M)

    return cdf, icdf

def choose_halo_masses(cdf, icdf, mpart, grid):

    tol = 2 * mpart
    i = 0
    masses = []
    for mcell in np.nditer(grid):

        if mcell < tol:
            print "Mass in cell ", i, " is negligible"
            masses.append([])
            continue

        diff = 3 * mpart
        m_in_cell = mcell
        j = 0

        while diff > tol:
            maxcum = cdf(np.log10(m_in_cell))
            mean_halo_mass = 10 ** maxcum / (m_in_cell - mpart)
            N_expected = 1.2 * m_in_cell / mean_halo_mass

            x = np.random.random(N_expected) * maxcum
            m = icdf(x)
            cumsum = np.cumsum(m)
            cross_ind = np.where(cumsum > m_in_cell)[0]
            over = abs(cumsum[cross_ind] - m_in_cell)
            under = abs(cumsum[cross_ind - 1] - m_in_cell)
            if over < tol and under < tol:
                if over < under:
                    m = m[:cross_ind - 1]
                else:
                    m = m[:cross_ind]
                diff = over
            elif over < tol:
                m = m[:cross_ind]
                diff = over
            else:
                m = m[:cross_ind - 1]
                diff = under
                m_in_cell -= cumsum[cross_ind - 1]

            if j > 0:
                cell_masses = np.concatenate((cell_masses, m))
            else:
                cell_masses = m
            j += 1

        masses.append(cell_masses)
        i += 1

    return masses

def place_halos(halomasses, grid, dm_pos,edges):
    
    for cell in halomasses:
        
