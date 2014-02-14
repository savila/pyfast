'''
Created on 12/02/2014

@author: Steven

This is supposed to run the whole idea
'''

MAKE_PLOTS = False
VERBOSE = True

from hmf import MassFunction
import numpy as np
from scipy.integrate import cumtrapz
from scipy.interpolate import InterpolatedUnivariateSpline as spline
import sys
#SA-----
from numpy.ctypeslib import ndpointer
lib = ctypes.cdll.LoadLibrary('./place_halos.so')
fun = lib.place_halos
fun.restype = None
fun.argtypes = [ctypes.long, ndpointer(ctypes.double), ctypes.int, ndpointer(ctypes.int),
ctypes.int, ndpointer(ctypes.double), ndpointer(ctypes.double), ndpointer(ctypes.double), ctypes.double,
ndpointer(ctypes.double), ndpointer(ctypes.double), ndpointer(ctypes.double))

#void place_halos(long NHalosTot, double *HaloMass, int NCells, int *FirstHaloInCell,  
#int NTotPart, double *PartX, double *PartY, double *PartZ, double L, 
#double *HaloX, double *HaloY, double *HaloZ) 
#-----SA




if MAKE_PLOTS:
    import matplotlib.pyplot as plt

def fast(filename, cell_size=1, frac_mass_in_halos=1.0, **mf_kwargs):

    # Read in the dm particle possitions from GADGET file
    if VERBOSE:
        print "READING GADGET FILE..."
    dm_pos, header = _get_gadget(filename)
    mpart = header['massarr'][1] * 1e10
    boxsize = header['boxsize']

    if type(cell_size) == int:
        cell_size = boxsize / cell_size
    # Assign particles to a grid. At this stage, don't get rid of original pos
    if VERBOSE:
        print "CREATING GRID CELLS..."
    grid, cell_size = _make_grid(dm_pos, boxsize, cell_size)
    grid *= mpart

    # Update mf_kwargs with the things we know are true
    mf_kwargs.update({"z":header['z'],
                      "omegav":header['omegav'],
                      'h':header['h'],
                      'omegam':header['omegam']})

    if VERBOSE:
        print "PREPARING THEORETICAL MASS FUNCTION..."
    cdf, icdf, M, dndm, m_outside_range = prepare_mf(mpart, grid, mf_kwargs)

    # Take off the mass that we'll never populate.
    grid -= m_outside_range
    if VERBOSE:
        print "ASSIGNING HALO MASSES"
    halomasses, offsets = choose_halo_masses(cdf, icdf, mpart, grid, boxsize)

    if MAKE_PLOTS:
        hist, edges = np.histogram(np.log(halomasses), bins=20)
        hist = hist.astype('float64')
        hist /= (edges[1] - edges[0])
        centres = (edges[:-1] + edges[1:]) / 2
        plt.clf()
        plt.scatter(np.exp(centres), hist)

        plt.plot(10 ** M, 10 ** M * dndm * boxsize ** 3)
        plt.yscale('log')
        plt.xscale('log')
        plt.savefig("dndm.pdf")

    sys.exit("stopping here for now")
    halopos = place_halos(halomasses, grid, dm_pos, offsets)

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

    return pos.view(np.float32).reshape(pos.shape + (-1,)) , header_dict


def _make_grid(dm_pos, boxsize, cell_size):
    bins = int(boxsize / cell_size)
    H = np.histogramdd(dm_pos, bins)[0]
    cell_size = boxsize / bins
    return H, cell_size

def prepare_mf(mpart, grid, mf_kwargs):
    M = np.linspace(np.log10(mpart), np.log10(grid.max()), 2000)
    mf_obj = MassFunction(M=M, **mf_kwargs)

    mf = mf_obj.dndm
    m_outside_range = mf_obj.mltm[0] + mf_obj.mgtm[-1]

    cumfunc = cumtrapz(10 ** M * mf, M, initial=0) * np.log(10)

    cdf = spline(M, cumfunc, k=3)
    icdf = spline(cumfunc, M, k=3)

    if MAKE_PLOTS:
        plt.clf()
        plt.plot(M, cumfunc)

        plt.plot(M, cdf(M))
        plt.savefig("cumfunc.pdf")

        plt.clf()
        mcumfunc = cumtrapz(10 ** (2 * M) * mf, dx=M[1] - M[0], initial=1e-20) * np.log(10)
        plt.plot(M, mcumfunc)
        plt.savefig("mcumfunc.pdf")

        # How much mass is above 10**12.5
        minvcumfunc = cumtrapz(10 ** (2 * M[::-1]) * mf[::-1], dx=M[1] - M[0]) * np.log(10)
        minvcumfunc = np.concatenate((np.array([minvcumfunc[0]]), minvcumfunc))
        minvcumfunc /= minvcumfunc[-1]
        plt.clf()
        plt.plot(M, minvcumfunc[::-1])
        plt.yscale('log')
        plt.grid(True)
        plt.savefig("minvcumfunc.pdf")

    return cdf, icdf, M, mf, m_outside_range

def choose_halo_masses(cdf, icdf, mpart, grid, boxsize):

    tol = 3 * mpart
    i = 0
    offsets = np.zeros(grid.size)
    for mcell in np.nditer(grid):

        if mcell < tol:
            print "Mass in cell ", i, " is negligible"
            continue

        diff = 4 * mpart
        m_in_cell = mcell.copy()
        j = 0

        if VERBOSE:
            print "CELL: ", i
        while diff > tol:
            # Figure out expected number of halos to get back mass in cell
            maxcum = cdf(np.log10(m_in_cell))
            # mean_halo_mass = 10 ** mcdf(np.log10(m_in_cell)) / (m_in_cell - mpart)
            N_expected = int(maxcum * boxsize ** 3 / grid.size)
            # Generate random variates from 0 to maxcum
            x = np.random.random(N_expected) * maxcum

            # Generate halo masses from mf distribution
            m = 10 ** icdf(x)

            # Make sure we don't have more or less mass than needed
            cumsum = np.cumsum(m)
            try:
                cross_ind = np.where(cumsum > m_in_cell)[0][0]
            except IndexError:
                cross_ind = len(cumsum) - 1

            if VERBOSE:
                print "  ITERATION: ", j, "Num. Vars: ", len(cumsum), ", Final # Halos: ", cross_ind + 1
            over = abs(cumsum[cross_ind] - m_in_cell)
            under = abs(cumsum[cross_ind - 1] - m_in_cell)
            if over < tol and under < tol:
                if over < under:
                    m = m[:cross_ind ]
                else:
                    m = m[:cross_ind + 1]
                diff = over
            elif over < tol:
                m = m[:cross_ind + 1]
                diff = over
            else:
                m = m[:cross_ind]
                diff = under
                m_in_cell -= cumsum[cross_ind - 1]

            # Save the halo masses in this cell (perhaps append them)
            if i == 0 and j == 0:
                cell_masses = m
            else:
                cell_masses = np.concatenate((cell_masses, m))

            j += 1

        cell_masses[offsets[i]:] = np.sort(cell_masses[offsets[i]:])[::-1]
        if i < grid.size - 1:
            offsets[i + 1] = len(cell_masses)
        i += 1

    return cell_masses, offsets

def place_halos(halomasses, grid, dm_pos, edges):

    for cell in halomasses:
        pass
