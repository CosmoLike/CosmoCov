# --->> This file contains the functions R_1(z,k) and R_K(z,k), which return the first-order
# isotropic and tidal power spectrum responses, respectively. This finds application in the 
# evaluation of the SSC covariance, which take these as input (see Barreira, Krause, Schmidt, arXiv:1711.07467).

# --->> The current implementation assumes existing tables of the nonlinear power spectrum, its derivatives and growth-only first order responses.

# ======================================================== #
# Be brute and load all stuff from libraries
# ======================================================== #
from numpy import *
from scipy import interpolate, integrate

lpath = '../cosmolike_core/theory/'

# ================================================================================ #
# Load spectra and growth-only response tables and create 2D interpolators
# ================================================================================ #
# Create the interpolators that are needed; it does trick of setting to zero (or linear result for responses) for k = 0 (to prevent crashes/problems when k<k_min)
# This needs to be adapted if the format of the files changes
def make_interpolator_forP(path1, path2, value_at_k0):
    filetmp  = loadtxt(path1, skiprows = 1)
    ztmp     = loadtxt(path2)
    ktmp     = append(0., filetmp[:,0])
    Ptmp     = vstack([value_at_k0*ones(len(ztmp)), filetmp[:,1:]])
    return interpolate.interp2d(ztmp, ktmp, Ptmp, kind='cubic')
# This interpolator add a correction for k>kmax
def make_interpolator_forG(path1, path2, value_at_k0, index_highk, asymptote, doextrapolation):
    filetmp  = loadtxt(path1, skiprows = 1)
    ztmp     = loadtxt(path2)
    # Set to linear theory for low-k
    ktmp     = append(0., filetmp[:,0])
    Gtmp     = vstack([value_at_k0*ones(len(ztmp)), filetmp[:,1:]])
    if(doextrapolation):
        # Set correction for high-k
        khig = 10.**linspace(log10(max(ktmp+0.000001)), 4, 500)
        kout = append(ktmp, khig)
        Ghig = zeros([len(kout), len(ztmp)])
        Ghig[0:len(ktmp), :] = Gtmp
        for i in range(len(ztmp)):
            G_at_kmax = Ghig[len(ktmp)-1, i]
            kmax      = ktmp[-1]
            Ghig[len(ktmp)::, i] = (kmax/khig)**index_highk * (G_at_kmax - asymptote) + asymptote
        return interpolate.interp2d(ztmp, kout, Ghig, kind='cubic')
    else:
        return interpolate.interp2d(ztmp, ktmp, Gtmp, kind='cubic')

# P_m(z, k)
Pnl_int   = make_interpolator_forP(lpath+'lookup_tables/P_m_runmode.dat', lpath+'lookup_tables/P_m_zvalues.dat', 0.0)
# dP_m/dk
dPnl_int  = make_interpolator_forP(lpath+'lookup_tables/dP_m_postrunmode.dat', lpath+'lookup_tables/P_m_zvalues.dat', 0.0)

# G_1(z, k)
G1_int   = make_interpolator_forG(lpath+'lookup_tables/Resp_G1_fromsims.dat', lpath+'lookup_tables/Resp_zvalues_fromsims.dat', 26./21, 0.5, -3./4, True)
# G_K(z, k)
GK_int   = make_interpolator_forG(lpath+'lookup_tables/Resp_GK_fromsims_Andreas.dat', lpath+'lookup_tables/Resp_zvalues_fromsims_Andreas.dat', 8./7, 0.5, -9./4, True)

# ================================================================================ #
# R_1(z, k), R_K(z, k) and function with all the combinations
# ================================================================================ #
def R_1_int(z, k):
    return (1. - (1./3)*k*dPnl_int(z, k)/Pnl_int(z, k) + G1_int(z, k))[0]
def R_K_int(z, k):
    return (GK_int(z, k) - k*dPnl_int(z, k)/Pnl_int(z, k))[0]

