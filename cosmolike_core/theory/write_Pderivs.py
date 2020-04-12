# --->> This file is prepares tables of dP/dk/P and d^2P/dk^2/P to be read in by covNG_resp.py

from numpy import *
from scipy import interpolate
import scipy.ndimage as ndimage

lpath = "../cosmolike_core/theory/"

# Define this to be used below
def Spline3(xx,yy,Order,X): ## Spline function (Order is the order of the derivative)
    ff=interpolate.splrep(xx,yy,k=3)
    ss=interpolate.splev(X,ff,Order)
    return ss

dosmoothing = False
if(dosmoothing):
    def smoother(ddP, kk, sigma, kmin, kmax):
        smoothed =  ndimage.filters.gaussian_filter(ddP, sigma)
        out = copy(ddP)
        out[where(kk<kmin)] = smoothed[where(kk<kmin)]
        out[where(kk>kmax)] = smoothed[where(kk>kmax)]
        return out
    print('')
    print('Smoothing the 2nd derivative to remove noise')
    print('')

# If the format of the files changes, this needs to change
# There are smarter ways of doing this, but this will do the job too
def write_dP_ddP(path_P, path_z, filename1, header1, filename2, header2):
    fin = loadtxt(path_P, skiprows = 1)
    kin = fin[:,0]
    Pin = fin[:,1:]
    zin = loadtxt(path_z)
    dP  = zeros(shape(Pin))
    ddP = zeros(shape(Pin))
    # Build derivatives with Spline function
    for j in range(len(zin)):
        Ptmp = Pin[:,j]
        for i in range(len(kin)):
            dP[i,j]  = Spline3(kin, Ptmp, 1, kin[i])
            ddP[i,j] = Spline3(kin, Ptmp, 2, kin[i])
    # Smooth the 2nd derivative to remove noise
    if(dosmoothing):
        for j in range(len(zin)):
            ddP[:,j] = smoother(ddP[:,j], kin, 0.75, 0.1, 0.1)
    # Write to file
    fout1 = open(filename1, 'w')
    fout2 = open(filename2, 'w')
    fout1.write('k[h/Mpc]'); fout1.write(' ');
    fout2.write('k[h/Mpc]'); fout2.write(' ');
    for k in range(len(zin)):
        fout1.write(header1+"(z="+"%0.2f" % zin[k]+")"); fout1.write(' ')
        fout2.write(header2+"(z="+"%0.2f" % zin[k]+")"); fout2.write(' ')
    fout1.write('\n')
    fout2.write('\n')
    for i in range(len(kin)):
        fout1.write(str(kin[i])); fout1.write(' ')
        fout2.write(str(kin[i])); fout2.write(' ')
        for j in range(len(zin)):
            fout1.write(str( dP[i,j])) ; fout1.write(' ')
            fout2.write(str(ddP[i,j])); fout2.write(' ')
        fout1.write('\n')
        fout2.write('\n')
    fout1.close()
    fout2.close()

write_dP_ddP(lpath+'lookup_tables/P_L_class.dat', lpath+'lookup_tables/P_L_zvalues.dat', lpath+'lookup_tables/dP_L_postclass.dat', "dP/dk", lpath+'lookup_tables/ddP_L_postclass.dat', "d^2P/dk^2")
write_dP_ddP(lpath+'lookup_tables/P_m_runmode.dat',lpath+ 'lookup_tables/P_m_zvalues.dat', lpath+'lookup_tables/dP_m_postrunmode.dat', "dP/dk", lpath+'lookup_tables/ddP_m_postrunmode.dat', "d^2P/dk^2")

