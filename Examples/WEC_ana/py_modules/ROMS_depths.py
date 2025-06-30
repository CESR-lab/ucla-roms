######################################
__title__          = "ROMS_depths.py"
__author__         = "Daniel Dauhajre, Cigdem Akan"
__date__           = "August 2017"
__email__          = "ddauhajre@atmos.ucla.edu"
__python_version__ = "2.7.9"

'''
PYTHON LIBRARY OF FUNCTIONS TO
CONVERT SIGMA LEVELS TO DEPTHS OF A ROMS
GRID BASED ON OUTPUT (i.e., free surface)
'''
######################################

import numpy as np
from netCDF4 import Dataset


def get_zr_zw_tind(nc_roms, nc_grd, tind, dim_bounds):
    '''
    GET DEPTHS FOR A SPECIFIC netcdf ROMS output file

    nc_roms --> single netcdf file with roms output
    nc_grd  --> single netcdf file of roms grid
    tind    --> time index
    dim_bounds   ---> [eta_0, eta_1, xi_0, xi_1] list of spatial bounds
    '''
    #######################
    #DECLARE DEPTH ARRAY
    ######################
    nt = len(nc_roms.variables['ocean_time'][:])
    N = len(nc_roms.variables['u'][0,:,0,0])
    [Ly_all,Lx_all] = nc_grd.variables['pm'].shape
    [eta_0, eta_1, xi_0, xi_1] = dim_bounds
    Ly = eta_1 - eta_0
    Lx = xi_1 - xi_0

    z_r = np.zeros([N,Ly_all,Lx_all])
    z_w = np.zeros([N+1,Ly_all,Lx_all])

    ####################################
    # SET DEPTH CALCULATION ATTRIBUTES
    #####################################

    # SEE CODE IN set_depth() for documentation on these values
    Vtrans  = 2
    Vstret = 4



    #ACCESS GLOBAL ATTRIBUES
    hc      = getattr(nc_roms, 'hc')
    theta_s = getattr(nc_roms, 'theta_s')
    theta_b = getattr(nc_roms, 'theta_b')
    N       = len(nc_roms.dimensions['s_rho'])

    #BOTTOM TOPOGRAPHY
    h = nc_grd.variables['h'][eta_0:eta_1,xi_0:xi_1]

    zeta = nc_roms.variables['zeta'][tind,eta_0:eta_1,xi_0:xi_1]
    print('Calculating z_r, z_w at time-step =', tind)
    z_r = set_depth(Vtrans, Vstret, theta_s, theta_b, hc, N, 1, h, zeta).T
    z_w = set_depth(Vtrans, Vstret, theta_s, theta_b, hc, N, 5, h, zeta).T


    return np.swapaxes(z_r,1,2), np.swapaxes(z_w,1,2)
    ##########################################################





def get_zs(nc_roms, nc_grd, levs, dim_bounds, igrid=1):
    '''
    GET DEPTHS FOR A SPECIFIC netcdf ROMS output file

    nc_roms --> single netcdf file with roms output
    nc_grd  --> single netcdf file of roms grid
    dim_bounds   ---> [eta_0, eta_1, xi_0, xi_1] list of spatial bounds
    '''
    #######################
    #DECLARE DEPTH ARRAY
    ######################
    [nt,N,Ly_all,Lx_all] = nc_roms.variables['temp'].shape
    [eta_0, eta_1, xi_0, xi_1] = dim_bounds
    Ly = eta_1 - eta_0
    Lx = xi_1 - xi_0
    nlevs = len(levs)


    z = np.zeros([nt,nlevs,Ly_all,Lx_all])

    ####################################
    # SET DEPTH CALCULATION ATTRIBUTES
    #####################################

    # SEE CODE IN set_depth.py for documentation on these values
    Vtrans  = 2
    Vstret = 4



    #ACCESS GLOBAL ATTRIBUES
    hc      = getattr(nc_roms, 'hc')
    theta_s = getattr(nc_roms, 'theta_s')
    theta_b = getattr(nc_roms, 'theta_b')
    N       = len(nc_roms.dimensions['s_rho'])

    #BOTTOM TOPOGRAPHY
    h = nc_grd.variables['h'][eta_0:eta_1,xi_0:xi_1]

    ####################################
    # CALCULATE DEPTHS AT TIME-STEPS
    ####################################

    for n in range(nt):
        zeta = nc_roms.variables['zeta'][n,eta_0:eta_1,xi_0:xi_1]
        print('Calculating depths at time-step =', n)
        print('     levs = ', levs)
        z_temp = set_depth(Vtrans, Vstret, theta_s, theta_b, hc, N, igrid, h, zeta)
        for k in range(nlevs):
            z[n,k,:,:] = z_temp[:,:,levs[k]]
    return z
    ##########################################################


    ######################################
    # BELOW: FUNCTIONS BY CIGDEM AKA
    # FOR DEPTH CONVERSION
    #######################################



"""
 Given a batymetry (h), free-surface (zeta) and terrain-following
 parameters, this function computes the 3D depths for the requested
 C-grid location. If the free-surface is not provided, a zero value
 is assumed resulting in unperturb depths.  This function can be
 used when generating initial conditions or climatology data for
 an application. Check the following link for details:

    https://www.myroms.org/wiki/index.php/Vertical_S-coordinate

 On Input:

    Vtransform    Vertical transformation equation:

                    Vtransform = 1,   original transformation

                    z(x,y,s,t)=Zo(x,y,s)+zeta(x,y,t)*[1+Zo(x,y,s)/h(x,y)]

                    Zo(x,y,s)=hc*s+[h(x,y)-hc]*C(s)

                    Vtransform = 2,   new transformation

                    z(x,y,s,t)=zeta(x,y,t)+[zeta(x,y,t)+h(x,y)]*Zo(x,y,s)

                    Zo(x,y,s)=[hc*s(k)+h(x,y)*C(k)]/[hc+h(x,y)]

    Vstretching   Vertical stretching function:
                    Vstretching = 1,  original (Song and Haidvogel, 1994)
                    Vstretching = 2,  A. Shchepetkin (UCLA-ROMS, 2005)
                    Vstretching = 3,  R. Geyer BBL refinement
                    Vstretching = 4,  A. Shchepetkin (UCLA-ROMS, 2010)

    theta_s       S-coordinate surface control parameter (scalar)

    theta_b       S-coordinate bottom control parameter (scalar)

    hc            Width (m) of surface or bottom boundary layer in which
                    higher vertical resolution is required during
                    stretching (scalar)

    N             Number of vertical levels (scalar)

    igrid         Staggered grid C-type (integer):
                    igrid=1  => density points
                    igrid=2  => streamfunction points
                    igrid=3  => u-velocity points
                    igrid=4  => v-velocity points
                    igrid=5  => w-velocity points

    h             Bottom depth, 2D array at RHO-points (m, positive),
                    h(1:Lp+1,1:Mp+1)

    zeta          Free-surface, 2D array at RHO-points (m), OPTIONAL,
                    zeta(1:Lp+1,1:Mp+1)

 On Output:

    z             Depths (m, negative), 3D array
"""
def set_depth( Vtr, Vstr, thts, thtb, hc, N, igrid, h, zeta ):
    Np      = N+1
    Lp,Mp   = np.shape(h)
    L       = Lp-1
    M       = Mp-1
    if (igrid==5):
        z   = np.empty((Lp,Mp,Np))
    else:
        z   = np.empty((Lp,Mp,N))

    hmin    = np.min(h)
    hmax    = np.max(h)

    if (igrid == 5):
        kgrid=1
    else:
        kgrid=0

    s,C = stretching(Vstr, thts, thtb, hc, N, kgrid);
    #-----------------------------------------------------------------------
    #  Average bathymetry and free-surface at requested C-grid type.
    #-----------------------------------------------------------------------

    if (igrid==1):
        hr    = h
        zetar = zeta
    elif (igrid==2):
        hp    = 0.25*(h[0:L,0:M]+h[1:Lp,0:M]+h[0:L,1:Mp]+h[1:Lp,1:Mp])
        zetap = 0.25*(zeta[0:L,0:M]+zeta[1:Lp,0:M]+zeta[0:L,1:Mp]+zeta[1:Lp,1:Mp])
    elif (igrid==3):
        hu    = 0.5*(h[0:L,0:Mp]+h[1:Lp,0:Mp])
        zetau = 0.5*(zeta[0:L,0:Mp]+zeta[1:Lp,0:Mp])
    elif (igrid==4):
        hv    = 0.5*(h[0:Lp,0:M]+h[0:Lp,1:Mp])
        zetav = 0.5*(zeta[0:Lp,0:M]+zeta[0:Lp,1:Mp])
    elif (igrid==5):
        hr    = h
        zetar = zeta

    #----------------------------------------------------------------------
    # Compute depths (m) at requested C-grid location.
    #----------------------------------------------------------------------
    if (Vtr == 1):
        if (igrid==1):
            for k in range (0,N):
                z0 = (s[k]-C[k])*hc + C[k]*hr
                z[:,:,k] = z0 + zetar*(1.0 + z0/hr)
        elif (igrid==2):
            for k in range (0,N):
                z0 = (s[k]-C[k])*hc + C[k]*hp
                z[:,:,k] = z0 + zetap*(1.0 + z0/hp)
        elif (igrid==3):
            for k in range (0,N):
                z0 = (s[k]-C[k])*hc + C[k]*hu
                z[:,:,k] = z0 + zetau*(1.0 + z0/hu)
        elif (igrid==4):
            for k in range (0,N):
                z0 = (s[k]-C[k])*hc + C[k]*hv
                z[:,:,k] = z0 + zetav*(1.0 + z0/hv)
        elif (igrid==5):
            z[:,:,0] = -hr
            for k in range (0,Np):
                z0 = (s[k]-C[k])*hc + C[k]*hr
                z[:,:,k] = z0 + zetar*(1.0 + z0/hr)
    elif (Vtr==2):
        if (igrid==1):
            for k in range (0,N):
                z0 = (hc*s[k]+C[k]*hr)/(hc+hr)
                z[:,:,k] = zetar+(zeta+hr)*z0
        elif (igrid==2):
            for k in range (0,N):
                z0 = (hc*s[k]+C[k]*hp)/(hc+hp)
                z[:,:,k] = zetap+(zetap+hp)*z0
        elif (igrid==3):
            for k in range (0,N):
                z0 = (hc*s[k]+C[k]*hu)/(hc+hu)
                z[:,:,k] = zetau+(zetau+hu)*z0
        elif (igrid==4):
            for k in range (0,N):
                z0 = (hc*s[k]+C[k]*hv)/(hc+hv)
                z[:,:,k] = zetav+(zetav+hv)*z0
        elif (igrid==5):
            for k in range (0,Np):
                z0 = (hc*s[k]+C[k]*hr)/(hc+hr)
                z[:,:,k] = zetar+(zetar+hr)*z0

    return z
    ###################################
"""

 STRETCHING:  Compute ROMS vertical coordinate stretching function

 [s,C]=stretching(Vstretching, theta_s, theta_b, hc, N, kgrid, report)

 Given vertical terrain-following vertical stretching parameters, this
 routine computes the vertical stretching function used in ROMS vertical
 coordinate transformation. Check the following link for details:

    https://www.myroms.org/wiki/index.php/Vertical_S-coordinate

 On Input:

    Vstretching   Vertical stretching function:
                    Vstretching = 1,  original (Song and Haidvogel, 1994)
                    Vstretching = 2,  A. Shchepetkin (UCLA-ROMS, 2005)
                    Vstretching = 3,  R. Geyer BBL refinement
                    Vstretching = 4,  A. Shchepetkin (UCLA-ROMS, 2010)
    theta_s       S-coordinate surface control parameter (scalar)
    theta_b       S-coordinate bottom control parameter (scalar)
    hc            Width (m) of surface or bottom boundary layer in which
                    higher vertical resolution is required during
                    stretching (scalar)
    N             Number of vertical levels (scalar)
    kgrid         Depth grid type logical switch:
                    kgrid = 0,        function at vertical RHO-points
                    kgrid = 1,        function at vertical W-points
 On Output:

    s             S-coordinate independent variable, [-1 <= s <= 0] at
                    vertical RHO- or W-points (vector)
    C             Nondimensional, monotonic, vertical stretching function,
                    C(s), 1D array, [-1 <= C(s) <= 0]

"""
import pylab as pl

def stretching(Vstr, thts, thtb, hc, N, kgrid):
    s=[]
    C=[]

    Np=N+1

    #-----------------------------------------------------------------
    # Compute ROMS S-coordinates vertical stretching function
    #-----------------------------------------------------------------

    # Original vertical stretching function (Song and Haidvogel, 1994).
    if (Vstr == 1):
        ds = 1.0/N

        if (kgrid == 1):
            Nlev = Np
            lev  = np.linspace(0.0,N,Np)
            s    = (lev-N)*ds
        else:
            Nlev = N
            lev  = np.linspace(1.0,N,Np)-0.5
            s    = (lev-N)*ds

        if (thts > 0):
            Ptheta = np.sinh(thts*s)/np.sinh(thts)
            Rtheta = np.tanh(thts*(s+0.5))/(2.0*np.tanh(0.5*thts))-0.5
            C      = (1.0-thtb)*Ptheta+thtb*Rtheta
        else:
            C=s

    # A. Shchepetkin (UCLA-ROMS, 2005) vertical stretching function.
    if (Vstr==2):
        alfa = 1.0
        beta = 1.0
        ds   = 1.0/N

        if (kgrid == 1):
            Nlev = Np
            lev  = np.linspace(0.0,N,Np)
            s    = (lev-N)*ds
        else:
            Nlev = N
            lev  = np.linspace(1.0,N,Np)-0.5
            s    = (lev-N)*ds

        if (thts > 0):
            Csur = (1.0-np.cosh(thts*s))/(np.cosh(thts)-1.0)
            if (thtb > 0):
                Cbot   = -1.0+np.sinh(thtb*(s+1.0))/np.sinh(thtb)
                weigth = (s+1.0)**alfa*(1.0+(alfa/beta)*(1.0-(s+1.0)**beta))
                C      = weigth*Csur+(1.0-weigth)*Cbot
            else:
                C=Csur
        else:
            C=s

    # R. Geyer BBL vertical stretching function.
    if (Vstr==3):
        ds   = 1.0/N

        if (kgrid == 1):
            Nlev = Np
            lev  = np.linspace(0.0,N,Np)
            s    = (lev-N)*ds
        else:
            Nlev = N
            lev  = np.linspace(1.0,N,Np)-0.5
            s    = (lev-N)*ds

        if (thts > 0):
            exp_s = thts   # surface stretching exponent
            exp_b = thtb   # bottom  stretching exponent
            alpha = 3      # scale factor for all hyperbolic functions
            Cbot  = np.log(np.cosh(alpha*(s+1.0)**exp_b))/np.log(np.cosh(alpha))-1.0
            Csur  = -np.log(cosh(alpha*abs(s)**exp_s))/log(cosh(alpha))
            weight= (1-np.tanh( alpha*(s+0.5)))/2.0
            C     = weight*Cbot+(1.0-weight)*Csur
        else:
            C=s

    # A. Shchepetkin (UCLA-ROMS, 2010) double vertical stretching function
    # with bottom refinement
    if (Vstr == 4):
        ds   = 1.0/N

        if (kgrid == 1):
            Nlev = Np
            lev  = np.linspace(0.0,N,Np)
            s    = (lev-N)*ds
        else:
            Nlev = N
            lev  = np.linspace(1.0,N,Np)-0.5
            s    = (lev-N)*ds

        if (thts > 0):
            Csur = (1.0-np.cosh(thts*s))/(np.cosh(thts)-1.0)
        else:
            Csur = -s**2

        if (thtb > 0):
            Cbot = (np.exp(thtb*Csur)-1.0)/(1.0-np.exp(-thtb))
            C    = Cbot
        else:
            C    = Csur

    return (s,C)

