'''
Create ROMS grid and initial conditions for 
an idealized problem

call the script like >> python setup_grid_init_bry.py params_file.py

where params_file.py is a separate script that initializes
parameters for grid, initial, and boundary conditions
'''
########################################
import os
import sys
from netCDF4 import Dataset as netcdf
import numpy as np
import ROMS_depths as RD
from datetime import date
import TTW_momentum as TTW_mom
import ROMS_tools as RT
from scipy import special
######################################

param_file = sys.argv[1]
#execfile('./'+param_file)
exec(open('./'+param_file).read())  # python3 version


################
#Horizontal grid
################
x = np.arange(-Lx/2,Lx/2,dx)
y = np.arange(-Ly/2,Ly/2,dy)
nx = len(x)
ny = len(y)
pm = np.zeros([ny,nx])
pn = np.zeros([ny,nx])
angle = np.zeros([ny,nx]) #no rotation
mask_rho = np.zeros([ny,nx]) #all water
mask_rho[:,:] = 1.
pm[:,:] = 1./dx
pn[:,:] = 1./dy
x_rho, y_rho = np.meshgrid(x,y)
############################
#Create idealized bathmetry
#############################
h = np.zeros([ny,nx])
i_Ld = abs(x-Ld).argmin()
for j in range(ny):
   h[j,0:i_Ld] = H
   h[j,i_Ld::] = H - s_h*(x[i_Ld::]-Ld) 
h[h<hmin] = hmin
##############################################

#########################################
# Write and save variables to netcdf file
#########################################
grd_nc = netcdf(sim_name + '_grd.nc', 'w')
print ('')
print ('Creating ROMS grid: ' + sim_name + '_grd.nc')
print ('')

# SET GLOBAL ATTRIBUTES
grd_nc.title = 'ROMS grid produced by setup_grid_init_bry.py'
grd_nc.date = date.today().strftime("%B %d, %Y")

#SET DIMENSIONS
grd_nc.createDimension('xi_rho', nx)
grd_nc.createDimension('xi_u', nx-1)
grd_nc.createDimension('eta_rho', ny)
grd_nc.createDimension('eta_v', ny-1)

#CREATE AND WRITE VARIABLES

spherical = grd_nc.createVariable('spherical', 'c')
setattr(spherical, 'long_name', 'Grid type logical switch')
spherical[0]='F'



#xl_nc = grd_nc.createVariable('xl','f4', ('eta_rho', 'xi_rho'))
xl_nc = grd_nc.createVariable('xl','d')
setattr(xl_nc, 'long_name', 'domain length xi-direction')
setattr(xl_nc, 'units', 'meter')
xl_nc[0]=Lx


el_nc = grd_nc.createVariable('el','d')
setattr(el_nc, 'long_name', 'domain length in eta-direction')
setattr(el_nc, 'units', 'm')
el_nc[0]=Ly


xr_nc = grd_nc.createVariable('x_rho','f4', ('eta_rho', 'xi_rho'))
setattr(xr_nc, 'long_name', 'physical dimension in xi-direction')
setattr(xr_nc, 'units', 'm')
xr_nc[:,:] = x_rho


yr_nc = grd_nc.createVariable('y_rho','f4', ('eta_rho', 'xi_rho'))
setattr(yr_nc, 'long_name', 'physical dimension in eta-direction')
setattr(yr_nc, 'units', 'm')
yr_nc[:,:] = y_rho

h_nc = grd_nc.createVariable('h','f4', ('eta_rho', 'xi_rho'))
setattr(h_nc, 'long_name', 'bottom topography')
setattr(h_nc, 'units', 'm')
h_nc[:,:] = h

pm_nc = grd_nc.createVariable('pm','f4', ('eta_rho', 'xi_rho'))
setattr(pm_nc, 'long_name', 'curvilinear coordinate metric in XI-direction')
setattr(pm_nc, 'units', '1/m')
pm_nc[:,:] = pm

pn_nc = grd_nc.createVariable('pn','f4', ('eta_rho', 'xi_rho'))
setattr(pn_nc, 'long_name', 'curvilinear coordinate metric in ETA-direction')
setattr(pn_nc, 'units', '1/m')
pn_nc[:,:] = pn

f_nc = grd_nc.createVariable('f','f4', ('eta_rho', 'xi_rho'))
setattr(f_nc, 'long_name', 'Coriolis parameter at RHO-points')
setattr(f_nc, 'units', '1/s')
f_nc[:,:] = f0

angle_nc = grd_nc.createVariable('angle','f4', ('eta_rho', 'xi_rho'))
setattr(angle_nc, 'long_name', 'angle between east and XI-directions')
setattr(angle_nc, 'units', 'degrees')
angle_nc[:,:] = grid_ang

mask_nc = grd_nc.createVariable('mask_rho','f4', ('eta_rho', 'xi_rho'))
setattr(mask_nc, 'long_name', 'land/sea mask at rho-points')
setattr(mask_nc, 'units', 'land/water (0/1)')
mask_nc[:,:] = mask_rho

grd_nc.close()

print ('Saved grid file')


################################################
#Initial Condition
###############################################

print (' ')
print ('         Computing Initial Condition ')

###################
#Make vertical grid
'''
Create approximate sigma level 
grid to mimic ROMS functionality
'''
####################
zeta = np.zeros([ny,nx])
Vtrans = 2
Vstret = 4
z_r = RD.set_depth(Vtrans,Vstret,theta_s,theta_b,hc,N,1,h,zeta)
z_w = RD.set_depth(Vtrans,Vstret,theta_s,theta_b,hc,N,5,h,zeta)
s,Cs_r = RD.stretching(Vstret,theta_s, theta_b, hc,N,0)
s,Cs_w = RD.stretching(Vstret,theta_s, theta_b, hc,N,1)
#Cs_r = Cs_r_temp[0:N]
################################################
#Idealized initial buoyancy and mixing profile
##############################################
h_sbl = np.zeros([ny,nx])
for i in range(nx):
    if structure=='fil':
       h_sbl[:,i] = h0 + dh * np.exp( -((x[i]) / L_fil)**2)
    elif structure=='front':
       h_sbl[:,i] = h0 - ((dh/2.) * special.erf(x[i] / L_fil)) 
    else:
        h_sbl[:,i] = h0


b = np.zeros([N,ny,nx])
for i in range(nx):
    for j in range(ny):
        b[:,j,i] = b0 + Nb * (z_r[j,i,:] + H) + (0.5*N0) * ( (1+ B) * z_r[j,i,:] - ( 1- B) *( h_sbl[j,i] + lambda_inv * np.log(np.cosh((1./lambda_inv) *(z_r[j,i,:] + h_sbl[j,i])))))
       
###################
#Make initial AKv
##################
Kv = np.zeros([N+1,ny,nx])
#Set background first
Kv[:,:,:] = K_bak

Km = np.zeros([ny,nx])
Km0 = (K0 - K_bak)  / (1 + dh/h0)
for i in range(nx):
    for j in range(ny):
        Km[j,i] = Km0 * h_sbl[j,i] / h0

ff = 4./27 * (1 + sig0)**2
for i in range(nx):
    for j in range(ny):
        for k in range(1,N):
            sig = -z_r[j,i,k] / h_sbl[j,i]
            if sig<=1:
               Kv[k,j,i] = Kv[k,j,i] + Km[j,i] * (sig + sig0) * (1-sig)**2/ff

######################################
#Initialize temperature, not buoyancy
'''
T = b / (alpha * g)
where alpha = thermal expansion coefficient
'''
####################################
temp = b / (2e-4 * 9.81)




#####################
#Flow initilization
####################
u = np.zeros([N,ny,nx-1])
v = np.zeros([N,ny-1,nx])
ubar = np.zeros([ny,nx-1])
vbar = np.zeros([ny-1,nx])
zeta = np.zeros([ny,nx])


###################################
#Solve TTW velocities
###################################
#Compute db/dx, db/dy on w-levels
#bx = np.gradient(b,axis=
bx = np.gradient(b,axis=2) / dx
by = np.gradient(b,axis=1) / dy
zr_swap = np.swapaxes(z_r.T,1,2)
zw_swap = np.swapaxes(z_w.T,1,2)

#Compute geostrophic along-filament flow
import scipy.integrate as integrate
bx_w = RT.rho2w(bx.T,zr_swap.T,zw_swap.T)
zw_in=np.swapaxes(z_w,0,1)
#Along-filament geostrophic velocity
vg = (integrate.cumtrapz(bx_w,zw_in,axis=2).T/f0).T  # DPD to me v_g is at w points?
vg_kji = np.swapaxes(vg,0,2)
#Compute zeta from this
#import matplotlib.pyplot as plt
#Integrate vg(z=0) in x to get free-surface
zeta_g = np.zeros([ny,nx])
for i in range(1,nx):       # DPD vg and zeta_g at the same horizontal coordinate, so i and i-1 create an
    cff = f0*dx / 9.81      # erroneous 1/2 grid cell shift
    zeta_g[:,i] = zeta_g[:,i-1] + (0.5 * (cff * vg_kji[-1,:,i] + cff * vg_kji[-1,:,i-1]) )
    #(f0 * vg_kji[-1,:,i] * dx/9.81)
     

############################################
#Now compute pressure and pressure gradients
###########################################
phi = np.zeros([N,ny,nx])
b_w = RT.rho2w(b.T,zr_swap.T,zw_swap.T).T
zw_zeta = np.copy(zw_swap)
zw_zeta[-1,:,:] = zw_swap[-1,:,:] + zeta_g
Hz = zw_zeta[1:,:,:] - zw_zeta[:-1,:,:]
phi = np.zeros([N,ny,nx])
phi[-1,:,:] = 9.81 * zeta_g
for k in range(N-1,0,-1):
    phi[k-1,:,:] = phi[k,:,:] + (0.5 * ((Hz[k,:,:] * (-b[k,:,:])) + (Hz[k-1,:,:] * (-b[k-1,:,:]))))


#Pressure gradients
dphidx = np.gradient(phi,axis=2) / dx
dphidy = np.gradient(phi,axis=1) / dy
##########################################


#SET UP ARRAYS FOR TTW module (need to swap axes)
#TTW solver takes conventions [i,j,k]
f_ij = np.zeros([nx,ny])
f_ij[:,:] = f0
sustr = np.zeros([nx,ny])
svstr = np.zeros([nx,ny])
Kv_swap =np.swapaxes(Kv,0,2) 
zw_in=np.swapaxes(z_w,0,1)
zr_in = np.swapaxes(z_r,0,1)
bot_stress_opt = -1 #no bottom stress

#Solve TTW system and also get geostrophic velocities
ut,vt,ug,vg = TTW_mom.solve_ttw_momentum_sig(dphidx.T,dphidy.T,Kv_swap,sustr,svstr,sustr,svstr,f_ij,zr_in,zw_in,bot_stress_opt)

#Put u at horizontal u-points
if flow =='TTW':
   u_in = ut
   v_in = vt
elif flow == 'geostrophic':
   u_in = ug
   v_in = vg
else:
    u_in = np.zeros([nx,ny,N])
    v_in = np.zeros([nx,ny,N])

#Put u at horizontal u-pionts
u_init = RT.rho2u(np.swapaxes(u_in,0,2))
#Put v at horizontal v-points
v_init = RT.rho2v(np.swapaxes(v_in,0,2))


##########################
#Write Variables
#########################
init_nc = netcdf(sim_name + '_init.nc', 'w')
print ('')
print ('Creating initial condition: ' + sim_name + '_init.nc')
print ('')


# SET GLOBAL ATTRIBUTES
init_nc.title = 'ROMS initial conditon produced by setup_grid_init_bry.py'
init_nc.date = date.today().strftime("%B %d, %Y")
init_nc.VertCoordType='SM09'
init_nc.theta_s = theta_s
init_nc.theta_b = theta_b
init_nc.hc = hc
init_nc.Cs_r = Cs_r
init_nc.Cs_w = Cs_w

#SET DIMENSIONS
init_nc.createDimension('xi_rho', nx)
init_nc.createDimension('xi_u', nx-1)
init_nc.createDimension('eta_rho', ny)
init_nc.createDimension('eta_v', ny-1)
init_nc.createDimension('s_rho', N)
init_nc.createDimension('s_w', N+1)



#CREATE AND WRITE VARIABLES
otime_nc = init_nc.createVariable('ocean_time', 'd')
setattr(otime_nc, 'long_name', 'averaged time since initialization')
setattr(otime_nc, 'units', 'second')
otime_nc[0] = 0.

temp_nc = init_nc.createVariable('temp','f4', ('s_rho','eta_rho', 'xi_rho'))
setattr(temp_nc, 'long_name', 'potential temperature')
setattr(temp_nc, 'units', 'deg C')
temp_nc[:,:,:] = temp

salt_nc = init_nc.createVariable('salt','f4', ('s_rho','eta_rho', 'xi_rho'))
setattr(salt_nc, 'long_name', 'salinity')
setattr(salt_nc, 'units', 'psu')
salt_nc[:,:,:] = S0

AKv_nc = init_nc.createVariable('AKv','f4', ('s_w','eta_rho', 'xi_rho'))
setattr(AKv_nc, 'long_name', 'vertical mixing')
setattr(AKv_nc, 'units', 'm^2 / s')
AKv_nc[:,:,:] = Kv



u_nc = init_nc.createVariable('u','f4', ('s_rho','eta_rho', 'xi_u'))
setattr(u_nc, 'long_name', 'XI-velocity component')
setattr(u_nc, 'units', 'm/s')
u_nc[:,:,:] = u_init

ubar_nc = init_nc.createVariable('ubar','f4', ('eta_rho', 'xi_u'))
setattr(ubar_nc, 'long_name', 'barotropic XI-velocity component')
setattr(ubar_nc, 'units', 'm/s')
ubar_nc[:,:] = np.nanmean(u_init,axis=0)  # DPD this looks wrong, can't just mean since vertical grid non-uniform

v_nc = init_nc.createVariable('v','f4', ('s_rho','eta_v', 'xi_rho'))
setattr(v_nc, 'long_name', 'ETA-velocity component')
setattr(v_nc, 'units', 'm/s')
v_nc[:,:,:] = v_init

vbar_nc = init_nc.createVariable('vbar','f4', ('eta_v', 'xi_rho'))
setattr(vbar_nc, 'long_name', 'barotropic ETA-velocity component')
setattr(vbar_nc, 'units', 'm/s')
vbar_nc[:,:] = np.nanmean(v_init,axis=0)  # DPD this looks wrong, can't just mean since vertical grid non-uniform


zeta_nc = init_nc.createVariable('zeta','f4', ('eta_rho', 'xi_rho'))
setattr(zeta_nc, 'long_name', 'free-surface elevation')
setattr(zeta_nc, 'units', 'm')
zeta_nc[:,:] = zeta_g


hbls_nc = init_nc.createVariable('hbls','f4', ('eta_rho', 'xi_rho'))
setattr(hbls_nc, 'long_name', 'surface boundary layer depth')
setattr(hbls_nc, 'units', 'm')
hbls_nc[:,:] = h_sbl


init_nc.close()
###############################################

