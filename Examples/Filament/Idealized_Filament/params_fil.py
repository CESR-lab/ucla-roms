'''
Parameters for grid, initial and boundary conditions
'''
sim_name = 'Fil1'
###########################
#Parameters
##########################
#Grid Paramters
dx = 100     #meters
dy = 100     #meters
Lx = 20200   #meters 
Ly = 1200    # meters
N  = 250     #grid levels (sigma)
theta_s = 6.
theta_b = 6.
hc = 25.

#Bathymetry Paramters
H = 250 #total tdepth

#below 3 don't matter for flat bottom
# but leave in in case want to do slope
s_h = 0 #bathymetric slope
Ld = 8e3 #length of flat shelf
hmin = 2. #minimum depth

#Other grid parameters
f0 = 7.81e-5 #coriolis parameters
grid_ang = 0. #grid rotation (default=0-->xi=east)

# surface stress
sustr0 = 5.0e-6    # 5e-6 produces 0.01 surface u current...
svstr0 = 0.0


####################################
#Initial Condition Parameters
'''
version of idealized filament
from  McWilliams 2017 (JFM)
'''
####################################
#Density structure
#structure = 'front'
structure = 'fil'
#structure = None


#Choose initial velocity
# either TTW or geostrophic
flow = 'TTW'
#flow = 'geostrophic'
#flow=None


##########################
#buoyancy constant (ms^-2)
b0 = 6.4e-3 
#Fractional reduction of surface stratification relative to interior
B =0.025
#Vertical scale transition from surface to interior (m)
lambda_inv = 10 #3 
#Background stratification (s^-2)
Nb =  1e-7
#interior stratification (s^-2)
N0 = 3.4e-5 

#half width of filament
L_fil = 2000#750. 
#base surface boundary layer depth (m)
h0 = 60
#deviation of surface boundary layer at front or filament (m)
dh = 15

#Vertical Mixing Formulation
sig0 = 5e-3
K0 = 0.01
K_bak = 1e-4 #background diffusivity

#Constant Salinity
S0 = 32
