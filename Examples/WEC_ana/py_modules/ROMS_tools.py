
#  Adapted from Jon Gula's code
#Daniel Dauhajre, UCLA 2013
#module containing various roms tools(no depths in here, depths.py is module for that)

##########################################################################################3

#load necessary modules
from pylab import *
from scipy.io.netcdf import netcdf_file as netcdf
import numpy as np
import os
def change_dir_py(dir_name):
    if not os.path.exists(dir_name):
       print('Creating directory: ' + dir_name)
       os.makedirs(dir_name)
    print('Moving into directory: ' + dir_name)
    os.chdir(dir_name)
    ####
def mask_2D_var(var_in,mask_rho):
    '''
    Multiply 2D variable by mask
    to nan-out masked regions
    '''
    return var_in*mask_rho
    #######################################


def mask_3D_var(var_in,mask_rho):
    '''
    Multiply 3D variable by mask
    to nan-out masked regions
    '''
    var_masked = np.zeros(var_in.shape)
    for k in range(var_in.shape[0]):
        var_masked[k,:,:] = var_in[k,:,:] * mask_rho
    return var_masked
    #######################################

def vort_2D(u,v,pm,pn, f):
    return psi2rho(curvilinear_vort(u,v,pm,pn)) / f
    ###################################

def div_2D(u,v,pm,pn, f):
    return curvilinear_div(u,v,pm,pn) / f
    ###################################


def vort_3D(u,v,pm,pn, f):
    [nz,ny,nx_u] = u.shape
    vort_out = np.zeros([nz,ny,nx_u+1])
    for k in range(nz):
        vort_out[k,:,:] = psi2rho(curvilinear_vort(u[k,:,:],v[k,:,:],pm,pn)) / f
    return vort_out
    ###################################

def div_3D(u,v,pm,pn, f):
    [nz,ny,nx_u] = u.shape
    div_out = np.zeros([nz,ny,nx_u+1])
    for k in range(nz):
        div_out[k,:,:] = curvilinear_div(u[k,:,:],v[k,:,:],pm,pn) / f
    return div_out
    ###################################

def curvilinear_div(ubar,vbar,pm,pn):
    """ Calculate divergence  of velocity
        in ROMS curvilinear coordinates

    div = (pm*pn) * ( d/dxi((1/pn)*u) + d/deta( (1/pm) * v)
      """
    # GET SIZE
    [neta,nxi] = ubar.shape
    div_horiz = np.zeros([neta,nxi])
    A = (pm * pn)
    temp1,dxi_div = np.gradient( (1/pn) * u2rho(ubar))
    deta_div,temp2 = np.gradient( (1/pm) * v2rho(vbar))
    div_horiz = A * (dxi_div + deta_div)
    return div_horiz


def curvilinear_grad(phi,pm,pn):
    """ Compute horizontal gradient
        of a scalar field in curvilinear coordinates
    
    grad(phi) = pm * dphi/dxi + pn * dphi/deta
   
    *** USES SIMPLE CENTERED DIFFERENCING ***

    returns 2-component gradient vector in [eta,xi] shape
    """
    [neta,nxi] = phi.shape

    dphi_e,dphi_x = np.gradient(phi)

    dphi_deta = dphi_e * pm
    dphi_dxi = dphi_x * pn

    return dphi_deta,dphi_dxi


def princ_strain_axis(ubar,vbar,pm,pn): 
    '''
    CALCULATE PRINCIPLE STRIAN AXIS
    tan(2*theta_p) = (v_x + u_y) / (u_x - v_y)
    '''

    v_eta, v_xi = curvilinear_grad(v2rho(vbar), pm, pn)
    u_eta, u_xi = curvilinear_grad(u2rho(ubar), pm, pn)

    theta_p = 0.5 * np.degrees(np.arctan2( (v_xi + u_eta) , (u_xi - v_eta)))

    return theta_p

def curvilinear_vort(ubar,vbar,pm,pn):
    """ Calculate vertical component of velocity curl
        in curvilinear coordinates
    
     
    vort = (pm*pn) * (d/dxi(((1/pn)*v) - d/deta((1/pm)*u))
             = (pm*pn) * (d1 - d2)
    """

    
    [neta,nxi] = pm.shape
   
    d1 = np.zeros([neta-1,nxi-1])
    d2 = np.zeros([neta-1,nxi-1])
    d1[:,:] = (rho2v(1/pn)[:,1:] *vbar[:,1:]) - (rho2v(1/pn)[:,:-1] * vbar[:,:-1])
    d2[:,:] = (rho2u(1/pm)[1:,:] *ubar[1:,:]) - (rho2u(1/pm)[:-1,:] * ubar[:-1,:])

    xi = rho2psi(pm*pn) * (d1 - d2)
    return xi
    

def curvilinear_div_strain(ubar,vbar,pm,pn):
    """ Calculate divergence and strain of velocity
        in ROMS curvilinear coordinates 
    
    div = (pm*pn) * ( d/dxi((1/pn)*u) + d/deta( (1/pm) * v) 
          
    S = \sqrt{ (pm*d/dxi(u) - pn * d/deta(v))^2 + (pm*d/dxi(v) + pn*d/deta(u))**2)
    """

    # GET SIZE
    #[neta,nxi] = ubar.shape
    [neta,nxi] = pm.shape

    div_horiz = np.zeros([neta,nxi])
    strain    = np.zeros([neta,nxi])


    A = (pm * pn)

    temp1,dxi_div = np.gradient( (1/pn) * u2rho(ubar))
    deta_div,temp2 = np.gradient( (1/pm) * v2rho(vbar))

    div_horiz = A * (dxi_div + deta_div)


    du_deta,du_dxi = np.gradient(u2rho(ubar))
    dv_deta,dv_dxi = np.gradient(v2rho(vbar))

    strain = np.sqrt( ((pm*du_dxi) - (pn * dv_deta))**2 + ((pm*dv_dxi) + (pn*du_deta))**2)

    return div_horiz,strain










     

##################################################
# vort_calc ==> compute relative vorticity 
############################################

'''
INPUTS
ubar - depth averaged U-velocity at single time step [eta_rho, xi_rho] shape
vbar - depth averaged V-velocity at single time step [eta_rho, xi_rho] shape
pm   - from grid file
pn   - form grid file
OUTPUTS
xi - vorticity with [eta_rho, xi_rho] shape, normalize by f if necessary outside of module

'''
def vort_calc(ubar, vbar, pm, pn):
    [Mp,Lp] = pm.shape
    L       = Lp - 1
    M       = Mp - 1
    Lm      = L - 1
    xi      = np.zeros([M,L])
    mn_p    = np.zeros([M,L])
    uom     = np.zeros([M,Lp])
    von     = np.zeros([Mp,L])
    
    
    uom = 2 * ubar / (pm[:,0:L] + pm[:,1:Lp])
    von = 2 * vbar / (pn[0:M,:] + pn[1:Mp,:])

    
   # uom = 2 * ubar / (pn[0:M,:] + pn[1:Mp,:])
    #von = 2* vbar /  (pm[:,0:L] + pm[:,1:Lp])

    mn   = pm * pn
    mn_p = ( mn[0:M,0:L] + mn[0:M,1:Lp] + mn[1:Mp, 1:Lp] + mn[1:Mp, 0:L]) / 4
    xi   = mn_p * ( von[:,1:Lp] - von[:,0:L]  - uom[1:Mp,:] + uom[0:M,:] )
    return xi

    

'''
matlab code
[Mp,Lp]=size(pm);
L=Lp-1;
M=Mp-1;
Lm=L-1;
Mm=M-1;
xi=zeros(M,L);
mn_p=zeros(M,L);
uom=zeros(M,Lp);
von=zeros(Mp,L);
uom=2*ubar./(pm(:,1:L)+pm(:,2:Lp));
von=2*vbar./(pn(1:M,:)+pn(2:Mp,:));
mn=pm.*pn;
mn_p=(mn(1:M,1:L)+mn(1:M,2:Lp)+...
      mn(2:Mp,2:Lp)+mn(2:Mp,1:L))/4;
xi=mn_p.*(von(:,2:Lp)-von(:,1:L)-uom(2:Mp,:)+uom(1:M,:));
'''







def div_strain_calc(ubar, vbar, pm, pn):
  
    [Mp,Lp] = pm.shape
    M = Mp-1
    L = Lp-1
    
    
    #ux,uy = np.gradient(ubar) * ( (rho2u(pm) + rho2u(pn)) / 2.)
    #vx,vy = np.gradient(vbar) * ( (rho2v(pm) + rho2v(pn)) / 2.)
    
    '''
    ASSUMES PYTHON CONVENTIONS OF READING IN ROMS NETCDF
    ubar.shape, vbar.shape = [eta, xi] = [y, x] 
    '''

    uy,ux = np.gradient(ubar) * ( (rho2u(pm) + rho2u(pn)) / 2.)
    vy,vx = np.gradient(vbar) * ( (rho2v(pm) + rho2v(pn)) / 2.)
    
    ux_r = u2rho(ux)
    uy_r = u2rho(uy)
    vx_r = v2rho(vx) 
    vy_r = v2rho(vy)


    div    = ux_r + vy_r
    strain = np.sqrt( (ux_r - vy_r) **2 + (vx_r + uy_r)**2)
    
    return div,strain






#######################################################
#Transfert a field at psi points to rho points
#######################################################

def psi2rho(var_psi):

    if np.ndim(var_psi)<3:
        var_rho = psi2rho_2d(var_psi)
    else:
        var_rho = psi2rho_3d(var_psi)

    return var_rho


##############################

def psi2rho_2d(var_psi):

    [M,L]=var_psi.shape
    Mp=M+1
    Lp=L+1
    Mm=M-1
    Lm=L-1

    var_rho=np.zeros((Mp,Lp))
    var_rho[1:M,1:L]=0.25*(var_psi[0:Mm,0:Lm]+var_psi[0:Mm,1:L]+var_psi[1:M,0:Lm]+var_psi[1:M,1:L])
    var_rho[0,:]=var_rho[1,:]
    var_rho[Mp-1,:]=var_rho[M-1,:]
    var_rho[:,0]=var_rho[:,1]
    var_rho[:,Lp-1]=var_rho[:,L-1]

    return var_rho

#############################

def psi2rho_3d(var_psi):


    [Nz,Mz,Lz]=var_psi.shape
    var_rho=np.zeros((Nz,Mz+1,Lz+1))

    for iz in range(0, Nz, 1):    
        var_rho[iz,:,:]=psi2rho_2d(var_psi[iz,:,:])


    return var_rho

#######################################################
#Transfert a field at rho points to psi points
#######################################################

def rho2psi(var_rho):

    if np.ndim(var_rho)<3:
        var_psi = rho2psi_2d(var_rho)
    else:
        var_psi = rho2psi_3d(var_rho)

    return var_psi


##############################

def rho2psi_2d(var_rho):

    var_psi = 0.25*(var_rho[1:,1:]+var_rho[1:,:-1]+var_rho[:-1,:-1]+var_rho[:-1,1:])

    return var_psi

#############################

def rho2psi_3d(var_rho):

    var_psi = 0.25*(var_rho[:,1:,1:]+var_rho[:,1:,:-1]+var_rho[:,:-1,:-1]+var_rho[:,:-1,1:])

    return var_psi


#######################################################
#Transfert a 2 or 3-D field at rho points to u points
#######################################################

def rho2u(var_rho):

    if np.ndim(var_rho)<3:
        var_u = rho2u_2d(var_rho)
    else:
        var_u = rho2u_3d(var_rho)

    return var_u

def rho2u_2d(var_rho):

    [Mp,Lp]=var_rho.shape
    L=Lp-1
    var_u=0.5*(var_rho[:,0:L]+var_rho[:,1:Lp])

    return var_u


def rho2u_3d(var_rho):

    [N,Mp,Lp]=var_rho.shape
    L=Lp-1
    var_u=0.5*(var_rho[:,:,0:L]+var_rho[:,:,1:Lp])

    return var_u

#######################################################
#Transfert a 3-D field at rho points to v points
#######################################################

def rho2v(var_rho):

    if np.ndim(var_rho)<3:
        var_v = rho2v_2d(var_rho)
    else:
        var_v = rho2v_3d(var_rho)

    return var_v

#######################################################

def rho2v_2d(var_rho):

    [Mp,Lp]=var_rho.shape
    M=Mp-1
    var_v=0.5*(var_rho[0:M,:]+var_rho[1:Mp,:]);

    return var_v

#######################################################

def rho2v_3d(var_rho):

    [N,Mp,Lp]=var_rho.shape
    M=Mp-1
    var_v=0.5*(var_rho[:,0:M,:]+var_rho[:,1:Mp,:]);

    return var_v

#############################################################
#Transfer a field at vertical rho-points to vertical w-points
'''
USE np.interp to interpolate instead of simple averaging
because of uneven sigma level spacing (dz)
'''

#############################################################
def rho2w(var_rho,z_rho,z_w):
    if np.ndim(var_rho) <3:
       var_w = rho2w_2d(var_rho,z_rho,z_w)
    else:
       var_w = rho2w_3d(var_rho,z_rho,z_w)

    return var_w


def rho2w_2d(var_rho,z_rho,z_w):
    [M,N] = var_rho.shape
    var_w = np.zeros([M,N+1])
    for j in range(M):
        var_w[j,:] = np.interp(z_w[j,:],z_rho[j,:],var_rho[j,:])
    return var_w
    
      
   
def rho2w_3d(var_rho,z_rho,z_w):
    [M,L,N] = var_rho.shape
    var_w = np.zeros([M,L,N+1])
    for j in range(M):
        for i in range(L):
            var_w[j,i,:] = np.interp(z_w[j,i,:],z_rho[j,i,:],var_rho[j,i,:]) 
    return var_w

###############################################################
# TRANSFER A FIELD AT VERTICAL W-POINTS TO VERTICAL RHO-POINTS
##############################################################
def w2rho(var_w,z_rho,z_w):
    if np.ndim(var_w) < 3:
       var_rho = w2rho_2d(var_w,z_rho,z_w)
    else:
        var_rho = w2rho_3d(var_w,z_rho,z_w)
    
    return var_rho


def w2rho_2d(var_w,z_rho,z_w):
    [M,N_w] = var_w.shape
    var_rho = np.zeros([M,N_w-1])
    for j in range(M):
        var_rho[j,:] = np.interp(z_rho[j,:],z_w[j,:],var_w[j,:])
    return var_rho

def w2rho_3d(var_w,z_rho,z_w):
    [M,L,N_w] = var_w.shape
    var_rho = np.zeros([M,L,N_w-1])
    for j in range(M):
        for i in range(L):
            var_rho[j,i,:] = np.interp(z_rho[j,i,:],z_w[j,i,:],var_w[j,i,:])

    return var_rho


#######################################################
#Transfert a 2-D field at u points to the rho points
#######################################################

def u2rho(var_u):


    if np.ndim(var_u)<3:
        var_rho = u2rho_2d(var_u)
    else:
        var_rho = u2rho_3d(var_u)

    return var_rho

#######################################################

def u2rho_2d(var_u):

    [Mp,L]=var_u.shape
    Lp=L+1
    Lm=L-1
    var_rho=np.zeros((Mp,Lp))
    var_rho[:,1:L]=0.5*(var_u[:,0:Lm]+var_u[:,1:L])
    var_rho[:,0]=var_rho[:,1]
    var_rho[:,Lp-1]=var_rho[:,L-1]
    return var_rho

#######################################################

def u2rho_3d(var_u):

    [N,Mp,L]=var_u.shape
    Lp=L+1
    Lm=L-1
    var_rho=np.zeros((N,Mp,Lp))
    var_rho[:,:,1:L]=0.5*(var_u[:,:,0:Lm]+var_u[:,:,1:L])
    var_rho[:,:,0]=var_rho[:,:,1]
    var_rho[:,:,Lp-1]=var_rho[:,:,L-1]
    return var_rho


#######################################################
#Transfert a 2 or 2-D field at v points to the rho points
#######################################################

def v2rho(var_v):

    if np.ndim(var_v)<3:
        var_rho = v2rho_2d(var_v)
    else:
        var_rho = v2rho_3d(var_v)

    return var_rho

#######################################################

def v2rho_2d(var_v):

    [M,Lp]=var_v.shape
    Mp=M+1
    Mm=M-1
    var_rho=np.zeros((Mp,Lp))
    var_rho[1:M,:]=0.5*(var_v[0:Mm,:]+var_v[1:M,:])
    var_rho[0,:]=var_rho[1,:]
    var_rho[Mp-1,:]=var_rho[M-1,:]

    return var_rho

#######################################################

def v2rho_3d(var_v):

    [N,M,Lp]=var_v.shape
    Mp=M+1
    Mm=M-1
    var_rho=np.zeros((N,Mp,Lp))
    var_rho[:,1:M,:]=0.5*(var_v[:,0:Mm,:]+var_v[:,1:M,:])
    var_rho[:,0,:]=var_rho[:,1,:]
    var_rho[:,Mp-1,:]=var_rho[:,M-1,:]

    return var_rho


#######################################################
#interpolate a 3D variable on horizontal levels of constant depths 
#######################################################

def vinterps(var,z,depths,topo, cubic=0,ground_interp=0):

    [N,Mp,Lp]=var.shape
    Nz=len(depths)

    #if var not on rho-grid: interpolate z and topo to the same grid than var (u,v,or psi)
    if var.shape!=z.shape:
        if (var.shape[1]==z.shape[1]-1) and (var.shape[2]==z.shape[2]-1):
            z = rho2psi(z); topo = rho2psi(topo)
        elif (var.shape[1]==z.shape[1]-1):
            z = rho2v(z); topo = rho2v(topo)
        elif (var.shape[2]==z.shape[2]-1):
            z = rho2u(z); topo = rho2u(topo)

    if len(depths)==1:
        vnew=vinterp(var,z,depths[0],topo,cubic,ground_interp)

    else:
        [N,Mp,Lp]=var.shape; Nz=len(depths); vnew=np.zeros((Nz, Mp,Lp))
        for iz in range(0, Nz, 1):
            vnew[iz,:,:]=vinterp(var,z,depths[iz],topo,cubic,ground_interp)

    return vnew

    
####################


def vinterp(var,z,depth,topo=None,cubic=0,ground_interp=0):

    [N,Mp,Lp]=z.shape


    if depth>0: 
    
        varz = np.nan

    #######################################################
    #Simple linear interpolation
    #######################################################

    elif cubic==0:

        levs2=sum(z<depth,0)-1
        levs2[levs2==N-1]=N-2
        levs2[levs2==-1]=0
        levs1=levs2+1

        X,Y=np.meshgrid(np.arange(0,Lp),np.arange(0,Mp))

        pos1=levs1,Y,X
        pos2=levs2,Y,X

        z1=z[pos1]
        z2=z[pos2]

        v1=var[pos1]
        v2=var[pos2]
    
        varz = (((v1-v2)*depth+v2*z1-v1*z2)/(z1-z2))
        if topo!=None: varz[depth<-1*topo]=np.nan

    #######################################################
    #Cubic interpolation (see ShchepetkinMcWilliams08.pdf)
    #######################################################

    elif cubic==1:

        print('cubic interpolation')

        #find the closest level BELOW depth
        levs2=sum(z<depth)-1
        levs1=copy(levs2)
        #levs2[levs1==N-1]=N-2

        #cubic interpolation will use 4 values of var and z in the vertical (2 below and 2 above depth)
        Nlev = 4 

        #prepare arrays for intermediate variables:
        X,Y=np.meshgrid(np.arange(0,Lp),np.arange(0,Mp))
        levs=np.zeros((Nlev,Mp,Lp),int); Xlev=np.zeros((Nlev,Mp,Lp),int); Ylev=np.zeros((Nlev,Mp,Lp),int)
        for ilev in range(Nlev):
            levs[ilev,:,:]=levs2+ilev-1
            Xlev[ilev,:,:]=X
            Ylev[ilev,:,:]=Y


        levs[levs>N-1]=N-1
        levs[levs<0]=0

        pos=levs,Y,X; zz=z[pos]; vark=var[pos]


        #######################################################

        test0=np.zeros((Mp,Lp)); test0[levs2==-1]=1; 
        test1=np.zeros((Mp,Lp)); test1[levs2==0]=1;
        testN1=np.zeros((Mp,Lp)); testN1[levs2==N-2]=1; 
        testN=np.zeros((Mp,Lp)); testN[levs2==N-1]=1;
 
        #######################################################

        zz[1:-1,:,:] = testN * zz[:-2,:,:] + test0 * zz[2:,:,:] + (1 - test0 - testN) * zz[1:-1,:,:]

        dzz = zz[1:,:,:]- zz[:-1,:,:]; 
        dzz[-1,:,:] = testN1 * dzz[1,:,:] + (1-testN1)* dzz[-1,:,:]
        dzz[0,:,:] = test1 * dzz[1,:,:] + (1-test1)* dzz[0,:,:]

        vark[1:-1,:,:] = testN * vark[:-2,:,:] + test0 * vark[2:,:,:] + (1 - test0 - testN) * vark[1:-1,:,:]

        dvark = vark[1:,:,:]-vark[:-1,:,:]; 
        dvark[-1,:,:] = testN1 * dvark[1,:,:] + (1-testN1)* dvark[-1,:,:]
        dvark[0,:,:] = test1 * dvark[1,:,:] + (1-test1)* dvark[0,:,:]

        FC0 = (dvark[1:,:,:]+dvark[:-1,:,:])*dzz[1:,:,:]*dzz[:-1,:,:]
        FC1 = (dzz[1:,:,:]+dzz[:-1,:,:])*dvark[1:,:,:]*dvark[:-1,:,:]
        val=dvark[1:,:,:]*dvark[:-1,:,:]

        FC0[val<=0]=1; FC1[val<=0]=0; FC = FC1/FC0


        #######################################################
        
        cff = 1/dzz[1,:,:]; p=depth-zz[1,:,:]; q=zz[2,:,:]-depth

        varz = cff*(q*vark[1,:,:]+p*vark[2,:,:]- (1-test0-testN) * cff*p*q*(cff*(q-p)*dvark[1,:,:]+p*FC[1,:,:]-q*FC[0,:,:]))     

        #######################################################

        # mask values below ground 
        if ground_interp==0: varz[depth<-1*topo]=np.nan

   
    return varz



###############################
# CALCULATE z_r and z_w for entire grid
# INPUTS
# h ---> from grid file, depths
# zeta --> free surface of entire grid
# hc --> global attr from output
# Cs_r --> global attr from output
# Cs_w --> global attr from output

#RETURNS
# z_r --> depths at rho-point vertical levels [N,M,L] N==> levels
# z_w --> depths at w-point vertical levels  [N+1,M,L]
###############################
def zlevs(h,zeta,hc,Cs_r,Cs_w):
    N = Cs_r.shape[0]
    z_r = np.zeros([N,h.shape[0], h.shape[1]])
    
    for k in range(N):
        cff        = hc * ((k+1-N) - 0.5)/N
        z_r[k,:,:] = zeta + (zeta+h) * (cff + Cs_r[k] * h) / (hc + h)

   
        z_w = np.zeros([N+1, h.shape[0], h.shape[1]])
        z_w[0,:,:]  = -h # set bottom to actual depth
        z_w[-1,:,:] = zeta # set top to free surface height 
    
    for k in range(1,N+1):
        cff        = hc * ((k-N)) / N
        z_w[k,:,:] = zeta + (zeta+h) * (cff + Cs_w[k] * h) / (hc+h)

   
    return z_r, z_w



##########################################
# calculate friction velocity based on log law of wall
# INPUTS
# u --> (np array), u velocity at lowest sigma level (MUST BE AT RHO-POINTS)
# v -->  (np array) v-velocity at lowest sigma level (MUST BE AT RHO-POINTS)
# z_sig --> (np array) depths of the lowest sigma level (MUST HAVE SAME DIMENSIONS AS u,v
#           z_sig is the depth at which each velocity measurement is taken
# z0 --> roughness height
# OUTPUTS:
# u_fric --> (np array), u component of friction velocity
# v_fric --> (np array), v component of friction velocity
#########################################
def calc_fric_veloc(u,v,z_sig, z0):
    # von Karman constant
    k = 0.41

    u_fric = (u*k) / (np.log(abs(z_sig)) - log(z0))
    v_fric = (v*k) / (np.log(abs(z_sig)) - log(z0))

    return u_fric, v_fric


