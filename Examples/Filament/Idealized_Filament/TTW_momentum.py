######## SOLVE TTW EQUATION FOR MOMENTUM (NOT SHEAR) #########


#############################################################
# Load necessary modules 
############################################################
import numpy as np
from copy import copy
import time as tm


# UTILIZE SPARSE MATRICES FOR FASTER COMPUTATION
from scipy.sparse import *
from scipy.sparse.linalg import spsolve


def solve_ttw_momentum_sig(px,py,AKv,sustr,svstr,Tbx,Tby,f,z_r,z_w,bot_stress,timing=False,debug=0,ekman=0):
    '''
    INPUTS:
    3D MATRICES:
    [nx,ny,nz] type shape, with nz for AKv being nz+1 (at w-levels, not rho)
    px -(vertical rho-levels) pressure gradient in x-direction (dp/dx)
    py -(vertical rho-levels) pressure gradient in y-direction (dp/dy)
    AKv - (vertical w-levels) vertical viscosity coefficient
    z_w - depths of w-levels (used for vertical discretization)

    2D matrices:
    sustr,svstr --> surface wind stress in x and y-direction
    Tbx,Tby -   --> bottom stress in x and y-direction
    f - coriolis force at i,j points
    pm, pn --> from ROMS grid
    
    
    OUTPUTS:
    u,v --> TTW velocities (at vertical rho-levels)
    ug,vg --> geostrophic velocities (at vertical rho-levels)


    THE TTW SYSTEM IS SOLVED USING A TRIDIAGONAL MATRIX FOR A COUPLED (u,v) system
    
    Solves an Ax=B system where x =  (u)
                                     (v)


    The conventions for the vertical are as follows:
    A[0,0] --> corresponds to bottom (i.e bottom sigma levels)
    A[nz-1,nz-1] --> correspond to surface

    The inputs should follow these conveitons as well!!!!!




    AKv indexed as 1/2...N+1/2 in formal discretization
    Thus AKv[:,:,0] --> AKv_1/2 and AKv[:,:-1] --> AKv_(N+1/2)
    Where N is the number of vertical rho-levels
    '''

    if timing: tstart = tm.time()

    ##################################
    # Create forcing
    ##################################
    new = np.zeros(px.shape)

    # GET DIMENSIONS OF SYSTEM
   
    # nz --> number of rho-levels
    [nx,ny,nz] = px.shape
    
    # vertical difference is of rho-levels, which places the dzs at w-levels 
    # dz vector will be [nx,ny,nz-1] in shape
    # dz should be positive everywhere (given normal rho-level placement)
    dz = z_r[:,:,1:] - z_r[:,:,:-1]
    Hz = z_w[:,:,1:] - z_w[:,:,:-1] 
    # Solutions
    ut,vt = copy(new), copy(new)

    # Create [A] matrix

    # Coupled system, thus A has shape[2*nz,2*nz]
    ndim = 2 *nz
    kv = 2 #2 variable coupled system
    for i in range(nx):
       
        if i%10==0: print ('Solving TTW momentum equation: ', round(100.*i/(nx-1)), ' %')

        for j in range(ny):
        
            A = lil_matrix((ndim,ndim))
            #A = np.zeros([ndim,ndim])
            R = np.zeros(ndim) 
            # Assign bottom and near-bottom
            idx = 0
            ############################################
            C5=AKv[i,j,1]/dz[i,j,0]
            #FILL FIRST ROW
            A[idx,idx+1]  = -f[i,j] * Hz[i,j,0]
            A[idx,idx+kv] = -C5
            A[idx,idx]    = C5
            
            #FILL SECOND ROW
            A[idx+1,idx]      = f[i,j] *Hz[i,j,0]
            A[idx+1,idx+1+kv] = -C5
            A[idx+1,idx+1]    = C5
            #####################
            #     INTERIOR
            ######################
            for k in range(1,nz-1):	    
                ############################################    
                # mapping AKv k-index from rho-level k-index
                #############################################
                # in this loop:
                # if formal indexing is: AKv[k-1/2] --> translates to AKv[k]
                # if formal indexing is: AKv[k+1/2] --> translates to AKv[k+1]
                # where k is the rho-level indexer (loop indexer)

                # dz is separately mapped from Akv since it is differences between  rho-levels
                # formal indexing: dz[k-1/2] --> dz[k-1]
                # formal indexing: dz[k+1/2] --> dz[k]
                
             
                C3 = AKv[i,j,k] / dz[i,j,k-1]
                C4 = AKv[i,j,k+1] / dz[i,j,k]
                C2 = C3 + C4

                idx = 2*k
                
                #FILL FIRST ROW
                A[idx,idx+kv] = -C4
                A[idx,idx]    = C2
                A[idx,idx-kv] = -C3
                A[idx,idx+1]  = -f[i,j] * Hz[i,j,k]

                #FILL SECOND ROW
                A[idx+1,idx+1+kv] = -C4
                A[idx+1,idx+1]    = C2
                A[idx+1,idx+1-kv] = -C3
                A[idx+1,idx]      = f[i,j] *Hz[i,j,k]

                ######################
                #        TOP
                #######################
                idx = (2*nz) - kv
                #C1 = AKv[N-1/2]/dz[N-1/2] ---> AKv[N-1/2] = AKv[nz-1] where nz is # of rho-levels
                #AKV[N-1/2] should be AKv at level below free surface (so it is non-zero)
                C1 = AKv[i,j,nz-1]/dz[i,j,-1]
               
                # TESTING NO STRESS CONDITION W/ GHOST POINTS
                # dz here is actually not physical, but it is a ghost cell so assume equal spacing 
                # as level below
                C0 = AKv[i,j,nz] / dz[i,j,-1]        

                #FILL FIRST ROW
                A[idx,idx-kv] = -C1
                A[idx,idx]    = C1
                A[idx,idx+1]  = -f[i,j] *Hz[i,j,-1]

                #FILL SECOND ROW
                A[idx+1,idx+1-kv] = -C1
                A[idx+1,idx+1]    = C1
                A[idx+1,idx]      = f[i,j]*Hz[i,j,-1]

            # FILL RHS MATRIX	   
            for k in range(nz):
                idx = 2 *k 
                R[idx] =  -px[i,j,k] * Hz[i,j,k]
                R[idx+1] = -py[i,j,k] * Hz[i,j,k]

            # BOUNDARY CONDITIONS
            #idx = (2*N_py) 

            # SURFACE B.C. --> WIND STRESS (m^2/s^2)
            R[-2] =  px[i,j,-1] * Hz[i,j,-1] + sustr[i,j]
            R[-1] = -py[i,j,-1] * Hz[i,j,-1] + svstr[i,j]
            # BOTTOM B.C. IF APPLICABLE, IF NOT R[0], R[1] are -px,-py
            if bot_stress == 1:
                R[0] = -Tbx[i,j] + px[i,j,0] * Hz[i,j,0]
                R[1] = -Tby[i,j] - -py[i,j,0]*Hz[i,j,0]


            if timing: print ('Matrix definition OK.....', tm.time()-tstart)
            if timing: tstart=tm.time()

            #########################################################
            # SOLVE MATRIX A
            #########################################################
            
            #if i == 45 and j == 0:
            #   return A,R
         
            A = A.tocsr()

            if timing: print ('Starting computation......', tm.time() - tstart)
            if timing: tstart = tm.time()
            #X = np.zeros(ndim)
            X =spsolve(A,R)
            #if i == 25 and j == 25:
            #   return A,R,X,i,j
            if timing: print ('computation OK..........', tm.time()-tstart)
            if timing: tstart = tm.time()

            ###########################################################

            # reorder results in (i,j,k)
            for k in range(nz):
                idx = 2 *k
                ut[i,j,k] = X[idx];
                vt[i,j,k] = X[idx+1];

            if timing: print ('allocation OK..........,', tm.time() - tstart)
            if timing: tstart = tm.time()

      

    #######################################################################
    
    # SOLVE FOR GEOSTROPHIC VELOCITIES
    ug,vg = copy(new), copy(new)
    

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                ug[i,j,k] = -py[i,j,k] / f[i,j]
                vg[i,j,k] =  px[i,j,k] / f[i,j]
    
    
    return ut,vt,ug,vg



















