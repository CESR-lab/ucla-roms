
#  Adapted from Jon Gula's code
#Daniel Dauhajre, UCLA 2013
#module containing various roms tools(no depths in here, depths.py is module for that)

##########################################################################################3

#load necessary modules
#from pylab import *
#from scipy.io.netcdf import netcdf_file as netcdf
import numpy as np
import os
#from numba import autojit, prange, jit, njit
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


