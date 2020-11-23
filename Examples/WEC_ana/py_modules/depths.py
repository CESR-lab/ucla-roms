#get_depths ----> get depths of sigma levels for UCLA ROMS output
#Daniel Dauhajre, UCLA AOS
#November 2012
#=============================================
#        INPUTS
#fname ==> string, netcdf file name, path included
#gname ==> string, netcdf grid name, path included
#tindex ==> integer, time index
#ctype ===> string, coordinate type  ('u' or 'v' or 'r')
#scoord ==> string,  S-coordinate
#============================================

#load necessary modules
from pylab import *
#from scipy.io.netcdf import netcdf_file as netcdf
from netCDF4 import Dataset
import numpy as np
import os

def get_depths(fname, gname, tindex, ctype, scoord):
    "This function gets the depths of the sigma levels of a single time step from a given ROMS netcdf file"
#write code to check for input var types for debugging

#read in grid file
    ncgrid=Dataset(gname, 'r')
    h=ncgrid.variables['h'][:,:]
#read in output file
    roms_out=Dataset(fname, 'r')
#get zeta
    zeta_orig=roms_out.variables['zeta'][:,:,:]
    zeta_nosqueeze=np.copy(zeta_orig) #size for L4PV [12, 562, 1602]
    zeta=zeta_nosqueeze[tindex,:,:] #size for L4PV [562, 1602]
    
#get theta b and theta_s
    if scoord=='new' or scoord=='old':
       theta_s=getattr(roms_out, 'theta_s') 
       theta_b=getattr(roms_out, 'theta_b')
       hc=getattr(roms_out, 'hc')
    if len(zeta)==0:
   	  zeta=0.*h
#get s_rho dimension (number of sigma levels)
    dim=roms_out.dimensions
    # added len() 9/24/14 for netcdf4 update
    N=len(dim['s_rho']) #number of sigma levels

   
    vtype=ctype
    if ctype=='u' or ctype=='v':
       vtype='r'
     
    if scoord=='new' or scoord=='old':
       

    #def zlevs(h,zeta,theta_s,theta_b,hc,N,vtype,scoord): THIS SOMEHOW WORKS W/OUT ZLEVS FUNCTION, SLOPPIER, FIX LATER
        alpha=0
        beta=1
   
        [L,M]=h.shape
    
#Set s-curves in domain [-1<sc<0] at vertical W- and RHO-points
    
        if vtype=='w':
           sc=(np.array(range(0,N+1),dtype='float64')-N)/N
           N=N+1
        else:
            sc=(np.array(range(1,N+1),dtype='float64')-N-0.5)/N

    

        if scoord=='new' or scoord=='NEW':
           if theta_s>0:
              csrf=(1.0-np.cosh(theta_s*sc))/(np.cosh(theta_s)-1.0)
           else:
               csrf=-sc**2
        if theta_b>0:
           Cs=(np.exp(theta_b*(csrf+1.0))-1)/(np.exp(theta_b)-1.0)-1.0
        else:
            Cs=csrf
        if theta_s<=0 and theta_b<=0: #uniform spacing
           Cs=sc
        elif scoord=='bot' or scoord=='BOT':
             if theta_b>0:
                 x=sc+1.0
                 wgt=(x**alpha)/beta*(alpha+beta-alpha*(x**beta))
                 csrf=(1.0-np.cosh(theta_s*sc))/(np.cosh(theta_s)-1.0)
                 cbot=np.sinh(theta_b*x)/np.sinh(theta_b)-1.0
                 Cs=wgt*csrf+(1.0-wgt)*cbot
             else:
                 Cs=(1.0-np.cosh(theta_s*sc))/(np.cosh(theta_s)-1.0)
        elif scoord=='old' or scoord=='OLD':
             cff1=1.0/np.sinh(theta_s)
             cff2=0.5/np.tanh(0.5*theta_s)
             Cs=(1.0-theta_b)*cff1*np.sinh(theta_s*sc)+theta_b*(cff2*np.tanh(theta_s*(sc+0.5))-0.5)

#Creat S-coordinate system: based on model topography h(i,j),
# fast-time-averaged free surface field and vertical coordinate
# transformation metrics compute evolving depths of the 3-dimensional
#model grid.


        z=np.zeros([N,L,M])
        if scoord=='old' or scoord=='OLD':
           hinv=1/h
           cff=hc*(sc-Cs)
           cff1=Cs
           cff2=sc+1
           for k in range(0,N):
               z0=cff[k]+cff1[k]*h
               z[k,:,:]=z0+zeta*(1.0+z0*hinv)
        else:
            hinv=1/(h+hc)
            cff=hc*sc
            cff1=Cs
            for k in range(0,N):
                z[k,:,:]=zeta+(zeta+h)*(cff[k]+cff1[k]*h)*hinv
        return [z, Cs]
    #return z


     

   
        


    
