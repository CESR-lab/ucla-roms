import os
import sys
sys.path.append( os.path.join(os.getcwd(),'py_modules') ) # adds relative path to these modules
sys.path.append( os.path.join(os.getcwd(),'py_modules','f2py_modules') ) # adds relative path to these modules
from pylab import *
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
###### DevinD - add for font size of plots
label_size = 6
axis_size = 6
medium_text = 8
noraml_text = 10
######
import ROMS_tools as RT
import seaborn as sns 
plt.ion()

A = 0.001
h=10. # bathymetry [m]
g =9.81
k=2.0*np.pi/1. #lambda=100m
#sigma=np.sqrt(g*k*np.tanh(k*h)) # oldcode
sigma=np.sqrt(g*k) # DevinD - got rid of tanh(k*h) - confirmed made no difference as equals 1.0.  
#print('tanh(k*h)= ', tanh(k*h)) # DevinD - tested = 1.0 so shouldn't have affected result 
#Cg=(sigma/(2.*k))*(1.+(2.*k*h/(np.sinh(2*h*k)))) # eq (5.23) # oldcode - different from roms but gives same answer so changed to below.
Cg = 0.5 * sqrt(g/k)
#Cphi=sigma/k
Clw= np.sqrt(g*h)
print('Cg = ', Cg, 'Clw = ', Clw, '   and u_lw to ust2d ratio: - Clw^2/(Clw^2-Cg^2) = ', -Clw**2/(Clw**2-Cg**2) )


################################
y = 1
# # NewCode # #
filename_results = os.path.join(os.getcwd(),'wpp_his.0000.nc')
roms_out_new = Dataset(filename_results,'r')
new_code_name = 'WEC_analytical' # text to append to filename of plot pdf 
otime    = roms_out_new.variables['ocean_time'][:] # DevinD checked - [0,1,...,200]
#ubar_old  = roms_out.variables['ubar'][:,y,:]
ubar_new  = roms_out_new.variables['ubar'][:,y,:]
#ust2d_old = roms_out.variables['ust2d'][:,y,:]
ust2d_new = roms_out_new.variables['ust2d'][:,y,:] 
#zeta_old  = roms_out.variables['zeta'][:,y,:] #- sup 
zeta_new  = roms_out_new.variables['zeta'][:,y,:]
#print('ust2d_new.shape: = ', ust2d_new.shape, '   zeta_new.shape: ', zeta_new.shape)


env=0.0001
#x = roms_grid.variables['x_rho'][y,:] # oldcode
x_u = roms_out_new.variables['x_rho'][y,:-1] + 0.5
x_rho = roms_out_new.variables['x_rho'][y,:]
print('x_u[0]: ' + str(x_u[0]) + '   x_u[-1]: ' + str(x_u[-1]),' x_u.shape: ', x_u.shape) #DevinD add
print('x_rho[0]: ' + str(x_rho[0]) + '   x_rho[-1]: ' + str(x_rho[-1]),' x_rho.shape: ', x_rho.shape) #DevinD add
pm = roms_out_new.variables['pm'][y,:] # 1/dx
dx = 1/pm[1] # pm = 1/dx
print('pm[0] = ', pm[0], '; dx = ', dx)
################################ 
window_x= 150 # DevinD only used for wave packet plotting 150m either side of centre of wave packet
LLm=pm.shape[0]-2
print('LLm: ' + str(LLm) )
init_shift=-LLm/2
f_u       = np.zeros([len(x_u)])
a2_u     = np.zeros([len(x_u)])
u_lw_a = np.zeros([len(x_u)])
u_lw_rho_a = np.zeros([len(x_rho)])
f_rho       = np.zeros([len(x_rho)])
a2_rho     = np.zeros([len(x_rho)])
#vel_lw_rho = np.zeros([len(x_rho)])
zeta_a  = np.zeros([len(x_rho)]) # DevinD
zeta_hat= np.zeros([len(x_rho)])
zeta_hat_u= np.zeros([len(x_u)]) # set-up at u-point coordinate
ust2d_a     = np.zeros([len(x_u)])
ust2d_rho_a     = np.zeros([len(x_rho)])
#zeta_a_shift = np.zeros([len(x_rho)])

# For Plotting
sns.set()
plt.figure() 

# create file to save error calculation values
file_error = open("WEC_error_analysis.txt","w")
file_error.write("WEC analytical wave packet - error analysis:\n")
# historical model error values
roms_benchmark_MAE_error  = 1.3384811426220066e-10 # time=60s - Mean Absolute Error for roms simulation to check if value changes. This result was found by running the model before.
roms_benchmark_RMSE_error = 1.6029348955464868e-10 # time=60s - Root Mean Squared Error for roms simulation to check if value changes. This result was found by running the model before.

# time-steps of interest
time_plt = 60

# Loop through time-steps (or just one time-step)
for tt in range(time_plt,time_plt+1):#len(otime)):
    t     = np.mod(otime[tt],LLm/Cg)

# calculate position of wave packet for plotting
    pos_p = np.mod(Cg*t-init_shift,LLm) # This is in meters
    # Get integer index for approximately middle of wave packet (int rounds down)
    pos_p_indx = int(pos_p * pm[0]) # where pm is number of grid points per meter
    pos_p_lwr_indx = pos_p_indx - int( window_x * pm[0] )
    pos_p_upr_indx = pos_p_indx + int( window_x * pm[0] )
    print('pos_p = ', pos_p, 'pos_p_int = ', pos_p_indx, 'pos_p_lwr_indx = ', pos_p_lwr_indx,'pos_p_upr_indx = ', pos_p_upr_indx)
    
# Calculate wave height formula in steps
#    xstar = x+init_shift-Cg*t
    xstar_u   = x_u  +init_shift-Cg*t
    xstar_rho = x_rho+init_shift-Cg*t # Think this is fine size xstar is a vector of numbers...
#    xstar2    = xstar + LLm
    xstar2_u    = xstar_u   + LLm
    xstar2_rho  = xstar_rho + LLm # Think this is fine size xstar is a vector of numbers...
#    f[:]  = A*(np.exp(-env*xstar**2)+np.exp(-env*xstar2**2))
    f_u[:]    = A*(np.exp(-env*xstar_u**2)  +np.exp(-env*xstar2_u**2))
    f_rho[:]  = A*(np.exp(-env*xstar_rho**2)+np.exp(-env*xstar2_rho**2))
#    a2[:] = f[:]*f[:]
    a2_u[:]   = f_u[:]  *f_u[:]
    a2_rho[:] = f_rho[:]*f_rho[:]  

##### Currently set to zero!
    zeta_hat[:] = 0.0#-a2[:]*k/(2.*np.sinh(2*k*h))
    zeta_hat_u[:] = 0.0

#    ust2d[:] = a2[:]*sigma/(2.*h)#*np.tanh(k*h)) # I divided by h to get ust2d not Tst
    ust2d_a[:] = a2_u[:]*sigma/(2.*h)#*np.tanh(k*h)) # I divided by h to get ust2d not Tst
    ust2d_rho_a[:] = a2_rho[:]*sigma/(2.*h)#*np.tanh(k*h)) # I divided by h to get ust2d not Tst

 
    #zeta_hat[:] = -a2/(4.*h) 
    u_lw_a[:]   = -((Clw**2)/((Clw**2)-(Cg**2)))*(ust2d_a-(Cg*zeta_hat_u/h))  #-(3./4)*a2*np.sqrt(g)/((k**2)*(h**(7./2)))
    u_lw_rho_a[:]   = -((Clw**2)/((Clw**2)-(Cg**2)))*(ust2d_rho_a-(Cg*zeta_hat/h))

    zeta_a[:]   = u_lw_rho_a*(Cg/g) # -(3./4)*a2/((k**2)*(h**3))


############# ERROR CALCS #############################

    file_error.writelines(["\n\nTime-step: ", str(tt), ", time: ", str(t), "secs"])

#### ZETA ERROR ####

    file_error.write("\n*** ZETA error analysis ***")

    # error = data (or analytical) - modelled
#    error_zeta_old = zeta_a - zeta_old[tt,:]
    error_zeta_new = zeta_a - zeta_new[tt,:]
    # range of entries window around wave packet - indices: pos_p_lwr_indx < i < pos_p_upr_indx
#    error_zeta_window_old = error_zeta_old[pos_p_lwr_indx:(pos_p_upr_indx+1)] 
    error_zeta_window_new = error_zeta_new[pos_p_lwr_indx:(pos_p_upr_indx+1)] # don't need -1 as python index from zero # +1 needed else won't include last value (python indexing)
#    x_rho_test = x_rho[pos_p_lwr_indx:(pos_p_upr_indx+1)] # +1 needed else won't include last value (python indexing)
#    print('x_rho_test[0]=',x_rho_test[0],'x_rho_test[-1]=',x_rho_test[-1])
    file_error.writelines(["\nDomain range (m): ", str(x_rho[pos_p_lwr_indx]), " - ", str(x_rho[pos_p_upr_indx]) ])
    file_error.writelines(["\nNumber of rho-points (m): ", str(len(error_zeta_window_new)) ])

    file_error.write("\n--- Benchmark Error Values (MAE) ---")
    file_error.writelines(["\nA) roms_benchmark_MAE_error = ", str(roms_benchmark_MAE_error) ])
    file_error.writelines(["\nB) roms_benchmark_RMSE_error = ", str(roms_benchmark_RMSE_error) ])


###### Mean Absolute Error (MAE) around wave packet!
#    MAErr_old = mean( abs( error_zeta_window_old ) ) # -1 as python indx from 0 not 1
    MAErr_new = mean( abs( error_zeta_window_new ) ) # -1 as python indx from 0 not 1
#    print('MAErr_new = ', MAErr_new, 'MAErr_old = ', MAErr_old) 


###### Root Mean Square Error (RMSE) around wave packet!
#    RMSErr_old = sqrt( mean( error_zeta_window_old**2 ) ) # -1 as python indx from 0 not 1
    RMSErr_new = sqrt( mean( error_zeta_window_new**2 ) ) # -1 as python indx from 0 not 1
#    print('RMSErr_new = ', RMSErr_new, 'RMSErr_old = ', RMSErr_old)


###### Values as a percentage of max zeta analytical
    file_error.write("\n\n--- Mean Absolute Error (MAE) ---")
    #MAE percent of max
#    # old code
#    MAE_old_precent_of_max = ( MAErr_old / max( abs(zeta_a) ) ) * 100 # zeta is negative hence max(abs(zeta))
#    text_err = ['\nMAErr_old = ', str(MAErr_old), '\nMAE_old_precent_of_max(%) = ', str(MAE_old_precent_of_max), '\nmax(abs(zeta_a)) = ', str(max(abs(zeta_a))) ] # Variables to print to terminal and file
#    print(text_err)
#    file_error.writelines(text_err) # write error values to text file
    # new code
    MAE_new_precent_of_max = ( MAErr_new / max( abs(zeta_a) ) ) * 100 # zeta is negative hence max(abs(zeta))
    text_err = ['\nA) MAErr_new = ', str(MAErr_new), '\nMAE_new_precent_of_max(%) = ', str(MAE_new_precent_of_max), '\nmax(abs(zeta_a)) = ', str(max(abs(zeta_a))) ] # Variables to print to terminal and file
    print(text_err)
    file_error.writelines(text_err) # write error values to text file


##### RMSE percent of max
    file_error.write("\n\n--- Root Mean Squared Error (RMSE) ---")
#    # old code
#    RMSE_old_precent_of_max = ( RMSErr_old / max( abs(zeta_a) ) ) * 100 # zeta is negative hence max(abs(zeta))
#    text_err = ['\nRMSErr_old = ', str(RMSErr_old), '\nRMSE_old_precent_of_max(%) = ', str(RMSE_old_precent_of_max), '\nmax(abs(zeta_a)) = ', str(max(abs(zeta_a))) ]
#    print(text_err)
#    file_error.writelines(text_err) # write error values to text file
    # new code
    RMSE_new_precent_of_max = ( RMSErr_new / max( abs(zeta_a) ) ) * 100 # zeta is negative hence max(abs(zeta))
    text_err = ['\nB) RMSErr_new = ', str(RMSErr_new), '\nRMSE_new_precent_of_max(%) = ', str(RMSE_new_precent_of_max), '\nmax(abs(zeta_a)) = ', str(max(abs(zeta_a))) ]
    print(text_err)
    file_error.writelines(text_err) # write error values to text file


##### Max single error
    file_error.write("\n\n--- Max error value ---")
    abs_error_zeta_window_new = abs(error_zeta_window_new)
    max_abs_err = max( abs_error_zeta_window_new ) # use window because of free wave creatign max error elsewhere
    max_abs_err_indx = np.argmax( abs_error_zeta_window_new )  # get the corresponding index of the max error
    text_err = ['\nMax zeta abs. error = ', str(max_abs_err), ' at rho-point: ', str(max_abs_err_indx + pos_p_lwr_indx), ' at coord (m) x = ', str( x_rho[max_abs_err_indx + pos_p_lwr_indx] ) ] # NO need to add +1 for rho node number as rho index starts from 0 as does python indexing from 0.
    print(text_err)
    file_error.writelines(text_err) # write error values to text file
    ratio_max_to_ana      = max_abs_err / abs(zeta_a[max_abs_err_indx + pos_p_lwr_indx])
    ratio_max_to_ana_max  = max_abs_err / max(abs(zeta_a))
    text_err = ['\nMax error % at point = ', str(ratio_max_to_ana*100), '\nMax error % of max = ', str(ratio_max_to_ana_max*100)]
    print(text_err)
    file_error.writelines(text_err) # write error values to text file


##### Automated confirm if error has not changed, and thus roms result still the same.  
    roms_change_RMSE = roms_benchmark_RMSE_error - RMSErr_new
    file_error.write("\n\n--- Check ROMS against repo benchmark RMSE ---")
    text_err = ['\nROMS change = ROMS benchmark RMSE - ROMS compiled RMSE = ', str( roms_change_RMSE )]
    print(text_err)
    file_error.writelines(text_err) # write error values to text file
    if roms_change_RMSE == 0.0:
	text_err = "\n\nSUCCESSFUL: ROMS build WEC results the same as benchmark ROMS repository version"
    	print(text_err)
    	file_error.writelines(text_err) # write error values to text file
    else:
	text_err = "\n\nFAILED -> WARNING: ROMS build WEC results are not the same as benchmark ROMS repository version. A bug may have been introduced! \nIf the error is less than 1e-14 then there may be no issue, it could just be compiled on a different machine or compiler."
    	print(text_err)
    	file_error.writelines(text_err) # write error values to text file

    print('See WEC_error_analysis.txt for summary of error.')

#### END ZETA ERROR ####	


#################### PLOTTING #############################


    ##### LW, Tst, Zeta ####

    plt.subplot(3,1,1)
    plt.cla()
    plt.xlim([pos_p-window_x,pos_p+window_x])#following the wave packet
    #plt.plot(x,ubar[tt,:]-ubarFW[tt,:],label='WEC-FW')
#    plt.plot(x_u,ubar_old[tt,:],label='old')
    #plt.plot(x,ubarFW[tt,:],label='FW')
    plt.plot(x_u,u_lw_a[:],'--',label='analytical')
    plt.plot(x_u,ubar_new[tt,:],':',label='new')
    #plt.plot(x,vel_lw[:]*h+ubarFW[tt,:],label='analytical+FW')
    #plt.plot(x,test,linestyle='--',label='test')
    plt.title('t= '+str(otime[tt])+' s')
    plt.ylabel(r'Long wave velocity, [ms$^{-1}$]', fontsize=label_size)
    plt.xticks(fontsize=axis_size)
    plt.yticks(fontsize=axis_size)
    plt.legend(fontsize=label_size)

    plt.subplot(3,1,2)
    plt.cla()
    plt.xlim([pos_p-window_x,pos_p+window_x])
#    plt.plot(x,h*ust2d[tt,:],label='old')
#    plt.plot(x_u,ust2d_old[tt,:],label='old')
#    plt.plot(x,Tst[:],linestyle='--',label='analytical')
    plt.plot(x_u,ust2d_a[:],linestyle='--',label='analytical')
    plt.plot(x_u,ust2d_new[tt,:],':',label='new') 
    plt.ylabel(r'ust2d, [m$^2$s$^{-1}$]', fontsize=label_size)
    plt.xticks(fontsize=axis_size)
    plt.yticks(fontsize=axis_size)
    plt.legend(fontsize=label_size) 

    plt.subplot(3,1,3)
    plt.cla()
    plt.ylim([-0.00000003,0])
    plt.xlim([pos_p-window_x,pos_p+window_x])
    #plt.plot(x,zeta[tt,:]-zetaFW[tt,:],label=r'WEC-FW') 
#    plt.plot(x_rho,zeta_old[tt,:],label=r'old') 
    #plt.plot(x,zetaFW[tt,:],label=r'FW') 
    plt.plot(x_rho,zeta_a[:],linestyle='--',label='analytical') 
    plt.plot(x_rho,zeta_new[tt,:],':',label='new')
    #plt.plot(x,zeta_a[:]+zetaFW[tt,:],label='analytical+FW') 
    #plt.plot(x,test2,linestyle='--',label='test')
    plt.gca()
    plt.ylabel(r'Elevation, [m]', fontsize=label_size)
    plt.xlabel(r'Distance $x$, [m]', fontsize=label_size)
    plt.xticks(fontsize=axis_size)
    plt.yticks(fontsize=axis_size)
    plt.legend(fontsize=label_size)

    plt.pause(1)
    plt.draw()
    
    # DevinD - save plot
    plt.savefig('wave_packet_offline_'+new_code_name+'.pdf')


    #### Zeta ONLY ####

#    plt.subplot(3,1,1)
#    plt.cla()
#    plt.ylim([-2e-7,2e-7])
#    plt.xlim([0,500])#following the wave packet
#    #plt.plot(x,ubar[tt,:]-ubarFW[tt,:],label='WEC-FW')
##    plt.plot(x,ubar[tt,:],label='ROMS')
#    #plt.plot(x,ubarFW[tt,:],label='FW')
##    plt.plot(x,vel_lw[:],linestyle='--',label='analytical')
#    plt.plot(x,zeta[tt,:],label=r'ROMS')  
#    plt.plot(x,zeta_a[:],linestyle='--',label='analytical') 
#    plt.plot(x,zeta_new[tt,:],':',label='new')
#    #plt.plot(x,vel_lw[:]*h+ubarFW[tt,:],label='analytical+FW')
#    #plt.plot(x,test,linestyle='--',label='test')
#    plt.title('t= '+str(otime[tt])+' s')
##    plt.ylabel(r'Long wave velocity, [ms$^{-1}$]', fontsize=label_size)
#    plt.ylabel(r'Elevation, [m]', fontsize=label_size)
#    plt.xticks(fontsize=axis_size)
#    plt.yticks(fontsize=axis_size)
#    plt.legend(fontsize=label_size)

#    plt.subplot(3,1,2)
#    plt.cla()
#    plt.ylim([-2e-7,2e-7])
#    plt.xlim([500,1000]) 
##    plt.plot(x,Tst[:],linestyle='--',label='analytical')
##    plt.ylabel(r'T$^{St}$, [m$^2$s$^{-1}$]', fontsize=label_size)
#    plt.plot(x,zeta[tt,:],label=r'ROMS')  
#    plt.plot(x,zeta_a[:],linestyle='--',label='analytical')
#    plt.plot(x,zeta_new[tt,:],':',label='new') 
#    plt.ylabel(r'Elevation, [m]', fontsize=label_size)
#    plt.xticks(fontsize=axis_size)
#    plt.yticks(fontsize=axis_size)
#    plt.legend(fontsize=label_size) 

#    plt.subplot(3,1,3)
#    plt.cla()
##    plt.ylim([-0.00000003,0])
#    plt.ylim([-2e-7,2e-7])
#    plt.xlim([1000,1500]) 
#    plt.plot(x,zeta[tt,:],label=r'ROMS')  
#    plt.plot(x,zeta_a[:],linestyle='--',label='analytical')
#    plt.plot(x,zeta_new[tt,:],':',label='new')  
#    plt.gca()
#    plt.ylabel(r'Elevation, [m]', fontsize=label_size)
#    plt.xlabel(r'Distance $x$, [m]', fontsize=label_size)
#    plt.xticks(fontsize=axis_size)
#    plt.yticks(fontsize=axis_size)
#    plt.legend(fontsize=label_size)

#    plt.pause(3)
#    plt.draw()
#    
#    # DevinD - save plot
#    plt.savefig('wave_xi_vs_zeta.pdf')


    #### Difference plots #######


#    plt.subplot(3,1,1)
#    plt.cla()
#    plt.xlim([pos_p-window_x,pos_p+window_x])#following the wave packet
#    #plt.plot(x,ubar[tt,:]-ubarFW[tt,:],label='WEC-FW')
#    plt.plot(x,ubar[tt,:],label='old')
#    #plt.plot(x,ubarFW[tt,:],label='FW')
#    plt.plot(x,vel_lw[:],'--',label='analytical')
#    plt.plot(x,ubar_new[tt,:],':',label='new')
#    #plt.plot(x,vel_lw[:]*h+ubarFW[tt,:],label='analytical+FW')
#    #plt.plot(x,test,linestyle='--',label='test')
#    plt.title('t= '+str(otime[tt])+' s')
#    plt.ylabel(r'Long wave velocity, [ms$^{-1}$]', fontsize=label_size)
#    plt.xticks(fontsize=axis_size)
#    plt.yticks(fontsize=axis_size)
#    plt.legend(fontsize=label_size)

#    plt.subplot(3,1,2)
#    plt.cla()
##    plt.ylim([-0.00000003,0])
#    plt.xlim([pos_p-window_x,pos_p+window_x]) 
##    plt.plot(x,zeta[tt,:]-zeta_a[:],label=r'old-error')  
#    plt.plot(x_rho,zeta[tt,:]-zeta_a[:],label=r'old-error')    
##    plt.plot(x,zeta_new[tt,:]-zeta_a[:],':',label='new-error') 
#    plt.plot(x_rho[0:-2],zeta_new[tt,0:-2]-zeta_a[1:-1],':',label='new-error')
##    plt.plot(x,zeta_new[tt,:]-zeta_a_shift[:],':',label='new-er-shft')
##    plt.plot(x_rho,zeta_new[tt,:]-zeta_a[:],':',label='new-er-shft')
#    plt.ylabel(r'Diff Elev, [m]', fontsize=label_size)
#    plt.xticks(fontsize=axis_size)
#    plt.yticks(fontsize=axis_size)
#    plt.legend(fontsize=label_size) 

#    plt.subplot(3,1,3)
#    plt.cla()
#    plt.ylim([-0.00000003,0])
#    plt.xlim([pos_p-window_x,pos_p+window_x]) 
##    plt.plot(x,zeta[tt,:],label=r'old')  
#    plt.plot(x_rho,zeta[tt,:],label=r'old') 
##    plt.plot(x,zeta_a[:],linestyle='--',label='analytical') 
#    plt.plot(x_rho,zeta_a[:],linestyle='--',label='analytical') 
##    plt.plot(x,zeta_new[tt,:],':',label='new') 
#    plt.plot(x_rho,zeta_new[tt,:],':',label='new') # DevinD changed for x_rho point
#    plt.gca()
#    plt.ylabel(r'Elevation, [m]', fontsize=label_size)
#    plt.xlabel(r'Distance $x$, [m]', fontsize=label_size)
#    plt.xticks(fontsize=axis_size)
#    plt.yticks(fontsize=axis_size)
#    plt.legend(fontsize=label_size)

#    plt.pause(3)
#    plt.draw()
#    
#    # DevinD - save plot
#    plt.savefig('difference_wave_packet_t_'+str(time_plt)+'.pdf')


    ###### xi vs error ##############

#    plt.subplot(3,1,1)
#    plt.cla()
##    plt.ylim([-2e-7,2e-7])
#    plt.xlim([0,500])#following the wave packet 
#    plt.plot(x,zeta[tt,:]-zeta_a[:],label=r'old-error')   
#    plt.plot(x,zeta_new[tt,:]-zeta_a[:],':',label='new-error') 
#    plt.title('t= '+str(otime[tt])+' s') 
#    plt.ylabel(r'Diff Elev, [m]', fontsize=label_size)
#    plt.xticks(fontsize=axis_size)
#    plt.yticks(fontsize=axis_size)
#    plt.legend(fontsize=label_size)

#    plt.subplot(3,1,2)
#    plt.cla()
##    plt.ylim([-2e-7,2e-7])
#    plt.xlim([500,1000]) 
#    plt.plot(x,zeta[tt,:]-zeta_a[:],label=r'old-error')   
#    plt.plot(x,zeta_new[tt,:]-zeta_a[:],':',label='new-error') 
#    plt.ylabel(r'Diff Elev, [m]', fontsize=label_size)
#    plt.xticks(fontsize=axis_size)
#    plt.yticks(fontsize=axis_size)
#    plt.legend(fontsize=label_size) 

#    plt.subplot(3,1,3)
#    plt.cla() 
##    plt.ylim([-2e-7,2e-7])
#    plt.xlim([1000,1500]) 
#    plt.plot(x,zeta[tt,:]-zeta_a[:],label=r'old-error')   
#    plt.plot(x,zeta_new[tt,:]-zeta_a[:],':',label='new-error')  
#    plt.gca()
#    plt.ylabel(r'Diff Elev, [m]', fontsize=label_size)
#    plt.xlabel(r'Distance $x$, [m]', fontsize=label_size)
#    plt.xticks(fontsize=axis_size)
#    plt.yticks(fontsize=axis_size)
#    plt.legend(fontsize=label_size)

#    plt.pause(3)
#    plt.draw()
#    
#    # DevinD - save plot
#    plt.savefig('difference_xi_domain_t_'+str(time_plt)+'.pdf')

file_error.close() 

