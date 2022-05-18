      module ncvars

#include "cppdefs.opt"

      use param
      implicit none

! This is "ncvars":   indices for character array "vname" to keep
!----- -- -----------   names of netCDF variables and their attributes
! indxTime      time
! indxZ         free-surface
! indxUb,Vb     vertically integrated 2D U,V-momentum components
! indxU,V       3D u- and v-velocity components
! indxT,S,...,Zoo tracers (temperature, salinity, biotracers).
! indxO         "omega" vertical velocity (flux)
! indxW         "true" vertical velocity
! indxR         density anomaly
! indxAkv,Akt,Aks vertical viscosity/diffusivity coefficients
! indxHbls      depth of planetary boundary layer in KPP model
! indxHbbl      depth of bottom boundary layer in BKPP model

! indxAi        fraction of cell covered by ice
! indxUi,Vi     U,V-components of sea ice velocity
! indxHi,HS     depth of ice cover and depth of snow cover
! indxTIsrf     temperature of ice surface

! indxBSD,BSS bottom sediment grain Density and Size.
! indxWWA,WWD,WWP  wind induced wave Amplitude, Direction,and Period

      integer, parameter :: indxTime=1, indxZ=2, indxUb=3, indxVb=4
#ifdef SOLVE3D
     &                    , indxU=5, indxV=6, indxO=7, indxW=8
     &                    , indxR=9, indxT=10
# ifdef SALINITY
     &                    , indxS=indxT+1
# endif
     &                    , indxAkv=indxT+NT,   indxAkt=indxAkv+1
# ifdef SALINITY
     &                    , indxAks=indxAkt+1
# endif
# ifdef LMD_KPP
#  ifdef SALINITY
     &                    , indxHbls=indxAks+1
#  else
     &                    , indxHbls=indxAkt+1
#  endif
# endif
# ifdef LMD_BKPP
     &                    , indxHbbl=indxHbls+1
# endif

# ifdef SG_BBL96
#  ifndef ANA_WWAVE
     &                    , indxWWA=???          indxWWD=indxWWA+1
     &                    , indxWWP=indxWWA+2
#  endif
# endif
#endif
#ifdef ICE
      integer, parameter :: indxAi=????,     indxUi=indxAi+1,
     &                    , indxVi=indxAi+2, indxHi=indxAi+3,
     &                      indxHS=indxAi+4, indxTIsrf=indxAi+5
#endif

! Length of netCDF variable "time_step"

      integer, parameter :: iaux=6


! Naming conventions for indices, variable IDs, etc...
!------- ----------- --- -------- -------- ---- ------
! prefix ncid_  means netCDF ID for netCDF file
!        nrec_  record number in netCDF file since initialization
!        nrpf_  maximum number of records per file  (output netCDF
!                                                       files only)
! prefix/ending rst_/_rst refers to restart  netCDF file
!               his_/_his           history
!               avg_/_avg           averages
!               stn_/_stn           stations
!                    _frc           forcing
!                    _clm           climatology
!
! endings refer to:  ___Time  time [in seconds]
!                    ___Tstep time step numbers and record numbers
! all objects with   ___Z     free-surface
! these endings are  ___Ub    vertically integrated 2D U-momentum
! either:            ___Vb    vertically integrated 2D V-momentum
!
!  netCDF IDs, if    ___U     3D U-momentum
!  occur with prefix ___V     3D V-momentum
!  rst/his/avg/stn   ___T(NT) tracers
!                    ___R     density anomaly
! or                 ___O     omega vertical velocity
!                    ___W     true vertical velocity
!  parameter indices
!  if combined with  ___Akv   vertical viscosity coefficient
!  prefix "indx"     ___Akt   vertical T-diffusion coefficient
!  (see above).      ___Aks   vertical S-diffusion coefficient
!                    ___Hbl   depth of KPP surface boundary layer
!                    ___Hbbl  thickness of KPP bottom boundary layer
!
! vname    character array for variable names and attributes;

      integer, parameter :: max_frc_files=360
      integer max_frc, ncfrc(max_frc_files), nrst,  ncrst,   nrecrst,
     &      nrrec, nrpfrst, ncidclm

! NetCFD IDs for model variables

      integer rstTime, rstTstep,  rstZ, rstUb, rstVb
#ifdef SOLVE3D
# ifdef EXACT_RESTART
#  ifdef EXTRAP_BAR_FLUXES
      integer rst_DU_avg2,    rst_DV_avg2,
     &        rst_DU_avg_bak, rst_DV_avg_bak
#  elif defined PRED_COUPLED_MODE
      integer rst_rufrc, rst_rvfrc
#  endif
# endif
      integer rstU, rstV, rstT(NT+1)

# ifdef LMD_KPP
      integer rstHbls
# endif
# ifdef LMD_BKPP
      integer rstHbbl
# endif
#endif


#ifdef SOLVE3D
! DevinD - 180 and 24 are values from old code
# define NWRTHIS 180
!# define NWRTHIS 16+NT-2
#else
# define NWRTHIS 24
!# define NWRTHIS 14
#endif

! Horizontal Grid Type Codes =  0,1,2,3 for RHO-, U-, V-, PSI-points

      integer, parameter :: r_var=0, u_var=1, v_var=2, q_var=3

!            Horizontal array dimensions in netCDF files. In the case
! xi_rho     of MPI code with PARALLEL_FILES activated these dimensions
! xi_u       depend on the corresponding sizes of individual subdomains
! eta_rho    rather than the whole physical grid, therefore they become
! eta_v      live variables placed in common block and assigned values
!            in "mpi_setup" rather tnan be parameters defined here.

#if defined MPI && defined PARALLEL_FILES
      integer xi_rho, xi_u, eta_rho, eta_v
#else
      integer, parameter :: xi_rho=LLm+2, eta_rho=MMm+2,
     &                     xi_u=xi_rho-1, eta_v=eta_rho-1
#endif

      integer, parameter :: max_name_size=256
      character date_str*44, title*80
      character(len=max_name_size) :: ininame, grdname, rstname,
     &                             hisname, frcfile(max_frc_files) ! DevinD hisname not needed due to output_root_name
#ifdef AVERAGES
      character(len=max_name_size) avgname
#endif
#if (defined TCLIMATOLOGY && !defined ANA_TCLIMA) || !defined ANA_SSH
      character(len=max_name_size) clm_file
#endif
#if defined T_FRC_BRY  || defined M2_FRC_BRY || \
    defined M3_FRC_BRY || defined Z_FRC_BRY
      integer, parameter :: max_bry_files=360
      character(len=max_name_size) bry_file(max_bry_files)
#endif

      character*42  vname(3, 180) ! DevinD 180 copied from old code only


      end module ncvars
