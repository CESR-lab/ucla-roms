      module ncvars

#include "cppdefs.opt"

      use param
      implicit none


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

#endif

! Length of netCDF variable "time_step"

      integer, parameter :: iaux=6

! vname    character array for variable names and attributes;

      integer, parameter :: max_frc_files=360
      integer max_frc, ncfrc(max_frc_files), nrrec



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
      character(len=max_name_size) :: ininame, grdname, frcfile(max_frc_files)

      character(len=max_name_size) clm_file  ! defined TCLIMATOLOGY && !defined ANA_TCLIMA) || !defined ANA_SSH

      integer, parameter :: max_bry_files=360
      character(len=max_name_size) bry_file(max_bry_files)

      character*42  vname(3, 180) ! DevinD 180 copied from old code only


      end module ncvars
