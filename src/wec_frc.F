      module wec_frc

      ! Wave Effect on Current (WEC) module
      !
      ! provides all WEC terms that
      ! enter in the momentum equations.
      !
      ! INFO  ![
      !
      ! coding: module format - Devin Dollery & Jeroen Molemaker (2020)
      !         original code - Leonel Romero
      !
      ! code essentially a copy of work already done under MRL_WCI cppflag
      ! of previous roms implementation, based on Delphine Hypolite's
      ! version of old (pre 2020) ROMS used for USWC simulations.
      !
      ! COMMENTS:
      !
      ! 2022/04:
      ! - JM & DPD wec_force should actually be called twice for the predictor step
      !   and the corrector step. Currently only done at the predictor step.
      ! - In set forces we normally interpolate forcing to t+0.5dt.
      !   However, this should likely formally be at t=n (pred) and t=n+1/2 (corr)
      !]

#include "cppdefs.opt"

      use param
      use grid
      use ocean_vars
      use scalars
      use netcdf
      use nc_read_write
      use roms_read_write
      use dimensions
      implicit none

      private                                              ! make all variable private to this module unless public specified
      save                                                 ! make all public variables global (might not be needed depending on compiler)

#include "wec_frc.opt"

! CPP flag WEC covering entire module so that the module is empty if
! WEC is unused, and thus module will not take up unnecessary size in
! the executable.
#ifdef WEC

      ! DECLARE ALL WEC VARIABLES ![

      ! module specification section:
      ! this part of the module (before the 'contains' section)
      ! is used to declare all variables that might be needed globally
      ! or if there are memory benefits to declaring only once, rather
      ! than for each call of a subroutine.

      ! variables as per DH's WAVE_PACKET forces.h for old MRL_WCI cppflag

!--------------------------------------------------------------------
!  WIND INDUCED WAVES: everything is defined at rho-point
!--------------------------------------------------------------------
! wfrq | BBL/MRL | wind-induced wave frequency [rad/s]
! uorb | BBL     | xi-component  of wave-induced bed orbital velocity [m/s]
! vorb | BBL     | eta-component of wave-induced bed orbital velocity [m/s]
! wdrx | MRL     | cosine of wave direction [non dimension]
! wdre | MRL     | sine of   wave direction [non dimension]
! whrm | MRL     | (RMS) wave height (twice the wave amplitude) [m]
! wdsp | MRL     | breaking dissipation rate (\epsilon_b term) [m3/s3]
! wdrg | MRL     | frictional dissipation rate (\epsilon_d term) [m3/s3]
! rdsp | ROLLER  | roller dissipation rate (\epsilon_r term) [m3/s3]
! wbst | MRL/BKPP| frictional dissipation stress (e_d k/sigma) [m2/s2]
!--------------------------------------------------------------------

#if defined BBL || defined WEC
      real,allocatable,dimension(:,:) :: wfrq
#endif

#ifdef BBL
      real,allocatable,dimension(:,:) :: uorb
      real,allocatable,dimension(:,:) :: vorb
#endif   /* BBL */


      real,allocatable,dimension(:,:) :: whrm
      real,allocatable,dimension(:,:) :: wdsp
      real,allocatable,dimension(:,:) :: wdrg
      real,allocatable,dimension(:,:) :: wbst
      real,allocatable,dimension(:,:) :: wdrx
      real,allocatable,dimension(:,:) :: wdre
# ifdef SURFACE_ROLLER
      real,allocatable,dimension(:,:) :: rdsp
# endif

!--------------------------------------------------------------------
!  WAVE AVEREAGED QUANTITIES AND TERMS
!--------------------------------------------------------------------
!  2D  |  brk2dx   |   xi-direciton 2D breaking dissipation (rho)
!  2D  |  brk2de   |  eta-direction 2D breaking dissipation (rho)
!  2D  |  frc2dx   |   xi-direciton 2D frictional dissipation (rho)
!  2D  |  frc2de   |  eta-direction 2D frictional dissipation (rho)
!  2D  |  ust2d    |   xi-direciton Stokes transport (u-point)
!  2D  |  vst2d    |  eta-direciton Stokes transport (v-point)
!  2D  |  ustr2d   |   xi-direciton Stokes transport (r-point)
!  2D  |  vstr2d   |  eta-direciton Stokes transport (r-point)
!  2D  |  ustr0    |   xi-direciton Surface Stokes drift (r-point)
!  2D  |  vstr0    |  eta-direciton Surface Stokes drift (r-point)
!  2D  |  sup      |  quasi-static wave set-up (rho-point)
!  2D  |  calP     |  pressure correction term (rho-point)
!  2D  |  Kapsrf   |  Bernoulli head terrm at the surface (rho-point)
!  2D  |  lm       |  mean wave length (rho-point)
!--------------------------------------------------------------------
!  3D  |  brk3dx   |   xi-direciton 3D breaking dissipation (rho)
!  3D  |  brk3de   |  eta-direction 3D breaking dissipation (rho)
!  3D  |  ust      |   xi-direciton 3D Stokes drift velocity (u-point)
!  3D  |  vst      |  eta-direciton 3D Stokes drift velocity (v-point)
!  3D  |  wst      |       vertical 3D Stokes drift velocity (rho-point)
!  3D  |  kvf      |  vertical vortex force term (K term, 3D, rho-point)
!  3D  |  Akb      |  breaking-wave-induced additional diffusivity (w-point)
!  3D  |  Akw      |  wave-induced additional diffusivity (rho-point)
!  3D  |  E_pre    |  previous time-step value for Akw estimation (rho)
!  3D  |  frc3dx   |   xi-direction 3D frictional dissipation (rho)
!  3D  |  frc3de   |  eta-direction 3D frictional dissipation (rho)
!--------------------------------------------------------------------

      real,allocatable,dimension(:,:),public :: brk2dx
      real,allocatable,dimension(:,:),public :: brk2de
      real,allocatable,dimension(:,:),public :: ust2d
      real,allocatable,dimension(:,:),public :: vst2d
      real,allocatable,dimension(:,:)        :: ustr2d ! At rho point
      real,allocatable,dimension(:,:)        :: vstr2d ! At rho point
      real,allocatable,dimension(:,:)        :: ustr0
      real,allocatable,dimension(:,:)        :: vstr0
      real,allocatable,dimension(:,:)        :: llmm
      real,allocatable,dimension(:,:),public :: frc2dx ! used in prestep, step2d, step3d
      real,allocatable,dimension(:,:),public :: frc2de ! used in prestep, step2d, step3d
      real,allocatable,dimension(:,:),public :: sup    ! used in step2d
# ifdef SOLVE3D
      real,allocatable,dimension(:,:),public :: calP   ! used in prsgrd32AC1.F
      real,allocatable,dimension(:,:),public :: Kapsrf ! used in prsgrd32AC1.F
#  ifndef SURFACE_BREAK
      real         :: brk3dx ! used in momentum eqn in wec_frc.F
      real         :: brk3de ! used in momentum eqn in wec_frc.F
#  endif
#  ifdef BODY_FRICTION
      real         :: frc3dx ! used in momentum eqn in wec_frc.F
      real         :: frc3de ! used in momentum eqn in wec_frc.F
#  endif
      real,allocatable,dimension(:,:,:),public :: ust
      real,allocatable,dimension(:,:,:),public :: vst
      real,allocatable,dimension(:,:,:),public :: wst    ! used in set_avg.F
      real,allocatable,dimension(:,:,:),public :: kvf    ! used in prsgrd32AC1.F
      real,allocatable,dimension(:,:,:),public :: Akb    ! used in lmd_kpp.F
      real,allocatable,dimension(:,:,:),public :: Akw    ! used in lmd_kpp.F
      real,allocatable,dimension(:,:,:)        :: E_pre
# endif  /* SOLVE3D */

      ! REALISTIC WEC VARIABLES (FROM READ IN NETCDF FILE)
      ! =======================

#if defined BBL || defined WEC
# if !defined ANA_WEC_FRC && !defined WKB_WWAVE
!--------------------------------------------------------------------
!wwv_time|                    | time of wind-induced waves
!--------------------------------------------------------------------
!  wwag  |                    | wave amplitude [m]
!  wwdg  |  (READ IN FROM     | wave direction [radians]
!  wwpg  |   INPUT FILE)      | wave period [s]
!  wwub  |                    | orbital velocity magnitude [m/s]
!  wwuob |  point data        | xi-orbital velocity [m/s]
!  wwvob |  for wind induced  | eta-orbital velocity [m/s]
!  wwdrx |  waves             | cosine wave direction [ND]
!  wwdre |                    | sine wave direction [ND]
!  wweb  |                    | breaking dissipation [m3/s3]
!  wwed  |                    | frictional dissipation [m3/s3]
!  wwUSg |                    | Depth averaged Stokes drift [m/s]
!  wwVSg |                    | Depth averaged Stokes drift [m/s]
!  wwus0g|                    | Surface Stokes drift [m/s]
!  wwvs0g|                    | Surface Stokes drift [m/s]
!  wwsup |                    | Set down [m]
!  wwlm  |                    | Mean wavelength [m] - for diffusivity approx.
!--------------------------------------------------------------------
!  wveb  |  interpolated      | breaking dissipation [m3/s3]
!  wved  |  for present time  | frictional dissipation [m3/s3]
!  wvqb  |  step              | fraction of breaking waves [ND]
!--------------------------------------------------------------------

      real,allocatable,dimension(:,:,:)        :: wwag
      real,allocatable,dimension(:,:,:)        :: wwdg
      real,allocatable,dimension(:,:,:)        :: wwpg
#  ifdef BBL
      real,allocatable,dimension(:,:,:)        :: wwub
      real,allocatable,dimension(:,:,:)        :: wwuob
      real,allocatable,dimension(:,:,:)        :: wwvob
#  endif
      real,allocatable,dimension(:,:,:)        :: wwdrx
      real,allocatable,dimension(:,:,:)        :: wwdre
      real,allocatable,dimension(:,:,:)        :: wweb
      real,allocatable,dimension(:,:),public   :: wveb
      real,allocatable,dimension(:,:),public   :: wved
      real,allocatable,dimension(:,:,:)        :: wwed
      real,allocatable,dimension(:,:,:)        :: wwUSg
      real,allocatable,dimension(:,:,:)        :: wwVSg
      real,allocatable,dimension(:,:,:)        :: wwus0g
      real,allocatable,dimension(:,:,:)        :: wwvs0g
#   if defined WAVE_OFFLINE && defined SURFACE_ROLLER
      real,allocatable,dimension(:,:),public   :: wvqb
      real,allocatable,dimension(:,:,:)        :: wwqb
#   endif
#   if defined WAVE_OFFLINE && defined SUP_OFF
      real,allocatable,dimension(:,:,:)        :: wwsup
#  endif
#   if defined WAVE_OFFLINE && defined WAVE_DIFF_FROM_LM
      real,allocatable,dimension(:,:,:),public :: wwlm
#   endif

      ! Grid rotation angles
      real,allocatable,dimension(:,:)          :: angler_cos
      real,allocatable,dimension(:,:)          :: angler_sin

!# elif defined ANA_WEC_FRC && defined WEC
!      real, public :: wveb(GLOBAL_2D_ARRAY)
!      real, public :: wved(GLOBAL_2D_ARRAY)

# endif /* !ANA_WEC_FRC && !WKB_WWAVE */
#endif /* BBL || WEC */


      !#### end variables added as per forces.h #####

      ! AVERAGED VARIABLES
      ! ==================

      ! from averages.h in MRL_WCI old roms code
#ifdef AVERAGES
      real,allocatable,dimension(:,:)   :: sup_avg
      real,allocatable,dimension(:,:)   :: ust2d_avg
      real,allocatable,dimension(:,:)   :: vst2d_avg
#  ifdef SOLVE3D
      real,allocatable,dimension(:,:,:) :: ust_avg
      real,allocatable,dimension(:,:,:) :: vst_avg
      real,allocatable,dimension(:,:,:) :: wst_avg
      real,allocatable,dimension(:,:,:) :: akb_avg
      real,allocatable,dimension(:,:,:) :: akw_avg
      real,allocatable,dimension(:,:,:) :: kvf_avg
      real,allocatable,dimension(:,:)   :: calp_avg
      real,allocatable,dimension(:,:)   :: kapsrf_avg
#  endif  /* SOLVE3D */
#endif /* AVERAGES */


      ! ========================
      ! NETCDF RELATED VARIABLES
      ! ========================

      ! Names of variables in WEC's netcdf forcing file or output results
      ! only short name is needed for read in variables
      ! Result variables for output file need
      ! short name, long name, and units.

      ! Wave parameters
      character(len=5) :: wec_awv_name    = 'Awave'    !! Amplitude
      character(len=5) :: wec_dwv_name    = 'Dwave'    !! Direction
      character(len=5) :: wec_pwv_name    = 'Pwave'    !! Period

      ! Vertical average of u_stokes
      character(len=5) :: ust2d_name      = 'ust2d'
      character(len=38):: ust2d_long_name =
     &                   'u depth-averaged Stokes drift velocity'
      character(len=14):: ust2d_units     = 'meter second-1'

      ! Vertical average of v_stokes
      character(len=5) :: vst2d_name      = 'vst2d'
      character(len=38):: vst2d_long_name =
     &                   'v depth-averaged Stokes drift velocity'
      character(len=14):: vst2d_units     = 'meter second-1'

# ifdef SOLVE3D

      ! u_stokes 3d
      character(len=3) :: ust_name        = 'ust'
      character(len=26):: ust_long_name   = 'u 3D Stokes drift velocity'
      character(len=14):: ust_units       = 'meter second-1'

      ! v_stokes 3d
      character(len=3) :: vst_name        = 'vst'
      character(len=26):: vst_long_name   = 'v 3D Stokes drift velocity'
      character(len=14):: vst_units       = 'meter second-1'

      ! Read in data surface u & v stokes drfit
      character(len=4) :: ust0_name  = 'ust0'
      character(len=4) :: vst0_name  = 'vst0'

      ! w_stokes 3d
      character(len=3) :: wst_name        = 'wst'
      character(len=35):: wst_long_name   =
     &                   'w vertical 3D Stokes drift velocity'
      character(len=14):: wst_units       = 'meter second-1'

      ! Akb
      character(len=3) :: Akb_name        = 'Akb'
      character(len=35):: Akb_long_name   =
     &                   'eddy viscosity due to wave breaking'
      character(len=15):: Akb_units       = 'meter2 second-1'

      ! Akw
      character(len=3) :: Akw_name        = 'Akw'
      character(len=37):: Akw_long_name   =
     &                   'eddy diffusivity due to primary waves'
      character(len=15):: Akw_units       = 'meter2 second-1'

      ! kvf
      character(len=3) :: kvf_name        = 'kvf'
      character(len=34):: kvf_long_name   =
     &                   'vertical vortex force (u^st du/dz)'
      character(len=14):: kvf_units       = 'meter second-2'

      ! calP
      character(len=4) :: calP_name        = 'calP'
      character(len=37):: calP_long_name   =
     &                   'surface pressure correction in prsgrd'
      character(len=15):: calP_units       = 'meter2 second-2'

      ! Kapsrf
      character(len=6) :: Kapsrf_name        = 'Kapsrf'
      character(len=37):: Kapsrf_long_name   =
     &                   'surface Bernoulli head in prsgrd term'
      character(len=15):: Kapsrf_units       = 'meter2 second-2'

# endif

      character(len=2) :: wec_ebwv_name = 'eb'      !! wweb
      character(len=2) :: wec_edwv_name = 'ed'      !! wweb
      character(len=3) :: wec_lwv_name = 'lmw'      !! Mean Wave length

      ! sup
      character(len=40):: sup_long_name =
     &                   'quasi-static sea-level response (set-up)'

      !! Wave forcing time (Also its dimension name in netcdf file!!)
      character(len=8) :: wec_tim_name = 'wwv_time'


      ! Add netcdf WEC terms previously contained in ncvars
      ! ===================================================

      ! Index numbers are arbitrary, they are just place holders for
      ! the true or false bool for outputing WEC variables stored in
      ! wrt_wec_vars
      integer, parameter :: indxSUP=1
     &                    , indxUST2D=indxSUP+1, indxVST2D=indxSUP+2
# ifdef SOLVE3D
     &                    , indxUST=indxSUP+3, indxVST=indxSUP+4
     &                    , indxAkb=indxSUP+6
     &                    , indxAkw=indxSUP+7, indxKVF=indxSUP+8
     &                    , indxCALP=indxSUP+9, indxKAPS=indxSUP+10
      ! indxWST used in set_avg.F
      integer, parameter, public :: indxWST=indxSUP+5
# endif


      ! ===============
      ! MISC. VARIABLES
      ! ===============

      ! coefficient variables for external calculations in prsgrd32AC1.F
      real, public :: cff_ext, cff_ext2

      ! timestep variables for read in WEC data
      integer              :: it1 = 1  ! Placeholder for read in time
      integer              :: it2 = 2  ! Placeholder for read in time
      real, dimension(2)   :: wec_tim = [-99,-99] ! used to store both time steps (n <-> n+1)
      integer    :: ifile = 0 ! read in force file number
      ! current record number of data point in read in data - irec
      integer    :: irec  = 0

      ! netcdf outputting:
      integer            :: output_time = 0         ! record number of output. 0 indicates we need new file
      integer            :: record = nrpf           ! to trigger the first file creation
      integer :: ncid=-1, prev_fill_mode
      integer :: navg_wec = 0                       ! number of samples in average

      public set_wec_frc_tile     ! called from set_forces.F
      public create_wec_vars      ! called from def_his.F
      public init_arrays_wec_tile ! called from init_arrays.F
      public compute_wec_rhs_uv_terms ! called from step3d_uv1.F & pre_step3d4S.F
      public wstokes_wec_tile     ! called from set_avg.F
      public wec_forces           ! called from main.F

      ! END OF MODULE PRE-AMBLE  !]

      contains

!-----------------------------------------------------------------------
#ifndef ANA_WEC_FRC
      subroutine read_wec_frc(fname,irec,idx,wec_in_time)  ![
      ! read wec variables from input file
      ! this version is intended to exclusively read WEC variables that are
      ! part of 'Leonel's method'
      use netcdf
      implicit none

      ! input/output
      character(len=*),intent(in) :: fname                           ! this isn't even used! Should be ifile
      integer         ,intent(in) :: irec,idx
      real,dimension(1)           :: wec_in_time                     ! set wec_time for read in variables

      ! local
      integer :: ierr,ncid
      integer :: varid
      integer :: coarse = 0                                          ! interpolate from coarse data. 0=don't interpolate
      integer,dimension(4) :: start

      ierr=nf90_open(frcfile(ifile), nf90_nowrite, ncid)

      call ncread(ncid,wec_tim_name,wec_in_time,(/irec/))            ! needs an array for last term hence syntax
      wec_in_time = wec_in_time * 86400                              ! covnert read in time (currently in days) to seconds

      start=1
      start(3)=irec
      call ncread(ncid,wec_awv_name, wwag(i0:i1,j0:j1,idx),start)    ! wave amplitude
      call exchange_tile(1,nx,1,ny,  wwag(:,:,idx),1)
      call ncread(ncid,wec_dwv_name, wwdg(i0:i1,j0:j1,idx),start)    ! wave direction
      call exchange_tile(1,nx,1,ny,  wwdg(:,:,idx),1)
      call ncread(ncid,wec_pwv_name, wwpg(i0:i1,j0:j1,idx),start)    ! wave period
      call exchange_tile(1,nx,1,ny,  wwpg(:,:,idx),1)
      call ncread(ncid,ust2d_name,  wwUSg(i0:i1,j0:j1,idx),start)    ! ust2d
      call exchange_tile(1,nx,1,ny, wwUSg(:,:,idx),1)
      call ncread(ncid,vst2d_name,  wwVSg(i0:i1,j0:j1,idx),start)    ! vst2d
      call exchange_tile(1,nx,1,ny, wwVSg(:,:,idx),1)
      call ncread(ncid,ust0_name,  wwUS0g(i0:i1,j0:j1,idx),start)    ! ustr0
      call exchange_tile(1,nx,1,ny,wwUS0g(:,:,idx),1)
      call ncread(ncid,vst0_name,  wwVS0g(i0:i1,j0:j1,idx),start)    ! vstr0
      call exchange_tile(1,nx,1,ny,wwVS0g(:,:,idx),1)
      call ncread(ncid,wec_ebwv_name,wweb(i0:i1,j0:j1,idx),start)    ! eb
      call exchange_tile(1,nx,1,ny,  wweb(:,:,idx),1)
      call ncread(ncid,wec_edwv_name,wwed(i0:i1,j0:j1,idx),start)    ! ed
      call exchange_tile(1,nx,1,ny,  wwed(:,:,idx),1)
      call ncread(ncid,wec_lwv_name, wwlm(i0:i1,j0:j1,idx),start)    ! lmv
      call exchange_tile(1,nx,1,ny,  wwlm(:,:,idx),1)
      call ncread(ncid,'sup',wwsup(i0:i1,j0:j1,idx),start)
      call exchange_tile(1,nx,1,ny, wwsup(:,:,idx),1)

      end subroutine read_wec_frc  !]
#endif /* !ANA_WEC_FRC */

! ----------------------------------------------------------------------
      subroutine create_wec_vars(ncid)  ![
      ! Add WEC variables to an open netcdf file
      implicit none

      ! input
      integer,intent(in) :: ncid
      ! local
      integer           :: ierr, varid

      if (wrt_SUP) then                                ! if = outputing variable yes or no
        varid = nccreate(ncid,'sup',(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',sup_long_name)
        ierr = nf90_put_att(ncid,varid,'units','m')
      endif
      if (wrt_UST2D) then                              ! ust2d: 2D, depth-averaged Stokes drift velocities.
        varid = nccreate(ncid,'ust2d',(/dn_xu,dn_yr,dn_tm/),(/xi_u,eta_rho,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',ust2d_long_name)
        ierr = nf90_put_att(ncid,varid,'units',ust2d_units)
      endif
      if (wrt_VST2D) then                              ! vst2d: 2D, depth-averaged Stokes drift velocities.
        varid = nccreate(ncid,vst2d_name,(/dn_xr,dn_yv,dn_tm/),(/xi_rho,eta_v,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',vst2d_long_name)
        ierr = nf90_put_att(ncid,varid,'units',vst2d_units)
      endif
# ifdef SOLVE3D
      if (wrt_UST) then                                ! ust: 3D Stokes drift velocities.
        varid = nccreate(ncid,ust_name,(/dn_xu,dn_yr,dn_zr,dn_tm/),(/xi_u,eta_rho,N,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',ust_long_name)
        ierr = nf90_put_att(ncid,varid,'units',ust_units)
      endif
      if (wrt_VST) then                                ! vst: 3D Stokes drift velocities.
        varid = nccreate(ncid,vst_name,(/dn_xr,dn_yv,dn_zr,dn_tm/),(/xi_rho,eta_v,N,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',vst_long_name)
        ierr = nf90_put_att(ncid,varid,'units',vst_units)
      endif
      if (wrt_WST) then                                ! wst: vertical Stokes drift velocity (m/s) at rho-point
        varid =
     &     nccreate(ncid,wst_name,(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',wst_long_name)
        ierr = nf90_put_att(ncid,varid,'units',wst_units)
      endif
      if (wrt_Akb) then                                ! Akb: Vertical eddy viscosity coef
        varid = nccreate(ncid,Akb_name,(/dn_xr,dn_yr,dn_zw,dn_tm/),(/xi_rho,eta_rho,N+1,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',Akb_long_name)
        ierr = nf90_put_att(ncid,varid,'units',Akb_units)
      endif
      if (wrt_Akw) then                                ! Akw: Vertical eddy diffusivity
        varid = nccreate(ncid,Akw_name,(/dn_xr,dn_yr,dn_zw,dn_tm/),(/xi_rho,eta_rho,N+1,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',Akw_long_name)
        ierr = nf90_put_att(ncid,varid,'units',Akw_units)
      endif
      if (wrt_KVF) then                                ! kvf: vertical vortex force term (u^St du/dz)
        varid =
     &    nccreate(ncid,kvf_name,(/dn_xr,dn_yr,dn_zr,dn_tm/),(/xi_rho,eta_rho,N,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',kvf_long_name)
        ierr = nf90_put_att(ncid,varid,'units',kvf_units)
      endif
      if (wrt_CALP) then                               ! calP: surface pressure correction
        varid = nccreate(ncid,calP_name,(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',calP_long_name)
        ierr = nf90_put_att(ncid,varid,'units',calP_units)
      endif
      if (wrt_KAPS) then                               ! Kapsrf: surface Bernoulli head
        varid = nccreate(ncid,Kapsrf_name,(/dn_xr,dn_yr,dn_tm/),(/xi_rho,eta_rho,0/))
        ierr = nf90_put_att(ncid,varid,'long_name',Kapsrf_long_name)
        ierr = nf90_put_att(ncid,varid,'units',Kapsrf_units)
      endif
# endif  /* SOLVE3D */

      end subroutine create_wec_vars  !]
!----------------------------------------------------------------------
      subroutine create_wec_file(fname)  ![
      implicit none

      !input/output
      character(len=99),intent(out) :: fname

      ! local
      integer :: ierr,varid
      character(len=10),dimension(4) :: dimnames           ! dimension names
      integer,          dimension(4) :: dimsizes

      if (wrt_file_avg) then
        call create_file('_wec_avg',fname)
      else
        call create_file('_wec',fname)
      endif

      ierr=nf90_open(fname,nf90_write,ncid)

      call create_wec_vars(ncid)

      if (wrt_file_avg) then
        ierr=nf90_put_att(ncid,nf90_global,'type','wec forcing average')
      else
        ierr=nf90_put_att(ncid,nf90_global,'type','wec forcing history')
      endif

      ierr = nf90_close(ncid)

      end subroutine create_wec_file !]
! ----------------------------------------------------------------------
      subroutine wrt_wec_data  ![
      ! write wec variables to output netcdf file
      ! don't include t=0 in averaging. This create 0.5dt error in averaging,
      ! but this 0.5dt error has always been in ROMS.
      ! for 2 steps. True avg would be 0.5*t0 + t1 + 0.5*t2, but we've never done that.
      implicit none

      ! local
      integer,dimension(4)   :: start
      logical,save           :: first_step=.true.
      character(len=99),save :: fname
      integer :: tile, ierr

      if (.not. first_step) then
        output_time = output_time + dt
      endif
      first_step=.false.

      if (wrt_file_avg) call calc_wec_avg

      if (output_time>=output_period) then ! time for an output
        navg_wec=0
        output_time = 0

        if (record==nrpf) then
          call create_wec_file(fname)
          record = 0
        endif
        record = record + 1

        ierr=nf90_open(fname,nf90_write,ncid)
        ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)

        call ncwrite(ncid,'ocean_time',(/time/),(/record/))

        if (wrt_sup) then
          if (wrt_file_avg) then
            call ncwrite(ncid,'sup',sup_avg(i0:i1,j0:j1),(/1,1,record/))
          else
            call ncwrite(ncid,'sup',sup(i0:i1,j0:j1),(/1,1,record/))
          endif
        endif
        if (wrt_ust2d) then
          if (wrt_file_avg) then
            call ncwrite(ncid,'ust2d',ust2d_avg(1:i1,j0:j1),(/1,1,record/))
          else
            call ncwrite(ncid,'ust2d',ust2d(1:i1,j0:j1),(/1,1,record/))
          endif
        endif
        if (wrt_ust2d) then
          if (wrt_file_avg) then
            call ncwrite(ncid,'vst2d',vst2d_avg(i0:i1,1:j1),(/1,1,record/))
          else
            call ncwrite(ncid,'vst2d',vst2d(i0:i1,1:j1),(/1,1,record/))
          endif
        endif
# ifdef SOLVE3D
        if (wrt_ust) then
          if (wrt_file_avg) then
            call ncwrite(ncid,'ust',ust_avg(1:i1,j0:j1,:),(/1,1,1,record/))
          else
            call ncwrite(ncid,'ust',ust(1:i1,j0:j1,:),(/1,1,1,record/))
          endif
        endif
        if (wrt_vst) then
          if (wrt_file_avg) then
            call ncwrite(ncid,'vst',vst_avg(i0:i1,1:j1,:),(/1,1,1,record/))
          else
            call ncwrite(ncid,'vst',vst(i0:i1,1:j1,:),(/1,1,1,record/))
          endif
        endif
        if (wrt_wst) then
          print *,'wrt_wst not functional'
          stop
          do tile=0,NSUB_X*NSUB_E-1
            call wstokes_wec(tile)                                   ! calculate wst
          enddo
!         call ncwrite(ncid,wst_name,wst_avg(i0:i1,j0:j1,:),start)
        endif
        if (wrt_Akb) then
          print *,'wrt_akb not ready'
          stop
!         call ncwrite(ncid,Akb_name,Akb_avg(i0:i1,j0:j1,:),start)
        endif
        if (wrt_Akw) then
          print *,'wrt_akw not ready'
          stop
!         call ncwrite(ncid,Akw_name,Akw_avg(i0:i1,j0:j1,:),start)
        endif
        if (wrt_KVF) then
          print *,'wrt_kvf not ready'
          stop
!         call ncwrite(ncid,kvf_name,kvf_avg(i0:i1,j0:j1,:),start)
        endif
        start(3)=record
        if (wrt_CALP) then
          print *,'wrt_calp not ready'
          stop
!         call ncwrite(ncid,calP_name,    calP_avg(i0:i1,j0:j1),start)
        endif
        if (wrt_KAPS) then
          print *,'wrt_kaps not ready'
          stop
!         call ncwrite(ncid,Kapsrf_name,Kapsrf_avg(i0:i1,j0:j1),start)
        endif
# endif  /* SOLVE3D */


        ierr=nf90_close(ncid)
        if (mynode == 0) then
            write(*,'(7x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')  ! confirm work completed
     &        'wec_frc :: wrote data, tdays =', tdays,
     &        'step =', iic-1, 'rec =', record
        endif
      endif

      end subroutine wrt_wec_data  !]
! ----------------------------------------------------------------------
#ifdef ANA_WEC_FRC
      subroutine ana_init_wec_frc_tile (istr,iend,jstr,jend)  ![
      use tracers, only: t
      implicit none

      ! inputs
      integer, intent(in) :: istr, iend, jstr, jend
      ! local
      integer :: i, j, k
      real    :: kkw,Clw,Cg

      ! compute_auxiliary_bounds.h is needed for istrR, jstrR, etc
      ! This had to appear after other variables declared else wont compile
#include "compute_auxiliary_bounds.h"

      kkw = 2.0*pi                                         ! DPD changed kw to kkw so as not to conflict part 2 also has kw
      Clw = sqrt(g*h(2,2))
      Cg  = 0.5*sqrt(g/kkw)

      do j=jstr-1,jend+1
        do i=istr-1,iend+1
          do k=1,N
            t(i,j,k,1,itemp)=18.0
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
          enddo

          if (i>istr-1) then

            ubar(i,j,1) = -((Clw**2)/((Clw**2)-(Cg**2)))*
     &                        ust2d(i,j)
            ubar(i,j,2) = ubar(i,j,1)                      ! JM not actually sure if that is needed

            ! zeta is rho-point, so averaging u-to-rho
            zeta(i,j,1) = -((Clw**2)/((Clw**2)-(Cg**2)))*
     &                         0.5*( ust2d(i,j)+ust2d(i+1,j) )  *Cg/g
            do k=1,N
              u(i,j,k,1) = ubar(i,j,1)                     ! actually need
              v(i,j,k,1) = vbar(i,j,1)                     ! to init u even if baro
            enddo
          end if
        enddo
      enddo

      end subroutine ana_init_wec_frc_tile  !]
#endif /* ANA_WEC_FRC */

! ----------------------------------------------------------------------
#ifdef ANA_WEC_FRC
      subroutine ana_wec_frc_tile (istr,iend,jstr,jend)  ![
      ! analytical WEC forcing
      implicit none

      ! inputs
      integer, intent(in) :: istr, iend, jstr, jend
      ! local
      integer :: i, j
      real    :: cff1,cff2
      real    :: kkw,Clw,Cg
      real    :: dx,x,timed
      real    :: amp,Lx


         ! get universal scalars, g and pi
! # include "scalars" - included in top of module
         ! get iSW_corn variable
! # include "param.h" - included in top of module
         ! get pm coordinate transformer from grid
! # include "grid" - included in top of module
         ! note haven't included DH's forces.h as all wave parameters declared in WEC

         ! For PART 2 below: the include ocean.h files for mrl_wci part
         ! included in top of module

            !*** PART 1 - from DH's WAVE_PACKET analytical.F ***!

          ! The following code is taken from analytical.F L995 from DH's
          ! WAVE_PACKET example code

! ######## NOTE: stokes variables calculated in wec_forces and take in
!                whrm as rho-point variable because that is how non-analytical
!                code works, so need to be consistent.


! Note: grid node i=1 is the western most node at x=0m <- account for this.
!       Hence also, u-point to rho conversion as index is at u-point
!       Hence: wrhm_rho(i,j) = 0.5*(tmp_whrm_u(i,j)+tmp_whrm_u(i+1,j))

      kkw=2.0*pi  ! Changed kw to kkw so as not to conflict part 2
      env=0.0001
      Cg=0.5*sqrt(g/kkw)
      amp = 0.001 ! wave amplitude
      Lx = gnx*dx ! Domain size in x direction

      dx = dm_u(1,1)  !! the grid has constant spacing dx

      do j=jstr-1,jend+1
        do i=istr-1,iend+1

              ! DevinD - whrm from Eq(29) of DH's 'A propagating wave packet.pdf

          timed = modulo(time, Lx / Cg )  ! t' from eq(29)

              ! ### whrm at rho-point=i ###

              ! Calculate x coordinate at rho-point=i.
              ! Since u-point i=1 should equal 0m for western tile.
              ! Since rho-point i=1 should equal 0.5m for western tile
              ! if grid spacing dx=1m. Hence -0.5*dx for rho-point coord.
          x = ( dble(i)+dble(iSW_corn)-0.5 ) * dx

              ! Position relative to wave centre.
              ! Cg*t' = Cg*cff1 gives metres of distance wave travelled,
              ! other variables in grid points not meters.
              ! -1 is needed because node number i=1 is at coord 0m

          ! gnx*dx = domain size
          cff2 = x - Cg*timed - Lx/2 ! initial shift mid domain
          whrm(i,j)= amp*2.0*( exp(-env * cff2**2)
     &                      + exp( -env * (cff2+ Lx)**2 ) )


          wfrq(i,j)= sqrt(g*2.0*pi)
          wdrx(i,j)= 1.0                             ! cosine wave direction (xi)
          wdre(i,j)= 0.0                             ! sine wave direction (eta)
          wveb(i,j)= 0.0
          wved(i,j)= 0.0
        enddo
      enddo

      ! NOW CALCULATE WEC FORCES USING WHRM, etc
      ! ========================================
        !  for now, in main.
!      call wec_wave2frc_method1(istr,iend,jstr,jend)

!      call wec_wv2frc_method2(istr,iend,jstr,jend)

      ! NOTE: step 0 need to calculate ubar, zeta, ust, vst
      ! ===================================================

        ! Once ust2d & vst2d are calculated in wec_wave2frc_method?
        ! for step 0 will need to call ana_init_wec_frc_tile from that
        ! point.

!      call ana_init_wec_frc_tile(istr,iend,jstr,jend)

      end subroutine ana_wec_frc_tile  !]
#endif /* ANA_WEC_FRC */

! ----------------------------------------------------------------------
      subroutine init_arrays_wec_tile (istr,iend,jstr,jend)  ![
      ! this is called from init_arrays_tile of init_arrays.F
      ! it is needed to "first touch" arrays to optimize
      ! NUMA shared memory allocation on linux.
      implicit none

      ! inputs
      integer,intent(in) :: istr, iend, jstr, jend
      ! local
      integer        :: i, j, k
      real,parameter :: init=0.

#include "compute_auxiliary_bounds.h"

#if defined BBL || defined WEC
      allocate( wfrq(GLOBAL_2D_ARRAY) )
#endif

#ifdef BBL
      allocate( uorb(GLOBAL_2D_ARRAY) )
      allocate( vorb(GLOBAL_2D_ARRAY) )
#endif   /* BBL */

      allocate( whrm(GLOBAL_2D_ARRAY) )
      allocate( wdsp(GLOBAL_2D_ARRAY) )
      allocate( wdrg(GLOBAL_2D_ARRAY) )
      allocate( wbst(GLOBAL_2D_ARRAY) )
      allocate( wdrx(GLOBAL_2D_ARRAY) )
      allocate( wdre(GLOBAL_2D_ARRAY) )
# ifdef SURFACE_ROLLER
      allocate( rdsp(GLOBAL_2D_ARRAY) )
# endif

      allocate( brk2dx(GLOBAL_2D_ARRAY) )
      allocate( brk2de(GLOBAL_2D_ARRAY) )
      allocate( ust2d(GLOBAL_2D_ARRAY) )
      allocate( vst2d(GLOBAL_2D_ARRAY) )
      allocate( ustr2d(GLOBAL_2D_ARRAY) ) ! At rho point
      allocate( vstr2d(GLOBAL_2D_ARRAY) ) ! At rho point
      allocate( ustr0(GLOBAL_2D_ARRAY) )
      allocate( vstr0(GLOBAL_2D_ARRAY) )
      allocate( llmm(GLOBAL_2D_ARRAY) )
      allocate( frc2dx(GLOBAL_2D_ARRAY) ) ! used in prestep, step2d, step3d
      allocate( frc2de(GLOBAL_2D_ARRAY) ) ! used in prestep, step2d, step3d
      allocate( sup(GLOBAL_2D_ARRAY) ) ! used in step2d
# ifdef SOLVE3D
      allocate( calP(GLOBAL_2D_ARRAY) ) ! used in prsgrd32AC1.F
      allocate( Kapsrf(GLOBAL_2D_ARRAY) ) ! used in prsgrd32AC1.F
#  ifndef SURFACE_BREAK
      allocate( brk3dx(GLOBAL_2D_ARRAY,N) ) ! used in momentum eqn in wec_frc.F
      allocate( brk3de(GLOBAL_2D_ARRAY,N) ) ! used in momentum eqn in wec_frc.F
#  endif
#  ifdef BODY_FRICTION
      allocate( frc3dx(GLOBAL_2D_ARRAY,N) ) ! used in momentum eqn in wec_frc.F
      allocate( frc3de(GLOBAL_2D_ARRAY,N) ) ! used in momentum eqn in wec_frc.F
#  endif
      allocate( ust(GLOBAL_2D_ARRAY,N) )
      allocate( vst(GLOBAL_2D_ARRAY,N) )
      allocate( wst(GLOBAL_2D_ARRAY,N) ) ! used in set_avg.F
      allocate( kvf(GLOBAL_2D_ARRAY,N) ) ! used in prsgrd32AC1.F
      allocate( Akb(GLOBAL_2D_ARRAY,0:N) ) ! used in lmd_kpp.F
      allocate( Akw(GLOBAL_2D_ARRAY,0:N) ) ! used in lmd_kpp.F
      allocate( E_pre(GLOBAL_2D_ARRAY,0:N) )
# endif  /* SOLVE3D */

      ! REALISTIC WEC VARIABLES (FROM READ IN NETCDF FILE)

#if defined BBL || defined WEC
# if !defined ANA_WEC_FRC && !defined WKB_WWAVE

      allocate( wwag(GLOBAL_2D_ARRAY,2) )
      allocate( wwdg(GLOBAL_2D_ARRAY,2) )
      allocate( wwpg(GLOBAL_2D_ARRAY,2) )
#  ifdef BBL
      allocate( wwub(GLOBAL_2D_ARRAY,2) )
      allocate( wwuob(GLOBAL_2D_ARRAY,2) )
      allocate( wwvob(GLOBAL_2D_ARRAY,2) )
#  endif
      allocate( wwdrx(GLOBAL_2D_ARRAY,2) )
      allocate( wwdre(GLOBAL_2D_ARRAY,2) )
      allocate( wweb(GLOBAL_2D_ARRAY,2) )
      allocate( wveb(GLOBAL_2D_ARRAY) )
      allocate( wved(GLOBAL_2D_ARRAY) )
      allocate( wwed(GLOBAL_2D_ARRAY,2) )
      allocate( wwUSg(GLOBAL_2D_ARRAY,2) )
      allocate( wwVSg(GLOBAL_2D_ARRAY,2) )
      allocate( wwus0g(GLOBAL_2D_ARRAY,2) )
      allocate( wwvs0g(GLOBAL_2D_ARRAY,2) )
#   if defined WAVE_OFFLINE && defined SURFACE_ROLLER
      allocate( wvqb(GLOBAL_2D_ARRAY) )
      allocate( wwqb(GLOBAL_2D_ARRAY,2) )
#   endif
#   if defined WAVE_OFFLINE && defined SUP_OFF
      allocate( wwsup(GLOBAL_2D_ARRAY,2) )
#  endif
#   if defined WAVE_OFFLINE && defined WAVE_DIFF_FROM_LM
      allocate( wwlm(GLOBAL_2D_ARRAY,2) )
#   endif

      ! Grid rotation angles
      allocate( angler_cos(GLOBAL_2D_ARRAY) )
      allocate( angler_sin(GLOBAL_2D_ARRAY) )

!# elif defined ANA_WEC_FRC && defined WEC
!      real, public :: wveb(GLOBAL_2D_ARRAY)
!      real, public :: wved(GLOBAL_2D_ARRAY)

# endif /* !ANA_WEC_FRC && !WKB_WWAVE */
#endif /* BBL || WEC */


      !#### end variables added as per forces.h #####

#ifdef AVERAGES
      if (wrt_file_avg) then                                         ! wrt_file_avg avoids any allocation
        if (wrt_SUP)   allocate(   sup_avg(GLOBAL_2D_ARRAY) )    ! if don't want avg file
        if (wrt_UST2D) allocate( ust2d_avg(GLOBAL_2D_ARRAY) )
        if (wrt_VST2D) allocate( vst2d_avg(GLOBAL_2D_ARRAY) )
#  ifdef SOLVE3D
        if (wrt_UST)   allocate(   ust_avg(GLOBAL_2D_ARRAY,N) )
        if (wrt_VST)   allocate(   vst_avg(GLOBAL_2D_ARRAY,N) )
        if (wrt_wST)   allocate(   wst_avg(GLOBAL_2D_ARRAY,N) )
        if (wrt_AKB)   allocate(   akb_avg(GLOBAL_2D_ARRAY,0:N) )
        if (wrt_AKW)   allocate(   akw_avg(GLOBAL_2D_ARRAY,0:N) )
        if (wrt_KVF)   allocate(   kvf_avg(GLOBAL_2D_ARRAY,N) )
        if (wrt_CALP)  allocate(  calp_avg(GLOBAL_2D_ARRAY) )
        if (wrt_KAPS)  allocate(kapsrf_avg(GLOBAL_2D_ARRAY) )
      endif
#  endif  /* SOLVE3D */
#endif /* AVERAGES */

      ! This section of code is taken from DH's WAVE_PACKET init_arrays.F

      do j=jstrR,jendR
        do i=istrR,iendR

! DevinD - NON-ANA old code L
#ifdef SOLVE3D
# if defined BBL || defined WEC
#  if !defined ANA_WWAVE && !defined WKB_WWAVE
          do k=1,2
            wwag(i,j,k)=init
            wwdg(i,j,k)=init
            wwpg(i,j,k)=init
#   ifdef BBL
            wwuob(i,j,k)=init
            wwvob(i,j,k)=init
#    ifdef WAVE_OFFLINE
            wwub(i,j,k) =init
#    endif
#   endif
#   ifdef WEC
            wwdrx(i,j,k)=init
            wwdre(i,j,k)=init
#    ifdef WAVE_OFFLINE
            wweb(i,j,k) =init
#     ifdef SURFACE_ROLLER
            wwqb(i,j,k) =init
#     endif
            wwUSg(i,j,k) = init
#    endif
#   endif /* WEC */
          enddo

          wfrq(i,j) =init
#   ifdef BBL
          uorb(i,j) =init
          vorb(i,j) =init
#   endif
#   ifdef WEC
          whrm(i,j) =init
          wdrx(i,j) =init
          wdre(i,j) =init
          wbst(i,j) =init
#    ifdef WAVE_OFFLINE
          wveb(i,j) =init
#     ifdef SURFACE_ROLLER
          wvqb(i,j) =init
#     endif
#    endif
#   endif /* WEC */
#  endif
# endif /* defined BBL || defined WEC */
#endif /* SOLVE3D */

      ! DevinD - even though they are split in the old code's init_arrays.F
      ! I don't see a need as it is the same loop so wave variables and
      ! stokes variables in the same loop

          sup(i,j)    =init
          ust2d(i,j)  =init
          vst2d(i,j)  =init
          brk2dx(i,j) =init
          brk2de(i,j) =init
          frc2dx(i,j) =init
          frc2de(i,j) =init
# ifdef SOLVE_3D
          calP(i,j)   =init
          Kapsrf(i,j) =init
          do k=1,N
#  ifndef SURFACE_BREAK
            brk3dx(i,j,k) =init
            brk3de(i,j,k) =init
#  endif
#  ifdef BODY_FRICTION
            frc3dx(i,j,k) =init
            frc3de(i,j,k) =init
#  endif
            ust(i,j,k) =init
            vst(i,j,k) =init
            wst(i,j,k) =init
            kvf(i,j,k) =init
          enddo
          do k=0,N
#  if defined LMD_KPP || defined LMD_BKPP
            Akb(i,j,k) =init
#  endif
            Akw(i,j,k) =init
            E_pre(i,j,k) =init
          enddo
# endif /* SOLVE_3D */
          wdsp(i,j) =init
#    ifdef SURFACE_ROLLER
          rdsp(i,j) =init
#    endif
        enddo
      enddo

      end subroutine init_arrays_wec_tile  !]

! ----------------------------------------------------------------------
#ifdef SOLVE3D
      subroutine wstokes_wec (tile)  ![
      ! vertical Stokes drift velocity at rho-point
      use hidden_mpi_vars
      use private_scratch
      implicit none
      integer :: tile

# include "compute_tile_bounds.h"
      call wstokes_wec_tile(istr,iend,jstr,jend)

      end subroutine wstokes_wec  !]

!-----------------------------------------------------------------------
      subroutine wstokes_wec_tile (istr,iend,jstr,jend)  ![

      implicit none
!! # include "param.h" ! already included at top of module
      integer istr,iend,jstr,jend, i,j,k
      real Wrk(0:N)
!! # include "grid" ! already included at top of module

      do j=jstr,jend
        do i=istr,iend
          Wrk(0)=0.D0
          do k=1,N,+1        !--> recursive
            Wrk(k) = Wrk(k-1) -Hz(i,j,k)
     &       *( pm(i,j)*(ust(i+1,j,k)-ust(i,j,k))
     &         +pn(i,j)*(vst(i,j+1,k)-vst(i,j,k)) )
          enddo

          wst(i,j,N)=+0.375*Wrk(N) +0.75*Wrk(N-1) -0.125*Wrk(N-2)
!          wst(i,j,N)=0. ! DevinD zeros
          do k=N-1,2,-1
            wst(i,j,k)=+0.5625*(Wrk(k  )+Wrk(k-1))
     &                 -0.0625*(Wrk(k+1)+Wrk(k-2))
!            wst(i,j,k)=0. ! DevinD zeros
          enddo
          wst(i,j,1)= -0.125*Wrk(2) +0.75*Wrk(1) +0.375*Wrk(0)
!          wst(i,j,1)= 0. ! DevinD zeros
        enddo
      enddo
# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then                       ! Set lateral
        do k=1,N                                   ! boundary
          do j=jstr,jend                           ! conditions
            wst(istr-1,j,k)=wst(istr,j,k)
          enddo
        enddo
      endif
      if (EASTERN_EDGE) then
        do k=1,N
          do j=jstr,jend
            wst(iend+1,j,k)=wst(iend,j,k)
          enddo
        enddo
      endif
# endif
# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do k=1,N
          do i=istr,iend
            wst(i,jstr-1,k)=wst(i,jstr,k)
          enddo
        enddo
      endif
      if (NORTHERN_EDGE) then
        do k=1,N
          do i=istr,iend
            wst(i,jend+1,k)=wst(i,jend,k)
          enddo
        enddo
      endif
#  ifndef EW_PERIODIC
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          wst(istr-1,jstr-1,k)=wst(istr,jstr,k)
        enddo
      endif
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          wst(istr-1,jend+1,k)=wst(istr,jend,k)
        enddo
      endif
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          wst(iend+1, jstr-1,k)=wst(iend,jstr,k)
        enddo
      endif
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          wst(iend+1,jend+1,k)=wst(iend,jend,k)
        enddo
      endif
#  endif
# endif

      end subroutine wstokes_wec_tile  !]
#endif  /* SOLVE3D */

! ----------------------------------------------------------------------
      subroutine compute_wec_rhs_uv_terms (istr,iend,jstr,jend,ru,rv,wrkone,wrktwo)  ![
      ! add some of the WEC variables to the RHS
      ! of the equation system
      implicit none

      ! inputs
      integer,intent(in) :: istr, iend, jstr, jend
      real,dimension(PRIVATE_2D_SCRATCH_ARRAY,N) :: ru,rv
      real,dimension(PRIVATE_2D_SCRATCH_ARRAY)   :: wrkone,wrktwo

      ! local
      integer :: i, j, k
      real    :: cff, cff1, gam1, gam2

# include "compute_auxiliary_bounds.h" /* Taken from top of rhs3d35S.F */
! VF_ROBUST is not undefined at end of rhs3d35S.F, hence not undefined here.
# define VF_ROBUST /* Taken from top of rhs3d35S.F */


      ! Add combined vortex-force and advection terms & breaking terms
      ! 2nd- and 4th-order centered schemes, or 2nd with local weighted
      ! filter (with VF_ROBUST) are available.

# if defined  WEC && defined UV_ADV
#  define utmp wrkone
#  define vtmp wrktwo
#  ifdef VF_ROBUST
      gam1=0.1        ! local weighted filter (gam1 may be 0-1/3)
      gam2=1.-3.*gam1  ! gam1=0 to revert back to 2nd-order centered
c      gam1=-1./24.     ! 4th-order centered scheme
c      gam2=9./8.
#  endif
      do k=1,N                       ! <-- start k-loop
#  ifdef VF_ROBUST
        do j=jstrV-2,jend+1
          do i=istrU-2,iend+1
            vtmp(i,j) =v(i,j,k,nrhs)
          enddo
        enddo
        do j=jstrV-2,jend+1
          do i=istrU-2,iend+1
            utmp(i,j) =u(i,j,k,nrhs)
          enddo
        enddo
#   ifndef EW_PERIODIC
        if (WESTERN_EDGE) then
          do j=jstrV-2,jend+1
            utmp(istrU-2,j) =utmp(istrU-1,j)
          enddo
          do j=jstr,jend+1
            vtmp(istrU-2,j) =vtmp(istrU-1,j)
          enddo
        endif
#   endif
#   ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do i=istr,iend+1
            utmp(i,jstrV-2) =utmp(i,jstrV-1)
          enddo
          do i=istrU-2,iend+1
            vtmp(i,jstrV-2) =vtmp(i,jstrV-1)
          enddo
        endif
#   endif
#  endif
        do j=jstr,jend
          do i=istrU,iend
            cff  =0.5*dn_u(i,j)*(Hz(i-1,j,k)+Hz(i,j,k))
            cff1 = 0.25*(
#  ifdef VF_ROBUST
     &                    4.*ust(i,j,k)*( gam1*(utmp(i+1,j)-utmp(i-2,j))
     &                                   +gam2*(utmp(i,j)-utmp(i-1,j)) )
     &        +(vst(i,j,k)+vst(i-1,j,k))
     &                                 *( gam1*(vtmp(i+1,j)-vtmp(i-2,j))
     &                                   +gam2*(vtmp(i,j)-vtmp(i-1,j)) )
     &        +(vst(i,j+1,k)+vst(i-1,j+1,k))
     &                             *( gam1*(vtmp(i+1,j+1)-vtmp(i-2,j+1))
     &                               +gam2*(vtmp(i,j+1)-vtmp(i-1,j+1)) )
#  else
     &         (ust(i+1,j,k)+ust(i,j,k))*(u(i+1,j,k,nrhs)-u(i,j,k,nrhs))
     &        +(ust(i,j,k)+ust(i-1,j,k))*(u(i,j,k,nrhs)-u(i-1,j,k,nrhs))
     &        +(vst(i,j,k)+vst(i-1,j,k))*(v(i,j,k,nrhs)-v(i-1,j,k,nrhs))
     &        +(vst(i,j+1,k)+vst(i-1,j+1,k))
     &                              *(v(i,j+1,k,nrhs)-v(i-1,j+1,k,nrhs))
#  endif

     &                                                                 )
            ru(i,j,k)=ru(i,j,k) + cff*( cff1
#  ifndef SURFACE_BREAK
     &                                          +dm_u(i,j)*brk3dx(i,j,k)
#  endif
#  if defined WAVE_FRICTION && defined BODY_FRICTION
     &                                          +dm_u(i,j)*frc3dx(i,j,k)
#  endif
     &                                                                 )
          enddo
        enddo

        do j=jstrV,jend
          do i=istr,iend
            cff  =0.5*dm_v(i,j)*(Hz(i,j-1,k)+Hz(i,j,k))
            cff1 = 0.25*(
#  ifdef VF_ROBUST
     &         (ust(i,j,k)+ust(i,j-1,k))
     &                                 *( gam1*(utmp(i,j+1)-utmp(i,j-2))
     &                                   +gam2*(utmp(i,j)-utmp(i,j-1)) )
     &        +(ust(i+1,j,k)+ust(i+1,j-1,k))
     &                             *( gam1*(utmp(i+1,j+1)-utmp(i+1,j-2))
     &                               +gam2*(utmp(i+1,j)-utmp(i+1,j-1)) )
     &                   +4.*vst(i,j,k)*( gam1*(vtmp(i,j+1)-vtmp(i,j-2))
     &                                   +gam2*(vtmp(i,j)-vtmp(i,j-1)) )
#  else
     &         (ust(i,j,k)+ust(i,j-1,k))*(u(i,j,k,nrhs)-u(i,j-1,k,nrhs))
     &        +(ust(i+1,j,k)+ust(i+1,j-1,k))
     &                              *(u(i+1,j,k,nrhs)-u(i+1,j-1,k,nrhs))
     &        +(vst(i,j+1,k)+vst(i,j,k))*(v(i,j+1,k,nrhs)-v(i,j,k,nrhs))
     &        +(vst(i,j,k)+vst(i,j-1,k))*(v(i,j,k,nrhs)-v(i,j-1,k,nrhs))
#  endif
     &                                                                 )
            rv(i,j,k)=rv(i,j,k) + cff*( cff1
#  ifndef SURFACE_BREAK
     &                                          +dn_v(i,j)*brk3de(i,j,k)
#  endif
#  if defined WAVE_FRICTION && defined BODY_FRICTION
     &                                          +dn_v(i,j)*frc3de(i,j,k)
#  endif
     &                                                                 )
          enddo
        enddo
      enddo                ! <-- end k-loop
#  undef utmp
#  undef vtmp
# endif
# undef VF_ROBUST


      end subroutine compute_wec_rhs_uv_terms  !]

!-----------------------------------------------------------------------
      subroutine set_wec_frc_tile (istr,iend,jstr,jend,ierr)  ![
      ! set wec forces (realistic or analytical forcing)
      implicit none

      ! input/outputs
      integer,intent(in) :: istr,iend,jstr,jend
      integer,intent(out):: ierr
      ! local
      integer :: tmp,i,j
      real    :: cff1,cff2

! Needed for istrR,iendR, etc
# include "compute_extended_bounds.h"

      ! IF ANALYTICAL FORCING
      ! =====================
# ifdef ANA_WEC_FRC
      call ana_wec_frc_tile(istr,iend,jstr,jend)

# else /* REALISTIC MODEL */

      ! IF REALISTIC FORCING
      ! ====================

      ! note: wec_time(it2) > wec_time(it1), but values of it1
      ! and it2 change between each reading of values so that only one
      ! value is read in.

      ! UPDATE LATEST READ IN VARIABLE?
      ! -------------------------------
      ! Need wec_tim(it1) < time < wec_tim(it2), thus:
      if (wec_tim(it2) < time) then  ! We need to refresh data

        ! FIRST TIMESTEP
        ! --------------
        ! wec_tim(it1) & wec_tim(its) declared < 0 at top of module:
        if (wec_tim(it1) < 0) then ! First timestep

          ! ifile and irec should be in the module pre-amble, initial as 0
          ! Find the force file (ifile) and the record (irec) for
          ! the last record with wec_time < tday
          ! First time, we need to read an extra slice of data
          ! so this routine provides irec = irec - 1 for first time.
          ! Must make sure model time format (secs/days) matches data time
          ! Hence tdays as input data times are in days
          call find_new_record(wec_awv_name,wec_tim_name,tdays,ifile,irec,wec_tim(it1))

          ! Reading is only done by master processor hence C$OMP MASTER
C$OMP MASTER
          ! For the first time irec = irec - 1, as above, such that
          ! closest record to wec_tim(it1) < time.
          ! read in raw input variables and
          ! read and set input time in wec_tim(it1) to seconds
          call read_wec_frc(frcfile(ifile),irec,it1,wec_tim(it1))
C$OMP END MASTER
C$OMP BARRIER
          ! Barrier needed because master processor must have read in all
          ! WEC data before child processors try to access read in data!

          ! Confirm to terminal reading of variables and
          ! time for read in data
          if (mynode == 0) then
            write(*,'(4x,A,5x,A,G14.6,1x,A,I4)')
     &               'set_wec_frc :: read WEC inputs',
     &               'wec_input_time (s) =',wec_tim(it1),
     &               'rec =',irec MYID
          endif

          ! Convert read inputs into useable, rotated & rmask'ed variables
          ! for the first, earlier read in time variables(time) < time.
          call set_wec_inp2vars (istr,iend,jstr,jend,it1)

        ! ALL OTHER TIMESTEPS
        ! -------------------
        ! Every time step switch indices such that
        ! var(i,j,1) & wec_tim(1) -> earlier && var(i,j,2) & wec_tim(2) -> later: step t=x
        ! var(i,j,1) & wec_tim(1) -> later   && var(i,j,2) & wec_tim(2) ->  earlier: step t=x+1
        ! and so on. The benifit is only need to change one variable each time.
        else
          tmp = it1
          it1 = it2 ! If it1 = 1, it now equals 2, and vice-versa
          it2 = tmp ! If it2 = 2, it now equals 1, and vice-versa

        endif ! End of first step

        ! Find the force file (ifile) and the record (irec) for
        ! the first record with wec_time > tday
        ! Must make sure model time format (secs/days) matches data time
        ! Hence tdays as input data times are in days
        call find_new_record(wec_awv_name,wec_tim_name,tdays,ifile,irec,wec_tim(it2))

             ! Causes heisenbug
!             if (mynode==0) print *,'irec=',irec ! Debug

        ! Reading is only done by master processor hence C$OMP MASTER
C$OMP MASTER
        ! Read in the first record with wec_time > tday
        ! Remember, it2 alternates between index 1 & 2 for each reading
        ! read in raw input variables and
        ! read and set input time in wec_tim(it2) to second
        call read_wec_frc(frcfile(ifile),irec,it2,wec_tim(it2))
C$OMP END MASTER
C$OMP BARRIER

        ! Confirm to terminal reading of variables and
        ! time for read in data
        if (mynode == 0) then                  ! Text formatting:
          write(*,'(4x,A,5x,A,G14.6,1x,A,I4)') ! 4x is spaces, A is string
     &               'set_wec_frc :: read WEC inputs', ! 5x is 5 spaces, A is string
     &               'wec_input_time (s) =',wec_tim(it2), ! ! G is number for time
     &               'rec =',irec MYID ! I is integer
        endif

        ! Convert read inputs into useable, rotated & rmask'ed variables
        ! for the latter read in time variables(time) > time.
        call set_wec_inp2vars (istr,iend,jstr,jend,it2)

      endif

      ! Temporal interpolation
      ! ======================

      ! Set coefficients for interpolation.
      ! Time difference in seconds. Not yet ratios of 1.
      ! cff1 with it2 because if time=it1 then difference=0, but want cff1 = 1.
      cff1=( wec_tim(it2)-time ) / ( wec_tim(it2)-wec_tim(it1) )
      cff2=( time-wec_tim(it1) ) / ( wec_tim(it2)-wec_tim(it1) )

      ! Check the interpolation coefficients such that modeltime is
      ! bounded by read in variable times.
      if (cff1.ge.0. .and. cff2.ge.0.) then

        do j=jstrR,jendR
          do i=istrR,iendR

            ! WAVE PARAMETERS - at rho-points
            ! ---------------

            ! wwag already converted to 2 x wwag = whrm in set_wec_inp2vars
            whrm(i,j)=cff1*wwag(i,j,it1)+cff2*wwag(i,j,it2)
            wdrx(i,j)=cff1*wwdrx(i,j,it1)+cff2*wwdrx(i,j,it2)
            wdre(i,j)=cff1*wwdre(i,j,it1)+cff2*wwdre(i,j,it2)
            ! wwpg already converted to peak per. to freq. in set_wec_inp2vars
            wfrq(i,j)=cff1*wwpg(i,j,it1)+cff2*wwpg(i,j,it2)

            ! STOKES DRIFT VARIABLES - at rho-points
            ! ----------------------

            ! Get x & y components of wwUSg and wwVSg
            ! Rotation components angler_cos & angler_sin
            ! are precomputed in set_wec_inp2vars
            ustr2d(i,j)=(cff1*wwUSg(i,j,it1)+cff2*wwUSg(i,j,it2))
     &                                         *angler_cos(i,j)
     &                +(cff1*wwVSg(i,j,it1)+cff2*wwVSg(i,j,it2))
     &                                         *angler_sin(i,j) ! U-east to xi-dir depth ave. Stokes drift
            vstr2d(i,j)=-(cff1*wwUSg(i,j,it1)+cff2*wwUSg(i,j,it2))
     &                                         *angler_sin(i,j)
     &                +(cff1*wwVSg(i,j,it1)+cff2*wwVSg(i,j,it2))
     &                                         *angler_cos(i,j) ! V-north to eta-dir depth ave. Stokdes drift

            ustr0(i,j)=(cff1*wwus0g(i,j,it1)+cff2*wwus0g(i,j,it2))
     &                                         *angler_cos(i,j)
     &                +(cff1*wwvs0g(i,j,it1)+cff2*wwvs0g(i,j,it2))
     &                                         *angler_sin(i,j) ! U-east to xi-dir depth ave. Stokes drift
            vstr0(i,j)=-(cff1*wwus0g(i,j,it1)+cff2*wwus0g(i,j,it2))
     &                                         *angler_sin(i,j)
     &                +(cff1*wwvs0g(i,j,it1)+cff2*wwvs0g(i,j,it2))
     &                                         *angler_cos(i,j) ! V-north to eta-dir depth ave. Stokdes drift

            wveb(i,j)=(cff1*wweb(i,j,it1)+cff2*wweb(i,j,it2))
!#   ifdef WAVE_FRICTION
            wved(i,j)=(cff1*wwed(i,j,it1)+cff2*wwed(i,j,it2))
!#   else
!            wved(i,j)=0.0
!#   endif

!#  if defined WAVE_OFFLINE && defined SURFACE_ROLLER
!            wvqb(i,j)=cff1*wwqb(i,j,it1)+cff2*wwqb(i,j,it2)
!#  endif
!#  if defined WAVE_OFFLINE && defined SUP_OFF
            sup(i,j)= cff1*wwsup(i,j,it1)+cff2*wwsup(i,j,it2)
!#  endif
!#  if defined WAVE_OFFLINE && defined WAVE_DIFF_FROM_LM
            llmm(i,j)= cff1*wwlm(i,j,it1)+cff2*wwlm(i,j,it2)
!#  endif

          enddo
        enddo

      ! Throw warning if interpolation times do not bound model time.
      elseif (ZEROTH_TILE) then
        write(*,'(/1x,3A/3(1x,A,F16.10)/)')       '### WARNING: ',
     &  'set_wec_frc_tile :: Current model time is outside bounds of ',
     &  '''wec_time''.',  'start =', wec_tim(it1)*sec2day,
     &                    'tdays =',     tdays,
     &                    'end =',   wec_tim(it2)*sec2day
        ierr=ierr+1
      endif
# endif /* ANA_WEC_FRC or REALISTIC MODEL */

      end subroutine set_wec_frc_tile  !]

!-------------------------------------------------------------------------
#ifndef ANA_WEC_FRC
      subroutine set_wec_inp2vars (istr,iend,jstr,jend,idx)  ![
      ! convert and mask read in wec data
      ! partial copy of set_wwave_tile from get_wwave.F in old code.
      implicit none

      ! input
      integer :: istr,iend,jstr,jend
      integer :: idx                   ! index relating to read in variable data time set 1 or set 2
      ! local
      integer :: i,j
      real    :: wdir                  ! interim variable to transform wave direction

# include "compute_extended_bounds.h"


      ! Compute bed wave orbital velocity (uorb,vorb) and wave frequency,
      ! wfrq, used with BBL, and/or Stokes transport (wvstx,wvsty) and wave
      ! set-up, wvsup, for MRL_WCI.  Dispersion relation is approximated by
      ! an explicit polynomial function proposed by Dean & Dalrymple (1991).

      do j=jstrR,jendR
        do i=istrR,iendR

          ! WAVE PARAMETERS
          ! ---------------

          ! Wave amplitude is read in, so convert to hrms for later use.
          ! 2xa = hrms. Keep same variable name to prevent cost
          ! of extra variable
          wwag(i,j,idx) =2.0*wwag(i,j,idx)
# ifdef MASKING
     &                           *rmask(i,j)
# endif

              ! Calculate wdir once as used twice below
          wdir=-(deg2rad*wwdg(i,j,idx)-1.5*pi)  ! FROM N (WAVE_OFFLINE) to travel to 0deg --> x-dir
# ifdef CURVGRID
     &                    -angler(i,j)          ! angler must be in radians
# endif
              ! Direction (wwdg) is read in as one variable so can't be
              ! reused for dir-x and dir-y, need new interim variables
          wwdrx(i,j,idx) =cos(wdir)
# ifdef MASKING
     &                             *rmask(i,j)
# endif
          wwdre(i,j,idx) =sin(wdir)
# ifdef MASKING
     &                             *rmask(i,j)
# endif

          wwpg(i,j,idx)=2.0*pi/max(wwpg(i,j,idx),0.1D0) !  peak per. to freq.
# ifdef MASKING
     &                             *rmask(i,j)
# endif


# ifdef MASKING /* remaining variables not converted just rmask'ed */

          wwUSg(i,j,idx) = wwUSg(i,j,idx)*rmask(i,j)       ! DevinD this is likely unneccessary
          wwVSg(i,j,idx) = wwVSg(i,j,idx)*rmask(i,j)
          wwus0g(i,j,idx) = wwus0g(i,j,idx)*rmask(i,j)
          wwvs0g(i,j,idx) = wwvs0g(i,j,idx)*rmask(i,j)

          wweb(i,j,idx) = wweb(i,j,idx)*rmask(i,j)
!#   if defined WAVE_OFFLINE && defined WAVE_FRICTION
          wwed(i,j,idx) = wwed(i,j,idx)*rmask(i,j)
!#   endif
!#   if defined WAVE_OFFLINE && defined SUP_OFF
          wwsup(i,j,idx) = wwsup(i,j,idx)*rmask(i,j)
!#   endif
!#   if defined WAVE_OFFLINE && defined WAVE_DIFF_FROM_LM
          wwlm(i,j,idx) = wwlm(i,j,idx)*rmask(i,j)
!#   endif
# endif

          ! Set grid rotation angles (if first timestep)
          ! ------------------------

          ! These are used in set_wec_frc_tile to get
          ! x & y forcing components from read in data.
          ! Need once as doesn't change.
          if(FIRST_TIME_STEP) then

            angler_cos(i,j) = cos(angler(i,j))
            angler_sin(i,j) = sin(angler(i,j))

          endif

        enddo
      enddo

      end subroutine set_wec_inp2vars  !]
#endif /* !ANA_WEC_FRC */

! ----------------------------------------------------------------------
#ifdef AVERAGES
      subroutine calc_wec_avg  ![
      implicit none

      ! local
      real :: coef

      navg_wec = navg_wec +1
      coef = 1./navg_wec

      ! need i0:i1 indices because arrays still GLOBAL_2D therefore wasted margin
      if (wrt_SUP)
     &    sup_avg(i0:i1,j0:j1)  =   sup_avg(i0:i1,j0:j1)*(1-coef) +   sup(i0:i1,j0:j1)*coef
      if (wrt_UST2D)
     &  ust2d_avg( 1:i1,j0:j1)  = ust2d_avg( 1:i1,j0:j1)*(1-coef) + ust2d( 1:i1,j0:j1)*coef
      if (wrt_VST2D)
     &  vst2d_avg(i0:i1, 1:j1)  = vst2d_avg(i0:i1, 1:j1)*(1-coef) + vst2d(i0:i1, 1:j1)*coef
#  ifdef SOLVE3D
      if (wrt_UST)
     &    ust_avg( 1:i1,j0:j1,:)= ust_avg( 1:i1,j0:j1,:)*(1-coef) + ust( 1:i1,j0:j1,:)*coef
      if (wrt_VST)
     &    vst_avg(i0:i1, 1:j1,:)= vst_avg(i0:i1, 1:j1,:)*(1-coef) + vst(i0:i1, 1:j1,:)*coef
      if (wrt_WST)
     &    wst_avg(i0:i1,j0:j1,:)= wst_avg(i0:i1,j0:j1,:)*(1-coef) + wst(i0:i1,j0:j1,:)*coef
      if (wrt_Akb)
     &    akb_avg(i0:i1,j0:j1,:)= akb_avg(i0:i1,j0:j1,:)*(1-coef) + akb(i0:i1,j0:j1,:)*coef
      if (wrt_Akw)
     &    akw_avg(i0:i1,j0:j1,:)= akw_avg(i0:i1,j0:j1,:)*(1-coef) + akw(i0:i1,j0:j1,:)*coef
      if (wrt_KVF)
     &    kvf_avg(i0:i1,j0:j1,:)= kvf_avg(i0:i1,j0:j1,:)*(1-coef) + kvf(i0:i1,j0:j1,:)*coef
      if (wrt_CALP)
     &   calp_avg(i0:i1,j0:j1)  =calp_avg(i0:i1,j0:j1)  *(1-coef) +calp(i0:i1,j0:j1)  *coef
      if (wrt_KAPS)
     & kapsrf_avg(i0:i1,j0:j1)  =kapsrf_avg(i0:i1,j0:j1)*(1-coef) +kapsrf(i0:i1,j0:j1)*coef
#  endif  /* SOLVE3D */

      end subroutine calc_wec_avg  !]
#endif /* AVERAGES */

!-----------------------------------------------------------------------
      subroutine wec_forces  ![
      ! this is mrl_wci (tile) from mrl_wci.F in old code.
      use hidden_mpi_vars
      use private_scratch
      implicit none

      integer,save :: tile=0
# include "compute_tile_bounds.h"

      call wec_forces_tile (istr,iend,jstr,jend,
# ifdef SOLVE3D
     &                              A3d(1,1),A3d(1,2),A3d(1,3),A3d(1,4),
# endif
     &            A2d(1,1),A2d(1,2),A2d(1,3),A2d(1,4),A2d(1,5),A2d(1,6),
     &                   A2d(1,7),A2d(1,8),A2d(1,9),A2d(1,10),A2d(1,11),
     &                                              A2d(1,12),A2d(1,13))

      end subroutine wec_forces  !]
! ----------------------------------------------------------------------
      subroutine wec_forces_tile (istr,iend,jstr,jend,  ![
# ifdef SOLVE3D
     &                                              wrk1,wrk2,wrk3,wrk7,
# endif
     &           wh,fr,kw,brk,stk,Dstp,act,kD,inv_d,inv_f,frc,ebrk,erol)

      ! NOTE: input wave data is at rho-points, thus conversion from
      ! whrm at rho-points, is needed to get ust2d and u-points, etc.

      ! This is taken from mrl_wci.F of old code.

!   Evaluating wave-averaged terms and Stokes drift based on MRL04
!   (see: McWilliams, Restrepo & Lane, 2004, JFM, 511, pp.135-178)
!   inputs: R.M.S. wave height (m);
!           peak wave freqiuency (rad/s);
!           mean wave direction: wdrx & wdre (non dimensional);
!           breaking dissipation   (\epsilon_b / \rho, m3/s3);
!           roller dissipation     (\epsilon_r / \rho, m3/s3); and
!           frictional dissipation (\epsilon_d / \rho, m3/s3)
!
!   Note that if wkb_wwave.F is used, all the dissipation terms are
!   divided by wave frequency, sigma (2pi/T).
!

      use coupling

      implicit none

!# include "param"
!# include "grid"
!# include "scalars"
# ifdef SOLVE3D
!#  if defined WAVE_FRICTION && defined BODY_FRICTION
!#   include "mixing"
!#  endif
# endif
!# ifdef ANA_BRY
!#   include "boundary"
!# endif
!# ifdef WKB_WWAVE
!#   include "wkb_wwave.h"
!# endif
!----- Cigdem (Uorb => Uorb1)
      integer istr,iend,jstr,jend,i,j,imin,jmin
      real    cff,cff1,cff2,cff3,cff4,inv_g,khd,kh,umag,
     &   eps,wave_ramp,wramp2,kbrk,fb,fb0,fb1,fb2,inv_fbs,
     &   bconst,tauc,tauw,Uorb1,ka_f00,hz0,z_tide,
     %   dtinv,inv_zb,a_brk,a_kv,c1o3,c4o3,inv_k,dd,ust_ker,intfb,
     &   fbsrf,fn1,fn2,khmax,bz1,bz2,beta,abot,delta,vst_ker,
     &   a_frc,kfrc,
     &             wh(PRIVATE_2D_SCRATCH_ARRAY),
     &             fr(PRIVATE_2D_SCRATCH_ARRAY),
     &             kw(PRIVATE_2D_SCRATCH_ARRAY),
     &            brk(PRIVATE_2D_SCRATCH_ARRAY),
     &            stk(PRIVATE_2D_SCRATCH_ARRAY),
     &           Dstp(PRIVATE_2D_SCRATCH_ARRAY),
     &            act(PRIVATE_2D_SCRATCH_ARRAY),
     &             kD(PRIVATE_2D_SCRATCH_ARRAY),
     &          inv_d(PRIVATE_2D_SCRATCH_ARRAY),
     &          inv_f(PRIVATE_2D_SCRATCH_ARRAY),
     &            frc(PRIVATE_2D_SCRATCH_ARRAY),
     &           ebrk(PRIVATE_2D_SCRATCH_ARRAY),
     &           erol(PRIVATE_2D_SCRATCH_ARRAY)
# if defined WAVE_OFFLINE && !defined U2010
     &         ,usec,vsec,USEZ,VSEZ,ustm0w
     &         ,FAC,FACo,COUNTR, USTERR, ustm0
     &         ,fac2d(PRIVATE_2D_SCRATCH_ARRAY)
     &         ,inv_exef(PRIVATE_2D_SCRATCH_ARRAY)
     &         ,keff(PRIVATE_2D_SCRATCH_ARRAY)
# endif
# if defined SUP_OFF
     &         ,supm, Fsup
     &         ,actf(PRIVATE_2D_SCRATCH_ARRAY)
# endif
# if defined WAVE_DIFF_FROM_LM
     &         ,keffkw
# endif
# ifdef SOLVE3D
      integer  k,kk
      real     wrk1(PRIVATE_2D_SCRATCH_ARRAY,0:N),
     &         wrk2(PRIVATE_2D_SCRATCH_ARRAY,0:N),
     &         wrk3(PRIVATE_2D_SCRATCH_ARRAY,0:N),
     &         wrk4(0:N),wrk5(0:N),wrk6(0:N), kvsurf
     &         ,wrk7(PRIVATE_2D_SCRATCH_ARRAY,0:N)
# endif
      parameter ( eps=1.e-10,
!# ifdef STOKESEKMAN
!     &            bconst=0.D0,
!# else
     &            bconst=0.1D0, !0.03D0,               ! breaking contribution to KPP
!# endif
     &            a_kv=1.2D0, !1.2D0                 ! breaking scale for eddy visc.
     &            a_brk=0.2D0,                ! breaking scale for body force
     &            a_frc=3.D0,                 ! friction scale for body force
     &            khmax=20.D0,                ! deep-water limit for k x dep
     &            c1o3=0.3333333333333333,    ! 1/3
     &            c4o3=1.3333333333333333     ! 4/3
     &           )
# ifdef WAVE_OFFLINE
      real roller_fraction
!#  ifdef SURFACE_ROLLER
!      parameter ( roller_fraction=1.0)
!#  else
      parameter ( roller_fraction=0.0)
!#  endif
#  define wkb_roller roller_fraction
# endif
# include "compute_auxiliary_bounds.h"


! TRANSLATIONS AND MANIPULATIONS OF WEC VARIABLES FROM RAW VALUES IN
! FORCING FILE TO NEEDED VARIABLES HERE. AS PER SET_WWAVE_TILE (OLD) OR
! SET_WEC_FRC_TILE (NEW).

! Small section taken from old code's set_wwave_tile subroutine
! Convert offline 2d stokes drift at rho-points to u-points.
! ==========================================================

      do j=jstrR,jendR
        do i=istr,iendR
          ust2d(i,j)=0.5*(ustr2d(i-1,j)+ustr2d(i,j))
# ifdef MASKING
!    &                  *umask(i,j)
# endif
        enddo
      enddo
      do j=jstr,jendR
        do i=istrR,iendR
          vst2d(i,j)=0.5*(vstr2d(i,j-1)+vstr2d(i,j))
# ifdef MASKING
!    &                  *vmask(i,j)
# endif
        enddo
      enddo



! AS PER MRL_WCI_TILE below
! =========================

! explicit wavenumber estimator

# undef KH_SOULSBY
# define KH_HUNT

! vertical distribution function F_Kv(z) for eddy viscosity

c# define FKV_FUNC1
# define FKV_FUNC2
c# define FKV_FUNC3

! vertical distribution function F_B(z) for breaking accerelation

# undef FB_FUNC0
c# define FB_FUNC1
c# define FB_FUNC2
# define FB_FUNC3
!----- Cigdem
c# define FB_WSCALE
c# define BRK0
c# define KVF0
c# define AKB0

! vertical distribution function F_B(z) for bottom streaming

c# define FF_FUNC1
c# define FF_FUNC2
# define FF_FUNC3

! initial ramping coefficients

!# if defined WAVE_RAMP
!#  ifdef STOKESEKMAN
!      wave_ramp = tanh(3.0*dt*sec2day*float(iic-ntstart))
!#  elif defined DUCK94
!      wave_ramp = tanh(24.0*dt*sec2day*float(iic-ntstart))
!#  elif defined MVCO
!      wave_ramp = tanh(3.0*dt*sec2day*float(iic-ntstart))
!!----- Cigdem (added RIP_CURRENT option)
!#  elif defined RIP_CURRENT
!c      wave_ramp = tanh(144.0*dt*sec2day*float(iic-ntstart)) ! 30 min
!      wave_ramp = tanh(432.0*dt*sec2day*float(iic-ntstart))  ! 10 min
!#  else
!      wave_ramp = tanh(24.0*dt*sec2day*float(iic-ntstart))
!#  endif
!# else
      wave_ramp = 1.0
!# endif
      wramp2= wave_ramp**2
!
! Evaluate 2DH wave-current intreraction variable.
! ================================================
!
! 2DH Stokes velocities, breaking, roller and bottom-friction
! dissipation terms, defined at horizontal rho-points.
!
      inv_g = 1.0/g
      do j=jstr-1,jend+1
        do i=istr-1,iend+1
          Dstp(i,j)=max(h(i,j)+zeta(i,j,knew),eps)
          inv_d(i,j)=1.0/Dstp(i,j)

          fr(i,j) =wfrq(i,j)
          wh(i,j) =wave_ramp*whrm(i,j)
          inv_f(i,j) =1.0/max(fr(i,j),eps)
          khd = Dstp(i,j)*(fr(i,j)**2)*inv_g

          kh = sqrt( khd*khd + khd/(1.0 + khd*(0.6666666666
     &               +khd*(0.3555555555 + khd*(0.1608465608
     &               +khd*(0.0632098765 + khd*(0.0217540484
     &                            +khd*0.0065407983)))))) )

          kD(i,j) =kh
          kw(i,j) =kh*inv_d(i,j)
          ebrk(i,j) =wramp2*wveb(i,j)            ! epsilon_b (m3/s3)

          cff =ebrk(i,j)

          frc(i,j) =wramp2*wved(i,j)*kw(i,j)
     &                           *inv_f(i,j)     ! 2DH friction term
          brk(i,j) =cff*kw(i,j)*inv_f(i,j)       ! 2DH breaking term
          act(i,j) =inv_f(i,j)*(
     &                      0.125*g*(wh(i,j)**2) ) ! wave action density

          stk(i,j) =act(i,j)*kw(i,j)*inv_d(i,j)  ! depth-averaged Stokes
# ifdef MASKING
           kw(i,j) = kw(i,j)*rmask(i,j)
           kD(i,j) = kD(i,j)*rmask(i,j)
          act(i,j) =act(i,j)*rmask(i,j)
          stk(i,j) =stk(i,j)*rmask(i,j)
          brk(i,j) =brk(i,j)*rmask(i,j)
          frc(i,j) =frc(i,j)*rmask(i,j)
# endif   /* apply land mask */
        enddo
      enddo              ! <-- discard inv_d (ifndef BBL_F00)

# ifdef BRK0
      do j=jstr-1,jend+1
        do i=istr-1,iend+1
          brk(i,j)=0.D0
        enddo
      enddo
# endif

!
! sup : quasi-static sea-level response, set-up. no interaction.
! ==============================================================
!

# if !defined SUP_OFF
      do j=jstrR,jendR
        do i=istrR,iendR
          sup(i,j) =-0.125*(wh(i,j)**2)*kw(i,j)
     &              /max(sinh(2.*min(kD(i,j),khmax)),eps)
#  ifdef MASKING
     &                                        *rmask(i,j)
#  endif
        enddo
      enddo
# endif

# if defined SUP_OFF
      actf(:,:)=0.
! determine the effective filtered action density g a^2/w
      do j=jstrR,jendR
        do i=istrR,iendR
          supm=-0.125*(wh(i,j)**2)*kw(i,j)
     &              /max(sinh(2.*min(kw(i,j)*h(i,j),khmax)),eps)
! re calculate sup with \hat{\eta}
          if (kD(i,j).lt.6.and.supm.ne.0.) then
          Fsup=min(1.,sqrt(sup(i,j)/supm))
          else
          Fsup=1.
          endif
          sup(i,j) =-0.125*(wh(i,j)**2)*kw(i,j)*Fsup**2
     &              /max(sinh(2.*min(kD(i,j),khmax)),eps)
!          sup=0. ! DevinD zeros
          actf(i,j)=act(i,j)*Fsup**2
        enddo
      enddo
# endif

!
! Dissipation terms : copy into shared arrays
! ===========================================
!
      do j=jstrR,jendR
        do i=istrR,iendR
          wdsp(i,j)=ebrk(i,j)         ! epsilon_b (m3/s3)
!# elif !defined WKB_WWAVE & defined WAVE_OFFLINE
          cff = wramp2*wved(i,j)*inv_f(i,j)
          wdrg(i,j)=cff*fr(i,j)       ! epsilon_d (m3/s3)
          wbst(i,j)=cff*kw(i,j)       ! e_d stress (m2/s2) for KPP
# ifdef MASKING
          wdsp(i,j)=wdsp(i,j)*rmask(i,j)
          wbst(i,j)=wbst(i,j)*rmask(i,j)
          wdrg(i,j)=wdrg(i,j)*rmask(i,j)
# endif
        enddo
      enddo
!# ifndef WKB_WWAVE
#  undef inv_f
!# endif
!
! 2DH depth-averaged Stokes drift and breaking acceleration
! =========================================================
!
      do j=jstrR,jendR
        do i=istr,iendR
!# if defined WKB_WWAVE || defined U2010
!          ust2d(i,j)  =0.5*(stk(i-1,j)*wdrx(i-1,j)+stk(i,j)*wdrx(i,j))
!# endif
          brk2dx(i,j) =0.5*(brk(i-1,j)*wdrx(i-1,j)+brk(i,j)*wdrx(i,j))
          frc2dx(i,j) =0.5*(frc(i-1,j)*wdrx(i-1,j)+frc(i,j)*wdrx(i,j))
#  ifdef MASKING
!#   ifdef SHOREFACE
!          ust2d(i,j)  = ust2d(i,j)*umask(i,j)
!#   endif
          brk2dx(i,j) =brk2dx(i,j)*umask(i,j)
          frc2dx(i,j) =frc2dx(i,j)*umask(i,j)
#  endif
        enddo
      enddo
      do j=jstr,jendR
        do i=istrR,iendR
!# if defined WKB_WWAVE || defined U2010
!          vst2d(i,j)  =0.5*(stk(i,j-1)*wdre(i,j-1)+stk(i,j)*wdre(i,j))
!# endif
          brk2de(i,j) =0.5*(brk(i,j-1)*wdre(i,j-1)+brk(i,j)*wdre(i,j))
          frc2de(i,j) =0.5*(frc(i,j-1)*wdre(i,j-1)+frc(i,j)*wdre(i,j))
#  ifdef MASKING
!#   ifdef SHOREFACE
!          vst2d(i,j)  = vst2d(i,j)*vmask(i,j)
!#   endif
          brk2de(i,j) =brk2de(i,j)*vmask(i,j)
          frc2de(i,j) =frc2de(i,j)*vmask(i,j)
#  endif
        enddo
      enddo       ! <-- discard stk
!
! Combined wave-current bottom drag coefficient r_D at rho-points
! ===============================================================
!
# if !defined BBL && (defined BBL_S95 || defined BBL_F00 || defined BBL_QUAD2D)
      inv_zb=1.0/max(Zob,eps)
      do j=jstr,jend
        do i=istr,iend
!#  if defined SOLVE3D && !defined BBL_QUAD2D
          cff1 = 0.5*(u(i,j,1,nrhs)+u(i+1,j,1,nrhs))
          cff2 = 0.5*(v(i,j,1,nrhs)+v(i,j+1,1,nrhs))
!#  else
!          cff1 = 0.5*(ubar(i,j,knew)+ubar(i+1,j,knew))
!          cff2 = 0.5*(vbar(i,j,knew)+vbar(i,j+1,knew))
!#  endif
          umag = sqrt( cff1*cff1+cff2*cff2 )
#  if defined BBL_F00 || defined BBL_S95
          Uorb1 = 0.5*fr(i,j)*wh(i,j)
     &              /max(sinh(min(kD(i,j),khmax)),eps)
!#  if defined BBL_F00
!          ka_f00  =0.0125   ! apparent roughness (Ruessink et al., 2001)
!          r_D(i,j)=0.015*((ka_f00*inv_d(i,j))**c1o3)
!     &               *sqrt(1.3456*(Uorb1**2)+umag**2)
!#  else
#    ifdef SOLVE3D
          cff1=umag*(vonKar/log(1.+0.5*Hz(i,j,1)*inv_zb))**2
c**          hz0 =max(z_r(i,j,1)-z_w(i,j,0),1.1*Zob)
#    else
          cff1=umag*(vonKar/log(1.+Dstp(i,j)*inv_zb))**2
c**          hz0 =max(0.5*Dstp(i,j),1.1*Zob)
#    endif
c**          cff1=umag*(vonKar/log(hz0*inv_zb))**2
          tauc=cff1*umag
          tauw=0.695*(Uorb1**1.48)*((Zob*fr(i,j))**0.52)
          cff2=1.0 + 1.2*((tauw/max(tauw+tauc,eps))**3.2)
          r_D(i,j)=cff1*cff2
!#  endif   /* BBL_F00 */
!# elif defined BBL_QUAD2D
!          r_D(i,j)=rdrg2*umag
#  endif    /* BBL_F00 || BBL_S95 */
#  ifdef MASKING
          r_D(i,j)=r_D(i,j)*rmask(i,j)
#  endif
        enddo
      enddo

#  ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=jstr,jend
          r_D(istr-1,j)=r_D(istr,j)
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=jstr,jend
          r_D(iend+1,j)=r_D(iend,j)
        enddo
      endif
#  endif
#  ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=istr,iend
          r_D(i,jstr-1)=r_D(i,jstr)
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=istr,iend
          r_D(i,jend+1)=r_D(i,jend)
        enddo
      endif
#   ifndef EW_PERIODIC
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        r_D(istr-1,jstr-1)=r_D(istr,jstr)
      endif
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        r_D(istr-1,jend+1)=r_D(istr,jend)
      endif
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        r_D(iend+1,jstr-1)=r_D(iend,jstr)
      endif
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        r_D(iend+1,jend+1)=r_D(iend,jend)
      endif
#   endif
#  endif
# endif   /* BBL */

# ifdef SOLVE3D
!
! Conservative 3D wave-current interaction variables.
! ===================================================
!
! We estimate cell-averaged 3D Stokes drift velocity (ust_r) by
! presuming that roller Stokes drift has the same z-dependency as
! the primary waves.  Analytical volume-averaging operation is
! applied to 3D Stokes drift velocity such that:
!
!   u^st = 0.5*g*A^2/frq*k/dz/(1-exp(-4*k*D))*
!         ( exp(2*k*(z_up +h-D)) -exp(-2*k*(z_up +h+D))
!          -exp(2*k*(z_low+h-D)) +exp(-2*k*(z_low+h+D)) )
!
! in which sinh, cosh, and tanh functions in u^st are expanded to
! exponential function so as not to be singular when kD gets very
! large.  Here 0.5*g*A^2/frq is equivalent to wave action density.
! The similar rule is also adapted in computation of CalP, Kapsrf,
! and Akw terms (see the code below).
!
# define ust_rm wrk1
! careful wrk3 also used for breaker body force, not for WAVE_OFFLINE
# if defined WAVE_OFFLINE
# define vst_r wrk7
# define ust_r wrk3
# endif
# define inv_dz wrk2
# define inv_ex inv_f
!
! 3D depth-dependent Stokes drift velocities
! ==========================================
!
! act= 0.5 g A^2/w
! kw= k
! h vs Dstp, Dstp= h+zeta+hzeta

      do j=jstr-1,jend+1
        do i=istr-1,iend+1
          inv_ex(i,j)=1.0/max(1.0-exp(-4.*kD(i,j)),eps)
          ust_ker=act(i,j)*kw(i,j)*inv_ex(i,j)
          cff1 = exp( 2.*kw(i,j)*(z_w(i,j,0)+h(i,j)-Dstp(i,j)))
     &          -exp(-2.*kw(i,j)*(z_w(i,j,0)+h(i,j)+Dstp(i,j)))
          do k=1,N,+1              ! <-- irreversible
            inv_dz(i,j,k)=1./Hz(i,j,k)
            cff2 = exp( 2.*kw(i,j)*(z_w(i,j,k)+h(i,j)-Dstp(i,j)))
     &            -exp(-2.*kw(i,j)*(z_w(i,j,k)+h(i,j)+Dstp(i,j)))
            ust_rm(i,j,k) =ust_ker*inv_dz(i,j,k)*( cff2-cff1 )
            cff1=cff2
          enddo
c        enddo
c      enddo                 ! <-- keep ust_r, inv_dz, &inv_ex

!
# if defined WAVE_OFFLINE && !defined U2010
!
! keff
c
c      do j=jstr-1,jend+1
c        do i=istr-1,iend+1
          keff(i,j)=max(2.*pi/20./h(i,j),0.5*sqrt(ustr0(i,j)**2
     &               +vstr0(i,j)**2)/
     &          max(.001,Dstp(i,j)*sqrt(ustr2d(i,j)**2+vstr2d(i,j)**2)))
           usec=ustr0(i,j)*tanh(2*keff(i,j)*Dstp(i,j))
           vsec=vstr0(i,j)*tanh(2*keff(i,j)*Dstp(i,j))
! rescale kD/kw*keff
          inv_exef(i,j)=1./max(1.0-exp(-4.*kD(i,j)/kw(i,j)*keff(i,j)),eps)
! act--> .5.*usec/keff**2
          ust_ker=0.5*usec/keff(i,j)*inv_exef(i,j)
          vst_ker=0.5*vsec/keff(i,j)*inv_exef(i,j)

          FACo=0.
          FAC=2.
          COUNTR=0.
          USTERR=.01
          ustm0=ust_rm(i,j,N)
          ustm0w=sqrt(ustr0(i,j)**2+vstr0(i,j))
          do while (abs(1-FAC).gt.USTERR.and.
     &         abs(FAC-FACo)/FAC.gt.1e-6.and.COUNTR.lt.15.)
          ! reset cff1
            cff1 = exp( 2.*keff(i,j)*(z_w(i,j,0)+h(i,j)-Dstp(i,j)))
     &            -exp(-2.*keff(i,j)*(z_w(i,j,0)+h(i,j)+Dstp(i,j)))
            USEZ=0.
            VSEZ=0.
            COUNTR=COUNTR+1.
            if (FACo.ne.0.) then
               ust_rm(i,j,:)=ust_rm(i,j,:)/FAC
            endif
            FACo=FAC
            do k=1,N,+1              ! <-- irreversible
c          if (COUNTR.ne.1.) then
c            ust_rm(i,j,k)=ust_rm(i,j,k)/FAC
c          endif
              if (kw(i,j)*Dstp(i,j).gt.2.*pi/20.
     &            .or.abs(ustm0-ustm0w)/ustm0w.gt..1) then ! should avoid if inside loop
C  abs(ustm0-ustm0w)/ustm0w> .1 helps in conditions of mixed wind sea
C  and swell
c            inv_dz(i,j,k)=1./Hz(i,j,k)
                cff2 = exp( 2.*keff(i,j)*(z_w(i,j,k)+h(i,j)-Dstp(i,j)))
     &                -exp(-2.*keff(i,j)*(z_w(i,j,k)+h(i,j)+Dstp(i,j)))
                vst_r(i,j,k)=ust_rm(i,j,k)*wdre(i,j)*
     &                       sqrt(tanh(abs(z_r(i,j,k)*keff(i,j)/2.)))+
     &                   (1.-sqrt(tanh(abs(z_r(i,j,k)*keff(i,j)/2.))))*
     &                       vst_ker*inv_dz(i,j,k)*( cff2-cff1 )
c             vst_r(i,j,k)= vst_ker*inv_dz(i,j,k)*( cff2-cff1 ) ! test
c             vst_r(i,j,k)= ust_r(i,j,k)*wdre(i,j) ! test
                ust_r(i,j,k)=ust_rm(i,j,k)*wdrx(i,j)*
     &                       sqrt(tanh(abs(z_r(i,j,k)*keff(i,j)/2.)))+
     &                   (1.-sqrt(tanh(abs(z_r(i,j,k)*keff(i,j)/2.))))*
     &                       ust_ker*inv_dz(i,j,k)*( cff2-cff1 )
                cff1=cff2
              else
                vst_r(i,j,k)= ust_rm(i,j,k)*wdre(i,j) !
                ust_r(i,j,k)= ust_rm(i,j,k)*wdrx(i,j)  ! monochromatic in shallow water
              endif
              USEZ=USEZ+ust_r(i,j,k)*Hz(i,j,k)
              VSEZ=VSEZ+vst_r(i,j,k)*Hz(i,j,k)
            enddo
            if (sqrt(ustr2d(i,j)**2+vstr2d(i,j)**2).ne.0) then ! ensure consistency betwew ustk and the depth integraded
              FAC=sqrt(USEZ**2+VSEZ**2)/(
     &            sqrt(ustr2d(i,j)**2+vstr2d(i,j)**2)*Dstp(i,j))
            else
              FAC=1.
            endif
!             write(*,*) 'FAC=', FAC
          end do ! while
          fac2d(i,j)=FAC
          if (kw(i,j)*Dstp(i,j).lt.2.*pi/20.
     &        .or.abs(ustm0-ustm0w)/ustm0w.lt..1) then
            keff(i,j)=kw(i,j)  ! shallow water
          endif

# endif
        enddo
      enddo                 ! <-- keep ust_r, inv_dz, &inv_ex

!  depth-integrate

! Horizontal interpolation of ust_r to u- & v- points.
! *** Notice that mask should NOT be applied here.
!
! DD test times 10
      do k=1,N
        do j=jstrR,jendR
          do i=istr,iendR
!# if !defined WAVE_OFFLINE || defined U2010
!            ust(i,j,k)=0.5*( ust_rm(i-1,j,k)*wdrx(i-1,j)
!     &                        +ust_rm(i,j,k)*wdrx(i,j) )
!# else
            ust(i,j,k)=0.5*( ust_r(i-1,j,k)
     &                        +ust_r(i,j,k))/fac2d(i,j)
# ifdef MASKING
!    &                  *umask(i,j)
# endif
!# endif
          enddo
        enddo
        do j=jstr,jendR
          do i=istrR,iendR
!# if !defined WAVE_OFFLINE || defined U2010
!            vst(i,j,k)=0.5*( ust_rm(i,j-1,k)*wdre(i,j-1)
!     &                        +ust_rm(i,j,k)*wdre(i,j) )
!# else
            vst(i,j,k)=0.5*( vst_r(i,j-1,k)
     &                        +vst_r(i,j,k) )/fac2d(i,j)
# ifdef MASKING
!    &                  *vmask(i,j)
# endif
!# endif
          enddo
        enddo
      enddo

!
!
! calP & Kapsrf: surface pressure & Bernoulli head at rho-point
! =============================================================
!
# define kv wrk4
# define dkvdz wrk5
# define d2kv  wrk6
# define actp umag
      do j=jstr,jend        ! kv (k dot v) at rho-point
        do i=istr,iend      ! and its 1st & 2nd derivertives
          do k=1,N          ! at rho-point
            kv(k) =0.5*kw(i,j)*(
     &         wdrx(i,j)* ( u(i,j,k,nstp)+u(i+1,j,k,nstp) )
     &       + wdre(i,j)* ( v(i,j,k,nstp)+v(i,j+1,k,nstp) ) )
          enddo
!# if !defined WAVE_OFFLINE || defined U2010
!          kvsurf =1.5*kv(N)-0.5*kv(N-1) ! extrapolate to surface
!# else
          kvsurf = 0.5*(ust_r(i,j,N)*( u(i,j,N,nstp)+u(i+1,j,N,nstp))
     &              + vst_r(i,j,N)*(v(i,j,N,nstp)+v(i,j+1,N,nstp) ))
!# endif
          do k=1,N-1
            dkvdz(k) =2.0*(kv(k+1)-kv(k))/(Hz(i,j,k+1)+Hz(i,j,k))
          enddo
          dkvdz(0) = dkvdz(1)    !2.*dkvdz(1)-dkvdz(2) ! severe!
          dkvdz(N) = dkvdz(N-1)  !2.*dkvdz(N-1)-dkvdz(N-2)
          do k=1,N
            d2kv(k) =dkvdz(k)-dkvdz(k-1) ! d^2kv/dz^2 x Hz
          enddo
          cff3 = 0.0
          do k=1,N
            dd   = z_r(i,j,k)-z_w(i,j,N)
            cff3 = cff3 + d2kv(k)*(
     &             exp( 2.*kw(i,j)*(dd-Dstp(i,j)))
     &           + exp(-2.*kw(i,j)*(dd+Dstp(i,j))) )
          enddo
          cff1 =-2.0*exp(-2.*kD(i,j))*inv_ex(i,j)*dkvdz(N)
          cff2 = dkvdz(0)/max(tanh(2.*kD(i,j)),eps)
          cff3 = cff3*inv_ex(i,j)
!# if !defined WAVE_OFFLINE || defined U2010
!          cff4 =-2.0*kw(i,j)*kvsurf
!# else
          cff4 = -2.0*(sinh(min(khmax,keff(i,j)*Dstp(i,j))))**2/
     &                 cosh(min(2.*khmax,2.*keff(i,j)*h(i,j)))*kvsurf
!# endif
! min limits to avoid blowups
!ccccc          actp =0.125*g*(wh(i,j)**2)/max(fr(i,j),eps)
!ccccc          calP(i,j) = actp*inv_g*tanh(kD(i,j))
# if defined SUP_OFF
          actp=actf(i,j)
# else
          actp=act(i,j)
# endif
!# if !defined WAVE_OFFLINE || defined U2010
!          calP(i,j) = actp*tanh(kD(i,j))
!     &                    *( cff1+cff2+cff3+cff4 )
!# else
          calP(i,j) = (actp*tanh(kD(i,j))
     &                    *( cff1+cff2+cff3)+ cff4) ! cfff4 now computed from surface Ustk with u
!# endif
#  ifdef MASKING
     &                                 *rmask(i,j)
#  endif
          cff=0.0
          do k=1,N
            cff = cff +  d2kv(k)*
     &         ( exp( 2.*kw(i,j)*(z_w(i,j,N)-z_r(i,j,k)-Dstp(i,j)))
     &          -exp(-2.*kw(i,j)*(z_w(i,j,N)-z_r(i,j,k)+Dstp(i,j))) )
          enddo

#  undef KAPSRF_SAFE
!#  ifdef KAPSRF_SAFE
!c#   define FLT 0.98
!#   define FLT 0.5
!          Kapsrf(i,j) = ( FLT*Kapsrf(i,j)  ! time-filter to suppress noises
!     &                   +(1.0-FLT)*cff*actp*inv_ex(i,j) )
!#  else
          Kapsrf(i,j) = cff*actp*inv_ex(i,j)
!          Kapsrf(i,j) = 0. ! DevinD zeros
!#  endif
#  ifdef MASKING
     &                           *rmask(i,j)
#  endif
        enddo
      enddo      ! <-- discard kv, dkvdz, d2kv, & inv_ex

# undef kv
# undef dkvdz
# undef d2kv
# undef inv_ex
!
! kvf : vertical vortex force term (K term) at rho-point
! ======================================================
!
# define kvr wrk5
      do j=jstr,jend
        do i=istr,iend
          do k=1,N-1       ! K term at horz rho- & vert w-point
!# if defined WAVE_OFFLINE && !defined U2010
           !! kvr = ust*u_z + vst*v_z (at w-points )
            kvr(k) = 0.25*(
     &                  (ust_r(i,j,k)+ust_r(i,j,k+1))*
     &                (u(i,j,k+1,nstp)-u(i,j,k,nstp)
     &                        +u(i+1,j,k+1,nstp)-u(i+1,j,k,nstp))
     &                 + (vst_r(i,j,k)+vst_r(i,j,k+1))*
     &                (v(i,j,k+1,nstp)-v(i,j,k,nstp)
     &                        +v(i,j+1,k+1,nstp)-v(i,j+1,k,nstp)))
     &                                 /(z_r(i,j,k+1)-z_r(i,j,k))
!# else
!            kvr(k) = 0.25*(ust_rm(i,j,k)+ust_rm(i,j,k+1))*
!     &                ( wdrx(i,j)*(u(i,j,k+1,nstp)-u(i,j,k,nstp)
!     &                        +u(i+1,j,k+1,nstp)-u(i+1,j,k,nstp))
!     &                 +wdre(i,j)*(v(i,j,k+1,nstp)-v(i,j,k,nstp)
!     &                        +v(i,j+1,k+1,nstp)-v(i,j+1,k,nstp)))
!     &                                 /(z_r(i,j,k+1)-z_r(i,j,k))
!# endif
          enddo
          kvr(0)=2.*kvr(1)-kvr(2)    ! 0.D0    ! bottom & top B.Cs.
          kvr(N)=2.*kvr(N-1)-kvr(N-2)
          do k=1,N                             ! move to vert rho-point
            kvf(i,j,k)=0.5*( kvr(k)+kvr(k-1) ) ! copy into a shared array
!            kvf(i,j,k)=0.
#  ifdef MASKING
     &                             *rmask(i,j)
#  endif
          enddo
        enddo
      enddo           ! <-- discard ust_r, kvr
# undef ust_r
# undef ust_rm
# undef vst_r
# undef kvr

!
! Akw: primary wave-induced additional diffusivity for tracers at w-point
! =======================================================================
!
      cff=wave_ramp/dt
      do j=jstrR,jendR
        do i=istrR,iendR
! Akw from Lm (llmm), mean wave length
          keffkw=2*pi/max(1.,llmm(i,j))
          cff2 = 0.25*wh(i,j)/max(1.-exp(-2.*keffkw*Dstp(i,j)),eps)
          do k=0,N
            cff1 = cff2*( exp( keffkw*(z_w(i,j,k)+h(i,j)-Dstp(i,j)))
     &                   -exp(-keffkw*(z_w(i,j,k)+h(i,j)+Dstp(i,j))) )
            cff3 = cff1**2
            Akw(i,j,k) = max(0.,cff*(cff3-E_pre(i,j,k))) ! avoid negatives
!            Akw(i,j,k) = 0. ! DevinD zeros
#  ifdef MASKING
     &                               *rmask(i,j)
#  endif
            E_pre(i,j,k) = cff3   ! store this for the next time step.
          enddo
        enddo
      enddo
!
! Non conservative wave-current interaction variables.
! ====================================================
!
! We estimate cell-averaged 3D breaking acceleration (brk_r),
! breaking-enhaced vertical eddy viscosity (Akb), and bottom
! wave streaming-induced acceleration (frc_r).  brk_r, Akb, and
! frc_r are assumed to follow an arbitrary vertical shape function,
! fb (or fkv, fwd for Akb & frc_r) with a surface/bottom scale,
! where kb^-1 = a_brk Hrms (for fb), or kb^-1 = a_kv Hrms (for fkv),
! or kb^-1 = a_frc x delta (for fwd), where Hrms is rms wave height
! and delta is turbulent wave boundary layer thickness.  a_brk,
! a_kv and a_frc are considered O(1) parameters which determine
! vertical penetration of breaking-driven momenta, additional eddy
! kinetic energy by breaking, and streaming-induced momenta near
! the bed.  fb (fkv and fwd) is designed with three (plus one
! special case for fb) different ways:
!
!  FUNC1: gb~1-tanh(kb(\zeta-z))^4 : Warner et al (2008)
!  FUNC2: gb~1-tanh(kb(\zeta-z))^2 : modified from the above
!  FUNC3: gb~cosh(kb(z+h))         : analogous to primary wave
!  FUNC0: gb~1                     : vertically uniform (Fb only)
!
! where fb(z) = gb(z) / [\int^D gb(z') dz'] = 1, thus vertical
! integral of fb is normalized.  Then volume-averaging operation is
! applied.  Exactly the same procedure is taken for fkv (not show
! here), while slightly different consideration is taken for fwd
! (see below).  In practice, vertical integrals of gb should have
! the following forms:
!
!   1) FUNC1
!
!                      4      3e^{4kbz}+3e^{2kbz}+2
!   \int gb(z') dz' = ---*-------------------------------
!                     3kb  e^{6kbz}+3e^{4kbz}+3e^{2kbz}+1
!
!                      4  3e^{-2kbz}+3e^{-4kbz}+2e^{-6kbz}
!                   = ---*---------------------------------
!                     3kb 1+3e^{-2kbz}+3e^{-4kbz}+e^{-6kbz}
!
!   2) FUNC2
!
!                      2     1         2  e^{-2kbz}
!   \int gb(z') dz' = --*---------- = --*-----------
!                     kb e^{2kbz}+1   kb 1+e^{-2kbz}
!
!   3) FUNC3
!
!   \int gb(z') dz' = [sinh(kb*z)]/kb
!
!   and then for FUNC3
!
!           kb       1     z2
!   fb = ----------*--*\int  cosh[kb(z+h)] dz
!        sinh(kb*D) dz     z1
!
!                1                                       z2
!      = -----------------*[e^{kb(z+h-D)}-e^{-kb(z+h+D)}]
!        dz*(1-e^{-2kb*D})                               z1
!
! in which hyperbolic functions (sinh, cosh, and tanh) in fb
! are expanded to exponential function so as not to be singular
! when kD is very large.
!
! Additionally, we can introduce gb=1 (vertically uniform function)
! with FB_FUNC0 for experimental purposes.  This option can be
! chosen only for fb since it may be unappropriate to fkv and fwd.
!
! Note that the above expressions are introduced in order to
! avoid the integral being Inf or NaN when kb^-1 -> 0.
!
# if defined LMD_MIXING && !defined AKB0
#  define kb0 wrk4
#  define vs  cff1
#  define kvb cff2
!
! Akb: Eddy viscosity due to depth-induced wave breaking
! ======================================================
!
!yu      do j=jstr-1,jend+1
!yu        do i=istr-1,iend+1
      do j=jstrR,jendR
        do i=istrR,iendR
          inv_k =min(max(a_kv*wh(i,j),eps),Dstp(i,j))
          kbrk  =1.0/inv_k         ! vertical scale for breaking
!#  ifdef FKV_FUNC1
!          fn1   =3.*(exp(-2.*kbrk*Dstp(i,j))+exp(-4.*kbrk*Dstp(i,j)))
!          fn2   =exp(-6.*kbrk*Dstp(i,j))
!          fb1   =c4o3*inv_k*(fn1+2.*fn2)/(1.+fn1+fn2)
!          intfb =c4o3*inv_k - fb1  ! intfb = \int^D gb(z') dz'
!          fb0   =1.D0/intfb        ! 1/[\int^D gb(z') dz']
!#  elif defined FKV_FUNC2
          fn1   =exp(-2.*kbrk*Dstp(i,j))
          fb1   =2.*inv_k*fn1/(1.+fn1)
          intfb =inv_k - fb1
          fb0   =1.D0/intfb
!#  elif defined FKV_FUNC3
!          fb1   =0.0               ! fb1 = \int gb dz at z = -h
!          fb0   =1.D0/(1.-exp(-2.*kbrk*Dstp(i,j)))
!#  endif
          vs    =ebrk(i,j)**c1o3   ! representative velocity scale
          kvb   =bconst*vs*wh(i,j) ! depth-averaged Akb
          do k=1,N,+1              ! <-- irreversible
!#  ifdef FKV_FUNC1
!            dd  =zeta(i,j,knew)-z_w(i,j,k)
!            fn1 =3.*(exp(-2.*kbrk*dd)+exp(-4.*kbrk*dd))
!            fn2 =exp(-6.*kbrk*dd)
!            fb2 =c4o3*inv_k*(fn1+2.*fn2)/(1.+fn1+fn2)
!#  elif defined FKV_FUNC2
            fn1 =exp(-2.*kbrk*(zeta(i,j,knew)-z_w(i,j,k)))
            fb2 =2.*inv_k*fn1/(1.+fn1)
!#  elif defined FKV_FUNC3
!            fb2 =exp( kbrk*(z_w(i,j,k)+h(i,j)-Dstp(i,j)))
!     &          -exp(-kbrk*(z_w(i,j,k)+h(i,j)+Dstp(i,j)))
!#  endif
            fb =fb0*inv_dz(i,j,k)*(fb2-fb1)
            kb0(k) = kvb*fb*Dstp(i,j)
            fb1=fb2     ! recursive procedure
          enddo
          do k=1,N-1
            Akb(i,j,k)=0.5*(kb0(k)+kb0(k+1))
!            Akb(i,j,k)=0. ! DevinD zeros
          enddo
          Akb(i,j,0)=0.D0
          Akb(i,j,N)=max(1.5*kb0(N)-0.5*kb0(N-1),0.0)
!          Akb(i,j,N)= 0. ! DevinD zeros
#  ifdef MASKING
          do k=0,N
            Akb(i,j,k)=Akb(i,j,k)*rmask(i,j)
          enddo
#  endif
        enddo
      enddo        ! <-- discard kb0
#  undef kb0
# endif    /* LMD_MIXING */

!# ifndef SURFACE_BREAK
!#  define brk_r wrk3
!!
!! 3D breaking acceleration term defined as a body force.
!! =====================================================
!!
!      do j=jstr-1,jend+1
!        do i=istr-1,iend+1
!!#  ifdef FB_FUNC0
!!          do k=1,N                 ! vertically uniform case
!!            brk_r(i,j,k) = brk(i,j)*inv_d(i,j)
!!          enddo                    ! <-- discard inv_d
!!#  else
!!#   ifdef FB_WSCALE
!!          inv_k =1.0/max(2.0*kw(i,j),eps)  ! Stokes scale
!!#   else
!          inv_k =min(max(a_brk*wh(i,j),eps),Dstp(i,j))
!!#   endif
!          kbrk  =1.0/inv_k         ! vertical scale for breaking
!#   ifdef FB_FUNC1
!          fn1   =3.*(exp(-2.*kbrk*Dstp(i,j))+exp(-4.*kbrk*Dstp(i,j)))
!          fn2   =exp(-6.*kbrk*Dstp(i,j))
!          fb1   =c4o3*inv_k*(fn1+2.*fn2)/(1.+fn1+fn2)
!          intfb =c4o3*inv_k - fb1  ! intfb = \int^D gb(z') dz'
!          fb0   =1.D0/intfb        ! 1/[\int^D gb(z') dz']
!#   elif defined FB_FUNC2
!          fn1   =exp(-2.*kbrk*Dstp(i,j))
!          fb1   =2.*inv_k*fn1/(1.+fn1)
!          intfb =inv_k - fb1
!          fb0   =1.D0/intfb
!#   elif defined FB_FUNC3
!          fb1   =0.0               ! fb1 = \int gb dz at z = -h
!          fb0   =1.D0/(1.-exp(-2.*kbrk*Dstp(i,j)))
!#   endif
!          do k=1,N,+1              ! <-- irreversible
!#   ifdef FB_FUNC1
!            dd  =zeta(i,j,knew)-z_w(i,j,k)
!            fn1 =3.*(exp(-2.*kbrk*dd)+exp(-4.*kbrk*dd))
!            fn2 =exp(-6.*kbrk*dd)
!            fb2 =c4o3*inv_k*(fn1+2.*fn2)/(1.+fn1+fn2)
!#   elif defined FB_FUNC2
!            fn1 =exp(-2.*kbrk*(zeta(i,j,knew)-z_w(i,j,k)))
!            fb2 =2.*inv_k*fn1/(1.+fn1)
!#   elif defined FB_FUNC3
!            fb2 =exp( kbrk*(z_w(i,j,k)+h(i,j)-Dstp(i,j)))
!     &          -exp(-kbrk*(z_w(i,j,k)+h(i,j)+Dstp(i,j)))
!#   endif
!            fb =fb0*inv_dz(i,j,k)*(fb2-fb1)
!            brk_r(i,j,k) = fb*brk(i,j)
!            fb1=fb2     ! recursive procedure
!          enddo
!#  endif       /* ifdef FB_FUNC0 */
!        enddo
!      enddo
!
!      do k=1,N
!        do j=jstrR,jendR
!          do i=istr,iendR
!            brk3dx(i,j,k)=0.5*( brk_r(i-1,j,k)*wdrx(i-1,j)
!     &                           +brk_r(i,j,k)*wdrx(i,j) )
!!            brk3dx(i,j,k)=0. ! DevinD zeros
!#  ifdef MASKING
!     &                                         *umask(i,j)
!#  endif
!          enddo
!        enddo
!        do j=jstr,jendR
!          do i=istrR,iendR
!            brk3de(i,j,k)=0.5*( brk_r(i,j-1,k)*wdre(i,j-1)
!     &                           +brk_r(i,j,k)*wdre(i,j) )
!!            brk3de(i,j,k)=0. ! DevinD
!#  ifdef MASKING
!     &                                         *vmask(i,j)
!#  endif
!          enddo
!        enddo
!      enddo       ! <-- discard brk_r, but keep ust_r
!#  undef brk_r
!# endif        /* ifndef SURFACE_BREAK */
# undef inv_d
!
! 3D bottom streaming acceleration term defined as a body force.
! ==============================================================
!
!  Convert wave-induced bottom streaming as a grid-averaged body
!  force. Similarly to the breaking acceleration, we can invoke one
!  of the following three kernal functions:
!
!   FUNC1: fb~1-tanh(kf(h+z))^4 : Warner et al (2008)
!   FUNC2: fb~1-tanh(kf(h+z))^2 : modified from the above
!   FUNC3: fb~cosh(kf(zeta-z))  : analogous to primary wave
!
!  Notice that for now the z-coordinate is inverted, and thus defined
!  from the bottom.  The length scale, k_wd^-1, is referenced to the
!  turbulent wave bottom boundary layer thickness, delta, which is
!  evaluated by Eq (2.45) in Fredsoe & Deigaard (Mechanics of Caostal
!  Sediment Transport, 1992, p.25).  There is also another option to
!  identify the laminar bottom streaming by Longuet-Higgins (1953),
!  which is required to reproduce Xu & Bowen (1994) and Lentz et al.
!  (2008).
!
!# if defined WAVE_FRICTION && defined BODY_FRICTION
!#  undef LH53_STREAMING
!#  ifdef ANA_VMIX
!#   define LH53_STREAMING
!#  endif
!#  define frc_r wrk3
!#  define frc2d stk
!      do j=jstr-1,jend+1
!        do i=istr-1,iend+1
!#  ifndef LH53_STREAMING
!          abot  =wh(i,j)/max(2.D0*sinh(min(kD(i,j),khmax)),eps)
!          delta =0.160455*((abot/Zob)**0.83)*Zob  ! BBL thickness
!          inv_k =min(max(a_frc*delta,eps),Dstp(i,j))
!          kfrc  =1.0/inv_k         ! vertical scale for wave friction
!#   ifdef FF_FUNC1
!          fn1   =3.*(exp(-2.*kfrc*Dstp(i,j))+exp(-4.*kfrc*Dstp(i,j)))
!          fn2   =exp(-6.*kfrc*Dstp(i,j))
!          fb1   =c4o3*inv_k        ! fb1: at z=-h (z_w at k=0)
!          fb2   =c4o3*inv_k*(fn1+2.*fn2)/(1.+fn1+fn2)
!          intfb =fb2-fb1           ! intfb = \int^D gb(z') dz'
!          fb0   =1.D0/intfb        ! 1/[\int^D gb(z') dz']
!#   elif defined FF_FUNC2
!          fn1   =exp(-2.*kfrc*Dstp(i,j))
!          fb1   =inv_k
!          fb2   =2.*inv_k*fn1/(1.+fn1)
!          intfb =fb2-fb1
!          fb0   =1.D0/intfb
!#   elif defined FF_FUNC3
!          fb1   =1.0-exp(-2.*kfrc*Dstp(i,j))
!          fb0   =1.D0/(exp(-2.*kfrc*Dstp(i,j))-1.0)
!#   endif
!          do k=1,N,+1              ! <-- irreversible
!#   ifdef FF_FUNC1
!            dd  =h(i,j)+z_w(i,j,k)
!            fn1 =3.*(exp(-2.*kfrc*dd)+exp(-4.*kfrc*dd))
!            fn2 =exp(-6.*kfrc*dd)
!            fb2 =c4o3*inv_k*(fn1+2.*fn2)/(1.+fn1+fn2)
!#   elif defined FF_FUNC2
!            fn1 =exp(-2.*kfrc*(h(i,j)+z_w(i,j,k)))
!            fb2 =2.*inv_k*fn1/(1.+fn1)
!#   elif defined FF_FUNC3
!            fb2 =exp( kfrc*(zeta(i,j,knew)-z_w(i,j,k)-Dstp(i,j)))
!     &          -exp(-kfrc*(zeta(i,j,knew)-z_w(i,j,k)+Dstp(i,j)))
!#   endif
!            fb =fb0*inv_dz(i,j,k)*(fb2-fb1)
!            frc_r(i,j,k) = fb*frc(i,j)
!            fb1=fb2
!          enddo
!#  else  /* LH53_STREAMING */
!          frc2d(i,j)=0.0
!          beta=sqrt(fr(i,j)/(Akv(i,j,0)+Akv(i,j,1)))
!          cff =wramp2*(wh(i,j)**2)*(fr(i,j)**2)*kw(i,j)
!     &         /max(16.0*(sinh(min(kD(i,j),khmax))**2)*beta,eps)
!          bz1 =beta*(h(i,j)+z_w(i,j,0))
!          cff1=2.*exp(-bz1)*(bz1*sin(bz1)+cos(bz1))-1.-exp(-2.*bz1)
!          cff3=1.D0 !min(1.0, pi/beta/Hz(i,j,1)) ! <-- unresoluved viscous layer
!          do k=1,N,+1              ! <-- irreversible
!            bz2 =beta*(h(i,j)+z_w(i,j,k))
!            cff2=2.*exp(-bz2)*(bz2*sin(bz2)+cos(bz2))-1.-exp(-2.*bz2)
!            frc_r(i,j,k)=-cff*(cff2-cff1)*inv_dz(i,j,k)*cff3
!            cff1=cff2; cff3=1.0
!            frc2d(i,j) = frc2d(i,j) + frc_r(i,j,k)*Hz(i,j,k)
!          enddo
!#  endif  /* !LH53_STREAMING */
!        enddo
!      enddo
!      do j=jstrR,jendR
!        do i=istr,iendR
!#  ifdef LH53_STREAMING
!          frc2dx(i,j)=0.5*( frc2d(i-1,j)*wdrx(i-1,j)
!     &                       +frc2d(i,j)*wdrx(i,j) )
!#   ifdef MASKING
!     &                                   *umask(i,j)
!#   endif
!#  endif
!          do k=1,N
!            frc3dx(i,j,k)=0.5*( frc_r(i-1,j,k)*wdrx(i-1,j)
!     &                           +frc_r(i,j,k)*wdrx(i,j) )
!#  ifdef MASKING
!     &                                         *umask(i,j)
!#  endif
!          enddo
!        enddo
!      enddo
!      do j=jstr,jendR
!        do i=istrR,iendR
!#  ifdef LH53_STREAMING
!          frc2de(i,j)=0.5*( frc2d(i,j-1)*wdre(i,j-1)
!     &                       +frc2d(i,j)*wdre(i,j) )
!#   ifdef MASKING
!     &                                   *vmask(i,j)
!#   endif
!#  endif
!          do k=1,N
!            frc3de(i,j,k)=0.5*( frc_r(i,j-1,k)*wdre(i,j-1)
!     &                           +frc_r(i,j,k)*wdre(i,j) )
!#  ifdef MASKING
!     &                                         *vmask(i,j)
!#  endif
!          enddo
!        enddo
!      enddo       ! <-- discard frc_r, frc2d
!#  undef frc_r
!#  undef frc2d
!# endif   /* BODY_FRICTION */
# undef inv_dz

# endif       /* SOLVE3D */

! Open boundary condition (valid only for idealized situations)
! =============================================================

!#   ifdef ANA_BRY
!      z_tide=0.0
!#    ifdef WKB_WWAVE
!      z_tide=wkb_tide
!#    endif
!
!#    if defined Z_FRC_BRY || defined M2_FRC_BRY || defined M3_FRC_BRY
!
!#     ifdef OBC_WEST
!      if (WESTERN_EDGE) then
!        do j=jstrR,jendR
!#      ifdef Z_FRC_BRY
!          zeta_west(j)=z_tide + sup(istrR,j)
!#      endif
!#      ifdef M2_FRC_BRY
!#       ifdef SHOREFACE
!          ubar_west(j)=ubar(istr,j,knew)
!#       else
!          ubar_west(j)=-ust2d(istrR,j)    ! Stokes drift (xi)
!#       endif
!          vbar_west(j)=vbar(istr,j,knew) ! Stokes drift (eta)
!#      endif
!#      if defined SOLVE3D && defined M3_FRC_BRY
!          do k=1,N
!            u_west(j,k)=u(istr,j,k,nnew)
!c            u_west(j,k)=-ust(istrR,j,k)
!            v_west(j,k)=v(istr,j,k,nnew)
!          enddo
!#      endif
!        enddo
!      endif
!#     endif  /* OBC_WEST */
!!
!#     ifdef OBC_EAST
!      if (EASTERN_EDGE) then
!        do j=jstrR,jendR
!#      ifdef Z_FRC_BRY
!          zeta_east(j)=z_tide + sup(iendR,j)
!#      endif
!#      ifdef M2_FRC_BRY
!c#      ifdef STOKESEKMAN
!c          ubar_east(j)=ubar(iend,j,knew) ! Neumann
!c#      else
!          ubar_east(j)=-ust2d(iend,j)  ! Stokes drift (xi)
!c#      endif
!          vbar_east(j)=vbar(iend,j,knew) !-vst2d(iendR,j)  ! Stokes drift (eta)
!#      endif
!#      if defined SOLVE3D && defined M3_FRC_BRY
!          do k=1,N
!            u_east(j,k)=u(iend,j,k,nnew)
!c            u_east(j,k)=-ust(iendR,j,k)
!            v_east(j,k)=v(iend,j,k,nnew)
!          enddo
!#      endif
!        enddo
!      endif
!#     endif  /* OBC_EAST */
!!
!#     ifdef OBC_SOUTH
!      if (SOUTHERN_EDGE) then
!        do i=istrR,iendR
!#      ifdef Z_FRC_BRY
!          zeta_south(i)=z_tide + sup(i,jstrR)
!#      endif
!#      ifdef M2_FRC_BRY
!          ubar_south(i)=ubar(i,jstr,knew)  ! Neumann
!          vbar_south(i)=-vst2d(i,jstr)     ! Stokes drift
!#      endif
!#      if defined SOLVE3D && defined M3_FRC_BRY
!          do k=1,N
!            u_south(i,k)=u(i,jstr,k,nnew)
!            v_south(i,k)=v(i,jstr,k,nnew)
!          enddo
!#      endif
!        enddo
!      endif
!#     endif  /* OBC_SOUTH */
!!
!#     ifdef OBC_NORTH
!      if (NORTHERN_EDGE) then
!        do i=istrR,iendR
!#      ifdef Z_FRC_BRY
!          zeta_north(i)=z_tide + sup(i,jendR)
!#      endif
!#      ifdef M2_FRC_BRY
!          ubar_north(i)=ubar(i,jend,knew)  ! Neumann
!          vbar_north(i)=-vst2d(i,jend)     ! Stokes drift
!#      endif
!#      if defined SOLVE3D && defined M3_FRC_BRY
!          do k=1,N
!            u_north(i,k)=u(i,jend,k,nnew)
!            v_north(i,k)=v(i,jend,k,nnew)
!          enddo
!#      endif
!        enddo
!      endif
!#     endif  /* OBC_NORTH */
!!
!#    endif  /* Z_FRC_BRY || M2_FRC_BRY || M3_FRC_BRY */
!!
!#   endif  /* ANA_BRY */


# ifdef EXCHANGE
      call exchange2d_4_tile(istr,iend,jstr,jend, sup,wdsp,wdrg,wbst)
      call exchange2d_3_tile(istr,iend,jstr,jend, ust2d,brk2dx,frc2dx)
      call exchange2d_3_tile(istr,iend,jstr,jend, vst2d,brk2de,frc2de)
!#  ifdef SURFACE_ROLLER
!      call exchange2d_tile(istr,iend,jstr,jend,rdsp)
!#  endif
#  if defined BBL_S95 || defined BBL_F00 || defined BBL_QUAD2D
      call exchange2d_tile (istr,iend,jstr,jend, r_D)
#  endif
#  ifdef SOLVE3D
      call exchange2d_2_tile(istr,iend,jstr,jend,Kapsrf,calP)
      call exchange_3_tile(istr,iend,jstr,jend, ust,N, vst,N, kvf,N)
#   ifndef SURFACE_BREAK
      call exchange_2_tile (istr,iend,jstr,jend, brk3dx,N, brk3de,N)
#   endif
!#   if defined WAVE_FRICTION && defined BODY_FRICTION
!      call exchange_2_tile (istr,iend,jstr,jend, frc3dx,N, frc3de,N)
!#   endif
#  endif
# endif

      if (wrt_file_his.or.wrt_file_avg) call wrt_wec_data

      end subroutine wec_forces_tile  !]

! ---------------------------------------------------------------------

#endif /* WEC for whole module */

      end module wec_frc
