      module wec_frc

      ! ===================================
      ! Wave Effect on Current (WEC) module
      ! ===================================
      ! which provides all WEC terms that
      ! enter in the momentum equations.
      !
      ! coding: Devin Dollery & Jeroen Molemaker
      !
      ! code essentially a copy of work already done under MRL_WCI cppflag
      ! of previous roms implementation, based on Delphine Hypolite's
      ! version of old (pre 2020) ROMS used for USWC simulations.

      use read_write ! new roms netcdf functionality module

      implicit none

#include "cppdefs.h"
! CPP flag WEC covering entire module so that the module is empty if
! WEC is unused, and thus module will not take up unnecessary size in
! the executable.
#if defined WEC /* WEC for whole module! */


      ! ****************************
      ! MODULE SPECIFICATION SECTION
      !
      ! This part of the module (before the 'contains' section)
      ! is used to declare all variables that might be needed globally
      ! or if there are memory benefits to declaring only once, rather
      ! than for each call of a subroutine.

      ! ****************************

      ! Make all variable private to this module unless public specified
      private

      ! save -Make all public variables global
      ! (might not be needed depending on compiler)
      save

      ! ==================================================
      ! INCLUDE FILES FOR COMMON VARIABLES AS PER OLD ROMS
      ! ==================================================

#include "param.h"
#include "grid.h"

! For netcdf variables
#include "ncvars.h"

! get universal scalars, g and pi
# include "scalars.h"
! get non-wec ocean variables, e.g. u & v
# include "ocean2d.h"
# ifdef SOLVE3D
#   include "ocean3d.h"
# endif


      ! =========================
      ! DECLARE ALL WEC VARIABLES
      ! =========================

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
      real wfrq(GLOBAL_2D_ARRAY)
#endif

#ifdef BBL
      real uorb(GLOBAL_2D_ARRAY)
      real vorb(GLOBAL_2D_ARRAY)
#endif   /* BBL */


      real whrm(GLOBAL_2D_ARRAY)
      real wdsp(GLOBAL_2D_ARRAY)
      real wdrg(GLOBAL_2D_ARRAY)
      real wbst(GLOBAL_2D_ARRAY)
      real wdrx(GLOBAL_2D_ARRAY)
      real wdre(GLOBAL_2D_ARRAY)
# ifdef SURFACE_ROLLER
      real rdsp(GLOBAL_2D_ARRAY)
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
!  3D  |  Kappa    |  3D Bernoulli head term (rho-point)
!  3D  |  kvf      |  vertical vortex force term (K term, 3D, rho-point)
!  3D  |  Akb      |  breaking-wave-induced additional diffusivity (w-point)
!  3D  |  Akw      |  wave-induced additional diffusivity (rho-point)
!  3D  |  E_pre    |  previous time-step value for Akw estimation (rho)
!  3D  |  frc3dx   |   xi-direction 3D frictional dissipation (rho)
!  3D  |  frc3de   |  eta-direction 3D frictional dissipation (rho)
!--------------------------------------------------------------------

      real, public :: brk2dx(GLOBAL_2D_ARRAY)
      real, public :: brk2de(GLOBAL_2D_ARRAY)
      real, public :: ust2d(GLOBAL_2D_ARRAY)
      real, public :: vst2d(GLOBAL_2D_ARRAY)
      real         :: ustr2d(GLOBAL_2D_ARRAY)
      real         :: vstr2d(GLOBAL_2D_ARRAY)
      real         :: ustr0(GLOBAL_2D_ARRAY)
      real         :: vstr0(GLOBAL_2D_ARRAY)
      real         :: llmm(GLOBAL_2D_ARRAY)
      real, public :: frc2dx(GLOBAL_2D_ARRAY) ! used in prestep, step2d, step3d
      real, public :: frc2de(GLOBAL_2D_ARRAY) ! used in prestep, step2d, step3d
      real, public :: sup(GLOBAL_2D_ARRAY) ! used in step2d
# ifdef SOLVE3D
      real, public :: calP(GLOBAL_2D_ARRAY) ! used in prsgrd32AC1.F
      real, public :: Kapsrf(GLOBAL_2D_ARRAY) ! used in prsgrd32AC1.F
#  ifndef SURFACE_BREAK
      real         :: brk3dx(GLOBAL_2D_ARRAY,N)
      real         :: brk3de(GLOBAL_2D_ARRAY,N)
#  endif
#  ifdef BODY_FRICTION
      real         :: frc3dx(GLOBAL_2D_ARRAY,N)
      real         :: frc3de(GLOBAL_2D_ARRAY,N)
#  endif
      real, public :: ust(GLOBAL_2D_ARRAY,N)
      real, public :: vst(GLOBAL_2D_ARRAY,N)
      real, public :: wst(GLOBAL_2D_ARRAY,N) ! used in set_avg.F
      real, public :: kvf(GLOBAL_2D_ARRAY,N) ! used in prsgrd32AC1.F
      real, public :: Akb(GLOBAL_2D_ARRAY,0:N) ! used in lmd_kpp.F
      real, public :: Akw(GLOBAL_2D_ARRAY,0:N) ! used in lmd_kpp.F
      real         :: E_pre(GLOBAL_2D_ARRAY,0:N)
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

      real         :: wwag(GLOBAL_2D_ARRAY,2)
      real         :: wwdg(GLOBAL_2D_ARRAY,2)
      real         :: wwpg(GLOBAL_2D_ARRAY,2)
#  ifdef BBL
      real         :: wwub(GLOBAL_2D_ARRAY,2)
      real         :: wwuob(GLOBAL_2D_ARRAY,2)
      real         :: wwvob(GLOBAL_2D_ARRAY,2)
#  endif
      real, public :: wwdrx(GLOBAL_2D_ARRAY,2)
      real, public :: wwdre(GLOBAL_2D_ARRAY,2)
      real, public :: wweb(GLOBAL_2D_ARRAY,2)
      real, public :: wveb(GLOBAL_2D_ARRAY)
      real, public :: wved(GLOBAL_2D_ARRAY)
      real, public :: wwed(GLOBAL_2D_ARRAY,2)
      real, public :: wwUSg(GLOBAL_2D_ARRAY,2)
      real, public :: wwVSg(GLOBAL_2D_ARRAY,2)
      real, public :: wwus0g(GLOBAL_2D_ARRAY,2)
      real, public :: wwvs0g(GLOBAL_2D_ARRAY,2)
#   if defined WAVE_OFFLINE && defined SURFACE_ROLLER
      real, public :: wvqb(GLOBAL_2D_ARRAY)
      real, public :: wwqb(GLOBAL_2D_ARRAY,2)
#   endif
#   if defined WAVE_OFFLINE && defined SUP_OFF
      real, public :: wwsup(GLOBAL_2D_ARRAY,2)
#  endif
#   if defined WAVE_OFFLINE && defined WAVE_DIFF_FROM_LM
      real, public :: wwlm(GLOBAL_2D_ARRAY,2)
#   endif

      real, public :: ww_cycle, wwv_time(2)
      integer, public :: ww_ncycle,  ww_rec,  itww, ntww,
     &        ww_file_id, ww_tid,  wwaid, wwpid, wwdid
#   if defined BBL && defined WAVE_OFFLINE
     &       ,wwuid, wwvid
#   endif
#   if defined WEC && defined WAVE_OFFLINE
     &       ,wweid, wwqid, wwfid, wwusaid, wwvsaid
     &       ,wwus0id, wwvs0id
#    if defined SUP_OFF
     &       ,wwsupid
#    endif
#    if defined WAVE_DIFF_FROM_LM
     &       ,wwlmid
#    endif
#   endif
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
      real sup_avg(GLOBAL_2D_ARRAY)
      real ust2d_avg(GLOBAL_2D_ARRAY)
      real vst2d_avg(GLOBAL_2D_ARRAY)
#  ifdef SOLVE3D
      real ust_avg(GLOBAL_2D_ARRAY,N)
      real vst_avg(GLOBAL_2D_ARRAY,N)
      real wst_avg(GLOBAL_2D_ARRAY,N)
      real akb_avg(GLOBAL_2D_ARRAY,0:N)
      real akw_avg(GLOBAL_2D_ARRAY,0:N)
      real kvf_avg(GLOBAL_2D_ARRAY,N)
      real calp_avg(GLOBAL_2D_ARRAY)
      real kaps_avg(GLOBAL_2D_ARRAY)
#  endif  /* SOLVE3D */
#endif /* AVERAGES */


      ! ========================
      ! NETCDF RELATED VARIABLES
      ! ========================

      ! Names of variables in WEC's netcdf forcing file
      character(len=5) :: wec_awv_name = 'Awave'    !! Amplitude
      character(len=5) :: wec_dwv_name = 'Dwave'    !! Direction
      character(len=5) :: wec_pwv_name = 'Pwave'    !! Period
      character(len=5) :: wec_usa_name = 'ust2d'    !! Vertical average of u_stokes
      character(len=5) :: wec_vsa_name = 'vst2d'    !! Vertical average of v_stokes
      character(len=4) :: wec_ust_name = 'ust0'     !! Surface u_stokes
      character(len=4) :: wec_vst_name = 'vst0'     !! Surface v_stokes
      character(len=2) :: wec_ebwv_name = 'eb'      !! wweb
      character(len=2) :: wec_edwv_name = 'ed'      !! wweb
      character(len=3) :: wec_lwv_name = 'lmw'
      character(len=4) :: wec_uob_name = 'uorb'
      character(len=4) :: wec_vob_name = 'vorb'
      character(len=3) :: wec_sup_name = 'sup'

      !! Wave forcing time (Also its dimension name in netcdf file!!)
      character(len=8) :: wec_tim_name = 'wwv_time'

      ! Add netcdf WEC terms previously contained in ncvars.h
      ! =====================================================

! ** WEC wave-current interaction variables ***
! < indx???: use between 81 and 100 >
! indxSUP               quasi-static sea-level response (wave set-up/down),
! indxUST2D,indxVST2D      2D depth-averaged Stokes drift velocities
! indxUST,indxVST,indxWST  3D Stokes drift velocities
! indxAkb,indxAkw       breaking eddy viscosity & primary wave eddy diffusivity
! indxKVF               vertical vortex force (u^St du/dz)
! indxCALP,indxKAPS     surface pressure correction & surface Bernoulli head
! indxUs0,indxVs0       Surface Stokes drift velocities
! indxLM                Mean wave length (for diffusivity approximation)
!

#if defined BBL || defined WEC
# if !defined ANA_WWAVE && !defined WKB_WWAVE
      integer, public, parameter :: indxWAVE =71,       indxWWA=indxWAVE
     &                    , indxWWD=indxWAVE+1, indxWWP=indxWAVE+2
#  ifdef WAVE_OFFLINE
#    ifdef BBL
      integer, public, parameter :: indxWUB=indxWAVE+3, indxWVB=indxWAVE+4
#    endif
#    ifdef WEC
      integer, public, parameter :: indxWEB=indxWAVE+5, indxWQB=indxWAVE+6
     &                    , indxWED=indxWAVE+7
#    endif
#  endif
# elif defined WKB_WWAVE
      integer, public, parameter :: indxHRM=71
     &                    , indxFRQ=indxHRM+1, indxWAC=indxHRM+2
     &                    , indxWKX=indxHRM+3, indxWKE=indxHRM+4
     &                    , indxEPB=indxHRM+5, indxEPD=indxHRM+6
     &                    , indxWAR=indxHRM+7, indxEPR=indxHRM+8
# endif
#endif

#ifdef WEC
      integer, public, parameter :: indxSUP=81
     &                    , indxUST2D=indxSUP+1, indxVST2D=indxSUP+2
# ifdef SOLVE3D
     &                    , indxUST=indxSUP+3, indxVST=indxSUP+4
     &                    , indxWST=indxSUP+5, indxAkb=indxSUP+6
     &                    , indxAkw=indxSUP+7, indxKVF=indxSUP+8
     &                    , indxCALP=indxSUP+9, indxKAPS=indxSUP+10
     &                    , indxUs0=indxSUP+11, indxVs0=indxSUP+12
     &                    , indxLM=indxSUP+13
# endif
#endif

      !!!!!!!!

#ifdef WEC
      integer, public :: hisSUP, hisUST2D, hisVST2D, hisUSTR0
# ifdef SOLVE3D
      integer, public :: hisUST, hisVST, hisAkb, hisAkw, hisKVF,
     &        hisCALP, hisKAPS, hisWST
# endif
#endif

#ifdef AVERAGES
# ifdef WEC
      integer, public :: avgSUP, avgUST2D, avgVST2D
#  ifdef SOLVE3D
      integer, public :: avgUST, avgVST, avgAkb, avgAkw, avgKVF,
     &        avgCALP, avgKAPS, avgWST
#  endif
# endif
#endif /* AVERAGES */

!      end type wec_ncvars

      ! DevinD - end wave terms previously contained in ncvars.h

      ! ===============
      ! MISC. VARIABLES
      ! ===============

      ! coefficient variables for external calculations in prsgrd32AC1.F
      real, public :: cff_ext, cff_ext2

      ! Rename grid variables
      integer, parameter :: gnx = LLm
      integer, parameter :: nz  = N

       ! timestep variables for read in WEC data
      integer              :: it1 = 1  ! step n. Note: should be parameter as unchanging
      integer              :: it2 = 2  ! step n+1. Note: should be parameter as unchanging
      real, dimension(2)   :: wec_tim = [-99,-99] ! used to store both time steps (n, n+1)
      integer    :: ifile = 0 ! read in file number
      ! current record number of data point in read in data - irec
      integer    :: irec  = 0

      ! flags taken from def_his.F
      ! ==========================
      ! and used in def_his_def_var_wec & def_his_nf_inq_varid
#ifndef AVRH
# define ncid nchis
# define fname hisname
# define wrt wrthis
# ifdef WEC
#  define vidSUP hisSUP
#  define vidUST2D hisUST2D
#  define vidVST2D hisVST2D
#  define vidUSTR0 hisUSTR0
#  define vidUST hisUST
#  define vidVST hisVST
#  define vidWST hisWST
#  define vidAkb hisAkb
#  define vidAkw hisAkw
#  define vidKVF hisKVF
#  define vidCALP hisCALP
#  define vidKAPS hisKAPS
# endif
#else
# define ncid ncavg
# define fname avgname
# define wrt wrtavg
# ifdef WEC
#  define vidSUP avgSUP
#  define vidUST2D avgUST2D
#  define vidVST2D avgVST2D
#  define vidUST avgUST
#  define vidWST avgWST
#  define vidVST avgVST
#  define vidAkb avgAkb
#  define vidAkw avgAkw
#  define vidKVF avgKVF
#  define vidCALP avgCALP
#  define vidKAPS avgKAPS
# endif
#endif


      ! ================================
      ! WEC subroutines needed as public (accessed outside of module)
      ! ================================

      public init_scalars_wec     ! called from init_scalars.F
      public read_inp_wec         ! called from read_inp.F
      public set_wec_frc_tile     ! called from set_forces.F
      public def_his_def_var_wec  ! called from def_his.F
      public def_his_nf_inq_varid ! called from def_his.F
      public init_arrays_wec_tile ! called from init_arrays.F
      public compute_wec_rhs_uv_terms ! called from step3d_uv1.F & pre_step3d4S.F
      public wstokes_wec_tile     ! called from set_avg.F
      public wrt_his_wec          ! called from wrt_his.F
      public wec_forces           ! called from main.F
#ifdef AVERAGES
      public read_inp_avg_wec     ! called from read_inp.F
      public set_avg_wec_tile     ! called from set_avg.F
      public wrt_avg_wec          ! called from wrt_avg.F
#endif

      ! =======================
      ! END OF MODULE PRE-AMBLE
      ! =======================

      contains 

      ! --------------------------------------------------------------------
      ! --------------------------------------------------------------------
      ! --------------------------------------------------------------------


      subroutine init_scalars_wec

      ! Called from init_scalars.F for WEC specific scalars

      implicit none

      vname(1,indxSUP)='sup'
      vname(2,indxSUP)='quasi-static sea-level response (set-up)'
      vname(3,indxSUP)='meter'

      vname(1,indxUST2D)='ust2d'
      vname(2,indxUST2D)='u depth-averaged Stokes drift velocity'
      vname(3,indxUST2D)='meter second-1'

      vname(1,indxVST2D)='vst2d'
      vname(2,indxVST2D)='v depth-averaged Stokes drift velocity'
      vname(3,indxVST2D)='meter second-1'
# ifdef SOLVE3D
      vname(1,indxUST)  ='ust'
      vname(2,indxUST)  ='u 3D Stokes drift velocity'
      vname(3,indxUST)  ='meter second-1'

      vname(1,indxVST)  ='vst'
      vname(2,indxVST)  ='v 3D Stokes drift velocity'
      vname(3,indxVST)  ='meter second-1'

      vname(1,indxWST)  ='wst'
      vname(2,indxWST)  ='w veitical 3D Stokes drift velocity'
      vname(3,indxWST)  ='meter second-1'

      vname(1,indxAKB)  ='Akb'
      vname(2,indxAKB)  ='eddy viscosity due to wave breaking'
      vname(3,indxAKB)  ='meter2 second-1'

      vname(1,indxAKW)  ='Akw'
      vname(2,indxAKW)  ='eddy diffusivity due to primary waves'
      vname(3,indxAKW)  ='meter2 second-1'

      vname(1,indxKVF)  ='kvf'
      vname(2,indxKVF)  ='vertical vortex force (u^st du/dz)'
      vname(3,indxKVF)  ='meter second-2'

      vname(1,indxCALP) ='calP'
      vname(2,indxCALP) ='surface pressure correction in prsgrd'
      vname(3,indxCALP) ='meter2 second-2'

      vname(1,indxKAPS) ='Kapsrf'
      vname(2,indxKAPS) ='surface Bernoulli head in prsgrd term'
      vname(3,indxKAPS) ='meter2 second-2'

      vname(1,indxUs0)  ='ust0'
      vname(2,indxUs0)  ='Surface Stokes drift velocity'
      vname(3,indxUs0)  ='meter second-1'

      vname(1,indxVs0)  ='vst0'
      vname(2,indxVs0)  ='Surface Stokes drift velocity'
      vname(3,indxVs0)  ='meter second-1'

      vname(1,indxLM)   ='lmw'
      vname(2,indxLM)   ='Mean Wave length'
      vname(3,indxLM)   ='meter'
# endif

#if defined BBL || defined WEC
# if !defined ANA_WWAVE && !defined WKB_WWAVE
      vname(1,indxWWA)='Awave'
      vname(2,indxWWA)='wind induced wave amplitude'
      vname(3,indxWWA)='meter'

      vname(1,indxWWD)='Dwave'
      vname(2,indxWWD)='wind induced wave direction'
      vname(3,indxWWD)='degrees'

      vname(1,indxWWP)='Pwave'
      vname(2,indxWWP)='wind induced wave Period'
      vname(3,indxWWP)='second'
#  ifdef WAVE_OFFLINE
#    ifdef BBL
      vname(1,indxWUB)='uorb'
      vname(2,indxWUB)='Eastward bottom orbital velocity (SWAN)'
      vname(3,indxWUB)='meter second-1'
      vname(1,indxWVB)='vorb'
      vname(2,indxWVB)='Northward bottom orbital velocity (SWAN)'
      vname(3,indxWVB)='meter second-1'
#    endif
#    ifdef WEC
      vname(1,indxWEB)='eb'
      vname(2,indxWEB)='breaking energy dissipation (SWAN)'
      vname(3,indxWEB)='meter3 second-3'

#      ifdef WAVE_FRICTION
      vname(1,indxWED)='ed'
      vname(2,indxWED)='bottom frictional dissipation (SWAN)'
      vname(3,indxWED)='meter3 second-3'
#      endif

#      ifdef SURFACE_ROLLER
      vname(1,indxWQB)='qb'
      vname(2,indxWQB)='fraction of breaking waves (SWAN)'
      vname(3,indxWQB)='non dimension'
#      endif
#    endif
#  endif
# endif
#endif  /* BBL || WEC */

      end subroutine init_scalars_wec

      ! --------------------------------------------------------------------

      subroutine init_arrays_wec_tile (istr,iend,jstr,jend)

      ! This is called from init_arrays_tile of init_arrays.F
      ! It is needed to "first touch" arrays to optimize
      ! NUMA shared memory allocation on linux.

      implicit none

      integer, intent(in) :: istr, iend, jstr, jend
      ! Declare looping indices
      integer :: i, j, k

      real, parameter :: init=0.

      ! compute_auxiliary_bounds.h is needed for istrR, jstrR, etc
      ! This had to declared last else won't compile
#include "compute_auxiliary_bounds.h"

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


      end subroutine init_arrays_wec_tile


      ! --------------------------------------------------------------------


      subroutine read_inp_wec (keyword, ierr, kwlen)

      implicit none

      character(len=32) :: keyword
      integer ierr, kwlen
      integer, parameter :: input=15

        call cancel_kwd (keyword(1:kwlen), ierr)

!# ifndef AVERAGES
        read(input,*,err=95) wrthis(indxSUP),  wrthis(indxUST2D)
     &                                        ,  wrthis(indxVST2D)
#  ifdef SOLVE3D
     &         , wrthis(indxUST), wrthis(indxVST), wrthis(indxwST)
     &         , wrthis(indxAKB), wrthis(indxAKW), wrthis(indxKVF)
     &                        , wrthis(indxCALP), wrthis(indxKAPS)
#  endif
          mpi_master_only write(*,'(/1x,A,3(/8x,A,T16,L1,T20,A))')
     &               'fields to be saved in MRL-WEC 2D history: (T/F)'
     &             ,  'SUP',    wrthis(indxSUP),    vname(2,indxSUP)
     &             ,  'UST2D',  wrthis(indxUST2D),  vname(2,indxUST2D)
     &             ,  'VST2D',  wrthis(indxVST2D),  vname(2,indxVST2D)

#  ifdef SOLVE3D
          mpi_master_only write(*,'(/1x,A,8(/8x,A,T16,L1,T20,A))')
     &               'fields to be saved in MRL-WEC 3D history: (T/F)'
     &             ,  'UST',    wrthis(indxUST),    vname(2,indxUST)
     &             ,  'VST',    wrthis(indxVST),    vname(2,indxVST)
     &             ,  'WST',    wrthis(indxWST),    vname(2,indxWST)
     &             ,  'Akb',    wrthis(indxAkb),    vname(2,indxAkb)
     &             ,  'Akw',    wrthis(indxAkw),    vname(2,indxAkw)
     &             ,  'KVF',    wrthis(indxKVF),    vname(2,indxKVF)
     &             ,  'CALP',   wrthis(indxCALP),   vname(2,indxCALP)
     &             ,  'KAPS',   wrthis(indxKAPS),   vname(2,indxKAPS)
#  endif


      goto 100
  95  write(*,'(/1x,4A/)') '### ERROR: read_inp_wec :: Cannot read ',
     &                       'entry ''', keyword(1:kwlen), '''.'
      ierr=ierr+1 ! Possible source for error here?
!      goto 99
!  99  close (input)
 100  continue

      end subroutine read_inp_wec

      ! ---------------------------------------------------
      ! ---------------------------------------------------
      ! ---------------------------------------------------

#ifdef AVERAGES
      subroutine read_inp_avg_wec (keyword, ierr, kwlen)

      implicit none

      character(len=32) :: keyword
      integer ierr, kwlen
      integer, parameter :: input=15

        call cancel_kwd (keyword(1:kwlen), ierr)

!# else
!        elseif (keyword(1:kwlen).eq.'wci_average_fields') then
!          call cancel_kwd (keyword(1:kwlen), ierr)
        read(input,*,err=95) wrtavg(indxSUP),  wrtavg(indxUST2D)
     &                                        ,  wrtavg(indxVST2D)
#  ifdef SOLVE3D
     &         , wrtavg(indxUST), wrtavg(indxVST), wrtavg(indxwST)
     &         , wrtavg(indxAKB), wrtavg(indxAKW), wrtavg(indxKVF)
     &                        , wrtavg(indxCALP), wrtavg(indxKAPS)
#  endif
          mpi_master_only write(*,'(/1x,A,3(/8x,A,T16,L1,T20,A))')
     &               'fields to be saved in MRL-WEC 2D average: (T/F)'
     &             ,  'SUP',    wrtavg(indxSUP),    vname(2,indxSUP)
     &             ,  'UST2D',  wrtavg(indxUST2D),  vname(2,indxUST2D)
     &             ,  'VST2D',  wrtavg(indxVST2D),  vname(2,indxVST2D)

#  ifdef SOLVE3D
          mpi_master_only write(*,'(/1x,A,8(/8x,A,T16,L1,T20,A))')
     &               'fields to be saved in MRL-WEC 3D average: (T/F)'
     &             ,  'UST',    wrtavg(indxUST),    vname(2,indxUST)
     &             ,  'VST',    wrtavg(indxVST),    vname(2,indxVST)
     &             ,  'WST',    wrtavg(indxWST),    vname(2,indxWST)
     &             ,  'Akb',    wrtavg(indxAkb),    vname(2,indxAkb)
     &             ,  'Akw',    wrtavg(indxAkw),    vname(2,indxAkw)
     &             ,  'KVF',    wrtavg(indxKVF),    vname(2,indxKVF)
     &             ,  'CALP',   wrtavg(indxCALP),   vname(2,indxCALP)
     &             ,  'KAPS',   wrtavg(indxKAPS),   vname(2,indxKAPS)
#  endif
!# endif /* AVERAGES */

      goto 100
  95  write(*,'(/1x,4A/)') '### ERROR: read_inp_wec :: Cannot read ',
     &                       'entry ''', keyword(1:kwlen), '''.'
      ierr=ierr+1 ! Possible source for error here?
!      goto 99
!  99  close (input)
 100  continue

      end subroutine read_inp_avg_wec
#endif /* AVERAGES */

      ! ---------------------------------------------------
      ! ---------------------------------------------------
      ! ---------------------------------------------------



      subroutine def_his_def_var_wec (ierr, r2dgrd, u2dgrd, v2dgrd
#ifdef SOLVE3D
     &                              , r3dgrd, u3dgrd, v3dgrd, w3dgrd
#endif
     &                                                             )

!      use netcdf, only: nf_float ! Didn't work

      implicit none



      integer lenstr, lvar, ierr, r2dgrd(3), u2dgrd(3), v2dgrd(3)
#ifdef SOLVE3D
     &      , r3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
#endif
      integer, external :: nf_def_var, nf_put_att_text, nf_put_att_real ! Taken from netcdf.inc
      character(len=64) text
      integer, parameter :: nf_float = 5

! Taken from def_his.F
#ifdef MASK_LAND_DATA
# ifdef HIS_DOUBLE
      real*8, parameter :: spv_set=1.D+33
# else
      real*4, parameter :: spv_set=1.E+33
# endif
#endif

!! ATTEMPTED POINTER IMPLEMENTATION
!
!      type domainptr
!        integer, pointer :: p => null() ! For target hisUST2D above
!      end type domainptr
!      ! These variables are new and needed for looping of attributes
!      type(domainptr), dimension(3) :: dom ! array of pointers to hisUST2d etc
!      integer indx_loop ! looping index
!      integer loop_grd(3,3) ! grid type - can make second dimension (i,:) for length 3 and 4
!      integer :: i=0 ! additional looping integer
!      ! 3D loop variables - may need to deallocate pointers at end of subroutine if memory leak!
!      type(domainptr), dimension(8) :: dom3D ! array of pointers to hisUST2d etc
!      integer loop_grd3D(8,4) ! grid type - can make second dimension (i,:) for length 3 and 4
!      integer loop_dims(8) ! array of variable's dimensions
!      ! Put in type of variable here.
!      loop_grd(1,:) = r2dgrd; loop_grd(2,:) = u2dgrd; loop_grd(3,:) = v2dgrd
!      ! Put in pointer to varid here (e.g. hisUST2D - see flags above)
!!      allocate(dom(1)%p) ; allocate(dom(2)%p) ; allocate(dom(3)%p)
!      dom(1)%p => vidSUP; dom(2)%p => vidUST2D; dom(3)%p => vidVST2D
!      write(*,*) 'vidSUP = ', vidSUP
!      do indx_loop = indxSUP, indxVST2D ! Loop through the relevant indices for wrt
!        i=i+1 ! This is for loop index from 1
!        if (wrt(indx_loop)) then
!           write(*,*) 'wrt(indx_loop)'
!          lvar=lenstr(vname(1,indx_loop))
!          ierr=nf_def_var (ncid, vname(1,indx_loop)(1:lvar), NF_FOUT,
!     &                                      3, loop_grd(i,:), dom(i)%p)
!# ifdef AVRH
!          text='averaged '/ /vname(2,indx_loop)
!# else
!          text=vname(2,indx_loop)
!# endif
!          lvar=lenstr(text)
!          ierr=nf_put_att_text (ncid, dom(i)%p, 'long_name', lvar,
!     &                                               text(1:lvar))
!          lvar=lenstr(vname(3,indx_loop))
!          ierr=nf_put_att_text (ncid, dom(i)%p, 'units',     lvar,
!     &                                 vname(3,indx_loop)(1:lvar))
!#  ifdef MASK_LAND_DATA
!          ierr=nf_put_att_FOUT (ncid, dom(i)%p, '_FillValue',
!     &                                NF_FOUT, 1, spv_set)
!#  endif
!          nullify(dom(i)%p) ! disassociate pointer from target
!!          deallocate(dom(i)%p)
!        endif
!      end do

!!!  END ATTEMPTED POINTER IMPLEMENTATION

! THIS IS OLD CODE REPLACED BY LOOP ABOVE - LEFT IN UNTIL SURE NO MEMORY LEAK
! WHEN DOING LARGE SCALE NON-ANALYTICAL MODEL

! Wave-current interaction diagnostic variables.
! =============================================
!
! sup: quasi-static sea-level response (wave set-up/down)
!
        if (wrt(indxSUP)) then
          lvar=lenstr(vname(1,indxSUP))
          ierr=nf_def_var (ncid, vname(1,indxSUP)(1:lvar), NF_FOUT,
     &                                           3, r2dgrd, vidSUP)
# ifdef AVRH
          text='averaged '/ /vname(2,indxSUP)
# else
          text=vname(2,indxSUP)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidSUP, 'long_name', lvar,
     &                                             text(1:lvar))
          lvar=lenstr(vname(3,indxSUP))
          ierr=nf_put_att_text (ncid, vidSUP, 'units',     lvar,
     &                                 vname(3,indxSUP)(1:lvar))
#  ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidSUP, '_FillValue',
     &                                NF_FOUT, 1, spv_set)
#  endif
        endif
!
! ust2d & vst2d: 2D, depth-averaged Stokes drift velocities.
!
        if (wrt(indxUST2D)) then
          lvar=lenstr(vname(1,indxUST2D))
          ierr=nf_def_var (ncid, vname(1,indxUST2D)(1:lvar), NF_FOUT,
     &                                           3, u2dgrd, vidUST2D)
# ifdef AVRH
          text='averaged '/ /vname(2,indxUST2D)
# else
          text=vname(2,indxUST2D)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidUST2D, 'long_name', lvar,
     &                                               text(1:lvar))
          lvar=lenstr(vname(3,indxUST2D))
          ierr=nf_put_att_text (ncid, vidUST2D, 'units',     lvar,
     &                                 vname(3,indxUST2D)(1:lvar))
#  ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidUST2D, '_FillValue',
     &                                NF_FOUT, 1, spv_set)
#  endif

        endif
        if (wrt(indxVST2D)) then
          lvar=lenstr(vname(1,indxVST2D))
          ierr=nf_def_var (ncid, vname(1,indxVST2D)(1:lvar), NF_FOUT,
     &                                           3, v2dgrd, vidVST2D)
# ifdef AVRH
          text='averaged '/ /vname(2,indxVST2D)
# else
          text=vname(2,indxVST2D)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidVST2D, 'long_name', lvar,
     &                                               text(1:lvar))
          lvar=lenstr(vname(3,indxVST2D))
          ierr=nf_put_att_text (ncid, vidVST2D, 'units',     lvar,
     &                                  vname(3,indxVST2D)(1:lvar))
#  ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidVST2D, '_FillValue',
     &                                NF_FOUT, 1, spv_set)
#  endif
        endif
# ifdef SOLVE3D

!!! ATTEMPTED POINTER IMPLEMENTATION

!      ! Put in type of variable here.
!      loop_grd3D(1,:) = u3dgrd; loop_grd3D(2,:) = v3dgrd; loop_grd3D(3,:) = r3dgrd
!      loop_grd3D(4,:) = w3dgrd; loop_grd3D(5,:) = w3dgrd; loop_grd3D(6,:) = r3dgrd
!      loop_grd3D(7,1:3) = r2dgrd; loop_grd3D(8,1:3) = r2dgrd ! Potential bug - 3 not 4 length integers
!      ! Put in pointer to varid here (e.g. hisUST2D - see flags above)
!      dom3D(1)%p => vidUST; dom3D(2)%p => vidVST; dom3D(3)%p => vidWST
!      dom3D(4)%p => vidAkb; dom3D(5)%p => vidAkw; dom3D(6)%p => vidKVF
!      dom3D(7)%p => vidCALP; dom3D(8)%p => vidKAPS
!      ! The number of dimensions of variable (includes time)
!      loop_dims = [4,4,4,4,4,4,3,3]
!
!      ! Could avoid having whole extra loop by using if statement for loop over 3d
!      ! variables or not. I.e. add to sup, ust2d and vst2d loop above.
!      i=0 ! reset loop counter
!      do indx_loop = indxUST, indxKAPS ! Loop through relevant indices for wrt
!        i=i+1 ! This is for loop index from 1
!                if (wrt(indx_loop)) then
!          lvar=lenstr(vname(1,indx_loop))
!          ierr=nf_def_var (ncid, vname(1,indx_loop)(1:lvar), NF_FOUT,
!     &                        loop_dims(i), loop_grd3D(i,:), dom3D(i)%p)
!# ifdef AVRH
!          text='averaged '/ /vname(2,indx_loop)
!# else
!          text=vname(2,indx_loop)
!# endif
!          lvar=lenstr(text)
!          ierr=nf_put_att_text (ncid, dom3D(i)%p, 'long_name', lvar,
!     &                                             text(1:lvar))
!          lvar=lenstr(vname(3,indx_loop))
!          ierr=nf_put_att_text (ncid, dom3D(i)%p, 'units',     lvar,
!     &                                 vname(3,indx_loop)(1:lvar))
!#  ifdef MASK_LAND_DATA
!          ierr=nf_put_att_FOUT (ncid, dom3D(i)%p, '_FillValue',
!     &                                NF_FOUT, 1, spv_set)
!#  endif
!          nullify(dom3D(i)%p) ! disassociate pointer from target
!        endif
!
!      end do

!!! END ATTEMPTED POINTER IMPLEMENTATION

! THIS IS OLD CODE REPLACED BY LOOP ABOVE - LEFT IN UNTIL SURE NO MEMORY LEAK
! WHEN DOING LARGE SCALE NON-ANALYTICAL MODEL

!
! ust3d & vst3d: 3D Stokes drift velocities.
!
        if (wrt(indxUST)) then
          lvar=lenstr(vname(1,indxUST))
          ierr=nf_def_var (ncid, vname(1,indxUST)(1:lvar), NF_FOUT,
     &                                           4, u3dgrd, vidUST)
# ifdef AVRH
          text='averaged '/ /vname(2,indxUST)
# else
          text=vname(2,indxUST)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidUST, 'long_name', lvar,
     &                                             text(1:lvar))
          lvar=lenstr(vname(3,indxUST))
          ierr=nf_put_att_text (ncid, vidUST, 'units',     lvar,
     &                                 vname(3,indxUST)(1:lvar))
#  ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidUST, '_FillValue',
     &                                NF_FOUT, 1, spv_set)
#  endif
        endif
        if (wrt(indxVST)) then
          lvar=lenstr(vname(1,indxVST))
          ierr=nf_def_var (ncid, vname(1,indxVST)(1:lvar), NF_FOUT,
     &                                           4, v3dgrd, vidVST)
# ifdef AVRH
          text='averaged '/ /vname(2,indxVST)
# else
          text=vname(2,indxVST)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidVST, 'long_name', lvar,
     &                                             text(1:lvar))
          lvar=lenstr(vname(3,indxVST))
          ierr=nf_put_att_text (ncid, vidVST, 'units',     lvar,
     &                                 vname(3,indxVST)(1:lvar))
#  ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidVST, '_FillValue',
     &                                NF_FOUT, 1, spv_set)
#  endif
        endif
!
! wst: vertical Stokes drift velocity (m/s) at rho-point
!
        if (wrt(indxWST)) then
          lvar=lenstr(vname(1,indxWST))
          ierr=nf_def_var (ncid, vname(1,indxWST)(1:lvar), NF_FOUT,
     &                                           4, r3dgrd, vidWST)
# ifdef AVRH
          text='averaged '/ /vname(2,indxWST)
# else
          text=vname(2,indxWST)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidWST, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxWST))
          ierr=nf_put_att_text (ncid, vidWST, 'units',     lvar,
     &                                 vname(3,indxWST)(1:lvar))
#  ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidWST, '_FillValue',
     &                                NF_FOUT, 1, spv_set)
#  endif
        endif
!
! Akb: Vertical eddy viscosity coefficient due to depth-induced wave breaking.
!
        if (wrt(indxAkb)) then
          lvar=lenstr(vname(1,indxAkb))
          ierr=nf_def_var (ncid, vname(1,indxAkb)(1:lvar), NF_FOUT,
     &                                           4, w3dgrd, vidAkb)
# ifdef AVRH
          text='averaged '/ /vname(2,indxAkb)
# else
          text=vname(2,indxAkb)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidAkb, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxAkb))
          ierr=nf_put_att_text (ncid, vidAkb, 'units',     lvar,
     &                                  vname(3,indxAkb)(1:lvar))
#  ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidAkb, '_FillValue',
     &                                NF_FOUT, 1, spv_set)
#  endif
        endif
! Akw: Vertical eddy diffusivity coefficient due to primary waves.
!
        if (wrt(indxAkw)) then
          lvar=lenstr(vname(1,indxAkw))
          ierr=nf_def_var (ncid, vname(1,indxAkw)(1:lvar), NF_FOUT,
     &                                           4, w3dgrd, vidAkw)
# ifdef AVRH
          text='averaged '/ /vname(2,indxAkw)
# else
          text=vname(2,indxAkw)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidAkw, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxAkw))
          ierr=nf_put_att_text (ncid, vidAkw, 'units',     lvar,
     &                                  vname(3,indxAkw)(1:lvar))
#  ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidAkw, '_FillValue',
     &                                NF_FOUT, 1, spv_set)
#  endif
        endif
!
! kvf: vertical vortex force term (u^St du/dz)
!
        if (wrt(indxKVF)) then
          lvar=lenstr(vname(1,indxKVF))
          ierr=nf_def_var (ncid, vname(1,indxKVF)(1:lvar), NF_FOUT,
     &                                           4, r3dgrd, vidKVF)
# ifdef AVRH
          text='averaged '/ /vname(2,indxKVF)
# else
          text=vname(2,indxKVF)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidKVF, 'long_name', lvar,
     &                                             text(1:lvar))
          lvar=lenstr(vname(3,indxKVF))
          ierr=nf_put_att_text (ncid, vidKVF, 'units',     lvar,
     &                                 vname(3,indxKVF)(1:lvar))
#  ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidKVF, '_FillValue',
     &                                NF_FOUT, 1, spv_set)
#  endif
        endif

!
! calP: surface pressure correction term appeared in prsgrd term
!
        if (wrt(indxCALP)) then
          lvar=lenstr(vname(1,indxCALP))
          ierr=nf_def_var (ncid, vname(1,indxCALP)(1:lvar), NF_FOUT,
     &                                           3, r2dgrd, vidCALP)
# ifdef AVRH
          text='averaged '/ /vname(2,indxCALP)
# else
          text=vname(2,indxCALP)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidCALP, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxCALP))
          ierr=nf_put_att_text (ncid, vidCALP, 'units',     lvar,
     &                                 vname(3,indxCALP)(1:lvar))
#  ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidCALP, '_FillValue',
     &                                NF_FOUT, 1, spv_set)
#  endif
        endif
!
! Kapsrf: surface Bernoulli head term appeared in prsgrd term
!
        if (wrt(indxKAPS)) then
          lvar=lenstr(vname(1,indxKAPS))
          ierr=nf_def_var (ncid, vname(1,indxKAPS)(1:lvar), NF_FOUT,
     &                                           3, r2dgrd, vidKAPS)
# ifdef AVRH
          text='averaged '/ /vname(2,indxKAPS)
# else
          text=vname(2,indxKAPS)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, vidKAPS, 'long_name', lvar,
     &                                              text(1:lvar))
          lvar=lenstr(vname(3,indxKAPS))
          ierr=nf_put_att_text (ncid, vidKAPS, 'units',     lvar,
     &                                 vname(3,indxKAPS)(1:lvar))
#  ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, vidKAPS, '_FillValue',
     &                                NF_FOUT, 1, spv_set)
#  endif
        endif
# endif  /* SOLVE3D */

      end subroutine def_his_def_var_wec


      ! --------------------------------------------------------------------
      ! --------------------------------------------------------------------
      ! --------------------------------------------------------------------

      subroutine def_his_nf_inq_varid (ierr, lfnm)

      implicit none

      integer lenstr, lvar, ierr, lfnm
      ! Taken from netcdf.inc to avoid include. Used before consider top of module
      integer, parameter :: nf_noerr=0
      integer, external  :: nf_inq_varid ! unique to this method

!
! sup: quasi-static sea-level response (wave set-up/down)
!
        if (wrt(indxSUP)) then
          lvar=lenstr(vname(1,indxSUP))
          ierr=nf_inq_varid (ncid, vname(1,indxSUP)(1:lvar), vidSUP)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxSUP)(1:lvar), fname(1:lfnm)
            goto 99                                       !--> ERROR
          endif
        endif
!
! ust2d & vst2d: 2D, depth-averaged Stokes drift velocities.
!
        if (wrt(indxUST2D)) then
          lvar=lenstr(vname(1,indxUST2D))
          ierr=nf_inq_varid (ncid, vname(1,indxUST2D)(1:lvar), vidUST2D)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxUST2D)(1:lvar), fname(1:lfnm)
            goto 99                                       !--> ERROR
          endif
        endif
        if (wrt(indxVST2D)) then
          lvar=lenstr(vname(1,indxVST2D))
          ierr=nf_inq_varid (ncid, vname(1,indxVST2D)(1:lvar), vidVST2D)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxVST2D)(1:lvar), fname(1:lfnm)
            goto 99                                       !--> ERROR
          endif
        endif
# ifdef SOLVE3D
!
! ust3d & vst3d: 3D Stokes drift velocities.
!
        if (wrt(indxUST)) then
          lvar=lenstr(vname(1,indxUST))
          ierr=nf_inq_varid (ncid, vname(1,indxUST)(1:lvar), vidUST)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxUST)(1:lvar), fname(1:lfnm)
            goto 99                                       !--> ERROR
          endif
        endif

        if (wrt(indxVST)) then
          lvar=lenstr(vname(1,indxVST))
          ierr=nf_inq_varid (ncid, vname(1,indxVST)(1:lvar), vidVST)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxVST)(1:lvar), fname(1:lfnm)
            goto 99                                       !--> ERROR
          endif
        endif
!
! wst: vertical Stokes drift velocity (m/s) at rho-point
!
        if (wrthis(indxWST)) then
          lvar=lenstr(vname(1,indxWST))
          ierr=nf_inq_varid (ncid, vname(1,indxWST)(1:lvar), vidWST)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxWST)(1:lvar), fname(1:lfnm)
            goto 99                                       !--> ERROR
          endif
        endif
!
! Akb: Vertical eddy viscosity coefficient due to depth-induced wave breaking.
!
        if (wrt(indxAkb)) then
          lvar=lenstr(vname(1,indxAkb))
          ierr=nf_inq_varid (ncid, vname(1,indxAkb)(1:lvar), vidAkb)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxAkb)(1:lvar), fname(1:lfnm)
            goto 99                                       !--> ERROR
          endif
        endif
!
! Akw: Vertical eddy diffusivity coefficient due to primary waves.
!
        if (wrt(indxAkw)) then
          lvar=lenstr(vname(1,indxAkw))
          ierr=nf_inq_varid (ncid, vname(1,indxAkw)(1:lvar), vidAkw)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxAkw)(1:lvar), fname(1:lfnm)
            goto 99                                       !--> ERROR
          endif
        endif
!
! kvf: vertical vortex force term (u^St du/dz)
!
        if (wrt(indxKVF)) then
          lvar=lenstr(vname(1,indxKVF))
          ierr=nf_inq_varid (ncid, vname(1,indxKVF)(1:lvar), vidKVF)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxKVF)(1:lvar), fname(1:lfnm)
            goto 99                                       !--> ERROR
          endif
        endif
!
! calP: surface pressure correction term appeared in prsgrd term
!
        if (wrt(indxCALP)) then
          lvar=lenstr(vname(1,indxCALP))
          ierr=nf_inq_varid (ncid, vname(1,indxCALP)(1:lvar), vidCALP)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxCALP)(1:lvar), fname(1:lfnm)
            goto 99                                       !--> ERROR
          endif
        endif
!
! Kapsrf: surface Bernoulli head term appeared in prsgrd term
!
        if (wrt(indxKAPS)) then
          lvar=lenstr(vname(1,indxKAPS))
          ierr=nf_inq_varid (ncid, vname(1,indxKAPS)(1:lvar), vidKAPS)
          if (ierr .ne. nf_noerr) then
            write(*,1) vname(1,indxKAPS)(1:lvar), fname(1:lfnm)
            goto 99                                       !--> ERROR
          endif
        endif
# endif  /* SOLVE3D */

      goto 100
  99  return

#ifdef AVRH
   1  format(/1x,'### ERROR: def_avg :: Cannot find variable ''',
#else
   1  format(/1x,'### ERROR: def_his :: Cannot find variable ''',
#endif
     &                            A, ''' in file ''', A, '''.'/)

 100  continue

      end subroutine def_his_nf_inq_varid


! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

#ifdef ANA_WEC_FRC
      subroutine ana_init_wec_frc_tile (istr,iend,jstr,jend)
      implicit none

      integer, intent(in) :: istr, iend, jstr, jend
      integer :: i, j, k

      real :: kkw,Clw,Cg


      ! compute_auxiliary_bounds.h is needed for istrR, jstrR, etc
      ! This had to appear after other variables declared else wont compile
#include "compute_auxiliary_bounds.h"


      kkw = 2.0*pi   ! Changed kw to kkw so as not to conflict part 2 also has kw
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
            ubar(i,j,2) = ubar(i,j,1) !! not actually sure if that is needed

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

      end subroutine ana_init_wec_frc_tile
#endif /* ANA_WEC_FRC */

! --------------------------------------------------------------------
! --------------------------------------------------------------------

#ifdef ANA_WEC_FRC
      subroutine ana_wec_frc_tile (istr,iend,jstr,jend)
      ! Analytical WEC forcing
      implicit none

      integer, intent(in) :: istr, iend, jstr, jend
      integer :: i, j
      ! Declare coefficient variables
      real :: cff1,cff2

      real :: kkw,Clw,Cg

      real :: dx,x,timed
      real :: amp,Lx


         ! get universal scalars, g and pi
! # include "scalars.h" - included in top of module
         ! get iSW_corn variable
! # include "param.h" - included in top of module
         ! get pm coordinate transformer from grid.h
! # include "grid.h" - included in top of module
         ! note haven't included DH's forces.h as all wave parameters declared in WEC

         ! For PART 2 below: the include ocean.h files for mrl_wci part
         ! included in top of module
!# include "ocean2d.h"
!# ifdef SOLVE3D
!#   include "ocean3d.h"
!# endif

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

      end subroutine ana_wec_frc_tile
#endif /* ANA_WEC_FRC */

! --------------------------------------------------------------------
! --------------------------------------------------------------------

#ifdef SOLVE3D
!
! Vertical Stokes drift velocity at rho-point
! ===========================================

      subroutine wstokes_wec (tile)
      implicit none
! already included at top of module # include "param.h"
      integer tile
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call wstokes_wec_tile(istr,iend,jstr,jend)
!      return
      end subroutine wstokes_wec


      subroutine wstokes_wec_tile (istr,iend,jstr,jend)

      implicit none
!! # include "param.h"
      integer istr,iend,jstr,jend, i,j,k
      real Wrk(0:N) ! DevinD - move to top of module once code is working
!! # include "grid.h"
!! # include "ocean3d.h"
!! # include "forces.h"
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

      end subroutine wstokes_wec_tile

#endif  /* SOLVE3D */

       ! --------------------------------------------------------------------
       ! --------------------------------------------------------------------
       ! --------------------------------------------------------------------

      subroutine compute_wec_rhs_uv_terms (istr,iend,jstr,jend,ru,rv,
     &                                           wrkone, wrktwo)

      implicit none

        integer, intent(in) :: istr, iend, jstr, jend
         ! Declare looping indices
        real, dimension(PRIVATE_2D_SCRATCH_ARRAY,N) :: ru,rv
        real, dimension(PRIVATE_2D_SCRATCH_ARRAY) ::  wrkone,wrktwo


        integer :: i, j, k

        real :: cff, cff1

      real :: gam1, gam2

# include "compute_auxiliary_bounds.h" /* Taken from top of rhs3d35S.F */
# define VF_ROBUST /* Taken from top of rhs3d35S.F */

!
! Add combined vortex-force and advection terms & breaking terms
! 2nd- and 4th-order centered schemes, or 2nd with local weighted
! filter (with VF_ROBUST) are available.
!
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


       end subroutine compute_wec_rhs_uv_terms


       ! --------------------------------------------------------------------
       ! --------------------------------------------------------------------
       ! --------------------------------------------------------------------

      subroutine wrt_his_wec (record, ierr, may_day_flag)

      ! use netcdf, only: nf_noerr ! didn't work for nf_noerr

      implicit none

      integer :: ierr, record, lvar, ncdf_write, lenstr, may_day_flag, tile
      integer, parameter :: nf_noerr=0 ! Taken from netcdf.inc to avoid include

      if (wrthis(indxSUP)) then
!        write(*,*) 'wrt_his: wrtSUP - ncidhis: ', ncidhis,
!     &         ' *** hisSUP: ', hisSUP
!        ierr=ncdf_write (sup, ncidhis, hisSUP, record, r2dvar)
        ierr=ncdf_write (nchis, hisSUP, record, r_var, sup, 1)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxSUP))
          write(*,1) vname(1,indxSUP)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his: wrtSUP' MYID
      endif

      if (wrthis(indxUST2D)) then
!        ierr=nf_fwrite (ust2d, ncidhis, hisUST2D, record, u2dvar)
        ierr=ncdf_write (nchis, hisUST2D, record, u_var, ust2d, 1)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxUST2D))
          write(*,1) vname(1,indxUST2D)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his: wrtUST2D' MYID
      endif

      if (wrthis(indxVST2D)) then
!        ierr=nf_fwrite (vst2d, ncidhis, hisVST2D, record, v2dvar)
        ierr=ncdf_write (nchis, hisVST2D, record, v_var, vst2d, 1)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxVST2D))
          write(*,1) vname(1,indxVST2D)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his: wrtUST2D' MYID
      endif
# ifdef SOLVE3D
      if (wrthis(indxUST)) then
!        ierr=nf_fwrite (ust, ncidhis, hisUST,  record,  u3dvar)
        ierr=ncdf_write (nchis, hisUST,  record,  u_var, ust(START_2D_ARRAY,1), N)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxUST))
          write(*,1) vname(1,indxUST)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his: wrtUST' MYID
      endif
      if (wrthis(indxVST)) then
!        ierr=nf_fwrite (vst, ncidhis, hisVST,  record,  v3dvar)
        ierr=ncdf_write (nchis, hisVST,  record,  v_var, vst, N)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxVST))
          write(*,1) vname(1,indxVST)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his: wrtVST' MYID
      endif
      if (wrthis(indxWST)) then
        do tile=0,NSUB_X*NSUB_E-1
          call wstokes_wec(tile)
        enddo
!        ierr=nf_fwrite (wst, ncidhis, hisWST, record, r3dvar)
        ierr=ncdf_write (nchis, hisWST, record, r_var, wst, N)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxWST))
          write(*,1) vname(1,indxWST)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his: wrtWST' MYID
      endif
      if (wrthis(indxAkb)) then
!        ierr=nf_fwrite (Akb, ncidhis, hisAkb, record, w3dvar)
        ierr=ncdf_write (nchis, hisAkb, record, r_var, Akb, N+1)
! DevinD - left this uncommented as don't know what w3dvar is in new code.
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxAkb))
          write(*,1) vname(1,indxAkb)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his: wrtAkb' MYID
      endif
      if (wrthis(indxAkw)) then
!        ierr=nf_fwrite (Akw, ncidhis, hisAkw, record, w3dvar)
        ierr=ncdf_write (nchis, hisAkw, record, r_var, Akw, N+1)
! DevinD - left this uncommented as don't know what w3dvar is in new code.
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxAkw))
          write(*,1) vname(1,indxAkw)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his: wrtAkw' MYID
      endif
      if (wrthis(indxKVF)) then
!        ierr=nf_fwrite (kvf, ncidhis, hisKVF, record, r3dvar)
        ierr=ncdf_write (nchis, hisKVF, record, r_var, kvf, N)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxKVF))
          write(*,1) vname(1,indxKVF)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his: wrtKVF' MYID
      endif
      if (wrthis(indxCALP)) then
!        ierr=nf_fwrite (calP, ncidhis, hisCALP, record, r2dvar)
        ierr=ncdf_write (nchis, hisCALP, record, r_var, calP, 1)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxCALP))
          write(*,1) vname(1,indxCALP)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his: wrtCALP' MYID
      endif
      if (wrthis(indxKAPS)) then
!        ierr=nf_fwrite (Kapsrf, ncidhis, hisKAPS, record, r2dvar)
        ierr=ncdf_write (nchis, hisKAPS, record, r_var, Kapsrf, 1)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxKAPS))
          write(*,1) vname(1,indxKAPS)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his: wrtKAPS' MYID
      endif
# endif  /* SOLVE3D */


  1   format(/1x, '### ERROR: wrt_his :: Cannot write variable ''',
     &              A, ''' into history file, rec =', i6, 3x,A,i4)
      goto 100
  99  if (may_day_flag == 0) may_day_flag=3
 100  continue

      end subroutine wrt_his_wec


       ! --------------------------------------------------------------------
       ! --------------------------------------------------------------------
       ! --------------------------------------------------------------------

#ifdef AVERAGES
      subroutine wrt_avg_wec (record, ierr, may_day_flag)

      ! Taken from wrt_avg.F in old code. Called from wrt_avg.F

      implicit none

      integer :: ierr, record, lvar, ncdf_write, lenstr, may_day_flag
      integer, parameter :: nf_noerr=0 ! Taken from netcdf.inc to avoid include

      if (wrtavg(indxSUP)) then
!        ierr=nf_fwrite (sup_avg(START_2D_ARRAY), ncidavg, avgSUP,
!     &                                            record, r2dvar)
        ierr=ncdf_write(ncavg, avgSUP, record, r_var,
     &                            sup_avg(START_2D_ARRAY), 1)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxSUP))
          write(*,1) vname(1,indxSUP)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
      endif
      if (wrtavg(indxUST2D)) then
        ierr=ncdf_write(ncavg, avgUST2D, record, u_var,
     &                           ust2d_avg(START_2D_ARRAY), 1)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxUST2D))
          write(*,1) vname(1,indxUST2D)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
      endif
      if (wrtavg(indxVST2D)) then
        ierr=ncdf_write(ncavg, avgVST2D, record, v_var,
     &                           vst2d_avg(START_2D_ARRAY), 1)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxVST2D))
          write(*,1) vname(1,indxVST2D)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
      endif

# ifdef SOLVE3D
      if (wrtavg(indxUST)) then
        ierr=ncdf_write(ncavg, avgUST, record, u_var,
     &                           ust_avg(START_2D_ARRAY,1), N)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxUST))
          write(*,1) vname(1,indxUST)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
      endif
      if (wrtavg(indxVST)) then
        ierr=ncdf_write(ncavg, avgVST, record, v_var,
     &                           vst_avg(START_2D_ARRAY,1), N)
        if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxVST))
          write(*,1) vname(1,indxVST)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
      endif
!      if (wrtavg(indxWST)) then
!        ierr=nf_fwrite (wst_avg(START_2D_ARRAY,1), ncidavg, avgWST,
!     &                                              record, r3dvar)
!        if (ierr .ne. nf_noerr) then
!         lvar=lenstr(vname(1,indxWST))
!          write(*,1) vname(1,indxWST)(1:lvar), record MYID
!          goto 99                                         !--> ERROR
!        endif
!      endif
      if (wrtavg(indxAkb)) then
        ierr=ncdf_write(ncavg, avgAkb, record, r_var,
     &                           akb_avg(START_2D_ARRAY,0), N+1 )
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxAkb))
          write(*,1) vname(1,indxAkb)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
      endif
      if (wrtavg(indxAkw)) then
        ierr=ncdf_write(ncavg,avgAkw, record, r_var,
     &                          akw_avg(START_2D_ARRAY,0), N+1 )
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxAkw))
          write(*,1) vname(1,indxAkw)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
      endif
      if (wrtavg(indxKVF)) then
        ierr=ncdf_write(ncavg, avgKVF,record, r_var,
     &                           kvf_avg(START_2D_ARRAY,1), N)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxKVF))
          write(*,1) vname(1,indxKVF)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
      endif
      if (wrtavg(indxCALP)) then
        ierr=ncdf_write(ncavg, avgCALP,record, r_var,
     &                           calp_avg(START_2D_ARRAY), 1)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxCALP))
          write(*,1) vname(1,indxCALP)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
      endif
      if (wrtavg(indxKAPS)) then
        ierr=ncdf_write(ncavg, avgKAPS,record, r_var,
     &                           kaps_avg(START_2D_ARRAY), 1)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxKAPS))
          write(*,1) vname(1,indxKAPS)(1:lvar), record MYID
          goto 99                                         !--> ERROR
        endif
      endif
# endif  /* SOLVE3D */

  1   format(/1x, '### ERROR: wrt_avg :: Cannot write variable ''',
     &             A, ''' into averages file, rec =', i6, 3x,A,i4)
      goto 100
  99  if (may_day_flag == 0) may_day_flag=3
 100  continue

      end subroutine wrt_avg_wec
#endif /* AVERAGES */

       ! --------------------------------------------------------------------
       ! --------------------------------------------------------------------
       ! --------------------------------------------------------------------

#ifndef ANA_WEC_FRC
!-------------------------------------------------------------------------
      subroutine read_wec_frc(fname,irec,idx,wec_time)
      ! This version is intended to exclusively read WEC variables that are
      ! part of 'Leonel's method'
      use netcdf
      implicit none

      ! input/output
      character(len=*),intent(in) :: fname
      integer         ,intent(in) :: irec,idx
      ! local
      integer ierr,ncid
      integer varid 
      ! external function
      integer ncdf_read
      ! Set wec_time for read in variables
      real,dimension(1):: wec_time

!      integer,parameter :: r_var = 0 !! grid type of var (u_var=1,v_var=2,q_var=3)
        ! -> comes from ncvars.h -> read_write module.

        ! ifile and irec should be in the module pre-amble, initial as 0
        ! Look for a file that has the requested variable and search
        ! for the last record for which wec_tim < time.

      ! open the file (maybe check if the file is open already)
      ! =============
      ierr=nf90_open(frcfile(ifile), nf90_nowrite, ncid)

      ! read the time
      ! =============
      ierr = nf90_inq_varid(ncid, wec_tim_name, varid) ! Get time ID (varid)
      ! nf90_get_var needs an array for last 2 terms hence syntax
      ierr = nf90_get_var(ncid, varid, wec_time, (/irec/), (/1/))
      wec_time = wec_time * 86400 ! covnert time in days to seconds

      ! read the vars
      ! =============
      ierr = nf90_inq_varid (ncid,wec_awv_name,varid) ! Wave amplitude
      ierr = ncdf_read(ncid,varid,irec,r_var,wwag(START_2D_ARRAY,idx),1)
      if(ierr /= 0) stop 'ERROR read_wec_frc'

      ierr = nf90_inq_varid (ncid,wec_dwv_name,varid) ! Wave direction
      ierr = ncdf_read(ncid,varid,irec,r_var,wwdg(START_2D_ARRAY,idx),1)
      if(ierr /= 0) stop 'ERROR read_wec_frc'

      ierr = nf90_inq_varid (ncid,wec_pwv_name,varid) ! Wave period
      ierr = ncdf_read(ncid,varid,irec,r_var,wwpg(START_2D_ARRAY,idx),1)
      if(ierr /= 0) stop 'ERROR read_wec_frc'

      ierr = nf90_inq_varid (ncid,wec_usa_name,varid) ! ust2d
      ierr = ncdf_read(ncid,varid,irec,r_var,wwUSg(START_2D_ARRAY,idx),1)
      if(ierr /= 0) stop 'ERROR read_wec_frc'

      ierr = nf90_inq_varid (ncid,wec_vsa_name,varid) ! vst2d
      ierr = ncdf_read(ncid,varid,irec,r_var,wwVSg(START_2D_ARRAY,idx),1)
      if(ierr /= 0) stop 'ERROR read_wec_frc'

      ierr = nf90_inq_varid (ncid,wec_ust_name,varid) ! ustr0
      ierr = ncdf_read(ncid,varid,irec,r_var,wwUS0g(START_2D_ARRAY,idx),1)
      if(ierr /= 0) stop 'ERROR read_wec_frc'

      ierr = nf90_inq_varid (ncid,wec_vst_name,varid) ! vstr0
      ierr = ncdf_read(ncid,varid,irec,r_var,wwVS0g(START_2D_ARRAY,idx),1)
      if(ierr /= 0) stop 'ERROR read_wec_frc'

      ierr = nf90_inq_varid (ncid,wec_ebwv_name,varid) ! eb
      ierr = ncdf_read(ncid,varid,irec,r_var,wweb(START_2D_ARRAY,idx),1)
      if(ierr /= 0) stop 'ERROR read_wec_frc'

      ierr = nf90_inq_varid (ncid,wec_edwv_name,varid) ! ed
      ierr = ncdf_read(ncid,varid,irec,r_var,wwed(START_2D_ARRAY,idx),1)
      if(ierr /= 0) stop 'ERROR read_wec_frc'

      ierr = nf90_inq_varid (ncid,wec_lwv_name,varid) ! lmv
      ierr = ncdf_read(ncid,varid,irec,r_var,wwlm(START_2D_ARRAY,idx),1)
      if(ierr /= 0) stop 'ERROR read_wec_frc'

      ierr = nf90_inq_varid (ncid,wec_sup_name,varid) ! sup
      ierr = ncdf_read(ncid,varid,irec,r_var,wwsup(START_2D_ARRAY,idx),1)
      if(ierr /= 0) stop 'ERROR read_wec_frc'

      end subroutine read_wec_frc
#endif /* !ANA_WEC_FRC */

!-------------------------------------------------------------------------
! #ifndef ANA_WEC_FRC
      subroutine set_wec_frc_tile (istr,iend,jstr,jend,ierr)
      implicit none
!      ! import/export
      integer istr,iend,jstr,jend ! ,intent(in) ::
      integer ierr ! ,intent(out)::
!      ! local
      integer            :: tmp,i,j
      real cff1,cff2, cos_angle, sin_angle

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

          !! First time, we need to read an extra slice of data
          !! Read in the last record with wec_time < tday
          !! this routine provides irec = irec - 1 for first time.
          call find_rec_new_mod(time,wec_awv_name,wec_tim_name,
     &                               ifile,irec)

C$OMP MASTER
          ! For the first time irec = irec - 1, as above, such that
          ! wec_tim(it1) < time.
          call read_wec_frc(frcfile(ifile),irec,it1,wec_tim(it1)) !! reads raw input values
C$OMP END MASTER
C$OMP BARRIER

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

        endif ! End of if first step

        call find_rec_new_mod(time,wec_awv_name,wec_tim_name,
     &                             ifile,irec)

        ! Read in the first record with wec_time > tday
        ! Remember, its alternates between index 1 & 2 for each reading
        call read_wec_frc(frcfile(ifile),irec,it2,wec_tim(it2))

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

            ! pre calculate angles to save costs of cos/sin routine
            cos_angle = cos(deg2rad*angler(i,j))
            sin_angle = sin(deg2rad*angler(i,j))

            ! Get x & y components of wwUSg and wwVSg
            ustr2d(i,j)=(cff1*wwUSg(i,j,it1)+cff2*wwUSg(i,j,it2))
     &                                         *cos_angle
     &                +(cff1*wwVSg(i,j,it1)+cff2*wwVSg(i,j,it2))
     &                                         *sin_angle ! U-east to xi-dir depth ave. Stokes drift
            vstr2d(i,j)=-(cff1*wwUSg(i,j,it1)+cff2*wwUSg(i,j,it2))
     &                                         *sin_angle
     &                +(cff1*wwVSg(i,j,it1)+cff2*wwVSg(i,j,it2))
     &                                         *cos_angle! V-north to eta-dir depth ave. Stokdes drift

            ustr0(i,j)=(cff1*wwus0g(i,j,it1)+cff2*wwus0g(i,j,it2))
     &                                         *cos_angle
     &                +(cff1*wwvs0g(i,j,it1)+cff2*wwvs0g(i,j,it2))
     &                                         *sin_angle! U-east to xi-dir depth ave. Stokes drift
            vstr0(i,j)=-(cff1*wwus0g(i,j,it1)+cff2*wwus0g(i,j,it2))
     &                                         *sin_angle
     &                +(cff1*wwvs0g(i,j,it1)+cff2*wwvs0g(i,j,it2))
     &                                         *cos_angle! V-north to eta-dir depth ave. Stokdes drift

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

      end subroutine set_wec_frc_tile
! #endif /* !ANA_WEC_FRC */
!-------------------------------------------------------------------------

#ifndef ANA_WEC_FRC
      subroutine set_wec_inp2vars (istr,iend,jstr,jend,idx)

      ! This is a partial copy of set_wwave_tile from get_wwave.F in old code.

      implicit none

      !Input
      integer idx ! Index relating to read in variables set 1 or set 2
      integer istr,iend,jstr,jend,i,j

      real wdir ! interim variable to transfor wave direction

! Needed for istrR,iendR, etc
# include "compute_extended_bounds.h"


! Compute bed wave orbital velocity (uorb,vorb) and wave frequency,
! wfrq, used with BBL, and/or Stokes transport (wvstx,wvsty) and wave
! set-up, wvsup, for MRL_WCI.  Dispersion relation is approximated by
! an explicit polynomial function proposed by Dean & Dalrymple (1991).

      do j=jstrR,jendR
        do i=istrR,iendR

          ! WAVE PARAMETERS
          ! ---------------

          wwag(i,j,idx) =2.0*wwag(i,j,idx)  ! Awave should be a 1/2 of Hrms
# ifdef MASKING
     &                           *rmask(i,j)
# endif

              ! Calculate wdir once as used twice below
          wdir=-(deg2rad*wwdg(i,j,idx)-1.5*pi)  ! FROM N (WAVE_OFFLINE) to travel to 0deg --> x-dir
# ifdef CURVGRID
     &                    -deg2rad*angler(i,j)
# endif
              ! Direction (wwdg) is read in as one variable so can't be
              ! reused for dir-x and dir-y
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

          wwUSg(i,j,idx) = wwUSg(i,j,idx)*rmask(i,j)
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

        enddo
      enddo

      end subroutine set_wec_inp2vars

#endif /* !ANA_WEC_FRC */

      ! -------------------------------------------
      ! -------------------------------------------
      ! -------------------------------------------

#ifdef AVERAGES
      subroutine set_avg_wec_tile (istr,iend,jstr,jend,cff,cff1
     &                                 ,jstrR,jendR,istrR,iendR)

      ! Taken from set_avg.F in old code.

! Compute time-averaged fields within a tile.

      implicit none

      integer istr,iend,jstr,jend,i,j,jstrR,jendR,istrR,iendR
# ifdef SOLVE3D
     &                       , k
# endif
      real cff,cff1

        if (wrtavg(indxSUP)) then
          do j=jstrR,jendR
            do i=istrR,iendR
              sup_avg(i,j)= cff1*sup_avg(i,j) +cff*sup(i,j)
            enddo
          enddo
        endif
        if (wrtavg(indxUST2D)) then
          do j=jstrR,jendR
            do i=istrR,iendR
              ust2d_avg(i,j)= cff1*ust2d_avg(i,j) +cff*ust2d(i,j)
            enddo
          enddo
        endif
        if (wrtavg(indxVST2D)) then
          do j=jstrR,jendR
            do i=istrR,iendR
              vst2d_avg(i,j)= cff1*vst2d_avg(i,j) +cff*vst2d(i,j)
            enddo
          enddo
        endif

#  ifdef SOLVE3D
        if (wrtavg(indxUST)) then
          do k=1,N
            do j=jstrR,jendR
              do i=istrR,iendR
                ust_avg(i,j,k)= cff1*ust_avg(i,j,k) +cff*ust(i,j,k)
              enddo
            enddo
          enddo
        endif
        if (wrtavg(indxVST)) then
          do k=1,N
            do j=jstrR,jendR
              do i=istrR,iendR
                vst_avg(i,j,k)= cff1*vst_avg(i,j,k) +cff*vst(i,j,k)
              enddo
            enddo
          enddo
        endif
        if (wrtavg(indxWST)) then
          do k=1,N
            do j=jstrR,jendR
              do i=istrR,iendR
                wst_avg(i,j,k)= cff1*wst_avg(i,j,k) +cff*wst(i,j,k)
              enddo
            enddo
          enddo
        endif
        if (wrtavg(indxAkb)) then
          do k=0,N
            do j=jstrR,jendR
              do i=istrR,iendR
                akb_avg(i,j,k)= cff1*akb_avg(i,j,k) +cff*Akb(i,j,k)
              enddo
            enddo
          enddo
        endif
        if (wrtavg(indxAkw)) then
          do k=0,N
            do j=jstrR,jendR
              do i=istrR,iendR
                akw_avg(i,j,k)= cff1*akw_avg(i,j,k) +cff*Akw(i,j,k)
              enddo
            enddo
          enddo
        endif
        if (wrtavg(indxKVF)) then
          do k=1,N
            do j=jstrR,jendR
              do i=istrR,iendR
                kvf_avg(i,j,k)= cff1*kvf_avg(i,j,k) +cff*kvf(i,j,k)
              enddo
            enddo
          enddo
        endif
        if (wrtavg(indxCALP)) then
          do j=jstrR,jendR
            do i=istrR,iendR
              calp_avg(i,j)= cff1*calp_avg(i,j) +cff*calP(i,j)
            enddo
          enddo
        endif
        if (wrtavg(indxKAPS)) then
          do j=jstrR,jendR
            do i=istrR,iendR
              kaps_avg(i,j)= cff1*kaps_avg(i,j) +cff*Kapsrf(i,j)
            enddo
          enddo
        endif
#  endif  /* SOLVE3D */

      end subroutine set_avg_wec_tile
#endif /* AVERAGES */

      ! -------------------------------------------
      ! -------------------------------------------
      ! -------------------------------------------

      subroutine wec_forces (tile)

      ! This is mrl_wci (tile) from mrl_wci.F in old code.

      implicit none

      integer  tile
!# include "param.h" ! DevinD - already included in top of module
# include "private_scratch.h"
# include "compute_tile_bounds.h"

      call wec_forces_tile (istr,iend,jstr,jend,
#  ifdef SOLVE3D
     &                              A3d(1,1),A3d(1,2),A3d(1,3),A3d(1,4),
#  endif
     &            A2d(1,1),A2d(1,2),A2d(1,3),A2d(1,4),A2d(1,5),A2d(1,6),
     &                   A2d(1,7),A2d(1,8),A2d(1,9),A2d(1,10),A2d(1,11),
     &                                              A2d(1,12),A2d(1,13))

      end subroutine wec_forces


      ! -------------------------------------------
      ! -------------------------------------------
      ! -------------------------------------------


      subroutine wec_forces_tile (istr,iend,jstr,jend,
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
      implicit none
!# include "param.h"
!# include "forces.h"
!# include "grid.h"
!# include "scalars.h"
!# include "ocean2d.h"
# ifdef SOLVE3D
!#   include "ocean3d.h"
!#  if defined WAVE_FRICTION && defined BODY_FRICTION
!#   include "mixing.h"
!#  endif
# endif
!# ifdef ANA_BRY
!#   include "boundary.h"
!# endif
!# ifdef WKB_WWAVE
!#   include "wkb_wwave.h"
!# endif
# if defined BBL_S95 || defined BBL_F00 || defined BBL_QUAD2D
#  include "coupling.h"
# endif
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

! TRANSLATIONS AND MANIPULATIONS OF WEV VARIABLES FROM RAW VALUES IN
! FORCING FILE TO NEEDED VARIABLES HERE. AS PER SET_WWAVE_TILE (OLD) OR
! SET_WEC_FRC_TILE (NEW).

! Small section taken from old code's set_wwave_tile subroutine
! Convert offline 2d stokes drift at rho-points to u-points.
! ==========================================================

      do j=jstrR,jendR
        do i=istr,iendR
          ust2d(i,j)=0.5*(ustr2d(i-1,j)+ustr2d(i,j))
        enddo
      enddo
      do j=jstr,jendR
        do i=istrR,iendR
          vst2d(i,j)=0.5*(vstr2d(i,j-1)+vstr2d(i,j))
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
     &          -exp(-2.*keff(i,j)*(z_w(i,j,0)+h(i,j)+Dstp(i,j)))
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
     &            -exp(-2.*keff(i,j)*(z_w(i,j,k)+h(i,j)+Dstp(i,j)))
            vst_r(i,j,k)=ust_rm(i,j,k)*wdre(i,j)*
     &                sqrt(tanh(abs(z_r(i,j,k)*keff(i,j)/2.)))+
     &            (1.-sqrt(tanh(abs(z_r(i,j,k)*keff(i,j)/2.))))*
     &             vst_ker*inv_dz(i,j,k)*( cff2-cff1 )
c             vst_r(i,j,k)= vst_ker*inv_dz(i,j,k)*( cff2-cff1 ) ! test
c             vst_r(i,j,k)= ust_r(i,j,k)*wdre(i,j) ! test
            ust_r(i,j,k)=ust_rm(i,j,k)*wdrx(i,j)*
     &                sqrt(tanh(abs(z_r(i,j,k)*keff(i,j)/2.)))+
     &            (1.-sqrt(tanh(abs(z_r(i,j,k)*keff(i,j)/2.))))*
     &             ust_ker*inv_dz(i,j,k)*( cff2-cff1 )
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
     &           sqrt(ustr2d(i,j)**2+vstr2d(i,j)**2)*Dstp(i,j))
             else
             FAC=1.
             endif
!             write(*,*) 'FAC=', FAC
             end do ! while
             fac2d(i,j)=FAC
          if (kw(i,j)*Dstp(i,j).lt.2.*pi/20.
     &            .or.abs(ustm0-ustm0w)/ustm0w.lt..1) then
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


      end subroutine wec_forces_tile

#endif /* WEC for whole module */
        
      end module wec_frc
