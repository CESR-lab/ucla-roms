#include "cppdefs.h"
      module wec_frc
      ! Wave Effect on Current (WEC) module which provides all the terms that enter
      ! in the momentum equations

       use dimensions, nt1 => nt, mynode1 => mynode, 
     &     rmask1 => rmask, max_frc1 => max_frc, grdname1 => grdname,
     &     frcfile1 => frcfile


       implicit none
! DevinD - need the following includes for GLOBAL_2D_ARRAY to work

#if defined WEC /* WEC for whole module! */

#include "param.h"
#include "grid.h"

! DevinD - added 05/07 for init_scalars_wec subroutine
#include "ncvars.h"

! DevinD - get universal scalars, g and pi
# include "scalars.h"
! DevinD - get ocean variables
# include "ocean2d.h"
# ifdef SOLVE3D
#   include "ocean3d.h"
# endif

       private

       ! Make all public variables global
       save

      ! Variables usable by code that uses this module


      !!!
      !###### DevinD - variables added as per DH's WAVE_PACKET forces.h
      !                for old MRL_WCI flag ######
      !!!

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
!
!--------------------------------------------------------------------
!  PUBLIC WAVE AVERAGED QUANTITIES THAT ARE USED BY THE CODE
!--------------------------------------------------------------------
!  2D  |  ust2d    |   xi-direction Stokes transport (u-point)
!  2D  |  vst2d    |  eta-direction Stokes transport (v-point)
!  3D  |  ust      |   xi-direction 3D Stokes drift velocity (u-point)
!  3D  |  vst      |  eta-direction 3D Stokes drift velocity (v-point)
!--------------------------------------------------------------------
!  2D  |  sup      |  quasi-static wave set-up (rho-point) 
!  2D  |  calP     |  pressure correction term (rho-point)
!  2D  |  Kapsrf   |  Bernoulli head terrm at the surface (rho-point)
!--------------------------------------------------------------------
!  2D  |  brk2dx   |   xi-direction 2D breaking dissipation (rho)
!  2D  |  brk2de   |  eta-direction 2D breaking dissipation (rho)
!  2D  |  frc2dx   |   xi-direction 2D frictional dissipation (rho)
!  2D  |  frc2de   |  eta-direction 2D frictional dissipation (rho)
!--------------------------------------------------------------------
!  3D  |  Akb      |  breaking-wave-induced additional diffusivity (w-point)
!  3D  |  Akw      |  wave-induced additional diffusivity (rho-point)
!  3D  |  E_pre    |  previous time-step value for Akw estimation (rho)
!--------------------------------------------------------------------
! not used yet due to cpp choices
!--------------------------------------------------------------------
!  3D  |  brk3dx   |   xi-direciton 3D breaking dissipation (rho)
!  3D  |  brk3de   |  eta-direction 3D breaking dissipation (rho)
!  3D  |  frc3dx   |   xi-direciton 3D frictional dissipation (rho)
!  3D  |  frc3de   |  eta-direction 3D frictional dissipation (rho)
!--------------------------------------------------------------------
!

      real, public :: ust2d(GLOBAL_2D_ARRAY)
      real, public :: vst2d(GLOBAL_2D_ARRAY)
# ifdef SOLVE3D
      real, public :: ust(GLOBAL_2D_ARRAY,N)
      real, public :: vst(GLOBAL_2D_ARRAY,N)
# endif  /* SOLVE3D */
      real, public :: calP(GLOBAL_2D_ARRAY)
# ifdef SOLVE3D
      real, public :: Kapsrf(GLOBAL_2D_ARRAY)
      real, public :: sup(GLOBAL_2D_ARRAY)
# endif  /* SOLVE3D */

      real, public :: brk2dx(GLOBAL_2D_ARRAY)
      real, public :: brk2de(GLOBAL_2D_ARRAY)
      real, public :: frc2dx(GLOBAL_2D_ARRAY)
      real, public :: frc2de(GLOBAL_2D_ARRAY)

# ifdef SOLVE3D
      real, public :: Akb(GLOBAL_2D_ARRAY,0:N)
      real, public :: Akw(GLOBAL_2D_ARRAY,0:N)
      real, public :: E_pre(GLOBAL_2D_ARRAY,0:N)

      real, public :: kvf(GLOBAL_2D_ARRAY,N)
# endif  /* SOLVE3D */


# ifdef SOLVE3D
#  ifndef SURFACE_BREAK
      real, public :: brk3dx(GLOBAL_2D_ARRAY,N)
      real, public :: brk3de(GLOBAL_2D_ARRAY,N)
#  endif
#  ifdef BODY_FRICTION
      real, public :: frc3dx(GLOBAL_2D_ARRAY,N)
      real, public :: frc3de(GLOBAL_2D_ARRAY,N)
#  endif
# endif  /* SOLVE3D */


      real :: llmm(GLOBAL_2D_ARRAY)
      real :: ustr2d(GLOBAL_2D_ARRAY)
      real :: vstr2d(GLOBAL_2D_ARRAY)
      real :: ustr0(GLOBAL_2D_ARRAY)
      real :: vstr0(GLOBAL_2D_ARRAY)

      ! NON-ANALYTICAL VARIABLES

#if defined BBL || defined WEC
# if !defined ANA_WEC_FRC && !defined WKB_WWAVE
!--------------------------------------------------------------------
!wwv_time|                    | time of wind-induced waves
!--------------------------------------------------------------------
!  wwag  |                    | wave amplitude [m]
!  wwdg  |                    | wave direction [radians]
!  wwpg  |                    | wave period [s]
!  wwub  |                    | orbital velocity magnitude [m/s]
!  wwfrq |  Two-time-level    | wave frequency [rad/s]
!  wwuob |  point data        | xi-orbital velocity [m/s]
!  wwvob |  for wind induced  ! eta-orbital velocity [m/s]
!  wwdrx |                    | cosine wave direction [ND]
!  wwdre |                    ! sine wave direction [ND]
!  wwhrm |                    ! (RMS) wave height [m]
!  wweb  |                    ! breaking dissipation [m3/s3]
!  wwed  |                    ! frictional dissipation [m3/s3]
!  wwqb  |                    ! fraction of breaking waves [ND]
!  wwUSg |                    ! Depth averaged Stokes drift [m/s]
!  wwVSg |                    ! Depth averaged Stokes drift [m/s]
!  wwus0g|                    ! Surface Stokes drift [m/s]
!  wwvs0g|                    ! Surface Stokes drift [m/s]
!  wwsup |                    | Set down [m]
!  wwlm  |                    | Mean wavelength [m] - for diffusivity approx.
!--------------------------------------------------------------------
!  wveb  |                    | breaking dissipation [m3/s3]
!  wved  |  for present time  | frictional dissipation [m3/s3]
!  wvqb  |  step              | fraction of breaking waves [ND]
!--------------------------------------------------------------------

      real, public :: wwag(GLOBAL_2D_ARRAY,2)
      real, public :: wwdg(GLOBAL_2D_ARRAY,2)
      real, public :: wwpg(GLOBAL_2D_ARRAY,2)
      real, public :: wwfrq(GLOBAL_2D_ARRAY,2)
#  ifdef BBL
#   ifdef WAVE_OFFLINE
      real, public :: wwub(GLOBAL_2D_ARRAY,2)
      real, public :: wwuob(GLOBAL_2D_ARRAY,2)
      real, public :: wwvob(GLOBAL_2D_ARRAY,2)
#   endif
#  endif
#  ifdef WEC
      real, public :: wwhrm(GLOBAL_2D_ARRAY,2)
      real, public :: wwdrx(GLOBAL_2D_ARRAY,2)
      real, public :: wwdre(GLOBAL_2D_ARRAY,2)
      real, public :: wweb(GLOBAL_2D_ARRAY,2)
      real, public :: wveb(GLOBAL_2D_ARRAY)
#   ifdef WAVE_OFFLINE
      real, public :: wved(GLOBAL_2D_ARRAY)
      real, public :: wwed(GLOBAL_2D_ARRAY,2)
      real, public :: wwUSg(GLOBAL_2D_ARRAY,2)
      real, public :: wwVSg(GLOBAL_2D_ARRAY,2)
      real, public :: wwus0g(GLOBAL_2D_ARRAY,2)
      real, public :: wwvs0g(GLOBAL_2D_ARRAY,2)
#   endif
#   if defined WAVE_OFFLINE && defined SURFACE_ROLLER
      real, public :: wvqb(GLOBAL_2D_ARRAY)
      real, public :: wwqb(GLOBAL_2D_ARRAY,2)
#   endif
#   if defined WAVE_OFFLINE && defined SUP_OFF
      real, public :: wwsup(GLOBAL_2D_ARRAY,2)
#  endif
#  if defined WAVE_OFFLINE && defined WAVE_DIFF_FROM_LM
      real, public :: wwlm(GLOBAL_2D_ARRAY,2)
#  endif
#  endif

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
# elif defined ANA_WEC_FRC && defined WEC
      real, public :: wveb(GLOBAL_2D_ARRAY)
      real, public :: wved(GLOBAL_2D_ARRAY)

# endif /* !ANA_WEC_FRC && !WKB_WWAVE */
#endif /* BBL || WEC */


      !#### DevinD - end variables added as per forces.h #####

      ! DevinD - variables used in ana_wec_frc_tile and ana_init_wec_frc_tile

       real :: env

       ! DevinD - coefficient variable for external calculations
       real, public :: cff_ext, cff_ext2

       ! DevinD - timestep variables as per river_frc.F
!       integer              :: it1 = 1  ! step n. Note: should be parameter as unchanging
!       integer              :: it2 = 2  ! step n+1. Note: should be parameter as unchanging
!       real, dimension(2)   :: rtim = [-99,-99] ! used to store both time steps (n, n+1)


      ! ########################################################
      ! DevinD - Add wave terms previously contained in ncvars.h

!      type, public :: wec_ncvars

! ** MRL wave-current interaction variables ***
! < indx???: use between 81 and 100 >
! indxSUP          quasi-static sea-level response (wave set-up/down),
! indxUST2D,indxVST2D      2D depth-averaged Stokes drift velocities
! indxUST,indxVST,indxWST  3D Stokes drift velocities
! indxAkb,indxAkw  breaking eddy viscosity & primary wave eddy diffusivity
! indxKVF          vertical vortex force (u^St du/dz)
! indxCALP,indxKAPS  surface pressure correction & surface Bernoulli head
! indxUs0,indxVs0      Surface Stokes drift velocities
! indxLM               Mean wave length (for diffusivity approximation)
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

      !#### Taken from old init_scalars.F #########

        ! Put in init_scalar_wec subroutine for now.

      !#### END old init_scalars.F #########

! flags taken from def_his.F and used in def_his_def_var_wec & def_his_nf_inq_varid
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


! ###### Taken from averages.h of old code ######
#ifdef AVERAGES
      real sup_avg(GLOBAL_2D_ARRAY)
      real ust2d_avg(GLOBAL_2D_ARRAY)
      real vst2d_avg(GLOBAL_2D_ARRAY)
#  ifdef SOLVE3D
      real ust_avg(GLOBAL_2D_ARRAY,N)
      real vst_avg(GLOBAL_2D_ARRAY,N)
      real akb_avg(GLOBAL_2D_ARRAY,0:N)
      real akw_avg(GLOBAL_2D_ARRAY,0:N)
      real kvf_avg(GLOBAL_2D_ARRAY,N)
      real calp_avg(GLOBAL_2D_ARRAY)
      real kaps_avg(GLOBAL_2D_ARRAY)
#  endif  /* SOLVE3D */
#endif /* AVERAGES */
! ###### END from averages.h of old code ######

      ! strings needed for wrt_his.F if using ncio method
!      character*15, public :: mynode_text
!      character*20, public :: mynode_nc_file
      ! end strings needed for wrt_his.F

      public init_scalars_wec
      public read_inp_wec
#ifndef ANA_WEC_FRC
      public read_wec_frc
      public set_wec_frc_tile
#endif
      public def_his_def_var_wec
      public def_his_nf_inq_varid
      public ana_wec_frc_tile
      public init_arrays_wec_tile
      public ana_init_wec_frc_tile
      public compute_wec_rhs_uv_terms
      public wrt_his_wec
      public wec_forces
#ifdef AVERAGES
      public read_inp_avg_wec
      public set_avg_wec_tile
      public wrt_avg_wec
#endif



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
            wwfrq(i,j,k)=init
#   ifdef BBL
            wwuob(i,j,k)=init
            wwvob(i,j,k)=init
#    ifdef WAVE_OFFLINE
            wwub(i,j,k) =init
#    endif
#   endif
#   ifdef WEC
            wwhrm(i,j,k)=init
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

! --------------------------------------------------------------------
! --------------------------------------------------------------------

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

!      !###### NOW CALCULATE WEC FORCES USING WHRM, etc ####
!  for now, in main.
!      call wec_wave2frc_method1(istr,iend,jstr,jend)

!      call wec_wv2frc_method2(istr,iend,jstr,jend)

      end subroutine ana_wec_frc_tile
! --------------------------------------------------------------------


! --------------------------------------------------------------------
      subroutine compute_wec_rhs_uv_terms (istr,iend,jstr,jend,ru,rv,
     &                                           wrkone, wrktwo)

      implicit none

      integer, intent(in) :: istr, iend, jstr, jend

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
      subroutine read_wec_frc (ierr)  !(time)
      ! Reads WEC forcing data from input file
      ! This version is now exclusively intended to read variables that are
      ! part of 'Leonels method'

      ! Read gridded wind wave height, direction, period at appropriate
      ! time from forcing netCDF file

      implicit none

# define WWAVE_DATA
!# include "param.h"
!# include "scalars.h"
!# include "forces.h"
!# include "ncvars.h"
# include "netcdf.inc"
      real cff
      integer ncid,ifrc, ierr, lfrc,lvar, lenstr, ncdf_read
      character*(max_name_size) frcname

! Initialization: Check, whether forcing netCDF file is already
!===============  opened, an if not, open it. Find and save netCDF
! IDs for relevant variables. Determine whether there is cycling to
! reuse the input data and what is cycling period "ww_cycle", find
! initial cycling index "ww_ncycle" and record index "ww_rec".
! Set initial value for time index "itww" and set both time record
! bounds to large negative artificial values to trigger the logic
! in reading part below.
c#define VERBOSE
#ifdef VERBOSE
      write(*,*) 'enter get_wwave' MYID
#endif

      ifrc=max(ww_file_id,1)
 100  ncid=ncfrc(ifrc) ! Old code was ncidfrc
      frcname=frcfile(ifrc)
      lfrc=lenstr(frcname)
      ierr=nf_noerr

#ifdef VERBOSE
      write(*,*) 'get_wwave : ncid = ',ncid MYID
      write(*,*) 'get_wwave : iic = ',iic,' ntstart = ',ntstart MYID
#endif
      if (iic.eq.ntstart) then !i1
!        write(*,*) 'HERE i1!!!!!!!!!!' ! DevinD debug
        if (ncid.eq.-1) then !i2
          ierr=nf_open (frcname(1:lfrc), nf_nowrite, ncid)
#ifdef VERBOSE
      write(*,*) 'get_wwave : nf_noerr = ',nf_noerr,' ierr = ',ierr MYID
#endif
          if (ierr.eq.nf_noerr) then !i3
            ncfrc(ifrc)=ncid ! DevinD - old code was ncidfrc
#ifdef VERBOSE
      write(*,*) 'get_wwave : opened' MYID
#endif
          else
            write(*,'(/1x,4A/)') '### ERROR in get_wwave :: ',
     &      'Cannot open netCDF file ''', frcname(1:lfrc), '''.'
          endif
        endif
        if (ierr .eq. nf_noerr) then !i4
          lvar=lenstr(vname(1,indxWWA))
          ierr=nf_inq_varid (ncid, vname(1,indxWWA)(1:lvar), wwaid)
#ifdef VERBOSE
      write(*,*) 'get_wwave : check variable ',vname(1,indxWWA)(1:lvar) MYID
#endif
          if (ierr .eq. nf_noerr) then !i5
            ww_file_id=ifrc
            mpi_nonexit_warn write(*,'(1x,5A,I3)') 'Found ''',
     &        vname(1,indxWWA)(1:lvar),
     &        ''' in ''', frcname(1:lfrc), ''', ifrc =', ifrc
            lvar=lenstr(vname(1,indxWWD))
            ierr=nf_inq_varid (ncid,vname(1,indxWWD)(1:lvar),wwdid)
            if (ierr .eq. nf_noerr) then !i6
              lvar=lenstr(vname(1,indxWWP))
              ierr=nf_inq_varid (ncid,vname(1,indxWWP)(1:lvar),wwpid)

#  if defined WAVE_OFFLINE && defined BBL
              if (ierr .eq. nf_noerr) then !i7
                lvar=lenstr(vname(1,indxWUB))
                ierr=nf_inq_varid (ncid,vname(1,indxWUB)(1:lvar),wwuid)
              if (ierr .eq. nf_noerr) then !i8
                lvar=lenstr(vname(1,indxWVB))
                ierr=nf_inq_varid (ncid,vname(1,indxWVB)(1:lvar),wwvid)
#  endif
#  if defined WAVE_OFFLINE && defined WEC
#   if !defined U2010
                if (ierr .eq. nf_noerr) then !9
                  lvar=lenstr(vname(1,indxUST2D))
                  ierr=nf_inq_varid (ncid,vname(1,indxUST2d)(1:lvar),
     &                                                     wwusaid)
                if (ierr .eq. nf_noerr) then !10
                  lvar=lenstr(vname(1,indxVST2D))
                  ierr=nf_inq_varid (ncid,vname(1,indxVST2d)(1:lvar),
     &                                                     wwvsaid)
                if (ierr .eq. nf_noerr) then !11
                  lvar=lenstr(vname(1,indxUs0))
                  ierr=nf_inq_varid (ncid,vname(1,indxUs0)(1:lvar),
     &                                                     wwus0id)
                if (ierr .eq. nf_noerr) then !12
                  lvar=lenstr(vname(1,indxVs0))
                  ierr=nf_inq_varid (ncid,vname(1,indxVs0)(1:lvar),
     &                                                     wwvs0id)
# endif
                if (ierr .eq. nf_noerr) then !i13
                  lvar=lenstr(vname(1,indxWEB))
                  ierr=nf_inq_varid (ncid,vname(1,indxWEB)(1:lvar),
     &                                                     wweid)
#   ifdef WAVE_FRICTION
                  if (ierr .eq. nf_noerr) then !i14
                    lvar=lenstr(vname(1,indxWED))
                    ierr=nf_inq_varid (ncid,vname(1,indxWED)(1:lvar),
     &                                                           wwfid)
#   endif
#   ifdef SURFACE_ROLLER
                  if (ierr .eq. nf_noerr) then !i15
                    lvar=lenstr(vname(1,indxWQB))
                    ierr=nf_inq_varid (ncid,vname(1,indxWQB)(1:lvar),
     &                                                           wwqid)
#   endif
#   ifdef SUP_OFF
                  if (ierr .eq. nf_noerr) then !i16
                    lvar=lenstr(vname(1,indxSUP))
                    ierr=nf_inq_varid (ncid,vname(1,indxSUP)(1:lvar),
     &                                                          wwsupid)
#  endif
#   ifdef WAVE_DIFF_FROM_LM
                  if (ierr .eq. nf_noerr) then !i16
                    lvar=lenstr(vname(1,indxLM))
                    ierr=nf_inq_varid (ncid,vname(1,indxLM)(1:lvar),
     &                                                          wwlmid)
#  endif
#  endif
                    if (ierr .eq. nf_noerr) then !i18
                      ierr=nf_inq_varid (ncid, 'wwv_time', ww_tid)
                      if (ierr .eq. nf_noerr) then !i18
!                        call set_cycle (ncid, ww_tid, ntww, ww_cycle,
!     &                                         ww_ncycle, ww_rec, ierr)
                        call find_record(ncid, ww_tid,
     &                         ww_rec, ntww,
     &                      ww_ncycle, ww_cycle,
     &                                  ifrc, max_frc,  ierr)
      ! DevinD - Example taken from set_stflux.F
!              call find_record(         ncid, stf_tid(itrc),
!     &                         stf_rec(itrc), ntstf(itrc),
!     &                      stf_ncycle(itrc), stf_cycle(itrc),
!     &                                  ifrc, max_frc,  ierr)
                        itww=1
                        wwv_time(1)=-1.E+20
                        wwv_time(2)=-1.E+20
                      else
                        write(*,1) 'wwv_time', frcname(1:lfrc)
                      endif !ie12
#  if defined WAVE_OFFLINE && defined WEC
#   ifdef WAVE_DIFF_FROM_LM
                    else
                      write(*,1) vname(1,indxLM)(1:lvar),
     &                                                  frcname(1:lfrc)
                    endif
#   endif
#   ifdef SUP_OFF
                    else
                      write(*,1) vname(1,indxSUP)(1:lvar),
     &                                                  frcname(1:lfrc)
                    endif
#   endif
#   ifdef SURFACE_ROLLER
                    else
                      write(*,1) vname(1,indxWQB)(1:lvar),
     &                                                  frcname(1:lfrc)
                    endif
#   endif
#   ifdef WAVE_FRICTION
                    else
                      write(*,1) vname(1,indxWED)(1:lvar),
     &                                                  frcname(1:lfrc)
                    endif
#   endif
                  else
                    write(*,1) vname(1,indxWEB)(1:lvar),
     &                                                  frcname(1:lfrc)
                  endif
#   if !defined U2010
                   else
                    write(*,1) vname(1,indxVs0)(1:lvar),
     &                                                  frcname(1:lfrc)
                  endif
                   else
                    write(*,1) vname(1,indxUs0)(1:lvar),
     &                                                  frcname(1:lfrc)
                  endif
                   else
                    write(*,1) vname(1,indxVST2D)(1:lvar),
     &                                                  frcname(1:lfrc)
                  endif
                   else
                    write(*,1) vname(1,indxUST2D)(1:lvar),
     &                                                  frcname(1:lfrc)
                  endif
#  endif
#  endif
#  if defined WAVE_OFFLINE && defined BBL
                else
                  write(*,1) vname(1,indxWVB)(1:lvar),
     &                                                  frcname(1:lfrc)
                endif
                else
                  write(*,1) vname(1,indxWUB)(1:lvar),
     &                                                  frcname(1:lfrc)
                endif
#  endif

              else
                write(*,1) vname(1,indxWWP)(1:lvar),
     &                                   frcname(1:lfrc)
              endif
            else
              write(*,1) vname(1,indxWWD)(1:lvar),
     &                                   frcname(1:lfrc)
            endif
          else
            write(*,1) vname(1,indxWWA)(1:lvar),
     &                                 frcname(1:lfrc)
            if (ifrc .lt. max_frc) then
              ifrc=ifrc+1
              goto 100
            endif
          endif
        endif
      endif
  9   format(1x,a,a,a,a,a,1x,i2)
  1   format(/1x,'WARNING: get_wwave :: Cannot find variable ''',
     &                       A, ''' in netCDF file ''', A, '''.')

!      write(*,*)'set_cycle: ',ntww, ww_cycle, itww, wwv_time(itww),
!     &                    wwaid,wwpid,wwdid,wweid

!
! Read data from the file:  Check if model time is bounded by past
!===== ==== ==== === =====  and future data times: if not, increment
! record and cycling indices, flip time index and read a new portion
! of data. Repeat until model time falls between the two data times.
!
!      do while (wwv_time(itww).lt.time+0.5*dt .and. ierr.eq.nf_noerr) !
!      this results in ### WARNING: set_wwave_tile :: Current model time
!      is outside bounds of 'wwv_time'
      do while (wwv_time(itww).lt.time .and. ierr.eq.nf_noerr)
!        call advance_cycle (ww_cycle, ntww, ww_ncycle,
!     &                                        ww_rec, ierr)
        call advance_record(ww_rec, ntww,
     &                      ww_ncycle, ww_ncycle,
     &                                  ifrc, max_frc,  ierr)
      ! DevinD - Example taken from set_stflux.F of new code
!        call advance_record(   stf_rec(itrc), ntstf(itrc),
!     &                      stf_ncycle(itrc), stf_cycle(itrc),
!     &                                  ifrc, max_frc,  ierr)
        if (ierr.eq.nf_noerr) then
          ierr=nf_get_var1_FTYPE (ncid, ww_tid, ww_rec, cff)
          if (ierr .eq. nf_noerr) then
            itww=min(3-itww,ntww)
            wwv_time(itww)=ww_cycle*ww_ncycle + cff*day2sec
!            ierr=nf_fread (wwag(START_2D_ARRAY,itww), ncid,
!     &                                      wwaid, ww_rec, r2dvar)
            ierr=ncdf_read(ncid, wwaid,  ww_rec, r_var,
     &                                 wwag(START_2D_ARRAY,itww), 1)

!      ! DevinD debug
!      write(*,*) 'wwag(1,1,itww)=?,node', wwag(1,1,itww), mynode

      ! Example from get_grid.F:
!      ierr=ncdf_read(ncid, varid, 0, r_var, h(START_2D_ARRAY), 1)
      ! Example from get_stflux.F
!            ierr=ncdf_read(ncid, stf_id(itrc),  stf_rec(itrc),  r_var,
!     &                     stflxg(START_2D_ARRAY,itstf(itrc),itrc), 1)

            if (ierr .eq. nf_noerr) then
              ierr=ncdf_read(ncid, wwdid, ww_rec, r_var,
     &                                     wwdg(START_2D_ARRAY,itww), 1 )
              if (ierr .eq. nf_noerr) then
                ierr=ncdf_read(ncid,wwpid, ww_rec, r_var,
     &                                     wwpg(START_2D_ARRAY,itww), 1 )

#  if defined WAVE_OFFLINE && defined BBL
                if (ierr .eq. nf_noerr) then
                  ierr=ncdf_read(ncid, wwuid, ww_rec, r_var,
     &                                wwuob(START_2D_ARRAY,itww), 1 )
                if (ierr .eq. nf_noerr) then
                  ierr=ncdf_read(ncid,wwvid, ww_rec, r_var,
     &                                     wwvob(START_2D_ARRAY,itww), 1 )
#  endif
#  if defined WAVE_OFFLINE && defined WEC
#  if !defined U2010
                  if (ierr .eq. nf_noerr) then
                    ierr=ncdf_read(ncid,wwusaid, ww_rec, r_var,
     &                                   wwUSg(START_2D_ARRAY,itww), 1 )
                  if (ierr .eq. nf_noerr) then
                    ierr=ncdf_read(ncid,wwvsaid, ww_rec, r_var,
     &                                   wwVSg(START_2D_ARRAY,itww), 1 )
                  if (ierr .eq. nf_noerr) then
                    ierr=ncdf_read(ncid, wwus0id, ww_rec, r_var,
     &                                  wwus0g(START_2D_ARRAY,itww), 1 )
                  if (ierr .eq. nf_noerr) then
                    ierr=ncdf_read(ncid,wwvs0id, ww_rec, r_var,
     &                                  wwvs0g(START_2D_ARRAY,itww), 1 )
# endif
                  if (ierr .eq. nf_noerr) then
                    ierr=ncdf_read(ncid,wweid, ww_rec, r_var,
     &                                  wweb(START_2D_ARRAY,itww), 1 )
#   ifdef WAVE_FRICTION
                    if (ierr .eq. nf_noerr) then
                      ierr=ncdf_read(ncid,wwfid, ww_rec, r_var,
     &                                    wwed(START_2D_ARRAY,itww), 1 )
#   endif
#   ifdef SURFACE_ROLLER
                    if (ierr .eq. nf_noerr) then
                      ierr=ncdf_read(ncid,wwqid, ww_rec, r_var,
     &                                    wwqb(START_2D_ARRAY,itww), 1 )
#   endif
#   ifdef SUP_OFF
                    if (ierr .eq. nf_noerr) then
                      ierr=ncdf_read(ncid,wwsupid, ww_rec, r_var,
     &                                    wwsup(START_2D_ARRAY,itww), 1)
#   endif
#   ifdef WAVE_DIFF_FROM_LM
                    if (ierr .eq. nf_noerr) then
                      ierr=ncdf_read(ncid,wwlmid, ww_rec, r_var,
     &                                    wwlm(START_2D_ARRAY,itww), 1 )
#   endif
#  endif

                      if (ierr .eq. nf_noerr) then
                       mpi_nonexit_warn
!----- Cigdem (extra & is removed)
     &                  write(*,'(5x,A,28x,A,F12.4,2(1x,A,I4))')
     &                   'get_wwave :: read wwave', 'time =',cff,
     &                   'rec =', ww_rec MYID
                         if (ntww.eq.1) return

#  if defined WAVE_OFFLINE && defined MRL_WCI
#   ifdef WAVE_DIFF_FROM_LM
                      else
                        write(*,2) 'Mean wavelength ', ww_rec
                      endif
#   endif
#   ifdef SUP_OFF
                      else
                        write(*,2) 'Set sup (down) ', ww_rec
                      endif
#   endif
#   ifdef SURFACE_ROLLER
                      else
                        write(*,2) 'fraction of breaking', ww_rec
                      endif
#   endif
#   ifdef WAVE_FRICTION
                      else
                        write(*,2) 'frictional dissipation', ww_rec
                      endif
#   endif
! kaida
                    else
                      write(*,2) 'breaking dissipation', ww_rec
                    endif
#  if !defined U2010
                      else
                        write(*,2) 'North Surface Stokes Drift',ww_rec
                      endif
                      else
                        write(*,2) 'East Surface Stokes Drift',ww_rec
                      endif
                      else
                        write(*,2) 'North depth avg. Stokes Drift',ww_rec
                      endif
                      else
                        write(*,2) 'East depth avg. Stokes Drift',ww_rec
                      endif
#  endif
#  endif
#  if defined WAVE_OFFLINE && defined BBL
                  else
                    write(*,2) 'North wave orbital velocity', ww_rec
                  endif
                  else
                    write(*,2) 'East wave orbital velocity', ww_rec
                  endif
#  endif
                else
                  write(*,2) 'wave period', ww_rec
                endif
              else
                write(*,2) 'wave direction', ww_rec
              endif
            else
              write(*,2) 'wave amplitude', ww_rec
            endif
          else
            write(*,2) 'wwv_time', ww_rec
          endif
        else
          write(*,'(/1x,A,I4,1x,A,I4/7x,4A/7x,2(A,F12.4)/)')
     &      'ERROR in get_wwave :: requested time record ', ww_rec,
     &      'exeeds the last record', ntww,  'available in forcing ',
     &      'netCDF file ''',   frcname(1:lfrc),   '''',  'tdays = ',
     &      tdays,            '  but the last available  wwv_time =',
     &                                        wwv_time(itww)*sec2day
        endif
      enddo
  2   format(/1x,'### ERROR: get_wwave :: Cannot read variable ''', A,
     &                                          ''' for record ', I4/)
# ifdef VERBOSE
      write(*,*) 'return from get_wwave' MYID
# endif
!      return


      end subroutine read_wec_frc

#endif /* !ANA_WEC_FRC */

      ! -----------------------------------------------------------
      ! -----------------------------------------------------------
      ! -----------------------------------------------------------

#ifndef ANA_WEC_FRC
      subroutine set_wec_frc_tile (istr,iend,jstr,jend,ierr)

      ! This is a copy of set_wwave_tile from get_wwave.F in old code.

      implicit none

      integer istr,iend,jstr,jend, ierr, i,j,k, it1,it2,kt1,kt2
      real cff,cff1,cff2, Dstp,cfrq,wdir, kh,khd, eps,xdir,ydir,
     &                                           wramp,m2p,isr2
# ifdef BBL
     &                                            , ab, orb_vel
# endif
# if !defined WAVE_OFFLINE && defined WEC
     &                                          , Btg, Gtg, sbc
# endif
      parameter (eps = 1.0d-10, m2p=1.2804973111, isr2=0.70710678)
# define WWAVE_DATA
!# include "param.h"
!# include "grid.h"
!# include "forces.h"
!# include "scalars.h"
!# include "ocean2d.h"
!
# include "compute_extended_bounds.h"

!
! Set coefficients for interpolation. Check that for the next time
! step [when time=time+dt] both weights will still be positive, and
! if not, set "synchro flag" to signal that new data should be read
! from an appropriate netCDF input file.
! After that either load time-invariant data, or interpolate in time
! or complain about error and signal to quit, if interpolation is
! needed, but not possible.
!
      it1=3-itww
      it2=itww
      cff1=wwv_time(it2)-time
      cff2=time-wwv_time(it1)
      if (cff1.lt.dt) synchro_flag=.true.

! Compute bed wave orbital velocity (uorb,vorb) and wave frequency,
! wfrq, used with BBL, and/or Stokes transport (wvstx,wvsty) and wave
! set-up, wvsup, for MRL_WCI.  Dispersion relation is approximated by
! an explicit polynomial function proposed by Dean & Dalrymple (1991).

!c      if (FIRST_TIME_STEP) then
      if (iic.le.ntstart) then
        kt1=1
        kt2=2
      elseif (cff2.le.dt) then
        kt1=it2
        kt2=it2
      else
        kt1=2
        kt2=1
      endif

      if (kt1.le.kt2) then
# if !defined WAVE_OFFLINE && defined MRL_WCI
        Btg=0.8; Gtg=0.4       ! breaking model by Thornton & Guza (1986)
        sbc=3.0*sqrt(pi)/16.0*g*(Btg**3)/(Gtg**4)/(2.0*pi)
!        sbc=0 ! kaida
# endif
        do k=kt1,kt2,+1
          do j=jstrR,jendR
            do i=istrR,iendR
              Dstp=max(h(i,j),eps)
!              cfrq=2.0*pi/max(wwpg(i,j,k)*m2p,0.1D0) ! mean to peak period
              cfrq=2.0*pi/max(wwpg(i,j,k),0.1D0) !  peak per. to freq.
# ifdef SHOREFACE
              wdir=deg2rad*wwdg(i,j,k)-1.5*pi  ! wdir is dir waves travel FROM N
# else
#  if defined WAVE_OFFLINE
              wdir=-(deg2rad*wwdg(i,j,k)-1.5*pi)  ! FROM N (WAVE_OFFLINE) to travel to 0deg --> x-dir
# endif
# endif
# ifdef CURVGRID
c#  ifndef USWC_WEC
#  if defined WAVE_OFFLINE
     &                    -deg2rad*angler(i,j)  !
#  endif
# endif
              xdir=cos(wdir)
              ydir=sin(wdir)
# ifdef BBL
#  ifdef WAVE_OFFLINE
!             wwuob(i,j,k) and wwvob(i,j,k) read above
#  else
              khd = Dstp*cfrq*cfrq/g
              kh = sqrt( khd*khd + khd/(1.0 + khd*(0.6666666666
     &                   +khd*(0.3555555555 + khd*(0.1608465608
     &                   +khd*(0.0632098765 + khd*(0.0217540484
     &                                   +khd*0.0065407983)))))) )
              ab =wwag(i,j,k)/max(sinh(kh),eps)
              orb_vel =cfrq*ab
              wwuob(i,j,k)=orb_vel*xdir
              wwvob(i,j,k)=orb_vel*ydir
#  endif
#  ifdef MASKING
              wwuob(i,j,k)=wwuob(i,j,k)*rmask(i,j)
              wwvob(i,j,k)=wwuob(i,j,k)*rmask(i,j)
#  endif
# endif

# if defined BBL || defined MRL_WCI
              wwfrq(i,j,k)=cfrq
#  ifdef MASKING
     &                             *rmask(i,j)
#  endif
# endif

# ifdef MRL_WCI
cccc              wwhrm(i,j,k) =2.0*wwag(i,j,k)*isr2 ! Hsig to Hrms conversion
              wwhrm(i,j,k) =2.0*wwag(i,j,k)  ! Awave should be a 1/2 of Hrms
              wwdrx(i,j,k) =xdir
              wwdre(i,j,k) =ydir
#  ifndef WAVE_OFFLINE
!              wweb(i,j,k)  =0.0 ! kaida
              wweb(i,j,k)  =sbc*cfrq*((0.7071*wwhrm(i,j,k))**7)/
     &                      max(Dstp**5,eps)   ! 0.7071 = 1/sqrt(2): Hsig -> Hrms
#  endif
#  ifdef MASKING
              wwhrm(i,j,k) =wwhrm(i,j,k)*rmask(i,j)
              wwdrx(i,j,k) =wwdrx(i,j,k)*rmask(i,j)
              wwdre(i,j,k) =wwdre(i,j,k)*rmask(i,j)
              wweb(i,j,k)  = wweb(i,j,k)*rmask(i,j)
#   if defined WAVE_OFFLINE && defined WAVE_FRICTION
              wwed(i,j,k)  = wwed(i,j,k)*rmask(i,j)
#   endif
#   if defined WAVE_OFFLINE && defined SURFACE_ROLLER
              wwqb(i,j,k)  = wwqb(i,j,k)*rmask(i,j)
#   endif
#   if defined WAVE_OFFLINE && defined SUP_OFF
              wwsup(i,j,k)  = wwsup(i,j,k)*rmask(i,j)
#   endif
#   if defined WAVE_OFFLINE && defined WAVE_DIFF_FROM_LM
              wwlm(i,j,k)  = wwlm(i,j,k)*rmask(i,j)
#   endif
#  endif
# endif
            enddo
          enddo
        enddo

c        if (istr.eq.1.and.jstr.eq.1) write(*,*) 'finish compute wave'
      endif    !<-- kt1 < kt2

# if !defined MRL_WCI && defined BBL && defined WAVE_RAMP
      wramp =tanh(3.0*dt*sec2day*float(iic-ntstart))
# else
      wramp =1.0
# endif

      if (ntww.eq.1) then                 ! Load time-invariant
        if (iic.eq.ntstart) then          ! wind-wave quantities
          do j=jstrR,jendR
            do i=IstR,iendR
# ifdef BBL
# ifdef WAVE_OFFLINE
              uorb(i,j)=wwuob(i,j,itww)*cos(deg2rad*angler(i,j))
     &                +wwvob(i,j,itww)*sin(deg2rad*angler(i,j))! U-east to xi-dir orbital velocity
              vorb(i,j)=-wwuob(i,j,itww)*sin(deg2rad*angler(i,j))
     &                +wwvob(i,j,itww)*cos(deg2rad*angler(i,j))! V-north to eta-dir orbital velocity
              uorb(i,j)=wramp*uorb(i,j)
              vorb(i,j)=wramp*vorb(i,j)
# else
              uorb(i,j)=wwuob(i,j,itww)   !  xi-dir orbital velocity
              vorb(i,j)=wwvob(i,j,itww)   ! eta-dir orbital velocity
# endif
# endif
# if defined BBL || defined MRL_WCI
              wfrq(i,j)=wwfrq(i,j,itww)
# endif
# ifdef MRL_WCI
              whrm(i,j)=wwhrm(i,j,itww)
              wdrx(i,j)=wwdrx(i,j,itww)
              wdre(i,j)=wwdre(i,j,itww)
              wveb(i,j)= wweb(i,j,itww)
!              wveb(i,j)= 0.0 ! kaida
#  if defined WAVE_OFFLINE && defined WAVE_FRICTION
              wved(i,j)= wwed(i,j,itww)
#  endif
#  if defined WAVE_OFFLINE && defined SURFACE_ROLLER
              wvqb(i,j)= wwqb(i,j,itww)
#  endif
#  if defined WAVE_OFFLINE && defined SUP_OFF
              sup(i,j)= wwsup(i,j,itww)
#  endif
#  if defined WAVE_OFFLINE && defined WAVE_DIFF_FROM_LM
              llmm(i,j)= wwlm(i,j,itww)
#  endif
#  if defined WAVE_OFFLINE && !defined U2010
c rotate
              ustr2d(i,j)=wwUSg(i,j,itww)*cos(deg2rad*angler(i,j))
     &                +wwVSg(i,j,itww)*sin(deg2rad*angler(i,j))! U-east to xi-dir depth ave. Stokes drift
              vstr2d(i,j)=-wwUSg(i,j,itww)*sin(deg2rad*angler(i,j))
     &                +wwVSg(i,j,itww)*cos(deg2rad*angler(i,j))! V-north to eta-dir depth ave. Stokdes drift

              ustr2d(i,j)=wramp*ustr2d(i,j)
              vstr2d(i,j)=wramp*vstr2d(i,j)

              ustr0(i,j)=wwus0g(i,j,itww)*cos(deg2rad*angler(i,j))
     &                +wwvs0g(i,j,itww)*sin(deg2rad*angler(i,j))! U-east to xi-dir depth ave. Stokes drift
              vstr0(i,j)=-wwus0g(i,j,itww)*sin(deg2rad*angler(i,j))
     &                +wwvs0g(i,j,itww)*cos(deg2rad*angler(i,j))! V-north to eta-dir depth ave. Stokdes drift
              ustr0(i,j)=wramp*ustr0(i,j)
              vstr0(i,j)=wramp*vstr0(i,j)
#  endif
# endif
            enddo
          enddo
# if defined WEC
#  if defined WAVE_OFFLINE && !defined U2010
cc interpolate ust and ut2d on u,v points
          do j=jstrR,jendR
            do i=istr,iendR
                 ust2d(i,j)=0.5*(ustr2d(i-1,j)+ustr2d(i,j))
c                 ust0(i,j)=0.5*(ustr0(i-1,j)+ustr0(i,j))
        enddo
        enddo
          do j=jstr,jendR
            do i=istrR,iendR
                 vst2d(i,j)=0.5*(vstr2d(i,j-1)+vstr2d(i,j))
c                 vst0(i,j)=0.5*(vstr0(i,j-1)+vstr0(i,j))
        enddo
        enddo

#  endif
# endif
        endif
      elseif (cff1.ge.0. .and. cff2.ge.0.) then
!        write(*,*) 'set_wwave_tile::cff1.ge.0.' ! DevinD debug
        cff=1./(cff1+cff2)
        cff1=cff1*cff
        cff2=cff2*cff
        do j=jstrR,jendR
          do i=istrR,iendR
# ifdef BBL
# ifdef WAVE_OFFLINE
              uorb(i,j)=(cff1*wwuob(i,j,it1)+cff2*wwuob(i,j,it2))
     &                                         *cos(deg2rad*angler(i,j))
     &                +(cff1*wwvob(i,j,it1)+cff2*wwvob(i,j,it2))
     &                                         *sin(deg2rad*angler(i,j))! U-east to xi-dir orbital velocity
              vorb(i,j)=-(cff1*wwuob(i,j,it1)+cff2*wwuob(i,j,it2))
     &                                         *sin(deg2rad*angler(i,j))
     &                +(cff1*wwvob(i,j,it1)+cff2*wwvob(i,j,it2))
     &                 *cos(deg2rad*angler(i,j))! V-north to eta-dir orbital velocity
              uorb(i,j)=wramp*uorb(i,j)
              vorb(i,j)=wramp*vorb(i,j)
# else
            uorb(i,j)=wramp*(cff1*wwuob(i,j,it1)+cff2*wwuob(i,j,it2))
            vorb(i,j)=wramp*(cff1*wwvob(i,j,it1)+cff2*wwvob(i,j,it2))
# endif
# endif

# if defined BBL || defined MRL_WCI
            wfrq(i,j)=cff1*wwfrq(i,j,it1)+cff2*wwfrq(i,j,it2)
# endif
# ifdef MRL_WCI
            whrm(i,j)=wramp*(cff1*wwhrm(i,j,it1)+cff2*wwhrm(i,j,it2))
            wdrx(i,j)=cff1*wwdrx(i,j,it1)+cff2*wwdrx(i,j,it2)
            wdre(i,j)=cff1*wwdre(i,j,it1)+cff2*wwdre(i,j,it2)
            wveb(i,j)=wramp*(cff1*wweb(i,j,it1)+cff2*wweb(i,j,it2))
!            wveb(i,j)=0.0 ! kaida
#  ifdef WAVE_OFFLINE
#   ifdef WAVE_FRICTION
            wved(i,j)=wramp*(cff1*wwed(i,j,it1)+cff2*wwed(i,j,it2))
#   else
            wved(i,j)=0.0
#   endif
# if !defined U2010
c rotate and interpolate in time
              ustr2d(i,j)=(cff1*wwUSg(i,j,it1)+cff2*wwUSg(i,j,it2))
     &                                         *cos(deg2rad*angler(i,j))
     &                +(cff1*wwVSg(i,j,it1)+cff2*wwVSg(i,j,it2))
     &                                         *sin(deg2rad*angler(i,j))! U-east to xi-dir depth ave. Stokes drift
              vstr2d(i,j)=-(cff1*wwUSg(i,j,it1)+cff2*wwUSg(i,j,it2))
     &                                         *sin(deg2rad*angler(i,j))
     &                +(cff1*wwVSg(i,j,it1)+cff2*wwVSg(i,j,it2))
     &                                         *cos(deg2rad*angler(i,j))! V-north to eta-dir depth ave. Stokdes drift
              ustr2d(i,j)=wramp*ustr2d(i,j)
              vstr2d(i,j)=wramp*vstr2d(i,j)

              ustr0(i,j)=(cff1*wwus0g(i,j,it1)+cff2*wwus0g(i,j,it2))
     &                                         *cos(deg2rad*angler(i,j))
     &                +(cff1*wwvs0g(i,j,it1)+cff2*wwvs0g(i,j,it2))
     &                                         *sin(deg2rad*angler(i,j))! U-east to xi-dir depth ave. Stokes drift
              vstr0(i,j)=-(cff1*wwus0g(i,j,it1)+cff2*wwus0g(i,j,it2))
     &                                         *sin(deg2rad*angler(i,j))
     &                +(cff1*wwvs0g(i,j,it1)+cff2*wwvs0g(i,j,it2))
     &                                         *cos(deg2rad*angler(i,j))! V-north to eta-dir depth ave. Stokdes drift

              ustr0(i,j)=wramp*ustr0(i,j)
              vstr0(i,j)=wramp*vstr0(i,j)
# endif
#  endif
#  if defined WAVE_OFFLINE && defined SURFACE_ROLLER
            wvqb(i,j)=cff1*wwqb(i,j,it1)+cff2*wwqb(i,j,it2)
#  endif
#  if defined WAVE_OFFLINE && defined SUP_OFF
            sup(i,j)= cff1*wwsup(i,j,it1)+cff2*wwsup(i,j,it2)
#  endif
#  if defined WAVE_OFFLINE && defined WAVE_DIFF_FROM_LM
            llmm(i,j)= cff1*wwlm(i,j,it1)+cff2*wwlm(i,j,it2)
#  endif
# endif
          enddo
        enddo
#  if defined MRL_WCI
#  if defined WAVE_OFFLINE && !defined U2010
cc interpolate ust and ut2d on u,v points
          do j=jstrR,jendR
            do i=istr,iendR
                 ust2d(i,j)=0.5*(ustr2d(i-1,j)+ustr2d(i,j))
c                 ust0(i,j)=0.5*(ustr0(i-1,j)+ustr0(i,j))
        enddo
        enddo
          do j=jstr,jendR
            do i=istrR,iendR
                 vst2d(i,j)=0.5*(vstr2d(i,j-1)+vstr2d(i,j))
c                 vst0(i,j)=0.5*(vstr0(i,j-1)+vstr0(i,j))
        enddo
        enddo
# endif
# endif
      elseif (ZEROTH_TILE) then
        write(*,'(/1x,3A/3(1x,A,F16.10)/)')       '### WARNING: ',
     &  'set_wwave_tile :: Current model time is outside bounds of ',
     &  '''wwv_time''.',  'ww_tstart =', wwv_time(it1)*sec2day,
     &                    'tdays =',     tdays,
     &                    'ww_tend =',   wwv_time(it2)*sec2day
        ierr=ierr+1
      endif


      !*** END PART 1 ***!

!      !###### NOW CALCULATE WEC FORCES USING WHRM, etc ####
!      call wec_forces(istr,iend,jstr,jend)



      end subroutine set_wec_frc_tile

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
# ifdef SOLVE3D
     &                              A3d(1,1),A3d(1,2),A3d(1,3),A3d(1,4),
# endif
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
#  if defined WAVE_FRICTION && defined BODY_FRICTION
#   include "mixing.h"
#  endif
# endif
# ifdef ANA_BRY
#   include "boundary.h"
# endif
# ifdef WKB_WWAVE
#   include "wkb_wwave.h"
# endif
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
# ifdef STOKESEKMAN
     &            bconst=0.D0,
# else
     &            bconst=0.1D0, !0.03D0,               ! breaking contribution to KPP
# endif
     &            a_kv=1.2D0, !1.2D0                 ! breaking scale for eddy visc.
     &            a_brk=0.2D0,                ! breaking scale for body force
     &            a_frc=3.D0,                 ! friction scale for body force
     &            khmax=20.D0,                ! deep-water limit for k x dep
     &            c1o3=0.3333333333333333,    ! 1/3
     &            c4o3=1.3333333333333333     ! 4/3
     &           )
# ifdef WAVE_OFFLINE
      real roller_fraction
#  ifdef SURFACE_ROLLER
      parameter ( roller_fraction=1.0)
#  else
      parameter ( roller_fraction=0.0)
#  endif
#  define wkb_roller roller_fraction
# endif
#include "compute_auxiliary_bounds.h"

!      write(*,*) 'WEC_FORCES:::' ! DevinD debug

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

!      if (time/86400.0.ge.4352.0.and.time/86400.0.le.4352.0001) then
!        write(*,*)'hello'
!# ifdef FB_WSCALE
!        write(*,*)'kb=2k'
!        write(*,*)'a_brk=',a_brk
!        write(*,*)'a_kv=',a_kv
!        write(*,*)'c_b=',bconst
!# else
!        write(*,*)'kb=1/(a_brk*Hs)'
!        write(*,*)'a_brk=',a_brk
!        write(*,*)'a_kv=',a_kv
!        write(*,*)'c_b=',bconst
!# endif
!# ifdef BRK0
!        write(*,*)'BRK0'
!# endif
!# ifdef KVF0
!        write(*,*)'KVF0'
!# endif
!# ifdef AKB0
!        write(*,*)'AKB0'
!# endif
!      endif

! vertical distribution function F_B(z) for bottom streaming

c# define FF_FUNC1
c# define FF_FUNC2
# define FF_FUNC3

! initial ramping coefficients

# if defined WAVE_RAMP
#  ifdef STOKESEKMAN
      wave_ramp = tanh(3.0*dt*sec2day*float(iic-ntstart))
#  elif defined DUCK94
      wave_ramp = tanh(24.0*dt*sec2day*float(iic-ntstart))
#  elif defined MVCO
      wave_ramp = tanh(3.0*dt*sec2day*float(iic-ntstart))
!----- Cigdem (added RIP_CURRENT option)
#  elif defined RIP_CURRENT
c      wave_ramp = tanh(144.0*dt*sec2day*float(iic-ntstart)) ! 30 min
      wave_ramp = tanh(432.0*dt*sec2day*float(iic-ntstart))  ! 10 min
#  else
      wave_ramp = tanh(24.0*dt*sec2day*float(iic-ntstart))
#  endif
# else
      wave_ramp = 1.0
# endif
      wramp2= wave_ramp**2
!
! Evaluate 2DH wave-current intreraction variable.
! ================================================
!
! 2DH Stokes velocities, breaking, roller and bottom-friction
! dissipation terms, defined at horizontal rho-points.
iled
make: *** [wec_frc.o] Error 1
!
      inv_g = 1.0/g
      do j=jstr-1,jend+1
        do i=istr-1,iend+1
          Dstp(i,j)=max(h(i,j)+zeta(i,j,knew),eps)
          inv_d(i,j)=1.0/Dstp(i,j)
# ifdef WKB_WWAVE
          fr(i,j) =frq(i,j,wnew)
          wh(i,j) =wave_ramp*hrm(i,j,wnew)
          kw(i,j) =max(wvn(i,j,wnew),eps)
#  ifndef SURFACE_ROLLER
          cff = wramp2*wsb(i,j,wnew)
          ebrk(i,j)=cff*fr(i,j)                  ! epsilon_b (m3/s3)
!          ebrk(i,j)=0                  ! epsilon_b (m3/s3) ! kaida
#  else
          ebrk(i,j)=wramp2*wsb(i,j,wnew)*fr(i,j) ! epsilon_b (m3/s3)
!          ebrk(i,j)=0 ! epsilon_b (m3/s3) ! kaida
          erol(i,j)=wamp2*wsr(i,j,wnew)*fr(i,j) ! epsilon_r (m3/s3)
!          erol(i,j)=0 !kaida
          cff =wramp2*( wsr(i,j,wnew)+           ! for breaking term
     &           (1.-wkb_roller)*wsb(i,j,wnew) ) ! primary + roller
#  endif
          brk(i,j) =cff*kw(i,j)                  ! 2DH breaking term
          frc(i,j) =wramp2*wfc(i,j,wnew)*kw(i,j) ! 2DH friction term
#  ifdef MRL_CEW
          kD(i,j) =kw(i,j)*Dstp(i,j)
#  else
          kD(i,j) =kw(i,j)*(h(i,j)+wkb_tide)
#  endif
#  ifndef SURFACE_ROLLER
          act(i,j) =wramp2*wac(i,j,wnew)         ! wave action density
#  else
          act(i,j) =wramp2*(wac(i,j,wnew)+war(i,j,wnew))
#  endif
# else    /* ifndef WKB_WWAVE */
          fr(i,j) =wfrq(i,j)
          wh(i,j) =wave_ramp*whrm(i,j)
          inv_f(i,j) =1.0/max(fr(i,j),eps)
          khd = Dstp(i,j)*(fr(i,j)**2)*inv_g
#  ifdef KH_SOULSBY
          if (khd.ge.1.0D0) then        ! explicit wavenumber estimator
           kh=khd
          else
           kh=sqrt(khd)
          endif
          do k=1,3
            cff = tanh(kh)
            kh = kh-(kh*cff-khd)/max(cff+kh*(1.0-cff**2),eps)
          enddo
#  elif defined KH_HUNT
          kh = sqrt( khd*khd + khd/(1.0 + khd*(0.6666666666
     &               +khd*(0.3555555555 + khd*(0.1608465608
     &               +khd*(0.0632098765 + khd*(0.0217540484
     &                            +khd*0.0065407983)))))) )
#  endif
          kD(i,j) =kh
          kw(i,j) =kh*inv_d(i,j)
          ebrk(i,j) =wramp2*wveb(i,j)            ! epsilon_b (m3/s3)
!          ebrk(i,j) =0            ! epsilon_b (m3/s3) kaida

#  ifdef SURFACE_ROLLER
          erol(i,j) =wramp2*0.5762*Dstp(i,j)*    ! epsilon_r (m3/s3)
     &               whrm(i,j)*wvqb(i,j)*        ! 0.5762=0.06*sinb*g^2
     &               kw(i,j)*inv_f(i,j)          ! Nairn et al (1991)
          cff =(1.-wkb_roller)*ebrk(i,j)+erol(i,j)
#  else
          cff =ebrk(i,j)
#  endif
          frc(i,j) =wramp2*wved(i,j)*kw(i,j)
     &                           *inv_f(i,j)     ! 2DH friction term
          brk(i,j) =cff*kw(i,j)*inv_f(i,j)       ! 2DH breaking term
          act(i,j) =inv_f(i,j)*(
     &                      0.125*g*(wh(i,j)**2) ! wave action density
#  ifdef SURFACE_ROLLER
     &          +wramp2*0.03*g*Dstp(i,j)*wh(i,j)*wvqb(i,j)
#  endif
     &                                                   )
# endif   /* ifdef WKB_WWAVE */
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

!      ! DevinD debug
!      if (mynode==5) then
!        write(*,*) 'tile-5: act(14,39)', act(14,39) ! DevinD debug
!      endif

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
# if defined WAVE_OFF_WTIDES
          supm=-0.125*(wh(i,j)**2)*kw(i,j)
     &              /max(sinh(2.*min(kD(i,j),khmax)),eps)
          if (kD(i,j).lt.6.and.supm.ne.0.) then
          Fsup=min(1.,sqrt(sup(i,j)/supm))
          else
          Fsup=1.
          endif
# else
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
# endif
          actf(i,j)=act(i,j)*Fsup**2
        enddo
      enddo
# endif

!      ! DevinD debug
!      if (mynode==5) then
!        write(*,*) 'tile-5: actf(14,39)', actf(14,39) ! DevinD debug
!        write(*,*) 'tile-5: Fsup', Fsup ! DevinD debug
!      endif

!
! Dissipation terms : copy into shared arrays
! ===========================================
!
      do j=jstrR,jendR
        do i=istrR,iendR
          wdsp(i,j)=ebrk(i,j)         ! epsilon_b (m3/s3)
# ifdef SURFACE_ROLLER
          rdsp(i,j)=erol(i,j)         ! epsilon_r (m3/s3)
# endif
# ifdef WKB_WWAVE
          cff = wramp2*wfc(i,j,wnew)
# else
!# elif !defined WKB_WWAVE & defined WAVE_OFFLINE
          cff = wramp2*wved(i,j)*inv_f(i,j)
# endif
          wdrg(i,j)=cff*fr(i,j)       ! epsilon_d (m3/s3)
          wbst(i,j)=cff*kw(i,j)       ! e_d stress (m2/s2) for KPP
# ifdef MASKING
          wdsp(i,j)=wdsp(i,j)*rmask(i,j)
          wbst(i,j)=wbst(i,j)*rmask(i,j)
          wdrg(i,j)=wdrg(i,j)*rmask(i,j)
#  ifdef SURFACE_ROLLER
          rdsp(i,j)=rdsp(i,j)*rmask(i,j)
#  endif
# endif
        enddo
      enddo
# ifndef WKB_WWAVE
#  undef inv_f
# endif
!
! 2DH depth-averaged Stokes drift and breaking acceleration
! =========================================================
!
      do j=jstrR,jendR
        do i=istr,iendR
# if defined WKB_WWAVE || defined U2010
          ust2d(i,j)  =0.5*(stk(i-1,j)*wdrx(i-1,j)+stk(i,j)*wdrx(i,j))
# endif
!          brk2dx(i,j) =0.5*(brk(i-1,j)*wdrx(i-1,j)+brk(i,j)*wdrx(i,j))
          brk2dx(i,j) =0.
!          frc2dx(i,j) =0.5*(frc(i-1,j)*wdrx(i-1,j)+frc(i,j)*wdrx(i,j))
          frc2dx(i,j) =0.
#  ifdef MASKING
#   ifdef SHOREFACE
          ust2d(i,j)  = ust2d(i,j)*umask(i,j)
#   endif
          brk2dx(i,j) =brk2dx(i,j)*umask(i,j)
          frc2dx(i,j) =frc2dx(i,j)*umask(i,j)
#  endif
        enddo
      enddo
      do j=jstr,jendR
        do i=istrR,iendR
# if defined WKB_WWAVE || defined U2010
          vst2d(i,j)  =0.5*(stk(i,j-1)*wdre(i,j-1)+stk(i,j)*wdre(i,j))
# endif
!          brk2de(i,j) =0.5*(brk(i,j-1)*wdre(i,j-1)+brk(i,j)*wdre(i,j))
          brk2de(i,j) =0.
!          frc2de(i,j) =0.5*(frc(i,j-1)*wdre(i,j-1)+frc(i,j)*wdre(i,j))
          frc2de(i,j) =0.
#  ifdef MASKING
#   ifdef SHOREFACE
          vst2d(i,j)  = vst2d(i,j)*vmask(i,j)
#   endif
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
#  if defined SOLVE3D && !defined BBL_QUAD2D
          cff1 = 0.5*(u(i,j,1,nrhs)+u(i+1,j,1,nrhs))
          cff2 = 0.5*(v(i,j,1,nrhs)+v(i,j+1,1,nrhs))
#  else
          cff1 = 0.5*(ubar(i,j,knew)+ubar(i+1,j,knew))
          cff2 = 0.5*(vbar(i,j,knew)+vbar(i,j+1,knew))
#  endif
          umag = sqrt( cff1*cff1+cff2*cff2 )
# if defined BBL_F00 || defined BBL_S95
          Uorb1 = 0.5*fr(i,j)*wh(i,j)
     &              /max(sinh(min(kD(i,j),khmax)),eps)
#  if defined BBL_F00
          ka_f00  =0.0125   ! apparent roughness (Ruessink et al., 2001)
          r_D(i,j)=0.015*((ka_f00*inv_d(i,j))**c1o3)
     &               *sqrt(1.3456*(Uorb1**2)+umag**2)
#  else
#   ifdef SOLVE3D
          cff1=umag*(vonKar/log(1.+0.5*Hz(i,j,1)*inv_zb))**2
c**          hz0 =max(z_r(i,j,1)-z_w(i,j,0),1.1*Zob)
#   else
          cff1=umag*(vonKar/log(1.+Dstp(i,j)*inv_zb))**2
c**          hz0 =max(0.5*Dstp(i,j),1.1*Zob)
#   endif
c**          cff1=umag*(vonKar/log(hz0*inv_zb))**2
          tauc=cff1*umag
          tauw=0.695*(Uorb1**1.48)*((Zob*fr(i,j))**0.52)
          cff2=1.0 + 1.2*((tauw/max(tauw+tauc,eps))**3.2)
          r_D(i,j)=cff1*cff2
#  endif   /* BBL_F00 */
# elif defined BBL_QUAD2D
          r_D(i,j)=rdrg2*umag
# endif    /* BBL_F00 || BBL_S95 */
# ifdef MASKING
          r_D(i,j)=r_D(i,j)*rmask(i,j)
# endif
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
!           inv_dz(i,j,k)=1./Hz(i,j,k)
            cff2 = exp( 2.*kw(i,j)*(z_w(i,j,k)+h(i,j)-Dstp(i,j)))
     &            -exp(-2.*kw(i,j)*(z_w(i,j,k)+h(i,j)+Dstp(i,j)))
!           ust_rm(i,j,k) =ust_ker*inv_dz(i,j,k)*( cff2-cff1 )
            ust_rm(i,j,k) =ust_ker*( cff2-cff1 )/Hz(i,j,k)
            if (mynode==1) then
              if (i==10.and.j==10.and.k<3) then
!               print *,'cff1/2:  ',cff2-cff1,inv_dz(i,j,k)
!               print *,'ust_rm:  ',ust_rm(i,j,k)
              endif
            endif
            cff1=cff2
          enddo
         enddo
       enddo                 ! <-- keep ust_r, inv_dz, &inv_ex
       if (mynode==1) then
c        if (i==10.and.j==10) then
           write(*,*) 'i,j:      ',i,j
           write(*,*) 'ust_ker:  ',ust_ker
           write(*,*) 'ust_rm(1): ',ust_rm(10,10,1)
           write(*,*) 'ust_rm(2): ',ust_rm(10,10,2)
           write(*,*) 'ust_rm(N): ',ust_rm(10,10,N)
c        endif
       endif
       stop

!
# if defined WAVE_OFFLINE && !defined U2010
!
! keff
c
       do j=jstr-1,jend+1
         do i=istr-1,iend+1

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

c             vst_r(i,j,:)= ust_rm(i,j,:)*wdre(i,j) !
c             ust_r(i,j,:)= ust_rm(i,j,:)*wdrx(i,j)  !
c             if (abs(1-FAC).gt.0.02) then
c             write(*,*) 'FAC=', FAC, 'count=',COUNTR,'i=',i,'j=',j
c!             write(*,*) 'FAC=', FAC, 'count=',COUNTR,'h=',h(i,j)
c             endif
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
# if !defined WAVE_OFFLINE || defined U2010
            ust(i,j,k)=0.5*( ust_rm(i-1,j,k)*wdrx(i-1,j)
     &                        +ust_rm(i,j,k)*wdrx(i,j) )
# else
            ust(i,j,k)=0.5*( ust_r(i-1,j,k)
     &                        +ust_r(i,j,k))/fac2d(i,j)
# endif
          enddo
        enddo
        do j=jstr,jendR
          do i=istrR,iendR
# if !defined WAVE_OFFLINE || defined U2010
            vst(i,j,k)=0.5*( ust_rm(i,j-1,k)*wdre(i,j-1)
     &                        +ust_rm(i,j,k)*wdre(i,j) )
# else
            vst(i,j,k)=0.5*( vst_r(i,j-1,k)
     &                        +vst_r(i,j,k) )/fac2d(i,j)
# endif
          enddo
        enddo
      enddo


!      ! DevinD debug
!     if (mynode==1) then
!        write(*,*) 'mynode=', mynode, 'ust(47,45,N)=', ust(47,45,N)
!        write(*,*) 'mynode=', mynode, 'ust(48,45,N)=', ust(48,45,N)
!        write(*,*) 'mynode=', mynode, 'ust(49,45,N)=', ust(49,45,N)
!        write(*,*) 'wdrx(50,45)=',wdrx(50,45)
!        write(*,*) 'ust_rm(50,45,N)=',ust_rm(50,45,N)
!        write(*,*) 'ust_r(49,45,N)=',ust_r(49,45,N)
!        write(*,*) 'ust_r(50,45,N)=',ust_r(50,45,N),'fac2d(50,45)=',fac2d(50,45)
!        write(*,*) 'mynode=', mynode, 'ust(50,45,N)=', ust(50,45,N)
!        write(*,*) 'mynode=', mynode, 'ust(51,45,N)=', ust(51,45,N)
!      endif
!      write(*,*) 'mynode=', mynode
!      write(*,*) 'ust(10,10,N)=', ust(10,10,N)
!!      write(*,*) 'ust_rm(9,10,N)=', ust_rm(9,10,N)
!!      write(*,*) 'wdre(9,10)=', wdre(9,10)
!!      write(*,*) 'ust_rm(10,10,N)=', ust_rm(10,10,N)
!!      write(*,*) 'wdre(10,10)=', wdre(10,10)
!      write(*,*) 'ustr0(10,10)=', ustr0(10,10)
!      write(*,*) 'vst(10,10,N)=', vst(10,10,N)
!      write(*,*) 'vstr0(10,10)=', vstr0(10,10)

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
# if !defined WAVE_OFFLINE || defined U2010
          kvsurf =1.5*kv(N)-0.5*kv(N-1) ! extrapolate to surface
# else
          kvsurf = 0.5*(ust_r(i,j,N)*( u(i,j,N,nstp)+u(i+1,j,N,nstp))
     &              + vst_r(i,j,N)*(v(i,j,N,nstp)+v(i,j+1,N,nstp) ))
# endif
          do k=1,N-1
            dkvdz(k) =2.0*(kv(k+1)-kv(k))/(Hz(i,j,k+1)+Hz(i,j,k))
          enddo
          dkvdz(0) = dkvdz(1)    !2.*dkvdz(1)-dkvdz(2) ! severe!
          dkvdz(N) = dkvdz(N-1)  !2.*dkvdz(N-1)-dkvdz(N-2)
!      ! DevinD debug
!      if (mynode==5 .and. i==14 .and. j==39) then
!        write(*,*) 'tile-5: Hz(14,39,49)', Hz(14,39,49) ! DevinD debug
!        write(*,*) 'tile-5: kv(48)', kv(48) ! DevinD debug
!        write(*,*) 'tile-5: kv(49)', kv(49) ! DevinD debug
!        write(*,*) 'tile-5: dkvdz(49)', dkvdz(49) ! DevinD debug
!      endif
          do k=1,N
            d2kv(k) =dkvdz(k)-dkvdz(k-1) ! d^2kv/dz^2 x Hz
          enddo
!      ! DevinD debug
!      if (mynode==5 .and. i==14 .and. j==39) then
!        write(*,*) 'tile-5: d2kv(49)', d2kv(49) ! DevinD debug
!      endif
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
# if !defined WAVE_OFFLINE || defined U2010
          cff4 =-2.0*kw(i,j)*kvsurf
# else
          cff4 = -2.0*(sinh(min(khmax,keff(i,j)*Dstp(i,j))))**2/
     &                 cosh(min(2.*khmax,2.*keff(i,j)*h(i,j)))*kvsurf
# endif
! min limits to avoid blowups
!ccccc          actp =0.125*g*(wh(i,j)**2)/max(fr(i,j),eps)
!ccccc          calP(i,j) = actp*inv_g*tanh(kD(i,j))
# if defined SUP_OFF
          actp=actf(i,j)
# else
          actp=act(i,j)
# endif
# if !defined WAVE_OFFLINE || defined U2010
          calP(i,j) = actp*tanh(kD(i,j))
     &                    *( cff1+cff2+cff3+cff4 )
# else
          calP(i,j) = (actp*tanh(kD(i,j))
     &                    *( cff1+cff2+cff3)+ cff4) ! cfff4 now computed from surface Ustk with u
!          calP(i,j) = 0. ! DevinD zeros
# endif
#  ifdef MASKING
     &                                 *rmask(i,j)
#  endif
          cff=0.0
          do k=1,N
            cff = cff +  d2kv(k)*
     &         ( exp( 2.*kw(i,j)*(z_w(i,j,N)-z_r(i,j,k)-Dstp(i,j)))
     &          -exp(-2.*kw(i,j)*(z_w(i,j,N)-z_r(i,j,k)+Dstp(i,j))) )
          enddo
!      ! DevinD debug
!      if (mynode==5 .and. i==14 .and. j==39) then
!        write(*,*) 'tile-5: actp', actp ! DevinD debug
!        write(*,*) 'tile-5: d2kv(48)', d2kv(48) ! DevinD debug
!        write(*,*) 'tile-5: d2kv(49)', d2kv(49) ! DevinD debug
!        write(*,*) 'tile-5: kw(14,39)', kw(14,39) ! DevinD debug
!        write(*,*) 'tile-5: z_w(14,39,N)', z_w(14,39,N) ! DevinD debug
!        write(*,*) 'tile-5: z_r(14,39,49)', z_r(14,39,49) ! DevinD debug
!        write(*,*) 'tile-5: Dstp(14,39)', Dstp(14,39) ! DevinD debug
!        write(*,*) 'tile-5: cff', cff ! DevinD debug
!      endif
#  undef KAPSRF_SAFE
#  ifdef KAPSRF_SAFE
c#   define FLT 0.98
#   define FLT 0.5
          Kapsrf(i,j) = ( FLT*Kapsrf(i,j)  ! time-filter to suppress noises
     &                   +(1.0-FLT)*cff*actp*inv_ex(i,j) )
#  else
          Kapsrf(i,j) = cff*actp*inv_ex(i,j)
!          Kapsrf(i,j) = 0. ! DevinD zeros
#  endif
#  ifdef MASKING
     &                           *rmask(i,j)
#  endif
        enddo
      enddo      ! <-- discard kv, dkvdz, d2kv, & inv_ex

!      ! DevinD debug
!      if (mynode==5) then
!        write(*,*) 'tile-5: inv_ex(14,39)', inv_ex(14,39) ! DevinD debug
!        write(*,*) 'tile-5: Kapsrf(14,39)', Kapsrf(14,39) ! DevinD debug
!      endif

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
# if defined WAVE_OFFLINE && !defined U2010
            kvr(k) = 0.25*(
     &                  (ust_r(i,j,k)+ust_r(i,j,k+1))*
     &                (u(i,j,k+1,nstp)-u(i,j,k,nstp)
     &                        +u(i+1,j,k+1,nstp)-u(i+1,j,k,nstp))
     &                 + (vst_r(i,j,k)+vst_r(i,j,k+1))*
     &                (v(i,j,k+1,nstp)-v(i,j,k,nstp)
     &                        +v(i,j+1,k+1,nstp)-v(i,j+1,k,nstp)))
     &                                 /(z_r(i,j,k+1)-z_r(i,j,k))
# else
            kvr(k) = 0.25*(ust_rm(i,j,k)+ust_rm(i,j,k+1))*
     &                ( wdrx(i,j)*(u(i,j,k+1,nstp)-u(i,j,k,nstp)
     &                        +u(i+1,j,k+1,nstp)-u(i+1,j,k,nstp))
     &                 +wdre(i,j)*(v(i,j,k+1,nstp)-v(i,j,k,nstp)
     &                        +v(i,j+1,k+1,nstp)-v(i,j+1,k,nstp)))
     &                                 /(z_r(i,j,k+1)-z_r(i,j,k))
# endif
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

# ifdef KVF0
      do k=1,N
        do j=jmin,jend
          do i=imin,iend
            kvf(i,j,k)=0.D0
          enddo
        enddo
      enddo
# endif


!
! Akw: primary wave-induced additional diffusivity for tracers at w-point
! =======================================================================
!
      cff=wave_ramp/dt
      do j=jstrR,jendR
        do i=istrR,iendR
# if !defined WAVE_DIFF_FROM_LM
          cff2 = 0.25*wh(i,j)/max(1.-exp(-2.*kD(i,j)),eps)
          do k=0,N
            cff1 = cff2*( exp( kw(i,j)*(z_w(i,j,k)+h(i,j)-Dstp(i,j)))
     &                   -exp(-kw(i,j)*(z_w(i,j,k)+h(i,j)+Dstp(i,j))) )
            cff3 = cff1**2
            Akw(i,j,k) = max(0.,cff*(cff3-E_pre(i,j,k))) ! avoid negatives
#  ifdef MASKING
     &                               *rmask(i,j)
#  endif
            E_pre(i,j,k) = cff3   ! store this for the next time step.
          enddo
# else
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
# endif
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
#  ifdef FKV_FUNC1
          fn1   =3.*(exp(-2.*kbrk*Dstp(i,j))+exp(-4.*kbrk*Dstp(i,j)))
          fn2   =exp(-6.*kbrk*Dstp(i,j))
          fb1   =c4o3*inv_k*(fn1+2.*fn2)/(1.+fn1+fn2)
          intfb =c4o3*inv_k - fb1  ! intfb = \int^D gb(z') dz'
          fb0   =1.D0/intfb        ! 1/[\int^D gb(z') dz']
#  elif defined FKV_FUNC2
          fn1   =exp(-2.*kbrk*Dstp(i,j))
          fb1   =2.*inv_k*fn1/(1.+fn1)
          intfb =inv_k - fb1
          fb0   =1.D0/intfb
#  elif defined FKV_FUNC3
          fb1   =0.0               ! fb1 = \int gb dz at z = -h
          fb0   =1.D0/(1.-exp(-2.*kbrk*Dstp(i,j)))
#  endif
#  ifndef SURFACE_ROLLER
          vs    =ebrk(i,j)**c1o3   ! representative velocity scale
#  else
          vs    =((1.-wkb_roller)*ebrk(i,j)+erol(i,j))**c1o3
#  endif
          kvb   =bconst*vs*wh(i,j) ! depth-averaged Akb
          do k=1,N,+1              ! <-- irreversible
#  ifdef FKV_FUNC1
            dd  =zeta(i,j,knew)-z_w(i,j,k)
            fn1 =3.*(exp(-2.*kbrk*dd)+exp(-4.*kbrk*dd))
            fn2 =exp(-6.*kbrk*dd)
            fb2 =c4o3*inv_k*(fn1+2.*fn2)/(1.+fn1+fn2)
#  elif defined FKV_FUNC2
            fn1 =exp(-2.*kbrk*(zeta(i,j,knew)-z_w(i,j,k)))
            fb2 =2.*inv_k*fn1/(1.+fn1)
#  elif defined FKV_FUNC3
            fb2 =exp( kbrk*(z_w(i,j,k)+h(i,j)-Dstp(i,j)))
     &          -exp(-kbrk*(z_w(i,j,k)+h(i,j)+Dstp(i,j)))
#  endif
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

# ifndef SURFACE_BREAK
#  define brk_r wrk3
!
! 3D breaking acceleration term defined as a body force.
! =====================================================
!
      do j=jstr-1,jend+1
        do i=istr-1,iend+1
#  ifdef FB_FUNC0
          do k=1,N                 ! vertically uniform case
            brk_r(i,j,k) = brk(i,j)*inv_d(i,j)
          enddo                    ! <-- discard inv_d
#  else
#   ifdef FB_WSCALE
          inv_k =1.0/max(2.0*kw(i,j),eps)  ! Stokes scale
#   else
          inv_k =min(max(a_brk*wh(i,j),eps),Dstp(i,j))
#   endif
          kbrk  =1.0/inv_k         ! vertical scale for breaking
#   ifdef FB_FUNC1
          fn1   =3.*(exp(-2.*kbrk*Dstp(i,j))+exp(-4.*kbrk*Dstp(i,j)))
          fn2   =exp(-6.*kbrk*Dstp(i,j))
          fb1   =c4o3*inv_k*(fn1+2.*fn2)/(1.+fn1+fn2)
          intfb =c4o3*inv_k - fb1  ! intfb = \int^D gb(z') dz'
          fb0   =1.D0/intfb        ! 1/[\int^D gb(z') dz']
#   elif defined FB_FUNC2
          fn1   =exp(-2.*kbrk*Dstp(i,j))
          fb1   =2.*inv_k*fn1/(1.+fn1)
          intfb =inv_k - fb1
          fb0   =1.D0/intfb
#   elif defined FB_FUNC3
          fb1   =0.0               ! fb1 = \int gb dz at z = -h
          fb0   =1.D0/(1.-exp(-2.*kbrk*Dstp(i,j)))
#   endif
          do k=1,N,+1              ! <-- irreversible
#   ifdef FB_FUNC1
            dd  =zeta(i,j,knew)-z_w(i,j,k)
            fn1 =3.*(exp(-2.*kbrk*dd)+exp(-4.*kbrk*dd))
            fn2 =exp(-6.*kbrk*dd)
            fb2 =c4o3*inv_k*(fn1+2.*fn2)/(1.+fn1+fn2)
#   elif defined FB_FUNC2
            fn1 =exp(-2.*kbrk*(zeta(i,j,knew)-z_w(i,j,k)))
            fb2 =2.*inv_k*fn1/(1.+fn1)
#   elif defined FB_FUNC3
            fb2 =exp( kbrk*(z_w(i,j,k)+h(i,j)-Dstp(i,j)))
     &          -exp(-kbrk*(z_w(i,j,k)+h(i,j)+Dstp(i,j)))
#   endif
            fb =fb0*inv_dz(i,j,k)*(fb2-fb1)
            brk_r(i,j,k) = fb*brk(i,j)
            fb1=fb2     ! recursive procedure
          enddo
#  endif       /* ifdef FB_FUNC0 */
        enddo
      enddo

      do k=1,N
        do j=jstrR,jendR
          do i=istr,iendR
            brk3dx(i,j,k)=0.5*( brk_r(i-1,j,k)*wdrx(i-1,j)
     &                           +brk_r(i,j,k)*wdrx(i,j) )
!            brk3dx(i,j,k)=0. ! DevinD zeros
#  ifdef MASKING
     &                                         *umask(i,j)
#  endif
          enddo
        enddo
        do j=jstr,jendR
          do i=istrR,iendR
            brk3de(i,j,k)=0.5*( brk_r(i,j-1,k)*wdre(i,j-1)
     &                           +brk_r(i,j,k)*wdre(i,j) )
!            brk3de(i,j,k)=0. ! DevinD
#  ifdef MASKING
     &                                         *vmask(i,j)
#  endif
          enddo
        enddo
      enddo       ! <-- discard brk_r, but keep ust_r
#  undef brk_r
# endif        /* ifndef SURFACE_BREAK */
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
# if defined WAVE_FRICTION && defined BODY_FRICTION
#  undef LH53_STREAMING
#  ifdef ANA_VMIX
#   define LH53_STREAMING
#  endif
#  define frc_r wrk3
#  define frc2d stk
      do j=jstr-1,jend+1
        do i=istr-1,iend+1
#  ifndef LH53_STREAMING
          abot  =wh(i,j)/max(2.D0*sinh(min(kD(i,j),khmax)),eps)
          delta =0.160455*((abot/Zob)**0.83)*Zob  ! BBL thickness
          inv_k =min(max(a_frc*delta,eps),Dstp(i,j))
          kfrc  =1.0/inv_k         ! vertical scale for wave friction
#   ifdef FF_FUNC1
          fn1   =3.*(exp(-2.*kfrc*Dstp(i,j))+exp(-4.*kfrc*Dstp(i,j)))
          fn2   =exp(-6.*kfrc*Dstp(i,j))
          fb1   =c4o3*inv_k        ! fb1: at z=-h (z_w at k=0)
          fb2   =c4o3*inv_k*(fn1+2.*fn2)/(1.+fn1+fn2)
          intfb =fb2-fb1           ! intfb = \int^D gb(z') dz'
          fb0   =1.D0/intfb        ! 1/[\int^D gb(z') dz']
#   elif defined FF_FUNC2
          fn1   =exp(-2.*kfrc*Dstp(i,j))
          fb1   =inv_k
          fb2   =2.*inv_k*fn1/(1.+fn1)
          intfb =fb2-fb1
          fb0   =1.D0/intfb
#   elif defined FF_FUNC3
          fb1   =1.0-exp(-2.*kfrc*Dstp(i,j))
          fb0   =1.D0/(exp(-2.*kfrc*Dstp(i,j))-1.0)
#   endif
          do k=1,N,+1              ! <-- irreversible
#   ifdef FF_FUNC1
            dd  =h(i,j)+z_w(i,j,k)
            fn1 =3.*(exp(-2.*kfrc*dd)+exp(-4.*kfrc*dd))
            fn2 =exp(-6.*kfrc*dd)
            fb2 =c4o3*inv_k*(fn1+2.*fn2)/(1.+fn1+fn2)
#   elif defined FF_FUNC2
            fn1 =exp(-2.*kfrc*(h(i,j)+z_w(i,j,k)))
            fb2 =2.*inv_k*fn1/(1.+fn1)
#   elif defined FF_FUNC3
            fb2 =exp( kfrc*(zeta(i,j,knew)-z_w(i,j,k)-Dstp(i,j)))
     &          -exp(-kfrc*(zeta(i,j,knew)-z_w(i,j,k)+Dstp(i,j)))
#   endif
            fb =fb0*inv_dz(i,j,k)*(fb2-fb1)
            frc_r(i,j,k) = fb*frc(i,j)
            fb1=fb2
          enddo
#  else  /* LH53_STREAMING */
          frc2d(i,j)=0.0
          beta=sqrt(fr(i,j)/(Akv(i,j,0)+Akv(i,j,1)))
          cff =wramp2*(wh(i,j)**2)*(fr(i,j)**2)*kw(i,j)
     &         /max(16.0*(sinh(min(kD(i,j),khmax))**2)*beta,eps)
          bz1 =beta*(h(i,j)+z_w(i,j,0))
          cff1=2.*exp(-bz1)*(bz1*sin(bz1)+cos(bz1))-1.-exp(-2.*bz1)
          cff3=1.D0 !min(1.0, pi/beta/Hz(i,j,1)) ! <-- unresoluved viscous layer
          do k=1,N,+1              ! <-- irreversible
            bz2 =beta*(h(i,j)+z_w(i,j,k))
            cff2=2.*exp(-bz2)*(bz2*sin(bz2)+cos(bz2))-1.-exp(-2.*bz2)
            frc_r(i,j,k)=-cff*(cff2-cff1)*inv_dz(i,j,k)*cff3
            cff1=cff2; cff3=1.0
            frc2d(i,j) = frc2d(i,j) + frc_r(i,j,k)*Hz(i,j,k)
          enddo
#  endif  /* !LH53_STREAMING */
        enddo
      enddo
      do j=jstrR,jendR
        do i=istr,iendR
#  ifdef LH53_STREAMING
          frc2dx(i,j)=0.5*( frc2d(i-1,j)*wdrx(i-1,j)
     &                       +frc2d(i,j)*wdrx(i,j) )
#   ifdef MASKING
     &                                   *umask(i,j)
#   endif
#  endif
          do k=1,N
            frc3dx(i,j,k)=0.5*( frc_r(i-1,j,k)*wdrx(i-1,j)
     &                           +frc_r(i,j,k)*wdrx(i,j) )
#  ifdef MASKING
     &                                         *umask(i,j)
#  endif
          enddo
        enddo
      enddo
      do j=jstr,jendR
        do i=istrR,iendR
#  ifdef LH53_STREAMING
          frc2de(i,j)=0.5*( frc2d(i,j-1)*wdre(i,j-1)
     &                       +frc2d(i,j)*wdre(i,j) )
#   ifdef MASKING
     &                                   *vmask(i,j)
#   endif
#  endif
          do k=1,N
            frc3de(i,j,k)=0.5*( frc_r(i,j-1,k)*wdre(i,j-1)
     &                           +frc_r(i,j,k)*wdre(i,j) )
#  ifdef MASKING
     &                                         *vmask(i,j)
#  endif
          enddo
        enddo
      enddo       ! <-- discard frc_r, frc2d
#  undef frc_r
#  undef frc2d
# endif   /* BODY_FRICTION */
# undef inv_dz

# endif       /* SOLVE3D */

! Open boundary condition (valid only for idealized situations)
! =============================================================

#   ifdef ANA_BRY
      z_tide=0.0
#    ifdef WKB_WWAVE
      z_tide=wkb_tide
#    endif

#    if defined Z_FRC_BRY || defined M2_FRC_BRY || defined M3_FRC_BRY

#     ifdef OBC_WEST
      if (WESTERN_EDGE) then
        do j=jstrR,jendR
#      ifdef Z_FRC_BRY
          zeta_west(j)=z_tide + sup(istrR,j)
#      endif
#      ifdef M2_FRC_BRY
#       ifdef SHOREFACE
          ubar_west(j)=ubar(istr,j,knew)
#       else
          ubar_west(j)=-ust2d(istrR,j)    ! Stokes drift (xi)
#       endif
          vbar_west(j)=vbar(istr,j,knew) ! Stokes drift (eta)
#      endif
#      if defined SOLVE3D && defined M3_FRC_BRY
          do k=1,N
            u_west(j,k)=u(istr,j,k,nnew)
c            u_west(j,k)=-ust(istrR,j,k)
            v_west(j,k)=v(istr,j,k,nnew)
          enddo
#      endif
        enddo
      endif
#     endif  /* OBC_WEST */
!
#     ifdef OBC_EAST
      if (EASTERN_EDGE) then
        do j=jstrR,jendR
#      ifdef Z_FRC_BRY
          zeta_east(j)=z_tide + sup(iendR,j)
#      endif
#      ifdef M2_FRC_BRY
c#      ifdef STOKESEKMAN
c          ubar_east(j)=ubar(iend,j,knew) ! Neumann
c#      else
          ubar_east(j)=-ust2d(iend,j)  ! Stokes drift (xi)
c#      endif
          vbar_east(j)=vbar(iend,j,knew) !-vst2d(iendR,j)  ! Stokes drift (eta)
#      endif
#      if defined SOLVE3D && defined M3_FRC_BRY
          do k=1,N
            u_east(j,k)=u(iend,j,k,nnew)
c            u_east(j,k)=-ust(iendR,j,k)
            v_east(j,k)=v(iend,j,k,nnew)
          enddo
#      endif
        enddo
      endif
#     endif  /* OBC_EAST */
!
#     ifdef OBC_SOUTH
      if (SOUTHERN_EDGE) then
        do i=istrR,iendR
#      ifdef Z_FRC_BRY
          zeta_south(i)=z_tide + sup(i,jstrR)
#      endif
#      ifdef M2_FRC_BRY
          ubar_south(i)=ubar(i,jstr,knew)  ! Neumann
          vbar_south(i)=-vst2d(i,jstr)     ! Stokes drift
#      endif
#      if defined SOLVE3D && defined M3_FRC_BRY
          do k=1,N
            u_south(i,k)=u(i,jstr,k,nnew)
            v_south(i,k)=v(i,jstr,k,nnew)
          enddo
#      endif
        enddo
      endif
#     endif  /* OBC_SOUTH */
!
#     ifdef OBC_NORTH
      if (NORTHERN_EDGE) then
        do i=istrR,iendR
#      ifdef Z_FRC_BRY
          zeta_north(i)=z_tide + sup(i,jendR)
#      endif
#      ifdef M2_FRC_BRY
          ubar_north(i)=ubar(i,jend,knew)  ! Neumann
          vbar_north(i)=-vst2d(i,jend)     ! Stokes drift
#      endif
#      if defined SOLVE3D && defined M3_FRC_BRY
          do k=1,N
            u_north(i,k)=u(i,jend,k,nnew)
            v_north(i,k)=v(i,jend,k,nnew)
          enddo
#      endif
        enddo
      endif
#     endif  /* OBC_NORTH */
!
#    endif  /* Z_FRC_BRY || M2_FRC_BRY || M3_FRC_BRY */
!
#   endif  /* ANA_BRY */


# ifdef EXCHANGE
      call exchange2d_4_tile(istr,iend,jstr,jend, sup,wdsp,wdrg,wbst)
      call exchange2d_3_tile(istr,iend,jstr,jend, ust2d,brk2dx,frc2dx)
      call exchange2d_3_tile(istr,iend,jstr,jend, vst2d,brk2de,frc2de)
#  ifdef SURFACE_ROLLER
      call exchange2d_tile(istr,iend,jstr,jend,rdsp)
#  endif
#  if defined BBL_S95 || defined BBL_F00 || defined BBL_QUAD2D
      call exchange2d_tile (istr,iend,jstr,jend, r_D)
#  endif
#  ifdef SOLVE3D
      call exchange2d_2_tile(istr,iend,jstr,jend,Kapsrf,calP)
      call exchange_3_tile(istr,iend,jstr,jend, ust,N, vst,N, kvf,N)
#   ifndef SURFACE_BREAK
      call exchange_2_tile (istr,iend,jstr,jend, brk3dx,N, brk3de,N)
#   endif
#   if defined WAVE_FRICTION && defined BODY_FRICTION
      call exchange_2_tile (istr,iend,jstr,jend, frc3dx,N, frc3de,N)
#   endif
#  endif
# endif


      end subroutine wec_forces_tile


      ! -------------------------------------------
      ! -------------------------------------------
      ! -------------------------------------------


!        subroutine init_wec_arrays ! DevinD added this routine.
!
!        ! Initiate arrays used by WEC (from DH's MRL init_arrays.F)
!        !      note: JM's river_frc module doesn't have an init step, it is done in main loop.
!
!            implicit none
!
!            ! beware ANA_SMFLUX and ANA_SRFLUX entry in init_arrays.F
!
!
!
!        end subroutine init_wec_arrays


#endif /* WEC for whole module */
        
      end module wec_frc
