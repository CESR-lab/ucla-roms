      module wec_frc
      ! Wave Effect on Current (WEC) module which provides all the terms that enter
      ! in the momentum equations

       use dimensions, nt_dimensions => nt, mynode_dims => mynode !! has dimensions, a list of forcing files and rmask

       implicit none

! DevinD - need the following includes for GLOBAL_2D_ARRAY to work
#include "cppdefs.h"
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
       real, public, dimension(GLOBAL_2D_ARRAY) :: brk2dx, brk2de !! xi/eta-direction 2D breaking dissipation (rho)
       real, public, dimension(GLOBAL_2D_ARRAY) :: frc2dx,frc2de  !! xi/eta-direction 2D frictional dissipation (rho)
       real, public, dimension(GLOBAL_2D_ARRAY) :: ust2d     !! xi-direction 2D Stokes transport (u-point)
       real, public, dimension(GLOBAL_2D_ARRAY) :: vst2d     !! eta-direction 2D Stokes transport (v-point)
       real, public, dimension(GLOBAL_2D_ARRAY) :: sup            !! quasi-static wave set-up (rho-point)
       real, public, dimension(GLOBAL_2D_ARRAY) :: calP           !! pressure correction term (rho-point)
       real, public, dimension(GLOBAL_2D_ARRAY) :: Kapsrf         !! Bernoulli head term at the surface (rho-point)
!      real, public, dimension(nx,ny) :: lm             !! Daniel's term: wavelength for diffusivity (not used)


       ! 3D terms
       real, public, dimension(GLOBAL_2D_ARRAY, nz) :: brk3dx, brk3de !! xi/eta-direction 3D breaking dissipation (rho)
       real, public, dimension(GLOBAL_2D_ARRAY, nz) :: ust       !! xi-direction 3D Stokes drift velocity (u-point)
       real, public, dimension(GLOBAL_2D_ARRAY, nz) :: vst       !! eta-direction 3D Stokes drift velocity (u-point)
       real, public, dimension(GLOBAL_2D_ARRAY, nz) :: wst,      !! vertical 3D Stokes drift velocity (rho-point)
     &                                                 Kappa,    !! 3D Bernoulli head term (rho-point)
     &                                                 kvf,      !! vertical vortex force term (K term, 3D, rho-point)
     &                                                 frc3dx, frc3de    !! xi/eta-direction 3D frictional dissipation (rho)


       ! terms with extra z-index
       real, public, dimension(GLOBAL_2D_ARRAY,0:nz) :: Akb,  !! breaking-wave-induced additional diffusivity (w-point
     &                                                  Akw,  !! wave-induced additional diffusivity (rho-point)
     &                                                  E_pre !! previous time-step value for Akw estimation (rho)

       ! DevinD added from set_forces.F
       real, public :: vstu,ustv,dudx,dvdx,dude,dvde


       logical :: analytical = .true.  ! analytical forcing

      !!!
      !!!!! DevinD - variables added as per DH's WAVE_PACKET forces.h for old MRL_WCI flag
      !!!

       ! wfrq | BBL/MRL | wind-induced wave frequency [rad/s]
       real, public, dimension(GLOBAL_2D_ARRAY) :: wfrq
       ! wdrx | MRL     | cosine of wave direction [non dimension]
       real, public, dimension(GLOBAL_2D_ARRAY) :: wdrx
       ! wdre | MRL     | sine of   wave direction [non dimension]
       real, public, dimension(GLOBAL_2D_ARRAY) :: wdre
       ! whrm | MRL     | (RMS) wave height (twice the wave amplitude) [m]
       real, public, dimension(GLOBAL_2D_ARRAY) :: whrm
       ! wdsp | MRL     | breaking dissipation rate (\epsilon_b term) [m3/s3]
       real, public, dimension(GLOBAL_2D_ARRAY) :: wdsp
       ! wdrg | MRL     | frictional dissipation rate (\epsilon_d term) [m3/s3]
       real, public, dimension(GLOBAL_2D_ARRAY) :: wdrg
       ! wbst | MRL/BKPP| frictional dissipation stress (e_d k/sigma) [m2/s2]
       real, public, dimension(GLOBAL_2D_ARRAY) :: wbst

       ! breaking dissipation [m3/s3]
       real, public, dimension(GLOBAL_2D_ARRAY) :: wveb
       ! frictional dissipation [m3/s3]
       real, public, dimension(GLOBAL_2D_ARRAY) :: wved


      !!!
      !!!!! DevinD - end variables added as per forces.h
      !!!

      ! DevinD - variables used in ana_wec_frc_tile and ana_init_wec_frc_tile

       ! For PART 1: Declare some variables from analytical.F
       real :: kkw, env, Clw, Cg

       ! DevinD - coefficient variable for external calculations
       real, public :: cff_ext, cff_ext2

       ! DevinD - timestep variables as per river_frc.F
       integer              :: it1 = 1  ! step n. Note: should be parameter as unchanging
       integer              :: it2 = 2  ! step n+1. Note: should be parameter as unchanging
       real, dimension(2)   :: rtim = [-99,-99] ! used to store both time steps (n, n+1)



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

# if !defined ANA_WWAVE && !defined WKB_WWAVE
        integer, public, parameter :: indxWAVE =71,       indxWWA=indxWAVE
     &                    , indxWWD=indxWAVE+1, indxWWP=indxWAVE+2
#  ifdef SWAN_OFFLINE
#    ifdef BBL
        integer, public, parameter :: indxWUB=indxWAVE+3
#    endif
#    ifdef MRL_WCI
        integer, public, parameter :: indxWEB=indxWAVE+4, indxWQB=indxWAVE+5
     &                    , indxWED=indxWAVE+6
#    endif
#  endif
# elif defined WKB_WWAVE
        integer, public, parameter :: indxHRM=71
     &                    , indxFRQ=indxHRM+1, indxWAC=indxHRM+2
     &                    , indxWKX=indxHRM+3, indxWKE=indxHRM+4
     &                    , indxEPB=indxHRM+5, indxEPD=indxHRM+6
     &                    , indxWAR=indxHRM+7, indxEPR=indxHRM+8
# endif



         integer, public, parameter :: indxSUP=81
     &                    , indxUST2D=indxSUP+1, indxVST2D=indxSUP+2
# ifdef SOLVE3D
     &                    , indxUST=indxSUP+3, indxVST=indxSUP+4
     &                    , indxWST=indxSUP+5, indxAkb=indxSUP+6
     &                    , indxAkw=indxSUP+7, indxKVF=indxSUP+8
     &                    , indxCALP=indxSUP+9, indxKAPS=indxSUP+10
# endif

      !!!!!!!!

        integer, public, target :: hisSUP, hisUST2D, hisVST2D ! Target for pointers if used in def_his_def_var_wec!!!
!      common /ncvars/ hisSUP, hisUST2D, hisVST2D

      type domainptr
        integer, pointer :: p ! For target hisUST2D above
      end type domainptr

# ifdef SOLVE3D
        integer, public, target :: hisUST, hisVST, hisAkb, hisAkw, hisKVF,
     &        hisCALP, hisKAPS, hisWST ! Target for pointers if used in def_his_def_var_wec!!!
!      common /ncvars/ hisUST, hisVST, hisAkb, hisAkw, hisKVF,
!     &        hisCALP, hisKAPS, hisWST
# endif

      integer, public :: avgSUP, avgUST2D, avgVST2D
!      common /ncvars/ avgSUP, avgUST2D, avgVST2D
#  ifdef SOLVE3D
      integer, public :: avgUST, avgVST, avgAkb, avgAkw, avgKVF,
     &        avgCALP, avgKAPS, avgWST
!      common /ncvars/ avgUST, avgVST, avgAkb, avgAkw, avgKVF,
!     &        avgCALP, avgKAPS, avgWST
#  endif

!      end type wec_ncvars

      ! DevinD - end wave terms previously contained in ncvars.h


! flags taken from def_his.F and used in def_his_def_var_wec & def_his_nf_inq_varid
#ifndef AVRH
# define ncid nchis
# define fname hisname
# define wrt wrthis
# ifdef WEC
#  define vidSUP hisSUP
#  define vidUST2D hisUST2D
#  define vidVST2D hisVST2D
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



      ! strings needed for wrt_his.F if using ncio method
!      character*15, public :: mynode_text
!      character*20, public :: mynode_nc_file
      ! end strings needed for wrt_his.F

      public init_scalars_wec
      public read_inp_wec
      public def_his_def_var_wec
      public def_his_nf_inq_varid
      public set_wec_frc
      public ana_wec_frc_tile
      public init_arrays_wec_tile
      public ana_init_wec_frc_tile
      public wstokes_wec_tile
      public compute_wec_rhs_uv_terms
      public wstokes_wec
      public wrt_his_wec



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
      vname(2,indxWST)  ='w vertical 3D Stokes drift velocity'
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
# endif

      end subroutine init_scalars_wec

      ! --------------------------------------------------------------------

        subroutine init_arrays_wec_tile (istr,iend,jstr,jend)

          ! This is called from init_arrays_tile of init_arrays.F
          ! It is needed to "first touch" arrays to optimize
          ! NUMA shared memory allocation on linux.

          implicit none

          integer, intent(in) :: istr, iend, jstr, jend
          ! Declare looping indices
          integer :: i, j

          real, parameter :: init=0.

          ! compute_auxiliary_bounds.h is needed for istrR, jstrR, etc
          ! This had to declared last else won't compile
#include "compute_auxiliary_bounds.h"

          ! This section of code is taken from DH's WAVE_PACKET init_arrays.F


              whrm(i,j) =init
              wdrx(i,j) =init
              wdre(i,j) =init
              wbst(i,j) =init

          do j=jstrR,jendR
            do i=istrR,iendR
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

# ifndef AVERAGES
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


# else
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
# endif /* AVERAGES */

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
      integer, external :: nf_def_var, nf_put_att_text ! Taken from netcdf.inc
      character(len=64) text
      integer, parameter :: nf_float = 5

      ! These variables are new and needed for looping of attributes
      type(domainptr), dimension(3) :: dom ! array of pointers to hisUST2d etc
      integer indx_loop ! looping index
      integer loop_grd(3,3) ! grid type - can make second dimension (i,:) for length 3 and 4
      integer :: i=0 ! additional looping integer

      ! 3D loop variables - may need to deallocate pointers at end of subroutine if memory leak!
      type(domainptr), dimension(8) :: dom3D ! array of pointers to hisUST2d etc
      integer loop_grd3D(8,4) ! grid type - can make second dimension (i,:) for length 3 and 4
      integer loop_dims(8) ! array of variable's dimensions

      ! Put in type of variable here.
      loop_grd(1,:) = r2dgrd; loop_grd(2,:) = u2dgrd; loop_grd(3,:) = v2dgrd
      ! Put in pointer to varid here (e.g. hisUST2D - see flags above)
      dom(1)%p => vidSUP; dom(2)%p => vidUST2D; dom(3)%p => vidVST2D

      do indx_loop = indxSUP, indxVST2D ! Loop through the relevant indices for wrt
        i=i+1 ! This is for loop index from 1
        if (wrt(indx_loop)) then
          lvar=lenstr(vname(1,indx_loop))
          ierr=nf_def_var (ncid, vname(1,indx_loop)(1:lvar), NF_FOUT,
     &                                      3, loop_grd(i,:), dom(i)%p)
# ifdef AVRH
          text='averaged '/ /vname(2,indx_loop)
# else
          text=vname(2,indx_loop)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, dom(i)%p, 'long_name', lvar,
     &                                               text(1:lvar))
          lvar=lenstr(vname(3,indx_loop))
          ierr=nf_put_att_text (ncid, dom(i)%p, 'units',     lvar,
     &                                 vname(3,indx_loop)(1:lvar))
#  ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, dom(i)%p, '_FillValue',
     &                                NF_FOUT, 1, spv_set)
#  endif

        endif
      end do

! THIS IS OLD CODE REPLACED BY LOOP ABOVE - LEFT IN UNTIL SURE NO MEMORY LEAK
! WHEN DOING LARGE SCALE NON-ANALYTICAL MODEL

! Wave-current interaction diagnostic variables.
! =============================================
!
! sup: quasi-static sea-level response (wave set-up/down)
!
!        if (wrt(indxSUP)) then
!          lvar=lenstr(vname(1,indxSUP))
!          ierr=nf_def_var (ncid, vname(1,indxSUP)(1:lvar), NF_FOUT,
!     &                                           3, r2dgrd, vidSUP)
!# ifdef AVRH
!          text='averaged '/ /vname(2,indxSUP)
!# else
!          text=vname(2,indxSUP)
!# endif
!          lvar=lenstr(text)
!          ierr=nf_put_att_text (ncid, vidSUP, 'long_name', lvar,
!     &                                             text(1:lvar))
!          lvar=lenstr(vname(3,indxSUP))
!          ierr=nf_put_att_text (ncid, vidSUP, 'units',     lvar,
!     &                                 vname(3,indxSUP)(1:lvar))
!#  ifdef MASK_LAND_DATA
!          ierr=nf_put_att_FOUT (ncid, vidiSUP, '_FillValue',
!     &                                NF_FOUT, 1, spv_set)
!#  endif
!        endif

!
! ust2d & vst2d: 2D, depth-averaged Stokes drift velocities.
!
!        if (wrt(indxUST2D)) then
!          lvar=lenstr(vname(1,indxUST2D))
!          ierr=nf_def_var (ncid, vname(1,indxUST2D)(1:lvar), NF_FOUT,
!     &                                           3, u2dgrd, vidUST2D)
!# ifdef AVRH
!          text='averaged '/ /vname(2,indxUST2D)
!# else
!          text=vname(2,indxUST2D)
!# endif
!          lvar=lenstr(text)
!          ierr=nf_put_att_text (ncid, vidUST2D, 'long_name', lvar,
!     &                                               text(1:lvar))
!          lvar=lenstr(vname(3,indxUST2D))
!          ierr=nf_put_att_text (ncid, vidUST2D, 'units',     lvar,
!     &                                 vname(3,indxUST2D)(1:lvar))
!#  ifdef MASK_LAND_DATA
!          ierr=nf_put_att_FOUT (ncid, vidUST2D, '_FillValue',
!     &                                NF_FOUT, 1, spv_set)
!#  endif
!
!        endif

!        if (wrt(indxVST2D)) then
!          lvar=lenstr(vname(1,indxVST2D))
!          ierr=nf_def_var (ncid, vname(1,indxVST2D)(1:lvar), NF_FOUT,
!     &                                           3, v2dgrd, vidVST2D)
!# ifdef AVRH
!          text='averaged '/ /vname(2,indxVST2D)
!# else
!          text=vname(2,indxVST2D)
!# endif
!          lvar=lenstr(text)
!          ierr=nf_put_att_text (ncid, vidVST2D, 'long_name', lvar,
!     &                                               text(1:lvar))
!          lvar=lenstr(vname(3,indxVST2D))
!          ierr=nf_put_att_text (ncid, vidVST2D, 'units',     lvar,
!     &                                  vname(3,indxVST2D)(1:lvar))
!#  ifdef MASK_LAND_DATA
!          ierr=nf_put_att_FOUT (ncid, vidVST2D, '_FillValue',
!     &                                NF_FOUT, 1, spv_set)
!#  endif
!        endif
# ifdef SOLVE3D

      ! Put in type of variable here.
      loop_grd3D(1,:) = u3dgrd; loop_grd3D(2,:) = v3dgrd; loop_grd3D(3,:) = r3dgrd
      loop_grd3D(4,:) = w3dgrd; loop_grd3D(5,:) = w3dgrd; loop_grd3D(6,:) = r3dgrd
      loop_grd3D(7,1:3) = r2dgrd; loop_grd3D(8,1:3) = r2dgrd ! Potential bug - 3 not 4 length integers
      ! Put in pointer to varid here (e.g. hisUST2D - see flags above)
      dom3D(1)%p => vidUST; dom3D(2)%p => vidVST; dom3D(3)%p => vidWST
      dom3D(4)%p => vidAkb; dom3D(5)%p => vidAkw; dom3D(6)%p => vidKVF
      dom3D(7)%p => vidCALP; dom3D(8)%p => vidKAPS
      ! The number of dimensions of variable (includes time)
      loop_dims = [4,4,4,4,4,4,3,3]

      ! Could avoid having whole extra loop by using if statement for loop over 3d
      ! variables or not. I.e. add to sup, ust2d and vst2d loop above.
      i=0 ! reset loop counter
      do indx_loop = indxUST, indxKAPS ! Loop through relevant indices for wrt
        i=i+1 ! This is for loop index from 1
                if (wrt(indx_loop)) then
          lvar=lenstr(vname(1,indx_loop))
          ierr=nf_def_var (ncid, vname(1,indx_loop)(1:lvar), NF_FOUT,
     &                        loop_dims(i), loop_grd3D(i,:), dom3D(i)%p)
# ifdef AVRH
          text='averaged '/ /vname(2,indx_loop)
# else
          text=vname(2,indx_loop)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, dom3D(i)%p, 'long_name', lvar,
     &                                             text(1:lvar))
          lvar=lenstr(vname(3,indx_loop))
          ierr=nf_put_att_text (ncid, dom3D(i)%p, 'units',     lvar,
     &                                 vname(3,indx_loop)(1:lvar))
#  ifdef MASK_LAND_DATA
          ierr=nf_put_att_FOUT (ncid, dom3D(i)%p, '_FillValue',
     &                                NF_FOUT, 1, spv_set)
#  endif
        endif

      end do

! THIS IS OLD CODE REPLACED BY LOOP ABOVE - LEFT IN UNTIL SURE NO MEMORY LEAK
! WHEN DOING LARGE SCALE NON-ANALYTICAL MODEL

!
! ust3d & vst3d: 3D Stokes drift velocities.
!
!        if (wrt(indxUST)) then
!          lvar=lenstr(vname(1,indxUST))
!          ierr=nf_def_var (ncid, vname(1,indxUST)(1:lvar), NF_FOUT,
!     &                                           4, u3dgrd, vidUST)
!# ifdef AVRH
!          text='averaged '/ /vname(2,indxUST)
!# else
!          text=vname(2,indxUST)
!# endif
!          lvar=lenstr(text)
!          ierr=nf_put_att_text (ncid, vidUST, 'long_name', lvar,
!     &                                             text(1:lvar))
!          lvar=lenstr(vname(3,indxUST))
!          ierr=nf_put_att_text (ncid, vidUST, 'units',     lvar,
!     &                                 vname(3,indxUST)(1:lvar))
!#  ifdef MASK_LAND_DATA
!          ierr=nf_put_att_FOUT (ncid, vidUST, '_FillValue',
!     &                                NF_FOUT, 1, spv_set)
!#  endif
!        endif

!        if (wrt(indxVST)) then
!          lvar=lenstr(vname(1,indxVST))
!          ierr=nf_def_var (ncid, vname(1,indxVST)(1:lvar), NF_FOUT,
!     &                                           4, v3dgrd, vidVST)
!# ifdef AVRH
!          text='averaged '/ /vname(2,indxVST)
!# else
!          text=vname(2,indxVST)
!# endif
!          lvar=lenstr(text)
!          ierr=nf_put_att_text (ncid, vidVST, 'long_name', lvar,
!     &                                             text(1:lvar))
!          lvar=lenstr(vname(3,indxVST))
!          ierr=nf_put_att_text (ncid, vidVST, 'units',     lvar,
!     &                                 vname(3,indxVST)(1:lvar))
!#  ifdef MASK_LAND_DATA
!          ierr=nf_put_att_FOUT (ncid, vidVST, '_FillValue',
!     &                                NF_FOUT, 1, spv_set)
!#  endif
!        endif
!
! wst: vertical Stokes drift velocity (m/s) at rho-point
!
!        if (wrt(indxWST)) then
!          lvar=lenstr(vname(1,indxWST))
!          ierr=nf_def_var (ncid, vname(1,indxWST)(1:lvar), NF_FOUT,
!     &                                           4, r3dgrd, vidWST)
!# ifdef AVRH
!          text='averaged '/ /vname(2,indxWST)
!# else
!          text=vname(2,indxWST)
!# endif
!          lvar=lenstr(text)
!          ierr=nf_put_att_text (ncid, vidWST, 'long_name', lvar,
!     &                                              text(1:lvar))
!          lvar=lenstr(vname(3,indxWST))
!          ierr=nf_put_att_text (ncid, vidWST, 'units',     lvar,
!     &                                 vname(3,indxWST)(1:lvar))
!#  ifdef MASK_LAND_DATA
!          ierr=nf_put_att_FOUT (ncid, vidWST, '_FillValue',
!     &                                NF_FOUT, 1, spv_set)
!#  endif
!        endif
!
! Akb: Vertical eddy viscosity coefficient due to depth-induced wave breaking.
!
!        if (wrt(indxAkb)) then
!          lvar=lenstr(vname(1,indxAkb))
!          ierr=nf_def_var (ncid, vname(1,indxAkb)(1:lvar), NF_FOUT,
!     &                                           4, w3dgrd, vidAkb)
!# ifdef AVRH
!          text='averaged '/ /vname(2,indxAkb)
!# else
!          text=vname(2,indxAkb)
!# endif
!          lvar=lenstr(text)
!          ierr=nf_put_att_text (ncid, vidAkb, 'long_name', lvar,
!     &                                              text(1:lvar))
!          lvar=lenstr(vname(3,indxAkb))
!          ierr=nf_put_att_text (ncid, vidAkb, 'units',     lvar,
!     &                                  vname(3,indxAkb)(1:lvar))
!#  ifdef MASK_LAND_DATA
!          ierr=nf_put_att_FOUT (ncid, vidAkb, '_FillValue',
!     &                                NF_FOUT, 1, spv_set)
!#  endif
!        endif
! Akw: Vertical eddy diffusivity coefficient due to primary waves.
!
!        if (wrt(indxAkw)) then
!          lvar=lenstr(vname(1,indxAkw))
!          ierr=nf_def_var (ncid, vname(1,indxAkw)(1:lvar), NF_FOUT,
!     &                                           4, w3dgrd, vidAkw)
!# ifdef AVRH
!          text='averaged '/ /vname(2,indxAkw)
!# else
!          text=vname(2,indxAkw)
!# endif
!          lvar=lenstr(text)
!          ierr=nf_put_att_text (ncid, vidAkw, 'long_name', lvar,
!     &                                              text(1:lvar))
!          lvar=lenstr(vname(3,indxAkw))
!          ierr=nf_put_att_text (ncid, vidAkw, 'units',     lvar,
!     &                                  vname(3,indxAkw)(1:lvar))
!#  ifdef MASK_LAND_DATA
!          ierr=nf_put_att_FOUT (ncid, vidAkw, '_FillValue',
!     &                                NF_FOUT, 1, spv_set)
!#  endif
!        endif
!
! kvf: vertical vortex force term (u^St du/dz)
!
!        if (wrt(indxKVF)) then
!          lvar=lenstr(vname(1,indxKVF))
!          ierr=nf_def_var (ncid, vname(1,indxKVF)(1:lvar), NF_FOUT,
!     &                                           4, r3dgrd, vidKVF)
!# ifdef AVRH
!          text='averaged '/ /vname(2,indxKVF)
!# else
!          text=vname(2,indxKVF)
!# endif
!          lvar=lenstr(text)
!          ierr=nf_put_att_text (ncid, vidKVF, 'long_name', lvar,
!     &                                             text(1:lvar))
!          lvar=lenstr(vname(3,indxKVF))
!          ierr=nf_put_att_text (ncid, vidKVF, 'units',     lvar,
!     &                                 vname(3,indxKVF)(1:lvar))
!#  ifdef MASK_LAND_DATA
!          ierr=nf_put_att_FOUT (ncid, vidKVF, '_FillValue',
!     &                                NF_FOUT, 1, spv_set)
!#  endif
!        endif

!
! calP: surface pressure correction term appeared in prsgrd term
!
!        if (wrt(indxCALP)) then
!          lvar=lenstr(vname(1,indxCALP))
!          ierr=nf_def_var (ncid, vname(1,indxCALP)(1:lvar), NF_FOUT,
!     &                                           3, r2dgrd, vidCALP)
!# ifdef AVRH
!          text='averaged '/ /vname(2,indxCALP)
!# else
!          text=vname(2,indxCALP)
!# endif
!          lvar=lenstr(text)
!          ierr=nf_put_att_text (ncid, vidCALP, 'long_name', lvar,
!     &                                              text(1:lvar))
!          lvar=lenstr(vname(3,indxCALP))
!          ierr=nf_put_att_text (ncid, vidCALP, 'units',     lvar,
!     &                                 vname(3,indxCALP)(1:lvar))
!#  ifdef MASK_LAND_DATA
!          ierr=nf_put_att_FOUT (ncid, vidCALP, '_FillValue',
!     &                                NF_FOUT, 1, spv_set)
!#  endif
!        endif
!
! Kapsrf: surface Bernoulli head term appeared in prsgrd term
!
!        if (wrt(indxKAPS)) then
!          lvar=lenstr(vname(1,indxKAPS))
!          ierr=nf_def_var (ncid, vname(1,indxKAPS)(1:lvar), NF_FOUT,
!     &                                           3, r2dgrd, vidKAPS)
!# ifdef AVRH
!          text='averaged '/ /vname(2,indxKAPS)
!# else
!          text=vname(2,indxKAPS)
!# endif
!          lvar=lenstr(text)
!          ierr=nf_put_att_text (ncid, vidKAPS, 'long_name', lvar,
!     &                                              text(1:lvar))
!          lvar=lenstr(vname(3,indxKAPS))
!          ierr=nf_put_att_text (ncid, vidKAPS, 'units',     lvar,
!     &                                 vname(3,indxKAPS)(1:lvar))
!#  ifdef MASK_LAND_DATA
!          ierr=nf_put_att_FOUT (ncid, vidKAPS, '_FillValue',
!     &                                NF_FOUT, 1, spv_set)
!#  endif
!        endif
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


      ! --------------------------------------------------------------------
      ! --------------------------------------------------------------------
      ! --------------------------------------------------------------------

        subroutine ana_init_wec_frc_tile (istr,iend,jstr,jend)

         ! Initialize wec forcing - chose to include because DH's WAVE_PACKET
         ! version is separate, and init version doesn't include time, whereas
         ! ana_wec_frc_tile does include time.
         ! Might be able to just add if statement to ana_wec_frc_tile instead?

         implicit none

         integer, intent(in) :: istr, iend, jstr, jend
         ! Declare looping indices
         integer :: i, j, k

         real :: cff1 ! This could be declared in module but not sure - here for now.
         real :: cff2 ! DevinD added for debug check

          ! compute_auxiliary_bounds.h is needed for istrR, jstrR, etc
          ! This had to appear after other variables declared else wont compile
#include "compute_auxiliary_bounds.h"


         ! Taken from ana_init.F of DH's WAVE_PACKET code

         real :: stk(PRIVATE_2D_SCRATCH_ARRAY)

          kkw=2.0*pi/1.0 ! Changed kw to kkw so as not to conflict part 2 also has kw
          env=0.0001
          Clw=sqrt(g*h(2,2))
          Cg=0.5*sqrt(g/kkw)
!           do j=jstrR,jendR
!            do i=istrR,iendR
          do j=jstr-1,jend+1
            do i=istr-1,iend+1
             do k=1,N
              t(i,j,k,1,itemp)=18.0
              t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
             enddo
             ! DevinD - whrm from Eq(29) of DH's 'A propagating wave packet.pdf
!             cff1=(dble(i)+dble(iSW_corn)-(LLm/2))/pm(2,2)
!             whrm(i,j)  = 0.001*2.0*exp(-env*cff1**2)       ! Hrms=2*A
             cff1=modulo(time,(gnx/Cg))        ! Changed LLm to gnx
             cff2=(dble(i)-(Cg*cff1)
     &                +dble(iSW_corn)-(gnx/2))/pm(2,2)    ! initial shift mid domain
             whrm(i,j)= 0.001*2.0*(exp(-env*cff2**2)
     &                      +exp(-env*(cff2+gnx/pm(2,2))**2))
             wfrq(i,j)   = sqrt(g*kkw)
             wdrx(i,j)   = 1.0                             ! cosine wave direction (xi)
             wdre(i,j)   = 0.0                             ! sine wave direction (eta)
             wdsp(i,j)   = 0.0
             sup(i,j)    = 0.0
    !         sup(i,j)    = -(0.25*whrm(i,j)**2)*kw/(2.*sinh(2*kw*h(i,j)))   ! A**2 = 0.25 Hrms**2
!             ust2d(i,j)  = (0.125*whrm(i,j)**2)*sqrt(g*kkw)/h(i,j)
             ! DevinD - stk & whrm appear to be at u-points because loop
             ! is through integers and i represents x. For this grid u-points
             ! are on the integer meters, and rho-points between on .5m
             stk(i,j)  = (0.125*whrm(i,j)**2)*sqrt(g*kkw)/h(i,j)
             if (i>istr-1) then
                ! DevinD - ust2d averaged to rho-point, but it seems this
                ! is the wrong conversion process for u-to-rho conversion.
                ust2d(i,j)= 0.5*(stk(i-1,j)*wdrx(i-1,j)+stk(i,j)*wdrx(i,j))
!                ubar(i,j,1) = -((Clw**2)/((Clw**2)-(Cg**2)))*
!     &                         ust2d(i,j)
                ! ubar should be at u-point - hence using stk
                ubar(i,j,1) = -((Clw**2)/((Clw**2)-(Cg**2)))*
     &                         stk(i,j)
!                ! DevinD - attempted hack to correct phase shift stk(i-1)
!                ubar(i,j,1) = -((Clw**2)/((Clw**2)-(Cg**2)))*
!     &                         stk(i-1,j)
!                zeta(i,j,1) = ubar(i,j,1)*Cg/g
                ! zeta should be at rho-point - hence using ust2d
                zeta(i,j,1) = -((Clw**2)/((Clw**2)-(Cg**2)))*
     &                         ust2d(i,j)*Cg/g
                brk2dx(i,j) = 0.0
                brk2de(i,j) = 0.0
                frc2dx(i,j) = 0.0
                frc2de(i,j) = 0.0

                vbar(i,j,1)= 0.0
                vst2d(i,j) = 0.0
                do k=1,N
                    u(i,j,k,1) = ubar(i,j,1)                     ! actually need
                    v(i,j,k,1) = vbar(i,j,1)                     ! to init u even if baro
                    ust(i,j,k) = ust2d(i,j)
                    vst(i,j,k) = vst2d(i,j)
                enddo
             end if
           enddo
          enddo

      ! DevinD  debug output
!      write(*,*) 'ust0(1,1,1)', ust(1,1,1)
!      write(*,*) 'ust0(1,1,2)', ust(1,1,2)

        end subroutine ana_init_wec_frc_tile

       ! --------------------------------------------------------------------




        subroutine ana_wec_frc_tile (istr,iend,jstr,jend)
         ! Analytical river forcing data

         implicit none

         integer, intent(in) :: istr, iend, jstr, jend
         ! Declare looping indices
         integer :: i, j

         ! get all cppdefs flags
! # include "cppdefs.h" - included in top of module (left here to remember where they were needed)
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

!          ! For PART 1: Declare some variables from analytical.F
!          real :: kkw, env, Cg - in top of module

          ! For PART 2: declare variables needed
          ! note - commented out SEEMINGLY unneeded variables
          real :: cff,cff1,cff2,cff3,cff4,
     &      inv_g,
     &      khd, kh,
     &      umag,
     &      eps,
     &      wave_ramp,
     &      wramp2,kbrk,fb,fb0,fb1,fb2,inv_fbs,
     &      bconst,tauc,tauw,Uorb1,ka_f00,hz0,z_tide,
     &      dtinv,inv_zb,a_brk,a_kv,c1o3,c4o3,inv_k,dd,ust_ker,intfb,
     &      fbsrf,fn1,fn2,khmax,bz1,bz2,beta,abot,delta,
     &      a_frc,kfrc

          real ::  wh(PRIVATE_2D_SCRATCH_ARRAY),
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



!          real, parameter :: eps=1.e-10,
          parameter ( eps=1.e-10,
# ifdef STOKESEKMAN
     &                bconst=0.D0,
# else
     &                bconst=0.1D0, !0.03D0,               ! breaking contribution to KPP
# endif
     &                a_kv=1.2D0, !1.2D0                 ! breaking scale for eddy visc.
     &                a_brk=0.2D0,                ! breaking scale for body force
     &                a_frc=3.D0,                 ! friction scale for body force
     &                khmax=20.D0,                ! deep-water limit for k x dep
     &                c1o3=0.3333333333333333,    ! 1/3
     &                c4o3=1.3333333333333333     ! 4/3
     &               )


# ifdef SOLVE3D
          integer  k,kk
          real     wrk1(PRIVATE_2D_SCRATCH_ARRAY,0:nz),
     &             wrk2(PRIVATE_2D_SCRATCH_ARRAY,0:nz),
     &             wrk3(PRIVATE_2D_SCRATCH_ARRAY,0:N),
     &             wrk4(0:N),wrk5(0:N),wrk6(0:N), kvsurf
# endif


          ! compute_auxiliary_bounds.h is needed for istrR, jstrR, etc
          ! This had to appear after eps=1.e-10 else wont compile
#include "compute_auxiliary_bounds.h"

      ! DevinD - this flag is needed for kh variable used for ust/vst
# define KH_HUNT

          !!!!!!!!!!
          ! If first step - Cancelled as put in ana_init.F instead
          !!!!!!!!!!        This is taken from river_frc.F read

!          if (rtim(it1) < 0) then !! first time
!
!
!            call ana_init_wec_frc_tile(istr,iend,jstr,jend)
!
!
!          !!!!!!!!!!
!          ! Else timestep > 0
!          !!!!!!!!!!
!
!          else


            !*** PART 1 - from DH's WAVE_PACKET analytical.F ***!

          ! The following code is taken from analytical.F L995 from DH's
          ! WAVE_PACKET example code



          kkw=2.0*pi/1.0 ! Changed kw to kkw so as not to conflict part 2
          env=0.0001
          Cg=0.5*sqrt(g/kkw)
          do j=jstr-1,jend+1
            do i=istr-1,iend+1
              cff1=modulo(time,(gnx/Cg))        ! Changed LLm to gnx
              cff2=(dble(i)-(Cg*cff1)
     &                +dble(iSW_corn)-(gnx/2))/pm(2,2)    ! initial shift mid domain
              whrm(i,j)= 0.001*2.0*(exp(-env*cff2**2)
     &                      +exp(-env*(cff2+gnx/pm(2,2))**2))
              wfrq(i,j)= sqrt(g*2.0*pi)
              wdrx(i,j)= 1.0                             ! cosine wave direction (xi)
              wdre(i,j)= 0.0                             ! sine wave direction (eta)
              wveb(i,j)= 0.0
              wved(i,j)= 0.0
            enddo
          enddo

            !*** END PART 1 ***!

          !!!!!
          ! Taken from mrl_wci.F
          !!!!!

            !*** PART 2 - from DH's WAVE_PACKET mrl_wci.F ***!

          wave_ramp = 1.0
!          wramp2= wave_ramp**2
          inv_g = 1.0/g
!
!        ! Evaluate 2DH wave-current interaction variable.
!        ! ================================================
!        !
!        ! 2DH Stokes velocities, breaking, roller and bottom-friction
!        ! dissipation terms, defined at horizontal rho-points.

          ! Note: this would normally involve A3d and A2d scratch
          ! variables, but I have ignored them for ease now.

          do j=jstr-1,jend+1
            do i=istr-1,iend+1
!              Dstp(i,j)=max(h(i,j)+zeta(i,j,knew),eps)
              Dstp(i,j)=h(i,j)
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
#  if defined KH_SOULSBY
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
          ebrk(i,j) =0            ! epsilon_b (m3/s3) kaida DD

#  ifdef SURFACE_ROLLER
          erol(i,j) =wramp2*0.5762*Dstp(i,j)*    ! epsilon_r (m3/s3)
     &               whrm(i,j)*wvqb(i,j)*        ! 0.5762=0.06*sinb*g^2
     &               kw(i,j)*inv_f(i,j)          ! Nairn et al (1991)
          erol(i,j) =0            ! DD
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
!          stk(i,j) =act(i,j)*kw(i,j)*inv_d(i,j)  ! depth-averaged Stokes
          stk(i,j) = 0.125*(whrm(i,j)**2)*wfrq(i,j)   ! DD
     &                      /Dstp(i,j)
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
    !      print *, 'depth-induced wave breaking deactivated'
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
      do j=jstrR,jendR
        do i=istrR,iendR
          sup(i,j) =0.0
!          sup(i,j) =-0.125*(wh(i,j)**2)*kw(i,j)
!     &              /sinh(2.*kD(i,j))
!     &              /max(sinh(2.*min(kD(i,j),khmax)),eps)
#  ifdef MASKING
     &                                        *rmask(i,j)
#  endif
        enddo
      enddo
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
!# elif !defined WKB_WWAVE & defined SWAN_OFFLINE
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
          ust2d(i,j)  =0.5*(stk(i-1,j)*wdrx(i-1,j)+stk(i,j)*wdrx(i,j))
          brk2dx(i,j) =0.5*(brk(i-1,j)*wdrx(i-1,j)+brk(i,j)*wdrx(i,j))
          frc2dx(i,j) =0.5*(frc(i-1,j)*wdrx(i-1,j)+frc(i,j)*wdrx(i,j))
#  ifdef MASKING
#   ifdef SHOREFACE
          ust2d(i,j)  = ust2d(i,j)*umask(i,j)
#   endif
          brk2dx(i,j) =brk2dx(i,j)*umask(i,j)
          frc2dx(i,j) =frc2dx(i,j)*umask(i,j)
#  endif
        enddo
      enddo

            ! DevinD - I don't know why vst2d is not included in loop
            ! above for efficiency, but I have kept it separate.

      do j=jstr,jendR
        do i=istrR,iendR
          vst2d(i,j)  =0.5*(stk(i,j-1)*wdre(i,j-1)+stk(i,j)*wdre(i,j))
          brk2de(i,j) =0.5*(brk(i,j-1)*wdre(i,j-1)+brk(i,j)*wdre(i,j))
          frc2de(i,j) =0.5*(frc(i,j-1)*wdre(i,j-1)+frc(i,j)*wdre(i,j))
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
!# if !defined BBL && (defined BBL_S95 || defined BBL_F00 || defined BBL_QUAD2D)
!      inv_zb=1.0/max(Zob,eps)
! continues here but have left out as does not seem relevant and a lot of code

#ifdef SOLVE3D
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
# define ust_r wrk1
# define inv_dz wrk2
# define inv_ex inv_f
!
! 3D depth-dependent Stokes drift velocities
! ==========================================
!

      do j=jstr-1,jend+1
        do i=istr-1,iend+1

          inv_ex(i,j)=1.0/(1.0-exp(-4.*kD(i,j))
     &                        -2.0*exp(-2.*kD(i,j)))
!          inv_ex(i,j)=1.0/max(1.0-exp(-4.*kD(i,j))
!     &                        -2.0*exp(-2.*kD(i,j)),eps)
          ust_ker=(act(i,j)*(sqrt(g*kw(i,j)*tanh(kD(i,j)))**2)/9.81)
     &                    *inv_ex(i,j)
!          ust_ker=act(i,j)*kw(i,j)*inv_ex(i,j)
          cff1 = exp( 2.*kw(i,j)*(z_w(i,j,0)+h(i,j)-Dstp(i,j)))
     &          -exp(-2.*kw(i,j)*(z_w(i,j,0)+h(i,j)+Dstp(i,j)))
          do k=1,N,+1              ! <-- irreversible
            inv_dz(i,j,k)=1./Hz(i,j,k)
            cff2 = exp( 2.*kw(i,j)*(z_w(i,j,k)+h(i,j)-Dstp(i,j)))
     &            -exp(-2.*kw(i,j)*(z_w(i,j,k)+h(i,j)+Dstp(i,j)))
            ust_r(i,j,k) =ust_ker*inv_dz(i,j,k)*( cff2-cff1 )
            cff1=cff2
!
! Horizontal interpolation of ust_r to u- & v- points.
! *** Notice that mask should NOT be applied here.
!
            if (i>istr-1) then
                    ust(i,j,k)=0.5*( ust_r(i-1,j,k)*wdrx(i-1,j)
     &                        +ust_r(i,j,k)*wdrx(i,j) )
            end if
            if (j>jstr-1) then
                    vst(i,j,k)=0.5*( ust_r(i,j-1,k)*wdre(i,j-1)
     &                        +ust_r(i,j,k)*wdre(i,j) )
            end if
          enddo
        enddo
      enddo                 ! <-- keep ust_r, inv_dz, inv_ex

      ! DevinD  debug output
!      write(*,*) 'act(1,1)', act(1,1)
!      write(*,*) 'inv_ex(1,1)', inv_ex(1,1)
!      write(*,*) 'kh', kh
!      write(*,*) 'kD(1,1)', kD(1,1)
!      write(*,*) 'tanh(kD(1,1))', tanh(kD(1,1))
!      write(*,*) 'ust_ker', ust_ker
!      write(*,*) 'inv_dz(1,1,2)', inv_dz(1,1,2)
!      write(*,*) 'ust_r(1,1,2)', ust_r(1,1,2)
!      write(*,*) 'ust(1,1,1)', ust(1,1,1)
!      write(*,*) 'ust(1,1,2)', ust(1,1,2)

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
          do k=1,N  ! DevinD change N to nz    ! at rho-point
            kv(k) =0.5*kw(i,j)*(
     &         wdrx(i,j)* ( u(i,j,k,nstp)+u(i+1,j,k,nstp) )
     &       + wdre(i,j)* ( v(i,j,k,nstp)+v(i,j+1,k,nstp) ) )
          enddo
          kvsurf =1.5*kv(N)-0.5*kv(N-1) ! extrapolate to surface
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
          cff4 =-2.0*kw(i,j)*kvsurf
!ccccc          actp =0.125*g*(wh(i,j)**2)/max(fr(i,j),eps)
!ccccc          calP(i,j) = actp*inv_g*tanh(kD(i,j))
          actp=act(i,j) ! DevinD - this is unnecessary as is define actp
          calP(i,j) = actp*tanh(kD(i,j))
     &                    *( cff1+cff2+cff3+cff4 )
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
#  ifdef KAPSRF_SAFE
c#   define FLT 0.98
#   define FLT 0.5
          Kapsrf(i,j) = ( FLT*Kapsrf(i,j)  ! time-filter to suppress noises
     &                   +(1.0-FLT)*cff*actp*inv_ex(i,j) )
#  else
          Kapsrf(i,j) = cff*actp*inv_ex(i,j)
#  endif
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
            kvr(k) = 0.25*(ust_r(i,j,k)+ust_r(i,j,k+1))*
     &                ( wdrx(i,j)*(u(i,j,k+1,nstp)-u(i,j,k,nstp)
     &                        +u(i+1,j,k+1,nstp)-u(i+1,j,k,nstp))
     &                 +wdre(i,j)*(v(i,j,k+1,nstp)-v(i,j,k,nstp)
     &                        +v(i,j+1,k+1,nstp)-v(i,j+1,k,nstp)))
     &                                 /(z_r(i,j,k+1)-z_r(i,j,k))
          enddo
          kvr(0)=2.*kvr(1)-kvr(2)    ! 0.D0    ! bottom & top B.Cs.
          kvr(N)=2.*kvr(N-1)-kvr(N-2)
          do k=1,N                             ! move to vert rho-point
            kvf(i,j,k)=0.5*( kvr(k)+kvr(k-1) ) ! copy into a shared array
#  ifdef MASKING
     &                             *rmask(i,j)
#  endif
          enddo
        enddo
      enddo           ! <-- discard ust_r, kvr
# undef ust_r
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
!      cff=wave_ramp/dt
      do j=jstrR,jendR
        do i=istrR,iendR
!          cff2 = 0.25*wh(i,j)/max(1.-exp(-2.*kD(i,j)),eps)
          do k=0,N
!            cff1 = cff2*( exp( kw(i,j)*(z_w(i,j,k)+h(i,j)-Dstp(i,j)))
!     &                   -exp(-kw(i,j)*(z_w(i,j,k)+h(i,j)+Dstp(i,j))) )
!            cff3 = cff1**2
            Akw(i,j,k) = 0.D0
!            Akw(i,j,k) = cff*(cff3-E_pre(i,j,k))
#  ifdef MASKING
     &                               *rmask(i,j)
#  endif
!            E_pre(i,j,k) = cff3   ! store this for the next time step.
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
          enddo
          Akb(i,j,0)=0.D0
          Akb(i,j,N)=max(1.5*kb0(N)-0.5*kb0(N-1),0.0)
!            if (mynode.eq.1) then
!              print*,' '
!              print*, Akb(10,10,k)
!              print*,' '
!            endif
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
#  ifdef MASKING
     &                                         *umask(i,j)
#  endif
          enddo
        enddo
        do j=jstr,jendR
          do i=istrR,iendR
            brk3de(i,j,k)=0.5*( brk_r(i,j-1,k)*wdre(i,j-1)
     &                           +brk_r(i,j,k)*wdre(i,j) )
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


#endif /* SOLVE3D */

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


        end subroutine ana_wec_frc_tile

       ! --------------------------------------------------------------------
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
          do k=N-1,2,-1
            wst(i,j,k)=+0.5625*(Wrk(k  )+Wrk(k-1))
     &                 -0.0625*(Wrk(k+1)+Wrk(k-2))
          enddo
          wst(i,j,1)= -0.125*Wrk(2) +0.75*Wrk(1) +0.375*Wrk(0)
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

!
! Add combined vortex-force and advection terms & breaking terms
! 2nd- and 4th-order centered schemes, or 2nd with local weighted
! filter (with VF_ROBUST) are available.
!
# if defined  WEC && defined UV_ADV
#  define utmp wrkone
#  define vtmp wrktwo
#  ifdef VF_ROBUST
      real :: gam1=0.1         ! local weighted filter (gam1 may be 0-1/3)
      real :: gam2=1.-3.*gam1  ! gam1=0 to revert back to 2nd-order centered
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

        subroutine read_wec_frc  !(time)
         ! Reads WEC forcing data from input file

         ! Not applicable for analytical version

        end subroutine read_wec_frc


               ! --------------------------------------------------------------------
       subroutine set_wec_frc(time)
         ! Use an analytical forcing or interpolate two instances
         ! of the forcing that are read from input files
         implicit none
         ! import/export
         real(kind=8), intent(in) :: time
         ! local
!         real    :: cf1,cf2
         real    :: tday

         tday = time/(24*3600) ! convert from time seconds to days

         if (analytical) then

!            call ana_wec_frc(tday)
        ! This is from river_frc but I don't believe it is relevant
        ! because ana_wec_frc_tile is called from set_forces.F

         else
           ! check whether we need to update forcing data
!           if (rtim(it2) < tday) then
!             call read_wec_frc(tday)
!           endif

!           ! linear interpolation in time
!           cf1 = ( rtim(it2) - tday )/(rtim(it2) - rtim(it1) )
!           cf2 = ( tday - rtim(it1) )/(rtim(it2) - rtim(it1) )
!           riv_vol = cf1*rvol(:,it1) + cf2 * rvol(:,it2)
!           riv_trc = cf1*rtrc(:,:,it1) + cf2 * rtrc(:,:,it2)
!          if (mynode==0) then
!            print *, 'volume: ',riv_vol(1:nriv_frc)
!          endif
!           ! tidal modulation of river volume
!           riv_vol = (cos(2*pi*tday) + 1)*riv_vol

           stop 'no analytical river data'

         endif
       end subroutine set_wec_frc


       ! --------------------------------------------------------------------

       subroutine find_wec_new  !(time,v_name,vt_name,ifile,irec)

        ! Get nc file model input data

       end subroutine find_wec_new

!        ! -------------------------------------
!        subroutine ana_wec_tile(istr,jstr,iend,jend)
!         !analytical wec forcing
!         ! Set everything to zero for testing
!         implicit none
!         real(kind=8), intent(in) :: istr,jstr,iend,jend
!         integer :: i, j, k
!
!          ! rho-point variables
!          do i=istr,iend
!            do j=jstr,jend
!               sup(i,j)    = 0.D0
!               brk2dx(i,j) = 0.D0
!               brk2de(i,j) = 0.D0
!               frc2dx(i,j) = 0.D0
!               frc2de(i,j) = 0.D0
!               calP(i,j)   = 0.D0
!               Kapsrf(i,j) = 0.D0
!!               lm(i,j)     = 0.D0
!               do k=1,nz
!                  brk3dx(i,j,k) = 0.D0
!                  brk3de(i,j,k) = 0.D0
!                  wst(i,j,k)    = 0.D0
!                  Kappa(i,j,k)  = 0.D0 !may not need this variable
!                  kvf(i,j,k)    = 0.D0
!                  Akw(i,j,k)    = 0.D0
!                  frc3dx(i,j,k) = 0.D0
!                  frc3de(i,j,k) = 0.D0
!                enddo
!                do k=0,nz
!                   Akb(i,j,k) = 0.D0
!                enddo
!             enddo
!          enddo
!
!          !u-point variables
!          do i=istr,iend-1
!           do j=jstr,jend
!              ust2d(i,j) = 0.D0
!              do k=1,nz
!                 ust(i,j,k) = 0.D0
!              enddo
!           enddo
!          enddo
!
!          !v-point variables
!          do i=istr,iend
!           do j=jstr,jend-1
!              vst2d(i,j) = 0.D0
!              do k=1,nz
!                 vst(i,j,k) =0.D0
!              enddo
!            enddo
!          enddo
!
!        end subroutine ana_wec_tile


        ! ----------------------------------------


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


        ! ----------------------------------------
                 
!        subroutine wstokes_tile(istr,jstr,iend,jend, pm, pn)
         ! Compute vertical stokes velocity         
!         implicit none
!         real(kind=8), intent(in) :: istr,jstr,iend,jend
!         integer :: i, j, k
!         real Wrk(0:nz)

!         do j=jstr,jend
!           do i=istr,iend
!              Wrk(0) = 0.D0
!              do k=1,nz,+1 ! --> recursive
!                  Wrk(k) = Wrk(k-1) -Hz(i,j,k)
!     &       *( pm(i,j)*(ust(i+1,j,k)-ust(i,j,k))
!     &         +pn(i,j)*(vst(i,j+1,k)-vst(i,j,k)) )
!             enddo 
!           wst(i,j,nz)=+0.375*Wrk(nz) +0.75*Wrk(nz-1) -0.125*Wrk(nz-2)
!           do k=nz-1,2,-1
!              wst(i,j,k)=+0.5625*(Wrk(k  )+Wrk(k-1))
!     &                 -0.0625*(Wrk(k+1)+Wrk(k-2))
!           enddo
!           wst(i,j,1)= -0.125*Wrk(2) +0.75*Wrk(1) +0.375*Wrk(0)
!          enddo
!         enddo
!# ifndef EW_PERIODIC
!         if (WESTERN_EDGE) then                       ! Set lateral
!           do k=1,nz                                   ! boundary
!             do j=jstr,jend                           ! conditions
!               wst(istr-1,j,k)=wst(istr,j,k)
!           enddo
!          enddo
!         endif
!        if (EASTERN_EDGE) then
!          do k=1,nz
!            do j=jstr,jend
!              wst(iend+1,j,k)=wst(iend,j,k)
!            enddo
!          enddo
!        endif
!# endif
!# ifndef NS_PERIODIC
!         if (SOUTHERN_EDGE) then
!           do k=1,nz
!             do i=istr,iend
!               wst(i,jstr-1,k)=wst(i,jstr,k)
!             enddo
!          enddo
!         endif
!         if (NORTHERN_EDGE) then
!           do k=1,nz
!             do i=istr,iend
!               wst(i,jend+1,k)=wst(i,jend,k)
!             enddo
!           enddo
!         endif
!#  ifndef EW_PERIODIC
!         if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
!           do k=1,nz
!             wst(istr-1,jstr-1,k)=wst(istr,jstr,k)
!           enddo
!         endif
!         if (WESTERN_EDGE .and. NORTHERN_EDGE) then
!           do k=1,nz
!             wst(istr-1,jend+1,k)=wst(istr,jend,k)
!          enddo
!         endif
!         if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
!           do k=1,nz
!             wst(iend+1, jstr-1,k)=wst(iend,jstr,k)
!           enddo
!         endif
!         if (EASTERN_EDGE .and. NORTHERN_EDGE) then
!           do k=1,nz
!             wst(iend+1,jend+1,k)=wst(iend,jend,k)
!           enddo
!         endif
!#  endif
!# endif
   
!        end subroutine wstokes_tile
        ! ----------------------------------------
        
      end module wec_frc
