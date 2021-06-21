      module diagnostics

      ! INFO: ![
      ! Online diagnostics for ROMS simulation.
      ! Initial coding DevinD and JeroenM - 2021/05

      ! To compute diagnostics use DIAGNOSTICS flag.
      ! To chose between history slices of diagnostics or averages,
      ! set diag_avg to .false. for history, or .true. for averaging.
      ! Averaging is much less efficient as diagnostics must be collected
      ! at everytime step to maintain a true average.

      ! The main subroutine is wrt_online_diags, which is called at the very
      ! beginning of main.F, in order to record the diagnostics of 'nnew' index,
      ! which was t=n+1 for the previous timestep, and is now t=n since roms_step
      ! has moved onto the next timestep. I.e. 'nnew' will become 'nstp' in next timestep.
      ! But wrt_online_diags is called before nnew and nstp indicies change in main.F.

      ! If you are needing averages then the diagnotics need to be computed at everytimestep.
      ! If you just need snapshots then the diagnostics only need to be calculated at
      ! your output frequency.
      ! Both of these are controlled by the logical 'calc_diag', which is set in
      ! wrt_online_diags at the beginning of a timestep, and then that logical is
      ! used in the equation system to trigger a calculation of diagnostic.

      ! Units:
      ! The roms history file give u/v in m/s, and 't' tracers in their
      ! concentrations. However, since the depth of the water column is also evolving
      ! during a timestep, it is easier to work in dz*u/v/w/t, since that is the unit
      ! of u/v/w/t when the depth is changed in the timestep as follows:
      ! step3d_uv1.F (depth not changed just to show position amongst code)
      ! step2d_tile -> set_depth_tile -> (Hz changes)
      ! step3d_uv2.F (depth not changed just to show position amongst code)
      !
      ! See Documentation/readme-specifics/README.diagnostics_u_v_value_changes
      ! for more details on how the units of 'u' change during a timestep.

      ! It is assumed that by using logicals that are set as parameters
      ! the if branch will be removed during compiling, since the outcome
      ! of the branch is always known and unchanging. If this is true it will allow
      ! for if statements to be placed within loops, which is usually not advised
      ! due to the cost of the if statment at every index iteration.

      ! We should move some parameters into .in file, as it is way more convenient,
      ! but a template of what should be in the .in file when using this module should be
      ! left at the top here, so that it can be copied into the .in, you don't need
      ! to find an old diagnostics .in file...

      ! It is my feeling that once we use averages, all diagnostic calculations should
      ! not be calculated in the step3d... loop, but in it's own subroutine to avoid
      ! having if statments in loops. This will regardless result in the cleanest code.

      ! Mdiag(2): hydrostatic pressure gradient  ! OLD JM NOTES
      ! Mdiag(3): cor(u)
      ! Mdiag(4): adv(u) (non-dissipative part)
      ! Mdiag(5): horizontal mixing (including numerical dissipation)
      ! Mdiag(6): vertical mixing   (including numerical dissipation)
      ! Mdiag(7): non-hydrostatic pressure gradient

      ! Tdiag(1): dt/dt  ! OLD JM NOTES
      ! Tdiag(2): adv(t) (non-dissipative part)
      ! Tdiag(3): adv(t) (dissipative part)
      ! Tdiag(4): vertical mixing
      ! Tdiag(5): sources/sinks
      !]

#include "cppdefs.h"
#ifdef DIAGNOSTICS

      use dimensions
      use read_write
      use netcdf
      use tracers

      implicit none
      private

#include "ocean2d.h"
#include "ocean3d.h"
! scalars.h for 'knew' var
#include "scalars.h"

      ! ***************************************************************
      ! User inputs
      logical, parameter :: write_file    = .true.  ! whether to write a file or not (probably remove this)
      logical, parameter :: diag_avg      = .true. ! compute history (F) or averages (T)
      integer            :: rec_rate      = 2       ! output record rate, move this to .in file rather?
      integer            :: recs_per_file = 20       ! total recs per file
      ! End of user inputs
      ! ***************************************************************

      ! Preamble:  ![

      logical, public :: calc_diag  ! flag to calculate diagnostics in equation system.
      logical, public :: set_dz_old ! set dz*u/v/w/t at previous timestep, needed for diag calculation when skipping steps if not averages.

      ! U & V MOMENTUM:
      ! ---------------

      ! NOTE: u & v are 1:N in z, but w is 0:N hence different variables needed

      integer,public,parameter                            :: nmd = 7    ! number of momentum diag terms

      real, public, dimension(GLOBAL_2D_ARRAY, N, 2, nmd) :: Mdiag=0    ! 2 for u & v
      real, public, dimension(GLOBAL_2D_ARRAY, N, 2)      :: uvdz_old=0 ! u from previous timestep, need to either store this or Hz to get old u*dz.
      real, public, dimension(GLOBAL_2D_ARRAY, N, 2)      :: uv_dif=0   ! gross u/v difference between current vs previous timestep u/v
      real, public, dimension(GLOBAL_2D_ARRAY, N)         :: Mdsum=0    ! sum of Mdiag
      real, public, dimension(GLOBAL_2D_ARRAY, N)         :: diff=0     ! = uv_dif - Mdsum
      real, public, dimension(GLOBAL_2D_ARRAY, N, 2)      :: uv_prev=0  ! store incremental u/v/w (nnew)

      real, allocatable, dimension(:, :, :, :, :)         :: Md_avg     ! averages of diagnostics. Allocatabled (in subroutine) hence memory only if diag_avg is true.
      real, allocatable, dimension(:, :, :, :)            :: Md_a_sum   ! sum of diag averages at timestep
      real, allocatable, dimension(:, :, :, :)            :: Md_a_ever  ! total of all diag averages for entire simulation. Should add up to change in u/v/t over simulation...
      real, allocatable, dimension(:, :, :, :)            :: uvdz_t0    ! original dz*u/v to compare total change in u/v over simulation.
      real, allocatable, dimension(:, :, :, :)            :: uvdz_ever  ! total change in dz*u/v over simulation.

      real, public                                        :: Mdtmp ! temporary scalar for calculations of Mdiag

      real, public, dimension(GLOBAL_1DY_ARRAY, N, 2)     :: uv_wbry    ! bry change for u and v

      integer, parameter, public ::
     &  iprsgr = 1, ! internal loop ranges only (i.e. excluding domain boundary)
     &  icori  = 2, ! internal loop ranges only
     &  iadv   = 3, ! internal loop ranges only
     &  ivisc  = 4, ! internal loop ranges only
     &  iuv2DC = 5, ! internal loop ranges only
     &  iu3d   = 6, ! bry only
     &  iuv2FC = 7  ! internal + bry (full domain) loop ranges

      character (len=20), dimension(nmd) :: vnames = (/
     &  'prsgrd',
     &  'coriolis',
     &  'advection',
     &  'visc3d',
     &  'vert_mix_2D3D',
     &  '3dbc',       ! should store 3dbc changes in 2D arrays not 3D, as just 0 in interior
     &  '2D3D_cor'
     &                                                /)

      character (len=70), dimension(nmd) :: lname =   (/               ! len=70 allows for 8 characters of 'Average ' if diag_avg
     &  'Hydrostatic pressure gradient - prsgrd.F',
     &  'Coriolis & curvilinear coords - compute_horiz_rhs_uv_term.h',
     &  'Advection - step3d_uv1.F',
     &  'Viscous terms - visc3d_GP.F',
     &  'Vertical mixing & 2D/3D coupling - step3d_uv2.F part 1',
     &  'Boundary changes - u/v3dbc_im.F',
     &  '2D-3D correction - step3d_uv2.F part 2'
     &                                                 /)

      character (len=20)              :: units = 'm^2/s (dz*(u or v))'

      character (len=2), dimension(3) :: uvw = (/ 'u_', 'v_', 'w_' /)

      ! TRACERS:
      ! --------
# ifdef DIAGS_TRC
      ! Since there could be many tracers, it is preferable not to carry any 't_old'
      ! or temporary arrays for each tracer to avoid the memory burden. This differs
      ! from the u,v,w implementation here currently.
      ! Need to consider pipe and river sources.

      integer, public, parameter                           :: ndt = 5  ! number of tracer diag terms
#  ifdef LMD_KPP
     &                                                            + 1
#  endif
      real, public, dimension(GLOBAL_2D_ARRAY, N, NT, ndt) :: Tdiag=0    ! consider ndt, NT order

      real, public, dimension(GLOBAL_2D_ARRAY, N)          :: T_dsum=0   ! sum of t diagnostics to compare with gross change
      real, public, dimension(GLOBAL_2D_ARRAY, N)          :: Hz_old=0   ! store previous Hz to calculate old dz*t
                                                                       ! better than storing dz*t for each tracer. Better to recalculate it?
      integer, parameter, public ::
     &  it_h_adv    = 1,  ! compute_horiz_tracer_fluxes.h
     &  it_v_adv    = 2,  ! compute_vert_tracer_fluxes.h
      ! Do we need surface forcing flux as a vertical? Shouldn't be a 3D array anyway
     &  it_imp_vdif = 3,  ! step3d_t_ISO.F
     &  it_t3dmix   = 4,
     &  it_lmdkpp   = 5,  ! not used if no LMD_KPP flag, hence at the end. ! step3d_t_ISO.F. This only change T & S so shouldn't be 'NT' in size, a waste.
     &  it_bc       = 6   ! BC's shouldn't have 3D array

      character (len=30), dimension(ndt) :: tdname = (/
     &  '_hori_adv',
     &  '_vert_adv',
     &  '_imp_vdif',
     &  '_t3dmix',
     &  '_lmdkpp',
     &  '_t3dbc'
     &                                               /)

      character (len=80), dimension(ndt) :: tlname = (/
     &  'Horizontal advective fluxes - compute_horiz_tracer_fluxes.h',
     &  'Vertical advective fluxes - compute_vert_tracer_fluxes.h.h',
     &  'Implicit time step for vertical diffusion - step3d_t_ISO.F',
     &  'Horizontal (along geopotential surfaces) Laplacian diffusion terms - t3dmix_GP.F',
     &  'Nonlocal transport flux for unstable (convective) forcing - step3d_t_ISO.F',
     &  'Boundary changes - t3dbc_im.F'
     &                                               /)

      character (len=20) :: tunits = 'C.m^2/s (dz*t)'  ! not sure on units here.

      public set_diags_t_at_bc
      public set_diags_t_at_t3dmix_GP
# endif /* DIAGS_TRC */

      ! Netcdf outputting:
      ! ------------------
      integer            :: output_rec = 0              ! record number of output. 0 indicates we need new file
      integer            :: total_output_rec = 0        ! total flux output recordings so far
      integer            :: tn                          ! roms timestep
      character (len=60) :: file_name = '_on_diags.nc'  ! suffix to append to simulation output rootname


      ! W MOMENTUM:
      ! -----------
# ifdef NHMG
      integer,public,parameter :: nwd = 5  !! number of 'w' momentum diag terms (for NHMG runs)

      real, public, dimension(GLOBAL_2D_ARRAY, 0:N, nwd) :: Wdiag=0    ! 3 for u,v,w
      real, public, dimension(GLOBAL_2D_ARRAY, 0:N)      :: wdz_old=0 ! u from previous timestep, need to either store this or Hz to get old u*dz.
      real, public, dimension(GLOBAL_2D_ARRAY, 0:N)      :: w_dif=0   ! gross u/v difference between current vs previous timestep u/v
      real, public, dimension(GLOBAL_2D_ARRAY, 0:N)      :: Wdsum=0    ! sum of Mdiag
      real, public, dimension(GLOBAL_2D_ARRAY, 0:N)      :: wdif=0     ! = w_dif - Mdsum
      real, public, dimension(GLOBAL_2D_ARRAY, 0:N)      :: w_prev=0  ! store incremental u/v/w (nnew)

      integer, parameter, public ::
     &  iwprsgr   = 1,
     &  iwhoriadv = 2,
     &  iwvertadv = 3,
     &  iwuv2     = 4,
     &  iwbc      = 5

      character (len=20), dimension(nwd) :: wvname = (/
     &  'w_prsgrd',
     &  'w_horiz_flx',
     &  'w_vert_adv',
     &  'w_2D_3D_coupling',
     &  'w_3dbc'
     &                                               /)

      character (len=60), dimension(nwd) :: wlname = (/
     &  'prsgrd.F',
     &  'compute_horiz_rhs_w_terms.h',
     &  'compute_vert_rhs_w_terms.h',
     &  'step3d_uv2.F',
     &  'w3dbc_im.F'
     &                                               /)

      character (len=20) :: wunits = 'm^2/s (dz*w)'  ! not sure on units here.

      ! Public subroutines:
      public set_diags_w_at_uv1
      public set_diags_w_at_uv2
      public set_diags_w_at_uv2_end
      public set_diags_w_at_bc
# endif /* NHMG */

      public init_diags
      public set_diags_at_visc3d
      public set_diags_at_uv1
      public set_diags_at_uv2_DC
      public set_diags_at_bc_u
      public set_diags_at_bc_v
      public set_diags_at_uv2_FCu
      public set_diags_at_uv2_FCv
      public set_diags_dz_old
      public wrt_diagnostics !]

      contains

! ----------------------------------------------------------------------
      subroutine init_diags(istr, jstr, istrR, iendR, jstrR, jendR)  ![
      ! Set original dz*u/v/w/t for initial timestep, hence only called once.
      ! Called from prestep before Hz has changed in the timestep.

      implicit none

      ! inputs
      integer, intent(in) :: istr, jstr, istrR, iendR, jstrR, jendR

      ! local
      integer :: i,j,k

      do k=1,N
        do j=jstrR,jendR   ! Full loop ranges since u-change over every point including bry.
          do i=istr,iendR  ! istr is the west boundary for u-points, istrR=E, jstrR=S, jendR=N
            uvdz_old(i,j,k, 1) = u(i,j,k,nstp) * (0.5*( Hz(i-1,j,k)+Hz(i,j,k) ))
          enddo
        enddo
        do j=jstr,jendR    ! Full loop ranges since u-change over every point including bry.
          do i=istrR,iendR ! istr is the west boundary for u-points, istrR=E, jstrR=S, jendR=N
            uvdz_old(i,j,k, 2) = v(i,j,k,nstp) * (0.5*( Hz(i,j-1,k)+Hz(i,j,k) ))
          enddo
        enddo
      enddo

      if (FIRST_TIME_STEP .and. diag_avg) uvdz_t0 = uvdz_old  ! store the initial condition for budget confirmation at simulation end

# ifdef DIAGS_TRC
      Hz_old = Hz ! store Hz to calculate old dz*t
# endif
# ifdef NHMG
      do k=1,N ! 0 is always 0?
        do j=jstrR,jendR   ! Full loop ranges since u-change over every point including bry.
          do i=istrR,iendR
            if(k==N) then
              wdz_old(i,j,k) = w(i,j,k,nstp) * Hz(i,j,k)
            else
              wdz_old(i,j,k) = w(i,j,k,nstp) * (0.5*( Hz(i,j,k+1)+Hz(i,j,k) ))
            endif
          enddo
        enddo
      enddo
# endif /* NHMG */

      if(mynode==0) print *, 'INIT HZ_OLD **************'  ! debug

      end subroutine init_diags  !]

! ----------------------------------------------------------------------
      subroutine set_diags_at_uv1( istr, iend, jstr, jend, istrA, j, dim, DC, ruv, suvstr, uv )  ![
      ! Set the diagnostic terms for prsgrd, coriolis, and advection.
      ! Works for both u and v.
      ! Prsgrd and coriolis need to also be convert to dz * u/v

      implicit none

      ! inputs
      integer, intent(in) :: istr, iend, jstr, jend, istrA, j, dim ! istr, iend, etc. 'A' for any as depends on inputs, which 'dim' dimension.
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N), intent(in) :: DC
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,N),   intent(in) :: ruv     ! ru or rv
      real, dimension(GLOBAL_2D_ARRAY),              intent(in) :: suvstr  ! sustr or svstr
      real, dimension(GLOBAL_2D_ARRAY,N,3),          intent(in) :: uv      ! u or v

      ! local
      integer :: i, k

      ! NOTE: currently missing WEC term

      do i=istrA,iend
        Mdiag(i,j,N,dim,iprsgr) = Mdiag(i,j,N,dim,iprsgr) * DC(i,0)  ! need to convert to dz*u now that we have DC
        Mdiag(i,j,1,dim,iprsgr) = Mdiag(i,j,1,dim,iprsgr) * DC(i,0)  ! consider taking this out of loop and doing seperately
        Mdiag(i,j,N,dim,icori)  = Mdiag(i,j,N,dim,icori)  * DC(i,0)  ! perhaps even in own subroutine to avoid clutter
        Mdiag(i,j,1,dim,icori)  = Mdiag(i,j,1,dim,icori)  * DC(i,0)
        Mdiag(i,j,N,dim,iadv)   = DC(i,0)*ruv(i,j,N) + dt*suvstr(i,j)                 ! all advection - includes 2nd part of compute_horiz_uv & compute_vert_uv.
     &                           - ( Mdiag(i,j,N,dim,iprsgr) + Mdiag(i,j,N,dim,icori) ) ! subtract previous 2 terms included in total
        uv_prev(i,j,N,dim)      = uv(i,j,N,nnew)     ! store for visc3d
        Mdiag(i,j,1,dim,iadv)   = DC(i,0)*ruv(i,j,1)
     &                           - ( Mdiag(i,j,1,dim,iprsgr) + Mdiag(i,j,1,dim,icori) )
        uv_prev(i,j,1,dim)      = uv(i,j,1,nnew)     ! store for visc3d
      enddo

      do k=2,N-1
        do i=istrA,iend
          Mdiag(i,j,k,dim,iprsgr) = Mdiag(i,j,k,dim,iprsgr) * DC(i,0)
          Mdiag(i,j,k,dim,icori)  = Mdiag(i,j,k,dim,icori)  * DC(i,0)
          uv_prev(i,j,k,dim)      = uv(i,j,k,nnew)
          Mdiag(i,j,k,dim,iadv)   = DC(i,0)*ruv(i,j,k) ! all advection - includes 2nd part of compute_horiz_uv & compute_vert_uv.
     &                            - ( Mdiag(i,j,k,dim,iprsgr) + Mdiag(i,j,k,dim,icori) )
        enddo
      enddo

      end subroutine set_diags_at_uv1  !]

! ----------------------------------------------------------------------
      subroutine set_diags_at_visc3d( istr, iend, jstr, jend, isA, jsA, dim, uv )  ![
      ! set diagnostics for u or v in visc3d_GP.F
      ! interior points only (i.e. excluding bry)

      implicit none

      ! inputs
      integer, intent(in) :: istr, iend, jstr, jend, isA, jsA, dim  ! isA jsA are starts depending on u or v variable
      real, dimension(GLOBAL_2D_ARRAY,N,3), intent(in) :: uv        ! u or v
      integer :: i, j, k  ! local

      do k=1,N
        do j=jsA,jend    ! loop ranges for interior of subdomains, exludes global domain boundaries
          do i=isA,iend
            Mdiag( i,j,k,dim,ivisc) = uv(i,j,k,nnew) - uv_prev(i,j,k,dim)  ! units already in dz*u
            uv_prev(i,j,k,dim)      = uv(i,j,k,nnew)
          enddo
        enddo
      enddo
!      Mdiag (istrU:iend,jstr:jend,:,1,ivisc) = u(istrU:iend,jstr:jend,:,nnew)
!      uv_prev(istrU:iend,jstr:jend,:)         = u(istrU:iend,jstr:jend,:,nnew) ! can remove

      end subroutine set_diags_at_visc3d  !]

! ----------------------------------------------------------------------
      subroutine set_diags_at_uv2_DC( istr, iend, jstr, jend, istrU, jstrV )  ![
      ! set diagnostics for u or v in step3d_uv2.F
      ! interior points only (i.e. excluding bry)

      implicit none

      ! inputs
      integer, intent(in) :: istr, iend, jstr, jend, istrU, jstrV
      integer :: i, j, k  ! local

      do k=1,N
        do j=jstr,jend    ! loop ranges for interior of subdomains, exludes global domain boundaries
          do i=istrU,iend
            Mdtmp                  = u(i,j,k,nnew) * (0.5*( Hz(i-1,j,k)+Hz(i,j,k) )) ! convert to dz*u
            Mdiag( i,j,k,1,iuv2DC) = Mdtmp - uv_prev(i,j,k,1)  ! avoids recomputing above
            uv_prev(i,j,k,1)       = Mdtmp
          enddo
        enddo
        do j=jstrV,jend    ! loop ranges for interior of subdomains, exludes global domain boundaries
          do i=istr,iend
            Mdtmp                  = v(i,j,k,nnew) * (0.5*( Hz(i,j-1,k)+Hz(i,j,k) )) ! convert to dz*v
            Mdiag( i,j,k,2,iuv2DC) = Mdtmp - uv_prev(i,j,k,2)
            uv_prev(i,j,k,2)       = Mdtmp
          enddo
        enddo
      enddo

      if(mynode==0) print *, 'uv2DC!!!!!!!!!'  ! debug only

      end subroutine set_diags_at_uv2_DC  !]

! ----------------------------------------------------------------------
      subroutine set_diags_at_bc_u( istr, iend, jstr, jend )  ![
      ! set the diagnostics terms for boundary changes in u
      ! called from step3d_uv2.F right after call to u3dbc_im.F

      implicit none

      ! inputs/outputs
      integer, intent(in) :: istr, iend, jstr, jend

      ! local
      integer :: i

      if (WESTERN_EDGE)  uv_prev(istr  ,:,:,1) = u(istr  ,:,:,nnew) * 0.5*(Hz(istr-1,:,:)+Hz(istr  ,:,:))
      if (EASTERN_EDGE)  uv_prev(iend+1,:,:,1) = u(iend+1,:,:,nnew) * 0.5*(Hz(iend  ,:,:)+Hz(iend+1,:,:))

      if (WESTERN_EDGE)  uv_wbry(:,:,1)           = uv_prev(istr  ,:,:,1) - uvdz_old(istr  ,:,:, 1)

      if (WESTERN_EDGE)  Mdiag(istr  ,:,:,1,iu3d) = uv_prev(istr  ,:,:,1) - uvdz_old(istr  ,:,:, 1) ! uv_prev is already u_current, can do this because have uvdz_old already
      if (EASTERN_EDGE)  Mdiag(iend+1,:,:,1,iu3d) = uv_prev(iend+1,:,:,1) - uvdz_old(iend+1,:,:, 1)

      do i=istr,iend+1 ! use loop here else indexing would be even messier
        if (SOUTHERN_EDGE) uv_prev(i,jstr-1,:,1)=u(i,jstr-1,:,nnew)*0.5*(Hz(i-1,jstr-1,:)+Hz(i,jstr-1,:))
        if (NORTHERN_EDGE) uv_prev(i,jend+1,:,1)=u(i,jend+1,:,nnew)*0.5*(Hz(i-1,jend+1,:)+Hz(i,jend+1,:))
      enddo
      if (SOUTHERN_EDGE) Mdiag(:,jstr-1,:,1,iu3d) = uv_prev(:,jstr-1,:,1) - uvdz_old(:,jstr-1,:, 1)
      if (NORTHERN_EDGE) Mdiag(:,jend+1,:,1,iu3d) = uv_prev(:,jend+1,:,1) - uvdz_old(:,jend+1,:, 1)

      end subroutine set_diags_at_bc_u  !]

! ----------------------------------------------------------------------
      subroutine set_diags_at_bc_v( istr, iend, jstr, jend )  ![
      ! set the diagnostics terms for boundary changes in v
      ! called from step3d_uv2.F right after call to v3dbc_im.F

      implicit none

      ! inputs/outputs
      integer, intent(in) :: istr, iend, jstr, jend

      ! local
      integer :: j

      do j=jstr,jend+1 ! use loop here else indexing would be even messier
        if (WESTERN_EDGE)  uv_prev(istr-1,j,:,2) = u(istr-1,j,:,nnew) * 0.5*(Hz(istr-1,j-1,:)+Hz(istr-1,j,:))
        if (EASTERN_EDGE)  uv_prev(iend+1,j,:,2) = u(iend+1,j,:,nnew) * 0.5*(Hz(iend+1,j-1,:)+Hz(iend+1,j,:))
      enddo
      if (WESTERN_EDGE)  Mdiag(istr-1,:,:,2,iu3d) = uv_prev(istr-1,:,:,2) - uvdz_old(istr-1,:,:,2) ! uv_prev is already u_current, can do this because have uvdz_old already
      if (EASTERN_EDGE)  Mdiag(iend+1,:,:,2,iu3d) = uv_prev(iend+1,:,:,2) - uvdz_old(iend+1,:,:,2)

      if (SOUTHERN_EDGE) uv_prev(:,jstr  ,:,2)=u(:,jstr  ,:,nnew)*0.5*(Hz(:,jstr-1,:)+Hz(:,jstr  ,:))
      if (NORTHERN_EDGE) uv_prev(:,jend+1,:,2)=u(:,jend+1,:,nnew)*0.5*(Hz(:,jend  ,:)+Hz(:,jend+1,:))
      if (SOUTHERN_EDGE) Mdiag(:,jstr  ,:,2,iu3d) = uv_prev(:,jstr  ,:,2) - uvdz_old(:,jstr  ,:,2)
      if (NORTHERN_EDGE) Mdiag(:,jend+1,:,2,iu3d) = uv_prev(:,jend+1,:,2) - uvdz_old(:,jend+1,:,2)

      end subroutine set_diags_at_bc_v  !]

! ----------------------------------------------------------------------
      subroutine set_diags_at_uv2_FCu( istr, iendR, jstrR, jendR )  ![
      implicit none

      ! inputs
      integer, intent(in) :: istr, iendR, jstrR, jendR ! changed the name to use for u or v
      integer :: i, j, k  ! local

      do k=1,N
        do j=jstrR,jendR  ! since this section uv2FC uses IU_RANGE & J_RANGE which includes boundary
          do i=istr,iendR
            Mdiag(i,j,k,1,iuv2FC) = u(i,j,k,nnew) * (0.5*( Hz(i-1,j,k)+Hz(i,j,k) )) - uv_prev(i,j,k,1)

            ! Full loop ranges since u-change over every point including bry.
            ! istr is the west boundary for u-points, istrR=E, jstrR=S, jendR=N
            ! can't just do (u(nnew)-u(nstp)*Hz here because actual need to old Hz for u(nstp) as Hz has evolved.
            uv_dif(i,j,k, 1)   = u(i,j,k,nnew) * (0.5*( Hz(i-1,j,k)+Hz(i,j,k) )) - uvdz_old(i,j,k, 1) ! do this here as units are correct
          enddo
        enddo
      enddo

      end subroutine set_diags_at_uv2_FCu  !]

! ----------------------------------------------------------------------
      subroutine set_diags_at_uv2_FCv( istrR, iendR, jstr, jendR )  ![
      implicit none

      ! inputs
      integer, intent(in) :: istrR, iendR, jstr, jendR ! changed the name to use for u or v
      integer :: i, j, k  ! local

      do k=1,N  ! COULD REUSE THIS BY LOOPING OVER a AND b AND SWITCHING ARGUEMENTS
        do j=jstr,jendR  ! since this section uv2FC uses IU_RANGE & J_RANGE which includes boundary
          do i=istrR,iendR
            Mdiag(i,j,k,2,iuv2FC) = v(i,j,k,nnew) * (0.5*( Hz(i,j-1,k)+Hz(i,j,k) )) - uv_prev(i,j,k,2)
            ! put = Mdiag below rather than same equation...
            ! Full loop ranges since u-change over every point including bry.
            ! istr is the west boundary for u-points, istrR=E, jstrR=S, jendR=N
            ! can't just do (u(nnew)-u(nstp)*Hz here because actual need to old Hz for u(nstp) as Hz has evolved.
            uv_dif(i,j,k, 2)   = v(i,j,k,nnew) * (0.5*( Hz(i,j-1,k)+Hz(i,j,k) )) - uvdz_old(i,j,k, 2) ! do this here as units are correct
          enddo
        enddo
      enddo

      end subroutine set_diags_at_uv2_FCv  !]

! ----------------------------------------------------------------------
      subroutine set_diags_dz_old(istr, jstr, istrR, iendR, jstrR, jendR)  ![
      implicit none

      ! inputs
      integer, intent(in) :: istr, jstr, istrR, iendR, jstrR, jendR

      ! local
      integer :: i, j, k

      ! can just loop both over istrR and jstrR and have one loop. No harm?
      ! This is currently used for bry calculation so is also needed for averaging. u_diff is other use but that's not needed for averaging.
      ! Should put an if statement here for diags_avg and only compute bry if averages...

      do k=1,N
        do j=jstrR,jendR
          do i=istr,iendR
            uvdz_old(i,j,k, 1) = u(i,j,k,nnew) * (0.5*( Hz(i-1,j,k)+Hz(i,j,k) ))
          enddo
        enddo
      enddo

      do k=1,N
        do j=jstr,jendR
          do i=istrR,iendR
            uvdz_old(i,j,k, 2) = v(i,j,k,nnew) * (0.5*( Hz(i,j-1,k)+Hz(i,j,k) ))
          enddo
        enddo
      enddo

      ! Don't forget 'w' here.

      if(mynode==0) print *, 'DZ_OLD  $$$$$$$$$$'  ! debug

      end subroutine set_diags_dz_old  !]

! ----------------------------------------------------------------------
# ifdef NHMG
      subroutine set_diags_w_at_uv1( istr, iend, jstr, jend, DC, rw )  ![
      ! Set the diagnostic terms for prsgrd, horiz & vert advection.
      ! Prsgrd and horiz adv still need to also be convert to dz * w

      implicit none

      ! inputs
      integer,                                       intent(in) :: istr, iend, jstr, jend
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N), intent(in) :: DC
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,0:N), intent(in) :: rw

      ! local
      integer :: i, j, k

      do k=1,N
        do j=jstr,jend
          do i=istr,iend
            Wdiag(i,j,k,iwprsgr)   = Wdiag(i,j,k,iwprsgr)   * DC(i,0)  ! need to convert to dz*w now that we have DC
            Wdiag(i,j,k,iwhoriadv) = Wdiag(i,j,k,iwhoriadv) * DC(i,0)
            Wdiag(i,j,k,iwvertadv) = DC(i,0)*rw(i,j,k)
     &                              - ( Wdiag(i,j,k,iwprsgr) + Wdiag(i,j,k,iwhoriadv) ) ! subtract previous 2 terms included in total
            w_prev(i,j,k)          = w(i,j,k,nnew)     ! store for next change
          enddo
        enddo
      enddo

      end subroutine set_diags_w_at_uv1  !]

! ----------------------------------------------------------------------
      subroutine set_diags_w_at_uv2( istr, iend, jstr, jend )  ![
      implicit none

      ! inputs
      integer, intent(in) :: istr, iend, jstr, jend

      ! local
      integer :: i, j, k

      ! This could well be wrong. Need to confirm with Jeroen.
      do k=1,N-1
        do j=jstr,jend  ! I changed istr to jstr...
          do i=Istr,Iend
            Wdiag(i,j,k,iwuv2) = w(i,j,k,nnew) * 0.5*( Hz(i,j,k+1)+Hz(i,j,k) ) - w_prev(i,j,k)
          enddo
        enddo
      enddo

      ! Need to do N seperately since can't do 0.5 * Hz(i,j,k+1)
      ! Do seperate loop rather than having if statement in loop? Not sure what is more efficient
      ! ifort compiler report could tell me if it fuses these loops...
      do j=jstr,jend  ! I changed istr to jstr...
        do i=Istr,Iend
          Wdiag(i,j,N,iwuv2) = w(i,j,N,nnew) * Hz(i,j,N) - w_prev(i,j,N)
        enddo
      enddo

      end subroutine set_diags_w_at_uv2  !]

! ----------------------------------------------------------------------
      subroutine set_diags_w_at_bc( istr, iend, jstr, jend )  ![
      ! set the diagnostics terms for boundary changes in w
      ! called from step3d_uv2.F right after call to w3dbc_im.F

      implicit none

      ! inputs/outputs
      integer, intent(in) :: istr, iend, jstr, jend

      ! local
      integer :: i, j, k

      if (WESTERN_EDGE) then
        i=istr-1
        do k=1,N-1            ! N-1 because of k+1
          do j=jstr-1,jend+1
            Wdiag(i,j,k,iwbc) = w(i,j,k,nnew) * 0.5*(Hz(i,j,k+1)+Hz(i,j,k)) - wdz_old(i,j,k) ! w_prev is already u_current, can do this because have uvdz_old already
          enddo
        enddo
        Wdiag(i,j,N,iwbc) = w(i,j,N,nnew) * Hz(i,j,N) - wdz_old(i,j,N)
      endif

      if (EASTERN_EDGE) then
        i=iend+1
        do k=1,N-1
          do j=jstr-1,jend+1
            Wdiag(i,j,k,iwbc) = w(i,j,k,nnew) * 0.5*(Hz(i,j,k+1)+Hz(i,j,k)) - wdz_old(i,j,k) ! w_prev is already u_current, can do this because have uvdz_old already
          enddo
        enddo
        Wdiag(i,j,N,iwbc) = w(i,j,N,nnew) * Hz(i,j,N) - wdz_old(i,j,N)
      endif

      if (SOUTHERN_EDGE) then
        j=jstr-1
        do k=1,N-1
          do i=istr-1,iend+1
            Wdiag(i,j,k,iwbc) = w(i,j,k,nnew) * 0.5*(Hz(i,j,k+1)+Hz(i,j,k)) - wdz_old(i,j,k) ! w_prev is already u_current, can do this because have uvdz_old already
          enddo
        enddo
        Wdiag(i,j,N,iwbc) = w(i,j,N,nnew) * Hz(i,j,N) - wdz_old(i,j,N)
      endif

      if (NORTHERN_EDGE) then
        j=jend+1
        do k=1,N-1
          do i=istr-1,iend+1
            Wdiag(i,j,k,iwbc) = w(i,j,k,nnew) * 0.5*(Hz(i,j,k+1)+Hz(i,j,k)) - wdz_old(i,j,k) ! w_prev is already u_current, can do this because have uvdz_old already
          enddo
        enddo
        Wdiag(i,j,N,iwbc) = w(i,j,N,nnew) * Hz(i,j,N) - wdz_old(i,j,N)
      endif

      end subroutine set_diags_w_at_bc  !]

! ----------------------------------------------------------------------
      subroutine set_diags_w_at_uv2_end( istrR, iendR, jstrR, jendR)  ![
      implicit none

      ! inputs
      integer, intent(in) :: istrR, iendR, jstrR, jendR ! changed the name to use for u or v
      integer :: i, j, k, tmp  ! local

      do k=1,N-1
        do j=jstrR,jendR
          do i=istrR,iendR
            ! Full loop ranges since u-change over every point including bry.
            tmp = w(i,j,k,nnew) * 0.5 * ( Hz(i,j,k+1)+Hz(i,j,k) )
            w_dif(i,j,k)   = tmp - wdz_old(i,j,k)
            wdz_old(i,j,k) = tmp
          enddo
        enddo
      enddo

      do j=jstrR,jendR
        do i=istrR,iendR
          tmp = w(i,j,N,nnew) * Hz(i,j,N)
          w_dif(i,j,N)   = tmp - wdz_old(i,j,N)
          wdz_old(i,j,N) = tmp
        enddo
      enddo

      end subroutine set_diags_w_at_uv2_end  !]
# endif /* NHMG */

! ----------------------------------------------------------------------
# ifdef DIAGS_TRC
      subroutine set_diags_t_at_bc(istr, iend, jstr, jend, itrc)  ![

      use tracers, only: t
      implicit none

      ! inputs/outputs
      integer, intent(in) :: istr, iend, jstr, jend, itrc

      ! local
      integer :: i, j

      ! PROBABLY DON'T NEED THIS CAN JUST USE T(NNEW) AT THE BOUNDARY SINCE ONLY UPDATE FROM HERE.

      if (WESTERN_EDGE) then
        i=istr-1
        Tdiag(i,:,:,itrc,it_bc) = t(i,:,:,nnew,itrc) *     Hz(i,:,:)
     &                          - t(i,:,:,nstp,itrc) * Hz_old(i,:,:)
      endif
      if (EASTERN_EDGE) then
        i=iend+1
        Tdiag(i,:,:,itrc,it_bc) = t(i,:,:,nnew,itrc) *     Hz(i,:,:)
     &                          - t(i,:,:,nstp,itrc) * Hz_old(i,:,:)
      endif
      if (SOUTHERN_EDGE) then
        j=jstr-1
        Tdiag(:,j,:,itrc,it_bc) = t(:,j,:,nnew,itrc) *     Hz(:,j,:)
     &                          - t(:,j,:,nstp,itrc) * Hz_old(:,j,:)
      endif
      if (NORTHERN_EDGE) then
        j=jend+1
        Tdiag(:,j,:,itrc,it_bc) = t(:,j,:,nnew,itrc) *     Hz(:,j,:)
     &                          - t(:,j,:,nstp,itrc) * Hz_old(:,j,:)
      endif

      end subroutine set_diags_t_at_bc  !]

! ----------------------------------------------------------------------
      subroutine set_diags_t_at_t3dmix_GP(istr,iend,jstr,jend)  ![

      use tracers, only: t
      implicit none           ! from t3dmix_GP.F

      ! inputs/outputs
      integer, intent(in) :: istr, iend, jstr, jend

      ! local
      integer :: i, j, k, itrc

      do itrc=1,NT ! CHANGE THIS TO ONLY tracers of interest!
        do k=1,N  ! K>0 in t3dmix
          do j=jstr,jend
            do i=istr,iend
              Tdiag(i,j,k,itrc,it_t3dmix) = t(i,j,k,nnew,itrc)
            enddo
          enddo
        enddo
      enddo

      end subroutine set_diags_t_at_t3dmix_GP  !]

! ----------------------------------------------------------------------
      subroutine set_diags_avg_tile(istr,iend,jstr,jend)  ![ set the averages of diagnostic variables
      implicit none

      ! MUST ONLY HAVE ONE TILE IN WRT_ONLINE_DIAGS THEN THIS WILL BECOME EASIER

      ! inputs
      integer, intent(in) :: istr,iend,jstr,jend

      ! local
      integer :: i, j, k, iMd

! Needed for istrU, istrR, etc
# include "compute_auxiliary_bounds.h"

      if(mynode==0) print *, '--- SETTING DIAGS AVERAGES' ! debug only

      ! This method is arguably more efficient than set_avg.F, however, it requires more
      ! lines of code. Since diagnostic averaging will be heavy, more efficient code is argueably
      ! more important.

      if (mod(tn,rec_rate)==0) then      ! timestep / rec_rate

!        do iMd=1,iuv2DC                 ! all diagnostic terms with internal loop ranges
!       ! for greatest efficiency should only do the loop ranges of the diag terms, but doing full
!       ! ranges for now just to get working.

        do iMd=1,nmd                     ! all diagnostic terms with internal loop ranges
          do k=1,N                       ! here we divide the average by the averaging rate
            do j=jstrR,jendR
              do i=istrR,iendR
                Md_avg(i,j,k,1,iMD) = ( Md_avg(i,j,k,1,iMD)+Mdiag(i,j,k,1,iMD) ) / rec_rate
              enddo
            enddo
          enddo
        enddo

        Md_a_sum = 0 ! improve later
        do iMd=1,nmd
          Md_a_sum(:,:,:,1) = Md_a_sum(:,:,:,1) + Md_avg(:,:,:,1,iMd)      ! sum all diag avgs at timestep
        enddo
        Md_a_ever = Md_a_ever + Md_a_sum * rec_rate                        ! keep an endless running total for budget closing
        call set_diags_uvdz_ever( istr, jstr, istrR, iendR, jstrR, jendR ) ! * rec_rate since tot = avg * t

        if(mynode==0) print *, '----> DIVIDE FOR AVERAGE'

      else

        if (mod(tn-1,rec_rate)==0) then   ! this is next step after we have averaged.

          ! or could have Md_avg = 0 at the top of routine but this is more efficient but more code?
          do iMd=1,nmd                    ! all diagnostic terms with internal loop ranges
            do k=1,N                      ! here we divide the average by the averaging rate
              do j=jstrR,jendR
                do i=istrR,iendR
                  Md_avg(i,j,k,1,iMD) = Mdiag(i,j,k,1,iMD)
                enddo
              enddo
            enddo
          enddo

        else                              ! all other timesteps

          do iMd=1,nmd                    ! all diagnostic terms with internal loop ranges
            do k=1,N                      ! here we divide the average by the averaging rate
              do j=jstrR,jendR
                do i=istrR,iendR
                  Md_avg(i,j,k,1,iMD) = Md_avg(i,j,k,1,iMD)+Mdiag(i,j,k,1,iMD)
                enddo
              enddo
            enddo
          enddo

        endif  ! <-- mod(tn-1,rec_rate)==0

      endif    ! <-- mod(tn  ,rec_rate)==0

      end subroutine set_diags_avg_tile  !]

! ----------------------------------------------------------------------
      subroutine set_diags_uvdz_ever( istr, jstr, istrR, iendR, jstrR, jendR )  ![
      ! need current dz*u/v to compare gross changes in u with sum of averages
      implicit none

      ! inputs
      integer, intent(in) :: istr, jstr, istrR, iendR, jstrR, jendR

      ! local
      integer :: i, j, k

      do k=1,N
        do j=jstrR,jendR
          do i=istr,iendR
            uvdz_ever(i,j,k, 1) = u(i,j,k,nnew) * (0.5*( Hz(i-1,j,k)+Hz(i,j,k) ))
     &                          - uvdz_t0(i,j,k,1)
          enddo
        enddo
      enddo

      do k=1,N
        do j=jstr,jendR
          do i=istrR,iendR
            uvdz_ever(i,j,k, 2) = v(i,j,k,nnew) * (0.5*( Hz(i,j-1,k)+Hz(i,j,k) ))
     &                          - uvdz_t0(i,j,k,2)
          enddo
        enddo
      enddo

      end subroutine set_diags_uvdz_ever  !]

! ----------------------------------------------------------------------
      subroutine wrt_tracer_diags_tile(istr,iend,jstr,jend, ncid)  ![
      ! since the value of the tracer has been stored as dz*t, we need to find the
      ! differences for the contribution of each diagnostic term

      use tracers, only: t
      implicit none

      ! inputs
      integer, intent(in)                 :: istr, iend, jstr, jend, ncid

      ! local
      integer                             :: itrc, idt, i, j, k, ierr, it_prev
      character(len=30)                   :: text
      real, dimension(GLOBAL_2D_ARRAY, N) :: Td_tmp=0  ! Move this to top of module for permanent scratch

! Needed for istrU, istrR, etc
# include "compute_auxiliary_bounds.h"

      do itrc=1,NT ! NT  ! CHANGE THIS TO ONLY tracers of interest!

        T_dsum=0 ! could just do boundaries as the rest can be overwritten in first below
        Td_tmp=0 ! Do more efficient reset of boundaries here. Above changes bry values but below doesn't use them but can be seen in output if not set to 0 again.

        do k=1,N            ! first diagnostic calculated using previous step's value
          do j=jstr,jend
            do i=istr,iend
              Td_tmp(i,j,k)=Tdiag(i,j,k,itrc,it_h_adv) - t(i,j,k,nstp,itrc)*Hz_old(i,j,k)
              T_dsum(i,j,k)=Td_tmp(i,j,k) ! add to total
            enddo
          enddo
        enddo
        text=trim(t_vname(itrc)) / / tdname(it_h_adv)
        call nc_write_var(ncid, Td_tmp, N, text, rp_var, output_rec, ierr)

        do k=1,N
          do j=jstr,jend
            do i=istr,iend
              Td_tmp(i,j,k)=Tdiag(i,j,k,itrc,it_v_adv) - Tdiag(i,j,k,itrc,it_h_adv)
              T_dsum(i,j,k)=T_dsum(i,j,k)+Td_tmp(i,j,k) ! add to total
            enddo
          enddo
        enddo
        text=trim(t_vname(itrc)) / / tdname(it_v_adv)
        call nc_write_var(ncid, Td_tmp, N, text, rp_var, output_rec, ierr)
        it_prev=it_v_adv ! store index of previous diagnostic

# ifdef LMD_KPP
        do k=1,N
          do j=jstr,jend
            do i=istr,iend
              Td_tmp(i,j,k)=Tdiag(i,j,k,itrc,it_lmdkpp) - Tdiag(i,j,k,itrc,it_v_adv)
              T_dsum(i,j,k)=T_dsum(i,j,k)+Td_tmp(i,j,k) ! add to total
            enddo
          enddo
        enddo
        text=trim(t_vname(itrc)) / / tdname(it_lmdkpp)
        call nc_write_var(ncid, Td_tmp, N, text, rp_var, output_rec, ierr)
        it_prev=it_lmdkpp ! store index of previous diagnostic if lmd_kpp
# endif

        ! From here onwards 't' was no longer in dz*t units:
        do k=1,N
          do j=jstr,jend
            do i=istr,iend  ! it_imp_vdif in 't' units, but it_prev in dz*t units
              Td_tmp(i,j,k)=Tdiag(i,j,k,itrc,it_imp_vdif)*Hz(i,j,k) - Tdiag(i,j,k,itrc,it_prev)
              T_dsum(i,j,k)=T_dsum(i,j,k)+Td_tmp(i,j,k) ! add to total
            enddo
          enddo
        enddo
        text=trim(t_vname(itrc)) / / tdname(it_imp_vdif)
        call nc_write_var(ncid, Td_tmp, N, text, rp_var, output_rec, ierr)

        do k=1,N
          do j=jstr,jend
            do i=istr,iend  ! both in dz*t units
              Td_tmp(i,j,k)=(Tdiag(i,j,k,itrc,it_t3dmix)-Tdiag(i,j,k,itrc,it_imp_vdif))*Hz(i,j,k)
              T_dsum(i,j,k)=T_dsum(i,j,k)+Td_tmp(i,j,k) ! add to total
            enddo
          enddo
        enddo
        text=trim(t_vname(itrc)) / / tdname(it_imp_vdif)
        call nc_write_var(ncid, Td_tmp, N, text, rp_var, output_rec, ierr)

        ! Do boundary changes - already calculated due to need for Hz and Hz_old to differ (they don't when wrt_online_diags is called)
        text=trim(t_vname(itrc)) / / tdname(it_bc)
        call nc_write_var(ncid, Tdiag(START_2D_ARRAY,1,itrc,it_bc), N, text, rp_var, output_rec, ierr)
        T_dsum=T_dsum+Tdiag(:,:,:,itrc,it_bc)  ! re-do this just for boundaries not full domain!

        text=trim(t_vname(itrc)) / / '_dsum'
        call nc_write_var(ncid, T_dsum, N, text, rp_var, output_rec, ierr)

        Td_tmp=0 ! probably don't need this since loop below is full range
        ! write gross difference:
        do k=1,N                    ! first diagnostic calculated using previous step's value
          do j=jstrR,jendR          ! Full domain loop range incl. bry.
            do i=istrR,iendR
              Td_tmp(i,j,k)=t(i,j,k,nnew,itrc)*Hz(i,j,k) - t(i,j,k,nstp,itrc)*Hz_old(i,j,k)
            enddo
          enddo
        enddo
        text=trim(t_vname(itrc)) / / '_diff'
        call nc_write_var(ncid, Td_tmp, N, text, rp_var, output_rec, ierr)

        Td_tmp=Td_tmp-T_dsum ! difference between gross change and sum of diags.
        text=trim(t_vname(itrc)) / / '_d-d'
        call nc_write_var(ncid, Td_tmp, N, text, rp_var, output_rec, ierr)

      enddo ! <- itrc=1,NT

      Hz_old = Hz  ! update Hz_old

!      if(WESTERN_EDGE) then  ! just here for jstrR, shouldn't be in tracers
!        call wrt_1D_or_2D_array(ncid, output_rec, 'u_wbry', jstrR, jendR, 1, N, uv_wbry(jstrR:jendR,1:N,1) )
!      endif

      end subroutine wrt_tracer_diags_tile  !]
# endif /* DIAGS_TRC */

! ----------------------------------------------------------------------
      subroutine wrt_diagnostics( tile )  ![
      implicit none

      ! input
      integer, intent(in) :: tile

# include "compute_tile_bounds.h"

      call wrt_diagnostics_tile(istr, iend, jstr, jend)

      end subroutine wrt_diagnostics  !]

! ----------------------------------------------------------------------
      subroutine wrt_diagnostics_tile(istr, iend, jstr, jend)  ![ write online diagnostics
      implicit none

      ! inputs
      integer, intent(in)     :: istr, iend, jstr, jend

      ! local
      integer                 :: prev_fill_mode,
     &                           r2dgrd(3), u2dgrd(3), v2dgrd(3), auxil(2)   ! 2D grid types
#ifdef SOLVE3D
     &                         , r3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)  ! 3D grid types
#endif
      integer                 :: ierr = 0, ncid = -1, dim, diag, i, j, k, tile, idt, itrc

      character(len=64)       :: fname
      character(len=30)       :: text
      integer, dimension(3)   :: rwebrygrd(3)                  ! experimental for if using 2D bry outputs
      integer, dimension(4,2) :: grdtyp3d                      ! u and v type to using in loop
      integer, dimension(2)   :: gt_wrt = (/ up_var, vp_var /) ! grid type for writing

      if (FIRST_TIME_STEP .and. diag_avg) then                 ! only if diag averaging and only done once
        allocate(    Md_avg(GLOBAL_2D_ARRAY, N, 2, nmd) )      ! allocate memory for diag averages
        allocate(  Md_a_sum(GLOBAL_2D_ARRAY, N, 2)      )      ! sum of diag averages at timestep
        allocate( Md_a_ever(GLOBAL_2D_ARRAY, N, 2)      )      ! continuing total of sum of averages Md_a_ever += Md_a_sum
        allocate(   uvdz_t0(GLOBAL_2D_ARRAY, N, 2)      )
        allocate( uvdz_ever(GLOBAL_2D_ARRAY, N, 2)      )
        Md_avg = 0; Md_a_ever=0                                ! Md_a_sum set to 0 elsewhere, so not needed here
        do diag=1,nmd                                          ! loop through diag terms
          lname(diag)  = 'Average ' / / trim(lname(diag))      ! and prepend 'Average'
        enddo
      endif

      tn = iic-ntstart      ! model step, as timestep t=n. (icc = t(n+1), and ntstart=1 if not from restart). We write the previous step results during the start of the next step

      ! might be able to use diag_avg .or. calc_diag and compiler sees or is parameter!
      if (diag_avg .or. mod(tn+1,rec_rate) == 0) then ! +1 as calculate in the step before writing. so if tn=5 then calc_diag at tn=4, hence need +1 to trigger mod(tn...)
        calc_diag = .true.                            ! This triggers diags to be calculated in equation system. Always on for averages
      else                                            ! set here once rather than everywhere in equation system redo 'mod'
        calc_diag = .false.                           ! since this is called at the beginning of main.
      endif
      if (diag_avg .or. mod(tn+2,rec_rate) == 0) then ! +2 as need to set the previous values before calc_diag step above
        set_dz_old=.true.                             ! So if tn=5 then calc_diag at tn=4, and need set_dz_old tn=3 hence +2 to trigger mod(tn...)
      else
        set_dz_old=.false.
      endif

      if (iic>ntstart) then ! catch to prevent output of zeros at inital step since mod(0,anything) = 0...

        if (diag_avg) then                               ! avoid very first step because diags are only all set by the end of the
            call set_diags_avg_tile(istr,iend,jstr,jend) ! first step, and then added to averages and/or output at the beginning of the next step.
        endif

        if ( write_file==.true.     .and.
     &       mod(tn,rec_rate) == 0        ) then

          if(total_output_rec == 0) then                            ! Add root output name to file name (only done once)
            file_name = trim(output_root_name) / / trim(file_name)  ! Notice below space between '/ /', this is needed to avoid
          endif                                                     ! cpp preprocessor seeing // as a c++ comment and deleting it.


          if(output_rec == 0) then                                  ! if file doesn't yet exist or need new file!

            call ncdf_create_file(file_name, ncid, prev_fill_mode   ! lose the prev_fill_mode? r2dgrd, ...., auxil are id's to dimensions, needed to create vars in the file
     &                    ,recs_per_file, total_output_rec, auxil
#ifdef SOLVE3D
!     &                    ,r3dgrd, u3dgrd, v3dgrd, w3dgrd           ! 3D grid dimensions
     &                    ,r3dgrd, grdtyp3d(:,1), grdtyp3d(:,2), w3dgrd           ! 3D grid dimensions
#endif
     &                    ,r2dgrd, u2dgrd, v2dgrd )                 ! 2D grid dimensions


            ! Define remaining variables:

            do dim=1,2       ! u,v dims
              do diag=1,nmd  ! diag terms
                text  = uvw(dim) / / vnames(diag)     ! diag short name
                call nc_define_var(ncid, text, lname(diag), units, grdtyp3d(:,dim), ierr )
              enddo

              ! following to compare gross change in u/v and sum of diags. Possibly remove.
              text  = uvw(dim) / / 'diff'  ! gross change in dz*u or dz*v over 1 timestep
              call nc_define_var(ncid, text, text, text, grdtyp3d(:,dim), ierr )
              text  = uvw(dim) / / 'dsum'  ! sum of all diag terms over 1 timestep
              call nc_define_var(ncid, text, text, text, grdtyp3d(:,dim), ierr )
              text  = uvw(dim) / / 'd-d'   ! compare gross change against sum of diags. Should be approx zero.
              call nc_define_var(ncid, text, text, text, grdtyp3d(:,dim), ierr )

              if (diag_avg) then
                text  = uvw(dim) / / 'dz_ever'  ! gross change in u/v up until now since t0
                call nc_define_var(ncid, text, text, text, grdtyp3d(:,dim), ierr )
                text  = uvw(dim) / / 'avg_ever' ! total diag avg in u/v up until now since t0
                call nc_define_var(ncid, text, text, text, grdtyp3d(:,dim), ierr )
              endif

            enddo ! <-- dim

# ifdef DIAGS_TRC
            do itrc=1,NT  ! CHANGE THIS TO ONLY tracers of interest!
              do idt=1,ndt ! tracer diag terms
                text=trim(t_vname(itrc)) / / tdname(idt)
                call nc_define_var(ncid, text, tlname(idt), tunits, r3dgrd, ierr )
              enddo
              text  = trim(t_vname(itrc)) / / '_diff'  ! gross change in dz*t over 1 timestep
              call nc_define_var(ncid, text, text, text, r3dgrd, ierr )
              text  = trim(t_vname(itrc)) / / '_dsum'  ! sum of all diag terms over 1 timestep
              call nc_define_var(ncid, text, text, text, r3dgrd, ierr )
              text  = trim(t_vname(itrc)) / / '_d-d'   ! compare gross change against sum of diags. Should be approx zero.
              call nc_define_var(ncid, text, text, text, r3dgrd, ierr )
            enddo  ! <-- itrc
# endif

!            ! define bry 2D vars (experimental):
!            rwebrygrd(1)=r3dgrd(2) ! set eta_rho dimension
!            rwebrygrd(2)=r3dgrd(3) ! set s_rho dimension
!            rwebrygrd(3)=r3dgrd(4) ! set time    dimension
!!          if(WESTERN_EDGE) then
!            call nc_define_var(ncid, 'u_wbry', 'u_wbry', 'u_wbry', rwebrygrd, ierr )
!!          endif

            ierr=nf90_enddef(ncid) ! end definition to write changes to disk
            if (mynode == 0) then
              write(*,'(6x,2A)') 'diagnostics :: created new netCDF file '
     &                            ,file_name
            endif

          endif ! (output_rec == 0; needed to create a file)


          ! WRITE VARIABLES TO FILE
          ! =======================

          ierr=nf90_open(file_name, nf90_write, ncid)

          ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)   ! set fill value - nf90_nofill for optimized writing

                                                                  ! advance output record number
          output_rec = output_rec + 1                             ! current file record
          total_output_rec = total_output_rec + 1                 ! total number of recs written


          call nc_write_time(ncid, output_rec, total_output_rec)  ! write ocean time

          do dim=1,2                            ! over dims
            do diag=1,nmd                       ! over diag terms
              text = uvw(dim) / / vnames(diag)
              if (diag_avg) then
                call nc_write_var(ncid, Md_avg(:,:,:,dim,diag), N, text, gt_wrt(dim), output_rec, ierr)
                text  = uvw(dim) / / 'dz_ever'  ! gross change in u/v up until now since t0
                call nc_write_var(ncid, uvdz_ever(:,:,:,dim),   N, text, gt_wrt(dim), output_rec, ierr)
                text  = uvw(dim) / / 'avg_ever' ! total diag avg in u/v up until now since t0
                call nc_write_var(ncid, Md_a_ever(:,:,:,dim),   N, text, gt_wrt(dim), output_rec, ierr)
              else
                call nc_write_var(ncid,  Mdiag(:,:,:,dim,diag), N, text, gt_wrt(dim), output_rec, ierr)
              endif
            enddo

            ! do uv_dif calculation at end of step3d_uv2.F as have istr, iend there already.
            ! but use 'nnew' there as still in previou timestep!
            text = uvw(dim) / / 'diff'  ! gross change in dz*u or dz*v over 1 timestep
            call nc_write_var(ncid, uv_dif(:,:,:,dim), N, text, gt_wrt(dim), output_rec, ierr)

            Mdsum = 0  ! refresh for this step
            do diag=1,nmd
              Mdsum = Mdsum + Mdiag(:,:,:,dim,diag)  ! should make the ranges more efficient (e.g. bry term don't need interior points)
            enddo
            text = uvw(dim) / / 'dsum'               ! sum of diag terms over 1 timestep
            call nc_write_var(ncid, Mdsum, N, text, gt_wrt(dim), output_rec, ierr)

            diff = uv_dif(:,:,:,dim)-Mdsum           ! difference between gross dz*u (or v) change and sum of diags of dz*u
            text = uvw(dim) / / 'd-d'
            call nc_write_var(ncid, diff, N, text, gt_wrt(dim), output_rec, ierr)

          enddo  ! dim

# ifdef DIAGS_TRC
          call wrt_tracer_diags_tile(istr, iend, jstr, jend, ncid)
# endif

          ierr=nf90_close (ncid)  ! close netcdf file

          if (mynode == 0) then
            write(*,'(6x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')      ! confirm work completed
     &       'diagnostics :: wrote output, tdays =', tdays,
     &       'step =', iic-1, 'rec =', output_rec, '/', total_output_rec
     &        MYID
          endif

          if (output_rec == recs_per_file) then  ! reset output_rec to create new file
            output_rec = 0
          endif

        endif ! <-- write_file==.true.

      endif   ! <-- icc>ntstart

      end subroutine wrt_diagnostics_tile  !]

! ----------------------------------------------------------------------

#endif /* DIAGNOSTICS */

      end module diagnostics
