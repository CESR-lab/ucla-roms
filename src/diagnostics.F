      module diagnostics

      ! Online diagnostics for ROMS simulation.
      ! Initial coding DevinD and JeroenM - 2021/05

#include "cppdefs.h"
#ifdef DIAGNOSTICS

      use dimensions
      use read_write
      use netcdf

      implicit none
      private

#include "ocean2d.h"
#include "ocean3d.h"
! scalars.h for 'knew' var
#include "scalars.h"

      ! ***************************************************************
      ! User inputs
      logical :: write_file    = .true.  ! whether to write a file or not (probably remove this)
      integer :: rec_rate      = 1       ! output record rate,
      integer :: recs_per_file = 20       ! total recs per file

      ! End of user inputs
      ! ***************************************************************

      ! Preamble:  ![

      ! Mdiag(2): hydrostatic pressure gradient
      ! Mdiag(3): cor(u)
      ! Mdiag(4): adv(u) (non-dissipative part)
      ! Mdiag(5): horizontal mixing (including numerical dissipation)
      ! Mdiag(6): vertical mixing   (including numerical dissipation)
      ! Mdiag(7): non-hydrostatic pressure gradient

      ! Tdiag(1): dt/dt
      ! Tdiag(2): adv(t) (non-dissipative part)
      ! Tdiag(3): adv(t) (dissipative part)
      ! Tdiag(4): vertical mixing
      ! Tdiag(5): sources/sinks

      integer,public,parameter :: nmd = 7  !! number of momentum diag terms
      integer,public,parameter :: ndt = 5  !! number of tracer diag terms

      ! NOTE: u & v are 1:N in z, but w is 0:N hence different variables needed

      real, public, dimension(GLOBAL_2D_ARRAY, N, 2, nmd) :: Mdiag=0    ! 3 for u,v,w
!      real,public, dimension(nx,ny,nz,ndt,nt) :: Tdiag

      real, public, dimension(GLOBAL_2D_ARRAY, N, 2)      :: uvdz_old=0 ! u from previous timestep, need to either store this or Hz to get old u*dz.
      real, public, dimension(GLOBAL_2D_ARRAY, N, 2)      :: uv_dif=0   ! gross u/v difference between current vs previous timestep u/v
      real, public, dimension(GLOBAL_2D_ARRAY, N)         :: Mdsum=0    ! sum of Mdiag
      real, public, dimension(GLOBAL_2D_ARRAY, N)         :: diff=0     ! = uv_dif - Mdsum
      real, public, dimension(GLOBAL_2D_ARRAY, N, 2)      :: uv_prev=0  ! store incremental u/v/w (nnew)

      real, public :: Mdtmp ! temporary scalar for calculations of Mdiag

      integer, parameter, public :: iprsgr = 1, ! internal loop ranges only (i.e. excluding domain boundary)
     &                              icori  = 2, ! internal loop ranges only
     &                              iadv   = 3, ! internal loop ranges only
     &                              ivisc  = 4, ! internal loop ranges only
     &                              iuv2DC = 5, ! internal loop ranges only
     &                              iu3d   = 6, ! bry only
     &                              iuv2FC = 7  ! internal + bry (full domain) loop ranges


      character (len=20), dimension(nmd) :: vname = (/
     &  'prsgrd',
     &  'coriolis',
     &  'advection',
     &  'visc3d',
     &  'uv2_DC',
     &  '3dbc',       ! should store 3dbc changes in 2D arrays not 3D, as just 0 in interior
     &  'uv2_FC'
     &                                                                  /)

      character (len=60), dimension(nmd) :: lname = (/
     &  'Hydrostatic pressure gradient - prsgrd.F',
     &  'Coriolis & curvilinear coords - compute_horiz_rhs_uv_term.h',
     &  'Advection - step3d_uv1.F',
     &  'Viscous terms - visc3d_GP.F',
     &  'Vertical mixing & 2D/3D coupling - step3d_uv2.F part 1',
     &  'Boundary changes - u/v3dbc_im.F',
     &  '2D-3D correction - step3d_uv2.F part 2'
     &                                                                  /)

      character (len=20), dimension(nmd) :: units = (/
     &  'm^2/s (dz*u)',
     &  'm^2/s (dz*u)',
     &  'm^2/s (dz*u)',
     &  'm^2/s (dz*u)',
     &  'm^2/s (dz*u)',
     &  'm^2/s (dz*u)',
     &  'm^2/s (dz*u)'
     &                                                                  /)

      character (len=2), dimension(3) :: uvw = (/ 'u_', 'v_', 'w_' /)

      ! Netcdf outputting:
      integer            :: output_rec = 0              ! record number of output. 0 indicates we need new file
      integer            :: total_output_rec = 0        ! total flux output recordings so far
      character (len=60) :: file_name = '_on_diags.nc'  ! suffix to append to simulation output rootname

# ifdef NHMG
      integer,public,parameter :: nwd = 5  !! number of 'w' momentum diag terms (for NHMG runs)

      real, public, dimension(GLOBAL_2D_ARRAY, 0:N, nwd) :: Wdiag=0    ! 3 for u,v,w
      real, public, dimension(GLOBAL_2D_ARRAY, 0:N)      :: wdz_old=0 ! u from previous timestep, need to either store this or Hz to get old u*dz.
      real, public, dimension(GLOBAL_2D_ARRAY, 0:N)      :: w_dif=0   ! gross u/v difference between current vs previous timestep u/v
      real, public, dimension(GLOBAL_2D_ARRAY, 0:N)      :: Wdsum=0    ! sum of Mdiag
      real, public, dimension(GLOBAL_2D_ARRAY, 0:N)      :: wdif=0     ! = w_dif - Mdsum
      real, public, dimension(GLOBAL_2D_ARRAY, 0:N)      :: w_prev=0  ! store incremental u/v/w (nnew)

      integer, parameter, public :: iwprsgr   = 1,
     &                              iwhoriadv = 2,
     &                              iwvertadv = 3,
     &                              iwuv2     = 4,
     &                              iwbc      = 5

      character (len=20), dimension(nwd) :: wvname = (/
     &  'w_prsgrd',
     &  'w_horiz_flx',
     &  'w_vert_adv',
     &  'w_2D_3D_coupling',
     &  'w_3dbc'
     &                                                                  /)

      character (len=60), dimension(nwd) :: wlname = (/
     &  'prsgrd.F',
     &  'compute_horiz_rhs_w_terms.h',
     &  'compute_vert_rhs_w_terms.h',
     &  'step3d_uv2.F',
     &  'w3dbc_im.F'
     &                                                                  /)

      character (len=20), dimension(nwd) :: wunits = (/
     &  'm^2/s (dz*u)',
     &  'm^2/s (dz*u)',
     &  'm^2/s (dz*u)',
     &  'm^2/s (dz*u)',
     &  'm^2/s (dz*u)'
     &                                                                  /)

      public set_diags_w_at_uv1
      public set_diags_w_at_uv2
      public set_diags_w_at_uv2_end
      public set_diags_w_at_bc
# endif /* NHMG */

      public set_diags_at_visc3d
      public set_diags_at_uv1
      public set_diags_at_uv2_DC
      public set_diags_at_bc_u
      public set_diags_at_bc_v
      public set_diags_at_uv2_FCu
      public set_diags_at_uv2_FCv
      public wrt_online_diags !]

      contains

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

      do k=1,N  ! COULD REUSE THIS BY LOOPING OVER a AND b AND SWITCHING ARGUEMENTS
        do j=jstrR,jendR  ! since this section uv2FC uses IU_RANGE & J_RANGE which includes boundary
          do i=istr,iendR
            Mdiag(i,j,k,1,iuv2FC) = u(i,j,k,nnew) * (0.5*( Hz(i-1,j,k)+Hz(i,j,k) )) - uv_prev(i,j,k,1)

            ! Full loop ranges since u-change over every point including bry.
            ! istr is the west boundary for u-points, istrR=E, jstrR=S, jendR=N
            ! can't just do (u(nnew)-u(nstp)*Hz here because actual need to old Hz for u(nstp) as Hz has evolved.
            uv_dif(i,j,k, 1)   = u(i,j,k,nnew) * (0.5*( Hz(i-1,j,k)+Hz(i,j,k) )) - uvdz_old(i,j,k, 1) ! do this here as units are correct
            uvdz_old(i,j,k, 1) = u(i,j,k,nnew) * (0.5*( Hz(i-1,j,k)+Hz(i,j,k) ))
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
            uvdz_old(i,j,k, 2) = v(i,j,k,nnew) * (0.5*( Hz(i,j-1,k)+Hz(i,j,k) ))
          enddo
        enddo
      enddo

!      if(mynode==0) then  ! debug
!        print *, 'v(10,10,48,nnew)=',v(10,10,48,nnew)
!        print *, 'uvdz_old(10,10,48,2)=', uvdz_old(10,10,48,2)
!        print *, 'uv_dif  (10,10,48,2)=', uv_dif(10,10,48,2)
!        print *, 'Hz     (10,10,48)  =', Hz(10,10,48)
!      endif

      end subroutine set_diags_at_uv2_FCv  !]


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
      subroutine wrt_online_diags  ![ write online diagnostics

      implicit none

      ! local
      integer :: prev_fill_mode,
     &           r2dgrd(3), u2dgrd(3), v2dgrd(3), auxil(2)
#ifdef SOLVE3D
     &         , r3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
#endif
      integer :: ierr = 0, ncid = -1, dim, diag, i, j, k
      character(len=64) :: fname
      character(len=20) :: text

      integer, dimension(4,2) :: grdtyp3d ! u and v type to using in loop
      integer, dimension(2) :: gt_wrt = (/ up_var, vp_var /) ! grid type for writing

      if (write_file==.true. .and. mod(iic-1,rec_rate) == 0) then ! 1 = ntstart

        if(total_output_rec == 0) then                            ! Add root output name to file name (only done once)
          file_name = trim(output_root_name) / / trim(file_name)  ! Notice below space between '/ /', this is needed to avoid
        endif                                                     ! cpp preprocessor seeing // as a c++ comment and deleting it.


        if(output_rec == 0) then                                  ! if file doesn't yet exist or need new file!

          call ncdf_create_file(file_name, ncid, prev_fill_mode   ! lose the prev_fill_mode? r2dgrd, ...., auxil are id's to dimensions, needed to create vars in the file
     &                  ,recs_per_file, total_output_rec, auxil
#ifdef SOLVE3D
!     &                  ,r3dgrd, u3dgrd, v3dgrd, w3dgrd           ! 3D grid dimensions
     &                  ,r3dgrd, grdtyp3d(:,1), grdtyp3d(:,2), w3dgrd           ! 3D grid dimensions
#endif
     &                  ,r2dgrd, u2dgrd, v2dgrd )                 ! 2D grid dimensions


          ! Define remaining variables:
          do dim=1,2  ! loop through dims
            do diag=1,nmd

              text  = uvw(dim) / / vname(diag)  ! diag short name
              call nc_define_var(ncid, text, lname(diag), units(diag), grdtyp3d(:,dim), ierr )
            enddo

            text  = uvw(dim) / / 'diff'
            call nc_define_var(ncid, text, text, text, grdtyp3d(:,dim), ierr )
            text  = uvw(dim) / / 'dsum'
            call nc_define_var(ncid, text, text, text, grdtyp3d(:,dim), ierr )
            text  = uvw(dim) / / 'd-d'
            call nc_define_var(ncid, text, text, text, grdtyp3d(:,dim), ierr )
          enddo

          ierr=nf90_enddef(ncid) ! end definition to write changes to disk
          if (mynode == 0) then
            write(*,'(6x,2A)') 'diagnostics :: created new netCDF file '
     &                          ,file_name
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

        do dim=1,2  ! loop through dims
          do diag=1,nmd
            text = uvw(dim) / / vname(diag)
            call nc_write_var(ncid, Mdiag(:,:,:,dim,diag), N, text, gt_wrt(dim), output_rec, ierr)
          enddo

          ! do uv_dif calculation at end of step3d_uv2.F as have istr, iend there already.
          ! but use 'nnew' there as still in previou timestep!
          text = uvw(dim) / / 'diff'
          call nc_write_var(ncid, uv_dif(:,:,:,dim), N, text, gt_wrt(dim), output_rec, ierr)

          Mdsum = 0  ! refresh for this step
          do diag=1,nmd
            Mdsum = Mdsum + Mdiag(:,:,:,dim,diag)
          enddo
          text = uvw(dim) / / 'dsum'
          call nc_write_var(ncid, Mdsum, N, text, gt_wrt(dim), output_rec, ierr)

          diff=uv_dif(:,:,:,dim)-Mdsum ! difference between gross dz*u change and sum of diags of dz*u
          text  = uvw(dim) / / 'd-d'
          call nc_write_var(ncid, diff, N, text, gt_wrt(dim), output_rec, ierr)

        enddo  ! dim

        ierr=nf90_close (ncid)  ! close netcdf file

        if (mynode == 0) then
          write(*,'(6x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')      ! confirm work completed
     &     'diagnostics :: wrote output, tdays =', tdays,
     &     'step =', iic-1, 'rec =', output_rec, '/', total_output_rec
     &      MYID
        endif

        if (output_rec == recs_per_file) then  ! reset output_rec to create new file
          output_rec = 0
        endif

      endif ! (write_file==.true.)

      end subroutine wrt_online_diags  !]

! ----------------------------------------------------------------------

#endif /* DIAGNOSTICS */

      end module diagnostics
