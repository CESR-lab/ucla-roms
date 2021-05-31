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

      real, public, dimension(GLOBAL_2D_ARRAY, N, 3, nmd)  :: Mdiag=0   ! 3 for u,v,w
!      real,public, dimension(nx,ny,nz,ndt,nt) :: Tdiag

      real, public, dimension(GLOBAL_2D_ARRAY, N, 3)       :: udz_old=0 ! u from previous timestep, need to either store this or Hz to get old u*dz.
      real, public, dimension(GLOBAL_2D_ARRAY, N, 3)       :: u_dif=0   ! u difference current vs previous u
      real, public, dimension(GLOBAL_2D_ARRAY, N)          :: Mdsum=0   ! sum of Mdiag

      real, public, dimension(GLOBAL_2D_ARRAY, N)          :: diff=0    ! = u_dif - Mdsum

      real, public, dimension(GLOBAL_2D_ARRAY, N, 3)       :: u_prev=0  ! store incremental u/v/w (nnew)

      real, public :: Mdtmp ! temporary scalar for calculations of Mdiag

      integer, parameter, public :: iprsgr = 1, ! internal loop ranges only (i.e. excluding domain boundary)
     &                              icori  = 2, ! internal loop ranges only
     &                              iadv   = 3, ! internal loop ranges only
     &                              ivisc  = 4, ! internal loop ranges only
     &                              iuv2DC = 5, ! internal loop ranges only
     &                              iu3d   = 6, ! bry only
     &                              iuv2FC = 7  ! internal + bry (full domain) loop ranges


      character (len=20), dimension(nmd) :: vname = (/
     &                                                 'prsgrd',
     &                                                 'coriolis',
     &                                                 'advection',
     &                                                 'visc3d',
     &                                                 'uv2_DC',
     &                                                 'u3dbc',
     &                                                 'uv2_FC'
     &                                                                  /)

      character (len=60), dimension(nmd) :: lname = (/
     &  'prsgrd L275',
     &  'horiz_uv.h L35 coriolis & curvilinear coords',
     &  'advection',
     &  'visc3d L376',
     &  'uv2_DC L176',
     &  'u3dbc',
     &  'uv2_FC L463'
     &                                                                  /)

      character (len=20), dimension(nmd) :: units = (/
     &                                                 'm^2/s (dz*u)',
     &                                                 'm^2/s (dz*u)',
     &                                                 'm^2/s (dz*u)',
     &                                                 'm^2/s (dz*u)',
     &                                                 'm^2/s (dz*u)',
     &                                                 'm^2/s (dz*u)',
     &                                                 'm^2/s (dz*u)'
     &                                                                  /)

      character (len=2), dimension(3) :: uvw = (/ 'u_', 'v_', 'w_' /)

      ! Netcdf outputting:
      integer            :: output_rec = 0              ! record number of output. 0 indicates we need new file
      integer            :: total_output_rec = 0        ! total flux output recordings so far
      character (len=60) :: file_name = '_on_diags.nc'  ! suffix to append to simulation output rootname

      public set_diags_at_visc3d
      public set_diags_at_uv1
      public set_diags_at_uv2_DC
      public set_diags_at_bc
      public set_diags_at_uv2_FCu
      public set_diags_at_uv2_FCv
      public wrt_online_diags

      contains

! ----------------------------------------------------------------------
      subroutine set_diags_at_uv1( istr, iend, jstr, jend, istrU, j, dim, DC, ru, sustr )  ![

      implicit none

      ! inputs
      integer, intent(in) :: istr, iend, jstr, jend, istrU, j, dim ! istr, iend, etc. 'A' for any as depends on inputs, which 'dim' dimension.
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N),   intent(in) :: DC       ! (:,:)   prevents local copy
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,N), intent(in) :: ru       ! (:,:,:) prevents local copy - ru or rv
      real, dimension(GLOBAL_2D_ARRAY)  , intent(in) :: sustr

      ! local
      integer :: i, k

      do i=istrU,iend
        Mdiag(i,j,N,1,iprsgr) = Mdiag(i,j,N,1,iprsgr) * DC(i,0)  ! need to convert to dz*u now that we have DC
        Mdiag(i,j,1,1,iprsgr) = Mdiag(i,j,1,1,iprsgr) * DC(i,0)  ! consider taking this out of loop and doing seperately
        Mdiag(i,j,N,1,icori)  = Mdiag(i,j,N,1,icori)  * DC(i,0)  ! perhaps even in own subroutine to avoid clutter
        Mdiag(i,j,1,1,icori)  = Mdiag(i,j,1,1,icori)  * DC(i,0)
        Mdiag(i,j,N,1,iadv)   = DC(i,0)*ru(i,j,N) + dt*sustr(i,j)                 ! all advection - includes 2nd part of compute_horiz_uv & compute_vert_uv.
     &                           - ( Mdiag(i,j,N,1,iprsgr) + Mdiag(i,j,N,1,icori) ) ! subtract previous 2 terms included in total
        u_prev(i,j,N,1)       = u(i,j,N,nnew)     ! store for visc3d
        Mdiag(i,j,1,1,iadv)   = DC(i,0)*ru(i,j,1)
     &                           - ( Mdiag(i,j,1,1,iprsgr) + Mdiag(i,j,1,1,icori) )
        u_prev(i,j,1,1)       = u(i,j,1,nnew)     ! store for visc3d
      enddo

      do k=2,N-1
        do i=istrU,iend
          Mdiag(i,j,k,1,iprsgr) = Mdiag(i,j,k,1,iprsgr) * DC(i,0)
          Mdiag(i,j,k,1,icori)  = Mdiag(i,j,k,1,icori)  * DC(i,0)
          u_prev(i,j,k,1)       = u(i,j,k,nnew)
          Mdiag(i,j,k,1,iadv)   = DC(i,0)*ru(i,j,k) ! all advection - includes 2nd part of compute_horiz_uv & compute_vert_uv.
     &                            - ( Mdiag(i,j,k,1,iprsgr) + Mdiag(i,j,k,1,icori) )
        enddo
      enddo

      end subroutine set_diags_at_uv1  !]

! ----------------------------------------------------------------------
      subroutine set_diags_at_visc3d( isA, iend, jsA, jend, dim )  ![
      implicit none

      ! inputs
      integer, intent(in) :: isA, iend, jsA, jend, dim ! changed the name to use for u or v
      integer :: i, j, k  ! local

      do k=1,N
        do j=jsA,jend    ! loop ranges for interior of subdomains, exludes global domain boundaries
          do i=isA,iend
            Mdiag(i,j,k,dim,ivisc) = u(i,j,k,nnew) - u_prev(i,j,k,dim)  ! units already in dz*u
            u_prev(i,j,k,dim)      = u(i,j,k,nnew)
          enddo
        enddo
      enddo
!      Mdiag (istrU:iend,jstr:jend,:,1,ivisc) = u(istrU:iend,jstr:jend,:,nnew)
!      u_prev(istrU:iend,jstr:jend,:)         = u(istrU:iend,jstr:jend,:,nnew) ! can remove

      end subroutine set_diags_at_visc3d  !]

! ----------------------------------------------------------------------
      subroutine set_diags_at_uv2_DC( isA, iend, jsA, jend, dim )  ![
      implicit none

      ! inputs
      integer, intent(in) :: isA, iend, jsA, jend, dim ! changed the name to use for u or v
      integer :: i, j, k  ! local

      do k=1,N
        do j=jsA,jend    ! loop ranges for interior of subdomains, exludes global domain boundaries
          do i=isA,iend  ! could do vector operation
            Mdtmp                 = u(i,j,k,nnew) * (0.5*( Hz(i-1,j,k)+Hz(i,j,k) )) ! convert to dz*u
            Mdiag(i,j,k,1,iuv2DC) = Mdtmp - u_prev(i,j,k,1)  ! avoids recomputing above
            u_prev(i,j,k,1) = Mdtmp
          enddo
        enddo
      enddo

      end subroutine set_diags_at_uv2_DC  !]

! ----------------------------------------------------------------------
      subroutine set_diags_at_bc( istr, iend, jstr, jend )  ![
      ! set the diagnostics terms for boundary changes

      implicit none

      ! inputs/outputs
      integer, intent(in) :: istr, iend, jstr, jend

      ! local
      integer :: i

      if (WESTERN_EDGE)  u_prev(istr  ,:,:,1) = u(istr  ,:,:,nnew) * 0.5*(Hz(istr-1,:,:)+Hz(istr  ,:,:))
      if (EASTERN_EDGE)  u_prev(iend+1,:,:,1) = u(iend+1,:,:,nnew) * 0.5*(Hz(iend  ,:,:)+Hz(iend+1,:,:))
      if (WESTERN_EDGE)  Mdiag(istr  ,:,:,1,iu3d) = u_prev(istr  ,:,:,1) - udz_old(istr  ,:,:, 1) ! u_prev is already u_current, can do this because have udz_old already
      if (EASTERN_EDGE)  Mdiag(iend+1,:,:,1,iu3d) = u_prev(iend+1,:,:,1) - udz_old(iend+1,:,:, 1)
      do i=istr,iend+1 ! use loop here else indexing would be even messier
        if (SOUTHERN_EDGE) u_prev(i,jstr-1,:,1)=u(i,jstr-1,:,nnew)*0.5*(Hz(i-1,jstr-1,:)+Hz(i,jstr-1,:))
        if (NORTHERN_EDGE) u_prev(i,jend+1,:,1)=u(i,jend+1,:,nnew)*0.5*(Hz(i-1,jend+1,:)+Hz(i,jend+1,:))
      enddo
      if (SOUTHERN_EDGE) Mdiag(:,jstr-1,:,1,iu3d) = u_prev(:,jstr-1,:,1) - udz_old(:,jstr-1,:, 1)
      if (NORTHERN_EDGE) Mdiag(:,jend+1,:,1,iu3d) = u_prev(:,jend+1,:,1) - udz_old(:,jend+1,:, 1)

      end subroutine set_diags_at_bc

! ----------------------------------------------------------------------
      subroutine set_diags_at_uv2_FCu( istr, iendR, jstrR, jendR )  ![
      implicit none

      ! inputs
      integer, intent(in) :: istr, iendR, jstrR, jendR ! changed the name to use for u or v
      integer :: i, j, k  ! local

      do k=1,N  ! COULD REUSE THIS BY LOOPING OVER a AND b AND SWITCHING ARGUEMENTS
        do j=jstrR,jendR  ! since this section uv2FC uses IU_RANGE & J_RANGE which includes boundary
          do i=istr,iendR
            Mdiag(i,j,k,1,iuv2FC) = u(i,j,k,nnew) * (0.5*( Hz(i-1,j,k)+Hz(i,j,k) )) - u_prev(i,j,k,1)

            ! Full loop ranges since u-change over every point including bry.
            ! istr is the west boundary for u-points, istrR=E, jstrR=S, jendR=N
            ! can't just do (u(nnew)-u(nstp)*Hz here because actual need to old Hz for u(nstp) as Hz has evolved.
            u_dif(i,j,k, 1) = u(i,j,k,nnew) * (0.5*( Hz(i-1,j,k)+Hz(i,j,k) )) - udz_old(i,j,k, 1) ! do this here as units are correct
            udz_old(i,j,k, 1)  = u(i,j,k,nnew) * (0.5*( Hz(i-1,j,k)+Hz(i,j,k) ))
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
            Mdiag(i,j,k,2,iuv2FC) = v(i,j,k,nnew) * (0.5*( Hz(i,j-1,k)+Hz(i,j,k) )) - u_prev(i,j,k,2)
            ! put = Mdiag below rather than same equation...
            ! Full loop ranges since u-change over every point including bry.
            ! istr is the west boundary for u-points, istrR=E, jstrR=S, jendR=N
            ! can't just do (u(nnew)-u(nstp)*Hz here because actual need to old Hz for u(nstp) as Hz has evolved.
            u_dif(i,j,k, 2)   = v(i,j,k,nnew) * (0.5*( Hz(i,j-1,k)+Hz(i,j,k) )) - udz_old(i,j,k, 2) ! do this here as units are correct
            udz_old(i,j,k, 2) = v(i,j,k,nnew) * (0.5*( Hz(i,j-1,k)+Hz(i,j,k) ))
          enddo
        enddo
      enddo

      if(mynode==0) then  ! debug
        print *, 'v(10,10,48,nnew)=',v(10,10,48,nnew)
        print *, 'udz_old(10,10,48,2)=', udz_old(10,10,48,2)
        print *, 'u_dif  (10,10,48,2)=', u_dif(10,10,48,2)
        print *, 'Hz     (10,10,48)  =', Hz(10,10,48)
      endif

      end subroutine set_diags_at_uv2_FCv  !]

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

          ! do u_dif calculation at end of step3d_uv2.F as have istr, iend there already.
          ! but use 'nnew' there as still in previou timestep!
          text = uvw(dim) / / 'diff'
          call nc_write_var(ncid, u_dif(:,:,:,dim), N, text, gt_wrt(dim), output_rec, ierr)

          Mdsum = 0  ! refresh for this step
          do diag=1,nmd
            Mdsum = Mdsum + Mdiag(:,:,:,dim,diag)
          enddo
          text = uvw(dim) / / 'dsum'
          call nc_write_var(ncid, Mdsum, N, text, gt_wrt(dim), output_rec, ierr)

          diff=u_dif(:,:,:,dim)-Mdsum ! difference between gross dz*u change and sum of diags of dz*u
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
