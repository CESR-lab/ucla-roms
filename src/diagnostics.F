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

      integer,public,parameter :: ndm = 5  !! number of momentum diag terms
      integer,public,parameter :: ndt = 5  !! number of tracer diag terms

      real, public, dimension(GLOBAL_2D_ARRAY, N, 3, ndm)  :: Mdiag=0  ! 3 for u,v,w
!      real,public, dimension(nx,ny,nz,ndt,nt) :: Tdiag
      real, public, dimension(GLOBAL_2D_ARRAY, N, ndm)     :: u_b4=0   ! temporary u for debugging.
      real, public, dimension(GLOBAL_2D_ARRAY, N, ndm)     :: u_af=0   ! u after

      real, public, dimension(GLOBAL_2D_ARRAY, N)          :: u_old=0! u from previous timestep
      real, public, dimension(GLOBAL_2D_ARRAY, N)          :: u_dif=0! u difference current vs previous u
      real, public, dimension(GLOBAL_2D_ARRAY, N)          :: Mdsum=0! sum of Mdiag

      real, public, dimension(GLOBAL_2D_ARRAY, N)          :: diff=0

      real, public, dimension(GLOBAL_2D_ARRAY, N, 3)       :: rus=0  ! 3 contributions to ru
      real, public, dimension(GLOBAL_2D_ARRAY, N, 3)       :: drus=0 ! changes in 3 ru terms

      real, public, dimension(GLOBAL_2D_ARRAY, N)          :: Hz_old=0

      integer, parameter, public :: iuv1   = 1,  ! indices in Mdiag array
     &                              ivisc  = 2,
     &                              iuv2DC = 3,
     &                              iuv2FC = 4,
     &                              iu3d   = 5,
     &                              iprsgr = 6,
     &                              icori  = 7,
     &                              iadv   = 8


      character (len=20), dimension(ndm) :: vname = (/
!     &                                                 'hydrostatic_pg'
     &                                                 'uv1_ru',
     &                                                 'visc3d',
     &                                                 'uv2_DC',
     &                                                 'uv2_FC',
     &                                                 'u3dbc'
     &                                                                  /)

      character (len=60), dimension(ndm) :: lname = (/
     &  'uv1_ru',
     &  'visc3d L376',
     &  'uv2_DC L176',
     &  'uv2_FC L463',
     &  'u3dbc'
     &                                                                  /)

      character (len=20), dimension(ndm) :: units = (/
     &                                                 'm^2/s (dz*u)',
     &                                                 'm^2/s (dz*u)',
     &                                                 'm^2/s (dz*u)',
     &                                                 'm^2/s (dz*u)',
     &                                                 'm^2/s (dz*u)'
     &                                                                  /)

      character (len=20), dimension(3) :: ruvname = (/
     &                                                 'prsgrd',
     &                                                 'coriolis',
     &                                                 'advection'
     &                                                                  /)

      character (len=60), dimension(3) :: rulname = (/
     &  'prsgrd L275',
     &  'horiz_uv.h L35 coriolis & curvilinear coords',
     &  'advection'
     &                                                                  /)

      character (len=20), dimension(3) :: ruunits = (/
     &                                                 'm^2/s (dz*u)',
     &                                                 'm^2/s (dz*u)',
     &                                                 'm^2/s (dz*u)'
     &                                                                  /)

      character (len=2), dimension(3) :: uvw = (/ 'u_', 'v_', 'w_' /)

      ! Netcdf outputting:
      integer            :: output_rec = 0              ! record number of output. 0 indicates we need new file
      integer            :: total_output_rec = 0        ! total flux output recordings so far
      character (len=60) :: file_name = '_on_diags.nc'  ! suffix to append to simulation output rootname

      public wrt_online_diags

      contains

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

      if (write_file==.true. .and. mod(iic-1,rec_rate) == 0) then ! 1 = ntstart

        if(total_output_rec == 0) then                            ! Add root output name to file name (only done once)
          file_name = trim(output_root_name) / / trim(file_name)  ! Notice below space between '/ /', this is needed to avoid
!#ifndef HZU
!          u_old = u(:,:,:,nstp)  !debug only possibly
!#else
!!          u_old = u(:,:,:,nstp) * Hz ! Now done in pre_step3d since need loop ranges for dz*u calc
!#endif
        endif                                                     ! cpp preprocessor seeing // as a c++ comment and deleting it.


        if(output_rec == 0) then                                  ! if file doesn't yet exist or need new file!

          call ncdf_create_file(file_name, ncid, prev_fill_mode   ! lose the prev_fill_mode? r2dgrd, ...., auxil are id's to dimensions, needed to create vars in the file
     &                  ,recs_per_file, total_output_rec, auxil
#ifdef SOLVE3D
     &                  ,r3dgrd, u3dgrd, v3dgrd, w3dgrd           ! 3D grid dimensions
#endif
     &                  ,r2dgrd, u2dgrd, v2dgrd )                 ! 2D grid dimensions


          ! Define remaining variables:
!          do diag=1,ndm
!            do dim=1,1  ! loop through dims
!              text  = uvw(dim) / / vname(diag)  ! diag short name
!              call nc_define_var(ncid, text, lname(diag), units(diag), u3dgrd, ierr )
!            enddo
!          enddo

!          do diag=1,ndm
!            text  = 'ub4_' / / vname(diag)  ! diag short name
!            call nc_define_var(ncid, text, lname(diag), units(diag), u3dgrd, ierr )
!          enddo

!          do diag=1,ndm
!            text  = 'uaf_' / / vname(diag)  ! diag short name
!            call nc_define_var(ncid, text, lname(diag), units(diag), u3dgrd, ierr )
!          enddo

          do diag=1,ndm
            text  = 'Md_' / / vname(diag)  ! diag short name
            call nc_define_var(ncid, text, lname(diag), units(diag), u3dgrd, ierr )
          enddo

          do diag=1,3
            do dim=1,1  ! loop through dims
              text  = uvw(dim) / / ruvname(diag)  ! diag short name
              call nc_define_var(ncid, text, rulname(diag), ruunits(diag), u3dgrd, ierr )
            enddo
          enddo

          call nc_define_var(ncid, 'u_diff', 'u_diff', 'u_diff', u3dgrd, ierr )
          call nc_define_var(ncid, 'Mdsum',  'Mdsum',  'Mdsum',  u3dgrd, ierr )
          call nc_define_var(ncid, 'rusum',  'rusum',  'rusum',  u3dgrd, ierr )
          call nc_define_var(ncid, 'diff',    'diff',   'diff',  u3dgrd, ierr )


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

!        do diag=1,ndm
!          do dim=1,1  ! loop through dims
!            text = uvw(dim) / / vname(diag)
!            call nc_write_var(ncid, Mdiag(:,:,:,dim,diag), N, text, up_var, output_rec, ierr)
!          enddo
!        enddo

!        do diag=1,ndm
!          text = 'ub4_' / / vname(diag)
!          call nc_write_var(ncid, u_b4(:,:,:,diag), N, text, up_var, output_rec, ierr)
!        enddo
!
!        do diag=1,ndm
!          text = 'uaf_' / / vname(diag)
!          call nc_write_var(ncid, u_af(:,:,:,diag), N, text, up_var, output_rec, ierr)
!        enddo

        do diag=1,ndm
          text = 'Md_' / / vname(diag)
          call nc_write_var(ncid, Mdiag(:,:,:,1,diag), N, text, up_var, output_rec, ierr)
        enddo

        drus(:,:,:,1)=rus(:,:,:,1)
        drus(:,:,:,2)=rus(:,:,:,2)-rus(:,:,:,1)
        drus(:,:,:,3)=rus(:,:,:,3)-rus(:,:,:,2)
        do diag=1,3
          do dim=1,1  ! loop through dims
            text = uvw(dim) / / ruvname(diag)
            call nc_write_var(ncid, drus(:,:,:,diag), N, text, up_var, output_rec, ierr)
          enddo
        enddo
        call nc_write_var(ncid, rus(:,:,:,3), N, 'rusum', up_var, output_rec, ierr)

!#ifndef HZU
!        u_dif = u(:,:,:,nstp) - u_old
!        call nc_write_var(ncid, u_dif, N, 'u_diff', up_var, output_rec, ierr)
!        u_old = u(:,:,:,nstp)
!#else
        ! do u_dif calculation at end of step3d_uv2.F as have istr, iend there already.
        ! but use 'nnew' there as still in previou timestep!
!        u_dif = ( u(:,:,:,nstp) * Hz ) - u_old ! Now in units of dz*u!
        call nc_write_var(ncid, u_dif, N, 'u_diff', up_var, output_rec, ierr)
!        u_old = u(:,:,:,nstp) * Hz
!#endif

        Mdsum = 0  ! refresh for this step
        do diag=2,ndm ! 2 to ignore uv1_ru
          Mdsum = Mdsum + Mdiag(:,:,:,1,diag)
        enddo
        Mdsum = Mdsum + rus(:,:,:,3)
        call nc_write_var(ncid, Mdsum, N, 'Mdsum', up_var, output_rec, ierr)

        diff=u_dif-Mdsum
        call nc_write_var(ncid, diff, N, 'diff', up_var, output_rec, ierr)

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
