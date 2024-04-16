#include "cppdefs.opt"

      subroutine init_arrays

      use param
      use hidden_mpi_vars
      use private_scratch

      implicit none

      integer,save :: tile=1
      integer      :: i,j

#include "compute_tile_bounds.h"

      call init_arrays_private_scratch
#ifdef SOLVE3D
      do i=1,N2d
        iA2d(i,1)=0
        iA2d(i,2)=0
      enddo
#endif

      call init_arrays_tile(istr,iend,jstr,jend)
      end

      subroutine init_arrays_tile (istr,iend,jstr,jend)

! This routine initialize "first-touches" model shared arrays. Most
! of them are assigned to zeros, vertical mixing coefficients are set
! to their background values and will remain unchanged if no vertical
! mixing scheme is applied. The main point here is that because of the
! "first touch" default data placement policy on Linux operating
! system, this operation actually performs distribution of the shared
! arrays accross the nun-uniform-access memory (NUMA) computer (i.e.,
! within a cluster node) unless another distribution policy is
! specified to override the default.

      ! Required modules
      ! ----------------

      use param
      use wec_frc
      use bulk_frc
      use surf_flux ! for: sustr, svstr, stflx, sst_data, sss_data
      use tracers   ! iTandS, t, t_avg
      use coupling
      use eos_vars
      use grid
      use mixing
      use ocean_vars
      use basic_output
      use scalars
      use mpi
      use work_mod
      use boundary
      use buffer
      use mess_buffers
      use dimensions
      use flux_frc
      use pipe_frc
      use bgc_ecosys_vars
      use bgc_forces
      use bgc

      implicit none
      integer, intent(in) :: istr, iend, jstr, jend
      integer             :: i, j, k, itrc, itavg
      real, parameter     :: init0=0.    !!!!  0xFFFA5A5A ==> NaN
#define ALL_DATA
#undef ALL_DATA
#ifdef PRINT_TILE_RANGES
# ifdef MPI
      integer status(MPI_STATUS_SIZE), blank, ierr
# endif
#endif

#include "compute_extended_bounds.h"

#ifdef PRINT_TILE_RANGES
# ifdef MPI
      if (mynode>0) then
        call MPI_Recv (blank, 1, MPI_INTEGER, mynode-1,
     &                 1, ocean_grid_comm, status, ierr)
      endif
      i=mynode
# else
      i=proc(2)
# endif
      write(*,'(I4/2(6x,A6,I3,3x,A6,I3))') i, 'istr =',istr,
     &        'iend =',iend,   'jstr =',jstr, 'jend =',jend
      write(*,'(4x,2(6x,A6,I3,3x,A6,I3)/)')   'istrR=',istrR,
     &        'iendR=',iendR, 'jstrR=',jstrR, 'jendR=',jendR
# ifdef MPI
      if (mynode < NNODES) then
        call MPI_Send (blank, 1, MPI_INTEGER, mynode+1,
     &                        1, ocean_grid_comm,  ierr)
      endif
# endif
#endif

      call init_arrays_ocean

      call init_avg_arrays

      call init_arrays_coupling

      call init_arrays_grid     ! lots of extra arrays here initialized which might negatively affect first touch principle.
                                ! before only rmask was first-touched here...

      call init_arrays_mess_buffers  ! DevinD placed here. Assumed used often for higher up first touch

      call init_arrays_eos_vars

#ifdef SOLVE3D
      do k=0,N
        do j=jstrR,jendR
          do i=istrR,iendR
            We(i,j,k)=init0
            Wi(i,j,k)=init0
# ifdef NHMG
            w(i,j,k,1)=init0
            w(i,j,k,2)=init0
# endif
          enddo
        enddo
      enddo

#endif /* SOLVE3D */

! Initialize forcing arrays.

      call init_arrays_surf_flx

#ifdef BULK_FRC
      call init_arrays_bulk_frc
#endif

# ifdef WEC
      call init_arrays_wec_tile(istr,iend,jstr,jend)
# endif

#if defined M2NUDGING && !defined M2_FRC_BRY
      do j=jstrR,jendR
        do i=istrR,iendR
          ssh(i,j)=init0
# ifndef ANA_SSH
          sshg(i,j,1)=init0
          sshg(i,j,2)=init0
# endif
        enddo
      enddo
#endif

! Set variable horizontal viscosities and tracer diffusion
! coefficients (see "mixing") to their background values.

      call init_arrays_mixing
      call init_arrays_work_mod   ! DevinD added this here to allocate
      call init_arrays_boundary
      call init_arrays_buffer     ! DevinD added this here to allocate

#ifndef BULK_FRC
      call init_arrays_flux_frc   ! DevinD added this here to allocate
#endif

      if (pipe_source) call init_arrays_pipes

#if defined(BIOLOGY_BEC2) || defined(MARBL)
      call init_arrays_bgc_ecosys_vars
      call init_arrays_bgc_forces
      call init_arrays_bgc_frc
#endif

      end
