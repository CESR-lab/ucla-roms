#include "cppdefs.opt"

      subroutine init_arrays(tile,my_first)

      use param
      use hidden_mpi_vars
      use private_scratch

      implicit none

      integer, intent(in) :: tile, my_first
      integer             :: i, j

#include "compute_tile_bounds.h"

      call init_arrays_private_scratch
#ifdef SOLVE3D
      do i=1,N2d
        iA2d(i,1)=0
        iA2d(i,2)=0
      enddo
#endif

      call init_arrays_tile(istr,iend,jstr,jend,tile,my_first)
      end

      subroutine init_arrays_tile (istr,iend,jstr,jend,tile,my_first)

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
      use bulk_frc  !, only: srflx
      use surf_flux ! for: sustr, svstr, stflx, sst_data, sss_data
      use tracers   ! iTandS, t, t_avg
      use averages
      use climat
      use coupling
      use eos_vars
      use grid
      use mixing
      use ocean2d
      use ocean3d
      use scalars
      use mpi
      use work_mod
      use boundary
      use buffer
      use mess_buffers
      use dimensions
      use flux_frc

      implicit none
      integer, intent(in) :: istr, iend, jstr, jend, tile, my_first
      integer             :: i, j, k, itrc, itavg
      real, parameter     :: init=0.    !!!!  0xFFFA5A5A ==> NaN
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

      call init_arrays_ocean2d

      call init_arrays_coupling

      call init_arrays_grid     ! lots of extra arrays here initialized which might negatively affect first touch principle.
                                ! before only rmask was first-touched here...
      call init_arrays_ocean3d

      call init_arrays_mess_buffers  ! DevinD placed here. Assumed used often for higher up first touch

#ifdef SOLVE3D
      do k=1,N                       ! Initialize
        do j=jstrR,jendR             ! 3-D primitive
          do i=istrR,iendR
# ifdef AVERAGES
            rho_avg(i,j,k)=init
# endif
          enddo
        enddo
      enddo

      call init_arrays_eos_vars

      do k=0,N
        do j=jstrR,jendR
          do i=istrR,iendR
            We(i,j,k)=init
            Wi(i,j,k)=init
# ifdef NHMG
            w(i,j,k,1)=init
            w(i,j,k,2)=init
# endif
# ifdef AVERAGES
            w_avg(i,j,k)=init
            wvl_avg(i,j,k)=init
# endif
          enddo
        enddo
      enddo

#endif /* SOLVE3D */

! Initialize forcing arrays.

      call init_arrays_surf_flx

# ifdef WEC
      call init_arrays_wec_tile(istr,iend,jstr,jend)
# endif

! Initialize climatology arrays (see "climat").

#if defined M2NUDGING && !defined M2_FRC_BRY
      do j=jstrR,jendR
        do i=istrR,iendR
          ssh(i,j)=init
# ifndef ANA_SSH
          sshg(i,j,1)=init
          sshg(i,j,2)=init
# endif
        enddo
      enddo
#endif
#ifdef SOLVE3D
# ifdef TCLIMATOLOGY
      do itrc=1,NT
        do j=jstrR,jendR
          do i=istrR,iendR
            Tnudgcof(i,j,itrc)=init
          enddo
        enddo
      enddo
# endif
# if (defined TCLIMATOLOGY || defined TNUDGING) && !defined ANA_TCLIMA
      do itrc=1,NT
        do k=1,N
          do j=jstrR,jendR
            do i=istrR,iendR
              tclm(i,j,k,itrc)=init
              tclima(i,j,k,1,itrc)=init
              tclima(i,j,k,2,itrc)=init
            enddo
          enddo
        enddo
      enddo
# endif
#endif

#ifdef UCLIMATOLOGY
# ifndef ANA_UCLIMA
      do j=jstrR,jendR
        do i=istrR,iendR
          ubclm(i,j)=init
          vbclm(i,j)=init
          ubclima(i,j,1)=init
          ubclima(i,j,2)=init
          vbclima(i,j,1)=init
          vbclima(i,j,2)=init
        enddo
      enddo
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            uclm(i,j,k)=init
            vclm(i,j,k)=init
            uclima(i,j,k,1)=init
            uclima(i,j,k,2)=init
            vclima(i,j,k,1)=init
            vclima(i,j,k,2)=init
          enddo
        enddo
      enddo
#  endif
# endif
#endif

! Set variable horizontal viscosities and tracer diffusion
! coefficients (see "mixing") to their background values.

      call init_arrays_mixing
      call init_arrays_dimensions ! DevinD added this here to allocate
      call init_arrays_work_mod   ! DevinD added this here to allocate
      call init_arrays_boundary
      call init_arrays_buffer     ! DevinD added this here to allocate

#ifndef BULK_FRC
      call init_arrays_flux_frc   ! DevinD added this here to allocate
#endif

      end
