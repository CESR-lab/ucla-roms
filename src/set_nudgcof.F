#include "cppdefs.opt"
#if defined SPONGE || defined TCLIMATOLOGY

      subroutine set_nudgcof (tile)  ! Set up nudging-to-climatology

      use param
      use grid
      use hidden_mpi_vars
      use private_scratch
      use mpi_exchanges

      implicit none                  ! inverse time-scale coeffcients
      integer tile                   ! (units [1/sec])

# include "compute_tile_bounds.h"
      call set_nudgcof_tile (istr,iend,jstr,jend, A2d(1,1))
      end

      subroutine set_nudgcof_tile (istr,iend,jstr,jend, wrk)

      use param
      use mixing
      use scalars

      implicit none
      integer istr,iend,jstr,jend
      real wrk(PRIVATE_2D_SCRATCH_ARRAY)
      integer i,j, isp , I1,I0,I2,J1,J0,J2, ibnd,itrc
      real cff,cff1,cff2,cff3, frac, Tau

# include "compute_extended_bounds.h"

!
! The following code segment is rather a generic code go setup
! constant-slope, "flat-top" roof shaped nudging and sponge bands
! along open boundaries, if any. Adopted from Pierrick Penven.
!
      sponge_size = 15                      ! number of interior points in sponge layer (excludes buffer)
      isp = sponge_size+1                   ! ROMS uses isp which includes buffer point

      if (mynode==0) write(*,'(/1x,A,I3,A)')
     &  'Sponge layer = ',sponge_size,' interior grid points.'

                                            ! Note that i+iSW_corn and
      do j=max(-1,jstrR-1),jendR            ! j+jSW_corn are indices on
        do i=max(-1,istrR-1),iendR          ! the physical grid (a.k.a.
          ibnd=isp                          ! absolute indices).
#  ifdef MPI
#   ifdef OBC_WEST
          ibnd=min(ibnd,         i+iSW_corn )
#   endif
#   ifdef OBC_EAST
          ibnd=min(ibnd, LLm+1 -(i+iSW_corn))
#   endif
#   ifdef OBC_SOUTH
          ibnd=min(ibnd,         j+jSW_corn )
#   endif
#   ifdef OBC_NORTH
          ibnd=min(ibnd, MMm+1 -(j+jSW_corn))
#   endif
#  else
#   ifdef OBC_WEST
          ibnd=min(ibnd,      i)
#   endif
#   ifdef OBC_EAST
          ibnd=min(ibnd, Lm+1-i)
#   endif
#   ifdef OBC_SOUTH
          ibnd=min(ibnd,      j)
#   endif
#   ifdef OBC_NORTH
          ibnd=min(ibnd, Mm+1-j)
#   endif
#  endif
#  ifdef PACIFIC
          wrk(i,j)=0.5 + 0.5*cos(pi*float(ibnd)/float(isp))
#  else
          wrk(i,j)=float(isp-ibnd)/float(isp)  ! 1 in interior, <1 within sponge
#  endif
        enddo
      enddo

#  ifdef SPONGE
#   ifdef UV_VIS2
      do j=jstrR,jendR
        do i=istrR,iendR
          visc2_r(i,j)=visc2_r(i,j)+v_sponge*wrk(i,j)
        enddo
      enddo

      do j=jstr,jendR        ! viscosity at psi points
        do i=istr,iendR
          visc2_p(i,j)=visc2_p(i,j)+0.25*v_sponge*( wrk(i,j)
     &                  +wrk(i-1,j)+wrk(i,j-1)+wrk(i-1,j-1))
        enddo
      enddo
#   endif

#   if defined SOLVE3D && defined TS_DIF2
      do itrc=1,NT           ! diffusivity for the tracers
        do j=jstrR,jendR
          do i=istrR,iendR
            diff2(i,j,itrc)=diff2(i,j,itrc)+v_sponge*wrk(i,j)
          enddo
        enddo
      enddo
#   endif
#  endif
      end
#else
      subroutine set_nudgcof_empty
      end
#endif
