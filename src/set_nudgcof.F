#include "cppdefs.opt"
#if defined SPONGE || defined TCLIMATOLOGY

      subroutine set_nudgcof  ! Set up nudging-to-climatology

      use param
      use private_scratch

      implicit none                  ! inverse time-scale coeffcients
      integer tile                   ! (units [1/sec])

      call set_nudgcof_tile (A2d(1,1))
      end

      subroutine set_nudgcof_tile (wrk)

      use param
      use dimensions
      use mixing
      use scalars

      implicit none
      real wrk(PRIVATE_2D_SCRATCH_ARRAY)
      integer i,j, isp,ibnd,itrc
      real cff,cff1,cff2,cff3, frac, Tau

!
! The following code segment is rather a generic code go setup
! constant-slope, "flat-top" roof shaped nudging and sponge bands
! along open boundaries, if any. Adopted from Pierrick Penven.
!
      sponge_size = 15                      ! number of interior points in sponge layer (excludes buffer)
      isp = sponge_size+1                   ! ROMS uses isp which includes buffer point

      if (mynode==0) write(*,'(/1x,A,I3,A)')
     &  'Sponge layer = ',sponge_size,' interior grid points.'

      ! GS - for 2.0, these are set to +1 to be safe

                                            ! Note that i+iSW_corn and
      do j=-1,ny+1                          ! j+jSW_corn are indices on
        do i=-1,nx+1                        ! the physical grid (a.k.a.
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

      ! GS: set to 0 per usage of visc2_r in visc3d_GP and visc3d_S
      !     set to nx/ny+1 per usage of diff2 in t3dmix_S
      do j=0,ny+1
        do i=0,nx+1
#  ifdef SPONGE
#   ifdef UV_VIS2
          visc2_r(i,j)=visc2_r(i,j)+v_sponge*wrk(i,j)
#   endif
#   if defined SOLVE3D && defined TS_DIF2
#    ifdef PACIFIC
          diff2(i,j,itemp)=diff2(i,j,itemp)+0.5*v_sponge*wrk(i,j)
#    else
          diff2(i,j,itemp)=diff2(i,j,itemp)+v_sponge*wrk(i,j)
#    endif
#   endif
#  endif
        enddo
      enddo

#  ifdef SPONGE
#   ifdef UV_VIS2
      ! GS: set to 1 to nx/ny+1 per usage in visc2d
      do j=1,ny+1            ! viscosity at psi points
        do i=1,nx+1
          visc2_p(i,j)=visc2_p(i,j)+0.25*v_sponge*( wrk(i,j)
     &                  +wrk(i-1,j)+wrk(i,j-1)+wrk(i-1,j-1))
        enddo
      enddo
#   endif
#   if defined SOLVE3D && defined TS_DIF2
      do itrc=2,NT           ! diffusivity for the other tracers
      ! GS: set to 0 to nx/ny+1 per usage of diff2 in t3dmix_S
        do j=0,ny+1
          do i=0,nx+1
            diff2(i,j,itrc)=diff2(i,j,itemp)  ! Waste of space to have this 2D array copied...
          enddo
        enddo
      enddo
#   endif
#   ifdef UV_VIS2
#    ifdef EXCHANGE
!     call exchange_xxx(visc2_p)
#    endif
#   endif
#  endif
      end
#else
      subroutine set_nudgcof_empty
      end
#endif
