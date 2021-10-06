#include "cppdefs.opt"
#if defined SPONGE || defined TCLIMATOLOGY

      subroutine set_nudgcof (tile)  ! Set up nudging-to-climatology

      use param
      use grid
      use hidden_mpi_vars
      use private_scratch

      implicit none                  ! inverse time-scale coeffcients
      integer tile                   ! (units [1/sec])

# include "compute_tile_bounds.h"
      call set_nudgcof_tile (istr,iend,jstr,jend, A2d(1,1))
      end

      subroutine set_nudgcof_tile (istr,iend,jstr,jend, wrk)

      use param
      use climat
      use mixing
      use scalars

      implicit none
      integer istr,iend,jstr,jend
      real wrk(PRIVATE_2D_SCRATCH_ARRAY)
      integer i,j, isp , I1,I0,I2,J1,J0,J2, ibnd,itrc
      real cff,cff1,cff2,cff3, frac, Tau

# include "compute_extended_bounds.h"

# if defined DAMEE_B
      Tau=1./(5.*86400.)              ! Make two nudging bands along
      do j=jstrR,jendR                ! southern and northern borders
        do i=istrR,iendR              ! of the domain.
          Tnudgcof(i,j,itemp)=0.
        enddo
#  ifdef MPI
        frac=(float(j+jSW_corn)-0.5)/float(MMm)
#  else
        frac=(float(j)         -0.5)/float(MMm)
#  endif
        cff=Tau*(179./168. - frac*352./21.)      ! All these fancy
        if (cff>0.) then                      ! looking rational
          do i=istrR,iendR                       ! fractions are just
            Tnudgcof(i,j,itemp)=cff              ! for the reasons of
          enddo                                  ! compatibility with
        endif                                    ! legacy code and do
        cff=Tau*(179./168. - (1.-frac)*352./21.) ! carry any mystical
        if (cff>0.) then                      ! meaning.
          do i=istrR,iendR
            Tnudgcof(i,j,itemp)=cff
          enddo
        endif
      enddo

      cff2=Tau/12.
      if (LLm == 128) then            ! Make special semi-circular
        cff3=Tau*95./84.              ! nudjing area to simulate
        frac=0.125                    ! Mediterranean outflow.
        I1=102                        ! I0,J0 are index coordinates
        I0=109                        ! of the center located at
        I2=108                        ! Gibraltar, I1..I2, J1..J2
        J1=74                         ! cut out the needed portion
        J0=77                         ! of the circle.
        J2=80
      elseif (LLm==256) then
        cff3=Tau*107./96.
        frac=0.1
        I1=229
        I0=238
        I2=237
        J1=154
        J0=159
        J2=164
      endif

#  ifdef MPI
      I1=I1 -iSW_corn
      I0=I0 -iSW_corn
      I2=I2 -iSW_corn

      J1=J1 -jSW_corn
      J0=J0 -jSW_corn
      J2=J2 -jSW_corn
#  endif
      do j=max(jstrR,J1),min(J2,jendR)
        do i=max(istrR,I1),min(I2,iendR)
          cff=frac*sqrt(float((i-I0)*(i-I0)+(j-J0)*(j-J0)))
          Tnudgcof(i,j,itemp)=max(0., cff3+cff*(cff2-cff3) )
        enddo
      enddo

      do j=jstrR,jendR
        do i=istrR,iendR
          Tnudgcof(i,j,isalt)=Tnudgcof(i,j,itemp)
        enddo
      enddo

c--# elif defined PACIFIC
c      Tau=1./(5.*86400.)              ! Make nudging band along
c      do j=jstrR,jendR                ! southern boundary
c--#  ifdef MPI
c        frac=(float(j+jSW_corn)-0.5)/float(MMm)
c--#  else
c        frac=(float(j)         -0.5)/float(MMm)
c--#  endif
c        cff=Tau*(1.-frac/16.)
c        if (cff>0.) then
c          do i=istrR,iendR
c            Tnudgcof(i,j,itemp)=cff
c          enddo
c        else
c          do i=istrR,iendR
c            Tnudgcof(i,j,itemp)=0.
c          enddo
c        endif
c      enddo
# else
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

      do j=jstrR,jendR
        do i=istrR,iendR
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
c        if (diff2(i,j,itemp)/=500.D0) then
c        write(*,*) i,j,diff2(i,j,itemp)
c        endif
#   endif
#  endif
#  if defined SOLVE3D && defined TCLIMATOLOGY
          Tnudgcof(i,j,itemp)=tauT_out*wrk(i,j)
#  endif
        enddo
      enddo

#  if defined SOLVE3D && defined TCLIMATOLOGY
      do itrc=2,NT               ! includes BIOLOGY variables
        do j=jstrR,jendR
          do i=istrR,iendR
            Tnudgcof(i,j,itrc)=Tnudgcof(i,j,itemp)
          enddo
        enddo
      enddo
#   if defined ROBUST_DIAG
      do j=jstrR,jendR         ! nudging everywhere
        do i=istrR,iendR
          Tnudgcof(i,j,itemp)=1./(360.*86400.)
#     if defined SALINITY
          Tnudgcof(i,j,isalt)=1./(360.*86400.)
#     endif
        enddo
      enddo
#   endif
#  endif

#  ifdef SPONGE
#   ifdef UV_VIS2
      do j=jstr,jendR        ! viscosity at psi points
        do i=istr,iendR
          visc2_p(i,j)=visc2_p(i,j)+0.25*v_sponge*( wrk(i,j)
     &                  +wrk(i-1,j)+wrk(i,j-1)+wrk(i-1,j-1))
        enddo
      enddo
#   endif
#   if defined SOLVE3D && defined TS_DIF2
      do itrc=2,NT           ! diffusivity for the other tracers
        do j=jstrR,jendR
          do i=istrR,iendR
            diff2(i,j,itrc)=diff2(i,j,itemp)  ! Waste of space to have this 2D array copied...
          enddo
        enddo
      enddo
#   endif
#   ifdef UV_VIS2
#    ifdef EXCHANGE
        call exchange2d_tile (istr,iend,jstr,jend, visc2_p)
#    endif
#   endif
#  endif
# endif
      end
#else
      subroutine set_nudgcof_empty
      end
#endif
