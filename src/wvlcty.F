#include "cppdefs.opt"
#ifdef SOLVE3D

      subroutine wvlcty (tile, Wvlc)
      use param
      use hidden_mpi_vars
      use private_scratch
      implicit none
      integer tile
      real Wvlc(GLOBAL_2D_ARRAY,0:N)
# include "compute_tile_bounds.h"
      call wvlcty_tile (istr,iend,jstr,jend, Wvlc, A2d(1,1), A2d(1,1),   ! looks like a bug here repeated A2d(1,1)?
     &                                                       A2d(1,2))
      end

      subroutine wvlcty_tile (istr,iend,jstr,jend, Wvlc, Wrk,Wxi,Weta)

! Compute absolute vertical velocity, which consists of three
! components: S-coordinate vertical velocity  w*pm*pn; projection
! of (quasi-)horizontal motions along S=const surfaces; and vertical
! velocity of moving grid-box interfaces due to the motion of free
! surface. This computation is done solely for diagnostic/output
! purposes and does not have any feedback onto the model.
! Unlike W, absolute vertical velocity is defined at RHO-points.

      use param
      use grid
      use ocean_vars
      use scalars

      implicit none
      integer istr,iend,jstr,jend, imin,imax,jmin,jmax, i,j,k
      real Wvlc(GLOBAL_2D_ARRAY,0:N)
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N) :: Wrk
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: Wxi, Weta

# include "compute_extended_bounds.h"

# ifdef EW_PERIODIC
      imin=max(istrR, istr-1)        ! The logic here is somewhat
      imax=min(iendR, iend+1)        ! opposite to the standard
#else
      imin=istr                      ! "auxiliary" bounds: extend
      imax=iend                      ! the range of computation of
# endif
# ifdef NS_PERIODIC
      jmin=max(jstrR, jstr-1)        ! "Wvlc", by one row of points
      jmax=min(jendR, jend+1)        ! on the side, but only in the
# else
      jmin=jstr                      ! case of periodic direction.
      jmax=jend
# endif

! Compute "omega" vertical velocity by means of integration of mass
! divergence of mass fluxes from bottom up. In this computation,
! unlike that in omega.F, there is (1) immediate multiplication by
! pm*pn so that the result has meaning of velocity, rather than
! finite volume mass flux through vertical facet of tracer grid box;
! and (2, also unlike omega.F) no subtraction of vertical velocity
! of moving grid-box interface (the effect of "breething" of vertical
! grid system due to evolving free surface) is made now.
! Consequently, Wrk(:,N)/=0, unlike its counterpart W(:,:,N).eqv.0
! in omega.F. Once omega vertical velocity is computed, interpolate
! it to vertical RHO-points.

      do j=jmin,jmax
        do i=imin,imax
          Wrk(i,0)=0.
        enddo
        do k=1,N,+1
          do i=imin,imax
            Wrk(i,k)=Wrk(i,k-1)-pm(i,j)*pn(i,j)*(
     &                      FlxU(i+1,j,k)-FlxU(i,j,k)
     &                     +FlxV(i,j+1,k)-FlxV(i,j,k))
c**    Wrk(i,k)=0.!(uncomment to test the second part)
          enddo
        enddo
        do i=imin,imax
          Wvlc(i,j,N)=+0.375*Wrk(i,N) +0.75*Wrk(i,N-1)
     &                                -0.125*Wrk(i,N-2)
        enddo
        do k=N-1,2,-1
          do i=imin,imax
            Wvlc(i,j,k)=+0.5625*(Wrk(i,k  )+Wrk(i,k-1))
     &                  -0.0625*(Wrk(i,k+1)+Wrk(i,k-2))
          enddo
        enddo
        do i=imin,imax
          Wvlc(i,j,  1)= -0.125*Wrk(i,2) +0.75*Wrk(i,1)
     &                                  +0.375*Wrk(i,0)
        enddo
      enddo

! Compute and add contributions due to quasi-horizontal motions along
! S=const surfaces by multiplying horizontal velocity components by
! slops S-coordinate surfaces:

      do k=1,N
        do j=jmin,jmax
          do i=imin,imax+1
            Wxi(i,j)=u(i,j,k,nstp)*(z_r(i,j,k)-z_r(i-1,j,k))
     &                                  *(pm(i,j)+pm(i-1,j))
          enddo
        enddo
        do j=jmin,jmax+1
          do i=imin,imax
            Weta(i,j)=v(i,j,k,nstp)*(z_r(i,j,k)-z_r(i,j-1,k))
     &                                   *(pn(i,j)+pn(i,j-1))
          enddo
        enddo
        do j=jmin,jmax
          do i=imin,imax
            Wvlc(i,j,k)=Wvlc(i,j,k)+0.25*(Wxi(i,j)+Wxi(i+1,j)
     &                                +Weta(i,j)+Weta(i,j+1))
          enddo
        enddo
      enddo
# ifdef NHMG
       ! overwrite with w var
      do k=1,N  ! WASTED COMPUTING TO OVERWRITE SOMETHING USE IF ELSE RATHER
        do j=jmin,jmax
          do i=imin,imax
            Wvlc(i,j,k) =  w(i,j,k,nstp)
          enddo
        enddo
      enddo
# endif

!  Set lateral boundary conditions: gradient only.

# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do k=1,N
          do j=jmin,jmax
            Wvlc(imin-1,j,k)=Wvlc(imin,j,k)
          enddo
        enddo
      endif
      if (EASTERN_EDGE) then
        do k=1,N
          do j=jmin,jmax
            Wvlc(imax+1,j,k)=Wvlc(imax,j,k)
          enddo
        enddo
      endif
# endif
# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do k=1,N
          do i=imin,imax
            Wvlc(i,jmin-1,k)=Wvlc(i,jmin,k)
          enddo
        enddo
      endif
      if (NORTHERN_EDGE) then
        do k=1,N
          do i=imin,imax
            Wvlc(i,jmax+1,k)=Wvlc(i,jmax,k)
          enddo
        enddo
      endif
#  ifndef EW_PERIODIC
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          Wvlc(imin-1,jmin-1,k)=Wvlc(imin,jmin,k)
        enddo
      endif
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          Wvlc(imin-1,jmax+1,k)=Wvlc(imin,jmax,k)
        enddo
      endif
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          Wvlc(imax+1,jmin-1,k)=Wvlc(imax,jmin,k)
        enddo
      endif
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          Wvlc(imax+1,jmax+1,k)=Wvlc(imax,jmax,k)
        enddo
      endif
#  endif
# endif
      end
#else
      subroutine wvlcty_empty
      end
#endif /* SOLVE3D */
