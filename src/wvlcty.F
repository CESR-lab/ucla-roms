      subroutine wvlcty(Wvlc)

! Compute absolute vertical velocity, which consists of three
! components: S-coordinate vertical velocity  w*pm*pn; projection
! of (quasi-)horizontal motions along S=const surfaces; and vertical
! velocity of moving grid-box interfaces due to the motion of free
! surface. This computation is done solely for diagnostic/output
! purposes and does not have any feedback onto the model.
! Unlike omega, absolute vertical velocity is defined at RHO-points.

!     use param
      use dimensions
      use grid
      use ocean_vars
      use scalars
      use roms_mpi

      implicit none

      ! import/export
      real,dimension(1-bf:nx+bf,1-bf:ny+bf,nz) :: wvlc
      ! local
      integer :: i,j,k

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

      do j=j0,j1
!       do i=i0,i1
!         wvlc(i,j,0)=0.
!       enddo
        do i=i0,i1
          wvlc(i,j,1) = -pm(i,j)*pn(i,j)*(
     &             FlxU(i+1,j,1)-FlxU(i,j,1)
     &            +FlxV(i,j+1,1)-FlxV(i,j,1))
        enddo
        do k=2,nz
          do i=i0,i1
            wvlc(i,j,k) = wvlc(i,j,k-1)-pm(i,j)*pn(i,j)*(
     &                      FlxU(i+1,j,k)-FlxU(i,j,k)
     &                     +FlxV(i,j+1,k)-FlxV(i,j,k))
          enddo
        enddo
        do i=i0,i1
          wvlc(i,j,1) = 0.5*(wvlc(i,j,1) )
        enddo
        do k=2,nz
          do i=i0,i1
            wvlc(i,j,k) = 0.5*(wvlc(i,j,k)+wvlc(i,j,k-1))
          enddo
        enddo
      enddo ! j=loop

! Compute and add contributions due to quasi-horizontal motions along
! S=const surfaces by multiplying horizontal velocity components by
! slops S-coordinate surfaces:

      do j=j0,j1
        do k=1,n
          do i=i0,i1
            wvlc(i,j,k) = wvlc(i,j,k) 
     &         + 0.25*u(i  ,j,k,nstp)*(z_r(i,j,k)-z_r(i-1,j,k))
     &                                 *(pm(i,j)+pm(i-1,j))
     &         + 0.25*u(i+1,j,k,nstp)*(z_r(i,j,k)-z_r(i+1,j,k))
     &                                 *(pm(i,j)+pm(i+1,j))
     &         + 0.25*v(i,j  ,k,nstp)*(z_r(i,j,k)-z_r(i,j-1,k))
     &                                  *(pn(i,j)+pn(i,j-1))
     &         + 0.25*v(i,j+1,k,nstp)*(z_r(i,j,k)-z_r(i,j+1,k))
     &                                  *(pn(i,j)+pn(i,j+1))
           enddo
         enddo
       enddo ! j=loop

      if (nhmg) then
       ! overwrite with w var
      do k=1,nz  ! WASTED COMPUTING TO OVERWRITE SOMETHING USE IF ELSE RATHER
        do j=j0,j1
          do i=i0,i1
            Wvlc(i,j,k) = w(i,j,k,nstp)
          enddo
        enddo
      enddo
      endif

!  Set lateral boundary conditions: no gradient for wvlc only.

      if (inode.eq.0.and.(.not.ew_periodic)) then
        wvlc(0,j0:j1,:) = wvlc(1,j0:j1,:)
      endif
      if (inode.eq.npx-1.and.(.not.ew_periodic)) then
        wvlc(nx+1,j0:j1,:) = wvlc(nx,j0:j1,:)
      endif
      if (jnode.eq.npy-1.and.(.not.ns_periodic)) then
        wvlc(i0:i1,ny+1,:) = wvlc(i0:i1,ny,:)
      endif
      if (jnode.eq.0.and.(.not.ns_periodic)) then
        wvlc(i0:i1,0,:) = wvlc(i0:i1,1,:)
      endif

      end
