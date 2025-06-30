
      ! Everything after the implicit none
      ! Replace with something less trivial when needed

      integer :: i,j,k


      ! necessary for FIRST_TIME_STEP flag to work, as per set_global_definitions.h
      ! this had no value for analytical examples previously.
#ifdef EXACT_RESTART
      forw_start=ntstart
#endif


      do k=1,nz
        do j=-1,ny+2
          do i=-1,nx+2
            t(i,j,k,1,itemp)= 20
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
# ifdef SALINITY
            t(i,j,k,1,isalt)=36.
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
# endif
          enddo
        enddo
      enddo

      do j=0,ny+1
        do i=0,nx+1
          ubar(i,j,1)=0.
          ubar(i,j,2)=ubar(i,j,1)
          vbar(i,j,1)=0
          vbar(i,j,2)=vbar(i,j,1)
          zeta(i,j,1)=9.
          zeta(i,j,2)=zeta(i,j,1)
          do k=1,nz
            u(i,j,k,1)=0.
            u(i,j,k,2)=u(i,j,k,1)
            v(i,j,k,1)=0.
            v(i,j,k,2)=v(i,j,k,1)
          enddo
        enddo
      enddo
