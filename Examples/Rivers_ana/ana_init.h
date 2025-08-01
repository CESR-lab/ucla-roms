      integer :: i,j,k

      ! necessary for FIRST_TIME_STEP flag to work, as per set_global_definitions.h
      ! this had no value for analytical examples previously.
#ifdef EXACT_RESTART
      forw_start=ntstart
#endif


      do j=1-bf,ny+bf          ! Set everything (except temperature
        do i=1-bf,nx+bf        ! and salinity) to all-zero state, then
          zeta(i,j,1)=0.       ! modify some of the variables, if a
          ubar(i,j,1)=0.       ! non-trivial initialization required.
          vbar(i,j,1)=0.       ! Note: A code to initialize T [and S]
        enddo                  ! must always be supplied for 3D
      enddo                    ! applications.
# ifdef SOLVE3D
      do k=1,nz
        do j=1-bf,ny+bf      
          do i=1-bf,nx+bf    
            u(i,j,k,1)=0.
            u(i,j,k,2)=0.
            v(i,j,k,1)=0.
            v(i,j,k,2)=0.
          enddo
        enddo
      enddo
# endif

#  ifdef SOLVE3D
      do k=1,nz
        do j=1-bf,ny+bf      
          do i=1-bf,nx+bf    
            t(i,j,k,1,itemp)=4.+10.*exp(z_r(i,j,k)/50.)
            t(i,j,k,1,isalt)=36.
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
          enddo
        enddo
      enddo
#  endif
