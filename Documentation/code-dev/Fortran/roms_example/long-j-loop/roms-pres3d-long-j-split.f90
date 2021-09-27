! Loop taken from roms pre_step3d.F
program main

    implicit none
    integer, parameter :: ntests = 1e4
    real(kind=8) :: t1, t2, elapsed_time
    integer(kind=8) :: tclock1, tclock2, clock_rate
    integer :: itest
    ! 1D array
    real(kind=8), allocatable, dimension(:) :: dd 
    integer, parameter :: N=50,x=80,y=100  
    ! 2d arrays
    real(kind=8), allocatable, dimension(:,:) :: pm,pn,DC,FC,WC,CF
    ! 3d arrays
    real(kind=8), allocatable, dimension(:,:,:) :: Hz_bak,Hz_fwd,rw,FlxU,FlxV,Hz,Akv,Wi
    ! 4d arrays
    real(kind=8), allocatable, dimension(:,:,:,:) :: w 
    
    ! 1D array
    allocate(dd(0:ntests))    
    ! 2D arrays
    allocate(pm(-1:x+1,-1:y+1), pn(-1:x+1,-1:y+1))
    allocate(DC(-1:x+1,0:N),FC(-1:x+1,0:N),WC(-1:x+1,0:N),CF(-1:x+1,0:N))
    ! 3D arrays
    allocate(FlxU(-1:x+1,-1:y+1,0:N), FlxV(-1:x+1,-1:y+1,0:N))
    allocate(Hz_bak(-1:x+1,-1:y+1,0:N), rw(-1:x+1,-1:y+1,0:N))
    allocate(Hz(-1:x+1,-1:y+1,0:N),Akv(-1:x+1,-1:y+1,0:N))
    allocate(Hz_fwd(-1:x+1,-1:y+1,0:N),Wi(-1:x+1,-1:y+1,0:N))
    ! 4d arrays
    allocate(w(-1:x+1,-1:y+1,0:N,2))
    
    ! fill a and b with 1's just for demo purposes:
    pm = 1.d0
    pn = 2.3
    FlxU = 1.5
    FlxV = 3.1   
    Hz_bak = 2.8
    Hz = 1.12
    Hz_fwd = -0.3
    rw = 0.11
    Akv = 12.1
    Wi = 3.4    
    
    ! -----------------------------------------------------------------

    call system_clock(tclock1)

    print *, "pre_step3d loop:"

    call cpu_time(t1)   ! start cpu timer

    dd(0)=0.1
    do itest=1,ntests
      dd(itest)=dd(itest-1)+0.001
      call prestep3d_loop(dd(itest),pm,pn,DC,FC,WC,CF,Hz_bak,Hz_fwd, &
                          rw,FlxU,FlxV,Hz,Akv,Wi,w)
    enddo

    call cpu_time(t2)   ! end cpu timer
    print 10, ntests, t2-t1
 10 format("Performed ",i7, " loops: CPU time = ",f12.8, " seconds")

    
    call system_clock(tclock2, clock_rate)
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
    print 11, elapsed_time
 11 format("Elapsed time = ",f12.8, " seconds")


end program main


subroutine prestep3d_loop( dd, pm,pn,DC,FC,WC,CF,Hz_bak,Hz_fwd, &
                           rw,FlxU,FlxV,Hz,Akv,Wi,w )

    implicit none

    ! Inputs
    real(kind=8) dd ! Scalar to multiply loop to ensure values are changing.

    integer :: i,j,k,itest
    ! ROMS
    integer, parameter :: N=50,x=80,y=100,istr=1,jstr=1,iend=x,indx=1,nstp=2
    ! Scalars
    real(kind=8) :: FlxDiv=0.d0,cff,dtau=3.,cf_stp=0.4,cf_bak=0.21
    ! 2d arrays
    real(kind=8), dimension(-1:x+1,-1:y+1) :: pm,pn
    real(kind=8), dimension(-1:x+1,0:N)    :: DC,FC,WC,CF
    ! 3d arrays
    real(kind=8), dimension(-1:x+1,-1:y+1,0:N) :: Hz_bak,Hz_fwd,rw,FlxU,FlxV,Hz,Akv,Wi
    ! 4d arrays
    real(kind=8), dimension(-1:x+1,-1:y+1,0:N,2) :: w
    
   
! ----------------------------------------------------------------    
      do j=1,y  !! Start of the giant j-loop
      
        do i=istr,iend
          DC(i,0) =dtau*pm(i,j)*pn(i,j)
        enddo
        
      enddo ! j 
      do j=1,y  !! Start j-loop 
              
        !! w(indx),w(nstp) in m/s, 
        !! w(0) is always zero
        do k=1,N-1
          do i=istr,iend
            DC(i,k)=dd*0.5*(Hz_bak(i,j,k+1)+Hz_bak(i,j,k))*( &
                       cf_stp*w(i,j,k,nstp)+cf_bak*w(i,j,k,indx) ) &
                                                + DC(i,0)*rw(i,j,k)

            w(i,j,k,indx)=dd*0.5*(Hz(i,j,k+1)+Hz(i,j,k))*w(i,j,k,nstp)
          enddo
        enddo
        
      enddo ! j 
      do j=1,y  !! Start j-loop 
              
        k = N   
        !! here is the special volume weighting for w(N)
        !! Still, write out the volume integrated (dz*w) vertical mixing eq.
        !! Just to check :)
        do i=istr,iend 
          DC(i,k)=dd*0.5*(   Hz_bak(i,j,k))*(  &
                     cf_stp*w(i,j,k,nstp)+cf_bak*w(i,j,k,indx) )  &
                                              + DC(i,0)*rw(i,j,k)

          w(i,j,k,indx)=dd*0.5*(  Hz(i,j,k))*w(i,j,k,nstp)
        enddo
        
      enddo ! j 
      do j=1,y  !! Start j-loop 
      
        !! start of tri-diag solve
        !   FC(i,1)=0.5*dtau*(Akt(i,j,k+1)+Akt(i,j,k)/( Hz_fwd(i,j,k) 
        do i=istr,iend
          FC(i,N-1)= dd*0.5*dtau*(-Akv(i,j,N)+Akv(i,j,N-1)) & ! see notes, exception for k=N
                           /Hz_fwd(i,j,N)

          WC(i,N-1)= DC(i,0)*0.5*(Wi(i,j,N)+Wi(i,j,N-1))

          cff = 0.5*Hz_fwd(i,j,N) +FC(i,N-1)-min(WC(i,N-1),0.)

          CF(i,N-1)= ( FC(i,N-1)+max(WC(i,N-1),0.) )/cff ! gam(n-1)=c(n)/bet 

          DC(i,N)= dd* DC(i,N)/cff                          ! u(1) = r(1)/bet
        enddo
        
      enddo ! j 
      do j=1,y  !! Start j-loop 
      
        do k=N-1,2,-1      !--> forward elimination for w
          do i=istr,iend
            FC(i,k-1)= dd*0.5*dtau*(Akv(i,j,k)+Akv(i,j,k-1))&!! A(k-1)
                        /Hz_fwd(i,j,k)

            WC(i,k-1)= dd*DC(i,0)*0.5*(Wi(i,j,k)+Wi(i,j,k-1))

            cff=dd*1./( 0.5*(Hz_fwd(i,j,k)+Hz_fwd(i,j,k-1))&
                                 +FC(i,k-1)-min(WC(i,k-1),0.)&
                                   +FC(i,k)+max(WC(i,k),0.)&
                          -CF(i,k)*(FC(i,k)-min(WC(i,k),0.))&
                                                            )
            CF(i,k-1)=dd*cff*( FC(i,k-1)+max(WC(i,k-1),0.) )

            DC(i,k)=dd*cff*(DC(i,k)+DC(i,k+1)*(FC(i,k)-min(WC(i,k),0.)))
          enddo
        enddo
        
      enddo ! j 
      do j=1,y  !! Start j-loop 
              
!        !! Use the fact that w(0) is always 0
        do i=istr,iend
          FC(i,0)= dd*0.5*dtau*(Akv(i,j,1)+Akv(i,j,0))&!! A(0)
                        /Hz_fwd(i,j,1)

          WC(i,0)= dd*DC(i,0)*0.5*(Wi(i,j,1)+Wi(i,j,0))

          w(i,j,1,indx)=dd*( DC(i,1)  +DC(i,2)*(FC(i,1)-min(WC(i,1),0.))&
                              )/( 0.5*(Hz_fwd(i,j,2)+Hz_fwd(i,j,1))&
                                            +FC(i,0)-min(WC(i,0),0.)&
                                            +FC(i,1)+max(WC(i,1),0.)&
                                   -CF(i,1)*(FC(i,1)-min(WC(i,1),0.))&
                                                                    )
        enddo
        
      enddo ! j 
      do j=1,y  !! Start j-loop 
              
        do k=2,N,+1          !--> backsubstitution for w
          do i=istr,iend
            w(i,j,k,indx)=dd*DC(i,k) +CF(i,k-1)*w(i,j,k-1,indx)
          enddo
        enddo
        !------- end computing of w(:,:,:,nnew) -----      
      
      enddo !! giant j-loop
! -----------------------------------------------------------------



end subroutine prestep3d_loop
