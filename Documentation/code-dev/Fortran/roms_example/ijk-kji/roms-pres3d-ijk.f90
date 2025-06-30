! Loop taken from roms pre_step3d.F

program prestep3d_loop

    implicit none
    integer, parameter :: ntests = 10000
!    real(kind=8), allocatable, dimension(:,:,:) :: a,b,c
    real(kind=8) :: t1, t2, elapsed_time
    integer(kind=8) :: tclock1, tclock2, clock_rate
    integer :: i,j,k,itest
    ! ROMS
    integer, parameter :: N=50,x=100,y=80
    real(kind=8) :: FlxDiv, cff, dtau=3.
    real(kind=8), allocatable, dimension(:,:) :: pm, pn
    real(kind=8), allocatable, dimension(:,:,:) :: FlxU,FlxV,Hz_bak,Hz_fwd,Hz,We,Wi

    call system_clock(tclock1)

    print *, "pre_step3d loop, array size x,y,N=", x, y, N

    allocate(pm(-1:x+1,-1:y+1), pn(-1:x+1,-1:y+1))
    allocate(FlxU(-1:x+1,-1:y+1,0:N), FlxV(-1:x+1,-1:y+1,0:N))
    allocate(Hz_bak(-1:x+1,-1:y+1,0:N), Hz_fwd(-1:x+1,-1:y+1,0:N))
    allocate(We(-1:x+1,-1:y+1,0:N), Wi(-1:x+1,-1:y+1,0:N))
    allocate(Hz(-1:x+1,-1:y+1,0:N))
    ! fill a and b with 1's just for demo purposes:
    pm = 1.d0
    pn = 2.3
    FlxDiv = 0.d0
    FlxU = 1.5
    FlxV = 3.1
    We = 8.1
    Wi = 4.95
    Hz_bak = 2.8
    Hz_fwd = 5.01
    Hz = 1.12

    call cpu_time(t1)   ! start cpu timer

! ----------------------------------------------------------------
    do itest=1,ntests
      do i=1,x                       ! compressible predictor substep,
        cff=0.5*dtau                  ! Eq. (4.7) from SM2005.
        do j=1,y
          do k=1,N
            FlxDiv=cff*pm(i,j)*pn(i,j)*( FlxU(i+1,j,k)-FlxU(i,j,k) &
                                        +FlxV(i,j+1,k)-FlxV(i,j,k) &
                     +We(i,j,k)+Wi(i,j,k) -We(i,j,k-1)-Wi(i,j,k-1) &
                                                                 )
            Hz_bak(i,j,k)=Hz(i,j,k) +FlxDiv
            Hz_fwd(i,j,k)=Hz(i,j,k) -FlxDiv
          enddo
        enddo
      enddo
    enddo
! -----------------------------------------------------------------

    call cpu_time(t2)   ! end cpu timer
    print 10, ntests, t2-t1
 10 format("Performed ",i7, " loops: CPU time = ",f12.8, " seconds")


    call system_clock(tclock2, clock_rate)
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
    print 11, elapsed_time
 11 format("Elapsed time = ",f12.8, " seconds")

end program prestep3d_loop
