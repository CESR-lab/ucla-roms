! $UWHPSC/codes/fortran/optimize/timings.f90

! Illustrate timing utilities in Fortran.
!  system_clock can be used to compute elapsed time between
!      two calls (wall time)
!  cpu_time can be used to compute CPU time used between two calls.

! Try compiling with different levels of optimization, e.g. -O3


program timings

    implicit none
    integer, parameter :: ntests = 100
    integer :: n = 200
    real(kind=8), allocatable, dimension(:,:,:) :: a,b,c
    real(kind=8) :: t1, t2, elapsed_time
    integer(kind=8) :: tclock1, tclock2, clock_rate
    integer :: i,j,k,itest

    call system_clock(tclock1)

    print *, "Will multiply n by n matrices array statements, n= ", n

    allocate(a(n,n,n), b(n,n,n), c(n,n,n))

    ! fill a and b with 1's just for demo purposes:
    a = 1.d0
    b = 1.d0
    c = 1.d0

    call cpu_time(t1)   ! start cpu timer

    do itest=1,ntests
      c = c + a*b
    enddo

    call cpu_time(t2)   ! end cpu timer
    print 10, ntests, t2-t1
 10 format("Performed ",i4, " matrix multiplies: CPU time = ",f12.8, " seconds")

    
    call system_clock(tclock2, clock_rate)
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
    print 11, elapsed_time
 11 format("Elapsed time = ",f12.8, " seconds")

end program timings
