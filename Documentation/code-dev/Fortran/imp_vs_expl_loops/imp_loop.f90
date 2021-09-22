! $UWHPSC/codes/fortran/optimize/timings.f90

! Illustrate timing utilities in Fortran.
!  system_clock can be used to compute elapsed time between
!      two calls (wall time)
!  cpu_time can be used to compute CPU time used between two calls.

! Try compiling with different levels of optimization, e.g. -O3


program timings

    implicit none
    integer, parameter :: ntests = 20
    integer :: n 
    real(kind=8), allocatable, dimension(:,:,:) :: a,b,c
    real(kind=8) :: t1, t2, elapsed_time
    integer(kind=8) :: tclock1, tclock2, clock_rate
    integer :: i,j,k,itest
    integer :: dk=3, dj=2, di=4

    call system_clock(tclock1)

    print *, "Will multiply n by n matrices k-j-i, n= ", n
    print *, "implicit loop (array operations):"
    print *, "k -> dk:n-dk, where dk=", dk
    print *, "j -> dj:n-dj, where dj=", dj
    print *, "i -> di:n-di, where di=", di
!    read *, n
    n=780 ! 810 seemed to be from old results but was killed this time. Maybe since I was using other resources.

    allocate(a(n,n,n), b(n,n,n), c(n,n,n))

    ! fill a and b with 1's just for demo purposes:
    a = 1.d0
    b = 1.d0

    call cpu_time(t1)   ! start cpu timer
    do itest=1,ntests
!        do k = indk,n-indk
!            do j = indj,n-indj
!                do i=indi,n-indi
     c(di:n-di,dj:n-dj,dk:n-dk) = 1.d0
     c(di:n-di,dj:n-dj,dk:n-dk) = c(di:n-di,dj:n-dj,dk:n-dk) + a(di:n-di,dj:n-dj,dk:n-dk)*b(di:n-di,dj:n-dj,dk:n-dk)
!                enddo
!            enddo
!        enddo
    enddo

    call cpu_time(t2)   ! end cpu timer
    print 10, ntests, t2-t1
 10 format("Performed ",i4, " matrix multiplies: CPU time = ",f12.8, " seconds")

    
    call system_clock(tclock2, clock_rate)
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
    print 11, elapsed_time
 11 format("Elapsed time = ",f12.8, " seconds")

end program timings
