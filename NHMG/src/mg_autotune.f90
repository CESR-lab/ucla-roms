module mg_autotune

  !- NHMG AUTOTUNE -!
  ! If the namelist parameter autotune is .True.
  ! test different values of namelist parameters to solve an output of CROCO 
  ! after a number of time steps (100 by default).
  ! Results, in terms of performance, are ranked and the best set of values underligned.
  ! 
  ! Nov 2017 : test only ns_pre and ns_post
  !
  ! - NHMG AUTOTUNE -!

  use mg_cst
  use mg_mpi
  use mg_namelist
  use mg_tictoc
  use mg_projection
  use mg_solvers

  implicit none

  integer(kind=rp) :: pre_min = 1
  integer(kind=rp) :: pre_max = 4
  integer(kind=rp) :: post_min = 1
  integer(kind=rp) :: post_max = 4

contains

  subroutine sb_autotune()

    character(len=64) :: testname
    character(len=16) :: cpre, cpost

    integer(kind=ip) :: orig_pre
    integer(kind=ip) :: orig_post

    integer(kind=ip) :: ipre, ipost

    integer(kind=ip) :: nx,ny,nz

    real(kind=rp), dimension(:,:)  , allocatable :: res
    real(kind=rp), dimension(:,:,:), allocatable :: pres

    if (myrank==0) write(*,*)'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
    if (myrank==0) write(*,*)'>>>>>> STARTING AUTOTUNE <<<<<<<'
    if (myrank==0) write(*,*)'<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'

    orig_pre  = ns_pre
    orig_post = ns_post

    allocate(res(pre_max-pre_min+1,post_max-post_min+1))

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz
    allocate(pres(nz,0:ny+1,0:nx+1))

    pres = grid(1)%p

    do ipre = pre_min, pre_max

       do ipost = post_min, post_max

          write(unit=cpre ,fmt=*) ipre
          write(unit=cpost,fmt=*) ipost

          testname = 'autotune_' // trim(adjustl(cpre)) // trim(adjustl(cpost))

          if (myrank==0) write(*,*) trim(testname)

          ns_pre = ipre   !- ns_pre  is declared in namelist.f90
          ns_post = ipost !- ns_post is declared in namelist.f90

          grid(1)%p = 0._rp

          call tic(1,trim(testname))

          !- WARNING : 
          !- Array values used in solve_p have to be exactly the same at each test.
          !- List of arrays:
          !-      -
          !-      -
          call solve_p()

          call toc(1,trim(testname))

          call get_tictoc(res(ipre,ipost),1,trim(testname))

          if (myrank==0) write(*,*) trim(testname),' results:', res(ipre,ipost)

          grid(1)%p = pres

       enddo ! pre

    enddo ! post

    deallocate(pres)

    if (myrank==0) write(*,*)'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
    if (myrank==0) write(*,*)'>>>>>> FINISHING AUTOTUNE <<<<<<'
    if (myrank==0) write(*,*)'>>>> AND STOP THE PROGRAM  <<<<<'
    if (myrank==0) write(*,*)'<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'

    stop

  end subroutine sb_autotune

end module mg_autotune
