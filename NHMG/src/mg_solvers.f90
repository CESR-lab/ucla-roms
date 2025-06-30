module mg_solvers

  use mg_cst
  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_intergrids
  use mg_relax
  use mg_netcdf_out

  implicit none

contains

  !---------------------------------------------------------------------
  subroutine solve_p()

    integer(kind=ip) :: nite
    integer(kind=ip) :: nx,ny,nz

    real(kind=rp)    :: rnorm,bnorm,res0,conv,rnorm0
    real(kind=rp), dimension(:,:,:), pointer :: p,b,r

    real(kind=lg) :: tstart,tend,perf
    real(kind=rp) :: rnxg,rnyg,rnzg
    real(kind=rp) :: rnpxg,rnpyg,fac

    integer, save :: count = 1
    logical :: verbose

    verbose = .false.
    if (mod(count,output_freq)==0) verbose=.true.
    if (autotune)  verbose=.true.
    count = count+1
    if ((myrank==0).and.verbose) write(*,*)'     ---------------'

    ! if commented, we use the previous pressure as first guess for the current projection
    ! grid(1)%p(:,:,:) = zero

    p  => grid(1)%p
    b  => grid(1)%b
    r  => grid(1)%r

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    call tic(1,'solve')
    call cpu_time(tstart)

    nite=0

    res0 = sum(grid(1)%b(1:nz,1:ny,1:nx)**2)
    call global_sum(1,res0,bnorm)
    bnorm = sqrt(bnorm)
    fac = sqrt(1./(nx*ny*nz))

    ! residual returns both 'r' and its norm
    call compute_residual(1,rnorm)
    res0   = rnorm/bnorm
    rnorm0 = res0

!   JM July 2020
!   solver_prec indicated an error reduction rather than an absolute error.
!   This is a problem for small bnorm (small absolute incoming divergence)
!   I added an absolute error criterium
    do while ((nite < solver_maxiter).and.(res0 > solver_prec).and.(res0*bnorm*fac>1e-15) )

       call Fcycle()

       call compute_residual(1,rnorm)
       rnorm = rnorm/bnorm
       conv = res0/rnorm ! error reduction after this iteration
       res0 = rnorm

       nite = nite+1
       if ((myrank == 0).and.verbose) write(*,10) nite, rnorm, conv

    enddo

    call cpu_time(tend)
    call toc(1,'solve')

    if ((myrank == 0).and.verbose) then
       rnpxg=real(grid(1)%npx,kind=rp)
       rnpyg=real(grid(1)%npy,kind=rp)
       rnxg =real(grid(1)%nx ,kind=rp)*rnpxg
       rnyg =real(grid(1)%ny ,kind=rp)*rnpyg
       rnzg =real(grid(1)%nz ,kind=rp)
       ! the rescaled time should be expressed in terms of error reduction,
       ! therefore the ratio rnorm/rnorm0
       perf = (tend-tstart)*(rnpxg*rnpyg)/(-log(rnorm/rnorm0)/log(10._rp))/(rnxg*rnyg*rnzg)
       write(*,*)'     --- summary ---'
       write(*,'(A,F8.3,A)')"     time spent to solve :",tend-tstart," s"
       write(*,'(A,E10.3)') "     rescaled performance:",perf
       write(*,*)'     ---------------'
    end if

10  format("     ite = ",I2,": res = ",E10.3," / conv = ",F10.3)

  end subroutine solve_p

  !---------------------------------------------------------------------
  subroutine Fcycle()

    integer(kind=ip):: lev,maxlev

    call tic(1,'Fcycle')

    maxlev=nlevs

    do lev=1,maxlev-1
       call fine2coarse(lev)
       grid(lev+1)%r=grid(lev+1)%b
    enddo

    call relax(maxlev, ns_coarsest)

    do lev=maxlev-1,1,-1
       call coarse2fine(lev)
       call Vcycle(lev)
    enddo

    call toc(1,'Fcycle')

  end subroutine Fcycle

  !----------------------------------------
  subroutine Vcycle(lev1)

    integer(kind=ip), intent(in) :: lev1
    integer(kind=ip)             :: lev
    real(kind=rp)                :: rnorm

    do lev=lev1,nlevs-1
       call relax(lev,ns_pre)
       call compute_residual(lev,rnorm)
       call fine2coarse(lev)
    enddo

    call relax(nlevs,ns_coarsest)

    do lev=nlevs-1,lev1,-1
       call coarse2fine(lev)
       call relax(lev,ns_post)
    enddo

  end subroutine Vcycle

end module mg_solvers
