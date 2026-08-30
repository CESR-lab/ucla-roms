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

  ! work arrays for the CG-accelerated solver (fine level only)
  real(kind=rp), dimension(:,:,:), allocatable :: sol  ! solution accumulator
  real(kind=rp), dimension(:,:,:), allocatable :: rcg  ! CG residual
  real(kind=rp), dimension(:,:,:), allocatable :: dcg  ! search direction
  real(kind=rp), dimension(:,:,:), allocatable :: qcg  ! A * dcg

contains

  !---------------------------------------------------------------------
  subroutine solve_p()

    if (solver_cg .and. (grid(1)%nz > 1)) then
       call solve_p_cg()
    else
       call solve_p_mg()
    endif

  end subroutine solve_p

  !---------------------------------------------------------------------
  subroutine solve_p_mg()

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

  end subroutine solve_p_mg

  !---------------------------------------------------------------------
  subroutine solve_p_cg()

    ! Flexible preconditioned conjugate gradients (FCG(1), Notay-type
    ! beta), preconditioned by one F-cycle. The flexible beta
    ! (A-orthogonalization against the previous direction) tolerates the
    ! nonsymmetric V(ns_pre,ns_post) red-black cycle, so the cycle needs
    ! no symmetrization. The operator is applied through matvec_lev1,
    ! the exact stencil of compute_residual_3D_8.
    !
    ! grid(1)%p/b/r are used as scratch by the preconditioner
    ! (b is set fresh by the caller before every solve and only read
    ! afterwards by debug output); the CG state lives in the module
    ! arrays sol/rcg/dcg/qcg. On exit grid(1)%p holds the solution.
    !
    ! Halo discipline: relax ends with fill_halo(p), so z=grid(1)%p and
    ! hence dcg (full-array updates) always carry valid halos; rcg/qcg
    ! halos are never read. No halo exchange is needed in this loop.

    integer(kind=ip) :: nite
    integer(kind=ip) :: nx,ny,nz

    real(kind=rp)    :: rnorm,bnorm,res0,conv,rnorm0
    real(kind=rp)    :: rho,zq,dq,alpha,beta,sumloc
    real(kind=rp), dimension(2) :: sum2
    real(kind=rp), dimension(:,:,:), pointer :: p,b,r

    real(kind=lg) :: tstart,tend,perf
    real(kind=rp) :: rnxg,rnyg,rnzg
    real(kind=rp) :: rnpxg,rnpyg,fac

    logical :: have_dir

    integer, save :: count = 1
    logical :: verbose

    verbose = .false.
    if (mod(count,output_freq)==0) verbose=.true.
    if (autotune)  verbose=.true.
    count = count+1
    if ((myrank==0).and.verbose) write(*,*)'     ------ CG -----'

    p  => grid(1)%p
    b  => grid(1)%b
    r  => grid(1)%r

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    if (.not.allocated(sol)) then
       allocate(sol(nz,0:ny+1,0:nx+1)); sol = zero
       allocate(rcg(nz,0:ny+1,0:nx+1)); rcg = zero
       allocate(dcg(nz,0:ny+1,0:nx+1)); dcg = zero
       allocate(qcg(nz,0:ny+1,0:nx+1)); qcg = zero
    endif

    call tic(1,'solve')
    call cpu_time(tstart)

    nite=0

    res0 = sum(b(1:nz,1:ny,1:nx)**2)
    call global_sum(1,res0,bnorm)
    bnorm = sqrt(bnorm)
    fac = sqrt(1./(nx*ny*nz))

    ! initial residual of the warm-started p:  rcg = b - A p
    call compute_residual(1,rnorm)
    rcg = r
    sol = p
    res0   = rnorm/bnorm
    rnorm0 = res0

    have_dir = .false.
    dq = one

    do while ((nite < solver_maxiter).and.(res0 > solver_prec).and.(res0*bnorm*fac>1e-15) )

       ! z = M^{-1} rcg : one F-cycle on (p=0, b=r=rcg); z is grid(1)%p
       p = zero
       b = rcg
       r = rcg
       call Fcycle()

       ! rho = <z,r>, zq = <z, A d_old> in one reduction
       sum2(1) = sum( p(1:nz,1:ny,1:nx)*rcg(1:nz,1:ny,1:nx) )
       if (have_dir) then
          sum2(2) = sum( p(1:nz,1:ny,1:nx)*qcg(1:nz,1:ny,1:nx) )
       else
          sum2(2) = zero
       endif
       call global_sum2(sum2)
       rho = sum2(1)
       zq  = sum2(2)

       if (have_dir) then
          beta = - zq/dq
          dcg = p + beta*dcg
       else
          dcg = p
          have_dir = .true.
       endif

       call matvec_lev1(dcg,qcg,sumloc)
       call global_sum(1,sumloc,dq)

       alpha = rho/dq
       sol = sol + alpha*dcg
       rcg(:,1:ny,1:nx) = rcg(:,1:ny,1:nx) - alpha*qcg(:,1:ny,1:nx)

       sumloc = sum(rcg(1:nz,1:ny,1:nx)**2)
       call global_sum(1,sumloc,rnorm)
       rnorm = sqrt(rnorm)/bnorm
       conv = res0/rnorm ! error reduction after this iteration
       res0 = rnorm

       nite = nite+1
       if ((myrank == 0).and.verbose) write(*,10) nite, rnorm, conv

    enddo

    p = sol

    call cpu_time(tend)
    call toc(1,'solve')

    if ((myrank == 0).and.verbose) then
       rnpxg=real(grid(1)%npx,kind=rp)
       rnpyg=real(grid(1)%npy,kind=rp)
       rnxg =real(grid(1)%nx ,kind=rp)*rnpxg
       rnyg =real(grid(1)%ny ,kind=rp)*rnpyg
       rnzg =real(grid(1)%nz ,kind=rp)
       perf = (tend-tstart)*(rnpxg*rnpyg)/(-log(rnorm/rnorm0)/log(10._rp))/(rnxg*rnyg*rnzg)
       write(*,*)'     --- summary ---'
       write(*,'(A,F8.3,A)')"     time spent to solve :",tend-tstart," s"
       write(*,'(A,E10.3)') "     rescaled performance:",perf
       write(*,*)'     ---------------'
    end if

10  format("     ite = ",I2,": res = ",E10.3," / conv = ",F10.3)

  end subroutine solve_p_cg

  !---------------------------------------------------------------------
  subroutine matvec_lev1(d3,q3,dq)

    ! q3 = A d3 on the fine level, plus dq = local <d3,q3>.
    ! Same stencil as compute_residual_3D_8 with b dropped and the sign
    ! flipped; keep the two in sync if the operator ever changes.
    ! d3 must have a valid halo; q3's halo is left untouched.

    real(kind=rp), dimension(:,:,:), allocatable, intent(in)    :: d3
    real(kind=rp), dimension(:,:,:), allocatable, intent(inout) :: q3
    real(kind=rp)                               , intent(out)   :: dq

    real(kind=rp), dimension(:,:,:,:), pointer :: cA
    integer(kind=ip) :: nx,ny,nz
    integer(kind=ip) :: i,j,k

    ! allocatable dummies preserve the actual bounds (1:nz,0:ny+1,0:nx+1),
    ! so natural indexing applies, exactly as in compute_residual_3D_8
    real(kind=rp) :: z

    cA => grid(1)%cA
    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    call tic(1,'matvec')

    dq = zero

    do i = 1,nx
       do j = 1,ny

          k=1 !lower level
          z =  cA(1,k,j,i)*d3(k,j,i)                                     &
             + cA(2,k+1,j,i)*d3(k+1,j,i)                                 &
             + cA(3,k,j,i)*d3(k+1,j-1,i)                                 &
             + cA(4,k,j,i)*d3(k  ,j-1,i) + cA(4,k  ,j+1,i)*d3(k  ,j+1,i) &
             + cA(5,k+1,j+1,i)*d3(k+1,j+1,i)                             &
             + cA(6,k,j,i)*d3(k+1,j,i-1)                                 &
             + cA(7,k,j,i)*d3(k  ,j,i-1) + cA(7,k  ,j,i+1)*d3(k  ,j,i+1) &
             + cA(8,k+1,j,i+1)*d3(k+1,j,i+1)                             &
!!  Special cross terms
             + cA(5,k,j,i)*d3(k,j+1,i-1) + cA(5,k,j-1,i+1)*d3(k,j-1,i+1) &
             + cA(8,k,j,i)*d3(k,j-1,i-1) + cA(8,k,j+1,i+1)*d3(k,j+1,i+1)
          q3(k,j,i) = z
          dq = dq + z*d3(k,j,i)

          do k = 2,nz-1 !interior levels
             z =  cA(1,k,j,i)*d3(k,j,i)                                        &
                + cA(2,k,j,i)*d3(k-1,j,i)   + cA(2,k+1,j,i)*d3(k+1,j,i)        &
                + cA(3,k,j,i)*d3(k+1,j-1,i) + cA(3,k-1,j+1,i)*d3(k-1,j+1,i)    &
                + cA(4,k,j,i)*d3(k  ,j-1,i) + cA(4,k  ,j+1,i)*d3(k  ,j+1,i)    &
                + cA(5,k,j,i)*d3(k-1,j-1,i) + cA(5,k+1,j+1,i)*d3(k+1,j+1,i)    &
                + cA(6,k,j,i)*d3(k+1,j,i-1) + cA(6,k-1,j,i+1)*d3(k-1,j,i+1)    &
                + cA(7,k,j,i)*d3(k  ,j,i-1) + cA(7,k  ,j,i+1)*d3(k  ,j,i+1)    &
                + cA(8,k,j,i)*d3(k-1,j,i-1) + cA(8,k+1,j,i+1)*d3(k+1,j,i+1)
             q3(k,j,i) = z
             dq = dq + z*d3(k,j,i)
          enddo

          k=nz !upper level
          z =  cA(1,k,j,i)*d3(k,j,i)                                     &
             + cA(2,k,j,i)*d3(k-1,j,i)                                   &
             + cA(3,k-1,j+1,i)*d3(k-1,j+1,i)                             &
             + cA(4,k,j,i)*d3(k  ,j-1,i) + cA(4,k  ,j+1,i)*d3(k  ,j+1,i) &
             + cA(5,k,j,i)*d3(k-1,j-1,i)                                 &
             + cA(6,k-1,j,i+1)*d3(k-1,j,i+1)                             &
             + cA(7,k,j,i)*d3(k  ,j,i-1) + cA(7,k  ,j,i+1)*d3(k  ,j,i+1) &
             + cA(8,k,j,i)*d3(k-1,j,i-1)
          q3(k,j,i) = z
          dq = dq + z*d3(k,j,i)

       enddo
    enddo

    call toc(1,'matvec')

  end subroutine matvec_lev1

  !---------------------------------------------------------------------
  subroutine global_sum2(sum2)

    ! fused 2-component global sum (fine level: no gather rescaling)

    real(kind=rp), dimension(2), intent(inout) :: sum2

    real(kind=rp), dimension(2) :: sumglo
    integer(kind=ip) :: ierr

    call MPI_ALLREDUCE(sum2,sumglo,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

    sum2 = sumglo

  end subroutine global_sum2

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
