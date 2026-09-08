module mg_projection

  use mg_cst
  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_mpi_exchange
  use mg_gather
  use mg_netcdf_out

  implicit none

  ! Face masks on the fine grid: 1 = fluid face, 0 = no-flux face.
  ! Set by set_face_masks (domain perimeter, non-periodic sides); an
  ! interior land mask can be folded in later through the same arrays.
  ! umsk(j,i): u-face i (west face of cell i);  vmsk(j,i): v-face j.
  real(kind=rp), dimension(:,:), allocatable :: umsk, vmsk

contains
  !-----------------------------------------------------------------------------------
  subroutine set_face_masks()
    integer(kind=ip) :: nx, ny
    nx = grid(1)%nx
    ny = grid(1)%ny
    if (.not.allocated(umsk)) allocate(umsk(0:ny+1,0:nx+1))
    if (.not.allocated(vmsk)) allocate(vmsk(0:ny+1,0:nx+1))
    umsk = one
    vmsk = one
    if (grid(1)%neighb(4) == MPI_PROC_NULL) umsk(:,1   ) = zero   ! west wall
    if (grid(1)%neighb(2) == MPI_PROC_NULL) umsk(:,nx+1) = zero   ! east wall
    if (grid(1)%neighb(1) == MPI_PROC_NULL) vmsk(1   ,:) = zero   ! south wall
    if (grid(1)%neighb(3) == MPI_PROC_NULL) vmsk(ny+1,:) = zero   ! north wall
  end subroutine set_face_masks

  !-----------------------------------------------------------------------------------
  subroutine set_matrices(first)

    logical, optional, intent(in) :: first   ! stage prints when .true.

    ! Define matrix coefficients cA
    ! Coefficients are stored in order of diagonals
    ! cA(1,:,:,:)      -> p(k,j,i)
    ! cA(2,:,:,:)      -> p(k-1,j,i)
    ! cA(3,:,:,:)      -> p(k+1,j-1,i)
    ! cA(4,:,:,:)      -> p(k,j-1,i)
    ! cA(5,:,:,:)      -> p(k-1,j-1,i)
    ! cA(6,:,:,:)      -> p(k+1,j,i-1)
    ! cA(7,:,:,:)      -> p(k,j,i-1)
    ! cA(8,:,:,:)      -> p(k-1,j,i-1)

    integer(kind=ip) :: lev
    integer(kind=ip) :: k,j,i
    integer(kind=ip) :: nx,ny,nz

    real(kind=rp), dimension(:,:),     pointer :: dx,dy
    real(kind=rp), dimension(:,:),     pointer :: dxu,dyv
    real(kind=rp), dimension(:,:,:),   pointer :: dzw
    real(kind=rp), dimension(:,:,:),   pointer :: Arx,Ary
    real(kind=rp), dimension(:,:)  ,   pointer :: Arz
    real(kind=rp), dimension(:,:,:),   pointer :: zxdy,zydx
    real(kind=rp), dimension(:,:,:),   pointer :: alpha
    real(kind=rp), dimension(:,:),     pointer :: beta
    real(kind=rp), dimension(:,:),     pointer :: gamu
    real(kind=rp), dimension(:,:),     pointer :: gamv
    real(kind=rp), dimension(:,:,:,:), pointer :: cA

    integer(kind=ip) :: dirichlet_flag 
    real(kind=rp)    :: sig, fj, fjm, fi, fim
    

    if (surface_neumann) then
       dirichlet_flag = 0
    else
       dirichlet_flag = 1
    endif

    do lev = 1, nlevs

       nx = grid(lev)%nx
       ny = grid(lev)%ny
       nz = grid(lev)%nz

       dx    => grid(lev)%dx    !
       dy    => grid(lev)%dy    !
       dxu   => grid(lev)%dxu   !
       dyv   => grid(lev)%dyv   !
       dzw   => grid(lev)%dzw   !
       Arx   => grid(lev)%Arx   !
       Ary   => grid(lev)%Ary   !
       Arz   => grid(lev)%Arz   !
       alpha => grid(lev)%alpha !
       beta  => grid(lev)%beta  !
       gamu  => grid(lev)%gamu  !
       gamv  => grid(lev)%gamv  !
       zxdy  => grid(lev)%zxdy  !
       zydx  => grid(lev)%zydx  !
       cA    => grid(lev)%cA    !

       !! interaction coeff with neighbours !!

       !---------------!
       !- lower level -!
       !---------------!
       k = 1
       do i = 1,nx
          do j = 1,ny+1
             ! couples with k+1,j-1
             cA(3,k,j,i) = qrt * ( zydx(k+1,j,i) + zydx(k,j-1,i) )
             ! couples with j-1
             cA(4,k,j,i) = hlf * (gamv(j,i) + gamv(j-1,i)) * Ary(k,j,i) / dyv(j,i) &
                  - qrt * ( zydx(k,j-1,i) - zydx(k,j,i) )
          enddo
       enddo

       do i = 1,nx+1
          do j = 1,ny
             ! couples with k+1,i-1
             cA(6,k,j,i) = qrt * ( zxdy(k+1,j,i) + zxdy(k,j,i-1) )
             ! couples with i-1
             cA(7,k,j,i) = hlf * (gamu(j,i) + gamu(j,i-1))  * Arx(k,j,i) / dxu(j,i) &
                  - qrt * ( zxdy(k,j,i-1) - zxdy(k,j,i) )
          enddo
       enddo

       do i = 1,nx+1
          do j = 0,ny
             ! only for k==1, couples with j+1,i-1
             cA(5,k,j,i) = beta(j,i-1) + beta(j+1,i)
          enddo
       enddo

       do i = 1,nx+1
          do j = 1,ny+1
             ! only for k==1, couples with j-1,i-1
             cA(8,k,j,i) = - beta(j,i-1) - beta(j-1,i)
          enddo
       enddo

       !-------------------!
       !- interior levels -!
       !-------------------!
       do i = 1,nx
          do j = 1,ny
             do k = 2,nz-1 
                ! couples with k-1
                cA(2,k,j,i) = &
                     ( Arz(j,i) / dzw(k,j,i) ) * &
                     hlf * (alpha(k,j,i)+alpha(k-1,j,i))
             enddo
          enddo
       enddo

       do i = 1,nx
          do j = 1,ny+1
             do k = 2,nz-1 
                ! couples with k+1,j-1
                cA(3,k,j,i) = qrt * ( zydx(k+1,j,i) + zydx(k,j-1,i) )
                ! couples with j-1
                cA(4,k,j,i) =  Ary(k,j,i) / dyv(j,i) 
                ! couples with k-1,j-1
                cA(5,k,j,i) = - qrt * ( zydx(k-1,j,i) + zydx(k,j-1,i) )
             enddo
          enddo
       enddo

       do i = 1,nx+1
          do j = 1,ny 
             do k = 2,nz-1 
                ! couples with k+1,i-1
                cA(6,k,j,i) = qrt * ( zxdy(k+1,j,i) + zxdy(k,j,i-1) )
                ! couples with i-1
                cA(7,k,j,i) = Arx(k,j,i) / dxu(j,i) 
                ! couples with k-1,i-1
                cA(8,k,j,i) = - qrt * ( zxdy(k-1,j,i) + zxdy(k,j,i-1) )
             enddo
          enddo
       enddo

       !---------------!
       !- upper level -!
       !---------------!
       k = nz
       do i = 1,nx    
          do j = 1,ny 
             ! couples with k-1
             cA(2,k,j,i) = (Arz(j,i)/dzw(k,j,i)) * hlf * ( alpha(k,j,i) + alpha(k-1,j,i) )
          enddo
       enddo

       do i = 1,nx
          do j = 1,ny+1
             ! couples with j-1
             fj = 1._rp ; fjm = 1._rp
             if (lev == 1) then
                fj = sfcf(j,i) ; fjm = sfcf(j-1,i)
             endif
             cA(4,k,j,i) = Ary(k,j,i) / dyv(j,i) &
                  + qrt * ( - zydx(k,j-1,i)*(2*sigtop(dzw(k+1,j-1,i))*fjm-1) &
                            + zydx(k,j  ,i)*(2*sigtop(dzw(k+1,j  ,i))*fj -1) ) 
             ! couples with k-1,j-1
             cA(5,k,j,i) = - qrt * ( zydx(k-1,j,i) + zydx(k,j-1,i) )
          enddo
       enddo

       do i = 1,nx+1
          do j = 1,ny 
             ! with Neumann BC, the CA(7,:,:) has flip signed on the slope term
             ! this will be in the paper
             ! couples with i-1
             fi = 1._rp ; fim = 1._rp
             if (lev == 1) then
                fi = sfcf(j,i) ; fim = sfcf(j,i-1)
             endif
             cA(7,k,j,i) = Arx(k,j,i) / dxu(j,i) &
                  + qrt * ( -zxdy(k,j,i-1)*(2*sigtop(dzw(k+1,j,i-1))*fim-1) &
                            +zxdy(k,j,i  )*(2*sigtop(dzw(k+1,j,i  ))*fi -1) ) 
             ! couples with k-1,i-1
             cA(8,k,j,i) = - qrt * ( zxdy(k-1,j,i) + zxdy(k,j,i-1) )
          enddo
       enddo

       call fill_halo(lev,cA)

       !! self-interaction coeff !!

       do i = 1,nx
          do j = 1,ny

             k = 1 !lower level
             cA(1,k,j,i) = &
                  - Arz(j,i) / dzw(k+1,j,i) * hlf * (alpha(k+1,j,i) + alpha(k,j,i)) &
                  - Arx(k,j,i  )/dxu(j,i  ) * hlf * (gamu(j,i) + gamu(j  ,i-1)) &
                  - Arx(k,j,i+1)/dxu(j,i+1) * hlf * (gamu(j,i) + gamu(j  ,i+1)) &
                  - Ary(k,j  ,i)/dyv(j  ,i) * hlf * (gamv(j,i) + gamv(j-1,i  )) &
                  - Ary(k,j+1,i)/dyv(j+1,i) * hlf * (gamv(j,i) + gamv(j+1,i  ))

             do k = 2,nz-1 !interior levels
                cA(1,k,j,i) = &
                     - Arz(j,i) / dzw(k+1,j,i) * hlf * (alpha(k+1,j,i) + alpha(k,j,i)) &
                     - Arz(j,i) / dzw(k  ,j,i) * hlf * (alpha(k-1,j,i) + alpha(k,j,i)) &
                     - Arx(k,j,i  )/dxu(j,i  )  &
                     - Arx(k,j,i+1)/dxu(j,i+1)  &
                     - Ary(k,j  ,i)/dyv(j  ,i)  &
                     - Ary(k,j+1,i)/dyv(j+1,i)
             enddo

             k=nz ! upper level
             fi = 1._rp
             if (lev == 1) fi = sfcf(j,i)
             cA(1,k,j,i) = &
                  - Arz(j,i) / dzw(k+1,j,i) * alpha(k,j,i) * sigtop(dzw(k+1,j,i)) * fi &
                  - Arz(j,i) / dzw(k  ,j,i) * hlf * (alpha(k-1,j,i) + alpha(k,j,i)) &
                  - Arx(k,j,i  )/dxu(j,i  )  &
                  - Arx(k,j,i+1)/dxu(j,i+1)  &
                  - Ary(k,j  ,i)/dyv(j  ,i)  &
                  - Ary(k,j+1,i)/dyv(j+1,i)

          enddo
       enddo

       if (netcdf_output) then
          if (myrank==0) write(*,*)'       write cA in a netcdf file'
          call write_netcdf(grid(lev)%cA,vname='ca',netcdf_file_name='cA.nc',rank=myrank,iter=lev)
       endif

       ! Fine grid: replace the hand-assembled stencil by the exact
       ! stencil of the MASKED operator  A = -D (M T M) G, so that the
       ! no-flux faces of the domain perimeter (and, later, of an
       ! interior land mask) are represented by zero rows/columns of T
       ! instead of a mirrored p.  Mirroring is a zero-flux condition
       ! only for the Arx*px part of the U row; with slope cross terms
       ! it makes the effective operator non-symmetric, which under the
       ! singular Neumann surface condition leaves an irreducible
       ! residual (LOCK, 2026-09-04).  The stencil is read off 27
       ! colouring vectors through correction_uvw itself, hence exact by
       ! construction; the interior is checked against the formulas.
       if (lev == 1) then
          call stage('face masks')
          call set_face_masks()
          call stage('stencil extraction (27 colours)')
          call assemble_masked_stencil()
       endif
       if (lev == 1 .and. nlevs > 1) call stage('coarse levels')

    enddo

  contains
    subroutine stage(what)
      character(len=*), intent(in) :: what
      if (.not.present(first)) return
      if (first .and. myrank==0) then
         write(*,'(A,A)') '  set_matrices: ', what
         flush(6)
      endif
    end subroutine stage

  end subroutine set_matrices

  !-------------------------------------------------------------------------
  subroutine assemble_masked_stencil()

    integer(kind=ip) :: nx,ny,nz,i,j,k,a,bb,c,di,dj,dk,m,ig,jg,pi,pj
    integer(kind=ip) :: ierr
    real(kind=rp), dimension(:,:,:,:), allocatable :: A27, cAf
    real(kind=rp), dimension(:,:,:), pointer :: p,du,dv,dw
    real(kind=rp), dimension(:,:,:,:), pointer :: cA
    real(kind=rp) :: dint, dext, dsym, dnul, amax, gl(4), lc(4)

    nx = grid(1)%nx;  ny = grid(1)%ny;  nz = grid(1)%nz
    p  => grid(1)%p;  du => grid(1)%du;  dv => grid(1)%dv;  dw => grid(1)%dw
    cA => grid(1)%cA
    pj = myrank/grid(1)%npx
    pi = mod(myrank,grid(1)%npx)

    allocate(A27(27,nz,ny,nx));  A27 = zero
    allocate(cAf(8,nz,0:ny+1,0:nx+1));  cAf = cA        ! formula stencil, kept for the check

    ! 27 colours on GLOBAL indices (consistent across seams)
    do a = 0,2
      do bb = 0,2
        do c = 0,2
          p = zero
          do i = 1,nx
            ig = i + pi*nx
            if (mod(ig,3) /= a) cycle
            do j = 1,ny
              jg = j + pj*ny
              if (mod(jg,3) /= bb) cycle
              do k = 1,nz
                if (mod(k,3) == c) p(k,j,i) = one
              enddo
            enddo
          enddo
          call fill_halo(1,p)
          dw(1,:,:) = zero
          call correction_uvw()                  ! (du,dv,dw) = M T M G p
          do i = 1,nx
            ig = i + pi*nx
            di = -2                                ! offset with colour a
            do m = -1,1
              if (mod(ig+m+3,3) == a) di = m
            enddo
            do j = 1,ny
              jg = j + pj*ny
              dj = -2
              do m = -1,1
                if (mod(jg+m+3,3) == bb) dj = m
              enddo
              do k = 1,nz
                dk = -2
                do m = -1,1
                  if (mod(k+m+3,3) == c) dk = m
                enddo
                m = 14 + di + 3*dj + 9*dk           ! 1..27, centre = 14
                A27(m,k,j,i) = -( du(k,j,i+1) - du(k,j,i)     &
                                + dv(k,j+1,i) - dv(k,j,i)     &
                                + dw(k+1,j,i) - dw(k,j,i) )
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    p = zero

    ! stencil -> the 8 stored diagonals (interior cells)
    do i = 1,nx
      do j = 1,ny
        do k = 1,nz
          cA(1,k,j,i) = A27(14      ,k,j,i)            ! ( 0, 0, 0)
          cA(2,k,j,i) = A27(14-9    ,k,j,i)            ! (-1, 0, 0)  k-1
          cA(3,k,j,i) = A27(14+9-3  ,k,j,i)            ! (+1,-1, 0)
          cA(4,k,j,i) = A27(14-3    ,k,j,i)            ! ( 0,-1, 0)
          cA(6,k,j,i) = A27(14+9-1  ,k,j,i)            ! (+1, 0,-1)
          cA(7,k,j,i) = A27(14-1    ,k,j,i)            ! ( 0, 0,-1)
          if (k >= 2) then
            cA(5,k,j,i) = A27(14-9-3,k,j,i)            ! (-1,-1, 0)
            cA(8,k,j,i) = A27(14-9-1,k,j,i)            ! (-1, 0,-1)
          else
            cA(5,k,j,i) = A27(14+3-1,k,j,i)            ! ( 0,+1,-1)  bottom xy
            cA(8,k,j,i) = A27(14-3-1,k,j,i)            ! ( 0,-1,-1)  bottom xy
          endif
        enddo
      enddo
    enddo
    ! halo-stored diagonals read by the interior rows (east, north, south)
    do j = 1,ny
      do k = 1,nz
        cA(7,k,j,nx+1) = A27(14+1,k,j,nx)                     ! ( 0, 0,+1)
        if (k >= 2) cA(6,k-1,j,nx+1) = A27(14-9+1,k,j,nx)     ! (-1, 0,+1)
        if (k <= nz-1) cA(8,k+1,j,nx+1) = A27(14+9+1,k,j,nx)  ! (+1, 0,+1)
      enddo
      cA(5,1,j-1,nx+1) = A27(14-3+1,1,j,nx)                   ! ( 0,-1,+1)
      cA(8,1,j+1,nx+1) = A27(14+3+1,1,j,nx)                   ! ( 0,+1,+1)
    enddo
    do i = 1,nx
      do k = 1,nz
        cA(4,k,ny+1,i) = A27(14+3,k,ny,i)                     ! ( 0,+1, 0)
        if (k >= 2) cA(3,k-1,ny+1,i) = A27(14-9+3,k,ny,i)     ! (-1,+1, 0)
        if (k <= nz-1) cA(5,k+1,ny+1,i) = A27(14+9+3,k,ny,i)  ! (+1,+1, 0)
      enddo
      cA(8,1,ny+1,i+1) = A27(14+3+1,1,ny,i)                   ! ( 0,+1,+1)
      cA(5,1,0,i+1)    = A27(14-3+1,1,1,i)                    ! ( 0,-1,+1)
    enddo

    ! checks: (1) interior agreement with the formulas away from the
    ! perimeter, (2) size of the coefficients outside the 15-point
    ! pattern, (3) symmetry A(c,n) = A(n,c) for in-rank pairs
    amax = maxval(abs(cA(1,1:nz,1:ny,1:nx)))
    dint = zero;  dext = zero;  dsym = zero;  dnul = zero
    do i = 1,nx
      do j = 1,ny
        do k = 1,nz
          if (i>=3 .and. i<=nx-2 .and. j>=3 .and. j<=ny-2) then
            dint = max(dint, maxval(abs(cA(1:8,k,j,i)-cAf(1:8,k,j,i))))
          else
            dext = max(dext, maxval(abs(cA(1:8,k,j,i)-cAf(1:8,k,j,i))))
          endif
          ! pattern: allowed offsets (dk,dj,di)
          do m = 1,27
            dk = (m-1)/9 - 1;  dj = mod((m-1)/3,3) - 1;  di = mod(m-1,3) - 1
            if (dj/=0 .and. di/=0 .and. .not.(k==1 .and. dk==0)) dnul = max(dnul,abs(A27(m,k,j,i)))
            if (dj/=0 .and. di/=0 .and. k==1 .and. dk/=0)       dnul = max(dnul,abs(A27(m,k,j,i)))
            if (dk/=0 .and. dj/=0 .and. di/=0)                  dnul = max(dnul,abs(A27(m,k,j,i)))
            ! symmetry with the reverse offset at the neighbour (in-rank)
            if (k+dk>=1 .and. k+dk<=nz .and. j+dj>=1 .and. j+dj<=ny .and. i+di>=1 .and. i+di<=nx) &
              dsym = max(dsym, abs(A27(m,k,j,i) - A27(28-m,k+dk,j+dj,i+di)))
          enddo
        enddo
      enddo
    enddo
    lc = (/ dint, dext, dsym, dnul /)
    call MPI_Allreduce(lc,gl,4,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
    lc(1) = amax
    call MPI_Allreduce(lc(1),amax,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
    if (myrank==0) then
       write(*,'(A)') '     masked stencil (fine grid):'
       write(*,'(A,ES9.2,A,ES9.2,A)') '       max|num-formula| interior ', gl(1)/amax, &
            '   near perimeter ', gl(2)/amax, '   (rel. to max|diag|)'
       write(*,'(A,ES9.2,A,ES9.2)') '       symmetry defect ', gl(3)/amax, &
            '   outside 15-point pattern ', gl(4)/amax
    endif
    if (gl(1)/amax > 1e-10_rp) then
       if (robin_beta > 0._rp) then
          ! Robin surface: the hand formulas carry a single-column surface
          ! factor in the slope terms, the extracted stencil the exact
          ! two-column one; the fine grid uses the extracted stencil.
          if (myrank==0) write(*,*) 'assemble_masked_stencil: interior mismatch expected with Robin surface (slope terms)'
       else
          if (myrank==0) write(*,*) 'assemble_masked_stencil: interior mismatch -- stopping'
          call MPI_Abort(MPI_COMM_WORLD,1,ierr)
       endif
    endif
    deallocate(A27,cAf)

  end subroutine assemble_masked_stencil

  !-------------------------------------------------------------------------     
  subroutine correction_uvw()

    !! u,v,w are fluxes, the correction is T*grad(p)

    integer(kind=ip):: k, j, i
    integer(kind=ip):: nx, ny, nz

    real(kind=rp) :: gamma
    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
    real(kind=rp), dimension(:,:)  , pointer :: dxu,dyv
    real(kind=rp), dimension(:,:)  , pointer :: Arz
    real(kind=rp), dimension(:,:,:), pointer :: dzw
    real(kind=rp), dimension(:,:,:), pointer :: Arx,Ary
    real(kind=rp), dimension(:,:,:), pointer :: zxdy,zydx
    real(kind=rp), dimension(:,:,:), pointer :: alpha
    real(kind=rp), dimension(:,:)  , pointer :: beta
    real(kind=rp), dimension(:,:,:), pointer :: p
    real(kind=rp), dimension(:,:,:), pointer :: px,py,pz
    real(kind=rp), dimension(:,:,:), pointer :: du,dv,dw

    integer(kind=ip) :: dirichlet_flag
 
!    if (myrank==0) write(*,*)'   - compute pressure gradient and translate to fluxes'

    if (surface_neumann) then
       dirichlet_flag = 0
    else
       dirichlet_flag = 1
    endif

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    dx    => grid(1)%dx
    dy    => grid(1)%dy
    dxu   => grid(1)%dxu
    dyv   => grid(1)%dyv
    dzw   => grid(1)%dzw
    Arx   => grid(1)%Arx
    Ary   => grid(1)%Ary
    Arz   => grid(1)%Arz
    alpha => grid(1)%alpha
    beta  => grid(1)%beta
    zxdy  => grid(1)%zxdy
    zydx  => grid(1)%zydx
    p     => grid(1)%p

    !! Pressure gradient -

    px => grid(1)%px
    py => grid(1)%py
    pz => grid(1)%pz

    do i = 1,nx+1
        do j = 0,ny+1
          do k = 1,nz
             px(k,j,i) = -one / dxu(j,i) * (p(k,j,i)-p(k,j,i-1))
          enddo
       enddo
    enddo

    do i = 0,nx+1
       do j = 1,ny+1 
          do k = 1,nz
             py(k,j,i) = -one / dyv(j,i) * (p(k,j,i)-p(k,j-1,i))
          enddo
       enddo
    enddo

    ! column masking: no pressure gradient through no-flux faces
    if (allocated(umsk)) then
       do i = 1,nx+1
          do j = 0,ny+1
             px(:,j,i) = px(:,j,i)*umsk(j,i)
          enddo
       enddo
       do i = 0,nx+1
          do j = 1,ny+1
             py(:,j,i) = py(:,j,i)*vmsk(j,i)
          enddo
       enddo
    endif

    do i = 0,nx+1
       do j = 0,ny+1

          k = 1 !bottom pressure gradient is undefined
          pz(k,j,i) = 9999999999.

          do k = 2,nz !interior levels
             pz(k,j,i) = -one / dzw(k,j,i) * (p(k,j,i)-p(k-1,j,i))
          enddo

          k = nz+1 !surface: Dirichlet, Neumann or Robin (q + beta dq/dz = f); 0 where prescribed
          pz(k,j,i) =  -one / dzw(k,j,i) * (-p(k-1,j,i)) * sigtop(dzw(k,j,i)) * sfcf(j,i)

       enddo
    enddo

    !! Correction for U -

    du => grid(1)%du

    do i = 1,nx+1  
       do j = 1,ny 
          k = 1
          gamma = one - qrt * ( &
               (zxdy(k,j,i  )/dy(j,i  ))**2/alpha(k,j,i  ) + &
               (zxdy(k,j,i-1)/dy(j,i-1))**2/alpha(k,j,i-1) )
          du(k,j,i) = gamma * Arx(k,j,i) * px(k,j,i) &
               - qrt * ( &
               + zxdy(k,j,i  ) * dzw(k+1,j,i  ) * pz(k+1,j,i  ) &
               + zxdy(k,j,i-1) * dzw(k+1,j,i-1) * pz(k+1,j,i-1) )  &
               - beta(j,i-1)   * dyv(j  ,i-1)   * py(k,j  ,i-1) &
               - beta(j,i-1)   * dyv(j+1,i-1)   * py(k,j+1,i-1) &
               - beta(j,i  )   * dyv(j  ,i  )   * py(k,j  ,i  ) &
               - beta(j,i  )   * dyv(j+1,i  )   * py(k,j+1,i  )

          do k = 2,nz-1 
             du(k,j,i) = Arx(k,j,i) * px(k,j,i) &
                  - qrt * ( &
                  + zxdy(k,j,i  ) * dzw(k  ,j,i  ) * pz(k  ,j,i  ) &
                  + zxdy(k,j,i  ) * dzw(k+1,j,i  ) * pz(k+1,j,i  ) &
                  + zxdy(k,j,i-1) * dzw(k  ,j,i-1) * pz(k  ,j,i-1) &
                  + zxdy(k,j,i-1) * dzw(k+1,j,i-1) * pz(k+1,j,i-1) )
          enddo

          k = nz
          du(k,j,i) = Arx(k,j,i) * px(k,j,i) &
               - qrt * ( &
               + zxdy(k,j,i  ) *       dzw(k  ,j,i  ) * pz(k  ,j,i  ) &
               + zxdy(k,j,i  ) * two * dzw(k+1,j,i  ) * pz(k+1,j,i  ) &
               + zxdy(k,j,i-1) *       dzw(k  ,j,i-1) * pz(k  ,j,i-1) &
               + zxdy(k,j,i-1) * two * dzw(k+1,j,i-1) * pz(k+1,j,i-1) )
       enddo
    enddo

    !! Correction for V - 

    dv => grid(1)%dv

    do i = 1,nx
       do j = 1,ny+1
          k = 1
          gamma = one - qrt * (  &
               (zydx(k,j  ,i)/dx(j  ,i))**2/alpha(k,j  ,i  ) + &
               (zydx(k,j-1,i)/dx(j-1,i))**2/alpha(k,j-1,i) )
          dv(k,j,i) = gamma * Ary(k,j,i) * py(k,j,i) &
               - qrt * ( &
               + zydx(k,j  ,i) * dzw(k+1,j  ,i) * pz(k+1,j  ,i) &
               + zydx(k,j-1,i) * dzw(k+1,j-1,i) * pz(k+1,j-1,i) ) &
               - beta(j-1,i)   * dxu(j-1,i  )   * px(k,j-1,i  ) &
               - beta(j-1,i)   * dxu(j-1,i+1)   * px(k,j-1,i+1) &
               - beta(j  ,i)   * dxu(j  ,i  )   * px(k,j  ,i  ) &
               - beta(j  ,i)   * dxu(j  ,i+1)   * px(k,j  ,i+1)

          do k = 2,nz-1
             dv(k,j,i) =  Ary(k,j,i) * py(k,j,i) &
                  - qrt * ( &
                  + zydx(k,j  ,i) * dzw(k  ,j  ,i) * pz(k  ,j  ,i) &
                  + zydx(k,j  ,i) * dzw(k+1,j  ,i) * pz(k+1,j  ,i) &
                  + zydx(k,j-1,i) * dzw(k  ,j-1,i) * pz(k  ,j-1,i) &
                  + zydx(k,j-1,i) * dzw(k+1,j-1,i) * pz(k+1,j-1,i) )
          enddo

          k = nz
          dv(k,j,i) = Ary(k,j,i) * py(k,j,i) &
               - qrt * ( &
               + zydx(k,j  ,i)       * dzw(k  ,j  ,i) * pz(k  ,j  ,i) &
               + zydx(k,j  ,i) * two * dzw(k+1,j  ,i) * pz(k+1,j  ,i) &
               + zydx(k,j-1,i)       * dzw(k  ,j-1,i) * pz(k  ,j-1,i) &
               + zydx(k,j-1,i) * two * dzw(k+1,j-1,i) * pz(k+1,j-1,i) ) 

       enddo
    enddo

    !! Correction for W -

    dw => grid(1)%dw

    do i = 1,nx
       do j = 1,ny

          do k = 2,nz
             dw(k,j,i) =  hlf * (alpha(k-1,j,i) + alpha(k,j,i)) * Arz(j,i) * pz(k,j,i) &
                  - qrt * ( &
                  + zxdy(k  ,j,i) * dxu(j,i  ) * px(k  ,j,i  ) &
                  + zxdy(k  ,j,i) * dxu(j,i+1) * px(k  ,j,i+1) &
                  + zxdy(k-1,j,i) * dxu(j,i  ) * px(k-1,j,i  ) &
                  + zxdy(k-1,j,i) * dxu(j,i+1) * px(k-1,j,i+1) ) &
                  - qrt * ( &
                  + zydx(k  ,j,i) * dyv(j  ,i) * py(k  ,j  ,i) &
                  + zydx(k  ,j,i) * dyv(j+1,i) * py(k  ,j+1,i) &
                  + zydx(k-1,j,i) * dyv(j  ,i) * py(k-1,j  ,i) &
                  + zydx(k-1,j,i) * dyv(j+1,i) * py(k-1,j+1,i) )
          enddo

          k = nz+1 
          ! surface: vertical part and tilt terms both carry the surface
          ! factor (1 Dirichlet, 0 Neumann, dzw/(dzw+beta) Robin) so that the
          ! operator stays symmetric with the two*dzw*pz terms of the u,v rows
          dw(k,j,i) = alpha(k-1,j,i) * Arz(j,i) * pz(k,j,i) &
               + sigtop(dzw(k,j,i)) * sfcf(j,i) * ( &
               - hlf * ( &
               + zxdy(k-1,j,i) * dxu(j,i  ) * px(k-1,j,i  ) &
               + zxdy(k-1,j,i) * dxu(j,i+1) * px(k-1,j,i+1) ) &
               - hlf * ( &
               + zydx(k-1,j,i) * dyv(j  ,i) * py(k-1,j  ,i) &
               + zydx(k-1,j,i) * dyv(j+1,i) * py(k-1,j+1,i) ) )
       enddo
    enddo

    ! row masking: no flux correction through no-flux faces
    if (allocated(umsk)) then
       do i = 1,nx+1
          do j = 1,ny
             du(:,j,i) = du(:,j,i)*umsk(j,i)
          enddo
       enddo
       do i = 1,nx
          do j = 1,ny+1
             dv(:,j,i) = dv(:,j,i)*vmsk(j,i)
          enddo
       enddo
    endif

  end subroutine correction_uvw

end module mg_projection

