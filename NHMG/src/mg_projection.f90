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

    logical, optional, intent(in) :: first   ! full extraction + stage prints
    logical :: full

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
          full = .true.
          if (present(first)) full = first
          if (full) then
             call stage('face masks')
             call set_face_masks()
             call stage('stencil extraction (27 colours)')
             call assemble_masked_stencil()
          else
             ! only the wall-adjacent rows differ from the formulas
             call refresh_wall_rows()
          endif
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
    ! Read the fine-grid operator M T M G off 27 colouring vectors and
    ! store it in cA (15-point symmetric storage: 8 slots per cell plus
    ! the halo slots read by the east, north and south rows).  Each cell
    ! meets each offset exactly once over the 27 colours, so every value
    ! goes straight to its slot: no extra 3D storage (a 27+8 field
    ! scratch was 0.7 GB per rank on a 96x96x256 subdomain, 2026-09-08).
    ! Checks: agreement with the formulas (interior cells) and the size
    ! of entries outside the pattern are taken as running maxima; the
    ! symmetry of the true operator is tested globally with
    ! <y,A x> = <x,A y> on two deterministic vectors.

    integer(kind=ip) :: nx,ny,nz,i,j,k,a,bb,c,di,dj,dk,m,n,ig,jg,pi,pj
    integer(kind=ip) :: ierr
    real(kind=rp), dimension(:,:,:), pointer :: p,du,dv,dw,x,y
    real(kind=rp), dimension(:,:,:,:), pointer :: cA
    real(kind=rp) :: dint, dext, dnul, amax, v, gl(4), lc(4), gs(2), ls(2), sxy, syx

    nx = grid(1)%nx;  ny = grid(1)%ny;  nz = grid(1)%nz
    p  => grid(1)%p;  du => grid(1)%du;  dv => grid(1)%dv;  dw => grid(1)%dw
    cA => grid(1)%cA
    x  => grid(1)%b;  y  => grid(1)%r        ! free work arrays at init
    pj = myrank/grid(1)%npx
    pi = mod(myrank,grid(1)%npx)

    dint = zero;  dext = zero;  dnul = zero

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
                v = -( du(k,j,i+1) - du(k,j,i)     &
                     + dv(k,j+1,i) - dv(k,j,i)     &
                     + dw(k+1,j,i) - dw(k,j,i) )
                n = stencil_slot(dk,dj,di,k)
                if (n > 0) then
                  ! the formula value still in the slot is the reference
                  if (i>=3 .and. i<=nx-2 .and. j>=3 .and. j<=ny-2) then
                    dint = max(dint, abs(v - cA(n,k,j,i)))
                  else
                    dext = max(dext, abs(v - cA(n,k,j,i)))
                  endif
                elseif (dj/=0 .and. di/=0 .and. .not.(k==1 .and. dk==0)) then
                  ! outside the 15-point pattern (xy corners except the
                  ! bottom level's own pair, and all xyz corners)
                  dnul = max(dnul,abs(v))
                endif
                call put_entry(cA,nx,ny,nz,v,dk,dj,di,k,j,i)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo

    ! symmetry of the true operator: <y,Ax> - <x,Ay> on two vectors
    do i = 1,nx
      ig = i + pi*nx
      do j = 1,ny
        jg = j + pj*ny
        do k = 1,nz
          x(k,j,i) = sin(0.731_rp*ig + 1.213_rp*jg + 2.117_rp*k)
          y(k,j,i) = cos(1.379_rp*ig + 0.577_rp*jg + 1.911_rp*k)
        enddo
      enddo
    enddo
    do i = 0,nx+1                       ! loops: pointer targets may alias,
      do j = 0,ny+1                     ! a whole-array copy would build a
        do k = 1,nz                     ! full-field temporary on the stack
          p(k,j,i) = x(k,j,i)
        enddo
      enddo
    enddo
    call fill_halo(1,p);  dw(1,:,:) = zero;  call correction_uvw()
    syx = zero
    do i = 1,nx
      do j = 1,ny
        do k = 1,nz
          syx = syx - y(k,j,i) * ( du(k,j,i+1) - du(k,j,i)  &
                                 + dv(k,j+1,i) - dv(k,j,i)  &
                                 + dw(k+1,j,i) - dw(k,j,i) )
        enddo
      enddo
    enddo
    do i = 0,nx+1
      do j = 0,ny+1
        do k = 1,nz
          p(k,j,i) = y(k,j,i)
        enddo
      enddo
    enddo
    call fill_halo(1,p);  dw(1,:,:) = zero;  call correction_uvw()
    sxy = zero
    do i = 1,nx
      do j = 1,ny
        do k = 1,nz
          sxy = sxy - x(k,j,i) * ( du(k,j,i+1) - du(k,j,i)  &
                                 + dv(k,j+1,i) - dv(k,j,i)  &
                                 + dw(k+1,j,i) - dw(k,j,i) )
        enddo
      enddo
    enddo
    p = zero;  x = zero;  y = zero

    amax = zero
    do i = 1,nx
      do j = 1,ny
        do k = 1,nz
          amax = max(amax, abs(cA(1,k,j,i)))
        enddo
      enddo
    enddo
    lc = (/ dint, dext, dnul, amax /)
    call MPI_Allreduce(lc,gl,4,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
    amax = gl(4)
    ls = (/ syx, sxy /)
    call MPI_Allreduce(ls,gs,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    syx = gs(1);  sxy = gs(2)
    if (myrank==0) then
       write(*,'(A)') '     masked stencil (fine grid):'
       write(*,'(A,ES9.2,A,ES9.2,A)') '       max|num-formula| interior ', gl(1)/amax, &
            '   near perimeter ', gl(2)/amax, '   (rel. to max|diag|)'
       write(*,'(A,ES9.2,A,ES9.2)') '       symmetry |<y,Ax>-<x,Ay>|/|<y,Ax>| ', &
            abs(syx-sxy)/max(abs(syx),tiny(one)), '   outside 15-point pattern ', gl(3)/amax
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

  end subroutine assemble_masked_stencil

  !-------------------------------------------------------------------------
  function stencil_slot(dk,dj,di,k) result(n)
    ! slot (1..8) of the 15-point storage that holds, at a cell on level
    ! k, the coupling to the neighbour at offset (dk,dj,di); 0 if that
    ! coupling is stored at the neighbour instead (or is off-pattern)
    integer(kind=ip), intent(in) :: dk,dj,di,k
    integer(kind=ip) :: n
    n = 0
    if (dk==0 .and. dj==0 .and. di==0) then;        n = 1
    elseif (dk==-1 .and. dj==0 .and. di==0) then;   n = 2
    elseif (dk==1 .and. dj==-1 .and. di==0) then;   n = 3
    elseif (dk==0 .and. dj==-1 .and. di==0) then;   n = 4
    elseif (dk==1 .and. dj==0 .and. di==-1) then;   n = 6
    elseif (dk==0 .and. dj==0 .and. di==-1) then;   n = 7
    elseif (k>=2) then
      if (dk==-1 .and. dj==-1 .and. di==0) n = 5
      if (dk==-1 .and. dj==0 .and. di==-1) n = 8
    else                                              ! bottom level: xy pair
      if (dk==0 .and. dj==1 .and. di==-1) n = 5
      if (dk==0 .and. dj==-1 .and. di==-1) n = 8
    endif
  end function stencil_slot

  !-------------------------------------------------------------------------
  subroutine put_entry(cA,nx,ny,nz,v,dk,dj,di,k,j,i)
    ! store the operator value v of row (k,j,i), offset (dk,dj,di):
    ! its own slot, and the halo slots read by the rows east, north and
    ! south of the subdomain (their mirror couplings)
    real(kind=rp), dimension(:,:,:,:), pointer, intent(in) :: cA
    integer(kind=ip), intent(in) :: nx,ny,nz,dk,dj,di,k,j,i
    real(kind=rp), intent(in) :: v
    integer(kind=ip) :: n
    n = stencil_slot(dk,dj,di,k)
    if (n > 0) cA(n,k,j,i) = v
    if (i == nx) then
      if (dk==0 .and. dj==0 .and. di==1)               cA(7,k,j,nx+1)   = v
      if (k>=2 .and. dk==-1 .and. dj==0 .and. di==1)   cA(6,k-1,j,nx+1) = v
      if (k<=nz-1 .and. dk==1 .and. dj==0 .and. di==1) cA(8,k+1,j,nx+1) = v
      if (k==1 .and. dk==0 .and. dj==-1 .and. di==1)   cA(5,1,j-1,nx+1) = v
      if (k==1 .and. dk==0 .and. dj==1 .and. di==1)    cA(8,1,j+1,nx+1) = v
    endif
    if (j == ny) then
      if (dk==0 .and. dj==1 .and. di==0)               cA(4,k,ny+1,i)   = v
      if (k>=2 .and. dk==-1 .and. dj==1 .and. di==0)   cA(3,k-1,ny+1,i) = v
      if (k<=nz-1 .and. dk==1 .and. dj==1 .and. di==0) cA(5,k+1,ny+1,i) = v
      if (k==1 .and. dk==0 .and. dj==1 .and. di==1)    cA(8,1,ny+1,i+1) = v
    endif
    if (j == 1) then
      if (k==1 .and. dk==0 .and. dj==-1 .and. di==1)   cA(5,1,0,i+1)    = v
    endif
  end subroutine put_entry

  !-------------------------------------------------------------------------
  subroutine refresh_wall_rows()
    ! Matrix recompute: the formulas are exact except in the rows of the
    ! cells touching a physical wall (masked faces). Re-extract those rows
    ! with the 27 colourings restricted to a band along each wall of this
    ! rank: the operator is applied in the band only, and the colour
    ! pattern is a function of global indices, so the seam halos are
    ! filled without communication. Ranks without walls do nothing.
    integer(kind=ip) :: nx,ny,nz,i,j,k,a,bb,c,ig,jg,pi,pj
    real(kind=rp), dimension(:,:,:), pointer :: p,du,dv,dw
    real(kind=rp), dimension(:,:,:,:), pointer :: cA
    logical :: west,east,south,north

    west  = grid(1)%neighb(4) == MPI_PROC_NULL
    east  = grid(1)%neighb(2) == MPI_PROC_NULL
    south = grid(1)%neighb(1) == MPI_PROC_NULL
    north = grid(1)%neighb(3) == MPI_PROC_NULL
    if (.not.(west.or.east.or.south.or.north)) return

    nx = grid(1)%nx;  ny = grid(1)%ny;  nz = grid(1)%nz
    p  => grid(1)%p;  du => grid(1)%du;  dv => grid(1)%dv;  dw => grid(1)%dw
    cA => grid(1)%cA
    pj = myrank/grid(1)%npx
    pi = mod(myrank,grid(1)%npx)

    do a = 0,2
      do bb = 0,2
        do c = 0,2
          do i = 0,nx+1
            ig = i + pi*nx
            do j = 0,ny+1
              jg = j + pj*ny
              do k = 1,nz
                p(k,j,i) = merge(one, zero, mod(ig+3,3)==a .and. mod(jg+3,3)==bb .and. mod(k,3)==c)
              enddo
            enddo
          enddo
          ! wall halos as fill_halo sets them (not seen by the masked rows)
          if (west)  p(:,1:ny,0)    = p(:,1:ny,1)
          if (east)  p(:,1:ny,nx+1) = p(:,1:ny,nx)
          if (south) p(:,0,1:nx)    = p(:,1,1:nx)
          if (north) p(:,ny+1,1:nx) = p(:,ny,1:nx)
          dw(1,:,:) = zero
          if (west)  call correction_uvw(i0=0,    i1=3)
          if (east)  call correction_uvw(i0=nx-2, i1=nx+1)
          if (south) call correction_uvw(j0=0,    j1=3)
          if (north) call correction_uvw(j0=ny-2, j1=ny+1)
          if (west)  call rows(1, 1, 1,ny)
          if (east)  call rows(nx,nx,1,ny)
          if (south) call rows(1,nx,1, 1)
          if (north) call rows(1,nx,ny,ny)
        enddo
      enddo
    enddo
    p = zero

  contains

    subroutine rows(ia,ib,ja,jb)
      integer(kind=ip), intent(in) :: ia,ib,ja,jb
      integer(kind=ip) :: di,dj,dk,m
      real(kind=rp) :: v
      do i = ia,ib
        ig = i + pi*nx
        di = -2                                ! offset with colour a
        do m = -1,1
          if (mod(ig+m+3,3) == a) di = m
        enddo
        do j = ja,jb
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
            v = -( du(k,j,i+1) - du(k,j,i)     &
                 + dv(k,j+1,i) - dv(k,j,i)     &
                 + dw(k+1,j,i) - dw(k,j,i) )
            call put_entry(cA,nx,ny,nz,v,dk,dj,di,k,j,i)
          enddo
        enddo
      enddo
    end subroutine rows

  end subroutine refresh_wall_rows

  !-------------------------------------------------------------------------     
  subroutine correction_uvw(i0,i1,j0,j1)
    ! optional bounds restrict the i and j loops (band extraction of the
    ! wall rows at a matrix recompute); default = whole subdomain
    integer(kind=ip), optional, intent(in) :: i0,i1,j0,j1

    !! u,v,w are fluxes, the correction is T*grad(p)

    integer(kind=ip):: k, j, i
    integer(kind=ip):: nx, ny, nz, ia, ib, ja, jb

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
    ia = 0;  ib = nx+1;  ja = 0;  jb = ny+1
    if (present(i0)) ia = i0
    if (present(i1)) ib = i1
    if (present(j0)) ja = j0
    if (present(j1)) jb = j1

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

    do i = max(1,ia),min(nx+1,ib)
        do j = max(0,ja),min(ny+1,jb)
          do k = 1,nz
             px(k,j,i) = -one / dxu(j,i) * (p(k,j,i)-p(k,j,i-1))
          enddo
       enddo
    enddo

    do i = max(0,ia),min(nx+1,ib)
       do j = max(1,ja),min(ny+1,jb)
          do k = 1,nz
             py(k,j,i) = -one / dyv(j,i) * (p(k,j,i)-p(k,j-1,i))
          enddo
       enddo
    enddo

    ! column masking: no pressure gradient through no-flux faces
    if (allocated(umsk)) then
       do i = max(1,ia),min(nx+1,ib)
          do j = max(0,ja),min(ny+1,jb)
             px(:,j,i) = px(:,j,i)*umsk(j,i)
          enddo
       enddo
       do i = max(0,ia),min(nx+1,ib)
          do j = max(1,ja),min(ny+1,jb)
             py(:,j,i) = py(:,j,i)*vmsk(j,i)
          enddo
       enddo
    endif

    do i = max(0,ia),min(nx+1,ib)
       do j = max(0,ja),min(ny+1,jb)
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

    do i = max(1,ia),min(nx+1,ib)
       do j = max(1,ja),min(ny,jb)
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

    do i = max(1,ia),min(nx,ib)
       do j = max(1,ja),min(ny+1,jb)
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

    do i = max(1,ia),min(nx,ib)
       do j = max(1,ja),min(ny,jb)
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
       do i = max(1,ia),min(nx+1,ib)
          do j = max(1,ja),min(ny,jb)
             du(:,j,i) = du(:,j,i)*umsk(j,i)
          enddo
       enddo
       do i = max(1,ia),min(nx,ib)
          do j = max(1,ja),min(ny+1,jb)
             dv(:,j,i) = dv(:,j,i)*vmsk(j,i)
          enddo
       enddo
    endif

  end subroutine correction_uvw

end module mg_projection

