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
  real(kind=rp), dimension(:,:), allocatable, target :: umsk, vmsk
  ! rmsk(j,i): cell mask (1 fluid, 0 land); have_mask: an interior land
  ! mask exists somewhere (global), which switches the coarse operators
  ! to Galerkin aggregation and adds identity rows for land
  real(kind=rp), dimension(:,:), allocatable :: rmsk
  logical :: have_mask = .false.

contains
  !-----------------------------------------------------------------------------------
  subroutine set_face_masks()
    ! walls of the domain perimeter (zeroed on top of any land mask set
    ! by set_land_masks)
    integer(kind=ip) :: nx, ny
    nx = grid(1)%nx
    ny = grid(1)%ny
    if (.not.allocated(umsk)) then
       allocate(umsk(0:ny+1,0:nx+1));  umsk = one
    endif
    if (.not.allocated(vmsk)) then
       allocate(vmsk(0:ny+1,0:nx+1));  vmsk = one
    endif
    if (grid(1)%neighb(4) == MPI_PROC_NULL) umsk(:,1   ) = zero   ! west wall
    if (grid(1)%neighb(2) == MPI_PROC_NULL) umsk(:,nx+1) = zero   ! east wall
    if (grid(1)%neighb(1) == MPI_PROC_NULL) vmsk(1   ,:) = zero   ! south wall
    if (grid(1)%neighb(3) == MPI_PROC_NULL) vmsk(ny+1,:) = zero   ! north wall
  end subroutine set_face_masks

  !-----------------------------------------------------------------------------------
  subroutine set_land_masks(rm,um,vm)
    ! interior land mask from the model (fine grid, halos included):
    ! cells rm, u faces um (face i = west face of cell i), v faces vm
    real(kind=rp), dimension(0:,0:), intent(in) :: rm,um,vm
    integer(kind=ip) :: nx,ny,ierr
    logical :: l
    nx = grid(1)%nx;  ny = grid(1)%ny
    if (.not.allocated(umsk)) allocate(umsk(0:ny+1,0:nx+1))
    if (.not.allocated(vmsk)) allocate(vmsk(0:ny+1,0:nx+1))
    if (.not.allocated(rmsk)) allocate(rmsk(0:ny+1,0:nx+1))
    umsk = um;  vmsk = vm;  rmsk = rm
    l = any(rmsk(1:ny,1:nx) == zero)
    call MPI_Allreduce(l,have_mask,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)
    grid(1)%lmsk = rmsk
    if (myrank==0 .and. have_mask) write(*,*) '    - interior land mask: coarse operators by Galerkin aggregation'
  end subroutine set_land_masks

  !-----------------------------------------------------------------------------------
  subroutine set_level_masks()
    ! fluid fraction of the coarse columns: mean over the 2x2 block,
    ! with the gather protocol of dx (see set_obc_leak)
    integer(kind=ip) :: lev,nx,ny,nxf,nyf,nxc,nyc,i,j
    real(kind=rp), dimension(:,:), pointer :: lf,lc
    do lev = 2,nlevs
       nx = grid(lev)%nx;  ny = grid(lev)%ny
       nxf = grid(lev-1)%nx;  nyf = grid(lev-1)%ny
       lf => grid(lev-1)%lmsk
       if (grid(lev)%gather == 1) then
          nxc = nx/grid(lev)%ngx;  nyc = ny/grid(lev)%ngy
          allocate(lc(0:nyc+1,0:nxc+1))
       else
          nxc = nx;  nyc = ny
          lc => grid(lev)%lmsk
       endif
       lc = one
       do i = 1,nxc
          do j = 1,nyc
             lc(j,i) = qrt*(lf(2*j-1,2*i-1) + lf(2*j,2*i-1) + lf(2*j-1,2*i) + lf(2*j,2*i))
          enddo
       enddo
       if (grid(lev)%gather == 1) then
          call gather(lev,lc,grid(lev)%lmsk)
          deallocate(lc)
       endif
       call fill_halo(lev,grid(lev)%lmsk)
    enddo
  end subroutine set_level_masks

  !-----------------------------------------------------------------------------------
  subroutine land_rows(lev)
    ! identity rows for land columns: no couplings (the ones through
    ! masked faces are already zero), unit diagonal, so p = b = 0 there
    integer(kind=ip), intent(in) :: lev
    integer(kind=ip) :: i,j,k,nx,ny,nz
    nx = grid(lev)%nx;  ny = grid(lev)%ny;  nz = grid(lev)%nz
    do i = 1,nx
       do j = 1,ny
          if (grid(lev)%lmsk(j,i) == zero) then
             do k = 1,nz
                grid(lev)%cA(:,k,j,i) = zero
                grid(lev)%cA(1,k,j,i) = -one
             enddo
             grid(lev)%cAg(:,j,i) = zero
             grid(lev)%cAg(1,j,i) = -one
          endif
       enddo
    enddo
  end subroutine land_rows

  !-----------------------------------------------------------------------------------
  subroutine set_matrices(first)

    logical, optional, intent(in) :: first   ! operator check + stage prints
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
    real(kind=rp), dimension(:,:,:,:), pointer :: cA

    integer(kind=ip) :: dirichlet_flag 
    real(kind=rp)    :: sig, fj, fjm, fi, fim, dm, dn, sfi
    ! face masks of the level: the fine-grid arrays (walls, land) on
    ! level 1, the walls of the level above it
    real(kind=rp), dimension(:,:), pointer :: um, vm
    real(kind=rp), dimension(:,:), allocatable, target :: umw, vmw

    if (surface_neumann) then
       dirichlet_flag = 0
    else
       dirichlet_flag = 1
    endif

    if (allocated(obcrob_u)) call set_obc_leak()
    if (have_mask) call set_level_masks()
    call stage('face masks')
    call set_face_masks()

    do lev = 1, nlevs

       nx = grid(lev)%nx
       ny = grid(lev)%ny
       nz = grid(lev)%nz

       ! coarse operator by Galerkin aggregation of the finer one: always
       ! for a 2D level (nz = 1: the 3D formulas below do not apply to its
       ! 5-point storage) and for an interior land mask (partially open
       ! coarse faces), otherwise on request (coarse_galerkin)
       if (lev > 1 .and. (coarse_galerkin .or. have_mask .or. nz == 1)) then
          call set_galerkin(lev)
          if (have_mask) call land_rows(lev)
          full = .true.
          if (present(first)) full = first
          if (full) call level_check(lev)
          cycle
       endif

       dx    => grid(lev)%dx    !
       dy    => grid(lev)%dy    !
       dxu   => grid(lev)%dxu   !
       dyv   => grid(lev)%dyv   !
       dzw   => grid(lev)%dzw   !
       Arx   => grid(lev)%Arx   !
       Ary   => grid(lev)%Ary   !
       Arz   => grid(lev)%Arz   !
       alpha => grid(lev)%alpha !
       zxdy  => grid(lev)%zxdy  !
       zydx  => grid(lev)%zydx  !
       cA    => grid(lev)%cA    !

       !! interaction coeff with neighbours !!

       !---------------!
       !- lower level -!
       !---------------!
       k = 1
       ! bottom level (ghost pressure level below it): the k=1 rows are
       ! the mirror of the Dirichlet surface rows, with the ghost value
       ! pg in place of the zero surface pressure.  Derived from
       ! correction_uvw: bottom face = kinematic row with the half-cell
       ! weights 1/2, k=1 u,v rows = interior form with the factor-two
       ! half-cell term; the k-1 slots (2,5,8) couple to the ghost rows.
       ! Checked against the operator by operator_check on the first call.
       do i = 1,nx
          do j = 1,ny+1
             ! couples with k+1,j-1
             cA(3,k,j,i) = qrt * ( zydx(k+1,j,i) + zydx(k,j-1,i) )
             ! couples with j-1
             cA(4,k,j,i) = Ary(k,j,i) / dyv(j,i) + qrt * ( zydx(k,j-1,i) - zydx(k,j,i) )
             ! couples with the ghost row at j-1
             cA(5,k,j,i) = - hlf * zydx(k,j-1,i)
          enddo
       enddo
       do i = 1,nx+1
          do j = 1,ny
             ! couples with k+1,i-1
             cA(6,k,j,i) = qrt * ( zxdy(k+1,j,i) + zxdy(k,j,i-1) )
             ! couples with i-1
             cA(7,k,j,i) = Arx(k,j,i) / dxu(j,i) + qrt * ( zxdy(k,j,i-1) - zxdy(k,j,i) )
             ! couples with the ghost row at i-1
             cA(8,k,j,i) = - hlf * zxdy(k,j,i-1)
          enddo
       enddo
       do i = 1,nx
          do j = 1,ny
             ! couples with the ghost row below
             cA(2,k,j,i) = Arz(j,i) / dzw(k,j,i) * alpha(k,j,i)
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
             ! mirror of the surface diagonal: full vertical coupling
             ! to the ghost row through the half cell dzw(1)
             cA(1,k,j,i) = &
                  - Arz(j,i) / dzw(k+1,j,i) * hlf * (alpha(k+1,j,i) + alpha(k,j,i)) &
                  - Arz(j,i) / dzw(k  ,j,i) * alpha(k,j,i) &
                  - Arx(k,j,i  )/dxu(j,i  )  &
                  - Arx(k,j,i+1)/dxu(j,i+1)  &
                  - Ary(k,j  ,i)/dyv(j  ,i)  &
                  - Ary(k,j+1,i)/dyv(j+1,i)

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

       ! masked faces (walls, open faces, land): the operator is
       ! -D (M T M) G with M the face masks, so every lateral and
       ! cross coupling, which goes through one lateral face, carries
       ! that face's mask.  The slope terms of the vertical couplings
       ! and of the k=1 and k=nz diagonals cancel between the two
       ! faces of a cell only when both are open: their residual is
       ! proportional to the mask jump across the cell.  Open faces
       ! are masked, with the leak as a separate diagonal term.  The
       ! coarse levels (walls only: whole faces at every level) take
       ! the same formulas on their coarsened metrics.
       if (lev == 1) then
          um => umsk;  vm => vmsk
       else
          allocate(umw(0:ny+1,0:nx+1), vmw(0:ny+1,0:nx+1))
          umw = one;  vmw = one
          if (grid(lev)%neighb(4) == MPI_PROC_NULL) umw(:,1   ) = zero
          if (grid(lev)%neighb(2) == MPI_PROC_NULL) umw(:,nx+1) = zero
          if (grid(lev)%neighb(1) == MPI_PROC_NULL) vmw(1   ,:) = zero
          if (grid(lev)%neighb(3) == MPI_PROC_NULL) vmw(ny+1,:) = zero
          um => umw;  vm => vmw
       endif
       do i = 1,nx
          do j = 1,ny+1
             do k = 1,nz
                cA(3,k,j,i) = cA(3,k,j,i) * vm(j,i)
                cA(4,k,j,i) = cA(4,k,j,i) * vm(j,i)
                cA(5,k,j,i) = cA(5,k,j,i) * vm(j,i)
             enddo
          enddo
       enddo
       do i = 1,nx+1
          do j = 1,ny
             do k = 1,nz
                cA(6,k,j,i) = cA(6,k,j,i) * um(j,i)
                cA(7,k,j,i) = cA(7,k,j,i) * um(j,i)
                cA(8,k,j,i) = cA(8,k,j,i) * um(j,i)
             enddo
          enddo
       enddo
       do i = 1,nx
          do j = 1,ny
             dm = um(j,i+1) - um(j,i)
             dn = vm(j+1,i) - vm(j,i)
             sfi = sigtop(dzw(nz+1,j,i))
             if (lev == 1) sfi = sfi * sfcf(j,i)
             do k = 1,nz
                ! direct couplings of the diagonal through masked faces
                cA(1,k,j,i) = cA(1,k,j,i) + (one - um(j,i  )) * Arx(k,j,i  )/dxu(j,i  ) &
                                          + (one - um(j,i+1)) * Arx(k,j,i+1)/dxu(j,i+1) &
                                          + (one - vm(j  ,i)) * Ary(k,j  ,i)/dyv(j  ,i) &
                                          + (one - vm(j+1,i)) * Ary(k,j+1,i)/dyv(j+1,i)
             enddo
             cA(2,1,j,i) = cA(2,1,j,i) + hlf * ( zxdy(1,j,i)*dm + zydx(1,j,i)*dn )
             cA(1,1,j,i) = cA(1,1,j,i) - hlf * ( zxdy(1,j,i)*dm + zydx(1,j,i)*dn )
             do k = 2,nz
                cA(2,k,j,i) = cA(2,k,j,i) + qrt * ( (zxdy(k,j,i) - zxdy(k-1,j,i))*dm &
                                                  + (zydx(k,j,i) - zydx(k-1,j,i))*dn )
             enddo
             cA(1,nz,j,i) = cA(1,nz,j,i) + hlf * (two*sfi - one) * ( zxdy(nz,j,i)*dm + zydx(nz,j,i)*dn )
          enddo
       enddo
       if (lev == 1 .and. allocated(obcrob_u)) then
          ! radiating open faces: du = obcrob*Ar*p(adjacent cell)
          ! (coarse levels: the coarsened leak, below)
          if (grid(1)%neighb(4) == MPI_PROC_NULL) then
             do j = 1,ny
                do k = 1,nz
                   cA(1,k,j,1)  = cA(1,k,j,1)  + obcrob_u(j,1)   *Arx(k,j,1)
                enddo
             enddo
          endif
          if (grid(1)%neighb(2) == MPI_PROC_NULL) then
             do j = 1,ny
                do k = 1,nz
                   cA(1,k,j,nx) = cA(1,k,j,nx) - obcrob_u(j,nx+1)*Arx(k,j,nx+1)
                enddo
             enddo
          endif
          if (grid(1)%neighb(1) == MPI_PROC_NULL) then
             do i = 1,nx
                do k = 1,nz
                   cA(1,k,1,i)  = cA(1,k,1,i)  + obcrob_v(1,i)   *Ary(k,1,i)
                enddo
             enddo
          endif
          if (grid(1)%neighb(3) == MPI_PROC_NULL) then
             do i = 1,nx
                do k = 1,nz
                   cA(1,k,ny,i) = cA(1,k,ny,i) - obcrob_v(ny+1,i)*Ary(k,ny+1,i)
                enddo
             enddo
          endif
       endif

       ! ghost rows (bottom pressure level): the kinematic row -dw(1)
       ! of T, written over the halo as well, since the
       ! slots read across seams and by the rows next to the east and
       ! north edges are never halo-filled (a fill would mirror the
       ! interior couplings into a wall halo)
       do i = 0,nx+1
          do j = 0,ny+1
             grid(lev)%cAg(1,j,i) = - Arz(j,i) / dzw(1,j,i) * alpha(1,j,i)   ! diagonal
             grid(lev)%cAg(2,j,i) =   Arz(j,i) / dzw(1,j,i) * alpha(1,j,i)   ! p(1,j,i)
             grid(lev)%cAg(3,j,i) =   hlf * zydx(1,j,i)                      ! p(1,j-1,i)
             grid(lev)%cAg(4,j,i) = - hlf * zydx(1,j,i)                      ! p(1,j+1,i)
             grid(lev)%cAg(5,j,i) =   hlf * zxdy(1,j,i)                      ! p(1,j,i-1)
             grid(lev)%cAg(6,j,i) = - hlf * zxdy(1,j,i)                      ! p(1,j,i+1)
          enddo
       enddo
       ! masked faces, as for the cell rows (the halo slots read
       ! across seams and walls, cAg(3,ny+1,:) and cAg(5,:,nx+1),
       ! included; the other halo slots are never read)
       do i = 0,nx+1
          do j = 0,ny+1
             grid(lev)%cAg(3,j,i) = grid(lev)%cAg(3,j,i) * vm(j,i)
             grid(lev)%cAg(5,j,i) = grid(lev)%cAg(5,j,i) * um(j,i)
             if (j <= ny) grid(lev)%cAg(4,j,i) = grid(lev)%cAg(4,j,i) * vm(j+1,i)
             if (i <= nx) grid(lev)%cAg(6,j,i) = grid(lev)%cAg(6,j,i) * um(j,i+1)
          enddo
       enddo
       do i = 1,nx
          do j = 1,ny
             dm = um(j,i+1) - um(j,i)
             dn = vm(j+1,i) - vm(j,i)
             grid(lev)%cAg(2,j,i) = grid(lev)%cAg(2,j,i) + hlf * ( zxdy(1,j,i)*dm + zydx(1,j,i)*dn )
          enddo
       enddo
       if (lev > 1) deallocate(umw, vmw)

       ! coarse levels: the open-boundary leak (fine grid: in the extracted
       ! stencil through correction_uvw) as a diagonal term per cell; the
       ! cycle diverges when the coarse operators lack it (COOL, 2026-09-09)
       if (lev > 1 .and. allocated(obcrob_u)) then
          do i = 1,nx
             do j = 1,ny
                do k = 1,nz
                   cA(1,k,j,i) = cA(1,k,j,i) - grid(lev)%leak(j,i)*grid(lev)%dz(k,j,i)
                enddo
             enddo
          enddo
       endif

       if (netcdf_output) then
          if (myrank==0) write(*,*)'       write cA in a netcdf file'
          call write_netcdf(grid(lev)%cA,vname='ca',netcdf_file_name='cA.nc',rank=myrank,iter=lev)
       endif

       ! Fine grid: the stencil of the MASKED operator  A = -D (M T M) G,
       ! so that the no-flux faces of the domain perimeter and of an
       ! interior land mask are represented by zero rows/columns of T
       ! instead of a mirrored p.  Mirroring is a zero-flux condition
       ! only for the Arx*px part of the U row; with slope cross terms
       ! it makes the effective operator non-symmetric, which under the
       ! singular Neumann surface condition leaves an irreducible
       ! residual (LOCK, 2026-09-04).  The formulas cover the masked
       ! rows; on the first call the matrix is checked against the
       ! operator as correction_uvw applies it (operator_check).
       full = .true.
       if (present(first)) full = first
       if (full) call level_check(lev)
       if (lev == 1) then
          if (full) then
             call stage('operator check')
             call operator_check()
          endif
          if (have_mask) call land_rows(1)
       endif
       if (lev == 1 .and. nlevs > 1) call stage('coarse levels')


    enddo

  contains
    subroutine level_check(lev)
      ! first call: per-level sanity of the assembled matrix (diagonal
      ! range over interior cells, NaN count), so a bad coarse level or
      ! gather shows in the start-up log
      integer(kind=ip), intent(in) :: lev
      integer(kind=ip) :: i,j,k,ierr,nx,ny,nz
      real(kind=rp) :: dmin,dmax,nnan,l(3),g(3)
      real(kind=rp), dimension(:,:,:,:), pointer :: cA
      nx = grid(lev)%nx;  ny = grid(lev)%ny;  nz = grid(lev)%nz
      cA => grid(lev)%cA
      dmin = huge(dmin);  dmax = zero;  nnan = zero
      do i = 1,nx
        do j = 1,ny
          do k = 1,nz
            if (any(cA(:,k,j,i) /= cA(:,k,j,i))) nnan = nnan + one
            dmin = min(dmin, abs(cA(1,k,j,i)))
            dmax = max(dmax, abs(cA(1,k,j,i)))
          enddo
        enddo
      enddo
        ! ghost rows counted with the cells (their diagonal is a bottom
        ! half-cell vertical coupling)
        do i = 1,nx
       do j = 1,ny
         if (any(grid(lev)%cAg(:,j,i) /= grid(lev)%cAg(:,j,i))) nnan = nnan + one
         dmin = min(dmin, abs(grid(lev)%cAg(1,j,i)))
         dmax = max(dmax, abs(grid(lev)%cAg(1,j,i)))
       enddo
        enddo
      l = (/ -dmin, dmax, nnan /)
      call MPI_Allreduce(l,g,3,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(nnan,l(3),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      if (myrank==0) write(*,'(A,I2,A,I0,A,I0,A,I0,A,I0,A,ES9.2,A,ES9.2,A,I0)') &
           '     level ',lev,': ',nx,'x',ny,'x',nz,' gather ',grid(lev)%gather, &
           '  |diag| min ',-g(1),' max ',g(2),'  NaN cells ',int(l(3))
    end subroutine level_check

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
  subroutine set_galerkin(lev)
    ! Coarse operator of level lev by piecewise-constant aggregation of
    ! the operator of level lev-1: R = sum over the 2x2x2 block (the
    ! restriction in use), P = R^T.  The gradient couplings get the
    ! factor 1/2 that makes RAP equal the geometric coarse operator on a
    ! uniform grid (a block-constant function sees the fine spacing);
    ! the diagonal is rebuilt as minus the coupling sum plus the surface
    ! (Robin/Dirichlet) and open-boundary leak terms from the coarse
    ! metrics, which are area terms and carry no 1/2.  The result fits
    ! the same 15-point storage (5-point at a 2D level); masked or
    ! leaking fine faces coarsen consistently.  Each fine coupling is
    ! stored once, from the fine cell whose block keeps the coarse pair;
    ! the reverse visit feeds the coupling sum only.
    integer(kind=ip), intent(in) :: lev
    integer(kind=ip) :: nxf,nyf,nzf,nxc,nyc,nzc,ndc,i,j,k,dk,dj,di,s
    integer(kind=ip) :: kc,jc,ic,kn,jn,in,ierr,m
    real(kind=rp) :: a
    real(kind=rp), dimension(:,:,:,:), pointer :: cf,cc,cA
    real(kind=rp), dimension(:,:,:), pointer :: t3,cgf,cgc,cAg
    real(kind=rp), dimension(:,:), pointer :: t2
    integer(kind=ip) :: k0,g
    integer(kind=ip), parameter :: noff = 14
    integer(kind=ip), dimension(3,noff), parameter :: off = reshape( (/ &
         -1,0,0,  1,0,0,  0,-1,0,  0,1,0,  0,0,-1,  0,0,1, &
         -1,-1,0, -1,1,0,  1,-1,0,  1,1,0, -1,0,-1, -1,0,1,  1,0,-1,  1,0,1 /), (/3,noff/) )

    nxf = grid(lev-1)%nx;  nyf = grid(lev-1)%ny;  nzf = grid(lev-1)%nz
    cf => grid(lev-1)%cA
    cA => grid(lev)%cA;    ndc = size(cA,1)
    nzc = grid(lev)%nz
    if (nzf == 1) then
       if (myrank==0) write(*,*) 'set_galerkin: 2D -> 2D coarsening not implemented'
       call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    endif
    k0 = 0                                  ! ghost row = level 0
    if (grid(lev)%gather == 1) then
       nxc = grid(lev)%nx/grid(lev)%ngx;  nyc = grid(lev)%ny/grid(lev)%ngy
       allocate(cc(ndc,nzc,0:nyc+1,0:nxc+1))
       allocate(cgc(6,0:nyc+1,0:nxc+1))
    else
       nxc = grid(lev)%nx;  nyc = grid(lev)%ny
       cc => cA
       cgc => grid(lev)%cAg
    endif
    cc = zero
    cgf => grid(lev-1)%cAg;  cgc = zero

    call fill_halo(lev-1,cf)

    do i = 1,nxf
       ic = (i+1)/2
       do j = 1,nyf
          jc = (j+1)/2
          do k = k0,nzf
             kc = (k+1)/2;  if (nzc == 1) kc = 1;  if (k == 0) kc = 0
             do m = 1,noff
                dk = off(1,m);  dj = off(2,m);  di = off(3,m)
                kn = k+dk;  jn = j+dj;  in = i+di
                if (kn < k0 .or. kn > nzf) cycle
                if (jn < 0 .or. jn > nyf+1 .or. in < 0 .or. in > nxf+1) cycle
                ! the fine coupling, wherever the symmetric storage keeps it
                if (k == 0) then                          ! ghost row: its own slots
                   g = gslot(dk,dj,di)
                   if (g == 0) cycle
                   a = cgf(g,j,i)
                elseif (kn == 0) then                     ! cell row -> ghost
                   s = stencil_slot(dk,dj,di)
                   if (s > 1) then
                      a = cf(s,k,j,i)
                   else
                      g = gslot(-dk,-dj,-di)
                      if (g == 0) cycle
                      a = cgf(g,jn,in)
                   endif
                else
                   s = stencil_slot(dk,dj,di)
                   if (s > 1) then
                      a = cf(s,k,j,i)
                   else
                      s = stencil_slot(-dk,-dj,-di)
                      if (s <= 1) cycle
                      a = cf(s,kn,jn,in)
                   endif
                endif
                if (a == zero) cycle
                if (kn > 0) then
                   kn = (kn+1)/2;  if (nzc == 1) kn = 1
                endif
                jn = (jn+1)/2;  in = (in+1)/2
                if (kn == kc .and. jn == jc .and. in == ic) cycle   ! internal to the block
                if (kc == 0) then
                   cgc(1,jc,ic) = cgc(1,jc,ic) - hlf*a               ! ghost coupling sum
                   g = gslot(kn-kc,jn-jc,in-ic)
                   if (g > 1) cgc(g,jc,ic) = cgc(g,jc,ic) + hlf*a
                   cycle
                endif
                cc(1,kc,jc,ic) = cc(1,kc,jc,ic) - hlf*a           ! coupling sum
                if (nzc > 1) then
                   s = stencil_slot(kn-kc,jn-jc,in-ic)
                else
                   s = slot2d(jn-jc,in-ic)
                endif
                if (s > 1) cc(s,kc,jc,ic) = cc(s,kc,jc,ic) + hlf*a
             enddo
          enddo
       enddo
    enddo

    if (grid(lev)%gather == 1) then
       ! gather slot by slot; the gather fills interior points only, so
       ! the halo slots (walls of the gathered grid) are set to zero here
       t3 => grid(lev)%dummy3
       cA = zero
       do s = 1,ndc
          do i = 0,nxc+1;  do j = 0,nyc+1;  do k = 1,nzc
             t3(k,j,i) = cc(s,k,j,i)
          enddo;  enddo;  enddo
          call gather(lev,t3,grid(lev)%r)
          do i = 1,grid(lev)%nx;  do j = 1,grid(lev)%ny;  do k = 1,nzc
             cA(s,k,j,i) = grid(lev)%r(k,j,i)
          enddo;  enddo;  enddo
       enddo
       deallocate(cc)
       grid(lev)%r = zero
       cAg => grid(lev)%cAg;  cAg = zero
       t2 => grid(lev)%dummy2
       do s = 1,6
          do i = 0,nxc+1;  do j = 0,nyc+1
             t2(j,i) = cgc(s,j,i)
          enddo;  enddo
          call gather(lev,t2,grid(lev)%rg)
          do i = 1,grid(lev)%nx;  do j = 1,grid(lev)%ny
             cAg(s,j,i) = grid(lev)%rg(j,i)
          enddo;  enddo
       enddo
       deallocate(cgc)
       grid(lev)%rg = zero
    endif

    ! surface and open-boundary leak terms (area terms: no 1/2)
    do i = 1,grid(lev)%nx
       do j = 1,grid(lev)%ny
          cA(1,nzc,j,i) = cA(1,nzc,j,i) - grid(lev)%Arz(j,i)/grid(lev)%dzw(nzc+1,j,i) &
               * grid(lev)%alpha(nzc,j,i) * sigtop(grid(lev)%dzw(nzc+1,j,i)) * grid(lev)%lmsk(j,i)
          if (allocated(obcrob_u)) then
             do k = 1,nzc
                cA(1,k,j,i) = cA(1,k,j,i) - grid(lev)%leak(j,i)*grid(lev)%dz(k,j,i)*grid(lev)%lmsk(j,i)
             enddo
          endif
       enddo
    enddo

    call fill_halo(lev,cA)
    call fill_halo_cAg(lev)

  contains
    function gslot(dk,dj,di) result(n)
      ! slot of the ghost row's storage for the coupling to the cell at
      ! offset (dk,dj,di): 1 self, 2..6 the five k=1 cells above
      integer(kind=ip), intent(in) :: dk,dj,di
      integer(kind=ip) :: n
      n = 0
      if (dk==0 .and. dj==0 .and. di==0) n = 1
      if (dk==1) then
         if (dj==0  .and. di==0 ) n = 2
         if (dj==-1 .and. di==0 ) n = 3
         if (dj==1  .and. di==0 ) n = 4
         if (dj==0  .and. di==-1) n = 5
         if (dj==0  .and. di==1 ) n = 6
      endif
    end function gslot
    function slot2d(dj,di) result(n)
      ! 5-point storage of a 2D level: 2 (j-1), 3 (i-1), 4 (j-1,i-1), 5 (j+1,i-1)
      integer(kind=ip), intent(in) :: dj,di
      integer(kind=ip) :: n
      n = 0
      if (dj==-1 .and. di==0) n = 2
      if (dj==0 .and. di==-1) n = 3
      if (dj==-1 .and. di==-1) n = 4
      if (dj==1 .and. di==-1) n = 5
    end function slot2d
  end subroutine set_galerkin

  !-------------------------------------------------------------------------
  subroutine set_obc_leak()
    ! Per-cell leak density of the radiating open faces, |obcrob|*(face
    ! width), so that the diagonal term of a cell is leak*dz; coarsened
    ! as a 2x2 block SUM (the block's total leak per unit dz) with the
    ! same gather protocol as dx.
    integer(kind=ip) :: lev,nx,ny,nxf,nyf,nxc,nyc,i,j
    real(kind=rp), dimension(:,:), pointer :: lf,lc

    do lev = 1,nlevs
       nx = grid(lev)%nx;  ny = grid(lev)%ny
       if (lev == 1) then
          lc => grid(1)%leak
          lc = zero
          if (grid(1)%neighb(4) == MPI_PROC_NULL) then
             do j = 1,ny
                lc(j,1)  = lc(j,1)  + abs(obcrob_u(j,1))   *grid(1)%dy(j,1)
             enddo
          endif
          if (grid(1)%neighb(2) == MPI_PROC_NULL) then
             do j = 1,ny
                lc(j,nx) = lc(j,nx) + abs(obcrob_u(j,nx+1))*grid(1)%dy(j,nx)
             enddo
          endif
          if (grid(1)%neighb(1) == MPI_PROC_NULL) then
             do i = 1,nx
                lc(1,i)  = lc(1,i)  + abs(obcrob_v(1,i))   *grid(1)%dx(1,i)
             enddo
          endif
          if (grid(1)%neighb(3) == MPI_PROC_NULL) then
             do i = 1,nx
                lc(ny,i) = lc(ny,i) + abs(obcrob_v(ny+1,i))*grid(1)%dx(ny,i)
             enddo
          endif
       else
          nxf = grid(lev-1)%nx;  nyf = grid(lev-1)%ny
          lf => grid(lev-1)%leak
          if (grid(lev)%gather == 1) then
             nxc = nx/grid(lev)%ngx;  nyc = ny/grid(lev)%ngy
             allocate(lc(0:nyc+1,0:nxc+1))
          else
             nxc = nx;  nyc = ny
             lc => grid(lev)%leak
          endif
          lc = zero
          do i = 1,nxc
             do j = 1,nyc
                lc(j,i) = lf(2*j-1,2*i-1) + lf(2*j,2*i-1) + lf(2*j-1,2*i) + lf(2*j,2*i)
             enddo
          enddo
          if (grid(lev)%gather == 1) then
             call gather(lev,lc,grid(lev)%leak)
             deallocate(lc)
          endif
       endif
       call fill_halo(lev,grid(lev)%leak)
    enddo
  end subroutine set_obc_leak


  !-------------------------------------------------------------------------
  subroutine operator_check()
    ! First-call check of the assembled fine-grid matrix against the
    ! operator it represents, -D (M T M) G as applied by correction_uvw:
    ! |A y - div(T grad y)| per row type on a deterministic vector y,
    ! which exposes any formula or halo-slot error that is not
    ! orthogonal to y, and the symmetry <y,A x> = <x,A y> of the
    ! operator itself on two vectors.  Stops on a mismatch above 1e-10
    ! of max|diag|: the recompute relies on every row of the formulas.
    ! Work arrays: b and r are free at init (y and r alias: y is
    ! consumed before the residual is formed).
    use mg_relax, only: compute_residual
    integer(kind=ip) :: nx,ny,nz,i,j,k,ig,jg,pi,pj,ierr
    real(kind=rp), dimension(:,:,:), pointer :: p,du,dv,dw,x,y,r
    real(kind=rp), dimension(:,:),   pointer :: pg
    real(kind=rp), dimension(:,:), allocatable :: xg,yg
    real(kind=rp) :: amax,syx,sxy,d,res,l(4),g(4),ls(2),gs(2)

    nx = grid(1)%nx;  ny = grid(1)%ny;  nz = grid(1)%nz
    p  => grid(1)%p;  du => grid(1)%du;  dv => grid(1)%dv;  dw => grid(1)%dw
    x  => grid(1)%b;  y  => grid(1)%r;   r  => grid(1)%r
    pg => grid(1)%pg
    allocate(xg(0:ny+1,0:nx+1), yg(0:ny+1,0:nx+1));  xg = zero;  yg = zero
    pj = myrank/grid(1)%npx
    pi = mod(myrank,grid(1)%npx)

    do i = 1,nx
      ig = i + pi*nx
      do j = 1,ny
        jg = j + pj*ny
        do k = 1,nz
          x(k,j,i) = sin(0.731_rp*ig + 1.213_rp*jg + 2.117_rp*k)
          y(k,j,i) = cos(1.379_rp*ig + 0.577_rp*jg + 1.911_rp*k)
        enddo
        xg(j,i) = sin(0.531_rp*ig + 1.713_rp*jg)
        yg(j,i) = cos(1.179_rp*ig + 0.877_rp*jg)
      enddo
    enddo

    ! <y,A x>
    do i = 0,nx+1                       ! loops: pointer targets may alias,
      do j = 0,ny+1                     ! a whole-array copy would build a
        do k = 1,nz                     ! full-field temporary on the stack
          p(k,j,i) = x(k,j,i)
        enddo
      enddo
    enddo
    pg = xg
    call fill_halo(1,p);  dw(1,:,:) = zero;  call correction_uvw()
    syx = zero
    do i = 1,nx
      do j = 1,ny
        do k = 1,nz
          syx = syx - y(k,j,i) * ( du(k,j,i+1) - du(k,j,i)  &
                                 + dv(k,j+1,i) - dv(k,j,i)  &
                                 + dw(k+1,j,i) - dw(k,j,i) )
        enddo
        syx = syx - yg(j,i) * dw(1,j,i)
      enddo
    enddo
    ! <x,A y>; du,dv,dw then hold T grad y
    do i = 0,nx+1
      do j = 0,ny+1
        do k = 1,nz
          p(k,j,i) = y(k,j,i)
        enddo
      enddo
    enddo
    pg = yg
    call fill_halo(1,p);  dw(1,:,:) = zero;  call correction_uvw()
    sxy = zero
    do i = 1,nx
      do j = 1,ny
        do k = 1,nz
          sxy = sxy - x(k,j,i) * ( du(k,j,i+1) - du(k,j,i)  &
                                 + dv(k,j+1,i) - dv(k,j,i)  &
                                 + dw(k+1,j,i) - dw(k,j,i) )
        enddo
        sxy = sxy - xg(j,i) * dw(1,j,i)
      enddo
    enddo
    ! A y by the matrix: b = 0 -> r = -A p
    x = zero;  grid(1)%bg = zero
    call compute_residual(1,res)
    l = zero
    do i = 1,nx
      do j = 1,ny
        do k = 1,nz
          d = ( du(k,j,i+1) - du(k,j,i) + dv(k,j+1,i) - dv(k,j,i) + dw(k+1,j,i) - dw(k,j,i) )
          if (k == 1) then
            l(1) = max(l(1), abs(d - r(k,j,i)))
          else
            l(2) = max(l(2), abs(d - r(k,j,i)))
          endif
          l(4) = max(l(4), abs(grid(1)%cA(1,k,j,i)))
        enddo
        l(3) = max(l(3), abs(dw(1,j,i) - grid(1)%rg(j,i)))
      enddo
    enddo
    call MPI_Allreduce(l,g,4,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
    ls = (/ syx, sxy /)
    call MPI_Allreduce(ls,gs,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    amax = g(4)
    if (myrank==0) then
       write(*,'(A,3ES10.2)') '     operator check |A y - div(T grad y)|/max|diag|: k=1, k>1, ghost ', g(1:3)/amax
       write(*,'(A,ES9.2)')   '       operator symmetry |<y,Ax>-<x,Ay>|/|<y,Ax>| ', &
            abs(gs(1)-gs(2))/max(abs(gs(1)),tiny(one))
    endif
    p = zero;  x = zero;  r = zero
    pg = zero;  grid(1)%bg = zero;  grid(1)%rg = zero
    deallocate(xg, yg)
    if (maxval(g(1:3))/amax > 1e-10_rp) then
       if (myrank==0) write(*,*) 'operator_check: matrix does not match the operator -- stopping'
       call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    endif
  end subroutine operator_check

  !-------------------------------------------------------------------------
  function stencil_slot(dk,dj,di) result(n)
    ! slot (1..8) of the 15-point storage that holds, at a cell, the
    ! coupling to the neighbour at offset (dk,dj,di); 0 if that
    ! coupling is stored at the neighbour instead (or is off-pattern)
    integer(kind=ip), intent(in) :: dk,dj,di
    integer(kind=ip) :: n
    n = 0
    if (dk==0 .and. dj==0 .and. di==0) then;        n = 1
    elseif (dk==-1 .and. dj==0 .and. di==0) then;   n = 2
    elseif (dk==1 .and. dj==-1 .and. di==0) then;   n = 3
    elseif (dk==0 .and. dj==-1 .and. di==0) then;   n = 4
    elseif (dk==1 .and. dj==0 .and. di==-1) then;   n = 6
    elseif (dk==0 .and. dj==0 .and. di==-1) then;   n = 7
    elseif (dk==-1 .and. dj==-1 .and. di==0) then;  n = 5   ! (k = 1: k-1 is the ghost row)
    elseif (dk==-1 .and. dj==0 .and. di==-1) then;  n = 8
    endif
  end function stencil_slot




  !-------------------------------------------------------------------------     
  subroutine fill_halo_cAg(lev)
    ! halo exchange of the six ghost-row coefficients (2D fills)
    integer(kind=ip), intent(in) :: lev
    integer(kind=ip) :: s,i,j,nx,ny
    real(kind=rp), dimension(:,:), pointer :: t2
    nx = grid(lev)%nx;  ny = grid(lev)%ny
    allocate(t2(0:ny+1,0:nx+1))
    do s = 1,6
       do i = 0,nx+1;  do j = 0,ny+1
          t2(j,i) = grid(lev)%cAg(s,j,i)
       enddo;  enddo
       call fill_halo(lev,t2)
       do i = 0,nx+1;  do j = 0,ny+1
          grid(lev)%cAg(s,j,i) = t2(j,i)
       enddo;  enddo
    enddo
    deallocate(t2)
  end subroutine fill_halo_cAg
  subroutine correction_uvw()
    !! u,v,w are fluxes, the correction is T*grad(p)

    integer(kind=ip):: k, j, i
    integer(kind=ip):: nx, ny, nz

    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
    real(kind=rp), dimension(:,:)  , pointer :: dxu,dyv
    real(kind=rp), dimension(:,:)  , pointer :: Arz
    real(kind=rp), dimension(:,:,:), pointer :: dzw
    real(kind=rp), dimension(:,:,:), pointer :: Arx,Ary
    real(kind=rp), dimension(:,:,:), pointer :: zxdy,zydx
    real(kind=rp), dimension(:,:,:), pointer :: alpha
    real(kind=rp), dimension(:,:,:), pointer :: p
    real(kind=rp), dimension(:,:)  , pointer :: pg
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
    zxdy  => grid(1)%zxdy
    zydx  => grid(1)%zydx
    p     => grid(1)%p
    pg => grid(1)%pg

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
          k = 1 ! bottom: the mirror of the surface value, from the ghost
                ! pressure p(0) = pg of the column
          pz(k,j,i) = -one / dzw(k,j,i) * (p(k,j,i)-pg(j,i))

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
          ! interior form with the factor-two half-cell
          ! term at the bottom face (mirror of the surface row)
          du(k,j,i) = Arx(k,j,i) * px(k,j,i) &
               - qrt * ( &
               + zxdy(k,j,i  ) * two * dzw(k  ,j,i  ) * pz(k  ,j,i  ) &
               + zxdy(k,j,i  ) *       dzw(k+1,j,i  ) * pz(k+1,j,i  ) &
               + zxdy(k,j,i-1) * two * dzw(k  ,j,i-1) * pz(k  ,j,i-1) &
               + zxdy(k,j,i-1) *       dzw(k+1,j,i-1) * pz(k+1,j,i-1) )

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
          dv(k,j,i) = Ary(k,j,i) * py(k,j,i) &
               - qrt * ( &
               + zydx(k,j  ,i) * two * dzw(k  ,j  ,i) * pz(k  ,j  ,i) &
               + zydx(k,j  ,i) *       dzw(k+1,j  ,i) * pz(k+1,j  ,i) &
               + zydx(k,j-1,i) * two * dzw(k  ,j-1,i) * pz(k  ,j-1,i) &
               + zydx(k,j-1,i) *       dzw(k+1,j-1,i) * pz(k+1,j-1,i) )

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
          ! bottom face: the kinematic row of T (mirror
          ! of the Dirichlet surface row); W(1) = 0 is enforced by the
          ! projection through the ghost pressure
          k = 1
          dw(k,j,i) = alpha(k,j,i) * Arz(j,i) * pz(k,j,i) &
               - hlf * ( &
               + zxdy(k,j,i) * dxu(j,i  ) * px(k,j,i  ) &
               + zxdy(k,j,i) * dxu(j,i+1) * px(k,j,i+1) ) &
               - hlf * ( &
               + zydx(k,j,i) * dyv(j  ,i) * py(k,j  ,i) &
               + zydx(k,j,i) * dyv(j+1,i) * py(k,j+1,i) )
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

    ! radiating open faces (implicit Flather, masked otherwise): the face
    ! flux leaks in proportion to the adjacent cell's pressure,
    ! du = n Arx/(theta c dt) p, which summed over the column is Flather
    ! at the new time level with zeta(n+1)-zeta(n) = pbar/(theta g dt)
    if (allocated(obcrob_u)) then
       if (grid(1)%neighb(4) == MPI_PROC_NULL) then
          do j = 1,ny
             do k = 1,nz
                du(k,j,1) = du(k,j,1) + obcrob_u(j,1)*Arx(k,j,1)*p(k,j,1)
             enddo
          enddo
       endif
       if (grid(1)%neighb(2) == MPI_PROC_NULL) then
          do j = 1,ny
             do k = 1,nz
                du(k,j,nx+1) = du(k,j,nx+1) + obcrob_u(j,nx+1)*Arx(k,j,nx+1)*p(k,j,nx)
             enddo
          enddo
       endif
       if (grid(1)%neighb(1) == MPI_PROC_NULL) then
          do i = 1,nx
             do k = 1,nz
                dv(k,1,i) = dv(k,1,i) + obcrob_v(1,i)*Ary(k,1,i)*p(k,1,i)
             enddo
          enddo
       endif
       if (grid(1)%neighb(3) == MPI_PROC_NULL) then
          do i = 1,nx
             do k = 1,nz
                dv(k,ny+1,i) = dv(k,ny+1,i) + obcrob_v(ny+1,i)*Ary(k,ny+1,i)*p(k,ny,i)
             enddo
          enddo
       endif
    endif

  end subroutine correction_uvw

end module mg_projection

