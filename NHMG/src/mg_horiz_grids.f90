module mg_horiz_grids

  use mg_cst
  use mg_mpi
  use mg_tictoc
  use mg_grids
  use mg_mpi_exchange
  use mg_gather
  use mg_namelist

  implicit none

contains

  !----------------------------------------
  subroutine set_horiz_grids()

    integer(kind=ip) :: lev
    integer(kind=ip) :: nx,ny
    integer(kind=ip) :: nxf,nxc
    integer(kind=ip) :: nyf,nyc

    real(kind=rp), dimension(:,:), pointer :: dxf,dxc
    real(kind=rp), dimension(:,:), pointer :: dyf,dyc
    real(kind=rp), dimension(:,:), pointer :: dxuf,dxuc
    real(kind=rp), dimension(:,:), pointer :: dyvf,dyvc

    do lev = 1, nlevs

       nx=grid(lev)%nx
       ny=grid(lev)%ny

       if (lev == 1) then

          grid(lev)%dxu(0:ny+1,1:nx+1) = hlf*(grid(lev)%dx(0:ny+1,0:nx)+grid(lev)%dx(0:ny+1,1:nx+1))
          grid(lev)%dyv(1:ny+1,0:nx+1) = hlf*(grid(lev)%dy(0:ny,0:nx+1)+grid(lev)%dy(1:ny+1,0:nx+1))

       else               ! coarsen dx,dy
          ! (needed when directly discretizing on coarser grids)

          nxf =grid(lev-1)%nx
          nyf =grid(lev-1)%ny

          dxf => grid(lev-1)%dx
          dyf => grid(lev-1)%dy
          dxuf => grid(lev-1)%dxu
          dyvf => grid(lev-1)%dyv

          if (grid(lev)%gather == 1) then
             nxc= nx/grid(lev)%ngx
             nyc= ny/grid(lev)%ngy
             allocate(dxc(0:nyc+1,0:nxc+1))
             allocate(dyc(0:nyc+1,0:nxc+1))
             allocate(dxuc(0:nyc+1,0:nxc+1))
             allocate(dyvc(0:nyc+1,0:nxc+1))
          else
             nxc = nx
             nyc = ny
             dxc => grid(lev)%dx
             dyc => grid(lev)%dy
             dxuc => grid(lev)%dxu
             dyvc => grid(lev)%dyv
          endif

          dxc(1:nyc,1:nxc) = hlf      * ( & ! only interior points
               dxf(1:nyf  :2,1:nxf  :2) + &
               dxf(2:nyf+1:2,1:nxf  :2) + &
               dxf(1:nyf  :2,2:nxf+1:2) + &
               dxf(2:nyf+1:2,2:nxf+1:2) )

          dyc(1:nyc,1:nxc) = hlf      * ( &
               dyf(1:nyf  :2,1:nxf  :2) + &
               dyf(2:nyf+1:2,1:nxf  :2) + &
               dyf(1:nyf  :2,2:nxf+1:2) + &
               dyf(2:nyf+1:2,2:nxf+1:2) )

          dxuc(1:nyc,1:nxc) = hlf      * ( & ! only interior points
               dxuf(1:nyf  :2,1:nxf  :2) + &
               dxuf(2:nyf+1:2,1:nxf  :2) + &
               dxuf(1:nyf  :2,2:nxf+1:2) + &
               dxuf(2:nyf+1:2,2:nxf+1:2) )

          dyvc(1:nyc,1:nxc) = hlf      * ( &
               dyvf(1:nyf  :2,1:nxf  :2) + &
               dyvf(2:nyf+1:2,1:nxf  :2) + &
               dyvf(1:nyf  :2,2:nxf+1:2) + &
               dyvf(2:nyf+1:2,2:nxf+1:2) )

          if (grid(lev)%gather == 1) then
             call gather(lev,dxc,grid(lev)%dx)
             call gather(lev,dyc,grid(lev)%dy)
             call gather(lev,dxuc,grid(lev)%dxu)
             call gather(lev,dyvc,grid(lev)%dyv)
             deallocate(dxc)
             deallocate(dyc)
             deallocate(dxuc)
             deallocate(dyvc)
          endif

       endif

       call fill_halo(lev,grid(lev)%dx)
       call fill_halo(lev,grid(lev)%dy)
       call fill_halo(lev,grid(lev)%dxu)
       call fill_halo(lev,grid(lev)%dyv)

    enddo

  end subroutine set_horiz_grids

end module mg_horiz_grids
