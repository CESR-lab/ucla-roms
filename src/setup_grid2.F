      subroutine setup_grid2

      ! Find minimum and maximum depth of model topography; minimum and
      ! maximum horizontal XI- and ETA-grid intervals around h-points;
      ! Check minimum and maximum values of Courant number for barotropic
      ! mode (based on surface gravity wave speed, barotropic time step
      ! and local grid intervals); combined crossection of all open
      ! boundaries.

      use param
      use roms_mpi
      use mpi
      use grid
      use comm_vars

      implicit none

!     use scalars

      ! Local 
      real,allocatable,dimension(:,:) :: fld
      real :: crssec
!     real :: area,volume
!     real :: hmax,grdmax,Cg_max,Cu_cor
!     real :: hmin,grdmin,Cg_min
      real,dimension(4) :: local,global
      integer :: ierr

      allocate(fld(1-bf:nx+bf,1-bf:ny+bf) )

      hmin = minval(h(1:nx,1:ny),mask=rmask(1:nx,1:ny)>0.5)
      hmax = maxval(h(1:nx,1:ny),mask=rmask(1:nx,1:ny)>0.5)

      fld = 1./sqrt(pm*pn)
      grdmin = minval(fld(1:nx,1:ny),mask=rmask(1:nx,1:ny)>0.5)
      grdmax = maxval(fld(1:nx,1:ny),mask=rmask(1:nx,1:ny)>0.5)

      fld = dtfast*sqrt(g*h*(pm*pm+pn*pn))
      Cg_min = minval(fld(1:nx,1:ny),mask=rmask(1:nx,1:ny)>0.5)
      Cg_max = maxval(fld(1:nx,1:ny),mask=rmask(1:nx,1:ny)>0.5)

      Cu_Cor = dt*maxval(abs(f(1:nx,1:ny)))

      ! Unperturbed volume of model grid assuming zero zeta
      fld = rmask/(pm*pn)
      area   = sum(fld(1:nx,1:ny))
      volume = sum(fld(1:nx,1:ny)*h(1:nx,1:ny))

      deallocate(fld)

      crssec = 0 
      if (inode.eq.0.and.obc_west) then
        crssec = crssec + 0.5*sum((h(0,1:ny)+h(1,1:ny))*dn_u(1,1:ny)*umask(1,1:ny))
      endif
      if (inode.eq.npx-1.and.obc_east) then
        crssec = crssec + 0.5*sum((h(nx,1:ny)+h(nx+1,1:ny))*dn_u(nx+1,1:ny)*umask(nx+1,1:ny))
      endif
      if (jnode.eq.0.and.obc_south) then
        crssec = crssec + 0.5*sum((h(1:nx,0)+h(1:nx,1))*dm_v(1:nx,1)*vmask(1:nx,1))
      endif
      if (jnode.eq.npy-1.and.obc_north) then
        crssec = crssec + 0.5*sum((h(1:nx,ny)+h(1:nx,ny+1))*dm_v(1:nx,ny+1)*vmask(1:nx,ny+1))
      endif

      local(1) = hmax
      local(2) = grdmax
      local(3) = Cg_max
      local(4) = Cu_cor
      call MPI_Allreduce(local, global, 4, MPI_DOUBLE_PRECISION,
     &                   MPI_MAX, ocean_grid_comm, ierr)
      hmax   = global(1)
      grdmax = global(2)
      Cg_max = global(3)
      Cu_Cor = global(4)

      local(1) = hmin
      local(2) = grdmin
      local(3) = Cg_min
      call MPI_Allreduce(local, global, 3, MPI_DOUBLE_PRECISION,
     &                   MPI_MIN, ocean_grid_comm, ierr)
      hmin   = global(1)
      grdmin = global(2)
      Cg_min = global(3)

      local(1) = area
      local(2) = volume
      local(3) = crssec
      call MPI_Allreduce(local, global, 3, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, ocean_grid_comm, ierr)
      area   = global(1)
      volume = global(2)
      crssec = global(3)

      if (mynode==0) then
        write(*,'(1x,A,F12.6,3x,A,ES14.7,5x,A,ES23.16)')
     &    'hmin =' ,hmin,  'grdmin =', grdmin, 'area =',   area
        write(*,'(1x,A,F12.6,3x,A,ES14.7,3x,A,ES23.16)')
     &    'hmax =' ,hmax,  'grdmax =', grdmax, 'volume =', volume
        write(*,'(43x,A,ES23.16)')        'open_cross =', crssec
        write(*,'(1x,A,F10.7,3x,A,F10.7,3x,A,F10.7)')
     &       'Cg_max =',Cg_max, 'Cg_min =',Cg_min, 'Cu_Cor =',Cu_Cor
      endif

      end
