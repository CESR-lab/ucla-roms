
      integer :: i,j
      integer :: ierr,ncid,varid
      real  :: riv_cells,riv_east,riv_west

      riv_west=xl*0.4 ! River west bank at 40% from west
      riv_east=xl*0.6 ! River west bank at 60% from west

      ! pm is constant for this case
      riv_cells = nint( (riv_east - riv_west)*pm(1,1)) !number of cells in this river
      do j=0,ny+1
        do i=0,nx+1
          if (xr(i,j)>riv_west .and. xr(i,j)<riv_east) then
            ! find 'coastline' masked cells
# ifdef MASKING
            if (rmask(i,j)==0 .and. rmask(i,j+1)==1) then
              rflx(i,j) = 1.0+1.0/riv_cells
            endif
# endif
          endif


        enddo
      enddo

      ierr=nf90_open(ana_grdname,nf90_write,ncid)
      varid = nccreate(ncid,'river_flux',(/dn_xr,dn_yr/),(/xi_rho,eta_rho/), nf90_double)
      ierr=nf90_put_att(ncid, varid,'long_name','River volume flux') 
!     ierr=nf90_close(ncid)
!     print *,'added river_flux',mynode
!     ierr=nf90_open(ana_grdname,nf90_write,ncid)
      call ncwrite(ncid,'river_flux', rflx(i0:i1,j0:j1))
      ierr=nf90_close(ncid)
