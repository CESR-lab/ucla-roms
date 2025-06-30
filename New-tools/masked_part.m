! ----------------------------------------------------------------------
      subroutine partition_mask(npx,npy)   ![
      implicit none

      ! input/outputs
      integer,           intent(in) :: npx,npy

      ! local
      real,dimension(:,:),allocatable :: mask
      real    :: msk_mx
      integer :: nx,ny
      integer :: ierr,ncid,dimid
      integer :: surplus_x,surplus_y
      integer :: loc_x,loc_y
      integer :: npi,npj
      integer :: count


      ierr = nf90_open(trim(adjustl(grdfile)), nf90_nowrite, ncid)
       if (ierr/=0) call handle_ierr(ierr,'opening: ',trim(adjustl(grdfile)))

      ierr = nf90_inq_dimid(ncid,'xi_rho',dimid)
       if (ierr/=0) call handle_ierr(ierr,'getting dimid: ','xi_rho')
      ierr = nf90_inquire_dimension(ncid, dimid, len=gnx)
       if (ierr/=0) call handle_ierr(ierr,'getting dimension: ','xi_rho')

      ierr = nf90_inq_dimid(ncid,'eta_rho',dimid)
       if (ierr/=0) call handle_ierr(ierr,'getting dimid: ','eta_rho')
      ierr = nf90_inquire_dimension(ncid, dimid, len=gny)
       if (ierr/=0) call handle_ierr(ierr,'getting dimension: ','eta_rho')

      gnx = gnx-2
      gny = gny-2

      allocate(mask(0:gnx+1,0:gny+1))
      call ncread(ncid, 'mask_rho',mask)
      ierr = nf90_close(ncid)
      mask = 1

      ! Get interior subdomain size
      nx = ceiling(1.0*gnx/npx)
      ny = ceiling(1.0*gny/npy)

      surplus_x = nx*npx - gnx
      surplus_y = ny*npy - gny

      ! Array's to store the location and size of each subdomain
      allocate(iloc(npx,3)) ! lnx = iloc(:,1), is = iloc(:,2), ie = iloc(:,3)
      allocate(jloc(npy,3)) ! lny = jloc(:,1), js = jloc(:,2), je = jloc(:,3)

      loc_x = 1
      do npi = 1,npx
        if (npi==1) then ! left-most subdomain
          iloc(npi,1) = nx - surplus_x/2
        elseif (npi==npx) then ! East-most subdomain
          iloc(npi,1) = nx - (surplus_x+1)/2
        else
          iloc(npi,1) = nx
        endif
        iloc(npi,2) = loc_x
        iloc(npi,3) = loc_x + iloc(npi,1) -1
        loc_x = loc_x + iloc(npi,1)
      enddo

      loc_y = 1
      do npj = 1,npy
        if (npj==1) then ! left-most subdomain
          jloc(npj,1) = ny - surplus_y/2
        elseif (npj==npy) then ! North-most subdomain
          jloc(npj,1) = ny - (surplus_y+1)/2
        else
          jloc(npj,1) = ny
        endif
        jloc(npj,2) = loc_y
        jloc(npj,3) = loc_y + jloc(npj,1) -1
        loc_y = loc_y + jloc(npj,1)
      enddo

      allocate(node_map(0:npy+1,0:npx+1))
      node_map = -98
      allocate(npi_c(npx*npy))
      allocate(npj_c(npx*npy))

      count = 0
      do npj = 1,npy
        do npi = 1,npx
          msk_mx = maxval(mask(iloc(npi,2):iloc(npi,3),jloc(npj,2):jloc(npj,3)))
          if (msk_mx /= 0) then
            count = count+1
            node_map(npj,npi) = count
            npi_c(count) = npi
            npj_c(count) = npj
          endif
        enddo
      enddo
      nparts = count
      print *,'nparts: ',nparts

      deallocate(mask)

      ! adjust iloc/jloc to reflect rho-vars in netcdf files
      ! first and last subdomains have buffers, so 1 point larger
      iloc(1,1) = iloc(1,1) + 1
      iloc(npx,1) = iloc(npx,1) + 1
      iloc(2:npx,2) = iloc(2:npx,2) + 1
      iloc(:,3) = iloc(:,3) + 1
      iloc(npx,3) = iloc(npx,3) + 1

      jloc(1,1) = jloc(1,1) + 1
      jloc(npy,1) = jloc(npy,1) + 1
      jloc(2:npy,2) = jloc(2:npy,2) + 1
      jloc(:,3) = jloc(:,3) + 1
      jloc(npy,3) = jloc(npy,3) + 1

      ! u-, and v- variables differ only for the first subdomain
      allocate(ilcu(npx,3))
      allocate(jlcv(npy,3))
      ilcu = iloc
      jlcv = jloc
      ilcu(1,1) = iloc(1,1)-1
      jlcv(1,1) = jloc(1,1)-1
      ilcu(2:npx,2) = ilcu(2:npx,2)-1
      ilcu(:,3) = ilcu(:,3)-1
      jlcv(2:npy,2) = jlcv(2:npy,2)-1
      jlcv(:,3) = jlcv(:,3)-1

      end subroutine partition_mask  !]
! ----------------------------------------------------------------------

 function optimal_partition(ncores,grdfile)

  mask = ncread(grdfile,'mask_rho');
  mask = mask(2:end-1,2:end-1);
  [gnx,gny] = size(mask)

  ar = gnx/gny                   % aspect ratio grid:
  ocean = sum(mask(:))/(gnx*gny) % percentage of grid that is not masked

  % first guess of x,y partition that gives roughly square subdomains
  % and should exceed the number of available cores by a bit.
  npx1 = ceiling(sqrt(ncores)*ar/ocean)
  npy1 = ceiling(sqrt(ncores)/ar/ocean)

  ! 16 guesses ok?
  npx0 = max(npx1-4,1)
  npy0 = max(npy1-4,1)

  i = 0
  for npy = npy0:npy1
    for npx = npx0:npx1
      call partition_mask(npx,npy)
      i = i+1
      active(i,1) =
      active(i,2) = npx
      active(i,3) = npy
    end
  end

      sort(active)


      end program optimal_part

