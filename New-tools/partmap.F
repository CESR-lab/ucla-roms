      program partmap

      use netcdf
      use nc_read_write
      use mpi
      use roms_part

      implicit none

!      integer :: npartx,nparty

      ! local
      real,dimension(:,:),allocatable :: mask,rflx,partition,inodes,jnodes
      integer,dimension(:,:),allocatable :: iloc,jloc
      integer,dimension(:,:),allocatable :: ilcu,jlcv

      real    :: msk_mx,rfx_mx
      integer :: nx,ny,npx,npy
      integer :: ierr,ncid,npid,dimid
      integer :: surplus_x,surplus_y
      integer :: loc_x,loc_y
      integer :: npi,npj
      integer :: count

      character(len=7)              :: vname1,vname2,vname3
      character(len=7),dimension(2) :: dname
      integer,dimension(2)          :: dsize
      integer                       :: varid


      ! Handle system arguments
      call cmd_args_map(npx,npy)

      ! read in data for global domains
      ierr = nf90_open(trim(adjustl(grdfile)), nf90_nowrite, ncid)
       if (ierr/=0) call handle_ierr(ierr,'opening: ',trim(adjustl(grdfile)))

      ! global dimensions
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

      ! mask rho
      allocate(mask(0:gnx+1,0:gny+1)); mask=1
      call ncread(ncid,'mask_rho',mask)

      ierr = nf90_close(ncid)
       if (ierr/=0) call handle_ierr(ierr,'closing: ',trim(adjustl(grdfile)))

      ! make new netcdf to visalize parition
      ierr = nf90_create('part.nc',nf90_netcdf4,npid)
       if (ierr/=0) call handle_ierr(ierr,'creating: ','part.nc')

      ! add variables
      dname=(/'xi_rho','eta_rho'/)
      dsize=(/gnx+2,gny+2/)
      vname1='node'
      varid=nccreate(npid,vname1,dname,dsize)
      vname2='inode'
      varid=nccreate(npid,vname2,dname,dsize)
      vname3='jnode'
      varid=nccreate(npid,vname3,dname,dsize)
      
      ! Set up subdomain indices
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

      allocate(partition(0:gnx,0:gny)); partition=-1
      allocate(inodes(0:gnx,0:gny));    inodes=-1
      allocate(jnodes(0:gnx,0:gny));    jnodes=-1

      count = 0
      do npj = 1,npy
        do npi = 1,npx
          msk_mx = maxval(mask(iloc(npi,2):iloc(npi,3),jloc(npj,2):jloc(npj,3)))
          !rfx_mx = maxval(rflx(iloc(npi,2):iloc(npi,3),jloc(npj,2):jloc(npj,3)))
          if (msk_mx > 0) then
            count = count+1
            partition(iloc(npi,2):iloc(npi,3),jloc(npj,2):jloc(npj,3))=count
            inodes(iloc(npi,2):iloc(npi,3),jloc(npj,2):jloc(npj,3))   =npi
            jnodes(iloc(npi,2):iloc(npi,3),jloc(npj,2):jloc(npj,3))   =npj
          endif
        enddo
      enddo

      ! add partition counts netcdf
      call ncwrite(npid,vname1,partition)
      call ncwrite(npid,vname2,inodes)
      call ncwrite(npid,vname3,jnodes)


      end program partmap
