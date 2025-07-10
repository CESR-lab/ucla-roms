      module roms_part
      
      ! ROMS functionality to facilitate partioning and joining of files

      use netcdf
      use mpi
      use nc_read_write

      private
       
      ! mpi variables
      integer,public :: part_comm,mynode,nnodes  ! mpi vars

      ! global dimensins
      integer,public :: gnx,gny

      ! Files
      integer,public :: nfiles ! number of src files
      character(len=128),public :: grdfile,srcfile
      character(len=128),public :: fullsrcfile,joinfile,partfile
      character(len=128),public :: gpath,spath
      character(len=128),public,allocatable, dimension(:) :: srcfiles
      character(len=128) :: fname

      ! Only used by mynode=0
      integer,dimension(:,:),allocatable :: iloc,jloc
      integer,dimension(:,:),allocatable :: ilcu,jlcv
      integer,dimension(:,:),allocatable :: node_map
      integer,dimension(:)  ,allocatable :: npi_c,npj_c
      integer,dimension(:)  ,allocatable :: boundary_nodes

      ! Dimension names for in netcdf files
      character(len=7),public,parameter :: dn_xr = 'xi_rho'
      character(len=7),public,parameter :: dn_yr = 'eta_rho'
      character(len=7),public,parameter :: dn_xu = 'xi_u'
      character(len=7),public,parameter :: dn_yv = 'eta_v'
      character(len=7),public,parameter :: dn_zr = 's_rho'
      character(len=7),public,parameter :: dn_zw = 's_w'
      character(len=7),public,parameter :: dn_tm = 'time'
      character(len=10),public :: tdim

      ! module local variables
      integer,public       :: nparts
      
      type,public :: ncstructure
        integer                                     :: nvars     ! number of variables
        integer,dimension(3)                        :: gsize     ! dimension size of the global grid (rho)
        character(len=30),dimension(:),allocatable  :: varname   ! name, if u/v or rho
        character(len=1),dimension(:),allocatable   :: vartype   ! 'f/r/u f/r/v f/r/w f/t'
        character(len=2),dimension(:),allocatable   :: bryvar    ! is this a boundary variable
        logical,dimension(:),allocatable            :: varpart   ! partitionable or not
        integer,dimension(:),allocatable            :: xtype     ! variable type
        integer,dimension(:),allocatable            :: ndims  
        integer,dimension(:),allocatable            :: natts
        type(dimstruct),dimension(:),allocatable    :: dims
        type(attstruct),dimension(:),allocatable    :: atts
      end type ncstructure

      type,public :: dimstruct
        character(len=10),dimension(:),allocatable ::dimnames
        integer,         dimension(:),allocatable ::dimsizes
      end type dimstruct

      type,public :: attstruct
        character(len=20),dimension(:),allocatable  ::attnames
        character(len=200),dimension(:),allocatable ::attvalues
        real,dimension(:),allocatable               ::attparams
      end type attstruct
      
      public :: cmd_args_optimal
      public :: cmd_args_map
      public :: cmd_args_partit
      public :: partition_mask
      public :: check_mask
      public :: make_partial_files
      public :: write_partial_files

      contains
! ----------------------------------------------------------------------
      subroutine cmd_args_optimal(cores,grdfile)  ![
      implicit none

!     import/export
      real,intent(out)            :: cores 
      character(len=128),intent(out) :: grdfile 

!     local
      integer            :: nargs
      character(len=128) :: arg 

      nargs =  command_argument_count()

      if (nargs < 2) call show_use_optimal

!     available cpus/cores
      call get_command_argument(1,arg);
      read(arg, *) cores
!     grdfile
      call get_command_argument(2,arg);
      grdfile = arg

      end subroutine cmd_args_optimal !]
! ----------------------------------------------------------------------
      subroutine cmd_args_map(npartx,nparty)  ![

!     import/export
      integer,intent(out) :: npartx,nparty

      integer :: nargs,i
      character(len=128) :: arg

      nargs =  command_argument_count()

      if (nargs < 3) then
        print*,'Command Template:  part_map npartx nparty grdfile.nc'
      endif
      
      ! first 2 should be number of subdomains in x- or y- direction
      call get_command_argument(1,arg);
      read(arg,*) npartx
      call get_command_argument(2,arg);
      read(arg,*) nparty
      call get_command_argument(3,arg)
      read(arg,*) grdfile

      end subroutine cmd_args_map  !]
! ----------------------------------------------------------------------
      subroutine cmd_args_partit(npartx,nparty)  ![
      implicit none

!     import/export
      integer,intent(out) :: npartx,nparty

      integer :: nargs,i
      character(len=128) :: arg

      nargs =  command_argument_count()

      if (nargs < 4) call show_use_partit
      
      ! first 2 should be number of subdomains in x- or y- direction
      call get_command_argument(1,arg);
      read(arg,*) npartx
      call get_command_argument(2,arg);
      read(arg,*) nparty
      print *,nargs,npartx,nparty

      call get_command_argument(3,arg)
      grdfile = trim(adjustl(arg))
!     read(arg,*) grdfile
!     print *,'gridfile: ',grdfile
      
      nfiles = nargs - 3
      allocate(srcfiles(nfiles))
      do i = 4,nargs
        call get_command_argument(i,arg)
        srcfiles(i-3) = trim(adjustl(arg))
!       read(arg,*) srcfiles(i-3)
      enddo
      print *,grdfile,srcfiles(1)

      end subroutine cmd_args_partit !]
! ----------------------------------------------------------------------
      subroutine show_use_map  ![
      implicit none

       print *, 'Usage of partit should be: ',
     &       'partmap npx npy grdfile.nc'
     
      stop

      end subroutine show_use_map  !] 
! ----------------------------------------------------------------------
      subroutine show_use_partit  ![
      implicit none

      print *, 'Usage of partit should be: ',
     &       'partit npx npy grdfile.nc, file1.nc...'
     
      stop

      end subroutine show_use_partit !]
! ----------------------------------------------------------------------
      subroutine show_use_optimal ![
      implicit none

      print *, 'Usage of optimal should be: ',
     &         'optimal_part cores grdfile.nc'
       
      stop

      end subroutine show_use_optimal !]
! ----------------------------------------------------------------------
      subroutine check_mask(npx,npy,gnx,gny,mask)   ![
      implicit none

      ! input/outputs
      integer, intent(in)                        :: npx,npy,gnx,gny
      real,dimension(0:gnx+1,0:gny+1),intent(in) :: mask

      ! local
      real,dimension(:,:),allocatable :: rflx
      real    :: msk_mx,rfx_mx
      integer :: nx,ny
      integer :: surplus_x,surplus_y
      integer :: loc_x,loc_y
      integer :: npi,npj
      integer :: count
      
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
          rfx_mx = 0
!         rfx_mx = maxval(rflx(iloc(npi,2):iloc(npi,3),jloc(npj,2):jloc(npj,3)))
          if (msk_mx > 0.or.rfx_mx>0) then
            count = count+1
            node_map(npj,npi) = count
            npi_c(count) = npi
            npj_c(count) = npj
          endif
        enddo
      enddo
      nparts = count

      deallocate(iloc,jloc,node_map,npi_c,npj_c)

      end subroutine check_mask !]
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
       print *, grdfile
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
!      mask = 1

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
      allocate(boundary_nodes(npx*npy))
      boundary_nodes = 0

      count = 0
      do npj = 1,npy
        do npi = 1,npx
          msk_mx = maxval(mask(iloc(npi,2):iloc(npi,3),jloc(npj,2):jloc(npj,3)))
          if (msk_mx > 0) then
            count = count+1
            node_map(npj,npi) = count
            npi_c(count) = npi
            npj_c(count) = npj
          endif
          ! record if at boundary
          if (((npi==1).or.(npi==npx)).or.((npj==1).or.(npj==npy))) then
            boundary_nodes(count)=1
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
      subroutine make_partial_files  ![
      implicit none

      ! local
      character(len=128) :: fname
      character(len=128) :: att_name
      integer,dimension(4) :: ibuff
      integer,dimension(8) :: loc_neigh
      integer :: iatt,ngatt,part
      integer :: i,j,varid
      integer :: dimids(4),ndims
      integer :: ivar,nvars,xtype
      character(len=20) :: vname
      type(ncstructure)  :: ncinfo
      character(len=10),dimension(5)  ::dname
      integer         ,dimension(5)  ::dsize
      character(len=20),dimension(10) ::attname
      character(len=200),dimension(10)::attvalue
      real,dimension(10)              ::attparam
      integer :: ierr, ncid,npid
      integer :: natts,idm,vid,did
      integer,dimension(2) :: subdompos

      ! open up reference netcdf for global attributes
      ierr = nf90_open(trim(srcfile), nf90_nowrite, ncid)
      if (ierr/=nf90_noerr) 
     &  call handle_ierr(ierr,'opening file: ',trim(srcfile) )

      call get_nc_info(ncid, ncinfo) 

      do part = 1,nparts
        i = npi_c(part)
        j = npj_c(part)
    
        fname = srcfile
        call insert_nodes(fname,part-1,nparts)

       ierr=nf90_create(trim(fname),nf90_netcdf4,npid)
!        ierr=nf90_create(trim(fname),nf90_classic_model,npid)
        if (ierr/=nf90_noerr) 
     &    call handle_ierr(ierr,'creating file: ',trim(fname) )

        ! global rho-dimensions
        ierr=nf90_put_att(npid, nf90_global, 'globaldims', ncinfo%gsize)

        ! subdompos: i,j position of sub-domain
        subdompos = (/ i-1, j-1 /) ! zero-based!
        ierr=nf90_put_att(npid, nf90_global, 'subdompos', subdompos)

        ! partition: node, total nodes, llc x, llc y
        ibuff = (/part-1, nparts, iloc(i,2), jloc(j,2) /)
        ierr=nf90_put_att(npid, nf90_global, 'partition', ibuff)

        ! starting indices: i0, j0, i0u, j0v
        ibuff = (/iloc(i,2),jloc(j,2),ilcu(i,2),jlcv(j,2) /)
        ierr=nf90_put_att(npid, nf90_global, 'start_idx', ibuff)

        ! ending indices: i1, j1, i1u, j1v
        ibuff = (/iloc(i,3),jloc(j,3),ilcu(i,3),jlcv(j,3) /)
        ierr=nf90_put_att(npid, nf90_global, 'end_idx', ibuff)

        ! neighbors: N,NE,E,SE,S,SW,W,NW (READ IN THAT ORDER)
        loc_neigh = (/ node_map(j+1,i), node_map(j+1,i+1),
     &                 node_map(j,i+1), node_map(j-1,i+1),
     &                 node_map(j-1,i), node_map(j-1,i-1),
     &                 node_map(j,i-1), node_map(j+1,i-1) /)
        ierr=nf90_put_att(npid, nf90_global, 'neighbors',loc_neigh-1)
         
        ! No I/O is performed when NF90_INQUIRE is called, since the required information is 
        ! available in memory for each open netCDF dataset.

        ! read all global attributes from srcfile and write them
        ierr = nf90_inquire(ncid, nAttributes=ngatt)
        do iatt = 1,ngatt
          ierr = nf90_inq_attname(ncid,nf90_global,iatt,att_name)
          ! if present, these attributes would mess up the writes above
          if (index(att_name,'start_idx')) cycle
          if (index(att_name,'end_idx')) cycle
          if (index(att_name,'neighbors')) cycle
          if (index(att_name,'subdompos')) cycle
          if (index(att_name,'globaldims')) cycle
          ierr = nf90_copy_att(ncid,nf90_global,att_name,npid,nf90_global)
        enddo 

        do ivar = 1,ncinfo%nvars
          vname   = ncinfo%varname(ivar)
          if (index(vname,'coarse')>0) then
            if (part==1) print *,'Skipping coarse variable: ', vname
            cycle
          endif
          xtype   = ncinfo%xtype(ivar)
          ndims   = ncinfo%ndims(ivar)
          dsize(1:ndims) = ncinfo%dims(ivar)%dimsizes
          dname(1:ndims) = ncinfo%dims(ivar)%dimnames
          natts   = ncinfo%natts(ivar)
          attname(1:natts) = ncinfo%atts(ivar)%attnames
          attvalue(1:natts)= ncinfo%atts(ivar)%attvalues
          attparam(1:natts)= ncinfo%atts(ivar)%attparams

          do idm =1,ndims
            if (index(dname(idm),'xi_rho') > 0)  dsize(idm)=iloc(i,1)
            if (index(dname(idm),'xi_u') > 0)    dsize(idm)=ilcu(i,1)
            if (index(dname(idm),'eta_rho') > 0) dsize(idm)=jloc(j,1)
            if (index(dname(idm),'eta_v') > 0)   dsize(idm)=jlcv(j,1)
            if (index(dname(idm),'time') > 0)    dsize(idm)= 0 
          enddo

          if (ndims==0) then
            dname='N/A'
          endif

          if (ncinfo%vartype(ivar).eq.'s') then
            continue
          elseif (ndims==0) then
            ierr=nf90_def_var(npid,vname,xtype,varid)
            do iatt=1,natts
              ierr = nf90_put_att(npid,varid,attname(iatt),attvalue(iatt))
            enddo
          else
            vid = nccreate(npid,vname,dname(1:ndims),dsize(1:ndims),xtype)
            do iatt=1,natts
              ! input attribute based on type
              if (attvalue(iatt)=='n/a') then  ! parameter (type(real))
                ierr = nf90_put_att(npid,vid,attname(iatt),attparam(iatt))
              elseif (attparam(iatt)==0) then  ! string    (type(char))
                ierr = nf90_put_att(npid,vid,attname(iatt),attvalue(iatt))
              endif
            enddo
          endif

        enddo


        ierr = nf90_close(npid)

      enddo 
      
      ! save specific name of time dimension
      ndims = ncinfo%ndims(1)
      dname(1:ndims) = ncinfo%dims(1)%dimnames
      do idm =1,ndims
        if (index(dname(idm),'tim')>0) then
          tdim = dname(idm)
          exit
        endif
      enddo

      ierr = nf90_close(ncid)

      end subroutine make_partial_files !]
! ----------------------------------------------------------------------
      subroutine write_partial_files  ![
      ! get info from first partial file in range
      ! read block of data for each partitionable variable
      ! and write to partial files, for all records in srcfile
      implicit none

      ! local
      character(len=128) :: fname
      integer :: ierr, ncid,npid,dimid
      integer,dimension(4) :: llidx
      integer :: nz
      integer,dimension(:),allocatable :: i0,j0
      integer,dimension(:),allocatable :: i1,j1
      integer,dimension(:),allocatable :: i0u,i1u
      integer,dimension(:),allocatable :: j0v,j1v
      integer :: li0,li1,lj0,lj1,li0u,li1u,lj0v,lj1v
      real,dimension(:,:,:),allocatable :: var3d
      real,dimension(:,:),allocatable   :: var2d,t2d
      real,dimension(:),allocatable     :: var0d
      character(len=1) :: charvar
      character(len=20) :: vname
      character(len=1)  :: bry
      type(ncstructure)  :: ncinfo
      integer,dimension(:),allocatable :: is,ie,js,je,
     &                                    bs,be         ! boundary


      integer :: part0,part1,myparts,part
      integer :: varid,ivar
      integer :: record,nrecords
      integer,dimension(4) :: startr,startw,startwr
      real,dimension(:),allocatable :: var1d
      integer :: sz1,sz2,sz3,ndims
      integer :: RSS

      ! profiling
      integer,dimension(:),allocatable :: ncparts
      real :: startrec,endrec,
     &        startwri,endwri,
     &        startvar,endvar

      myparts = ceiling(nparts*1.0/nnodes)
      part0 = myparts*mynode
      part1 = min(myparts*(mynode+1) -1,nparts-1)

      ! get info about variables from first partial file in range

      ! Read partial files in range
      allocate(i0(part0:part1)); allocate(i0u(part0:part1))
      allocate(i1(part0:part1)); allocate(i1u(part0:part1))
      allocate(j0(part0:part1)); allocate(j0v(part0:part1))
      allocate(j1(part0:part1)); allocate(j1v(part0:part1))

      allocate(is(part0:part1))
      allocate(js(part0:part1))
      allocate(ie(part0:part1))
      allocate(je(part0:part1))
      allocate(bs(part0:part1))
      allocate(be(part0:part1))

      allocate(ncparts(part0:part1))

      do part = part0,part1  ! zero based ![
    
        fname = srcfile
        call insert_nodes(fname,part,nparts)

        ierr = nf90_open(trim(adjustl(fname)),nf90_write,npid)
        if (ierr/=nf90_noerr) 
     &    call handle_ierr(ierr,'opening file: ',trim(fname) )

        if (part.eq.part0) call get_nc_info(npid, ncinfo) 

        ierr = nf90_get_att(npid,nf90_global,'start_idx',llidx)
        i0(part) = llidx(1)
        j0(part) = llidx(2)
        i0u(part)= llidx(3)
        j0v(part)= llidx(4)

        ierr = nf90_get_att(npid,nf90_global,'end_idx',llidx)
        i1(part) = llidx(1)
        j1(part) = llidx(2)
        i1u(part)= llidx(3)
        j1v(part)= llidx(4)

!        ierr = nf90_close(npid)
        ncparts(part) = npid

      enddo ! part=part0,part1  !]

      li0 = minval(i0)
      li1 = maxval(i1)
      lj0 = minval(j0)
      lj1 = maxval(j1)
      li0u = minval(i0u)
      li1u = maxval(i1u)
      lj0v = minval(j0v)
      lj1v = maxval(j1v)

      call system_mem_usage(RSS)

      ierr = nf90_open(trim(srcfile),nf90_nowrite,ncid)
      if (ierr/=0) call handle_ierr(ierr,'opening: ',trim(srcfile))

      ierr = nf90_inq_dimid(ncid,'time',dimid)
      if (ierr /= nf90_noerr) then
        ierr = nf90_inq_dimid(ncid,'rec_time',dimid)    ! backwards compatibility for bry files
        if (ierr /= nf90_noerr) then
          ierr = nf90_inq_dimid(ncid,tdim,dimid)        ! backwards compatibility for sample input data
          if (ierr /= nf90_noerr) then
            print *,'no time dim'
            nrecords=1
          else
            ierr = nf90_inquire_dimension(ncid,dimid,len=nrecords)
            if (ierr /= nf90_noerr) call handle_ierr(ierr)
          endif
        else
            ierr = nf90_inquire_dimension(ncid,dimid,len=nrecords)
            if (ierr /= nf90_noerr) call handle_ierr(ierr)
        endif
      else
        ierr = nf90_inquire_dimension(ncid,dimid,len=nrecords)
        if (ierr /= nf90_noerr) call handle_ierr(ierr)
      endif
     
      print *,'nrecords: ',nrecords

      do record = 1,nrecords

        if (mynode==0) print *,'Processing record: ',record

      do ivar=1,ncinfo%nvars

        vname = ncinfo%varname(ivar)
        if (mynode==0) print *,'Processing variable: ',vname
!       if (mynode==0) print *,'Partitionable? ',ncinfo%varpart(ivar)

        if (ncinfo%xtype(ivar).eq.nf90_char) goto 99

        bry = ncinfo%bryvar(ivar)

        ndims = ncinfo%ndims(ivar)
        if (vname=='time_step') cycle
        if (index(vname,'coarse')>0) cycle

        if (ncinfo%varpart(ivar)) then
          select case(ndims)
          case(1)
            ! save for later
          case(2)
            ! ----- BOUNDARIES ONLY ----- ![
            if (index(ncinfo%dims(ivar)%dimnames(2),'tim')>0) then
              if ((bry=='w').or.(bry=='e')) then
                select case(ncinfo%vartype(ivar))
                case('r')
                  allocate(var2d(lj0:lj1,1))
                  startr = (/lj0,record,1,1/)
                  bs = j0; be = j1
                case('v')
                  allocate(var2d(lj0v:lj1v,1))
                  startr = (/lj0v,record,1,1/)
                  bs = j0v; be = j1v
                end select
              elseif ((bry=='s').or.(bry=='n')) then
                select case(ncinfo%vartype(ivar))
                case('r')
                  allocate(var2d(li0:li1,1))
                  startr = (/li0,record,1,1/)
                  bs = i0; be = i1
                case('u')
                  allocate(var2d(li0u:li1u,1))
                  startr = (/li0u,record,1,1/)
                  bs = i0u; be = i1u
                end select
              endif
              startw = (/1,record,1,1/)
            ! --------------------------- !]
            else
              if (record>1) continue
              select case(ncinfo%vartype(ivar))
              case('r')
                allocate(var2d(li0:li1,lj0:lj1))
                startr = (/li0,lj0,record,1/)
                is = i0; ie = i1
                js = j0; je = j1
              case('u')
                allocate(var2d(li0u:li1u,lj0:lj1))
                startr = (/li0u,lj0,record,1/)
                is = i0u; ie = i1u
                js = j0; je = j1
              case('v')
                allocate(var2d(li0:li1,lj0v:lj1v))
                startr = (/li0,lj0v,record,1/)
                is = i0; ie = i1
                js = j0v; je = j1v
              case('j')
                allocate(var2d(li0u:li1u,lj0v:lj1v))
                startr = (/li0u,lj0v,record,1/)
                is = i0u; ie = i1u
                js = j0v; je = j1v
              end select
              startwr = (/i0,js,record,1/)
              startw  = (/1,1,record,1/)
            endif
          case(3)
            ! ----- BOUNDARIES ONLY ----- ![
            if ((index(ncinfo%dims(ivar)%dimnames(3),'tim')>0).and.(bry/='i')) then
              if ((bry=='w').or.(bry=='e')) then
                nz = ncinfo%dims(ivar)%dimsizes(2)
                select case(ncinfo%vartype(ivar))
                case('r')
                  allocate(var3d(lj0:lj1,nz,1))
                  startr = (/lj0,1,record,1/)
                  bs = j0; be = j1
                case('v')
                  allocate(var3d(lj0:lj1v,nz,1))
                  startr = (/lj0v,1,record,1/)
                  bs = j0v; be = j1v
                end select
              elseif ((bry=='s').or.(bry=='n')) then
                select case(ncinfo%vartype(ivar))
                case('r')
                  allocate(var3d(li0:li1,nz,1))
                  startr = (/li0,1,record,1/)
                  bs = i0; be = i1
                case('u')
                  allocate(var3d(li0:li1u,nz,1))
                  startr = (/li0u,1,record,1/)
                  bs = i0u; be = i1u
                end select
              endif
              startw = (/1,1,record,1/)
            ! --------------------------- !]
            elseif (index(ncinfo%dims(ivar)%dimnames(3),'time')>0) then
              select case(ncinfo%vartype(ivar))
              case('r')
                allocate(var3d(li0:li1,lj0:lj1,1))
                startr = (/li0,lj0,record,1/)
                is = i0; ie = i1
                js = j0; je = j1
              case('u')
                allocate(var3d(li0u:li1u,lj0:lj1,1))
                startr = (/li0u,lj0,record,1/)
                is = i0u; ie = i1u
                js = j0; je = j1
              case('v')
                allocate(var3d(li0:li1,lj0v:lj1v,1))
                startr = (/li0,lj0v,record,1/)
                is = i0; ie = i1
                js = j0v; je = j1v
              end select
              startw = (/1,1,record,1/)
            else
              if (record>1) continue
              nz = ncinfo%dims(ivar)%dimsizes(3)
              select case(ncinfo%vartype(ivar))
              case('r')
                allocate(var3d(li0:li1,lj0:lj1,nz))
                startr = (/li0,lj0,1,1/)
                is = i0; ie = i1
                js = j0; je = j1
              case('u')
                allocate(var3d(li0u:li1u,lj0:lj1,nz))
                startr = (/li0u,lj0,1,1/)
                is = i0u; ie = i1u
                js = j0; je = j1
              case('v')
                allocate(var3d(li0:li1,lj0v:lj1v,nz))
                startr = (/li0,lj0v,1,1/)
                is = i0; ie = i1
                js = j0v; je = j1v
              end select
              startw = (/1,1,1,1/)
            ! save for later (not present in pachug_rst.nc)
            endif
          case(4)
            nz = ncinfo%dims(ivar)%dimsizes(3)
            select case(ncinfo%vartype(ivar))
            case('r')
              allocate(var3d(li0:li1,lj0:lj1,nz))
              startr = (/li0,lj0,1,record/)
              is = i0; ie = i1
              js = j0; je = j1
            case('u')
              allocate(var3d(li0u:li1u,lj0:lj1,nz))
              startr = (/li0u,lj0,1,record/)
              is = i0u; ie = i1u
              js = j0; je = j1
            case('v')
              allocate(var3d(li0:li1,lj0v:lj1v,nz))
              startr = (/li0,lj0v,1,record/)
              is = i0; ie = i1
              js = j0v; je = j1v
            end select
            startw = (/1,1,1,record/)
          end select
        else   ! not partitionable
          select case(ndims)
          case(0) ! needed to handle single value variables
            goto 98
          case(1)
            if (index(ncinfo%dims(ivar)%dimnames(1),'time')>0) then
              startr = (/record,1,1,1/)
              startw = (/record,1,1,1/)
              allocate(var3d(1,1,1))
              is = 1; ie = 1
              js = 1; je = 1
            else
              if (record>1) continue
              startr = (/1,1,1,1/)
              startw = (/1,1,1,1/)
              sz1 = ncinfo%dims(ivar)%dimsizes(1)
              allocate(var3d(sz1,1,1))
              is = 1; ie = sz1
              js = 1; je = 1
            endif
          case(2)
            if (index(ncinfo%dims(ivar)%dimnames(2),'time')>0) then
              startr = (/1,record,1,1/)
              startw = (/1,record,1,1/)
              sz1 = ncinfo%dims(ivar)%dimsizes(1)
              allocate(var2d(sz1,1))
              is = 1; ie = sz1
              js = 1; je = 1
            else
              if (record>1) continue
              startr = (/1,1,1,1/)
              startw = (/1,1,1,1/)
              sz1 = ncinfo%dims(ivar)%dimsizes(1)
              sz2 = ncinfo%dims(ivar)%dimsizes(2)
              allocate(var2d(sz1,sz2))
              is = 1; ie = sz1
              js = 1; je = sz2
            endif
          case(3)
            if (index(ncinfo%dims(ivar)%dimnames(3),'time')>0) then
              startr = (/1,1,record,1/)
              startw = (/1,1,record,1/)
              sz1 = ncinfo%dims(ivar)%dimsizes(1)
              sz2 = ncinfo%dims(ivar)%dimsizes(2)
              allocate(var3d(sz1,sz2,1))
              is = 1; ie = sz1
              js = 1; je = sz2
            else
              startr = (/1,1,1,1/)
              startw = (/1,1,1,1/)
              sz1 = ncinfo%dims(ivar)%dimsizes(1)
              sz2 = ncinfo%dims(ivar)%dimsizes(2)
              sz3 = ncinfo%dims(ivar)%dimsizes(3)
              allocate(var3d(sz1,sz2,sz3))
              is = 1; ie = sz1
              js = 1; je = sz2
            endif
          end select
        endif

        ! variables that need to be saved once
!        if ((index(ncinfo%dims(ivar)%dimnames(1),'nriver')>0)
!     &  .or.(index(ncinfo%dims(ivar)%dimnames(1),'npipe')>0)) then
!        if (record==1) then
!          if (ndims==2) then
!            call ncread(ncid,vname,var2d)
!            do part = part0,part1
!              call ncwrite(ncparts(part),vname,var2d)
!            enddo
!          endif
!          goto 99
!        else
!          cycle
!        endif
!        endif

        if (ndims==1) call ncread(ncid,vname,var3d,startr)
        if (ndims==2) call ncread(ncid,vname,var2d,startr)
        if (ndims==3) call ncread(ncid,vname,var3d,startr)
        if (ndims==4) call ncread(ncid,vname,var3d,startr)

        do part = part0,part1  ! zero based
          fname = srcfile

          npid = ncparts(part)

          if ((bry=='w').or.(bry=='e').or.(bry=='s').or.(bry=='n')) then  ! check if boundary var
            if (boundary_nodes(part+1)==1) then                           ! only write at boundary
              if (ndims==2) then
                call ncwrite(ncparts(part),vname,var2d(bs(part):be(part),:),startw)
              elseif (ndims==3) then
                call ncwrite(ncparts(part),vname,var3d(bs(part):be(part),:,:),startw)
              endif
            else
              cycle
            endif  ! <-- not at boundary node
          else     ! <-- not boundary var
            if (ndims==1) then
              call ncwrite(ncparts(part),vname,var3d(is(part):ie(part),js(part):je(part),:),startw)
            elseif (ndims==2) then
              call ncwrite(ncparts(part),vname,var2d(is(part):ie(part),js(part):je(part)),startw)
            elseif (ndims==3) then
              call ncwrite(ncparts(part),vname,var3d(is(part):ie(part),js(part):je(part),:),startw)
            elseif (ndims==4) then
              call ncwrite(ncparts(part),vname,var3d(is(part):ie(part),js(part):je(part),:),startw)
            endif
          endif

        enddo

        if (allocated(var2d)) deallocate(var2d)
        if (allocated(var3d)) deallocate(var3d)

        ! single value (double/real) variables
  98    if ((ncinfo%xtype(ivar).eq.nf90_double).and.(ndims<=1)) then
          allocate(var0d(1))
          ierr = nf90_inq_varid(ncid,vname,varid)
          ierr = nf90_get_var(ncid,varid,var0d)
          do part = part0,part1
            npid = ncparts(part)
            ierr = nf90_put_var(npid,varid,var0d)
             if (ierr/=nf90_noerr) call handle_ierr(ierr)
          enddo
          if (allocated(var0d)) deallocate(var0d)
        endif

        ! handle char variables
  99    if (ncinfo%xtype(ivar).eq.nf90_char) then        
          ierr = nf90_inq_varid(ncid,vname,varid)
          ierr = nf90_get_var(ncid,varid,charvar)
          do part = part0,part1
            npid = ncparts(part)
            ierr = nf90_put_var(npid,varid,charvar)
             if (ierr/=nf90_noerr) call handle_ierr(ierr)
          enddo
        endif
        
      enddo  ! ivars

      enddo  ! records

      ierr = nf90_close(ncid)

      deallocate(i0); deallocate(i0u);
      deallocate(i1); deallocate(i1u);
      deallocate(j0); deallocate(j0v);
      deallocate(j1); deallocate(j1v);

      end subroutine write_partial_files !]
! ----------------------------------------------------------------------
      subroutine get_nc_info(ncid, ncinfo)  ![

      implicit none        
        
      ! input/output
      integer,           intent(in)  :: ncid
      type(ncstructure), intent(out) :: ncinfo

      ! local
      character(len=20) :: attname
      character(len=200):: attvalue
      real              :: attparam
      character(len=10),allocatable,dimension(:) :: dnames
      character(len=10)  :: dname
      character(len=20) :: vname
      integer :: ierr, i, xtype,natts,dsize
      integer,dimension(5) :: dimids
      integer :: dimid,ndims,ivar,nvars
      integer :: nrecords

      ! local allocation
      ierr = nf90_inquire(ncid,nvariables=nvars)

      ncinfo%nvars = nvars
      allocate(ncinfo%varname(nvars))
      allocate(ncinfo%vartype(nvars))
      allocate(ncinfo%varpart(nvars))
      allocate(ncinfo%bryvar(nvars))
      allocate(ncinfo%xtype(nvars))
      allocate(ncinfo%ndims(nvars))
      allocate(ncinfo%dims(nvars))
      allocate(ncinfo%natts(nvars))
      allocate(ncinfo%atts(nvars))

      allocate(dnames(ndims))
      ncinfo%varpart = .false.

      ! append info to allocated arrays apart of derived type
      do ivar = 1,nvars
        
        ! Get information per variable
        ierr = nf90_inquire_variable(ncid,ivar,vname,xtype,ndims,dimids,natts)
        
        ! name of variable
        ncinfo%varname(ivar) = vname

        ! check variable to determine boundary 
        if (index(vname,'_west')>0) then
          ncinfo%bryvar(ivar) = 'w'
        elseif (index(vname,'_east')>0) then
          ncinfo%bryvar(ivar) = 'e'
        elseif (index(vname,'_south')>0) then
          ncinfo%bryvar(ivar) = 's'
        elseif (index(vname,'_north')>0) then
          ncinfo%bryvar(ivar) = 'n'
        else
          ncinfo%bryvar(ivar) = 'i'
        endif

        ! variable type
        ncinfo%xtype(ivar) = xtype

        ! r,u,v, or s (skip)
        ncinfo%vartype(ivar) = 'r'

        ! dimensions
        ncinfo%ndims(ivar) = ndims
        allocate(ncinfo%dims(ivar)%dimsizes(ndims))
        allocate(ncinfo%dims(ivar)%dimnames(ndims))
        do i = 1,ndims
          ierr = nf90_inquire_dimension(ncid,dimids(i),dname,dsize)
          ncinfo%dims(ivar)%dimnames(i) = dname
          ncinfo%dims(ivar)%dimsizes(i) = dsize
          if (index(dname,'xi_')>0 .or. index(dname,'eta_')>0) then
            ncinfo%varpart(ivar) = .true.
          endif
          if (index(dname,'eta_u')>0 ) then
            ncinfo%dims(ivar)%dimnames(i) = 'eta_rho'
          endif
          if (index(dname,'xi_v')>0 ) then
            ncinfo%dims(ivar)%dimnames(i) = 'xi_rho'
          endif
          if (index(dname,'xi_u')>0 ) then
            ncinfo%vartype(ivar) = 'u'
          endif
          if (index(dname,'eta_v')>0) then
            ncinfo%vartype(ivar) = 'v'
          endif
          if (index(dname,'xi_rho')>0) then
            ncinfo%gsize(1) = dsize
          endif
          if (index(dname,'eta_rho')>0) then
            ncinfo%gsize(2) = dsize
          endif
          if (index(dname,'s_rho')>0) then
            ncinfo%gsize(3) = dsize
          endif
          if (index(dname,'s_w')>0) then
            ncinfo%gsize(3) = dsize
          endif

        enddo

        ! check for xi_u/eta_v dims
        if (ndims>=2) then
          if ((index(ncinfo%dims(ivar)%dimnames(2),'eta_v')>0).and.
     &        (index(ncinfo%dims(ivar)%dimnames(1),'xi_u')>0)) then
              ncinfo%vartype(ivar) = 'j'
          endif
        endif

        ! attributes
        ncinfo%natts(ivar) = natts
        allocate(ncinfo%atts(ivar)%attnames(natts))
        allocate(ncinfo%atts(ivar)%attvalues(natts))
        allocate(ncinfo%atts(ivar)%attparams(natts))
        do i = 1,natts
          ierr = nf90_inq_attname(ncid,ivar,i,attname)
          ncinfo%atts(ivar)%attnames(i)  = attname
          ierr = nf90_get_att(ncid,ivar,attname,attvalue)
          if (ierr/=0) then
            ierr = nf90_get_att(ncid,ivar,attname,attparam)
            ncinfo%atts(ivar)%attparams(i) = attparam
            ncinfo%atts(ivar)%attvalues(i) = 'n/a'
          else 
            ncinfo%atts(ivar)%attvalues(i) = attvalue
            ncinfo%atts(ivar)%attparams(i) = 0
          endif
        enddo
        
      enddo

      end subroutine get_nc_info !]
! ----------------------------------------------------------------------
      subroutine insert_nodes(fname, mynode, nnodes) ![
      ! Insert date and node number string into filename
      ! If the nodate argument is present, only add the node number
      implicit none

      ! import/export
      character(len=*),intent(inout):: fname 
      integer, intent(in) :: mynode, nnodes

      ! Local
      character(len=6)   :: indxstr
      integer :: label
      logical :: ncf
      integer :: istart

      ncf =  findstr(fname,'.nc',istart) 
      if (ncf) then
        fname(istart:istart+2) = '   '
      else
        print *,'error insert_node. No .nc found in file name'
        stop
      endif


      if (nnodes<=10) then
        write(indxstr,'(A,I0.1)') '.',mynode
      elseif (nnodes<=100) then
        write(indxstr,'(A,I0.2)') '.',mynode
      elseif (nnodes<=1000) then
        write(indxstr,'(A,I0.3)') '.',mynode
      elseif (nnodes<=10000) then
        write(indxstr,'(A,I0.4)') '.',mynode
      elseif (nnodes<=100000) then
        write(indxstr,'(A,I0.5)') '.',mynode
      elseif (nnodes>100001) then
         print *, 'ROMS is not ready for this future'
         stop
      endif
      fname = trim(fname) //trim(indxstr)

      fname = trim(fname) // '.nc'

      end subroutine insert_nodes !]
! ----------------------------------------------------------------------
      subroutine system_mem_usage(valueRSS) ![

      use ifport !if on intel compiler

      implicit none
      integer, intent(out) :: valueRSS

      character(len=200) :: filename = ' '
      character(len=80)  :: line
      character(len=8)   :: pid_char = ' '
      integer            :: pid
      logical            :: ifxst
      integer            :: ios

      valueRSS = -1    ! return negative number if not found

      !--- get process ID
      pid = getpid()
      write(pid_char, '(I8)') pid
      filename = '/proc/' // trim(adjustl(pid_char)) // '/status'

      !--- read system file
      inquire(file = filename, exist = ifxst)
      if (.not. ifxst) then
        write (*,*) 'system file does not exist'
        return
      endif

      open(unit = 100, file = filename, action = 'read')
      do
        read(100, '(a)', iostat=ios) line
        if (ios /= 0) then
          print *, 'error in RSS calc.'
          exit
        endif
        if (line(1:6) .eq. 'VmRSS:') then
          read(line(7:), *) valueRSS
          exit  ! Exit loop when VmRSS line is found
        endif
      enddo
      close(100)

      end subroutine system_mem_usage !]
! ----------------------------------------------------------------------
      logical function findstr(string,pattern,istart) ![
      implicit none

      !input/output
      character(len=*),intent(in)  :: string        ! string
      character(len=*),intent(in)  :: pattern       ! desired pattern to find within string
      integer,optional,intent(out) :: istart

      !local
      integer :: nl,nlv,i

      nl  = len(trim(string))
      nlv = len(pattern)

      findstr = .false.
      do i = 1,nl-nlv+1
         if (string(i:i+nlv-1) == pattern) then
          findstr = .true.
          exit
         endif
      enddo

      if (present(istart)) then
        if (findstr) then
          istart=i                                  ! return string starting index
        else
          istart=0
        endif
      endif

      end function findstr !]
! ----------------------------------------------------------------------

      end module roms_part
