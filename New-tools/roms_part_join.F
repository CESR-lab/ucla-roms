      module roms_part_join
      
      ! ROMS functionality to facilitate partioning and joining of files

      use netcdf
      use mpi
      use nc_read_write

      private
       
      ! mpi variables
      integer,public :: part_comm,mynode,nnodes  ! mpi vars

      ! range of partial files for this node
      integer :: part0,part1,myparts,part

      ! global dimensions
      integer :: gnx,gny,nz,nzw

      ! Files
      integer,public :: nfiles ! number of src files
      integer,public :: nparts ! number of partial files
      character(len=128),public :: fullsrcfile,srcfile,joinfile,partfile
      character(len=128),public :: grdfile
      character(len=128),public :: gpath,spath
      character(len=128),public,allocatable, dimension(:) :: srcfiles
      character(len=128),public,allocatable, dimension(:) :: partfiles

      ! Only used by mynode=0
      integer,dimension(:,:),allocatable :: iloc,jloc
      integer,dimension(:,:),allocatable :: ilcu,jlcv
      integer,dimension(:,:),allocatable :: node_map
      integer,dimension(:)  ,allocatable :: npi_c,npj_c

      ! Dimension names for in netcdf files
      character(len=7),public,parameter :: dn_xr = 'xi_rho'
      character(len=7),public,parameter :: dn_yr = 'eta_rho'
      character(len=7),public,parameter :: dn_xu = 'xi_u'
      character(len=7),public,parameter :: dn_yv = 'eta_v'
      character(len=7),public,parameter :: dn_zr = 's_rho'
      character(len=7),public,parameter :: dn_zw = 's_w'
      character(len=7),public,parameter :: dn_tm = 'time'
      
      type,public :: ncstructure
        integer                                     :: nrecords  ! number of records
        integer                                     :: nvars     ! number of variables
        integer,dimension(3)                        :: gsize     ! dimension size of the global grid (rho)
        character(len=30),dimension(:),allocatable  :: varname   ! name, if u/v or rho
        character(len=1),dimension(:), allocatable  :: vartype   ! 'f/r/u f/r/v f/r/w f/t'
        logical,dimension(:),allocatable            :: varpart   ! partitionable or not
        integer,dimension(:),allocatable            :: xtype     ! variable type
        integer,dimension(:),allocatable            :: ndims  
        integer,dimension(:),allocatable            :: natts
        type(dimstruct),dimension(:),allocatable    :: dims
        type(attstruct),dimension(:),allocatable    :: atts
      end type ncstructure

      type,public :: dimstruct
        character(len=7),dimension(:),allocatable ::dimnames
        integer,         dimension(:),allocatable ::dimsizes
      end type dimstruct

      type,public :: attstruct
        character(len=20),dimension(:),allocatable ::attnames
        character(len=200),dimension(:),allocatable ::attvalues
      end type attstruct

      type(ncstructure),public :: ncinfo
      
      public :: cmd_args_partit
      public :: cmd_args_ncjoin
      public :: partition_mask
      public :: make_partial_files
      public :: write_partial_files
      public :: make_joined_file
      public :: write_joined_file

      contains
! ----------------------------------------------------------------------
      subroutine cmd_args_partit(npartx,nparty)  ![
      implicit none

!     import/export
      integer,intent(out) :: npartx,nparty

      integer :: nargs,i,idx
      character(len=128) :: arg
      character(len=128) :: gpath,fpath

      nargs =  command_argument_count()

      if (nargs < 4) call show_use_partit
      
      ! first 2 should be number of subdomains in x- or y- direction
      call get_command_argument(1,arg);
      read(arg,*) npartx
      call get_command_argument(2,arg);
      read(arg,*) nparty
      print *,nargs,npartx,nparty

      call get_command_argument(3,arg)
!     read(arg,*) grdfile

      idx  = index(arg,'/',back=.true.)
!     gpath = arg(1:idx)
!     grdfile = arg(idx+1:128)
!     print *, gpath
!     print *, grdfile
      grdfile = arg
      
      nfiles = nargs - 3
      allocate(srcfiles(nfiles))
      do i = 4,nargs
        call get_command_argument(i,arg)
!       read(arg,*) srcfiles(i-3)
        idx  = index(arg,'/',back=.true.)
        spath = arg(1:idx)
        srcfiles(i-3) = arg(idx+1:128)
      enddo
      print *,grdfile,srcfiles(1)

      end subroutine cmd_args_partit !]
! ----------------------------------------------------------------------
      subroutine cmd_args_ncjoin  ![
      implicit none

      integer :: nargs,i
      character(len=128) :: arg

      nargs =  command_argument_count()

      if (nargs < 1) call show_use_ncjoin
      
      nparts = nargs
      allocate(partfiles(0:nparts-1))
      do i = 1,nargs
        call get_command_argument(i,arg)
        read(arg,*) partfiles(i-1)
      enddo
      print *,partfiles(1)

      myparts = ceiling(nparts*1.0/nnodes)
      part0 = myparts*mynode
      part1 = min(myparts*(mynode+1) -1,nparts-1)

      end subroutine cmd_args_ncjoin !]
! ----------------------------------------------------------------------
      subroutine show_use_partit  ![
      implicit none

      print *, 'Usage of partit should be: ',
     &       'partit npx npy grdfile.nc, file1.nc...'
     
      stop

      end subroutine show_use_partit !]
! ----------------------------------------------------------------------
      subroutine show_use_ncjoin  ![
      implicit none

      print *, 'Usage of ncjoin should be: ',
     &       'ncjoin file1.??.nc ...'
     
      stop

      end subroutine show_use_ncjoin !]
! ----------------------------------------------------------------------
      subroutine partition_mask(npx,npy)   ![
      implicit none

      ! input/outputs
      integer,           intent(in) :: npx,npy

      ! local
      real,dimension(:,:),allocatable :: mask,rflx
      real    :: msk_mx,rfx_mx
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
      allocate(rflx(0:gnx+1,0:gny+1))
      call ncread(ncid, 'mask_rho',mask)
      call ncread(ncid, 'river_flux',rflx)
      ierr = nf90_close(ncid)
!     mask = 1

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
          rfx_mx = maxval(rflx(iloc(npi,2):iloc(npi,3),jloc(npj,2):jloc(npj,3)))
          if (msk_mx > 0.or.rfx_mx>0) then
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
      subroutine make_partial_files ![
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
      character(len=7),dimension(5)  ::dname
      integer         ,dimension(5)  ::dsize
      character(len=20),dimension(5) ::attname
      character(len=200),dimension(5)::attvalue
      integer :: ierr, ncid,npid
      integer :: natts,idm,vid
      integer,dimension(2) :: subdompos

      ! open up reference netcdf for global attributes
      ierr = nf90_open(trim(fullsrcfile), nf90_nowrite, ncid)
      if (ierr/=nf90_noerr) 
     &  call handle_ierr(ierr,'opening file: ',trim(fullsrcfile) )

      call get_nc_info(ncid)

      do part = 1,nparts
        i = npi_c(part)
        j = npj_c(part)
    
        fname = srcfile
        call insert_nodes(fname,part-1,nparts)

!       ierr=nf90_create(trim(fname),nf90_netcdf4,npid)
        ierr=nf90_create(trim(fname),nf90_classic_model,npid)
        if (ierr/=nf90_noerr) 
     &    call handle_ierr(ierr,'creating file: ',trim(fname) )

        ! global rho-dimensions
        ierr=nf90_put_att(npid, nf90_global, 'globaldims', ncinfo%gsize)

        ! subdompos: i,j position of sub-domain
        subdompos = (/ i, j /)
        ierr=nf90_put_att(npid, nf90_global, 'subdompos', subdompos)

        ! partition: node, total nodes, llc x, llc y
        ibuff = (/part-1, nparts, iloc(i,2), jloc(j,2) /)
        ierr=nf90_put_att(npid, nf90_global, 'partition', ibuff)

        ! starting indices: i0, j0, i0u, j0v
        ibuff = (/iloc(i,2),jloc(j,2),ilcu(i,2),jlcv(j,2) /)
        ierr=nf90_put_att(npid, nf90_global, 'start_idx', ibuff)

        ! starting indices: i0, j0, i0u, j0v
        ibuff = (/iloc(i,3),jloc(j,3),ilcu(i,3),jlcv(j,3) /)
        ierr=nf90_put_att(npid, nf90_global, 'end_idx', ibuff)

          ! Neigbors: N,NE,E,SE,S,SW,W,NW (READ IN THAT ORDER)
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
          ierr = nf90_copy_att(ncid,nf90_global,att_name,npid,nf90_global)
        enddo 

        do ivar = 1,ncinfo%nvars
          vname   = ncinfo%varname(ivar)
          xtype   = ncinfo%xtype(ivar)
          ndims   = ncinfo%ndims(ivar)
          dsize(1:ndims) = ncinfo%dims(ivar)%dimsizes
          dname(1:ndims) = ncinfo%dims(ivar)%dimnames
          natts   = ncinfo%natts(ivar)
          attname(1:natts) = ncinfo%atts(ivar)%attnames
          attvalue(1:natts)= ncinfo%atts(ivar)%attvalues

          do idm =1,ndims
            if (index(dname(idm),'xi_rho') > 0)  dsize(idm)=iloc(i,1)
            if (index(dname(idm),'xi_u') > 0)    dsize(idm)=ilcu(i,1)
            if (index(dname(idm),'eta_rho') > 0) dsize(idm)=jloc(j,1)
            if (index(dname(idm),'eta_v') > 0)   dsize(idm)=jlcv(j,1)
            if (index(dname(idm),'time') > 0)    dsize(idm)= 0 
          enddo

          if (ncinfo%vartype(ivar).eq.'s') then
            continue
          else
            vid = nccreate(npid,vname,dname(1:ndims),dsize(1:ndims),xtype)
            do iatt=1,natts
              ierr = nf90_put_att(npid,vid,attname(iatt),attvalue(iatt))
            enddo
          endif

        enddo
        ierr = nf90_close(npid)

      enddo 

      ierr = nf90_close(ncid)

      end subroutine make_partial_files !]
! ----------------------------------------------------------------------
      subroutine write_partial_files ![
      ! get info from first partial file in range
      ! read block of data for each partitionable variable
      ! and write to partial files, for all records in srcfile
      implicit none

      ! local
      character(len=128) :: fname
      integer :: ierr, ncid,npid,dimid
      integer,dimension(4) :: llidx
      integer,dimension(:),allocatable :: i0,j0
      integer,dimension(:),allocatable :: i1,j1
      integer,dimension(:),allocatable :: i0u,i1u
      integer,dimension(:),allocatable :: j0v,j1v
!     integer :: lnx,lny
      integer :: li0,li1,lj0,lj1,li0u,li1u,lj0v,lj1v
      real,dimension(:,:,:),allocatable :: var3d
      character(len=20) :: vname
      integer,dimension(:),allocatable :: is,ie,js,je

      integer :: varid,ivar
      integer :: record,nrecords
      integer,dimension(4) :: startr,startw
      real,dimension(:),allocatable :: var1d
      integer :: sz1,sz2,sz3,ndims
      integer :: RSS

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

      print *, 'scanning files'
      do part = part0,part1
    
        fname = srcfile
        call insert_nodes(fname,part,nparts)

        ierr = nf90_open(trim(fname),nf90_write,npid)
        if (ierr/=nf90_noerr) 
     &    call handle_ierr(ierr,'opening file: ',trim(fname) )

        if (part.eq.part0) call get_nc_info(npid)

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

        ierr = nf90_close(npid)

      enddo 

      li0 = minval(i0)
      li1 = maxval(i1)
      lj0 = minval(j0)
      lj1 = maxval(j1)
      li0u = minval(i0u)
      li1u = maxval(i1u)
      lj0v = minval(j0v)
      lj1v = maxval(j1v)

      nz = ncinfo%gsize(3)

      call system_mem_usage(RSS)
      print * ,'Mem usage: ', RSS/1000

!     ierr = nf90_open_par(trim(fullsrcfile),nf90_nowrite,part_comm,MPI_Info_Null, ncid)
      ierr = nf90_open(trim(fullsrcfile),nf90_nowrite,ncid)
      if (ierr/=0) call handle_ierr(ierr,'opening: ',trim(fullsrcfile))

      ierr = nf90_inq_dimid(ncid,'time',dimid)
       if (ierr /= nf90_noerr) call handle_ierr(ierr)
      ierr = nf90_inquire_dimension(ncid,dimid,len=nrecords)
       if (ierr /= nf90_noerr) call handle_ierr(ierr)
     
      print *,'nrecords: ',nrecords


      do record = 1,nrecords

        if (mynode==0) print *,'Processing record: ',record
      do ivar=1,ncinfo%nvars

        vname = ncinfo%varname(ivar)
        if (mynode==0) print *,'Processing: ',trim(vname)

        ndims = ncinfo%ndims(ivar)

        if (ncinfo%varpart(ivar)) then
          select case(ndims)
          case(1)
            ! save for later
          case(2)
            ! save for later (not present in pachug_rst.nc)
          case(3)
            if (index(ncinfo%dims(ivar)%dimnames(3),'time')>0) then
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
            ! save for later (not present in pachug_rst.nc)
            endif
          case(4)
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
        else
          select case(ndims)
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
              allocate(var3d(sz1,1,1))
              is = 1; ie = sz1
              js = 1; je = 1
            else
              if (record>1) continue
              startr = (/1,1,1,1/)
              startw = (/1,1,1,1/)
              sz1 = ncinfo%dims(ivar)%dimsizes(1)
              sz2 = ncinfo%dims(ivar)%dimsizes(2)
              allocate(var3d(sz1,sz2,1))
              is = 1; ie = sz1
              js = 1; je = sz2
            endif
          end select

        endif

        call system_mem_usage(RSS)
        print * ,'Mem usage 1: ', vname,'  ',RSS/1000

        call ncread(ncid,vname,var3d,startr)

        call system_mem_usage(RSS)
        print * ,'Mem usage 2: ', vname,'  ',RSS/1000

        do part = part0,part1  ! zero based
          fname = srcfile
          call insert_nodes(fname,part,nparts)
          ierr = nf90_open(trim(fname), nf90_write, npid)
          call ncwrite(npid,vname,var3d(is(part):ie(part),js(part):je(part),:),startw)
          ierr = nf90_close(npid)
        enddo

        call system_mem_usage(RSS)
        print * ,'Mem usage 3: ', vname,'  ',RSS/1000

        deallocate(var3d)

        call system_mem_usage(RSS)
        print * ,'Mem usage 4: ', vname,'  ',RSS/1000

!       print *,'mynode: ', mynode
      enddo  ! ivars

      enddo ! records
      ierr = nf90_close(ncid)

      deallocate(i0); deallocate(i0u);
      deallocate(i1); deallocate(i1u);
      deallocate(j0); deallocate(j0v);
      deallocate(j1); deallocate(j1v);

      end subroutine write_partial_files !]
! ----------------------------------------------------------------------
      subroutine make_joined_file ![
      ! single node only
      implicit none

      ! local
      character(len=128) :: partfile
      character(len=128) :: att_name
      integer :: iatt,ngatt,part
      integer :: i,j,varid
      integer :: dimids(4),ndims
      integer :: ivar,nvars,xtype
      character(len=20) :: vname
      character(len=7),dimension(5)  ::dname
      integer         ,dimension(5)  ::dsize
      character(len=20),dimension(5) ::attname
      character(len=200),dimension(5)::attvalue
      integer :: ierr, ncid,npid
      integer :: natts,idm,vid
      integer :: idx1
      integer,dimension(3)  ::gdims
      real,dimension(:,:,:),allocatable :: var3d
      integer,dimension(4) :: startr

      partfile = trim(partfiles(part0))
      idx1 = index(partfile,'.nc',back=.true.)
      joinfile = partfile(1:idx1-1)
      idx1 = index(joinfile,'.',back=.true.)
      joinfile = joinfile(1:idx1-1) //'.nc'
      print *, trim(partfile),' ',trim(joinfile)

      ! Open first partial file in range 
      ierr = nf90_open(trim(partfile), nf90_nowrite, npid)
      if (ierr/=nf90_noerr) 
     &    call handle_ierr(ierr,'opening file: ',trim(partfile) )

      print *,'getting info from: ',trim(partfile)
      call get_nc_info(npid)

      ierr = nf90_get_att(npid,nf90_global,'globaldims',gdims)
      if (ierr/=nf90_noerr) 
        ! Global dimensions attribute not found, use old method ![

        partfile = partfiles(nparts)
        ierr = nf90_open(trim(partfile), nf90_nowrite, npid)
        if (ierr/=nf90_noerr) 
     &    call handle_ierr(ierr,'opening file: ',trim(partfile) )

        ierr=nf90_get_att(npid, nf90_global, 'partition', ibuff)

        nf90_inquire(npid,ndimensions=ndims)
        do i = 1,ndims
          ierr = nf90_inquire_dimension(ncid,i,dname,dsize)
          if (index(dname,'xi_rho')>0)  gdims(1) = ibuff(3) + dsize -1
          if (index(dname,'eta_rho')>0) gdims(2) = ibuff(4) + dsize -1
          if (index(dname,'s_rho')>0)   gdims(3) = dsize
        enddo
        ierr = nf90_close(npid)

      endif
      gnx = gdims(1)
      gny = gdims(2)
      nz  = gdims(3)
      
      if (mynode==0) then ! single node to create joined file

        print *,'Creating joined file'

        ierr=nf90_create(trim(joinfile),nf90_netcdf4,ncid)
!       ierr=nf90_create(trim(joinfile),nf90_classic_model,ncid)
        if (ierr/=nf90_noerr) 
     &    call handle_ierr(ierr,'creating file: ',trim(joinfile) )

        ! read all global attributes from the partial file and write them
        ierr = nf90_inquire(npid, nAttributes=ngatt)
        if (ierr/=nf90_noerr) 
     &    call handle_ierr(ierr,'getting atts: ',trim(partfile) )
        do iatt = 1,ngatt
          ierr = nf90_inq_attname(npid,nf90_global,iatt,att_name)
          if (index(att_name,'globaldims')>0) cycle
          if (index(att_name,'subdompos')>0) cycle
          if (index(att_name,'partition')>0) cycle
          if (index(att_name,'start_idx')>0) cycle
          if (index(att_name,'neighbors')>0) cycle
          ierr = nf90_copy_att(npid,nf90_global,att_name,ncid,nf90_global)
        enddo 

        do ivar = 1,ncinfo%nvars
          vname   = ncinfo%varname(ivar)
          xtype   = ncinfo%xtype(ivar)
          ndims   = ncinfo%ndims(ivar)
          dsize(1:ndims) = ncinfo%dims(ivar)%dimsizes
          dname(1:ndims) = ncinfo%dims(ivar)%dimnames
          natts   = ncinfo%natts(ivar)
          attname(1:natts) = ncinfo%atts(ivar)%attnames
          attvalue(1:natts)= ncinfo%atts(ivar)%attvalues

          do idm =1,ndims
            if (index(dname(idm),'xi_rho') > 0)  dsize(idm)=gnx
            if (index(dname(idm),'xi_u') > 0)    dsize(idm)=gnx-1
            if (index(dname(idm),'eta_rho') > 0) dsize(idm)=gny
            if (index(dname(idm),'eta_v') > 0)   dsize(idm)=gny-1
            if (index(dname(idm),'time') > 0)    dsize(idm)= 0 
          enddo

          vid = nccreate(ncid,vname,dname(1:ndims),dsize(1:ndims),xtype)
          ierr = nf90_var_par_access(ncid, varid, nf90_collective)
          do iatt=1,natts
            ierr = nf90_put_att(ncid,vid,attname(iatt),attvalue(iatt))
          enddo
        enddo

        ierr = nf90_close(ncid)

      endif ! mynode==0

      ierr = nf90_close(npid)

      if (gdims(1)==0) then 
    
      endif 
  

      end subroutine make_joined_file !]
! ----------------------------------------------------------------------
      subroutine write_joined_file ![
      ! read data for each variable in partial files
      ! and write to joinedl file, for all records
      implicit none

      ! local
      integer :: ierr, ncid,npid,dimid
      integer,dimension(4) :: llidx
      integer :: nx,ny,nxu=0,nyv=0
      integer,dimension(:),allocatable :: i0,j0
      integer,dimension(:),allocatable :: i1,j1
      integer,dimension(:),allocatable :: i0u,i1u
      integer,dimension(:),allocatable :: j0v,j1v
      integer :: lnx,lny,li0,li1,lj0,lj1,li0u,li1u,lj0v,lj1v
      real,dimension(:,:,:),allocatable :: var3d
      character(len=30) :: vname
      integer,dimension(:),allocatable :: is,ie,js,je


      integer :: varid,ivar
      integer :: record,nrecords
      integer,dimension(4) :: startr,startw
      real,dimension(:),allocatable :: var1d
      integer :: sz1,sz2,sz3,ndims
      integer :: RSS


      ! get info about variables from first partial file in range

      ! Read partial files in range
      allocate(i0(part0:part1))
      allocate(j0(part0:part1))
      allocate(i1(part0:part1))
      allocate(j1(part0:part1))
      allocate(i0u(part0:part1))
      allocate(i1u(part0:part1))
      allocate(j0v(part0:part1))
      allocate(j1v(part0:part1))

      allocate(is(part0:part1))
      allocate(js(part0:part1))
      allocate(ie(part0:part1))
      allocate(je(part0:part1))

      print *,'scanning partial files'
      do part = part0,part1 
    
        partfile = partfiles(part)

        ierr = nf90_open(trim(partfile),nf90_write,npid)
        if (ierr/=nf90_noerr) 
     &    call handle_ierr(ierr,'opening file: ',trim(partfile) )


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

        ierr = nf90_get_att(npid,nf90_global,'partition',ibuff)
        i0(part) = ibuff(3)
        j0(part) = ibuff(4)
        i0u(part)= max(i0(part)-1,1)
        j0v(part)= max(j0(part)-1,1)

        nf90_inquire(npid,ndimensions=ndims)
        do i = 1,ndims
          ierr = nf90_inquire_dimension(ncid,i,dname,dsize)
          if (index(dname,'xi_rho')>0)  nx = dsize
          if (index(dname,'eta_rho')>0) ny = dsize
          if (index(dname,'xi_u')>0)    nxu = dsize
          if (index(dname,'eta_v')>0)   nyv = dsize
        enddo
        i1(part) = i0(part) + nx-1
        j1(part) = j0(part) + ny-1
        i1u(part)= i0u(part)+ nxu-1
        j1v(part)= j0v(part)+ nyv-1

      enddo

      li0 = minval(i0)
      li1 = maxval(i1)
      lj0 = minval(j0)
      lj1 = maxval(j1)
      li0u = minval(i0u)
      li1u = maxval(i1u)
      lj0v = minval(j0v)
      lj1v = maxval(j1v)
      lnx = li1 - li0 + 1
      lny = lj1 - lj0 + 1

      print *,'do all records: ',nrecords

      do record = 1,ncinfo%nrecords

        if (mynode==0) print *,'Processing record: ',record

        do ivar=1,ncinfo%nvars

          vname = ncinfo%varname(ivar)
          if (mynode==0) print *,'Processing: ',trim(vname)

            ndims = ncinfo%ndims(ivar)

        if (ncinfo%varpart(ivar)) then
          select case(ndims)
          case(1)
            ! save for later
          case(2)
            ! save for later (not present in pachug_rst.nc)
          case(3)
            if (index(ncinfo%dims(ivar)%dimnames(3),'time')>0) then
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
            ! save for later (not present in pachug_rst.nc)
            endif
          case(4)
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
        else
          select case(ndims)
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
              allocate(var3d(sz1,1,1))
              is = 1; ie = sz1
              js = 1; je = 1
            else
              if (record>1) continue
              startr = (/1,1,1,1/)
              startw = (/1,1,1,1/)
              sz1 = ncinfo%dims(ivar)%dimsizes(1)
              sz2 = ncinfo%dims(ivar)%dimsizes(2)
              allocate(var3d(sz1,sz2,1))
              is = 1; ie = sz1
              js = 1; je = sz2
            endif
          end select

        endif

        call system_mem_usage(RSS)
        print * ,'Mem usage 1: ', vname,'  ',RSS/1000

!       print *, li0,li1,lj0,lj1
        var3d = 0

        do part = part0,part1  ! zero based
          partfile = partfiles(part)
          ierr = nf90_open(trim(partfile), nf90_nowrite, npid)
          call ncread(npid,vname,var3d(is(part):ie(part),js(part):je(part),:),startw)
          ierr = nf90_close(npid)
        enddo

        call system_mem_usage(RSS)
        print * ,'Mem usage 2: ', vname,'  ',RSS/1000

        print *,'maxval: ',maxval(abs(var3d))
        ierr = nf90_open_par(trim(joinfile),nf90_write,part_comm,MPI_Info_Null, ncid)
        ierr = nf90_var_par_access(ncid, ivar, nf90_collective)
!       ierr = nf90_open(trim(joinfile),nf90_write,ncid)
        if (ierr/=nf90_noerr) 
     &    call handle_ierr(ierr,'opening file: ',trim(joinfile) )
        call ncwrite(ncid,vname,var3d,startr)
        ierr = nf90_close(ncid)


        call system_mem_usage(RSS)
        print * ,'Mem usage 3: ', vname,'  ',RSS/1000

        deallocate(var3d)

        call system_mem_usage(RSS)
        print * ,'Mem usage 4: ', vname,'  ',RSS/1000

!       print *,'mynode: ', mynode
      enddo  ! ivars

      enddo ! records
      ierr = nf90_close(ncid)

      deallocate(i0)
      deallocate(i1)
      deallocate(j0)
      deallocate(j1)

      end subroutine write_joined_file !]
! ----------------------------------------------------------------------
      subroutine get_nc_info(ncid)  ![

      implicit none        
        
      ! input/output
      integer,           intent(in)  :: ncid

      ! local
      character(len=20) :: attname
      character(len=200):: attvalue
      character(len=7),allocatable,dimension(:) :: dnames
      character(len=7)  :: dname
      character(len=20) :: vname
      integer :: ierr, i, xtype,natts,dsize
      integer,dimension(5) :: dimids
      integer :: dimid,ndims,ivar,nvars


      print *,'Getting info'
      ! local allocation
      ierr = nf90_inquire(ncid,nvariables=nvars)

      ncinfo%nvars = nvars
      allocate(ncinfo%varname(nvars))
      allocate(ncinfo%vartype(nvars))
      allocate(ncinfo%varpart(nvars))
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

        ! variable type
        ncinfo%xtype(ivar) = xtype

        ! r,u,v, or s (skip)
        ncinfo%vartype(ivar) = 'r'

        ! dimensions
        ncinfo%gsize = 1
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
          if (index(dname,'time')>0) then
            ncinfo%nrecords = dsize
          endif
        enddo

        ! attributes
        ncinfo%natts(ivar) = natts
        allocate(ncinfo%atts(ivar)%attnames(natts))
        allocate(ncinfo%atts(ivar)%attvalues(natts))
        do i = 1,natts
          ierr = nf90_inq_attname(ncid,ivar,i,attname)
          ierr = nf90_get_att(ncid,ivar,attname,attvalue)
          ncinfo%atts(ivar)%attnames(i) = attname
          ncinfo%atts(ivar)%attvalues(i) = attvalue
        enddo
        
      enddo  ! 1,nvars

      print *,'nrecords: ',ncinfo%nrecords

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

      end module roms_part_join
