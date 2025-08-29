      module roms_part
      
      ! ROMS functionality to facilitate partioning and joining of files

      use netcdf
      use mpi
      use nc_read_write
      use netcdf

      private
       
      ! usage flags
      logical,public :: del_part_files,classic_part

      ! partition/joining dimensions
      public :: npartx,nparty

      ! mpi variables
      integer,public :: part_comm,mynode,nnodes  ! mpi vars

      ! range of partial files for this node
      integer :: myparts,part
      integer,public :: part0,part1

      ! global dimensions
      integer,public :: gnx,gny,nz,nzw

      ! Files
      integer,public :: nfiles ! number of src files
      character(len=128),public :: grdfile,srcfile
      character(len=128),public :: fullsrcfile,joinfile,partfile
      character(len=128),public :: gpath,spath
      character(len=128),public,allocatable, dimension(:) :: srcfiles
      character(len=128) :: fname
      character(len=128),public,allocatable, dimension(:) :: partfiles
      character(len=128),public,dimension(:),allocatable:: allfiles
      integer,public,dimension(:),allocatable :: file_indices

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
        integer                                     :: nrecords  ! number of records
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
        character(len=10),dimension(:),allocatable  ::dimnames
        integer,         dimension(:),allocatable   ::dimsizes
      end type dimstruct

      type,public :: attstruct
        character(len=20),dimension(:),allocatable  ::attnames
        character(len=200),dimension(:),allocatable ::attvalues
        real,dimension(:),allocatable               ::attparams
      end type attstruct
 
      type,public :: filestruct
        integer,dimension(:),allocatable            :: numfiles
        character(len=128),dimension(:),allocatable :: filelist
      end type filestruct
      
      public :: cmd_args_optimal
      public :: cmd_args_map
      public :: cmd_args_partit
      public :: cmd_args_ncjoin
      public :: partition_mask
      public :: check_mask
      public :: make_partial_files
      public :: write_partial_files
      public :: make_joined_file
      public :: write_joined_file

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
      integer :: nflags = 0  ! assume no flags present
      character(len=128) :: arg

      nargs =  command_argument_count()

      ! first 2 should be number of subdomains in x- or y- direction
      call get_command_argument(1,arg);
      ! if instead we have the classic flag present, all indices are shifted by 1
      if (index(arg,'-classic')>0) then
        if (nargs < 5) call show_use_partit
        print*,'Partition set to classic mode (ignoring mpi_masking)'
        classic_part = .true.
        nflags=1
        ! get npartx
        call get_command_argument(1+nflags,arg);
        read(arg,*) npartx
      else
        if (nargs < 4) call show_use_partit
        read(arg,*) npartx
      endif

      call get_command_argument(2+nflags,arg);
      read(arg,*) nparty
      print *,nargs,npartx,nparty

      call get_command_argument(3+nflags,arg)
      grdfile = trim(adjustl(arg))
      
      nfiles = nargs - 3 - nflags
      allocate(srcfiles(nfiles))
      do i = (4+nflags),nargs
        call get_command_argument(i,arg)
        srcfiles(i-3-nflags) = trim(adjustl(arg))
      enddo
      print *,grdfile,srcfiles(1)

      end subroutine cmd_args_partit !]
! ----------------------------------------------------------------------
      subroutine cmd_args_ncjoin  ![
      implicit none

      integer :: nargs,nfiles,nflags,idx,total_files,i,count
      character(len=128) :: arg,file,newfile

      ! baseline options for flags
      del_part_files = .false.

      nargs =  command_argument_count()

      if (nargs < 1) call show_use_ncjoin
      
      nflags = 0
      nfiles = 0
      do i = 1,nargs
        call get_command_argument(i,arg)
        ! look for flags ('--' or '-')
        ! "delete" flag: removes all partial files
        if (index(arg,'-delete')>0) then
          if (mynode==0) print*,'Flag to delete partial files raised'
          del_part_files = .true.
          nflags         = nflags+1
        ! Per 'show_use_ncjoin' flags go in front. Any following input must be
        ! a netcdf file. This convention allows for allocation of an array whose
        ! length is exactly as many files being read in.
        else
          if (.not.(allocated(allfiles))) allocate(allfiles(0:nargs-nflags-1))  ! -1 for zero-based
          read(arg,*) allfiles(i-nflags-1)                                      ! file naming convention
          nfiles=nfiles+1
        endif
      enddo

      !print*,'allfiles: ',allfiles

      ! check for multiple file submissions - generate queue of files if applicable
      partfile = trim(allfiles(0))
      idx = index(partfile,'.nc',back=.true.)
      idx = index(partfile(1:idx-1),'.',back=.true.)

      ! starting file
      file = partfile(1:idx-1)
      
      ! search for other files that wished to be partitioned - recording total amount
      total_files=1
      if (.not.(allocated(file_indices))) allocate(file_indices(nfiles))
      file_indices=-1  ! initialize with -1 as a form of mask
      file_indices(1)=0
      count=2          ! count to add properly file indices
      do i=0,nfiles-1  ! zero-based
        if (index(allfiles(i),trim(file))==0) then
          ! save indices where this file ends and next begins
          file_indices(count)=i-1
          file_indices(count+1)=i
          count=count+2  ! advance by two, such that each unique file has two indices (start and end)
          ! move on to search for next unique set of files to be joined
          newfile = allfiles(i)
          idx = index(newfile ,'.nc',back=.true.)
          idx = index(newfile(1:idx-1),'.',back=.true.)
          file = newfile(1:idx-1)   ! fill file with new name to update search
          total_files=total_files+1 ! and record presence of additional file
        endif
        ! close file_indices
        if (i==(nfiles-1)) file_indices(count)=i
      enddo

      !print*,'total files: ',total_files

      file_indices=pack(file_indices, file_indices>-1)

      !print*,'file_indices: ',(file_indices)

!      do i=1,size(file_indices)
!        if (MOD(i,2)/=0) then
!          stridx=file_indices(i); endidx=file_indices(i+1)
!          allocate(partfiles(endidx-stridx+1))
!          print*,'number of files: ', size(partfiles)
!          partfiles=allfiles(stridx:endidx)
!          print*,'File Ranges: ',partfiles
!          deallocate(partfiles)
!        endif
!      enddo

!      if (total_files==1) then
!        nparts = nfiles
!        part0 = 0
!        part1 = nparts-1
!      else
!      endif

      end subroutine cmd_args_ncjoin !]
! ----------------------------------------------------------------------
      subroutine show_use_map  ![
      implicit none

       print *, 'Usage of partmap should be: ',
     &       'partmap npx npy grdfile.nc'
     
      stop

      end subroutine show_use_map  !] 
! ----------------------------------------------------------------------
      subroutine show_use_partit  ![
      implicit none

      write(*,'(/1x,A//10x,2A)') 'Correct usage of partit:',
     &      'mpiexec  -n  np  partit  ',
     &      'np_x  np_y  grdfile.nc file1.nc...'
          write(*,'(1x,A/10x,2A/)') 'or',
     &      'mpiexec  -n  np  partit  ',
     &      'np_x  np_y --classic  grdfile.nc file1.nc...'
          write (*,'(/1x,2A/1x,2A//1x,2A/1x,A//1x,A//1x,A//)')
     &      'np_x and np_y are the ',
     &      'number of sub-domains you choose in x and y...',
     &      '-classic will include masked domains in your partition'

     
      stop

      end subroutine show_use_partit !]
! ----------------------------------------------------------------------
      subroutine show_use_ncjoin  ![
      implicit none

      write(*,'(/1x,A//10x,2A)') 'Correct usage of ncjoin:',
     &      'mpiexec  -n  np  ncjoin file1.nc...  '
          write(*,'(1x,A/10x,2A/)') 'or',
     &      'mpiexec  -n  np  ncjoin',
     &      '--delete file1.nc...'
          write (*,'(/1x,2A/1x,2A//1x,2A/1x,A//1x,A//1x,A//)')
     &      '-delete will remove all partital files once joined'

      stop

      end subroutine show_use_ncjoin !]
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

      ! Classic modes presence means all masks are set to 1 - cannot be ignored
      if (classic_part) mask = 1

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
        print*,jloc(npj,2),jloc(npj,3)
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
          if (msk_mx /= 0) then
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
      character(len=50) :: vname
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
        if (ierr/=nf90_noerr) 
     &    call handle_ierr(ierr,'creating file: ',trim(fname) )

        ! global rho-dimensions
        ierr=nf90_put_att(npid, nf90_global, 'globaldims', 
     &                   (/ncinfo%gsize(1)-2,ncinfo%gsize(2)-2/))

        ! NOTE: These are only saved when not in classical format such that
        ! ncjoin and roms_mpi.F are properly informed of how files were partitioned
        if (.not.classic_part) then
        ! i and j nodes (how it was partitioned)
        ierr=nf90_put_att(npid, nf90_global, 'ij_nodes', (/npartx,nparty/))
        ! subdompos: i,j position of sub-domain
        subdompos = (/ i-1, j-1 /) ! zero-based!
        ierr=nf90_put_att(npid, nf90_global, 'subdompos', subdompos)
        ! neighbors: N,NE,E,SE,S,SW,W,NW (READ IN THAT ORDER)
        loc_neigh = (/ node_map(j+1,i), node_map(j+1,i+1),
     &                 node_map(j,i+1), node_map(j-1,i+1),
     &                 node_map(j-1,i), node_map(j-1,i-1),
     &                 node_map(j,i-1), node_map(j+1,i-1) /)
        ierr=nf90_put_att(npid, nf90_global, 'neighbors',loc_neigh-1)
        endif

        ! partition: node, total nodes, llc x, llc y
        ibuff = (/part-1, nparts, iloc(i,2), jloc(j,2) /)
        ierr=nf90_put_att(npid, nf90_global, 'partition', ibuff)

        ! starting indices: i0, j0, i0u, j0v
        ibuff = (/iloc(i,2),jloc(j,2),ilcu(i,2),jlcv(j,2) /)
        ierr=nf90_put_att(npid, nf90_global, 'start_idx', ibuff)

        ! ending indices: i1, j1, i1u, j1v
        ibuff = (/iloc(i,3),jloc(j,3),ilcu(i,3),jlcv(j,3) /)
        ierr=nf90_put_att(npid, nf90_global, 'end_idx', ibuff)

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
      character(len=50) :: vname
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

        if (mynode==0) then
          if (record==1) print *,'Processing variable: ',vname
        endif

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
              nz = ncinfo%dims(ivar)%dimsizes(2)
              if ((bry=='w').or.(bry=='e')) then
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
              case('w')
                allocate(var3d(li0:li1,lj0:lj1,nz))
                startr = (/li0,lj0,1,1/)
                is = i0; ie = i1
                js = j0; je = j1
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
            case('w')
              allocate(var3d(li0:li1,lj0:lj1,nz))
              startr = (/li0,lj0,1,record/)
              is = i0; ie = i1
              js = j0; je = j1
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
      integer,dimension(2) :: ij
      type(ncstructure)  :: ncinfo
      character(len=50) :: vname
      character(len=7),dimension(5)  :: dname
      integer         ,dimension(5)  :: dsize
      integer         ,dimension(5)  :: did
      integer         ,dimension(5)  :: csizes
      character(len=20),dimension(5) :: attname
      character(len=200),dimension(8):: attvalue
      real,dimension(10)             ::attparam
      integer :: ierr, ncid,npid
      integer :: natts,idm,vid
      integer :: idx1
      integer,dimension(3)  ::gdims
      real,dimension(:,:,:),allocatable :: var3d
      integer,dimension(4) :: startr,ibuff

      partfile = trim(partfiles(part0))
      idx1 = index(partfile,'.nc',back=.true.)
      joinfile = partfile(1:idx1-1)
      idx1 = index(joinfile,'.',back=.true.)
      joinfile = joinfile(1:idx1-1) //'.nc'

      if (mynode==0) then
        print *,'Joining:';
        print *, trim(joinfile)
      endif

      ! Open first partial file in range 
      ierr = nf90_open(trim(partfile), nf90_nowrite, npid)
      if (ierr/=nf90_noerr) 
     &    call handle_ierr(ierr,'opening file: ',trim(partfile) )

      call get_nc_info(npid,ncinfo)

      ! get global dimensions - look for attribute
      ierr = nf90_get_att(npid,nf90_global,'globaldims',gdims)
      ! if file has attribute, mpi masking is present
      if (ierr==nf90_noerr) then
        ierr = nf90_open(trim(partfile), nf90_nowrite, npid)
        if (ierr/=nf90_noerr)
     &    call handle_ierr(ierr,'opening file: ',trim(partfile) )
        ierr=nf90_inquire(npid,ndimensions=ndims)
        do i = 1,ndims
          ierr = nf90_inquire_dimension(npid,i,dname(i),dsize(i))
          if (index(dname(i),'s_rho')>0)   gdims(3) = dsize(i)
        enddo
      ! if not present, then mpi_masking is not present, so simply grab
      ! final partial file to calculate global domains
      elseif (ierr/=nf90_noerr) then
        print*,'WARNING: globaldims attribute not found'
        print*,'Assuming final partial file is the uppermost right part'
        partfile = partfiles(nparts-1)
        print*,'final partfile: ',partfile
        ierr = nf90_open(trim(partfile), nf90_nowrite, npid)
        if (ierr/=nf90_noerr) 
     &    call handle_ierr(ierr,'opening file: ',trim(partfile) )
        ierr=nf90_get_att(npid, nf90_global, 'partition', ibuff)
        ierr=nf90_inquire(npid,ndimensions=ndims)
        do i = 1,ndims
          ierr = nf90_inquire_dimension(npid,i,dname(i),dsize(i))
          if (index(dname(i),'xi_rho')>0)  gdims(1) = ibuff(3) + dsize(i) -3
          if (index(dname(i),'eta_rho')>0) gdims(2) = ibuff(4) + dsize(i) -3
          if (index(dname(i),'s_rho')>0)   gdims(3) = dsize(i)
        enddo
        ierr = nf90_close(npid)
      endif

      gnx = gdims(1)+2
      gny = gdims(2)+2
      nz  = gdims(3)
      
      if (mynode==0) then ! single node to create joined file

        ierr=nf90_create(trim(joinfile),nf90_netcdf4,ncid)
        if (ierr/=nf90_noerr) 
     &    call handle_ierr(ierr,'creating file: ',trim(joinfile) )

        ! Adding dmimensions to joined file
        ierr = nf90_open(trim(partfiles(nparts-1)), nf90_nowrite, npid)
        if (ierr/=nf90_noerr) 
     &    call handle_ierr(ierr,'opening file: ',trim(partfile) )
        ierr=nf90_inquire(npid,ndimensions=ndims)
        do i = 1,ndims
          ierr = nf90_inquire_dimension(npid,i,dname(i),dsize(i))
          if (index(dname(i),'xi_rho')>0)   dsize(i) = gnx
          if (index(dname(i),'xi_u') > 0) dsize(i) = gnx-1
          if (index(dname(i),'eta_rho')>0)  dsize(i) = gny
          if (index(dname(i),'eta_v') > 0)dsize(i) = gny-1
          if (index(dname(i),'time') > 0) dsize(i) = 0 
          if (index(dname(i),'s_rho') > 0) dsize(i) = nz
          if (index(dname(i),'s_w') > 0) dsize(i) = nz+1
          ierr=nf90_def_dim(ncid,dname(i),dsize(i),did(i))
        enddo

        ! Read all global attributes from partial file and write to joined
        ierr = nf90_inquire(npid, nAttributes=ngatt)
        if (ierr/=nf90_noerr) 
     &    call handle_ierr(ierr,'getting global atts: ',trim(partfile) )
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
          natts   = ncinfo%natts(ivar)
          attname(1:natts) = ncinfo%atts(ivar)%attnames
          attvalue(1:natts)= ncinfo%atts(ivar)%attvalues
          attparam(1:natts)= ncinfo%atts(ivar)%attparams
          dsize(1:ndims) = ncinfo%dims(ivar)%dimsizes
          dname(1:ndims) = ncinfo%dims(ivar)%dimnames

          if (xtype.eq.nf90_char) goto 98

          do idm =1,ndims
            if (index(dname(idm),'xi_rho') > 0)  dsize(idm)=gnx
            if (index(dname(idm),'xi_u') > 0)    dsize(idm)=gnx-1
            if (index(dname(idm),'eta_rho') > 0) dsize(idm)=gny
            if (index(dname(idm),'eta_v') > 0)   dsize(idm)=gny-1
            if (index(dname(idm),'time') > 0)    dsize(idm)=0
            ierr = nf90_inq_dimid(ncid,dname(idm),did(idm))
          enddo

          if (ndims<=1) goto 99

          ierr=nf90_def_var(ncid,vname,xtype,did(1:ndims),vid,
     &            chunksizes = csizes(1:ndims),
     &            deflate_level=1, shuffle=.true.)
          if (ierr/=nf90_noerr) 
     &      call handle_ierr(ierr,'creating joined var: ',trim(vname) )
 

  98      if (ncinfo%xtype(ivar).eq.nf90_char) then
            ierr=nf90_def_var(ncid,vname,xtype,vid)
            if (ierr/=nf90_noerr) 
     &       call handle_ierr(ierr,'creating joined var: ',trim(vname) )
          endif

  99      if (((ncinfo%xtype(ivar).eq.nf90_double).and.(ndims<=1))
     &       .or.((ncinfo%xtype(ivar).eq.nf90_float).and.(ndims<=1)))then
            if (ndims==0) then
              ierr=nf90_def_var(ncid,vname,xtype,vid)
            elseif (ndims==1) then
              ierr=nf90_def_var(ncid,vname,xtype,did(1:ndims),vid)
            endif
            if (ierr/=nf90_noerr) 
     &       call handle_ierr(ierr,'creating joined var: ',trim(vname) )
          endif

          ierr = nf90_var_par_access(ncid, vid, nf90_collective)
          do iatt=1,natts
            ! input attribute based on type
            if (attvalue(iatt)=='n/a') then  ! parameter (type(real))
              ierr = nf90_put_att(ncid,vid,attname(iatt),attparam(iatt))
            elseif (attparam(iatt)==0) then  ! string    (type(char))
              ierr = nf90_put_att(ncid,vid,attname(iatt),attvalue(iatt))
            endif
          enddo

        enddo

        ierr = nf90_close(ncid)
        ierr = nf90_close(npid)

      endif ! mynode==0

  
      end subroutine make_joined_file !]
! ----------------------------------------------------------------------
      subroutine write_joined_file ![
      ! read data for each variable in partial files
      ! and write to joinedl file, for all records
      implicit none

      ! local
      integer :: ierr, ncid,npid,dimid
      integer,dimension(4) :: llidx,ibuff
      integer :: nx,ny,nxu=0,nyv=0
      integer,dimension(:),allocatable :: i0,j0
      integer,dimension(:),allocatable :: i1,j1
      integer,dimension(:),allocatable :: i0u,i1u
      integer,dimension(:),allocatable :: j0v,j1v
      integer :: lnx,lny,li0,li1,lj0,lj1,li0u,li1u,lj0v,lj1v
      real,dimension(:,:,:),allocatable :: var3d,tmp3d
      character(len=50) :: vname
      integer,dimension(:),allocatable :: is,ie,js,je
      type(ncstructure)  :: ncinfo
      character(len=7) :: ldname
      integer :: ldsize

      integer :: varid,ivar
      integer :: record,nrecords
      integer,dimension(4) :: startr,startw
      real,dimension(:),allocatable :: var1d
      integer :: sz1,sz2,sz3,ndims
      integer :: RSS,i

      integer,dimension(:),allocatable :: jparts
      integer                          :: jsize
      integer                          :: jrow,nrows,row,rowpart,rp0,rp1
      integer                          :: inode,jnode
      integer,dimension(:),allocatable :: srow
      integer,dimension(:),allocatable :: partitions
      integer,dimension(3)             :: csizes
      integer,dimension(2)             :: ij
      integer                          :: mpi_role,write_comm,counter,total_nodes
      character(len=1)                 :: bryvariable
      logical                          :: newfile

      if (mynode==0) then
      ! get number of rows (jnodes)
      ierr = nf90_open(trim(partfiles(part0)),nf90_nowrite,npid)
      ierr = nf90_get_att(npid,nf90_global,'ij_nodes',ij)
      ! first check if provided in attributes
      if (ierr==nf90_noerr) then
        ierr=nf90_close(npid)  ! file closed within conditionals as prev. ierr is needed
        nrows=ij(2)
        ! make note of mpi_masking:
        newfile=.true.
      elseif (ierr/=nf90_noerr) then
        newfile=.false.
        ierr=nf90_close(npid)
        allocate(j0(part0:part1))
        allocate(j1(part0:part1))
        allocate(jparts(part0:part1))
        nrows = 1
        do part=part0,part1
          partfile=partfiles(part)
          ierr = nf90_open(trim(partfile),nf90_nowrite,npid)
          if (ierr/=nf90_noerr) 
     &      call handle_ierr(ierr,'opening file: ',trim(partfile) )
          ierr = nf90_get_att(npid,nf90_global,'partition',ibuff)
          j0(part) = ibuff(4)
          ierr=nf90_inquire(npid,ndimensions=ndims)
          do i = 1,ndims
            ierr = nf90_inquire_dimension(npid,i,ldname,ldsize)
            if (index(ldname,'eta_rho')>0) ny = ldsize
          enddo
          ! get nrows(jnodes)
          j1(part) = j0(part) + ny-1
          if (part==0) jrow=j1(part)
          if (j1(part)/=jrow) then
            nrows=nrows+1; jrow=j1(part)
          endif
          jparts(part)=nrows
          ierr=nf90_close(npid)
        enddo
        deallocate(j0)
        deallocate(j1)
        jsize=size(jparts)
      endif  ! attribute version
      endif  ! mynode==0

      call MPI_BCAST(newfile,1,MPI_LOGICAL,0,part_comm,ierr)

      ! broadcast size of jparts (only needed for no mpi masking)
      if (.not.newfile) call MPI_BCAST(jsize,1,MPI_INTEGER,0,part_comm,ierr)
      if (mynode/=0) then
        allocate(jparts(0:jsize))
      endif
      if (.not.newfile) call MPI_BCAST(jparts,jsize,MPI_INTEGER,0,part_comm,ierr)

      call MPI_BCAST(nrows,1,MPI_INTEGER,0,part_comm,ierr)
 
      call MPI_Barrier(part_comm,ierr)

      ! get j partitions per row
      rowpart = ceiling(1.0*nrows/nnodes)
      !if (mynode==0) print*,'Jnodes per MPI Row: ',rowpart

      if ((mynode+1)*rowpart<=nrows) then
        mpi_role=1
      else
        mpi_role=MPI_UNDEFINED
      endif

      call MPI_Comm_split(part_comm,mpi_role,mynode,write_comm,ierr)

      ! only the mpi nodes needed will continue
      if (mpi_role==1) then

      call MPI_Comm_size(write_comm,total_nodes,ierr)

      ! on each node, grab partial files in its row
      allocate(partitions(nparts))
      partitions = -1
      counter=1
      do part=part0,part1
        ! get jnodes
        if (newfile) then
          ierr=nf90_open_par(trim(partfiles(part)),nf90_nowrite,write_comm,MPI_Info_Null,npid)
          ierr = nf90_get_att(npid,nf90_global,'subdompos',ij)
          jnode=ij(2)
        else
          jnode=jparts(part)-1
        endif
        ! check if current partitioned files jnode is within desired mpinode(row)
        if ((jnode >= (mynode*rowpart)).and.(jnode < ((mynode+1)*rowpart))) then
          partitions(counter)=part
          counter=counter+1  
        ! the final row, which is last mpi node, will overflow with remainding nodes if needed
        elseif (((mynode+1)==total_nodes).and.(jnode > (mynode)*rowpart)) then
          partitions(counter)=part
          counter=counter+1  
        endif
        ierr=nf90_close(npid)
      enddo

      ! partition ranges
      partitions=pack(partitions, partitions>-1)
      rp0 = partitions(1)
      rp1 = partitions(size(partitions))
  
      call MPI_Barrier(write_comm,ierr)

      ! Read partial files in range
      allocate(i0(rp0:rp1))
      allocate(j0(rp0:rp1))
      allocate(i1(rp0:rp1))
      allocate(j1(rp0:rp1))
      allocate(i0u(rp0:rp1))
      allocate(i1u(rp0:rp1))
      allocate(j0v(rp0:rp1))
      allocate(j1v(rp0:rp1))

      allocate(is(rp0:rp1))
      allocate(js(rp0:rp1))
      allocate(ie(rp0:rp1))
      allocate(je(rp0:rp1))

      !if (mynode==0) print *,'scanning partial files'
      do part = rp0,rp1 
    
        partfile = partfiles(part)

        ierr = nf90_open(trim(partfile),nf90_nowrite,npid)

        if (part.eq.rp0) call get_nc_info(npid,ncinfo)

        if (ierr/=nf90_noerr) 
     &    call handle_ierr(ierr,'opening file: ',trim(partfile) )

        ierr = nf90_get_att(npid,nf90_global,'partition',ibuff)
        i0(part) = ibuff(3)
        j0(part) = ibuff(4)
        i0u(part)= max(i0(part)-1,1)
        j0v(part)= max(j0(part)-1,1)

        ierr=nf90_inquire(npid,ndimensions=ndims)
        do i = 1,ndims
          ierr = nf90_inquire_dimension(npid,i,ldname,ldsize)
          if (index(ldname,'xi_rho')>0)  nx = ldsize
          if (index(ldname,'eta_rho')>0) ny = ldsize
          if (index(ldname,'xi_u')>0)    nxu = ldsize
          if (index(ldname,'eta_v')>0)   nyv = ldsize
        enddo
        i1(part) = i0(part) + nx-1
        j1(part) = j0(part) + ny-1
        i1u(part)= i0u(part)+ nxu-1
        j1v(part)= j0v(part)+ nyv-1

        ierr = nf90_close(npid)
      enddo

      call MPI_Barrier(write_comm,ierr)

      if (mynode==0) print *,'Total records: ',ncinfo%nrecords

      do record = 1,ncinfo%nrecords

        if (mynode==0.and.ncinfo%nrecords/=1) then
          print *,'Processing record: ',record,' of', ncinfo%nrecords
        endif

        do ivar=1,ncinfo%nvars

          vname = ncinfo%varname(ivar)

          ! check if bry var
          bryvariable=ncinfo%bryvar(ivar)

          if (mynode==0) print *,'Processing: ',trim(vname)

          ndims = ncinfo%ndims(ivar)

        if (ncinfo%varpart(ivar)) then

          if (ndims>=2) then
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
          endif

          select case(ndims)
          case(1)
            ! save for later
          case(2)
            if (index(ncinfo%dims(ivar)%dimnames(2),'time')>0) then
              select case(ncinfo%vartype(ivar))
              case('r')
                if ((bryvariable=='w').or.(bryvariable=='e')) then
                  allocate(var3d(lj0:lj1,1,1)); allocate(tmp3d(lj0:lj1,1,1))
                  is = j0; ie = j1
                  startw = (/lj0,record,1,1/)
                elseif ((bryvariable=='s').or.(bryvariable=='n')) then
                  allocate(var3d(gnx,1,1)); allocate(tmp3d(gnx,1,1))
                  is = i0; ie = i1
                  startw = (/1,record,1,1/)
                endif
              end select
              startr = (/1,record,1,1/)
            else
              select case(ncinfo%vartype(ivar))
              case('r')
                allocate(var3d(gnx,lj0:lj1,1))
                is = i0; ie = i1
                js = j0; je = j1
                startw = (/1,lj0,record,1/)
              case('u')
                allocate(var3d(gnx-1,lj0:lj1,1))
                is = i0; ie = i1
                js = j0; je = j1
                startw = (/1,lj0,record,1/)
              case('v')
                allocate(var3d(gnx,lj0v:lj1v,1))
                is = i0; ie = i1
                js = j0; je = j1
                startw = (/1,lj0,record,1/)
              case('j')
                print*,vname,'psi variables skipped'; cycle
              end select
              startr = (/1,1,record,1/)
            endif
          case(3)
            if (index(ncinfo%dims(ivar)%dimnames(3),'time')>0) then
              select case(ncinfo%vartype(ivar))
              case('r')
                allocate(var3d(gnx,lj0:lj1,1))
                is = i0; ie = i1
                js = j0; je = j1
                startw = (/1,lj0,record,1/)
              case('u')
                allocate(var3d(gnx-1,lj0:lj1,1))
                is = i0u; ie = i1u
                js = j0; je = j1
                startw = (/1,lj0,record,1/)
              case('v')
                allocate(var3d(gnx,lj0v:lj1v,1))
                is = i0; ie = i1
                js = j0v; je = j1v
                startw = (/1,lj0v,record,1/)
              end select
              startr = (/1,1,record,1/)
            else
              if (record>1) continue
              select case(ncinfo%vartype(ivar))
              case('r')
                allocate(var3d(gnx,lj0:lj1,nz))
                is = i0; ie = i1
                js = j0; je = j1
                startw = (/1,lj0,1,record/)
              case('u')
                allocate(var3d(gnx-1,lj0:lj1,nz))
                is = i0u; ie = i1u
                js = j0; je = j1
                startw = (/1,lj0,1,record/)
              case('v')
                allocate(var3d(gnx,lj0v:lj1v,nz))
                is = i0; ie = i1
                js = j0v; je = j1v
                startw = (/1,lj0v,1,record/)
              case('w') ! vertical components
                allocate(var3d(gnx,lj0:lj1,nz+1))
                is = i0; ie = i1
                js = j0; je = j1
                startw = (/1,lj0,1,record/)
              end select
              startr = (/1,1,1,record/)
            ! save for later (not present in pachug_rst.nc)
            endif
          case(4)
              select case(ncinfo%vartype(ivar))
              case('r')
                allocate(var3d(gnx,lj0:lj1,nz))
                is = i0; ie = i1
                js = j0; je = j1
                startw = (/1,lj0,1,record/)
              case('u')
                allocate(var3d(gnx-1,lj0:lj1,nz))
                is = i0u; ie = i1u
                js = j0; je = j1
                startw = (/1,lj0,1,record/)
              case('v')
                allocate(var3d(gnx,lj0v:lj1v,nz))
                is = i0; ie = i1
                js = j0v; je = j1v
                startw = (/1,lj0v,1,record/)
              case('w') ! vertical components
                allocate(var3d(gnx,lj0:lj1,nz+1))
                is = i0; ie = i1
                js = j0; je = j1
                startw = (/1,lj0,1,record/)
              end select
              startr = (/1,1,1,record/)
          end select

          var3d = 0
          ! ----- Boundary Variables ----- ![
          ! conditional splits up loops, as boundary variables are dependent on current parition
          ! unlike other partitionable variables. Conditional outside of loop results in two
          ! seperate ones, but is slighlty faster if a large amount of partitions is present.
          if (bryvariable/='i') then
          tmp3d = 0
          do part = rp0,rp1
            partfile = partfiles(part)
            ierr = nf90_open(trim(partfile),nf90_nowrite,npid)
            call ncread(npid,vname,tmp3d(is(part):ie(part),:,:),startr)
            ! since var3d is added to for ncwrite, we use a temporary array to see if this 
            ! bry variable is present in current partition. If not, then we skip and do not
            ! overwrite var3d with empty data.
            if (minval(tmp3d(is(part):ie(part),:,:))==maxval(tmp3d(is(part):ie(part),:,:))) cycle
            var3d(is(part):ie(part),:,:)=tmp3d(is(part):ie(part),:,:)
          enddo
          else
          ! ------------------------------ !]
          do part = rp0,rp1
            partfile = partfiles(part)
            ierr = nf90_open(trim(partfile),nf90_nowrite,npid)
            call ncread(npid,vname,var3d(is(part):ie(part),js(part):je(part),:),startr)
          enddo
          endif

          call MPI_Barrier(write_comm,ierr)

          ierr = nf90_open_par(trim(joinfile),nf90_write,write_comm,MPI_Info_Null,ncid)
          if (ierr/=nf90_noerr) 
     &      call handle_ierr(ierr,'opening file: ',trim(joinfile) )
          ierr = nf90_var_par_access(ncid,ivar,nf90_collective)

          call ncwrite(ncid,vname,var3d,startw)

          ierr = nf90_close(ncid)

          if (allocated(var3d)) deallocate(var3d)
          if (allocated(tmp3d)) deallocate(tmp3d)

        else !not-partitionable

          if (ncinfo%xtype(ivar).eq.nf90_char) cycle
 
          if (mynode==0) then
            select case(ndims)
            case(0)
              cycle
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

            var3d = 0
            partfile = partfiles(rp0)
            ierr = nf90_open(trim(partfile), nf90_nowrite, npid)
            if (ierr/=nf90_noerr) 
     &        call handle_ierr(ierr,'opening file: ',trim(partfile))
            call ncread(npid,trim(vname),var3d(is(rp0):ie(rp0),js(rp0):je(rp0),:),startr)
            ierr = nf90_close(npid)

            ierr = nf90_open(trim(joinfile),nf90_write,ncid)
            if (ierr/=nf90_noerr) 
     &        call handle_ierr(ierr,'opening file: ',trim(joinfile) )
            call ncwrite(ncid,vname,var3d,startw)
            ierr = nf90_close(ncid)

            if (allocated(var3d)) deallocate(var3d)

          endif !mynode=0

        endif

      enddo  ! ivars

      enddo ! records

      deallocate(i0)
      deallocate(i1)
      deallocate(j0)
      deallocate(j1)

      endif  ! <-- mpi conditional


      end subroutine write_joined_file !]
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
      character(len=50) :: vname
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
            ncinfo%vartype(ivar) = 'w'
!            ncinfo%gsize(3) = dsize
          endif
          if (index(dname,'time')>0) then
            ncinfo%nrecords = dsize
          else
            ncinfo%nrecords = 1
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
