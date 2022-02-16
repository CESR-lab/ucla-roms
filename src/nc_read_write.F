      module nc_read_write

      use netcdf
!     use param

      implicit none
      private

      ! ****************************************************************
      ! START USER INPUT
      integer :: deflate_level=1 ! Compression level for joined file
      logical :: shuffle=.true.  ! Shuffle on for extra conmpression
      ! END USER INPUT
      ! ****************************************************************

      interface ncwrite
        module procedure  ncwrite_1D, ncwrite_2D, ncwrite_3D
      end interface

      interface ncread
        module procedure  ncread_1D, ncread_2D
      end interface

      public :: ncread
      public :: ncwrite
      public :: nccreate

      contains

!----------------------------------------------------------
      subroutine ncread_1D(ncid,vname,dat,start)
        implicit none

        ! Input
        integer           :: ncid
        character(len=*)  :: vname
        integer,dimension(:),optional :: start

        real,dimension(:) :: dat

        ! local
        integer              :: varid,ierr
        integer,dimension(1) :: dims
        integer,dimension(2) :: sizes

        dims = shape(dat)
        ierr = nf90_inq_varid(ncid,vname,varid)

        if (present(start)) then
          sizes = (/dims, 1/)
          ierr = nf90_get_var(ncid,varid,dat,start,sizes) 
        else
          ierr = nf90_get_var(ncid,varid,dat)
        endif

      end subroutine ncread_1D

!----------------------------------------------------------
      subroutine ncread_2D(ncid,vname,dat,start)
        implicit none

        ! Input
        integer           :: ncid
        character(len=*)  :: vname
        integer,dimension(:),optional :: start

        real,dimension(:,:) :: dat

        ! local
        integer              :: varid,ierr
        integer,dimension(2) :: dims
        integer,dimension(3) :: sizes

        dims = shape(dat)
        ierr = nf90_inq_varid(ncid,vname,varid)

        if (present(start)) then
          sizes = (/dims, 1/)
          ierr = nf90_get_var(ncid,varid,dat,start,sizes) 
        else
          ierr = nf90_get_var(ncid,varid,dat)
        endif

      end subroutine ncread_2D

!----------------------------------------------------------
      subroutine ncwrite_1D(ncid,vname,dat,start)
      implicit none

      ! Input
      integer          :: ncid
      character(len=*) :: vname
      integer,dimension(:),optional :: start

      real,dimension(:) :: dat

      ! local
      integer              :: varid,ierr
      integer,dimension(1) :: dims
      integer,dimension(2) :: sizes

      dims = shape(dat)
      ierr = nf90_inq_varid(ncid,vname,varid)
      if (ierr/=nf90_noerr) call handle_err(ierr)

      if (present(start)) then
        sizes = (/dims, 1/)
        ierr = nf90_put_var(ncid,varid,dat,start,sizes) 
      else
        ierr = nf90_put_var(ncid,varid,dat)
      endif
      if (ierr/=nf90_noerr) call handle_err(ierr)

      end subroutine ncwrite_1D

!----------------------------------------------------------
      subroutine ncwrite_2D(ncid,vname,dat,start)
      implicit none

      ! Input
      integer          :: ncid
      character(len=*) :: vname
      integer,dimension(:),optional :: start

      real,dimension(:,:) :: dat


      ! local
      integer              :: varid,ierr
      integer,dimension(2) :: dims
      integer,dimension(3) :: sizes

      dims = shape(dat)
      ierr = nf90_inq_varid(ncid,vname,varid)
      if (ierr/=nf90_noerr) then
        call handle_err(ierr,'ncwrite: '/ /vname)
      endif

      if (present(start)) then
        sizes = (/dims, 1/)
        ierr = nf90_put_var(ncid,varid,dat,start,sizes) 
      else
        ierr = nf90_put_var(ncid,varid,dat)
      endif
      if (ierr/=nf90_noerr) call handle_err(ierr)

      end subroutine ncwrite_2D

!----------------------------------------------------------
      subroutine ncwrite_3D(ncid,vname,dat,start)
      implicit none

      ! Input
      integer          :: ncid
      character(len=*) :: vname
      integer,dimension(:),optional :: start

      real,dimension(:,:,:) :: dat

      ! local
      integer              :: varid,ierr
      integer,dimension(3) :: dims
      integer,dimension(4) :: sizes

      dims = shape(dat)
      ierr = nf90_inq_varid(ncid,vname,varid)
      if (ierr/=nf90_noerr) then
        call handle_err(ierr,'ncwrite :: '/ /vname)
      endif

      if (present(start)) then
        sizes = (/dims, 1/)
        ierr = nf90_put_var(ncid,varid,dat,start,sizes) 
      else
        ierr = nf90_put_var(ncid,varid,dat)
      endif

      if (ierr/=nf90_noerr) then
        call handle_err(ierr,'ncwrite :: '/ /vname)
      endif

      end subroutine ncwrite_3D

!----------------------------------------------------------
      integer function nccreate(ncid,varname,dimname,dimsize,vartype)
      ! Create a variable with dimensions in an existing file
      ! Returns the varid of the variable

      ! ncid:    ID of an opened netcdf file
      ! varname: name of the variable
      ! dimname: Names of dimensions of the variable
      ! dimsize: Dimension length (used if dimension is not yet defined)
      ! vartype: netcdf data type
      implicit none

      ! import/export
      integer, intent(in) :: ncid
      character(len=*),             intent(in) :: varname
      character(len=*),dimension(:),intent(in) :: dimname
      integer,         dimension(:),intent(in) :: dimsize
      integer,         optional,    intent(in) :: vartype 
      ! local
      integer :: i,ndim,varid,ierr,did
      integer,allocatable,dimension(:) :: dimid

      ndim = size(dimsize)
      allocate(dimid(ndim))

      ! Get dimension ids. Create if needed.
      do i = 1,ndim 
        ierr = nf90_inq_dimid(ncid,dimname(i),did)
        if (ierr.ne.nf90_noerr) then
          ierr=nf90_def_dim(ncid,dimname(i),dimsize(i),did)
        endif
        dimid(i) = did
      enddo
      ! Add compression and shuffle
      ierr=nf90_def_var(ncid,varname,vartype,dimid,varid)
      if (ierr/=nf90_noerr) call handle_err(ierr)

      nccreate = varid

      end function nccreate

!-------------------------------------------
      subroutine handle_err(status,txt)
      implicit none
      integer,         intent(in)          :: status
      character(len=*),intent(in),optional :: txt


      if(status /= nf90_noerr) then
        if (present(txt)) then
         print *, txt,', ',trim(nf90_strerror(status))
        else
         print *, trim(nf90_strerror(status))
        endif
         stop "Stopped"
      end if
      end subroutine handle_err


      end module nc_read_write
