      module nc_read_write

#include "cppdefs.opt"

      use netcdf

      implicit none
      private

#include "nc_read_write.opt"

      interface ncwrite
        module procedure  ncwrite_1D, ncwrite_2D, ncwrite_3D
      end interface

      interface ncread
        module procedure  ncread_1D, ncread_2D, ncread_3D
      end interface

      public :: ncread
      public :: ncwrite
      public :: nccreate
      public :: handle_ierr

      contains

!----------------------------------------------------------------------
      subroutine ncread_1D(ncid,vname,dat,start)  ![
      implicit none

      ! input
      integer           :: ncid
      character(len=*)  :: vname
      real,dimension(:) :: dat

      integer,dimension(:),optional :: start

      ! local
      integer              :: varid,ierr
      integer,dimension(1) :: dims
      integer,dimension(2) :: count

      dims = shape(dat)
      ierr = nf90_inq_varid(ncid,vname,varid)
      if (ierr/=0) call handle_ierr(ierr,'ncread inq_varid var=',vname)

      if (present(start)) then
        count = (/dims, 1/)
        ierr = nf90_get_var(ncid,varid,dat,start,count)
      else
        ierr = nf90_get_var(ncid,varid,dat)
      endif
      if (ierr/=0) call handle_ierr(ierr,'ncread read var=',vname)

      end subroutine ncread_1D  !]
!----------------------------------------------------------------------
      subroutine ncread_2D(ncid,vname,dat,start)  ![
      implicit none

      ! input
      integer             :: ncid
      character(len=*)    :: vname
      real,dimension(:,:) :: dat

      integer,dimension(:),optional :: start

      ! local
      integer              :: varid,ierr
      integer,dimension(2) :: dims
      integer,dimension(3) :: count

      dims = shape(dat)
      ierr = nf90_inq_varid(ncid,vname,varid)
      if (ierr/=0) call handle_ierr(ierr,'ncread inq_varid var=',vname)

      if (present(start)) then
        count = (/dims, 1/)
        ierr = nf90_get_var(ncid,varid,dat,start,count)
      else
        ierr = nf90_get_var(ncid,varid,dat)
      endif
      if (ierr/=0) call handle_ierr(ierr,'ncread read var=',vname)

      end subroutine ncread_2D  !]
!----------------------------------------------------------------------
      subroutine ncread_3D(ncid,vname,dat,start)  ![
      implicit none

      ! input
      integer               :: ncid
      character(len=*)      :: vname
      real,dimension(:,:,:) :: dat

      integer,dimension(:),optional :: start

      ! local
      integer              :: varid,ierr
      integer,dimension(3) :: dims
      integer,dimension(4) :: count

      dims = shape(dat)
      ierr = nf90_inq_varid(ncid,vname,varid)
      if (ierr/=0) call handle_ierr(ierr,'ncread inq_varid var=',vname)

      if (present(start)) then
        count = (/dims, 1/)
        ierr = nf90_get_var(ncid,varid,dat,start,count)
      else
        ierr = nf90_get_var(ncid,varid,dat)
      endif
      if (ierr/=0) call handle_ierr(ierr,'ncread read var=',vname)

      end subroutine ncread_3D  !]
!----------------------------------------------------------------------
      subroutine ncwrite_1D(ncid,vname,dat,start)  ![
      implicit none

      ! input
      integer           :: ncid
      character(len=*)  :: vname
      real,dimension(:) :: dat

      integer,dimension(:),optional :: start

      ! local
      integer              :: varid,ierr
      integer,dimension(1) :: dims
      integer,dimension(2) :: count

      dims = shape(dat)
      ierr = nf90_inq_varid(ncid,vname,varid)
      if (ierr/=0) call handle_ierr(ierr,'ncwrite inq_varid var=',vname)

      if (present(start)) then
        count = (/dims, 1/)
        ierr = nf90_put_var(ncid,varid,dat,start,count)
      else
        ierr = nf90_put_var(ncid,varid,dat)
      endif
      if (ierr/=0) call handle_ierr(ierr,'ncwrite write var=',vname)

      end subroutine ncwrite_1D  !]
!----------------------------------------------------------------------
      subroutine ncwrite_2D(ncid,vname,dat,start)  ![
      implicit none

      ! input
      integer             :: ncid
      character(len=*)    :: vname
      real,dimension(:,:) :: dat

      integer,dimension(:),optional :: start

      ! local
      integer              :: varid,ierr
      integer,dimension(2) :: dims
      integer,dimension(3) :: count

      dims = shape(dat)
      ierr = nf90_inq_varid(ncid,vname,varid)
      if (ierr/=0) call handle_ierr(ierr,'ncwrite inq_varid var=',vname)

      if (present(start)) then
        count = (/dims, 1/)
        ierr = nf90_put_var(ncid,varid,dat,start,count)
      else
        ierr = nf90_put_var(ncid,varid,dat)
      endif
      if (ierr/=0) call handle_ierr(ierr,'ncwrite write var=',vname)

      end subroutine ncwrite_2D  !]
!----------------------------------------------------------------------
      subroutine ncwrite_3D(ncid,vname,dat,start)  ![
      implicit none

      ! input
      integer               :: ncid
      character(len=*)      :: vname
      real,dimension(:,:,:) :: dat

      integer,dimension(:),optional :: start

      ! local
      integer              :: varid,ierr
      integer,dimension(3) :: dims
      integer,dimension(4) :: count

      dims = shape(dat)
      ierr = nf90_inq_varid(ncid,vname,varid)
      if (ierr/=0) call handle_ierr(ierr,'ncwrite inq_varid var=',vname)

      if (present(start)) then
        count = (/dims, 1/)
        ierr = nf90_put_var(ncid,varid,dat,start,count)
      else
        ierr = nf90_put_var(ncid,varid,dat)
      endif
      if (ierr/=0) call handle_ierr(ierr,'ncwrite write var=',vname)

      end subroutine ncwrite_3D  !]
!----------------------------------------------------------------------
      integer function nccreate(ncid,varname,dimname,dimsize,vartype)  ![
      ! Create a variable with dimensions in an existing file
      ! Returns the varid of the variable

      ! ncid:    ID of an opened netcdf file
      ! varname: name of the variable
      ! dimname: Names of dimensions of the variable
      ! dimsize: Dimension length (used if dimension is not yet defined)
      ! vartype: netcdf data type
      implicit none

      ! import/export
      integer,                      intent(in) :: ncid
      character(len=*),             intent(in) :: varname
      character(len=*),dimension(:),intent(in) :: dimname
      integer,dimension(:),optional,intent(in) :: dimsize
      integer,optional,             intent(in) :: vartype
      ! local
      integer :: i,ndim,varid,ierr,did,xtype
      integer,allocatable,dimension(:) :: dimid

      if (present(vartype)) then                           ! handle optional arguement
        xtype = vartype
      else
        xtype = default_prec 
      endif

      ndim = size(dimname)
      allocate(dimid(ndim))

      do i = 1,ndim                                        ! get dimension ids. Create if needed.
        ierr = nf90_inq_dimid(ncid,dimname(i),did)
        if (ierr/=nf90_noerr) then
          if (.not. present(dimsize)) then                 ! only an issue if dimension doesn't exist yet
            print *, 'ERROR: nccreate, no dimsize given for var=',varname
            error stop
          endif
          ierr=nf90_def_dim(ncid,dimname(i),dimsize(i),did)
          if (ierr/=0) call handle_ierr(ierr,'nccreate :: def_dim var=',varname)
        endif
        dimid(i) = did
      enddo

      ierr=nf90_def_var(ncid,varname,xtype,dimid,varid,
     &                  deflate_level=deflate_level, shuffle=shuffle)
      if (ierr/=nf90_noerr) then
        print *,'ERROR creating: ',varname
        call handle_ierr(ierr,'nccreate :: def_var=',varname)
      endif

      nccreate = varid

      end function nccreate  !]
! ---------------------------------------------------------------------
      subroutine handle_ierr(ierr, err_msg1, err_msg2)  ![
      ! handle fatal ierr with error messages.
      ! err_msg1 & 2 are optional messages
      ! get this from roms_read_write module rather at some point
      use param
      implicit none

      ! inputs
      integer ierr
      character(len=*), optional :: err_msg1, err_msg2
      ! local
      character(:), allocatable :: f_msg, f_msg1, f_msg2

      ! Check if error message 1 is available
      if(present(err_msg1))then
          f_msg1=err_msg1
      else
          f_msg1=''
      endif
      ! Check if error message 2 is available
      if(present(err_msg2))then
          f_msg2=err_msg2
      else
          f_msg2=''
      endif
      ! Combine optional messages:
      f_msg = trim(f_msg1) / / ' ' / / trim(f_msg2)

      write(*,'(/3x,A/5x,I4,A/12x,A/)') 'ERROR:: ',mynode,' ',
     &       trim(f_msg), nf90_strerror(ierr)

      error stop ! Stop simulation assumed fatal error

      end subroutine handle_ierr  !]
! ---------------------------------------------------------------------

      end module nc_read_write
