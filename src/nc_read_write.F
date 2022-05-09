      module nc_read_write

      use netcdf
!     use param

      implicit none
      private

#include "nc_read_write.opt"

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
        real,dimension(:) :: dat

        integer,dimension(:),optional :: start

        ! local
        integer              :: varid,ierr
        integer,dimension(1) :: dims
        integer,dimension(2) :: count

        dims = shape(dat)
        ierr = nf90_inq_varid(ncid,vname,varid)

        if (present(start)) then
          count = (/dims, 1/)
          ierr = nf90_get_var(ncid,varid,dat,start,count)
        else
          ierr = nf90_get_var(ncid,varid,dat)
        endif

      end subroutine ncread_1D

!----------------------------------------------------------
      subroutine ncread_2D(ncid,vname,dat,start)
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

        if (present(start)) then
          count = (/dims, 1/)
          ierr = nf90_get_var(ncid,varid,dat,start,count)
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
      integer,dimension(2) :: count

      dims = shape(dat)
      ierr = nf90_inq_varid(ncid,vname,varid)
      if (ierr/=nf90_noerr) call handle_err(ierr)

      if (present(start)) then
        count = (/dims, 1/)
        ierr = nf90_put_var(ncid,varid,dat,start,count)
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
      integer,dimension(3) :: count

      dims = shape(dat)
      ierr = nf90_inq_varid(ncid,vname,varid)
      if (ierr/=nf90_noerr) then
        call handle_err(ierr,'ncwrite: '/ /vname)
      endif

      if (present(start)) then
        count = (/dims, 1/)
        ierr = nf90_put_var(ncid,varid,dat,start,count)
      else
        ierr = nf90_put_var(ncid,varid,dat)
      endif
      if (ierr/=nf90_noerr) then
        call handle_err(ierr,'ncwrite: '/ /vname)
      endif

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
      integer,dimension(4) :: count

      dims = shape(dat)
      ierr = nf90_inq_varid(ncid,vname,varid)
      if (ierr/=nf90_noerr) then
        call handle_err(ierr,'ncwrite :: '/ /vname)
      endif

      if (present(start)) then
        count = (/dims, 1/)
        ierr = nf90_put_var(ncid,varid,dat,start,count)
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
      integer,                      intent(in) :: ncid
      character(len=*),             intent(in) :: varname
      character(len=*),dimension(:),intent(in) :: dimname
      integer,dimension(:),optional,intent(in) :: dimsize
      integer,optional,             intent(in) :: vartype
      ! local
      integer :: i,ndim,varid,ierr,did,xtype
      integer,allocatable,dimension(:) :: dimid
      real(kind=8),parameter :: spv_set_dbl=1.D+33
      real(kind=4),parameter :: spv_set_flt=1.E+33

      if (present(vartype)) then                           ! handle optional arguement
        xtype = vartype
!        if(xtype/=nf90_float .and. xtype/=nf90_double) then
!          write(*,*)
!     &     'ERROR: nc_define_var - optional arguement must be nf90_float or nf90_double'
!           ierr=1; goto 2                                 ! can't use error stop as OMP_MASTER only operation?
!        endif
      else
        xtype = nf90_float                                 ! default output number type is float
      endif

      ndim = size(dimname)
      allocate(dimid(ndim))

      do i = 1,ndim                                        ! get dimension ids. Create if needed.
        ierr = nf90_inq_dimid(ncid,dimname(i),did)
        if (ierr.ne.nf90_noerr) then
          if (.not. present(dimsize)) then                 ! only an issue if dimension doesn't exist yet
            print *, 'ERROR: nccreate, no dimsize given for var=',varname
            error stop
          endif
          ierr=nf90_def_dim(ncid,dimname(i),dimsize(i),did)
          if (ierr/=nf90_noerr) call handle_err(ierr)
        endif
        dimid(i) = did
      enddo

      if (xtype==nf90_float .or. xtype==nf90_double) then  ! catch for non-numeric type
#ifdef HIS_DOUBLE
        xtype = nf90_double                                ! default output as double
#endif
      endif
      ierr=nf90_def_var(ncid,varname,xtype,dimid,varid,
     &                  deflate_level=deflate_level, shuffle=shuffle)
      if (ierr/=nf90_noerr) then
        print *,'ERROR creating: ',varname
        call handle_err(ierr)
      endif

      if     (xtype==nf90_float)  then                               ! fill value to blank out mask
        ierr = nf90_put_att(ncid,varid,'_FillValue',spv_set_flt)     ! in ncview
      elseif (xtype==nf90_double) then                               ! if we set the mask to 0 not e+33
        ierr = nf90_put_att(ncid,varid,'_FillValue',spv_set_dbl)     ! then we can remove this
      endif

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
