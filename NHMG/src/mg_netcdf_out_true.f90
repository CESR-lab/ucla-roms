module mg_netcdf_out

  !===================!
  !- USE association -!
  !===================!
  use mg_cst
  use netcdf

  !===============================!
  !- DECLARATIONS of GLOBAL DATA -!
  !===============================!
  implicit none

  !====================================!
  !- DECLARATIONS of INTERFACE BLOCKS -!
  !====================================!
  interface write_netcdf
     module procedure                &
          sub_netcdf_write_fast_r2D, &
          sub_netcdf_write_fast_r3D, &
!!$       sub_netcdf_write_fast_r3D_p, &
          sub_netcdf_write_fast_r4D
  end interface write_netcdf

contains
  !********************************************************************!
  !
  ! Create a NetCDF File and store in it a 2D array of reals
  !
  !********************************************************************!
  subroutine sub_netcdf_write_fast_r2D( &
       array_r2D                      , & ! 2D array of reals.
       vname                          , & ! NetCDF variable name (opt).
       netcdf_file_name               , & ! NetCDF file name (opt).
       rank                           , & ! Add my rank in the NetCDF file name (opt)
       iter                           )   ! Add a number in the NetCDF 
    !                                     ! file name (opt).

    !=============================!
    !- DECLARATIONS of ARGUMENTS -!
    !=============================!
    real(kind = rp), dimension(:,:), intent(in) :: array_r2D
    character(len = *)                , optional, intent(in) :: vname
    character(len = *)                , optional, intent(in) :: netcdf_file_name
    integer(kind = ip)             , optional, intent(in) :: rank
    integer(kind = ip)             , optional, intent(in) :: iter

    !===================================!
    !- DECLARATIONS of LOCAL VARIABLES -!
    !===================================!
    integer(kind = ip), parameter :: nb_digit = 5

    character(len=48)     :: cmft
    integer(kind = ip) :: size_string

    character(len = 64)   :: final_netcdf_file_name
    character(len = 32)   :: final_vname
    integer(kind = ip) :: fiter
    integer(kind = ip) :: is_err
    integer(kind = ip) :: add_rank

    integer(kind = ip) :: dim1
    integer(kind = ip) :: dim2

    integer(kind = ip) :: nc_id
    integer(kind = ip) :: var_id
    integer(kind = ip) :: dim1_id
    integer(kind = ip) :: dim2_id

    !======================================!
    !- SECTION 1: test optional arguments -!
    !======================================!

    if (present(netcdf_file_name)) then
       final_netcdf_file_name = trim(netcdf_file_name)
    else
       final_netcdf_file_name = 'netcdf_file.nc'
    end if

    !- Test if the end of the netcdf file is ".nc".   -!
    !- If not stop the program with an error message. -!
    size_string = len_trim(final_netcdf_file_name)
    if (final_netcdf_file_name(size_string-2:size_string) /= '.nc') then
       write(*,*)'Error: mg_netcdf_out.f90'
       write(*,*)'NetCDF file name has to finish with .nc'
       write(*,*)'final_netcdf_file_name(size_string-2:size_string) =', &
            final_netcdf_file_name(size_string-2:size_string)
       stop -1
    end if

    if (present(vname)) then
       final_vname = trim(vname)
    else
       final_vname = 'array_r3D'
    end if

    !- Add var name in the NetCDF file name -!
    !- just before ".nc"
    if (final_vname /= 'array_r3D') then
       write(cmft, 1001) size_string-3, len_trim(final_vname)
       write(final_netcdf_file_name,FMT=trim(cmft)) &
            final_netcdf_file_name(1:size_string-3),'_',final_vname,'.nc'
    end if

1001 format('(A',I2,',A1,A',I2,',A3)')

    if ( present(rank)) then
       size_string = len_trim(final_netcdf_file_name)
       add_rank = rank
    else
       add_rank = -1
    end if

    !- Add rank in the NetCDF file name -!
    !- just before ".nc"
    if (add_rank /= -1) then
       write(cmft, 1002) size_string-3,3,3
       write(final_netcdf_file_name,FMT=trim(cmft)) &
            final_netcdf_file_name(1:size_string-3),'_',add_rank,'.nc'
    end if

1002 format('(A',I2,',A1,I',I2,'.',I2,',A3)')

    if ( present(iter)) then
       size_string = len_trim(final_netcdf_file_name)
       fiter = iter
    else
       fiter = -1
    end if

    !- Add a digit number in the NetCDF file name -!
    !- just before ".nc"
    if (fiter /= -1) then
       write(cmft, 1003) size_string-3, nb_digit, nb_digit
       write(final_netcdf_file_name, FMT=trim(cmft))   &
            final_netcdf_file_name(1:size_string-3), '_', &
            fiter,'.nc'
    end if

1003 format('(A',I2,',A1,I',I2,'.',I2,',A3)')

    !======================!
    !- CREATE NetCDF FILE -!
    !======================!
    is_err = nf90_create                     ( &
         path  = trim(final_netcdf_file_name), & ! NetCDF filename created.
         cmode = NF90_64BIT_OFFSET           , & ! 64_bits record, Large File
         ncid  = nc_id                       )

    if (is_err /= nf90_noerr) then
       stop
    end if

    !------------------!
    !- Set Dimensions -!
    !------------------!
    dim1 = size(array_r2D,1)
    dim2 = size(array_r2D,2)

    is_err=NF90_def_dim ( &
         ncid  = nc_id  , &
         name  = 'dim1' , &
         len   = dim1   , &
         dimid = dim1_id)

    if (is_err /= nf90_noerr) then
       stop
    end if

    is_err=NF90_def_dim ( &
         ncid  = nc_id  , &
         name  = 'dim2' , &
         len   = dim2   , &
         dimid = dim2_id)

    if (is_err /= nf90_noerr) then
       stop
    end if

    !-----------------------!
    !- Variable Definition -!
    !-----------------------!
    is_err=NF90_def_var              ( &
         ncid   = nc_id              , &
         name   = trim(final_vname)  , &
         xtype  = NF90_DOUBLE        , &
         dimids = (/dim1_id,dim2_id/), &
         varid  = var_id             )

    if (is_err /= nf90_noerr) then
       stop
    end if

    !------------------------------!
    !- Close the header defintion -!
    !------------------------------!
    is_err =  NF90_enddef(ncid = nc_id)

    if (is_err /= nf90_noerr) then
       stop
    end if

    !-------------------------------!
    !- Put data in the NetCDF file -!
    !-------------------------------!
    is_err = NF90_put_var   ( &
         ncid   = nc_id     , &
         varid  = var_id    , &
         values = array_r2D )

    if (is_err /= nf90_noerr) then
       stop
    end if

    !-------------------------!
    !- Close the netcdf file -!
    !-------------------------!
    is_err = nf90_close(ncid = nc_id)

    if (is_err /= nf90_noerr) then
       stop
    end if

  end subroutine sub_netcdf_write_fast_r2D

  !********************************************************************!
  !
  ! Create a NetCDF File and store in it a 3D array of reals
  !
  !********************************************************************!
  subroutine sub_netcdf_write_fast_r3D( &
       array_r3D                      , & ! 3D array of reals.
       vname                          , & ! NetCDF variable name (opt).
       netcdf_file_name               , & ! NetCDF file name (opt).
       rank                           , & ! Add my rank in the NetCDF file name (opt)
       iter                           )   ! Add a number in the NetCDF 
    !                                     ! file name (opt).

    !=============================!
    !- DECLARATIONS of ARGUMENTS -!
    !=============================!
    real(kind = rp), dimension(:,:,:)          , intent(in) :: array_r3D
    character(len = *)               , optional, intent(in) :: vname
    character(len = *)               , optional, intent(in) :: netcdf_file_name
    integer(kind = ip)               , optional, intent(in) :: rank
    integer(kind = ip)               , optional, intent(in) :: iter

    !===================================!
    !- DECLARATIONS of LOCAL VARIABLES -!
    !===================================!
    integer(kind = ip), parameter :: nb_digit = 5

    character(len=48)     :: cmft
    integer(kind = ip) :: size_string

    character(len = 32)   :: final_netcdf_file_name
    character(len = 32)   :: final_vname
    integer(kind = ip) :: fiter
    integer(kind = ip) :: is_err
    integer(kind = ip) :: add_rank

    integer(kind = ip) :: dim1
    integer(kind = ip) :: dim2
    integer(kind = ip) :: dim3

    integer(kind = ip) :: nc_id
    integer(kind = ip) :: var_id
    integer(kind = ip) :: dim1_id
    integer(kind = ip) :: dim2_id
    integer(kind = ip) :: dim3_id

    !======================================!
    !- SECTION 1: test optional arguments -!
    !======================================!

    if (present(netcdf_file_name)) then
       final_netcdf_file_name = trim(netcdf_file_name)
    else
       final_netcdf_file_name = 'netcdf_file.nc'
    end if

    !- Test if the end of the netcdf file is ".nc".   -!
    !- If not stop the program with an error message. -!
    size_string = len_trim(final_netcdf_file_name)
    if (final_netcdf_file_name(size_string-2:size_string) /= '.nc') then
       write(*,*)'Error: mg_netcdf_out.f90'
       write(*,*)'NetCDF file name has to finish with .nc'
       write(*,*)'final_netcdf_file_name(size_string-2:size_string) =', &
            final_netcdf_file_name(size_string-2:size_string)
       stop -1
    end if

    if (present(vname)) then
       final_vname = trim(vname)
    else
       final_vname = 'array_r3D'
    end if

    !- Add var name in the NetCDF file name -!
    !- just before ".nc"
    if (final_vname /= 'array_r3D') then
       write(cmft, 1001) size_string-3, len_trim(final_vname)
       write(final_netcdf_file_name,FMT=trim(cmft)) &
            final_netcdf_file_name(1:size_string-3),'_',final_vname,'.nc'
    end if

1001 format('(A',I2,',A1,A',I2,',A3)')

    if ( present(rank)) then
       size_string = len_trim(final_netcdf_file_name)
       add_rank = rank
    else
       add_rank = -1
    end if

    !- Add rank in the NetCDF file name -!
    !- just before ".nc"
    if (add_rank /= -1) then
       write(cmft, 1002) size_string-3,3,3
       write(final_netcdf_file_name,FMT=trim(cmft)) &
            final_netcdf_file_name(1:size_string-3),'_',add_rank,'.nc'
    end if

1002 format('(A',I2,',A1,I',I2,'.',I2,',A3)')

    if ( present(iter)) then
       size_string = len_trim(final_netcdf_file_name)
       fiter = iter
    else
       fiter = -1
    end if

    !- Add a digit number in the NetCDF file name -!
    !- just before ".nc"
    if (fiter /= -1) then
       write(cmft, 1003) size_string-3, nb_digit, nb_digit
       write(final_netcdf_file_name, FMT=trim(cmft))   &
            final_netcdf_file_name(1:size_string-3), '_', &
            fiter,'.nc'
    end if

1003 format('(A',I2,',A1,I',I2,'.',I2,',A3)')

    !======================!
    !- CREATE NetCDF FILE -!
    !======================!
    is_err = nf90_create                     ( &
         path  = trim(final_netcdf_file_name), & ! NetCDF filename created.
         cmode = NF90_64BIT_OFFSET           , & ! 64_bits record, Large File
         ncid  = nc_id                       )

    if (is_err /= nf90_noerr) then
       stop
    end if

    !------------------!
    !- Set Dimensions -!
    !------------------!
    dim1 = size(array_r3D,1)
    dim2 = size(array_r3D,2)
    dim3 = size(array_r3D,3)

    is_err=NF90_def_dim ( &
         ncid  = nc_id  , &
         name  = 'dim1' , &
         len   = dim1   , &
         dimid = dim1_id)

    if (is_err /= nf90_noerr) then
       stop
    end if

    is_err=NF90_def_dim ( &
         ncid  = nc_id  , &
         name  = 'dim2' , &
         len   = dim2   , &
         dimid = dim2_id)

    if (is_err /= nf90_noerr) then
       stop
    end if

    is_err=NF90_def_dim ( &
         ncid  = nc_id  , &
         name  = 'dim3' , &
         len   = dim3   , &
         dimid = dim3_id)

    if (is_err /= nf90_noerr) then
       stop
    end if

    !-----------------------!
    !- Variable Definition -!
    !-----------------------!
    is_err=NF90_def_var                      ( &
         ncid   = nc_id                      , &
         name   = trim(final_vname)          , &
         xtype  = NF90_DOUBLE                , &
         dimids = (/dim1_id,dim2_id,dim3_id/), &
         varid  = var_id                     )

    if (is_err /= nf90_noerr) then
       stop
    end if

    !------------------------------!
    !- Close the header defintion -!
    !------------------------------!
    is_err =  NF90_enddef(ncid = nc_id)

    if (is_err /= nf90_noerr) then
       stop
    end if

    !-------------------------------!
    !- Put data in the NetCDF file -!
    !-------------------------------!
    is_err = NF90_put_var   ( &
         ncid   = nc_id     , &
         varid  = var_id    , &
         values = array_r3D )

    if (is_err /= nf90_noerr) then
       stop
    end if

    !-------------------------!
    !- Close the netcdf file -!
    !-------------------------!
    is_err = nf90_close(ncid = nc_id)

    if (is_err /= nf90_noerr) then
       stop
    end if

  end subroutine sub_netcdf_write_fast_r3D

!!$  !********************************************************************!
!!$  !
!!$  ! Create a NetCDF File and store in it a 3D array of reals
!!$  !
!!$  !********************************************************************!
!!$  subroutine sub_netcdf_write_fast_r3D_p( &
!!$       array_r3D                      , & ! 3D array of reals.
!!$       vname                          , & ! NetCDF variable name (opt).
!!$       netcdf_file_name               , & ! NetCDF file name (opt).
!!$       rank                           , & ! Add my rank in the NetCDF file name (opt)
!!$       iter                           )   ! Add a number in the NetCDF 
!!$    !                                     ! file name (opt).
!!$
!!$    !=============================!
!!$    !- DECLARATIONS of ARGUMENTS -!
!!$    !=============================!
!!$    real(kind = rp), dimension(:,:,:), pointer , intent(in) :: array_r3D
!!$    character(len = *)                  , optional, intent(in) :: vname
!!$    character(len = *)                  , optional, intent(in) :: netcdf_file_name
!!$    integer(kind = ip)               , optional, intent(in) :: rank
!!$    integer(kind = ip)               , optional, intent(in) :: iter
!!$
!!$    !===================================!
!!$    !- DECLARATIONS of LOCAL VARIABLES -!
!!$    !===================================!
!!$    integer(kind = ip), parameter :: nb_digit = 5
!!$
!!$    character(len=48)     :: cmft
!!$    integer(kind = ip) :: size_string
!!$
!!$    character(len = 32)   :: final_netcdf_file_name
!!$    character(len = 32)   :: final_vname
!!$    integer(kind = ip) :: fiter
!!$    integer(kind = ip) :: is_err
!!$    integer(kind = ip) :: add_rank
!!$
!!$    integer(kind = ip) :: dim1
!!$    integer(kind = ip) :: dim2
!!$    integer(kind = ip) :: dim3
!!$
!!$    integer(kind = ip) :: nc_id
!!$    integer(kind = ip) :: var_id
!!$    integer(kind = ip) :: dim1_id
!!$    integer(kind = ip) :: dim2_id
!!$    integer(kind = ip) :: dim3_id
!!$
!!$    !======================================!
!!$    !- SECTION 1: test optional arguments -!
!!$    !======================================!
!!$
!!$    if (present(netcdf_file_name)) then
!!$       final_netcdf_file_name = trim(netcdf_file_name)
!!$    else
!!$       final_netcdf_file_name = 'netcdf_file.nc'
!!$    end if
!!$
!!$    !- Test if the end of the netcdf file is ".nc".   -!
!!$    !- If not stop the program with an error message. -!
!!$    size_string = len_trim(final_netcdf_file_name)
!!$    if (final_netcdf_file_name(size_string-2:size_string) /= '.nc') then
!!$       write(*,*)'Error: mg_netcdf_out.f90'
!!$       write(*,*)'NetCDF file name has to finish with .nc'
!!$       write(*,*)'final_netcdf_file_name(size_string-2:size_string) =', &
!!$            final_netcdf_file_name(size_string-2:size_string)
!!$       stop -1
!!$    end if
!!$
!!$    if (present(vname)) then
!!$       final_vname = trim(vname)
!!$    else
!!$       final_vname = 'array_r3D'
!!$    end if
!!$
!!$    !- Add var name in the NetCDF file name -!
!!$    !- just before ".nc"
!!$    if (final_vname /= 'array_r3D') then
!!$       write(cmft, 1001) size_string-3, len_trim(final_vname)
!!$       write(final_netcdf_file_name,FMT=trim(cmft)) &
!!$            final_netcdf_file_name(1:size_string-3),'_',final_vname,'.nc'
!!$    end if
!!$
!!$1001 format('(A',I2,',A1,A',I2,',A3)')
!!$
!!$    if ( present(rank)) then
!!$       size_string = len_trim(final_netcdf_file_name)
!!$       add_rank = rank
!!$    else
!!$       add_rank = -1
!!$    end if
!!$
!!$    !- Add rank in the NetCDF file name -!
!!$    !- just before ".nc"
!!$    if (add_rank /= -1) then
!!$       write(cmft, 1002) size_string-3,3,3
!!$       write(final_netcdf_file_name,FMT=trim(cmft)) &
!!$            final_netcdf_file_name(1:size_string-3),'_',add_rank,'.nc'
!!$    end if
!!$
!!$1002 format('(A',I2,',A1,I',I2,'.',I2,',A3)')
!!$
!!$    if ( present(iter)) then
!!$       size_string = len_trim(final_netcdf_file_name)
!!$       fiter = iter
!!$    else
!!$       fiter = -1
!!$    end if
!!$
!!$    !- Add a digit number in the NetCDF file name -!
!!$    !- just before ".nc"
!!$    if (fiter /= -1) then
!!$       write(cmft, 1003) size_string-3, nb_digit, nb_digit
!!$       write(final_netcdf_file_name, FMT=trim(cmft))   &
!!$            final_netcdf_file_name(1:size_string-3), '_', &
!!$            fiter,'.nc'
!!$    end if
!!$
!!$1003 format('(A',I2,',A1,I',I2,'.',I2,',A3)')
!!$
!!$    !======================!
!!$    !- CREATE NetCDF FILE -!
!!$    !======================!
!!$    is_err = nf90_create                     ( &
!!$         path  = trim(final_netcdf_file_name), & ! NetCDF filename created.
!!$         cmode = NF90_64BIT_OFFSET           , & ! 64_bits record, Large File
!!$         ncid  = nc_id                       )
!!$
!!$    if (is_err /= nf90_noerr) then
!!$       stop
!!$    end if
!!$
!!$    !------------------!
!!$    !- Set Dimensions -!
!!$    !------------------!
!!$    dim1 = size(array_r3D,1)
!!$    dim2 = size(array_r3D,2)
!!$    dim3 = size(array_r3D,3)
!!$
!!$    is_err=NF90_def_dim ( &
!!$         ncid  = nc_id  , &
!!$         name  = 'dim1' , &
!!$         len   = dim1   , &
!!$         dimid = dim1_id)
!!$
!!$    if (is_err /= nf90_noerr) then
!!$       stop
!!$    end if
!!$
!!$    is_err=NF90_def_dim ( &
!!$         ncid  = nc_id  , &
!!$         name  = 'dim2' , &
!!$         len   = dim2   , &
!!$         dimid = dim2_id)
!!$
!!$    if (is_err /= nf90_noerr) then
!!$       stop
!!$    end if
!!$
!!$    is_err=NF90_def_dim ( &
!!$         ncid  = nc_id  , &
!!$         name  = 'dim3' , &
!!$         len   = dim3   , &
!!$         dimid = dim3_id)
!!$
!!$    if (is_err /= nf90_noerr) then
!!$       stop
!!$    end if
!!$
!!$    !-----------------------!
!!$    !- Variable Definition -!
!!$    !-----------------------!
!!$    is_err=NF90_def_var                      ( &
!!$         ncid   = nc_id                      , &
!!$         name   = trim(final_vname)          , &
!!$         xtype  = NF90_DOUBLE                , &
!!$         dimids = (/dim1_id,dim2_id,dim3_id/), &
!!$         varid  = var_id                     )
!!$
!!$    if (is_err /= nf90_noerr) then
!!$       stop
!!$    end if
!!$
!!$    !------------------------------!
!!$    !- Close the header defintion -!
!!$    !------------------------------!
!!$    is_err =  NF90_enddef(ncid = nc_id)
!!$
!!$    if (is_err /= nf90_noerr) then
!!$       stop
!!$    end if
!!$
!!$    !-------------------------------!
!!$    !- Put data in the NetCDF file -!
!!$    !-------------------------------!
!!$    is_err = NF90_put_var   ( &
!!$         ncid   = nc_id     , &
!!$         varid  = var_id    , &
!!$         values = array_r3D )
!!$
!!$    if (is_err /= nf90_noerr) then
!!$       stop
!!$    end if
!!$
!!$    !-------------------------!
!!$    !- Close the netcdf file -!
!!$    !-------------------------!
!!$    is_err = nf90_close(ncid = nc_id)
!!$
!!$    if (is_err /= nf90_noerr) then
!!$       stop
!!$    end if
!!$
!!$  end subroutine sub_netcdf_write_fast_r3D_p
  !********************************************************************!
  !
  ! Create a NetCDF File and store in it a 4D array of reals
  !
  !********************************************************************!
  subroutine sub_netcdf_write_fast_r4D( &
       array_r4D                      , & ! 4D array of reals.
       vname                          , & ! NetCDF variable name (opt).
       netcdf_file_name               , & ! NetCDF file name (opt).
       rank                           , & ! Add my rank in the NetCDF file name (opt)
       iter                           )   ! Add a number in the NetCDF 
    !                                     ! file name (opt).

    !=============================!
    !- DECLARATIONS of ARGUMENTS -!
    !=============================!
    real(kind = rp), dimension(:,:,:,:), pointer , intent(in) :: array_r4D
    character(len = *)                    , optional, intent(in) :: vname
    character(len = *)                    , optional, intent(in) :: netcdf_file_name
    integer(kind = ip)                 , optional, intent(in) :: rank
    integer(kind = ip)                 , optional, intent(in) :: iter

    !===================================!
    !- DECLARATIONS of LOCAL VARIABLES -!
    !===================================!
    integer(kind = ip), parameter :: nb_digit = 5

    character(len=48)     :: cmft
    integer(kind = ip) :: size_string

    character(len = 32)   :: final_netcdf_file_name
    character(len = 32)   :: final_vname
    integer(kind = ip) :: fiter
    integer(kind = ip) :: is_err
    integer(kind = ip) :: add_rank

    integer(kind = ip) :: dim1
    integer(kind = ip) :: dim2
    integer(kind = ip) :: dim3
    integer(kind = ip) :: dim4

    integer(kind = ip) :: nc_id
    integer(kind = ip) :: var_id
    integer(kind = ip) :: dim1_id
    integer(kind = ip) :: dim2_id
    integer(kind = ip) :: dim3_id
    integer(kind = ip) :: dim4_id

    !======================================!
    !- SECTION 1: test optional arguments -!
    !======================================!

    if (present(netcdf_file_name)) then
       final_netcdf_file_name = trim(netcdf_file_name)
    else
       final_netcdf_file_name = 'netcdf_file.nc'
    end if

    !- Test if the end of the netcdf file is ".nc".   -!
    !- If not stop the program with an error message. -!
    size_string = len_trim(final_netcdf_file_name)
    if (final_netcdf_file_name(size_string-2:size_string) /= '.nc') then
       write(*,*)'Error: mg_netcdf_out.f90'
       write(*,*)'NetCDF file name has to finish with .nc'
       write(*,*)'final_netcdf_file_name(size_string-2:size_string) =', &
            final_netcdf_file_name(size_string-2:size_string)
       stop -1
    end if

    if (present(vname)) then
       final_vname = trim(vname)
    else
       final_vname = 'array_r4D'
    end if

    !- Add var name in the NetCDF file name -!
    !- just before ".nc"
    if (final_vname /= 'array_r4D') then
       write(cmft, 1001) size_string-3, len_trim(final_vname)
       write(final_netcdf_file_name,FMT=trim(cmft)) &
            final_netcdf_file_name(1:size_string-3),'_',final_vname,'.nc'
    end if

1001 format('(A',I2,',A1,A',I2,',A3)')

    if ( present(rank)) then
       size_string = len_trim(final_netcdf_file_name)
       add_rank = rank
    else
       add_rank = -1
    end if

    !- Add rank in the NetCDF file name -!
    !- just before ".nc"
    if (add_rank /= -1) then
       write(cmft, 1002) size_string-3,3,3
       write(final_netcdf_file_name,FMT=trim(cmft)) &
            final_netcdf_file_name(1:size_string-3),'_',add_rank,'.nc'
    end if

1002 format('(A',I2,',A1,I',I2,'.',I2,',A3)')

    if ( present(iter)) then
       size_string = len_trim(final_netcdf_file_name)
       fiter = iter
    else
       fiter = -1
    end if

    !- Add a digit number in the NetCDF file name -!
    !- just before ".nc"
    if (fiter /= -1) then
       write(cmft, 1003) size_string-3, nb_digit, nb_digit
       write(final_netcdf_file_name, FMT=trim(cmft))   &
            final_netcdf_file_name(1:size_string-3), '_', &
            fiter,'.nc'
    end if

1003 format('(A',I2,',A1,I',I2,'.',I2,',A3)')

    !======================!
    !- CREATE NetCDF FILE -!
    !======================!
    is_err = nf90_create                     ( &
         path  = trim(final_netcdf_file_name), & ! NetCDF filename created.
         cmode = NF90_64BIT_OFFSET           , & ! 64_bits record, Large File
         ncid  = nc_id                       )

    if (is_err /= nf90_noerr) then
       stop
    end if

    !------------------!
    !- Set Dimensions -!
    !------------------!
    dim1 = size(array_r4D,1)
    dim2 = size(array_r4D,2)
    dim3 = size(array_r4D,3)
    dim4 = size(array_r4D,4)

    is_err=NF90_def_dim ( &
         ncid  = nc_id  , &
         name  = 'dim1' , &
         len   = dim1   , &
         dimid = dim1_id)

    if (is_err /= nf90_noerr) then
       stop
    end if

    is_err=NF90_def_dim ( &
         ncid  = nc_id  , &
         name  = 'dim2' , &
         len   = dim2   , &
         dimid = dim2_id)

    if (is_err /= nf90_noerr) then
       stop
    end if

    is_err=NF90_def_dim ( &
         ncid  = nc_id  , &
         name  = 'dim3' , &
         len   = dim3   , &
         dimid = dim3_id)

    if (is_err /= nf90_noerr) then
       stop
    end if

    is_err=NF90_def_dim ( &
         ncid  = nc_id  , &
         name  = 'dim4' , &
         len   = dim4   , &
         dimid = dim4_id)

    if (is_err /= nf90_noerr) then
       stop
    end if

    !-----------------------!
    !- Variable Definition -!
    !-----------------------!
    is_err=NF90_def_var                              ( &
         ncid   = nc_id                              , &
         name   = trim(final_vname)                  , &
         xtype  = NF90_DOUBLE                        , &
         dimids = (/dim1_id,dim2_id,dim3_id,dim4_id/), &
         varid  = var_id                             )

    if (is_err /= nf90_noerr) then
       stop
    end if

    !------------------------------!
    !- Close the header defintion -!
    !------------------------------!
    is_err =  NF90_enddef(ncid = nc_id)

    if (is_err /= nf90_noerr) then
       stop
    end if

    !-------------------------------!
    !- Put data in the NetCDF file -!
    !-------------------------------!
    is_err = NF90_put_var   ( &
         ncid   = nc_id     , &
         varid  = var_id    , &
         values = array_r4D )

    if (is_err /= nf90_noerr) then
       stop
    end if

    !-------------------------!
    !- Close the netcdf file -!
    !-------------------------!
    is_err = nf90_close(ncid = nc_id)

    if (is_err /= nf90_noerr) then
       stop
    end if

  end subroutine sub_netcdf_write_fast_r4D

end module mg_netcdf_out
