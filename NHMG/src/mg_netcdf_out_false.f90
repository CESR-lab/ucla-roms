module mg_netcdf_out

  !===================!
  !- USE association -!
  !===================!
  use mg_cst

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

  end subroutine sub_netcdf_write_fast_r3D


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

  end subroutine sub_netcdf_write_fast_r4D

end module mg_netcdf_out
