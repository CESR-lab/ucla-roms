module mg_namelist

  use mg_cst
  use mg_tictoc

  implicit none

  integer(kind=ip) :: nsmall      =   8        ! smallest dimension ever for a subdomain; triggers a gather

  integer(kind=ip) :: ns_coarsest =  40        ! Number of relax sweeps for the coarsest grid level
  integer(kind=ip) :: ns_pre      =   3        ! Number of relax sweeps before coarsening  (going down)
  integer(kind=ip) :: ns_post     =   2        ! Number of relax sweeps after interpolation(going up)

  real(kind=rp)    :: solver_prec    = 1.d-6   !- solver precision
  integer(kind=ip) :: solver_maxiter = 10      !- maximum of solver iterations

  logical          :: autotune    =.false.     !- tuning test after a number of time steps (by default 100)
  integer(kind=ip) :: autotune_ts = 100        !- tuning test time step (if autotune=.true.)

  character(len=16) :: relax_method ='RB'      !- 'Gauss-Seidel', 'GS',
  !                                            !- 'Red-Black'   , 'RB',
  !                                            !- 'Four-Color'  , 'FC'

  logical           :: netcdf_output = .false. !- .false. or .true.

  integer(kind=ip)  :: output_freq = 100000000 ! Number of iterations between output of statistics

  logical           :: surface_neumann  = .true.

  logical           :: east_west_perio = .false.
  logical           :: north_south_perio = .false.

  namelist/nhparam/    &
       solver_prec   , &
       solver_maxiter, &
       nsmall        , &
       ns_coarsest   , &
       ns_pre        , &
       ns_post       , &
       autotune      , &
       autotune_ts   , &
       relax_method  , &
       netcdf_output , &
       output_freq   , &
       surface_neumann, &
       east_west_perio, &
       north_south_perio

contains

  !--------------------------------------------------------------------
  subroutine read_nhnamelist(filename, verbose, vbrank)

    character(len=*), optional, intent(in) :: filename
    logical         , optional, intent(in) :: verbose
    integer(kind=4) , optional, intent(in) :: vbrank

    character(len=64) :: fn_nml
    logical           :: vb
    integer(kind=ip)  :: lun_nml = 4
    integer(kind=4)   :: rank

    logical :: exist=.false.

    !- Namelist file name, by default 'nhmg_namelist'
    if (present(filename)) then
       fn_nml = filename
    else
       fn_nml = 'nhmg_namelist'
    endif

    !- Check if a namelist file exist
    inquire(file=fn_nml, exist=exist)

    !- Read namelist file if it is present, else use default values
    if (exist) then

       if (vbrank == 0) then
          write(*,*)"  Opening and Reading namelist file:", TRIM(fn_nml)
       endif

       open(unit=lun_nml, File=fn_nml, ACTION='READ')

       rewind(unit=lun_nml)
       read(unit=lun_nml, nml=nhparam)

    endif

    !- Print parameters or not !
    if (present(verbose)) then
       vb = verbose
    else
       vb = .true.
    endif

    if (vb) then

       if (present(vbrank)) then
          rank = vbrank
       else
          rank = 0
       endif

       if (rank == 0) then
          write(*,*)'  Non hydrostatic parameters:'
          write(*,*)'  - solver_prec   : ', solver_prec
          write(*,*)'  - solver_maxiter: ', solver_maxiter
          write(*,*)'  - nsmall        : ', nsmall
          write(*,*)'  - ns_coarsest   : ', ns_coarsest
          write(*,*)'  - ns_pre        : ', ns_pre
          write(*,*)'  - ns_post       : ', ns_post
          write(*,*)'  - autotune      : ', autotune
          write(*,*)'  - autotune_ts   : ', autotune_ts
          write(*,*)'  - relax_method  : ', trim(relax_method)
          write(*,*)'  - netcdf_output : ', netcdf_output
          write(*,*)'  - output freq   : ', output_freq
          write(*,*)'  - surf neumann  : ', surface_neumann
          write(*,*)'  - E/W periodic  : ', east_west_perio
          write(*,*)'  - N/S periodic  : ', north_south_perio
          write(*,*)'  '
       endif
    endif

  end subroutine read_nhnamelist

end module mg_namelist
