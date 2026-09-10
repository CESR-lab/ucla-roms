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
  logical           :: coarse_galerkin  = .false. !- coarse operators by aggregation (Galerkin) of the
                                                  !- fine one (always for a 2D level); .false.: geometric
  real(kind=rp)     :: robin_beta = 0._rp  !- implicit free surface: q + robin_beta dq/dz = f at the
                                          !- surface (0 = Dirichlet); set at run time by the caller
  real(kind=rp), allocatable :: sfcfac(:,:) !- fine-grid per-column multiplier of the surface factor
                                          !- (j,i): 1 Robin/Dirichlet, 0 prescribed surface flux
                                          !- (clamped open-boundary cells); allocated by the caller
  real(kind=rp), allocatable :: obcrob_u(:,:), obcrob_v(:,:)
                                          !- fine-grid radiating (implicit Flather) open faces, (j,i)
                                          !- on u faces 1..nx+1 / v faces 1..ny+1: n/(theta c dt) with
                                          !- n the outward normal sign (-1 west/south, +1 east/north),
                                          !- c = sqrt(g h); 0 elsewhere. The face flux correction is
                                          !- du = obcrob_u*Arx*p(adjacent cell): a diagonal-only,
                                          !- symmetric leak term. Allocated by the caller (implicit_fs).

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
       coarse_galerkin, &
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
          write(*,*)'  - coarse galerk : ', coarse_galerkin
          write(*,*)'  - E/W periodic  : ', east_west_perio
          write(*,*)'  - N/S periodic  : ', north_south_perio
          write(*,*)'  '
       endif
    endif

  end subroutine read_nhnamelist

  !----------------------------------------------------------------
  function sigtop(dzwtop) result(sig)
    ! surface factor of the top vertical coupling: 1 Dirichlet (q=0),
    ! 0 Neumann (dq/dz=0), dzw/(dzw+beta) for the Robin condition
    ! q + beta dq/dz = f (implicit free surface, beta = theta g dt^2)
    real(kind=rp), intent(in) :: dzwtop
    real(kind=rp) :: sig
    if (surface_neumann) then
       sig = 0._rp
    else
       sig = dzwtop / (dzwtop + robin_beta)
    endif
  end function sigtop

  !----------------------------------------------------------------
  function sfcf(j,i) result(f)
    ! fine-grid column multiplier of the surface factor (1 if unset)
    integer(kind=ip), intent(in) :: j,i
    real(kind=rp) :: f
    if (allocated(sfcfac)) then
       f = sfcfac(j,i)
    else
       f = 1._rp
    endif
  end function sfcf

end module mg_namelist
