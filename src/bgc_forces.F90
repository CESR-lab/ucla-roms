module bgc_forces
  ! formerly bgc_forces.h

#include "cppdefs.opt"
#if defined(BIOLOGY_BEC2) || defined(MARBL)

  use param, only: lm, mm
  implicit none

#ifdef PCO2AIR_FORCING
! xCO2air concentration
! ------- -------------
!     xCO2air: xCO2air concentraion [ppm]
  real(kind=8),allocatable,dimension(:,:) :: xco2air
#ifdef MARBL
  real(kind=8),allocatable,dimension(:,:) :: xco2air_alt
#endif
!SDISTRIBUTE_RESHAPE  xCO2air(BLOCK_PATTERN,*) BLOCK_CLAUSE
# if defined PCO2AIR_DATA || defined ALL_DATA
# ifndef SET_SMTH
#  undef PCO2AIR_DATA
# endif
!SDISTRIBUTE_RESHAPE  xco2airg(BLOCK_PATTERN,*) BLOCK_CLAUSE
  real(kind=8) :: xco2air_cycle, xco2air_time(2)
  integer(kind=4) :: xco2air_ncycle,  xco2air_rec, itxco2air, ntxco2air,&
  &xco2air_file_id, xco2air_id,  xco2air_tid
# endif
#endif /* PCO2AIR_FORCING */

! Atmospheric xCO2 [ppm] applied when PCO2AIR_FORCING is off (not read from a forcing
! file). Set via the BGC_SETTINGS namelist; read in bgc_shared_vars::read_nml_bgc.
  real(kind=8) :: xco2air_default = 284.7   ! pre-industrial [micromol/mol]

#if defined DAILYPAR_PHOTOINHIBITION || defined DAILYPAR_BEC
! daily avg swrad
! ------- -------------
  real(kind=8),allocatable,dimension(:,:) :: swrad_avg
!SDISTRIBUTE_RESHAPE  swrad_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
# if defined SWRAD_AVG_DATA || defined ALL_DATA
# ifndef SET_SMTH
#  undef SWRAD_AVG_DATA
# endif
!SDISTRIBUTE_RESHAPE  swrad_avgg(BLOCK_PATTERN,*) BLOCK_CLAUSE
  real(kind=8) :: swrad_avg_cycle, swrad_avg_time(2)
  integer(kind=4) :: swrad_avg_ncycle,  swrad_avg_rec, itswrad_avg, ntswrad_avg,&
  &swrad_avg_file_id, swrad_avg_id,  swrad_avg_tid
# endif
#endif /* DAILYPAR_PHOTOINHIBITION || defined DAILYPAR_BEC */


#ifdef NHY_FORCING
! NHY flux
! --- ----
  real(kind=8),allocatable,dimension(:,:) :: nhy
!SDISTRIBUTE_RESHAPE  nhy(BLOCK_PATTERN,*) BLOCK_CLAUSE
# if defined NHY_DATA || defined ALL_DATA
# ifndef SET_SMTH
#  undef NHY_DATA
# endif
  real(kind=8),allocatable,dimension(:,:,:) :: nhyg
!SDISTRIBUTE_RESHAPE  nhyg(BLOCK_PATTERN,*) BLOCK_CLAUSE
  real(kind=8)    :: nhy_cycle, nhy_time(2)
  integer(kind=4) :: nhy_ncycle,  nhy_rec, itnhy, ntnhy,&
  &nhy_file_id, nhy_id, nhy_tid
# endif
#endif /* NHY_FORCING */

#ifdef NOX_FORCING
! NOX flux
! --- ----
  real(kind=8),allocatable,dimension(:,:) :: nox
!SDISTRIBUTE_RESHAPE  nox(BLOCK_PATTERN,*) BLOCK_CLAUSE
# if defined NOX_DATA || defined ALL_DATA
# ifndef SET_SMTH
#  undef NOX_DATA
# endif
  real(kind=8),allocatable,dimension(:,:,:) :: noxg
!SDISTRIBUTE_RESHAPE  noxg(BLOCK_PATTERN,*) BLOCK_CLAUSE
  real(kind=8)    :: nox_cycle, nox_time(2)
  integer(kind=4) :: nox_ncycle,  nox_rec, itnox, ntnox,&
  &nox_file_id, nox_id, nox_tid
# endif
#endif /* NOX_FORCING */

#if defined BIOLOGY_BEC || defined BIOLOGY_BEC2 || defined(MARBL)
! dust flux
! ---- ----
!      dust: dust flux [kg m-2 s-1]
  real(kind=8),allocatable,dimension(:,:) :: dust
!SDISTRIBUTE_RESHAPE  dust(BLOCK_PATTERN) BLOCK_CLAUSE
# if defined DUST_DATA || defined ALL_DATA
#  ifndef SET_SMTH
#   undef DUST_DATA
#  endif
!SDISTRIBUTE_RESHAPE  dustg(BLOCK_PATTERN,*) BLOCK_CLAUSE
  real(kind=8) :: dust_cycle, dust_time(2)
  integer(kind=4) :: dust_ncycle,  dust_rec, itdust, ntdust,&
  &dust_file_id, dust_id,  dust_tid
# endif /* defined DUST_DATA || defined ALL_DATA */

! iron flux
! ---- ----
!     iron: iron flux [nmol cm-2 s-1]
  real(kind=8),allocatable,dimension(:,:) :: iron
!SDISTRIBUTE_RESHAPE  iron(BLOCK_PATTERN) BLOCK_CLAUSE
# if defined IRON_DATA || defined ALL_DATA
#  ifndef SET_SMTH
#   undef IRON_DATA
#  endif
!SDISTRIBUTE_RESHAPE  irong(BLOCK_PATTERN,*) BLOCK_CLAUSE
  real(kind=8) :: iron_cycle, iron_time(2)
  integer(kind=4) :: iron_ncycle,  iron_rec, itiron, ntiron,&
  &iron_file_id, iron_id,  iron_tid
# endif /* defined IRON_DATA || defined ALL_DATA */

#endif /* BIOLOG_BEC || MARBL */

contains

!----------------------------------------------------------------------
  subroutine init_arrays_bgc_forces  ![
    implicit none

#ifdef PCO2AIR_FORCING
    allocate( xco2air(GLOBAL_2D_ARRAY) )
#ifdef MARBL
    allocate( xco2air_alt(GLOBAL_2D_ARRAY) )
#endif
#endif /* PCO2AIR_FORCING */

#if defined DAILYPAR_PHOTOINHIBITION || defined DAILYPAR_BEC
    allocate( swrad_avg(GLOBAL_2D_ARRAY) )  ! daily average short-wave rad
#endif

#ifdef NHY_FORCING
    allocate( nhy(GLOBAL_2D_ARRAY) )
# if defined NHY_DATA || defined ALL_DATA
# ifndef SET_SMTH
#  undef NHY_DATA
# endif
    allocate( nhyg(GLOBAL_2D_ARRAY,2) )
# endif
#endif /* NHY_FORCING */

#ifdef NOX_FORCING
    allocate( nox(GLOBAL_2D_ARRAY) )
# if defined NOX_DATA || defined ALL_DATA
# ifndef SET_SMTH
#  undef NOX_DATA
# endif
    allocate( noxg(GLOBAL_2D_ARRAY,2) )
# endif
#endif /* NOX_FORCING */

#if defined BIOLOGY_BEC || defined BIOLOGY_BEC2 || defined MARBL
    allocate( dust(GLOBAL_2D_ARRAY) )
# if defined DUST_DATA || defined ALL_DATA
#  ifndef SET_SMTH
#   undef DUST_DATA
#  endif
# endif /* defined DUST_DATA || defined ALL_DATA */

    allocate( iron(GLOBAL_2D_ARRAY) )
# if defined IRON_DATA || defined ALL_DATA
#  ifndef SET_SMTH
#   undef IRON_DATA
#  endif
# endif /* defined IRON_DATA || defined ALL_DATA */

#endif /* BIOLOG_BEC || BIOLOG_BEC2 */

  end subroutine init_arrays_bgc_forces  !]

!----------------------------------------------------------------------

#endif /* BIOLOGY_BEC2  || MARBL */
end module bgc_forces
