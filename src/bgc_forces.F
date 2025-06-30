      module bgc_forces
      ! formerly bgc_forces.h

#include "cppdefs.opt"
#if defined(BIOLOGY_BEC2) || defined(MARBL)

      use param
      implicit none

#ifdef PCO2AIR_FORCING
! pCO2air concentration
! ------- -------------
!     pCO2air: pCO2air concentraion [ppm]
      real,allocatable,dimension(:,:) :: pco2air
#ifdef MARBL
      real,allocatable,dimension(:,:) :: pco2air_alt
#endif
CSDISTRIBUTE_RESHAPE  pCO2air(BLOCK_PATTERN,*) BLOCK_CLAUSE
# if defined PCO2AIR_DATA || defined ALL_DATA
# ifndef SET_SMTH
#  undef PCO2AIR_DATA
# endif
CSDISTRIBUTE_RESHAPE  pco2airg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real :: pco2air_cycle, pco2air_time(2)
      integer :: pco2air_ncycle,  pco2air_rec, itpco2air, ntpco2air,
     &           pco2air_file_id, pco2air_id,  pco2air_tid
# endif
#endif /* PCO2AIR_FORCING */

#if defined DAILYPAR_PHOTOINHIBITION || defined DAILYPAR_BEC
! daily avg swrad
! ------- -------------
      real,allocatable,dimension(:,:) :: swrad_avg
CSDISTRIBUTE_RESHAPE  swrad_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
# if defined SWRAD_AVG_DATA || defined ALL_DATA
# ifndef SET_SMTH
#  undef SWRAD_AVG_DATA
# endif
CSDISTRIBUTE_RESHAPE  swrad_avgg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real :: swrad_avg_cycle, swrad_avg_time(2)
      integer :: swrad_avg_ncycle,  swrad_avg_rec, itswrad_avg, ntswrad_avg,
     &           swrad_avg_file_id, swrad_avg_id,  swrad_avg_tid
# endif
#endif /* DAILYPAR_PHOTOINHIBITION || defined DAILYPAR_BEC */


#ifdef NHY_FORCING
! NHY flux
! --- ----
      real,allocatable,dimension(:,:) :: nhy
CSDISTRIBUTE_RESHAPE  nhy(BLOCK_PATTERN,*) BLOCK_CLAUSE
# if defined NHY_DATA || defined ALL_DATA
# ifndef SET_SMTH
#  undef NHY_DATA
# endif
      real,allocatable,dimension(:,:,:) :: nhyg
CSDISTRIBUTE_RESHAPE  nhyg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real    :: nhy_cycle, nhy_time(2)
      integer :: nhy_ncycle,  nhy_rec, itnhy, ntnhy,
     &           nhy_file_id, nhy_id, nhy_tid
# endif
#endif /* NHY_FORCING */

#ifdef NOX_FORCING
! NOX flux
! --- ----
      real,allocatable,dimension(:,:) :: nox
CSDISTRIBUTE_RESHAPE  nox(BLOCK_PATTERN,*) BLOCK_CLAUSE
# if defined NOX_DATA || defined ALL_DATA
# ifndef SET_SMTH
#  undef NOX_DATA
# endif
      real,allocatable,dimension(:,:,:) :: noxg
CSDISTRIBUTE_RESHAPE  noxg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real    :: nox_cycle, nox_time(2)
      integer :: nox_ncycle,  nox_rec, itnox, ntnox,
     &           nox_file_id, nox_id, nox_tid
# endif
#endif /* NOX_FORCING */

#if defined BIOLOGY_BEC || defined BIOLOGY_BEC2 || defined(MARBL)
! dust flux
! ---- ----
!      dust: dust flux [kg m-2 s-1]
      real,allocatable,dimension(:,:) :: dust
CSDISTRIBUTE_RESHAPE  dust(BLOCK_PATTERN) BLOCK_CLAUSE
# if defined DUST_DATA || defined ALL_DATA
#  ifndef SET_SMTH
#   undef DUST_DATA
#  endif
CSDISTRIBUTE_RESHAPE  dustg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real :: dust_cycle, dust_time(2)
      integer :: dust_ncycle,  dust_rec, itdust, ntdust,
     &           dust_file_id, dust_id,  dust_tid
# endif /* defined DUST_DATA || defined ALL_DATA */

! iron flux
! ---- ----
!     iron: iron flux [nmol cm-2 s-1]
      real,allocatable,dimension(:,:) :: iron
CSDISTRIBUTE_RESHAPE  iron(BLOCK_PATTERN) BLOCK_CLAUSE
# if defined IRON_DATA || defined ALL_DATA
#  ifndef SET_SMTH
#   undef IRON_DATA
#  endif
CSDISTRIBUTE_RESHAPE  irong(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real :: iron_cycle, iron_time(2)
      integer :: iron_ncycle,  iron_rec, itiron, ntiron,
     &           iron_file_id, iron_id,  iron_tid
# endif /* defined IRON_DATA || defined ALL_DATA */

#endif /* BIOLOG_BEC || MARBL */

      contains

!----------------------------------------------------------------------
      subroutine init_arrays_bgc_forces  ![
      implicit none

#ifdef PCO2AIR_FORCING
      allocate( pco2air(GLOBAL_2D_ARRAY) )
#ifdef MARBL
      allocate( pco2air_alt(GLOBAL_2D_ARRAY) )
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
