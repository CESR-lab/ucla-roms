      module bgc_forces
      ! formerly bgc_forces.h

#include "cppdefs.opt"
#ifdef BIOLOGY_BEC2

      use param
      implicit none

#ifdef PCO2AIR_FORCING
! pCO2air concentration
! ------- -------------
!     pCO2air: pCO2air concentraion [ppm]

      real,allocatable,dimension(:,:) :: pco2air
!      common /frc_pco2air/ pco2air
CSDISTRIBUTE_RESHAPE  pCO2air(BLOCK_PATTERN,*) BLOCK_CLAUSE
# if defined PCO2AIR_DATA || defined ALL_DATA
# ifndef SET_SMTH
#  undef PCO2AIR_DATA
# endif
! DevinD: pco2airg not needed, flx_data in bgc.F instead
!      real pco2airg(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE  pco2airg(BLOCK_PATTERN,*) BLOCK_CLAUSE
!      common /pco2airg_dat/ pco2airg
      real :: pco2air_cycle, pco2air_time(2)
      integer :: pco2air_ncycle,  pco2air_rec, itpco2air, ntpco2air,
     &           pco2air_file_id, pco2air_id,  pco2air_tid
!      common /pco2airdat/ pco2air_cycle,    pco2air_time,
!     &        pco2air_ncycle,  pco2air_rec, itpco2air, ntpco2air,
!     &        pco2air_file_id, pco2air_id,  pco2air_tid
# endif
#endif /* PCO2AIR_FORCING */

!--------- NHY_FORCING_START

#ifdef NHY_FORCING
! NHY flux
! --- ----

      real,allocatable,dimension(:,:) :: nhy
!      common /frc_nhy/ nhy
CSDISTRIBUTE_RESHAPE  nhy(BLOCK_PATTERN,*) BLOCK_CLAUSE
# if defined NHY_DATA || defined ALL_DATA
# ifndef SET_SMTH
#  undef NHY_DATA
# endif
      real,allocatable,dimension(:,:,:) :: nhyg
CSDISTRIBUTE_RESHAPE  nhyg(BLOCK_PATTERN,*) BLOCK_CLAUSE
!      common /nhyg_dat/ nhyg
      real    :: nhy_cycle, nhy_time(2)
      integer :: nhy_ncycle,  nhy_rec, itnhy, ntnhy,
     &        nhy_file_id, nhy_id, nhy_tid
!      common /nhydat/ nhy_cycle, nhy_time,
!     &        nhy_ncycle, nhy_rec, itnhy, ntnhy,
!     &        nhy_file_id, nhy_id, nhy_tid
# endif
#endif /* NHY_FORCING */
!--------- NHY_FORCING_END

#ifdef NOX_FORCING
! NOX flux
! --- ----

      real,allocatable,dimension(:,:) :: nox
!      common /frc_nox/ nox
CSDISTRIBUTE_RESHAPE  nox(BLOCK_PATTERN,*) BLOCK_CLAUSE
# if defined NOX_DATA || defined ALL_DATA
# ifndef SET_SMTH
#  undef NOX_DATA
# endif
      real,allocatable,dimension(:,:,:) :: noxg
CSDISTRIBUTE_RESHAPE  noxg(BLOCK_PATTERN,*) BLOCK_CLAUSE
!      common /noxg_dat/ noxg
      real    :: nox_cycle, nox_time(2)
      integer :: nox_ncycle,  nox_rec, itnox, ntnox,
     &           nox_file_id, nox_id, nox_tid
!      common /noxdat/ nox_cycle, nox_time,
!     &        nox_ncycle, nox_rec, itnox, ntnox,
!     &        nox_file_id, nox_id, nox_tid
# endif
#endif /* NOX_FORCING */

#if defined BIOLOGY_BEC || defined BIOLOGY_BEC2
! dust flux
! ---- ----
!      dust: dust flux [kg m-2 s-1]

      real,allocatable,dimension(:,:) :: dust
CSDISTRIBUTE_RESHAPE  dust(BLOCK_PATTERN) BLOCK_CLAUSE
!           common /frc_dust/dust
# if defined DUST_DATA || defined ALL_DATA
#  ifndef SET_SMTH
#   undef DUST_DATA
#  endif
! DevinD: dustg not needed, flx_data in bgc.F instead
!      real dustg(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE  dustg(BLOCK_PATTERN,*) BLOCK_CLAUSE
!      common /dustg_dat/dustg
      real :: dust_cycle, dust_time(2)
      integer :: dust_ncycle,  dust_rec, itdust, ntdust,
     &           dust_file_id, dust_id,  dust_tid
!      common /dustdat/ dust_cycle,    dust_time,
!     &        dust_ncycle,  dust_rec, itdust, ntdust,
!     &        dust_file_id, dust_id,  dust_tid
# endif /* defined DUST_DATA || defined ALL_DATA */

! iron flux
! ---- ----
!     iron: iron flux [nmol cm-2 s-1]

      real,allocatable,dimension(:,:) :: iron
CSDISTRIBUTE_RESHAPE  iron(BLOCK_PATTERN) BLOCK_CLAUSE
!           common /frc_iron/iron
# if defined IRON_DATA || defined ALL_DATA
#  ifndef SET_SMTH
#   undef IRON_DATA
#  endif
! DevinD: irong not needed, flx_data in bgc.F instead
!      real irong(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE  irong(BLOCK_PATTERN,*) BLOCK_CLAUSE
!      common /irong_dat/irong
      real :: iron_cycle, iron_time(2)
      integer :: iron_ncycle,  iron_rec, itiron, ntiron,
     &           iron_file_id, iron_id,  iron_tid
!      common /irondat/ iron_cycle,    iron_time,
!     &        iron_ncycle,  iron_rec, itiron, ntiron,
!     &        iron_file_id, iron_id,  iron_tid
# endif /* defined IRON_DATA || defined ALL_DATA */

# if defined  RIVER_LOAD_N || defined RIVER_LOAD_P
! river input: DIN, DIP from rivers as surface flux
!        din_river, dip_river
      real,allocatable,dimension(:,:) :: din_river, dip_river
!      common /frc_river/ din_river, dip_river
CSDISTRIBUTE_RESHAPE  din_river(BLOCK_PATTERN) BLOCK_CLAUSE
CSDISTRIBUTE_RESHAPE  dip_river(BLOCK_PATTERN) BLOCK_CLAUSE
#  if defined RIVER_DATA || defined ALL_DATA
#   ifndef SET_SMTH
#    undef RIVER_DATA
#   endif
      real,allocatable,dimension(:,:,:) :: din_riverg, dip_riverg
CSDISTRIBUTE_RESHAPE  din_riverg(BLOCK_PATTERN,*) BLOCK_CLAUSE
CSDISTRIBUTE_RESHAPE  dip_riverg(BLOCK_PATTERN,*) BLOCK_CLAUSE
!      common /riverg_dat/din_riverg,dip_riverg
      real :: din_river_cycle, din_river_time(2),
     &        dip_river_cycle, dip_river_time(2)
      integer :: din_river_ncycle, din_river_rec, itdin_river, ntdin_river,
     &        din_river_file_id, din_river_id, din_river_tid,
     &        dip_river_ncycle, dip_river_rec, itdip_river, ntdip_river,
     &        dip_river_file_id, dip_river_id,  dip_river_tid
!      common/riverdat/
!     &        din_river_cycle, din_river_time,
!     &        dip_river_cycle, dip_river_time,
!     &        din_river_ncycle, din_river_rec, itdin_river, ntdin_river,
!     &        din_river_file_id, din_river_id, din_river_tid,
!     &        dip_river_ncycle, dip_river_rec, itdip_river, ntdip_river,
!     &        dip_river_file_id, dip_river_id, dip_river_tid
#  endif /* defined RIVER_DATA || defined ALL_DATA */
# endif /* RIVER_LOAD_N || RIVER_LOAD_P */
#endif /* BIOLOG_BEC || BIOLOG_BEC2 */

      contains

!----------------------------------------------------------------------
      subroutine init_arrays_bgc_forces  ![
      implicit none

#ifdef PCO2AIR_FORCING
      allocate( pco2air(GLOBAL_2D_ARRAY) )
#endif /* PCO2AIR_FORCING */

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

#if defined BIOLOGY_BEC || defined BIOLOGY_BEC2
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

# if defined  RIVER_LOAD_N || defined RIVER_LOAD_P
      allocate( din_river(GLOBAL_2D_ARRAY), dip_river(GLOBAL_2D_ARRAY) )
#  if defined RIVER_DATA || defined ALL_DATA
#   ifndef SET_SMTH
#    undef RIVER_DATA
#   endif
      allocate( din_riverg(GLOBAL_2D_ARRAY), dip_riverg(GLOBAL_2D_ARRAY) )
#  endif /* defined RIVER_DATA || defined ALL_DATA */
# endif /* RIVER_LOAD_N || RIVER_LOAD_P */
#endif /* BIOLOG_BEC || BIOLOG_BEC2 */

      end subroutine init_arrays_bgc_forces  !]

!----------------------------------------------------------------------

#endif /* BIOLOGY_BEC2 */
      end module bgc_forces
