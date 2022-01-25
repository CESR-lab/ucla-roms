      module boundary

! Meaning of variables in "bry_time_vars": 
!-------- -- --------- -- -----------------
! itbry - array index corresponding to "new" data, either 1 or 2,
!         so model time is bounded by data time as
!
!         bry_time(3-itbry) + bry_ncycle*bry_cycle <= time
!                   time < bry_time(itbry) + bry_ncycle*bry_cycle
!
! bry_rec - record within netCDF file corresponding to bry_time(itbry);
! ntbry - total number of records in the file;
! ibry  - file index within the sequence of files, 1<= ibry <= max_bry;
! max_bry - total number of files in the sequence;
! bry_id - netCDF ID of current file;
! bry_time_id - netCDF variable ID for time variable;
! All others are netCDF variable IDs for the corresponding variables.

#include "cppdefs.opt"

      use param
      implicit none

#ifndef ANA_BRY
      real(kind=8) bry_cycle, bry_time(2)
      integer itbry, ntbry, bry_ncycle, ibry, max_bry,
     &                  bry_rec,  bry_id, bry_time_id

# ifdef OBC_WEST
#  ifdef Z_FRC_BRY
      integer zeta_west_id
#  endif
#  ifdef M2_FRC_BRY
      integer ubar_west_id, vbar_west_id
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      integer u_west_id, v_west_id, w_west_id
#   endif
#   ifdef T_FRC_BRY
      integer t_west_id(NT)
#   endif
#  endif
# endif

# ifdef OBC_EAST
#  ifdef Z_FRC_BRY
      integer zeta_east_id
#  endif
#  ifdef M2_FRC_BRY
      integer ubar_east_id, vbar_east_id
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      integer u_east_id, v_east_id, w_east_id
#   endif
#   ifdef T_FRC_BRY
      integer t_east_id(NT)
#   endif
#  endif
# endif

# ifdef OBC_SOUTH
#  ifdef Z_FRC_BRY
      integer zeta_south_id
#  endif
#  ifdef M2_FRC_BRY
      integer ubar_south_id, vbar_south_id
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      integer u_south_id, v_south_id, w_south_id
#   endif
#   ifdef T_FRC_BRY
      integer t_south_id(NT)
#   endif
#  endif
# endif

# ifdef OBC_NORTH
#  ifdef Z_FRC_BRY
      integer zeta_north_id
#  endif
#  ifdef M2_FRC_BRY
      integer ubar_north_id, vbar_north_id
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      integer u_north_id, v_north_id, w_north_id
#   endif
#   ifdef T_FRC_BRY
      integer t_north_id(NT)
#   endif
#  endif
# endif
#endif  /* ~ANA_BRY */


# ifdef OBC_WEST
#  ifdef Z_FRC_BRY
      real,allocatable,dimension(:)   :: zeta_west
      real,allocatable,dimension(:,:) :: zeta_west_dt
#  endif
#  ifdef M2_FRC_BRY
      real,allocatable,dimension(:)   :: ubar_west, vbar_west
      real,allocatable,dimension(:,:) :: ubar_west_dt, vbar_west_dt
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      real,allocatable,dimension(:,:)   :: u_west, v_west, w_west
      real,allocatable,dimension(:,:,:) :: u_west_dt, v_west_dt, w_west_dt
#   endif
#   ifdef T_FRC_BRY
      real,allocatable,dimension(:,:,:)   :: t_west
      real,allocatable,dimension(:,:,:,:) :: t_west_dt
#   endif
#  endif
# endif

# ifdef OBC_EAST
#  ifdef Z_FRC_BRY
      real,allocatable,dimension(:)   :: zeta_east
      real,allocatable,dimension(:,:) :: zeta_east_dt
#  endif
#  ifdef M2_FRC_BRY
      real,allocatable,dimension(:)   :: ubar_east, vbar_east
      real,allocatable,dimension(:,:) :: ubar_east_dt, vbar_east_dt
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      real,allocatable,dimension(:,:)   :: u_east, v_east, w_east
      real,allocatable,dimension(:,:,:) :: u_east_dt, v_east_dt, w_east_dt
#   endif
#   ifdef T_FRC_BRY
      real,allocatable,dimension(:,:,:)   :: t_east
      real,allocatable,dimension(:,:,:,:) :: t_east_dt
#   endif
#  endif
# endif

# ifdef OBC_SOUTH
#  ifdef Z_FRC_BRY
      real,allocatable,dimension(:)   :: zeta_south
      real,allocatable,dimension(:,:) :: zeta_south_dt
#  endif
#  ifdef M2_FRC_BRY
      real,allocatable,dimension(:)   :: ubar_south, vbar_south
      real,allocatable,dimension(:,:) :: ubar_south_dt, vbar_south_dt
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      real,allocatable,dimension(:,:)   :: u_south, v_south, w_south
      real,allocatable,dimension(:,:,:) :: u_south_dt, v_south_dt, w_south_dt
#   endif
#   ifdef T_FRC_BRY
      real,allocatable,dimension(:,:,:)   :: t_south
      real,allocatable,dimension(:,:,:,:) :: t_south_dt
#   endif
#  endif
# endif

# ifdef OBC_NORTH
#  ifdef Z_FRC_BRY
      real,allocatable,dimension(:)   :: zeta_north
      real,allocatable,dimension(:,:) :: zeta_north_dt
#  endif
#  ifdef M2_FRC_BRY
      real,allocatable,dimension(:)   :: ubar_north, vbar_north
      real,allocatable,dimension(:,:) :: ubar_north_dt, vbar_north_dt
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      real,allocatable,dimension(:,:)   :: u_north, v_north, w_north
      real,allocatable,dimension(:,:,:) :: u_north_dt, v_north_dt, w_north_dt
#   endif
#   ifdef T_FRC_BRY
      real,allocatable,dimension(:,:,:)   :: t_north
      real,allocatable,dimension(:,:,:,:) :: t_north_dt
#   endif
#  endif
# endif

      contains

      !----------------------------------------------------------------
      subroutine init_arrays_boundary  ![
      implicit none

# ifdef OBC_WEST
#  ifdef Z_FRC_BRY
      allocate( zeta_west(0:Mm+1), zeta_west_dt(0:Mm+1,2) )
#  endif
#  ifdef M2_FRC_BRY
      allocate( ubar_west(0:Mm+1), ubar_west_dt(0:Mm+1,2),
     &     vbar_west(0:Mm+1), vbar_west_dt(0:Mm+1,2) )
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      allocate( u_west(0:Mm+1,N), u_west_dt(0:Mm+1,N,2),
     &     v_west(0:Mm+1,N), v_west_dt(0:Mm+1,N,2),
     &     w_west(0:Mm+1,N), w_west_dt(0:Mm+1,N,2) )
#   endif
#   ifdef T_FRC_BRY
      allocate( t_west(0:Mm+1,N,NT), t_west_dt(0:Mm+1,N,2,NT) )
#   endif
#  endif
# endif

# ifdef OBC_EAST
#  ifdef Z_FRC_BRY
      allocate( zeta_east(0:Mm+1), zeta_east_dt(0:Mm+1,2) )
#  endif
#  ifdef M2_FRC_BRY
      allocate( ubar_east(0:Mm+1), ubar_east_dt(0:Mm+1,2),
     &     vbar_east(0:Mm+1), vbar_east_dt(0:Mm+1,2) )
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      allocate( u_east(0:Mm+1,N), u_east_dt(0:Mm+1,N,2),
     &     v_east(0:Mm+1,N), v_east_dt(0:Mm+1,N,2),
     &     w_east(0:Mm+1,N), w_east_dt(0:Mm+1,N,2) )
#   endif
#   ifdef T_FRC_BRY
      allocate( t_east(0:Mm+1,N,NT), t_east_dt(0:Mm+1,N,2,NT) )
#   endif
#  endif
# endif

# ifdef OBC_SOUTH
#  ifdef Z_FRC_BRY
      allocate( zeta_south(0:Lm+1), zeta_south_dt(0:Lm+1,2) )
#  endif
#  ifdef M2_FRC_BRY
      allocate( ubar_south(0:Lm+1), ubar_south_dt(0:Lm+1,2),
     &     vbar_south(0:Lm+1), vbar_south_dt(0:Lm+1,2) )
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      allocate( u_south(0:Lm+1,N), u_south_dt(0:Lm+1,N,2),
     &     v_south(0:Lm+1,N), v_south_dt(0:Lm+1,N,2),
     &     w_south(0:Lm+1,N), w_south_dt(0:Lm+1,N,2) )
#   endif
#   ifdef T_FRC_BRY
      allocate( t_south(0:Lm+1,N,NT), t_south_dt(0:Lm+1,N,2,NT) )
#   endif
#  endif
# endif

# ifdef OBC_NORTH
#  ifdef Z_FRC_BRY
      allocate( zeta_north(0:Lm+1), zeta_north_dt(0:Lm+1,2) )
#  endif
#  ifdef M2_FRC_BRY
      allocate( ubar_north(0:Lm+1), ubar_north_dt(0:Lm+1,2),
     &     vbar_north(0:Lm+1), vbar_north_dt(0:Lm+1,2) )
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      allocate( u_north(0:Lm+1,N), u_north_dt(0:Lm+1,N,2),
     &     v_north(0:Lm+1,N), v_north_dt(0:Lm+1,N,2),
     &     w_north(0:Lm+1,N), w_north_dt(0:Lm+1,N,2) )
#   endif
#   ifdef T_FRC_BRY
      allocate( t_north(0:Lm+1,N,NT), t_north_dt(0:Lm+1,N,2,NT) )
#   endif
#  endif
# endif

      end subroutine init_arrays_boundary !]

      !----------------------------------------------------------------

      end module boundary
