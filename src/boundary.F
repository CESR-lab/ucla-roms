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
      real zeta_west(0:Mm+1), zeta_west_dt(0:Mm+1,2)
#  endif
#  ifdef M2_FRC_BRY
      real ubar_west(0:Mm+1), ubar_west_dt(0:Mm+1,2),
     &     vbar_west(0:Mm+1), vbar_west_dt(0:Mm+1,2)
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      real u_west(0:Mm+1,N), u_west_dt(0:Mm+1,N,2),
     &     v_west(0:Mm+1,N), v_west_dt(0:Mm+1,N,2),
     &     w_west(0:Mm+1,N), w_west_dt(0:Mm+1,N,2)
#   endif
#   ifdef T_FRC_BRY
      real t_west(0:Mm+1,N,NT), t_west_dt(0:Mm+1,N,2,NT)
#   endif
#  endif
# endif

# ifdef OBC_EAST
#  ifdef Z_FRC_BRY
      real zeta_east(0:Mm+1), zeta_east_dt(0:Mm+1,2)
#  endif
#  ifdef M2_FRC_BRY
      real ubar_east(0:Mm+1), ubar_east_dt(0:Mm+1,2),
     &     vbar_east(0:Mm+1), vbar_east_dt(0:Mm+1,2)
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      real u_east(0:Mm+1,N), u_east_dt(0:Mm+1,N,2),
     &     v_east(0:Mm+1,N), v_east_dt(0:Mm+1,N,2),
     &     w_east(0:Mm+1,N), w_east_dt(0:Mm+1,N,2)
#   endif
#   ifdef T_FRC_BRY
      real t_east(0:Mm+1,N,NT), t_east_dt(0:Mm+1,N,2,NT)
#   endif
#  endif
# endif

# ifdef OBC_SOUTH
#  ifdef Z_FRC_BRY
      real zeta_south(0:Lm+1), zeta_south_dt(0:Lm+1,2)
#  endif
#  ifdef M2_FRC_BRY
      real ubar_south(0:Lm+1), ubar_south_dt(0:Lm+1,2),
     &     vbar_south(0:Lm+1), vbar_south_dt(0:Lm+1,2)
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      real u_south(0:Lm+1,N), u_south_dt(0:Lm+1,N,2),
     &     v_south(0:Lm+1,N), v_south_dt(0:Lm+1,N,2),
     &     w_south(0:Lm+1,N), w_south_dt(0:Lm+1,N,2)
#   endif
#   ifdef T_FRC_BRY
      real t_south(0:Lm+1,N,NT), t_south_dt(0:Lm+1,N,2,NT)
#   endif
#  endif
# endif

# ifdef OBC_NORTH
#  ifdef Z_FRC_BRY
      real zeta_north(0:Lm+1), zeta_north_dt(0:Lm+1,2)
#  endif
#  ifdef M2_FRC_BRY
      real ubar_north(0:Lm+1), ubar_north_dt(0:Lm+1,2),
     &     vbar_north(0:Lm+1), vbar_north_dt(0:Lm+1,2)
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      real u_north(0:Lm+1,N), u_north_dt(0:Lm+1,N,2),
     &     v_north(0:Lm+1,N), v_north_dt(0:Lm+1,N,2),
     &     w_north(0:Lm+1,N), w_north_dt(0:Lm+1,N,2)
#   endif
#   ifdef T_FRC_BRY
      real t_north(0:Lm+1,N,NT), t_north_dt(0:Lm+1,N,2,NT)
#   endif
#  endif
# endif

      end module boundary
