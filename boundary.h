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

#ifndef ANA_BRY
      real(kind=8) bry_cycle, bry_time(2)
      integer itbry, ntbry, bry_ncycle, ibry, max_bry,
     &                  bry_rec,  bry_id, bry_time_id
      common /bry_time_vars/ bry_cycle, bry_time,
     &        itbry, ntbry, bry_ncycle, ibry, max_bry, 
     &                  bry_rec,  bry_id, bry_time_id

# ifdef OBC_WEST
#  ifdef Z_FRC_BRY
      integer zeta_west_id
      common /bry_time_vars/ zeta_west_id
#  endif
#  ifdef M2_FRC_BRY
      integer ubar_west_id, vbar_west_id
      common /bry_time_vars/ ubar_west_id, vbar_west_id
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      integer u_west_id, v_west_id, w_west_id
      common /bry_time_vars/ u_west_id, v_west_id, w_west_id
#   endif
#   ifdef T_FRC_BRY
      integer t_west_id(NT)
      common /bry_time_vars/ t_west_id
#   endif
#  endif
# endif

# ifdef OBC_EAST
#  ifdef Z_FRC_BRY
      integer zeta_east_id
      common /bry_time_vars/ zeta_east_id
#  endif
#  ifdef M2_FRC_BRY
      integer ubar_east_id, vbar_east_id
      common /bry_time_vars/ ubar_east_id, vbar_east_id
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      integer u_east_id, v_east_id, w_east_id
      common /bry_time_vars/ u_east_id, v_east_id, w_east_id
#   endif
#   ifdef T_FRC_BRY
      integer t_east_id(NT)
      common /bry_time_vars/ t_east_id
#   endif
#  endif
# endif

# ifdef OBC_SOUTH
#  ifdef Z_FRC_BRY
      integer zeta_south_id
      common /bry_time_vars/ zeta_south_id
#  endif
#  ifdef M2_FRC_BRY
      integer ubar_south_id, vbar_south_id
      common /bry_time_vars/ ubar_south_id, vbar_south_id
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      integer u_south_id, v_south_id, w_south_id
      common /bry_time_vars/ u_south_id, v_south_id, w_south_id
#   endif
#   ifdef T_FRC_BRY
      integer t_south_id(NT)
      common /bry_time_vars/ t_south_id
#   endif
#  endif
# endif

# ifdef OBC_NORTH
#  ifdef Z_FRC_BRY
      integer zeta_north_id
      common /bry_time_vars/ zeta_north_id
#  endif
#  ifdef M2_FRC_BRY
      integer ubar_north_id, vbar_north_id
      common /bry_time_vars/ ubar_north_id, vbar_north_id
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      integer u_north_id, v_north_id, w_north_id
      common /bry_time_vars/ u_north_id, v_north_id, w_north_id
#   endif
#   ifdef T_FRC_BRY
      integer t_north_id(NT)
      common /bry_time_vars/ t_north_id
#   endif
#  endif
# endif
#endif  /* ~ANA_BRY */





# ifdef OBC_WEST
#  ifdef Z_FRC_BRY
      real zeta_west(0:Mm+1), zeta_west_dt(0:Mm+1,2)
      common /bry_west/ zeta_west, zeta_west_dt
#  endif
#  ifdef M2_FRC_BRY
      real ubar_west(0:Mm+1), ubar_west_dt(0:Mm+1,2),
     &     vbar_west(0:Mm+1), vbar_west_dt(0:Mm+1,2)
      common /bry_west/ ubar_west, ubar_west_dt,
     &                  vbar_west, vbar_west_dt
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      real u_west(0:Mm+1,N), u_west_dt(0:Mm+1,N,2),
     &     v_west(0:Mm+1,N), v_west_dt(0:Mm+1,N,2),
     &     w_west(0:Mm+1,N), w_west_dt(0:Mm+1,N,2)
      common /bry_west/ u_west, u_west_dt,
     &                  v_west, v_west_dt,
     &                  w_west, w_west_dt
#   endif
#   ifdef T_FRC_BRY
      real t_west(0:Mm+1,N,NT), t_west_dt(0:Mm+1,N,2,NT)
      common /bry_west/ t_west, t_west_dt
#   endif
#  endif
# endif

# ifdef OBC_EAST
#  ifdef Z_FRC_BRY
      real zeta_east(0:Mm+1), zeta_east_dt(0:Mm+1,2)
      common /bry_east/ zeta_east, zeta_east_dt
#  endif
#  ifdef M2_FRC_BRY
      real ubar_east(0:Mm+1), ubar_east_dt(0:Mm+1,2),
     &     vbar_east(0:Mm+1), vbar_east_dt(0:Mm+1,2)
      common /bry_east/ ubar_east, ubar_east_dt,
     &                  vbar_east, vbar_east_dt
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      real u_east(0:Mm+1,N), u_east_dt(0:Mm+1,N,2),
     &     v_east(0:Mm+1,N), v_east_dt(0:Mm+1,N,2),
     &     w_east(0:Mm+1,N), w_east_dt(0:Mm+1,N,2)
      common /bry_east/ u_east, u_east_dt,
     &                  v_east, v_east_dt,
     &                  w_east, w_east_dt
#   endif
#   ifdef T_FRC_BRY
      real t_east(0:Mm+1,N,NT), t_east_dt(0:Mm+1,N,2,NT)
      common /bry_east/ t_east, t_east_dt
#   endif
#  endif
# endif

# ifdef OBC_SOUTH
#  ifdef Z_FRC_BRY
      real zeta_south(0:Lm+1), zeta_south_dt(0:Lm+1,2)
      common /bry_south/ zeta_south, zeta_south_dt
#  endif
#  ifdef M2_FRC_BRY
      real ubar_south(0:Lm+1), ubar_south_dt(0:Lm+1,2),
     &     vbar_south(0:Lm+1), vbar_south_dt(0:Lm+1,2)
      common /bry_south/ ubar_south, ubar_south_dt,
     &                   vbar_south, vbar_south_dt
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      real u_south(0:Lm+1,N), u_south_dt(0:Lm+1,N,2),
     &     v_south(0:Lm+1,N), v_south_dt(0:Lm+1,N,2),
     &     w_south(0:Lm+1,N), w_south_dt(0:Lm+1,N,2)
      common /bry_south/ u_south, u_south_dt,
     &                   v_south, v_south_dt,
     &                   w_south, w_south_dt
#   endif
#   ifdef T_FRC_BRY
      real t_south(0:Lm+1,N,NT), t_south_dt(0:Lm+1,N,2,NT)
      common /bry_south/ t_south, t_south_dt
#   endif
#  endif
# endif

# ifdef OBC_NORTH
#  ifdef Z_FRC_BRY
      real zeta_north(0:Lm+1), zeta_north_dt(0:Lm+1,2)
      common /bry_north/ zeta_north, zeta_north_dt
#  endif
#  ifdef M2_FRC_BRY
      real ubar_north(0:Lm+1), ubar_north_dt(0:Lm+1,2),
     &     vbar_north(0:Lm+1), vbar_north_dt(0:Lm+1,2)
      common /bry_north/ ubar_north, ubar_north_dt,
     &                   vbar_north, vbar_north_dt
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      real u_north(0:Lm+1,N), u_north_dt(0:Lm+1,N,2),
     &     v_north(0:Lm+1,N), v_north_dt(0:Lm+1,N,2),
     &     w_north(0:Lm+1,N), w_north_dt(0:Lm+1,N,2)
      common /bry_north/ u_north, u_north_dt,
     &                   v_north, v_north_dt,
     &                   w_north, w_north_dt
#   endif
#   ifdef T_FRC_BRY
      real t_north(0:Lm+1,N,NT), t_north_dt(0:Lm+1,N,2,NT)
      common /bry_north/ t_north, t_north_dt
#   endif
#  endif
# endif

