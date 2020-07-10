#include "cppdefs.h"
#if !defined ANA_BRY && ( defined OBC_WEST   || defined OBC_EAST  || \
                          defined OBC_SOUTH  || defined OBC_NORTH )  \
                     && ( defined M2_FRC_BRY || defined Z_FRC_BRY || \
                          defined M3_FRC_BRY || defined T_FRC_BRY )

      subroutine get_bry_all(ierr_ext)     ! Read side boundary forcing
      implicit none                        ! fields from boundary file
      integer ierr_ext
# include "param.h"
# include "scalars.h"
# include "ncvars.h"
# include "netcdf.inc"
# include "boundary.h"
      real dtime
      integer ierr, ierr_all, itrc, lsrc, lenstr
# if defined OBC_WEST  || defined OBC_EAST
     &                          , read_bry_EW
# endif
# if defined OBC_SOUTH || defined OBC_NORTH
     &                          , read_bry_NS
# endif
      logical need_bry
# ifdef OBC_WEST
     &      , west_bry_active
# endif
# ifdef OBC_EAST
     &      , east_bry_active
# endif
# ifdef OBC_SOUTH
     &      , south_bry_active
# endif
# ifdef OBC_NORTH
     &      , north_bry_active
# endif
# ifdef VERBOSE
      write(*,*) 'enter get_bry_all' MYID
# endif

! Determine which side(s) require boundary conditions and check
! whether it is needed at all. Logical variable "need_bry" is used
! as a trigger here: if it remains .false. for the current MPI
! subdomain it does not have edge adjacent to an open boundary.

      need_bry=.false.
# ifdef OBC_WEST
#  ifdef MPI
      west_bry_active=.not.west_exchng
#  else
      west_bry_active=.true.
#  endif
      need_bry=need_bry .or. west_bry_active
# endif
# ifdef OBC_EAST
#  ifdef MPI
      east_bry_active=.not.east_exchng
#  else
      east_bry_active=.true.
#  endif
      need_bry=need_bry .or. east_bry_active
# endif
# ifdef OBC_SOUTH
#  ifdef MPI
      south_bry_active=.not.south_exchng
#  else
      south_bry_active=.true.
#  endif
      need_bry=need_bry .or. south_bry_active
# endif
# ifdef OBC_NORTH
#  ifdef MPI
      north_bry_active=.not.north_exchng
#  else
      north_bry_active=.true.
#  endif
      need_bry=need_bry .or. north_bry_active
# endif

! Initialization:  The overall logic is similar to "get_smth.F" with
!----------------  the exception that presence of all needed variables
! in netCDF file is mandatory, so there is no search for variables
! across multiple input files: if more than one, they are considered
! solely as time sequence. Other than that, as in "get_smth.F" -- check
! whether boundary forcing file is already in open state, open it if
! not. Find and save netCDF IDs for relevant variables.  Determine
! whether there is cycling to reuse the input data and find cycling
! period "bry_cycle", set initial cycling index "bry_ncycle" and
! record index "bry_rec", time index "itbry" and both time record
! bounds to large negative artificial values, so that it will trigger
! the logic in reading part below.


! Error signalling is as follows:
! "ierr_ext" is IN/OUT external error flag which remains unchanged
!            and untouched if no errors occur within this subroutine;
!            it will be changed by this routine only if
!                  (i) its incoming value is no-error,  AND
!                 (ii) at least one error occurs inside this routine;
! "ierr_all" is cumulative "OR" status of all errors occurring here;
!            essentially this variable "locks" on the first error
!            occurring solely within this routine;


      if (iic == ntstart) then
        ibry=1 ; itbry=1 ; bry_time(1)=-1.E+33
        bry_rec=-1       ; bry_time(2)=-1.E+33
      endif

 100  lsrc=lenstr(bry_file(ibry))
      ierr=nf_noerr ; ierr_all=nf_noerr

      if (bry_rec < 0) then  !<-- initialization or file switch
        if (bry_id == -1) then
          ierr=nf_open(bry_file(ibry)(1:lsrc), nf_nowrite, bry_id)
          if (ierr /= nf_noerr) write(*,'(/1x,4A/12x,A/)')
     &          '### ERROR: get_bry_all :: Cannot open file ''',
     &              bry_file(ibry)(1:lsrc),   ''' for reading.',
     &                                      nf_strerror(ierr)

        endif
        if (ierr == nf_noerr) then
          ierr=nf_inq_varid(bry_id, 'bry_time',  bry_time_id)
          if (ierr == nf_noerr) then
            call find_record(bry_id, bry_time_id,  bry_rec, ntbry,
     &              bry_ncycle, bry_cycle,  ibry, max_bry,  ierr)

            if (bry_rec < 0) then
              ierr=nf_close(bry_id) ; bry_id=-1
              ibry=ibry+1 ;  goto 100  !--> proceed to the next file
            endif

          else
            write(*,1) 'bry_time', bry_file(ibry)(1:lsrc)
          endif
          if (ierr_all == nf_noerr) ierr_all=ierr



# ifdef OBC_WEST
          if (west_bry_active) then
#  ifdef Z_FRC_BRY
            ierr=nf_inq_varid(bry_id, 'zeta_west', zeta_west_id)
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1) 'zeta_west', bry_file(ibry)(1:lsrc)
            endif
#  endif
#  ifdef M2_FRC_BRY
            ierr=nf_inq_varid(bry_id, 'ubar_west', ubar_west_id)
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1) 'ubar_west', bry_file(ibry)(1:lsrc)
            endif

            ierr=nf_inq_varid(bry_id, 'vbar_west', vbar_west_id)
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1) 'vbar_west', bry_file(ibry)(1:lsrc)
            endif
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
            ierr=nf_inq_varid(bry_id,    'u_west', u_west_id)
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1) 'u_west', bry_file(ibry)(1:lsrc)
            endif

            ierr=nf_inq_varid(bry_id,    'v_west', v_west_id)
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1) 'v_west', bry_file(ibry)(1:lsrc)
            endif
#    ifdef NHMG
            ierr=nf_inq_varid(bry_id,    'w_west', w_west_id)
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1) 'w_west', bry_file(ibry)(1:lsrc)
            endif
#    endif
#   endif
#   ifdef T_FRC_BRY
            ierr=nf_inq_varid(bry_id,'temp_west', t_west_id(itemp))
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1) 'temp_west', bry_file(ibry)(1:lsrc)
            endif
#    ifdef SALINITY
            ierr=nf_inq_varid(bry_id, 'salt_west', t_west_id(isalt))
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1) 'salt_west', bry_file(ibry)(1:lsrc)
            endif
#    endif
#   endif
#  endif
          endif !<-- west_bry_active
# endif


# ifdef OBC_EAST
          if (east_bry_active) then
#  ifdef Z_FRC_BRY
            ierr=nf_inq_varid(bry_id, 'zeta_east', zeta_east_id)
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1) 'zeta_east', bry_file(ibry)(1:lsrc)
            endif
#  endif
#  ifdef M2_FRC_BRY
            ierr=nf_inq_varid(bry_id, 'ubar_east', ubar_east_id)
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1) 'ubar_east', bry_file(ibry)(1:lsrc)
            endif

            ierr=nf_inq_varid (bry_id, 'vbar_east', vbar_east_id)
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1) 'vbar_east', bry_file(ibry)(1:lsrc)
            endif
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
            ierr=nf_inq_varid(bry_id,    'u_east',   u_east_id)
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1)   'u_east', bry_file(ibry)(1:lsrc)
            endif

            ierr=nf_inq_varid(bry_id,    'v_east',   v_east_id)
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1)   'v_east', bry_file(ibry)(1:lsrc)
            endif
#    ifdef NHMG
            ierr=nf_inq_varid(bry_id,    'w_east', w_east_id)
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1) 'w_east', bry_file(ibry)(1:lsrc)
            endif
#    endif
#   endif
#   ifdef T_FRC_BRY
            ierr=nf_inq_varid(bry_id, 'temp_east', t_east_id(itemp))
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1) 'temp_east', bry_file(ibry)(1:lsrc)
            endif
#    ifdef SALINITY
            ierr=nf_inq_varid(bry_id, 'salt_east', t_east_id(isalt))
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1) 'salt_east', bry_file(ibry)(1:lsrc)
            endif
#    endif
#   endif
#  endif
          endif  !<-- east_bry_active
# endif


# ifdef OBC_SOUTH
          if (south_bry_active) then
#  ifdef Z_FRC_BRY
            ierr=nf_inq_varid(bry_id, 'zeta_south', zeta_south_id)
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1) 'zeta_south', bry_file(ibry)(1:lsrc)
            endif
#  endif
#  ifdef M2_FRC_BRY
            ierr=nf_inq_varid(bry_id, 'ubar_south', ubar_south_id)
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1) 'ubar_south', bry_file(ibry)(1:lsrc)
            endif

            ierr=nf_inq_varid(bry_id, 'vbar_south', vbar_south_id)
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1) 'vbar_south', bry_file(ibry)(1:lsrc)
            endif
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
            ierr=nf_inq_varid(bry_id,    'u_south',  u_south_id)
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1)    'u_south', bry_file(ibry)(1:lsrc)
            endif

            ierr=nf_inq_varid(bry_id,    'v_south',  v_south_id)
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1)    'v_south', bry_file(ibry)(1:lsrc)
            endif
#    ifdef NHMG
            ierr=nf_inq_varid(bry_id,    'w_south', w_south_id)
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1) 'w_south', bry_file(ibry)(1:lsrc)
            endif
#    endif
#   endif
#   ifdef T_FRC_BRY
            ierr=nf_inq_varid(bry_id,'temp_south',t_south_id(itemp))
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1) 'temp_south', bry_file(ibry)(1:lsrc)
            endif
#    ifdef SALINITY
            ierr=nf_inq_varid(bry_id,'salt_south', t_south_id(isalt))
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1) 'salt_south', bry_file(ibry)(1:lsrc)
            endif
#    endif
#   endif
#  endif
          endif  !<-- south_bry_active
# endif


# ifdef OBC_NORTH
          if (north_bry_active) then
#  ifdef Z_FRC_BRY
            ierr=nf_inq_varid(bry_id, 'zeta_north', zeta_north_id)
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1) 'zeta_north', bry_file(ibry)(1:lsrc)
            endif
#  endif
#  ifdef M2_FRC_BRY
            ierr=nf_inq_varid(bry_id, 'ubar_north', ubar_north_id)
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1) 'ubar_north', bry_file(ibry)(1:lsrc)
            endif

            ierr=nf_inq_varid(bry_id, 'vbar_north', vbar_north_id)
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1) 'vbar_north', bry_file(ibry)(1:lsrc)
            endif
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
            ierr=nf_inq_varid(bry_id,    'u_north',  u_north_id)
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1)    'u_north', bry_file(ibry)(1:lsrc)
            endif

            ierr=nf_inq_varid(bry_id,    'v_north',  v_north_id)
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1)  'v_north', bry_file(ibry)(1:lsrc)
            endif
#    ifdef NHMG
            ierr=nf_inq_varid(bry_id,    'w_north', w_north_id)
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1) 'w_north', bry_file(ibry)(1:lsrc)
            endif
#    endif
#   endif
#   ifdef T_FRC_BRY
            ierr=nf_inq_varid(bry_id,'temp_north', t_north_id(itemp))
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1) 'temp_north', bry_file(ibry)(1:lsrc)
            endif
#    ifdef SALINITY
            ierr=nf_inq_varid(bry_id,'salt_north', t_north_id(isalt))
            if (ierr /= nf_noerr) then
              if (ierr_all == nf_noerr) ierr_all=ierr
              write(*,1) 'salt_north', bry_file(ibry)(1:lsrc)
            endif
#    endif
#   endif
#  endif
          endif !<-- north_bry_active
# endif
        endif
      endif   !<-- bry_rec < 0

  1   format(1x,'### ERROR: get_bry_all :: Cannot find variable ''',
     &                                A, ''' in file ''', A, '''.')

! Read data from the file:  Check if model time is bounded by past and
!----- ---- ---- --- -----  future data times: if not, increment record
! and cycling indices, flip time index and read a new portion of data.
! Repeat until model time falls between the two data times.

! IMPORTANT! Always read "bry_time" even if no actual data needs to be
! read: this is necessary to make sure that running out of data records
! condition is the same for all MPI nodes.

      do while(bry_time(itbry) < time+0.5*dt .and. ierr == nf_noerr)
        call advance_record(bry_rec, ntbry,  bry_ncycle, bry_cycle,
     &                                       ibry, max_bry,  ierr)
        if (bry_rec < 0) then
          ierr=nf_close(bry_id) ; bry_id=-1
          ibry=ibry+1 ;  goto 100  !--> proceed to the next file
        elseif (ierr /= nf_noerr) then
          if (ierr_all == nf_noerr) ierr_all=ierr
        endif

        if (ierr == nf_noerr) then
          ierr=nf_get_var1_double(bry_id, bry_time_id, bry_rec, dtime)
          if (ierr == nf_noerr) then
            itbry=min(3-itbry,ntbry)
            bry_time(itbry)=dtime*day2sec + bry_cycle*bry_ncycle

            if (need_bry) then
# ifdef OBC_WEST
              if (west_bry_active) then
#  ifdef Z_FRC_BRY
                ierr=read_bry_EW(bry_id, zeta_west_id, bry_rec, r_var,
     &                                   zeta_west_dt(0,itbry), 1)
                if (ierr_all == nf_noerr) ierr_all=ierr
#  endif
#  ifdef M2_FRC_BRY
                ierr=read_bry_EW(bry_id, ubar_west_id, bry_rec, u_var,
     &                                   ubar_west_dt(0,itbry), 1)
                if (ierr_all == nf_noerr) ierr_all=ierr

                ierr=read_bry_EW(bry_id, vbar_west_id, bry_rec, v_var,
     &                                   vbar_west_dt(0,itbry), 1)
                if (ierr_all == nf_noerr) ierr_all=ierr
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
                ierr=read_bry_EW(bry_id, u_west_id, bry_rec, u_var,
     &                                   u_west_dt(0,1,itbry), N)
                if (ierr_all == nf_noerr) ierr_all=ierr

                ierr=read_bry_EW(bry_id, v_west_id, bry_rec, v_var,
     &                                   v_west_dt(0,1,itbry), N)
                if (ierr_all == nf_noerr) ierr_all=ierr
#    ifdef NHMG
                ierr=read_bry_EW(bry_id, w_west_id, bry_rec, r_var,
     &                                   w_west_dt(0,1,itbry), N)
                if (ierr_all == nf_noerr) ierr_all=ierr
#    endif
#   endif
#   ifdef T_FRC_BRY
                do itrc=1,NT
                  ierr=read_bry_EW(bry_id, t_west_id(itrc),  bry_rec,
     &                           r_var, t_west_dt(0,1,itbry,itrc), N)
                  if (ierr_all == nf_noerr) ierr_all=ierr
                enddo
#   endif
#  endif
              endif
# endif

# ifdef OBC_EAST
              if (east_bry_active) then
#  ifdef Z_FRC_BRY
                ierr=read_bry_EW(bry_id, zeta_east_id, bry_rec, r_var,
     &                                   zeta_east_dt(0,itbry), 1)
                if (ierr_all == nf_noerr) ierr_all=ierr
#  endif
#  ifdef M2_FRC_BRY
                ierr=read_bry_EW(bry_id, ubar_east_id, bry_rec, u_var,
     &                                   ubar_east_dt(0,itbry), 1)
                if (ierr_all == nf_noerr) ierr_all=ierr

                ierr=read_bry_EW(bry_id, vbar_east_id, bry_rec, v_var,
     &                                   vbar_east_dt(0,itbry), 1)
                if (ierr_all == nf_noerr) ierr_all=ierr
#    ifdef NHMG
                ierr=read_bry_EW(bry_id, w_east_id, bry_rec, r_var,
     &                                   w_east_dt(0,1,itbry), N)
                if (ierr_all == nf_noerr) ierr_all=ierr
#    endif
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
                ierr=read_bry_EW(bry_id, u_east_id, bry_rec, u_var,
     &                                   u_east_dt(0,1,itbry), N)
                if (ierr_all == nf_noerr) ierr_all=ierr

                ierr=read_bry_EW(bry_id, v_east_id, bry_rec, v_var,
     &                                   v_east_dt(0,1,itbry), N)
                if (ierr_all == nf_noerr) ierr_all=ierr
#    ifdef NHMG
                ierr=read_bry_EW(bry_id, w_east_id, bry_rec, r_var,
     &                                   w_east_dt(0,1,itbry), N)
                if (ierr_all == nf_noerr) ierr_all=ierr
#    endif
#   endif
#   ifdef T_FRC_BRY
                do itrc=1,NT
                  ierr=read_bry_EW(bry_id, t_east_id(itrc),  bry_rec,
     &                           r_var, t_east_dt(0,1,itbry,itrc), N)
                  if (ierr_all == nf_noerr) ierr_all=ierr
                enddo
#   endif
#  endif
              endif
# endif

# ifdef OBC_SOUTH
              if (south_bry_active) then
#  ifdef Z_FRC_BRY
                ierr=read_bry_NS(bry_id, zeta_south_id, bry_rec, r_var,
     &                                   zeta_south_dt(0,itbry), 1)
                if (ierr_all == nf_noerr) ierr_all=ierr
#  endif
#  ifdef M2_FRC_BRY
                ierr=read_bry_NS(bry_id, ubar_south_id, bry_rec, u_var,
     &                                   ubar_south_dt(0,itbry), 1)
                if (ierr_all == nf_noerr) ierr_all=ierr

                ierr=read_bry_NS(bry_id, vbar_south_id, bry_rec, v_var,
     &                                   vbar_south_dt(0,itbry), 1)
                if (ierr_all == nf_noerr) ierr_all=ierr
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
                ierr=read_bry_NS(bry_id, u_south_id,  bry_rec, u_var,
     &                                   u_south_dt(0,1,itbry), N)
                if (ierr_all == nf_noerr) ierr_all=ierr

                ierr=read_bry_NS(bry_id, v_south_id,  bry_rec, v_var,
     &                                   v_south_dt(0,1,itbry), N)
                if (ierr_all == nf_noerr) ierr_all=ierr
#    ifdef NHMG
                ierr=read_bry_NS(bry_id, w_south_id, bry_rec, r_var,
     &                                   w_south_dt(0,1,itbry), N)
                if (ierr_all == nf_noerr) ierr_all=ierr
#    endif
#   endif
#   ifdef T_FRC_BRY
                do itrc=1,NT
                  ierr=read_bry_NS(bry_id, t_south_id(itrc),  bry_rec,
     &                           r_var, t_south_dt(0,1,itbry,itrc), N)
                  if (ierr_all == nf_noerr) ierr_all=ierr
                enddo
#   endif
#  endif
              endif
# endif

# ifdef OBC_NORTH
              if (north_bry_active) then
#  ifdef Z_FRC_BRY
                ierr=read_bry_NS(bry_id, zeta_north_id, bry_rec, r_var,
     &                                   zeta_north_dt(0,itbry), 1)
                if (ierr_all == nf_noerr) ierr_all=ierr
#  endif
#  ifdef M2_FRC_BRY
                ierr=read_bry_NS(bry_id, ubar_north_id, bry_rec, u_var,
     &                                   ubar_north_dt(0,itbry), 1)
                if (ierr_all == nf_noerr) ierr_all=ierr

                ierr=read_bry_NS(bry_id, vbar_north_id, bry_rec, v_var,
     &                                   vbar_north_dt(0,itbry), 1)
                if (ierr_all == nf_noerr) ierr_all=ierr
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
                ierr=read_bry_NS(bry_id, u_north_id, bry_rec, u_var,
     &                                   u_north_dt(0,1,itbry), N)
                if (ierr_all == nf_noerr) ierr_all=ierr

                ierr=read_bry_NS(bry_id, v_north_id, bry_rec, v_var,
     &                                   v_north_dt(0,1,itbry), N)
                if (ierr_all == nf_noerr) ierr_all=ierr
#    ifdef NHMG
                ierr=read_bry_NS(bry_id, w_north_id, bry_rec, r_var,
     &                                   w_north_dt(0,1,itbry), N)
                if (ierr_all == nf_noerr) ierr_all=ierr
#    endif
#   endif
#   ifdef T_FRC_BRY
                do itrc=1,NT
                  ierr=read_bry_NS(bry_id, t_north_id(itrc),  bry_rec,
     &                           r_var, t_north_dt(0,1,itbry,itrc), N)
                  if (ierr_all == nf_noerr) ierr_all=ierr
                enddo
#   endif
#  endif
              endif
# endif
              ierr=ierr_all
              if (ierr == nf_noerr) then
# ifdef MPI_SILENT_MODE
                if (mynode == 0) then
# endif
                  write(*,'(9x,A,3x,A,F12.4,2(1x,A,I5))')
     &           'get_bry_all :: read all boundary fields',
     &           'bry_time =', dtime, 'rec =', bry_rec MYID
# ifdef MPI_SILENT_MODE
                endif
# endif
              else
                write(*,'(1x,2A,I4,1x,3A,1x,A,I4)') '### ERROR: ',
     &             'get_bry_all :: Cannot read rec =',    bry_rec,
     &             'of ''',  bry_file(ibry)(1:lsrc),  '''.'  MYID
              endif
            endif !<-- need_bry
          else
            write(*,'(8x,2A,I4)')  '### ERROR: get_bry_all :: ',
     &             'Cannot read variable ''bry_time''' MYID
          endif
        else
           write(*,'(/1x,A,I4,1x,A,I4/7x,3A/7x,2(A,G12.4)/)')
     &   '### ERROR: get_bry_all :: requested time record ', bry_rec,
     &   'exeeds the last record',   ntbry,   'available in file ''',
     &    bry_file(ibry)(1:lsrc),  '''',  'tdays = ', tdays,
     &   '  but the last available  bry_time =',
     &                                   bry_time(itbry)*sec2day
        endif
      enddo
      if (ierr_ext == nf_noerr) ierr_ext=ierr_all
# ifdef VERBOSE
      write(*,*) 'leaving get_bry_all' MYID
# endif
      end




      subroutine set_bry_all_tile (istr,iend,jstr,jend, ierr)
      implicit none
      integer istr,iend,jstr,jend, ierr,
     &        imin,imax,jmin,jmax, it1,it2
# ifdef SOLVE3D
     &                           , k, itrc
# endif
# if defined OBC_WEST  || defined OBC_EAST
     &                           , j
# endif
# if defined OBC_SOUTH || defined OBC_NORTH
     &                           , i
# endif

      logical need_bry
      real tmid, cff,cff1,cff2
# include "param.h"
# include "scalars.h"
# include  "boundary.h"

      need_bry=.false.
                                  ! NOTE: Temporal interpolation of
      imin=istr                   ! BRY forcing variables is done over
# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then      ! non-standard ranges of spatial
#  ifdef OBC_WEST
        need_bry=.true.           ! indices which is motivated by the
#  endif
        imin=istr-1               ! desire to avoid using MPI calls
      endif                       ! to exchange ghost points for BRY
# elif defined EXCH_EAST_WEST
      if (WEST_EXCHNG) then       ! arrays, while, at the same time,
        imin=istr-1               ! having only a single row of ghost
      endif                       ! points (so the standard "extended"
# endif
                                  ! bounds are not applicable here).
      imax=iend
# ifndef EW_PERIODIC
      if (EASTERN_EDGE) then
#  ifdef OBC_EAST
        need_bry=.true.
#  endif
        imax=iend+1
      endif
# elif defined EXCH_EAST_WEST
      if (EAST_EXCHNG) then
        imax=iend+1
      endif
# endif

      jmin=jstr
# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
#  ifdef OBC_SOUTH
        need_bry=.true.
#  endif
        jmin=jstr-1
      endif
# elif defined EXCH_NORTH_SOUTH
      if (SOUTH_EXCHNG) then
        jmin=jstr-1
      endif
# endif

      jmax=jend
# ifndef NS_PERIODIC
      if (NORTHERN_EDGE) then
#  ifdef OBC_NORTH
        need_bry=.true.
#  endif
        jmax=jend+1
      endif
# elif defined EXCH_NORTH_SOUTH
      if (NORTH_EXCHNG) then
        jmax=jend+1
      endif
# endif
                                            ! Always set synchro_flag
      tmid=time+0.5*dt                      ! regardless whether the
      it2=itbry ; cff1=bry_time(it2)-tmid   ! boundary is needed or
      it1=3-it2 ; cff2=tmid-bry_time(it1)   ! not: this is necessary
      if (cff1 < dt) synchro_flag=.true.    ! to make sure that running
      cff=1./(cff1+cff2)                    ! out of data records in
      cff1=cff1*cff ; cff2=cff2*cff         ! the file is detected by 
                                            ! all MPI nodes.
      if (need_bry) then
# ifdef OBC_WEST
        if (WESTERN_EDGE) then
#  ifdef Z_FRC_BRY
          do j=jmin,jmax
            zeta_west(j)=cff1*zeta_west_dt(j,it1)
     &                  +cff2*zeta_west_dt(j,it2)
          enddo
#  endif
#  ifdef M2_FRC_BRY
          do j=jmin,jmax
            ubar_west(j)=cff1*ubar_west_dt(j,it1)
     &                  +cff2*ubar_west_dt(j,it2)
            vbar_west(j)=cff1*vbar_west_dt(j,it1)
     &                  +cff2*vbar_west_dt(j,it2)
          enddo
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
          do k=1,N
            do j=jmin,jmax
              u_west(j,k)=cff1*u_west_dt(j,k,it1)
     &                   +cff2*u_west_dt(j,k,it2)
              v_west(j,k)=cff1*v_west_dt(j,k,it1)
     &                   +cff2*v_west_dt(j,k,it2)
#     ifdef NHMG
              w_west(j,k)= 0.0
#     endif
            enddo
          enddo
#   endif
#   ifdef T_FRC_BRY
          do itrc=1,NT
            do k=1,N
              do j=jmin,jmax
                t_west(j,k,itrc)=cff1*t_west_dt(j,k,it1,itrc)
     &                          +cff2*t_west_dt(j,k,it2,itrc)
              enddo
            enddo
          enddo
#   endif
#  endif
        endif
# endif

# ifdef OBC_EAST
        if (EASTERN_EDGE) then
#  ifdef Z_FRC_BRY
          do j=jmin,jmax
            zeta_east(j)=cff1*zeta_east_dt(j,it1)
     &                  +cff2*zeta_east_dt(j,it2)
          enddo
#  endif
#  ifdef M2_FRC_BRY
          do j=jmin,jmax
            ubar_east(j)=cff1*ubar_east_dt(j,it1)
     &                  +cff2*ubar_east_dt(j,it2)
            vbar_east(j)=cff1*vbar_east_dt(j,it1)
     &                  +cff2*vbar_east_dt(j,it2)
          enddo
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
          do k=1,N
            do j=jmin,jmax
              u_east(j,k)=cff1*u_east_dt(j,k,it1)
     &                   +cff2*u_east_dt(j,k,it2)
              v_east(j,k)=cff1*v_east_dt(j,k,it1)
     &                   +cff2*v_east_dt(j,k,it2)
#     ifdef NHMG
              w_east(j,k)= 0.0
#     endif
            enddo
          enddo
#   endif
#   ifdef T_FRC_BRY
          do itrc=1,NT
            do k=1,N
              do j=jmin,jmax
                t_east(j,k,itrc)=cff1*t_east_dt(j,k,it1,itrc)
     &                          +cff2*t_east_dt(j,k,it2,itrc)
              enddo
            enddo
          enddo
#   endif
#  endif
        endif
# endif

# ifdef OBC_SOUTH
        if (SOUTHERN_EDGE) then
#  ifdef Z_FRC_BRY
          do i=imin,imax
            zeta_south(i)=cff1*zeta_south_dt(i,it1)
     &                   +cff2*zeta_south_dt(i,it2)
          enddo
#  endif
#  ifdef M2_FRC_BRY
          do i=imin,imax
            ubar_south(i)=cff1*ubar_south_dt(i,it1)
     &                   +cff2*ubar_south_dt(i,it2)
            vbar_south(i)=cff1*vbar_south_dt(i,it1)
     &                   +cff2*vbar_south_dt(i,it2)
          enddo
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
          do k=1,N
            do i=imin,imax
              u_south(i,k)=cff1*u_south_dt(i,k,it1)
     &                    +cff2*u_south_dt(i,k,it2)
              v_south(i,k)=cff1*v_south_dt(i,k,it1)
     &                    +cff2*v_south_dt(i,k,it2)
#     ifdef NHMG
              w_south(i,k)= 0.0
#     endif
            enddo
          enddo
#   endif
#   ifdef T_FRC_BRY
          do itrc=1,NT
            do k=1,N
              do i=imin,imax
                t_south(i,k,itrc)=cff1*t_south_dt(i,k,it1,itrc)
     &                           +cff2*t_south_dt(i,k,it2,itrc)
              enddo
            enddo
          enddo
#   endif
#  endif
        endif
# endif

# ifdef OBC_NORTH
        if (NORTHERN_EDGE) then
#  ifdef Z_FRC_BRY
          do i=imin,imax
            zeta_north(i)=cff1*zeta_north_dt(i,it1)
     &                   +cff2*zeta_north_dt(i,it2)
          enddo
#  endif
#  ifdef M2_FRC_BRY
          do i=imin,imax
            ubar_north(i)=cff1*ubar_north_dt(i,it1)
     &                   +cff2*ubar_north_dt(i,it2)
            vbar_north(i)=cff1*vbar_north_dt(i,it1)
     &                   +cff2*vbar_north_dt(i,it2)
          enddo
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
          do k=1,N
            do i=imin,imax
              u_north(i,k)=cff1*u_north_dt(i,k,it1)
     &                    +cff2*u_north_dt(i,k,it2)
              v_north(i,k)=cff1*v_north_dt(i,k,it1)
     &                    +cff2*v_north_dt(i,k,it2)
#     ifdef NHMG
              w_north(i,k)= 0.0
#     endif
            enddo
          enddo
#   endif
#   ifdef T_FRC_BRY
          do itrc=1,NT
            do k=1,N
              do i=imin,imax
                t_north(i,k,itrc)=cff1*t_north_dt(i,k,it1,itrc)
     &                           +cff2*t_north_dt(i,k,it2,itrc)
              enddo
            enddo
          enddo
#   endif
#  endif
        endif
# endif
        if (cff1 < 0. .or. cff2 < 0.) then
           write(*,'(/3A/3(1x,A,F16.10)/)')            '### WARNING: ',
     &      'set_bry_all_tile :: current model time is out of bounds ',
     &      'of ''bry_time''.', 'bry_tstart =', bry_time(it1)*sec2day,
     &      'tdays =',  tdays,  'bry_tend =',   bry_time(it2)*sec2day
          ierr=ierr+1
        endif
      endif  ! <-- need_bry
      end
#else
      subroutine get_bry_empty
      end
#endif

