#include "cppdefs.opt"
#ifdef AVERAGES

      subroutine wrt_avg             ! Write time-averaged fields
                                     ! into averages netCDF file.

      use param
      use wec_frc    ! , only: wrt_avg_wec
      use bgc        ! wrt_avg_bgc_diag
      use tracers    ! , only: wrt_avg_trc
      use read_write ! nc_write_var, rp_var, etc
      use netcdf     ! nf90 subroutines
      use ocean_vars
      use eos_vars
      use mixing
      use ncvars
      use scalars

      implicit none

      integer record, i, ierr, lstr, lvar, lenstr, ncdf_write, var_id_tmp
# ifdef SOLVE3D
     &             , itrc
# endif
      character(len=18) tstring

# include "netcdf.inc"
# if defined MPI && !defined PARALLEL_FILES
      integer status(MPI_STATUS_SIZE), blank
      if (mynode > 0) then
        call MPI_Recv (blank, 1, MPI_INTEGER, mynode-1,  3,
     &                       ocean_grid_comm, status, ierr)
      endif
# endif

! Create/open averages file; write grid arrays, if so needed,
! check for errors, then set record within the file "record".

      call def_avg (nrecavg, ierr)
      lstr=lenstr(avgname)
      if (ierr /= nf90_noerr) goto 99         !!! WARNING: Here it is
      nrecavg=max(nrecavg,1)                !!! assumed that global
      if (nrpfavg == 0) then                !!! averages record index
        record=nrecavg                      !!! "nrecavg" is already
      else                                  !!! advanced by main.
        record=1+mod(nrecavg-1, nrpfavg)
      endif

! Write out time-averaged variables:
! ----- --- ------------- ----------

      ! write ocean time
      ! ----------------
      ierr=nf90_inq_varid(ncavg, 'ocean_time', var_id_tmp)
      ierr=nf90_put_var(ncavg, var_id_tmp, time_avg,(/record/))

! Barotropic mode variables: free-surface and 2D momentum components
! in XI-,ETA-directions.

      if (wrtavg(indxZ) .and. ierr == nf90_noerr) then
        call nc_write_var(ncavg, zeta_avg, 1, vname(1,indxZ), rp_var,
     &                    record, ierr)
      endif

      if (wrtavg(indxUb) .and. ierr == nf90_noerr) then
        call nc_write_var(ncavg, ubar_avg, 1, vname(1,indxUb), up_var,
     &                    record, ierr)
      endif

      if (wrtavg(indxVb) .and. ierr == nf90_noerr) then
        call nc_write_var(ncavg, vbar_avg, 1, vname(1,indxVb), vp_var,
     &                    record, ierr)
      endif
# ifdef SOLVE3D

! 3D momentum components in XI- and ETA-directions.

      if (wrtavg(indxU) .and. ierr == nf90_noerr) then
        call nc_write_var(ncavg, u_avg, N, vname(1,indxU), up_var,
     &                    record, ierr)
      endif

      if (wrtavg(indxV) .and. ierr == nf90_noerr) then
        call nc_write_var(ncavg, v_avg, N, vname(1,indxV), vp_var,
     &                    record, ierr)
      endif

! Tracer variables.

      call wrt_avg_trc(record, ierr)

! Density anomaly.

      if (wrtavg(indxR) .and. ierr == nf90_noerr) then
        call nc_write_var(ncavg, rho_avg, N, vname(1,indxR), rp_var,
     &                    record, ierr)
      endif

! S-coordinate omega vertical velocity (m/s).

      if (wrtavg(indxO) .and. ierr == nf90_noerr) then
        call nc_write_var(ncavg, w_avg, N+1, vname(1,indxO), rp_var,
     &                    record, ierr)
      endif

! True vertical velocity (m/s).

      if (wrtavg(indxW) .and. ierr == nf90_noerr) then
        call nc_write_var(ncavg, wvl_avg, N, vname(1,indxW), rp_var,
     &                    record, ierr)
      endif

! Vertical eddy viscosity coefficient for momentum.

      if (wrtavg(indxAkv) .and. ierr == nf90_noerr) then
        call nc_write_var(ncavg, akv_avg(START_2D_ARRAY,0), N+1,
     &                    vname(1,indxAkv), rp_var, record, ierr)
      endif

! Vertical heat conductivity

      if (wrtavg(indxAkt) .and. ierr == nf90_noerr) then
        call nc_write_var(ncavg, akt_avg(START_2D_ARRAY,0), N+1,
     &                    vname(1,indxAkt), rp_var, record, ierr)
      endif
#  ifdef SALINITY

! Vertical diffusion coefficient for salinity.

      if (wrtavg(indxAks) .and. ierr == nf90_noerr) then
        call nc_write_var(ncavg, Aks_avg(START_2D_ARRAY,0), N+1,
     &                    vname(1,indxAks), rp_var, record, ierr)
      endif
#  endif
#  ifdef LMD_KPP

! Depth of planetary boundary layer [m].

      if (wrtavg(indxHbls) .and. ierr == nf90_noerr) then
        call nc_write_var(ncavg, hbl_avg, 1, vname(1,indxHbls), rp_var,
     &                    record, ierr)
      endif
#  endif
#  ifdef LMD_BKPP

! Thickness of KPP bottom boundary layer [m].

      if (wrtavg(indxHbbl) .and. ierr == nf90_noerr) then
        call nc_write_var(ncavg, hbbl_avg, 1, vname(1,indxHbbl), rp_var,
     &                    record, ierr)
      endif
#  endif
# endif /* SOLVE3D */

#ifdef WEC
      call wrt_avg_wec(record, ierr, may_day_flag)
#endif

#if defined BIOLOGY_BEC2 && defined BEC2_DIAG
      call wrt_avg_bgc_diag(record, ierr)
#endif

  1   format(/1x, '### ERROR: wrt_avg :: Cannot write variable ''',
     &             A, ''' into averages file, rec =', i6, 3x,A,i4)
      goto 100
  99  if (may_day_flag == 0) may_day_flag=3
 100  continue

! Close or synchronize netCDF file to disk to allow other
! processes to access data immediately after it is written.

# if defined MPI && !defined PARALLEL_FILES
      ierr=nf90_close(ncavg)
      if (nrpfavg>0 .and. record>=nrpfavg) ncavg=-1
# else
      if (nrpfavg>0 .and. record>=nrpfavg) then
        ierr=nf90_close(ncavg)
        ncavg=-1
      else
        ierr=nf90_sync(ncavg)
      endif
# endif
      if (ierr == nf90_noerr) then
# ifdef MPI_SILENT_MODE
        if (mynode==0) then
# endif
          write(tstring,'(F18.8)') time_avg*sec2day
          i=1
          do while (i<18 .and. tstring(i:i)==' ')
            i=i+1
          enddo
          write(*,'(7x,2A,1x,A,1x,A,I4,A,I4,1x,A,I3)') 'wrt_avg :: ',
     &      'wrote averages for time_avg =', tstring(i:i+8), 'rec =',
     &       record, '/', nrecavg  MYID
# ifdef MPI_SILENT_MODE
        endif
# endif
      else
        write(*,'(/1x,2A/)')      '### ERROR: wrt_avg :: Cannot ',
     &                         'synchronize/close averages file.'
        if (may_day_flag==0) may_day_flag=3
      endif
# if defined MPI && !defined PARALLEL_FILES
      if (mynode < NNODES-1) call MPI_Send (blank, 1, MPI_INTEGER,
     &                           mynode+1, 3, ocean_grid_comm, ierr)
# endif
      end
#else
      subroutine wrt_avg_empty
      end
#endif /* AVERAGES */
