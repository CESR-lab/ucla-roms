#include "cppdefs.opt"

      subroutine get_init(req_rec, tindx)    ! Reads initial conditions
                                             ! from netCDF file.

      use param
      use tracers, only: t_vname, t, iTandS
      use bgc ! get_init_bgc_diag
#ifdef MARBL
      use marbl_driver, only: marbldrv_read_ss_vars_from_rst
#endif
      use coupling
      use grid
      use hidden_mpi_vars
      use analytical
      use mixing
      use ocean_vars
      use basic_output
      use scalars
      use netcdf
      use bgc_ecosys_vars
      use diagnostics
      use sponge_tune
      use nc_read_write, only: handle_ierr, ncread
      use roms_read_write
      use dimensions

      implicit none

      integer :: req_rec, tindx, max_rec, record, ncid, varid, i, j, k,
     &           ierr,  start(4), count(2), ibuff(iaux),
     &           checkdims , itrc
      real    :: time_scale
      integer :: init_type, tile=0 ! set tile to zero for compute_tile_bounds.h below. Needed to fixed closed boundaries.
      integer, parameter :: init_run=1,  rst_run=2
#ifdef EXACT_RESTART
     &                    , apprx_rst=3, exact_rst=4
      real time_bak
# ifdef NHMG
!     integer i,j,k
# endif
#endif
#define time illegal
#define tdays illegal
#define nrrec illegal

! needed for istr, iend, etc in order to fix closed boundaries to u(istr)=0, etc
# include "compute_tile_bounds.h"

c--#define VERBOSE
#ifdef VERBOSE
      write(*,'(3(2x,A,I3))') 'enter get_init: req_rec =', req_rec,
     &            'tindx =', tindx, 'may_day_flag =', may_day_flag
#endif
      if (may_day_flag /= 0) return      !-->  EXIT

! Open conditions netCDF file for reading.  Check that all spatial
! dimensions in that file are consistent with the model, determine how
! many time records are available in the file and determine the number
! of the record from which the data will be read. The record is set as
! follows:
!         (1) if there is only one time record available in the file,
!             then use that record REGARDLESS of the value of "nrrec"
!             supplied in the input parameter file;
!
!         (2) if the file has multiple records and
!
!             (2a) nrrec > 0 then record number "nrrec" is read,
!                  provided that "nrrec" is within the available
!                  records (error message is printed otherwise);
!
!             (2b) nrrec < 0, THE LAST available record is used.

      init_type=0
      ierr=nf90_open(ininame, nf90_nowrite, ncid)
      if (ierr == nf90_noerr) then
        ierr=checkdims (ncid, ininame, max_rec)
        if (ierr == nf90_noerr) then
          if (max_rec > 0) then
            if (req_rec > 0) then
              if (req_rec <= max_rec) then
                record=req_rec
              else
                write(*,'(/1x,2A,I4,1x,A/12x,A,I4,1x,3A/)')
     &              '### ERROR: get_init :: requested restart time ',
     &              'record',req_rec, 'exceeds number', 'of records',
     &                       max_rec,  'available in netCDF file ''',
     &                                       trim(ininame), '''.'
              endif
            else
              record=max_rec
            endif
          else
            record=1
          endif
        endif
      else
        write(*,'(/1x,4A/12x,A/)')  '### ERROR: get_init :: Cannot ',
     &                'open netCDF file ''', trim(ininame), '''.',
     &                                             nf90_strerror(ierr)
      endif
      if (ierr /= nf90_noerr) goto 99                     !--> ERROR

! Read in evolving model variables:
!----- -- -------- ----- ----------
! Time: find netCDF id, read value, read attribute 'units' and set
! starting time index and time clock in days.  Note that time units
! read below also saved as vname(3,indxTime) and thereafter used to 
! control output time units literally copying it from the initial
! condition to restart/history/averages output files and writing time
! is seconds or days accordingly.
!
! Note that if EXACT_RESTART CPP-switch is defined, make a "soft"
! attempt to do exact restart, where "soft" means that exact restart
! is done only when file of initial conditions contains sufficient
! data, i.e. two consecutive time records of evolving fields one time
! step apart from each other.  Thus, in order to accept the file for
! exact restart, it must pass two consecutive checks: (i) ocean_time
! values in two consecutive values must differ by "dt" of current run,
! and (ii) "ntstart" from two consecutive records of netCDF structure
! "time_step" must differ by one.  If either check fails, forward step
! is used as the initial time step.  "get_init" is expected to be
! called twice consecutively as
!
!           call get_init(req_rec=rec-1, tindx=2)
!           call get_init(req_rec=rec,   tindx=1)
!
! where "rec" is record number in netCDF file which contains fields
! corresponding to time step "n" while "rec-1" corresponds to "n-1"
! (hence, making it possible to start time stepping with regular LF
! predictor step rather than forward), both checks are performed
! during the first call, tindx=2. If either check fails, the exact
! restart is cancelled and no reading of 2D and 3D fields will be
! performed for tindx=2.
!
! The possibility of exact restart is communicated with the rest of
! the code via integer variable "forw_start" which is set exclussively
! by this routine and is used as part of CPP-macro
!
!             FIRST_TIME_STEP iic==forw_start
!
! where the possibilities are as follows:
!
!       forw_start=1   means that "exact" restart, i.e., forward
!                      step is to be performed only during absolutely
!                      first time step, iic=1, and
!
!       forw_start=ntstart, means that restarted time stepping should
!                           also begin with forward step (approximate
!                           restart).
!
! This mechanism of exact restart is designed to handle essentially
! three situations: (1) initial run using a 3rd-party file which
! contain initial time, but does not contain "time_step". In this
! case ntstart is set to 1, and forward step is assumed at the first
! step; (2) restart from restart file generated by this code, but
! with deactivated CPP-switch EXACT_RESTART.  This file contains both
! both time variable and "time_step" structure, but only one
! consecutive record. This situation is identified automatically and
! approximate restart is assumed. This is compatibility mode. This
! also incldes restart from a history file generated by this code.
! (3) restart from a file created by this code with activated
! EXACT_RESTART.

      ierr=nf90_inq_varid(ncid, vname(1,indxTime), varid)

! The following is done for backward compatibility: normally time
! variable is named "ocean_time", but legacy startup files may name
! it either "roms_time" or "scrum_time".

      if (ierr /= nf90_noerr) then
        ierr=nf90_inq_varid(ncid, 'roms_time', varid)
      endif
      if (ierr /= nf90_noerr) then
        ierr=nf90_inq_varid(ncid, 'scrum_time', varid)
      endif
      if (ierr/=0) call handle_ierr(ierr,'get_init: time var not found')

      ierr=nf90_get_var(ncid, varid, start_time, start=(/record/) )
      if (ierr/=0) call handle_ierr(ierr,'get_init: get time')
      ierr=nf90_get_att(ncid, varid, 'units', vname(3,indxTime))
      if (ierr/=0) call handle_ierr(ierr,'get_init: get time units')

      if (vname(3,indxTime)(1:6) == 'second') then
        time_scale=1.D0
      elseif (vname(3,indxTime)(1:3) == 'day') then
        time_scale=day2sec
      else
        write(*,'(/1x,4A/12x,3A/)')    '### ERROR: get_init :: ',
     &  'unknown units for variable ''',trim(vname(1,indxTime)),
     &         '''', 'in netCDF file ''', trim(ininame), '''.'
        goto 99                                    !--> ERROR
      endif
      start_time=start_time*time_scale
!     start_time=start_time*time_scale + 288

#ifdef EXACT_RESTART
      if (tindx == 2) then
        forw_start=0
        if (record < max_rec) then
          time_bak=start_time
          ierr=nf90_get_var(ncid, varid, start_time, start=(/record+1/))
          if (ierr == nf90_noerr) then
            start_time=start_time*time_scale
# ifdef VERBOSE
            write(*,'(3(1x,A,F16.6))') 'time_bak =', time_bak,
     &               'start_time =', start_time, 'dt =', dt
# endif

! Note that expression "abs(start_time-time_bak-dt) < 0.001*dt" below
! is a roundoff-error tolerant version of "start_time == time_bak+dt".

            if (abs(start_time-time_bak-dt) < 0.01*dt) then
               forw_start=1
            else
              mpi_nonexit_warn write(*,'(1x,2A,2I4/10x,4A/10x,A/)')
     &                'WARNING: Exact restart is requested, but ',
     &                'is not possible: records', record,record+1,
     &                'in ''', trim(ininame),  ''' are not ',
     &                'consecutive time steps ==> proceeding ',
     &                          'with forward initial time step.'
            endif
          else
            write(*,2) trim(vname(1,indxTime)), record,
     &                 trim(ininame), nf90_strerror(ierr)
            goto 99
          endif
        else
           mpi_nonexit_warn write(*,'(1x,2A/10x,4A)')
     &           'WARNING: Exact restart is requested, but is not ',
     &           'possible: initial',  'file ''', trim(ininame),
     &                    ''' does not contain sufficient records.'
        endif
        if (forw_start /= 1) return
        forw_start=0
      endif
#endif
      if (ierr /= nf90_noerr)  goto 99                    !--> ERROR

! Check whether variable 'time_step' is present, which can be a
! structure of four to up to eight integer numbers storing time step
! number and the corresponding record numbers for output files.
! If present, use them to restart the time step number and record
! counters (i.e., technically  this is "restart" as opposite to
! "initial run");  otherwise initialise all the counters to zeroes.

      ierr=nf90_inq_varid(ncid, 'time_step', varid)
      if (ierr == nf90_noerr) then
        ierr=nf90_inquire_variable(ncid, varid, dimids=ibuff)
        if (ierr/=0) 
     &    call handle_ierr(ierr,'get_init: cannot get time_step dimids')
        ierr=nf90_inquire_dimension(ncid, ibuff(1), len=count(1))
        if (ierr/=0) 
     &    call handle_ierr(ierr,'get_init: cannot get time dim size')
        start(1)=1 ; start(2)=record ; count(2)=1
        ibuff(1:iaux)=0
        ierr=nf90_get_var(ncid, varid, ibuff, start, count)
        if (ierr/=0) 
     &  call handle_ierr(ierr,'get_init: cannot read time_step')
        ntstart=ibuff(1)+1

#ifdef EXACT_RESTART
        if (tindx == 2 .and. record < max_rec) then
          start(2)=record+1
          ierr=nf90_get_var(ncid, varid, ibuff, start, count)
          if (ierr/=0) 
     &  call handle_ierr(ierr,'get_init: cannot read 2nd time_step')
# ifdef VERBOSE
          write(*,*) 'ibuff(1),ntstart =', ibuff(1), ntstart
# endif
          if (ibuff(1) == ntstart) then
            forw_start=1
          else
          mpi_nonexit_warn write(*,'(1x,2A,2I4/10x,4A/10x,A)')
     &             'WARNING: Exact restart is requested, but is not ',
     &             'possible: records',  record,   record+1,  'in ''',
     &              trim(ininame),  ''' are not consecutive time ',
     &                                    'steps ==> proceeding with',
     &                                   'forward initial time step.'
            return   !--> no need to read preliminary record
          endif
        elseif (tindx == 1) then
          if (forw_start == 1) then
            init_type=exact_rst
          else
            init_type=apprx_rst
          endif
        endif
#else
        init_type=rst_run
#endif
        if (ierr /= nf90_noerr)  goto 99                  !--> ERROR
      else
        init_type=init_run
        ntstart=1               ! netCDF variable "time_step" not
!        nrecrst=0               ! found: proceed with initializing
!        nrechis=0               ! all counters to zero (initial run).

        ierr=0                                                       ! DevinD added to not trigger break later
      endif
#ifdef EXACT_RESTART
      if (tindx == 1 .and. forw_start == 0) forw_start=ntstart
#endif
#ifdef VERBOSE
      write(*,'(1x,2A,F12.4,1x,A,I4)')   'get_init: reading initial ',
     &   'fields for time =', start_time*sec2day, 'record =', record
#endif

! Read initial fields:
!---------------------

! Free-surface and barotropic 2D momentuma, XI- and ETA-components

      start=1; start(3)=record                                       ! 2D vars

      call ncread(ncid,'zeta',zeta(i0:i1,j0:j1,1),start)
      zeta(i0:i1,j0:j1,1)=zeta(i0:i1,j0:j1,1)*rmask(i0:i1,j0:j1)
      call exchange_xxx(zeta(:,:,1))                              

      call ncread(ncid,'ubar',ubar( 1:i1,j0:j1,1),start)
      ubar(1:i1,j0:j1,1) = ubar(1:i1,j0:j1,1)*umask(1:i1,j0:j1)  
      call exchange_xxx(ubar(:,:,1))

      call ncread(ncid,'vbar',vbar(i0:i1, 1:j1,1),start)
      vbar(i0:i1,1:j1,1) = vbar(i0:i1,1:j1,1)*vmask(i0:i1,1:j1)
      call exchange_xxx(vbar(:,:,1))

! Two sets of fast-time-averaged barotropic fluxes needed for exact
! restart in the case when using Adams-Bashforth-like extrapolation of
! vertically-integrated 3D velocities for computing momentum advection
! and Coriolis terms for 3D --> 2D forcing of barotropic mode.  Once
! again, adopting "soft policy" with respect to their presense/absence
! in the file: if not found use forward step instead of exact restart.

#ifdef SOLVE3D
# ifdef EXACT_RESTART
#  ifdef EXTRAP_BAR_FLUXES
      ierr=nf90_inq_varid(ncid, 'DU_avg2', varid)  ! change these to varid
      if (ierr == nf90_noerr) then
        ierr=nf90_inq_varid(ncid, 'DV_avg2', varid)
        if (ierr == nf90_noerr) then
          ierr=nf90_inq_varid(ncid, 'DU_avg_bak', varid)
          if (ierr == nf90_noerr) then
            ierr=nf90_inq_varid(ncid, 'DV_avg_bak', varid)
            if (ierr == nf90_noerr) then
              call ncread(ncid,'DU_avg2',      DU_avg2( 1:i1,j0:j1),start)
              call ncread(ncid,'DV_avg2',      DV_avg2(i0:i1, 1:j1),start)
              call ncread(ncid,'DU_avg_bak',DU_avg_bak( 1:i1,j0:j1),start)
              call ncread(ncid,'DV_avg_bak',DV_avg_bak(i0:i1, 1:j1),start)
              if (ierr /= nf90_noerr) goto 99                  !--> ERROR
            else
              forw_start=ntstart      !--> cancel exact restart
            endif
          else
            forw_start=ntstart      !--> cancel exact restart
          endif
        else
          forw_start=ntstart    !--> cancel exact restart
        endif
      else
        forw_start=ntstart    !--> cancel exact restart
      endif

#  elif defined PRED_COUPLED_MODE
      ierr=nf90_inq_varid(ncid, 'rufrc_bak', varid)
      if (ierr == nf90_noerr) then
        ierr=nf90_inq_varid(ncid,  'rvfrc_bak', varid)
        if (ierr == nf90_noerr) then

          call ncread(ncid,'rufrc_bak', rufrc_bak( 1:i1,j0:j1,tindx),start)    ! DevinD this update was untested
          call ncread(ncid,'rvfrc_bak', rvfrc_bak(i0:i1, 1:j1,tindx),start)
          if (ierr /= nf90_noerr) goto 99                     !--> ERROR

        else
          forw_start=ntstart    !--> cancel exact restart
        endif
      else
        forw_start=ntstart    !--> cancel exact restart
      endif
#  endif
# endif /*EXACT_RESTART*/

! 3D momentum components in XI- and ETA-directions

      start=1; start(4)=record                                       ! 3D vars

      call ncread(ncid,vname(1,indxU), u(1:i1,j0:j1,:,tindx),start) 
      do k=1,nz
        u(1:i1,j0:j1,k,tindx) = u(1:i1,j0:j1,k,tindx)*umask(1:i1,j0:j1)
      enddo
      call exchange_xxx(u(:,:,:,tindx) )

      call ncread(ncid,vname(1,indxV), v(i0:i1, 1:j1,:,tindx),start)           ! v requires masking and
       do k=1,nz
         v(i0:i1,1:j1,k,tindx) = v(i0:i1,1:j1,k,tindx)*vmask(i0:i1,1:j1)
       enddo
      call exchange_xxx(v(:,:,:,tindx) )

# if defined NHMG
      ! since ncdf_read already does an exchange to fill buffers, we're
      ! ready to initialy w using no-div
!     do j=1,Lm
!       do k=1,N
!         do i=1,Mm
!           DC(i,0)=pm(i,j)*pn(i,j)
!           w(i,j,k,nnew)=w(i,j,k-1,nnew)
!    &      - DC(i,0)*(
!    &       +0.5*(Hz(i,j,k)+Hz(i+1,j,k))*dn_u(i+1,j)*u(i+1,j,k,nnew)
!    &       -0.5*(Hz(i-1,j,k)+Hz(i,j,k))*dn_u(i,j)*u(i,j,k,nnew)
!    &       +0.5*(Hz(i,j,k)+Hz(i,j+1,k))*dm_v(i,j+1)*v(i,j+1,k,nnew)
!    &       -0.5*(Hz(i,j-1,k)+Hz(i,j,k))*dm_v(i,j)*v(i,j,k,nnew) )
!         enddo
!       enddo
!     enddo
# endif


! Tracer variables.

      do itrc=1,nt

        ierr=nf90_inq_varid(ncid, t_vname(itrc), varid)
        if (ierr == nf90_noerr) then
          call ncread(ncid, t_vname(itrc), t(i0:i1,j0:j1,:,tindx,itrc), start)
          do k=1,N
            t(i0:i1,j0:j1,k,tindx,itrc)=t(i0:i1,j0:j1,k,tindx,itrc)*rmask(i0:i1,j0:j1)
          enddo
          call exchange_xxx(t(:,:,:,tindx,itrc) ) 
        else

          if (itrc <= iTandS) then ! temperature and salt(s) always require inital condition
            write(*,1) trim(t_vname(itrc)), trim(ininame),
     &                  nf90_strerror(ierr)
            goto 99                                         !--> ERROR
          else
            t(:,:,:,tindx,itrc) = 0.0
            if(mynode==0) write(*,*) ' --- WARNING: '
     &         , trim(t_vname(itrc))
     &         , ' not in initial file.  Initialized to 0.0'
            ierr=nf90_noerr
          endif

        endif
      enddo

      start=1; start(3)=record                                       ! 2D vars
# ifdef LMD_KPP
      ierr=nf90_inq_varid(ncid, vname(1,indxHbls), varid)
      if (ierr == nf90_noerr) then
        call ncread(ncid,vname(1,indxHbls),hbls(i0:i1,j0:j1),start)      ! mask doesn't change answer but makes
        hbls(i0:i1,j0:j1)=hbls(i0:i1,j0:j1)*rmask(i0:i1,j0:j1)     ! output mask 0, rather than random numbers
        call exchange_xxx(hbls(:,:) )
      else
        mpi_nonexit_warn write(*,'(1x,6A)')        'WARNING: netCDF ',
     & 'variable ''', trim(vname(1,indxHbls)), ''' not found in ''',
     &           trim(ininame), ''' ==> initialized to zero state.'

      endif
# endif
# ifdef LMD_BKPP
      ierr=nf90_inq_varid(ncid, vname(1,indxHbbl), varid)
      if (ierr == nf90_noerr) then
        call ncread(ncid,vname(1,indxHbbl),hbbl(i0:i1,j0:j1),start)      ! mask doesn't change answer but makes
        hbbl(i0:i1,j0:j1)=hbbl(i0:i1,j0:j1)*rmask(i0:i1,j0:j1)     ! output mask 0, rather than random numbers
        call exchange_xxx(hbbl(:,:) )            ! confirmed needs an exchange, not sure why
      else
        mpi_nonexit_warn write(*,'(1x,6A)')        'WARNING: netCDF ',
     & 'variable ''', trim(vname(1,indxHbbl)), ''' not found in ''',
     &      trim(ininame), ''' ==> initialized to zero state.'

      endif
# endif

#ifdef MARBL
!     Read in MARBL saved_state variables here
      call marbldrv_read_ss_vars_from_rst(ncid,record)
#endif /* MARBL */      
      
#endif /* SOLVE3D */


#ifdef DIAGNOSTICS
      if (diag_pflx) call get_init_slow(ncid,record,tindx)
      if (ub_tune)   call get_init_ub(ncid,record)
#endif      

! Close input NetCDF file and  write greeting message depending
! on the the type of initial/restart procedure performed above.

      ierr=nf90_close(ncid)
#ifdef MPI_SILENT_MODE
      if (mynode == 0) then
#endif
       if (tindx == 1) then
         if (init_type == init_run) then
          write(*,'(6x,2A,F12.4,1x,A,I4)') 'get_init :: Read initial ',
     &   'conditions for day =', start_time*sec2day, 'record =',record
#ifdef EXACT_RESTART
         elseif (init_type == exact_rst) then
           write(*,'(6x,A,F12.4,1x,A,I4,A)')
     &    'get_init :: Exact restart from day =',  start_time*sec2day,
     &    'rec =', record, '.'
         elseif (init_type == apprx_rst) then
           write(*,'(6x,A,F12.4,1x,A,I4,A)')
     &      'get_init :: Approximate, single-step restart from day =',
     &       start_time*sec2day,   'rec =', record, '.'
#else
         elseif (init_type == rst_run) then
           write(*,'(6x,A,F12.4,1x,A,I4,A)')
     &          'get_init: Restarted from day =', start_time*sec2day,
     &     'rec =', record, '.'
#endif
         else
           write(*,'(/1x,A/)') '### ERROR: get_init :: Unknown Error.'
           goto 99
         endif
       endif  !<-- tindex==1
#ifdef MPI_SILENT_MODE
      endif
#endif

      ! Here we catch an initial file that does not have the boundary values of u/v
      ! set to zero, even though we prescribed a closed boundary. This is necessary
      ! for conservation tests.
#if !defined OBC_WEST && !defined EW_PERIODIC
      if(WESTERN_EDGE) then
        if(mynode==0) write (*,'(6x,A)')        ! S&W bry always contained in start node (==0)
     &  'get_init :: set closed west  boundary row to zero (u/ub/v/vb)'
        u(   istr  ,:     ,:,tindx)=0.0         ! ensure closed boundary u on boundary is zero
        ubar(istr  ,:     ,1)      =0.0         ! ubar read into index 1 above
      endif
#endif
#if !defined OBC_EAST && !defined EW_PERIODIC
      if(EASTERN_EDGE) then
        if(mynode==NNODES-1) write (*,'(6x,A)') ! N&E bry always contained in final node (==NNODES-1)
     &  'get_init :: set closed east  boundary row to zero (u/ub/v/vb)'
        u(   iend+1,:     ,:,tindx)=0.0
        ubar(iend+1,:     ,1)      =0.0
      endif
#endif
#if !defined OBC_SOUTH && !defined NS_PERIODIC
      if(SOUTHERN_EDGE) then
        if(mynode==0) write (*,'(6x,A)')
     &  'get_init :: set closed south boundary row to zero (u/ub/v/vb)'
        v(   :     ,jstr  ,:,tindx)=0.0
        vbar(:     ,jstr  ,1)      =0.0
      endif
#endif
#if !defined OBC_NORTH && !defined NS_PERIODIC
      if(NORTHERN_EDGE) then
        if(mynode==NNODES-1) write (*,'(6x,A)')
     &  'get_init :: set closed north boundary row to zero (u/ub/v/vb)'
        v(   :     ,jend+1,:,tindx)=0.0
        vbar(:     ,jend+1,1)      =0.0
      endif
#endif

#if defined ANA_BRY
      call ana_init_generic(1,nx,1,ny)
#endif


#ifdef VERBOSE
      write(*,'(1x,3(1x,A,I6))') 'return from get_init, ntstart =',
     &                                                  ntstart
# ifdef EXACT_RESTART
     &                              ,   'forw_start =', forw_start
# endif
#endif
      return
  99  may_day_flag=2
      return
  1   format(/1x,'### ERROR: get_init :: Cannot find variable ''', A,
     &                 '''' /12x, 'in netCDF file ''', A, ''''/12x,A)
  2   format(/1x,'### ERROR: get_init :: Cannot read variable ''', A,
     &        '''' /12x, 'rec =', I5, 1x, 'file = ''', A, ''''/12x,A)
      end
