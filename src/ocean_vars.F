
      module ocean_vars
      ! ocean variables and output routines for those variables

#include "cppdefs.opt"

      use param
      use netcdf
      use nc_read_write
      use roms_read_write
      implicit none
      private

      ! module preamble:  ![
      ! 2D - taken from what was ocean2d.F
      real,public,allocatable,dimension(:,:,:) :: zeta            ! free surface elevation [m] and barotropic
      real,public,allocatable,dimension(:,:,:) :: ubar            ! velocity components in XI-directions
      real,public,allocatable,dimension(:,:,:) :: vbar
      real,public,allocatable,dimension(:,:) :: zeta_avg
      real,public,allocatable,dimension(:,:) :: ubar_avg
      real,public,allocatable,dimension(:,:) :: vbar_avg
#ifdef NHMG
      real,public,allocatable,dimension(:,:) :: nh_ubar
      real,public,allocatable,dimension(:,:) :: nh_vbar
      real,public,allocatable,dimension(:,:) :: nh_wcor
#endif

      ! 3D - taken from what was ocean2d.F
#ifdef SOLVE3D
      real,public,allocatable,dimension(:,:,:,:) :: u
      real,public,allocatable,dimension(:,:,:,:) :: v
      real,public,allocatable,dimension(:,:,:)   :: u_avg
      real,public,allocatable,dimension(:,:,:)   :: v_avg
# if defined NHMG
      real,public,allocatable,dimension(:,:,:,:) :: w
      real,public,allocatable,dimension(:,:,:,:) :: nhdu
      real,public,allocatable,dimension(:,:,:,:) :: nhdv
      real,public,allocatable,dimension(:,:,:,:) :: nhdw
# endif

      real,public,allocatable,dimension(:,:,:) :: FlxU
      real,public,allocatable,dimension(:,:,:) :: FlxV
      real,public,allocatable,dimension(:,:,:) :: We      ! explicit
      real,public,allocatable,dimension(:,:,:) :: Wi      ! implicit

      real,public,allocatable,dimension(:,:,:) :: w_avg
      real,public,allocatable,dimension(:,:,:) :: wvl_avg

      real,public,allocatable,dimension(:,:,:) :: Hz      ! height of rho-cell
      real,public,allocatable,dimension(:,:,:) :: z_r     ! depth at rho-points
      real,public,allocatable,dimension(:,:,:) :: z_w     ! depth at   w-points
# if defined NHMG || defined NONTRAD_COR
      real,public,allocatable,dimension(:,:,:) :: dzdxi
      real,public,allocatable,dimension(:,:,:) :: dzdeta
# endif
#endif  /* SOLVE3D */

      integer,public :: rechis                            ! record in history file, needed by bgc and wec

      public :: init_arrays_ocean2d
      public :: init_arrays_ocean3d
      public :: wrt_his, wrt_rst

      contains  !]

!----------------------------------------------------------------------
      subroutine init_arrays_ocean2d  ![
      use ncvars
      use scalars

      implicit none

      allocate( zeta(GLOBAL_2D_ARRAY,4) ); zeta=0.         ! zeta(:,:,knew) needs to be =0. for set_depth_tile
      allocate( ubar(GLOBAL_2D_ARRAY,4) ); ubar=0.         ! since knew can change if exact restart or not
      allocate( vbar(GLOBAL_2D_ARRAY,4) ); vbar=0.         ! set all zeta = 0.

      if (wrtavg(indxZ)) then
        allocate( zeta_avg(GLOBAL_2D_ARRAY) )
        zeta_avg=0.                                        ! avg needs to be 0. because Nan x 0. = Nan in set_avg.F
      endif
      if (wrtavg(indxUb)) then
        allocate( ubar_avg(GLOBAL_2D_ARRAY) )
        ubar_avg=0.
      endif
      if (wrtavg(indxVb)) then
        allocate( vbar_avg(GLOBAL_2D_ARRAY) )
        vbar_avg=0.
      endif

#ifdef NHMG
      allocate( nh_ubar(GLOBAL_2D_ARRAY) )
      allocate( nh_vbar(GLOBAL_2D_ARRAY) )
      allocate( nh_wcor(GLOBAL_2D_ARRAY) )
#endif

      end subroutine init_arrays_ocean2d  !]

!----------------------------------------------------------------------
      subroutine init_arrays_ocean3d  ![
      use ncvars
      use scalars
      implicit none

      ! DevinD - again lots of extra first touched variables that weren't in init_arrays originally?

#ifdef SOLVE3D
      allocate( u(GLOBAL_2D_ARRAY,N,3) )
      u(:,:,:,1)=init
      u(:,:,:,2)=0.0                             ! index 2 used on rhs u(indx) in pre_step for DC.
      u(:,:,:,3)=init                            ! multiplied by zero but can't be a nan.
      allocate( v(GLOBAL_2D_ARRAY,N,3) )
      v(:,:,:,1)=init
      v(:,:,:,2)=0.0                             ! index 2 used on rhs v(indx) in pre_step for DC.
      v(:,:,:,3)=init                            ! multiplied by zero but can't be a nan.

      if (wrtavg(indxU)) then
        allocate( u_avg(GLOBAL_2D_ARRAY,N) )
        u_avg=0.
      endif
      if (wrtavg(indxV)) then
        allocate( v_avg(GLOBAL_2D_ARRAY,N) )
        v_avg=0.
      endif

# if defined NHMG
      allocate( w(GLOBAL_2D_ARRAY,0:N,3) )
      allocate( nhdu(GLOBAL_2D_ARRAY,1:N,2) )
      allocate( nhdv(GLOBAL_2D_ARRAY,1:N,2) )
      allocate( nhdw(GLOBAL_2D_ARRAY,0:N,2) )
# endif

      allocate( FlxU(GLOBAL_2D_ARRAY,N) ) ; FlxU=init
      allocate( FlxV(GLOBAL_2D_ARRAY,N) ) ; FlxV=init
      allocate( We(GLOBAL_2D_ARRAY,0:N) ) ; We=init        ! explicit
      allocate( Wi(GLOBAL_2D_ARRAY,0:N) ) ; Wi=init        ! implicit

      if (wrtavg(indxO)) then
        allocate( w_avg(GLOBAL_2D_ARRAY,0:N) )
        w_avg=0
      endif
      if (wrtavg(indxW)) then
        allocate( wvl_avg(GLOBAL_2D_ARRAY,0:N) )
        wvl_avg=0.
      endif


      allocate( Hz(GLOBAL_2D_ARRAY,N) )    ; Hz=init       ! height of rho-cell
      allocate( z_r(GLOBAL_2D_ARRAY,N) )   ; z_r=init      ! depth at rho-points
      allocate( z_w(GLOBAL_2D_ARRAY,0:N) ) ; z_w=init      ! depth at   w-points
# if defined NHMG || defined NONTRAD_COR
      allocate( dzdxi(GLOBAL_2D_ARRAY,1:N)  )
      allocate( dzdeta(GLOBAL_2D_ARRAY,1:N) )
# endif
#endif  /* SOLVE3D */

      end subroutine init_arrays_ocean3d  !]

!----------------------------------------------------------------------
      ! def_his calls itself to get def_avg. So need to fix that first before adding here

!----------------------------------------------------------------------
      subroutine def_rst(total_rec, ierr)  ![
      ! create or openrestart netCDF file.

      use tracers,         only: t_vname, t_lname, t_units, wrt_t
      use roms_read_write, only: output_root_name
      use ncvars
      use scalars
      use grid                                                       ! def_grid

      implicit none

      ! input
      integer :: total_rec, ierr
      ! local
      logical :: create_new_file
      integer :: lenstr, checkdims, my_nf_def_dim, varid
      integer :: rec, lrst, lvar, timedim, prev_fill_mode, itrc
#ifdef MASK_LAND_DATA
# ifdef DBLEPREC
      real(kind=8),parameter :: spv_set=1.D+33
# else
      real(kind=4),parameter :: spv_set=1.E+33
# endif
#endif

! Put time record index into file name. In  the case when model
! output is to be arranged into sequence of named files, the naming
! convention is as follows: 'rst_root.INDEX.[MPI_node.]nc', where
! INDEX is an integer number such that (i) it is divisible by the
! specified number of records per file; and (ii)
!
!      INDEX + record_within_the_file = total_record
!
! where, 1 =< record_within_the_file =< records_per_file, so that
! total_record changes continuously throughout the sequence of files.

      ierr=0
      rstname=trim(output_root_name) / / '_rst.nc'
      lrst=lenstr(rstname)
#if defined MPI && defined PARALLEL_FILES
        ! Insert MPI node numbers to file name
        call insert_node(rstname, lrst, mynode, NNODES, ierr)  ! previously done in ncdf_create_file:
#endif
      lrst=lenstr(rstname)
      if (nrpfrst > 0) then
        lvar=total_rec - (1+mod(total_rec-1, nrpfrst))
        call insert_time_index (rstname, lrst, lvar, ierr)
        if (ierr /= 0) goto 99
      endif

! Decide whether to create a new file, or open existing one.
! Overall the whole code below is organized into 3-way switch,
!
! 10  if (create_new_file) then
!        .... create new file, save netCDF ids for all variables;
!     elseif (ncrst==-1) then
!        .... try to open existing file and check its dimensions
!       if (cannot be opened or rejected) then
!         create_new_file=.true.
!         goto 10
!       endif   and prepare
!        .... prepare the file for adding new data,
!        .... find and save netCDF ids for all variables
!     else
!        .... just reopen, no checking, all ids are assumed
!        .... to be already known (this branch is only for
!        .... output into single file from MPI code).
!     endif
!
! which is designed to implement a flexible opening policy: if
!  (i) ldefhis=.true., it forces creation of a new file [if the file
!      already exists, it will be overwritten]; on the other hand, if
! (ii) ldefhis=.false., it is assumed that the file already exists
!      and an attempt to open it is made first; if the attempt is
!      successful, the file is prepared for appending hew data;
!      if it fails, a new file is created.

      create_new_file=ldefhis
      if (ncrst /= -1) create_new_file=.false.
#if defined MPI && !defined PARALLEL_FILES
      if (mynode > 0) create_new_file=.false.
#endif

  10  if (create_new_file) then
        ierr=nf90_create(rstname(1:lrst), nf90_netcdf4, ncrst)
        if (ierr == nf90_noerr) then
          ierr=nf90_set_fill(ncrst, nf90_nofill, prev_fill_mode)
          if (ierr /= nf90_noerr) then
            write(*,'(1x,4A,I4)') '### ERROR: def_rst :: Cannot ',
     &        'switch to ''nf_nofill'' mode.', nf90_strerror(ierr)
     &         MYID
            goto 99
          endif
        else
          write(*,'(/1x,4A/12x,A/)')  '### ERROR: def_rst :: Cannot ',
     &               'create restart NetCDF file ''', rstname(1:lrst),
     &                                       '''',  nf90_strerror(ierr)
          goto 99              !--> ERROR
        endif

        call put_global_atts (ncrst, ierr)

#ifdef PUT_GRID_INTO_RESTART
        if (total_rec == 1) call def_grid(ncrst)
#endif

! Define evolving model variables:
! ------ -------- ----- ----------
! Time step number and time record numbers:
        varid = nccreate(ncrst,'time_step',(/'auxil','time'/),(/iaux,0/),nf90_int)
        ierr=nf90_put_att (ncrst, varid, 'long_name',
     &             'time step and record numbers from initialization')

! Time.
        varid = nccreate(ncrst,vname(1,indxTime),(/'time'/),(/0/),nf90_double)
        ierr=nf90_put_att(ncrst, varid, 'long_name', vname(2,indxTime))
        ierr=nf90_put_att(ncrst, varid, 'units', vname(3,indxTime))

! Free-surface.
        varid = nccreate(ncrst,vname(1,indxZ),(/'xi_rho','eta_rho','time'/),(/xi_rho,eta_rho,0/), NF_FTYPE)
        ierr=nf90_put_att (ncrst, varid, 'long_name', vname(2,indxZ))
        ierr=nf90_put_att (ncrst, varid, 'units', vname(3,indxZ))
        ierr = nf90_put_att(ncrst,varid,'_FillValue',spv_set)

! 2D momenta in XI- and ETA-components
        varid = nccreate(ncrst,vname(1,indxUb),(/'xi_u','eta_rho','time'/),(/xi_u,eta_rho,0/), NF_FTYPE)
        ierr=nf90_put_att(ncrst, varid, 'long_name',
     &                                   vname(2,indxUb))
        ierr=nf90_put_att(ncrst, varid, 'units', vname(3,indxUb))
        ierr = nf90_put_att(ncrst,varid,'_FillValue',spv_set)
        varid = nccreate(ncrst,vname(1,indxVb),(/'xi_rho','eta_v','time'/),(/xi_rho,eta_v,0/), NF_FTYPE)
        ierr=nf90_put_att(ncrst, varid, 'long_name',
     &                                   vname(2,indxVb))
        ierr=nf90_put_att(ncrst, varid, 'units', vname(3,indxVb))
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)

#ifdef SOLVE3D
# ifdef EXACT_RESTART
#  ifdef EXTRAP_BAR_FLUXES
        varid = nccreate(ncrst,'DU_avg2',(/'xi_u','eta_rho','time'/),(/xi_u,eta_rho,0/), NF_FTYPE)
        ierr=nf90_put_att(ncrst, varid, 'long_name',
     &           '<<fast-time averaged ubar(:,:,n+1/2)>>')
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)

        varid = nccreate(ncrst,'DV_avg2',(/'xi_rho','eta_v','time'/),(/xi_rho,eta_v,0/), NF_FTYPE)
        ierr=nf90_put_att(ncrst, varid, 'long_name',
     &           '<<fast-time-averaged vbar(:,:,n+1/2)>>')
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)

        varid = nccreate(ncrst,'DU_avg_bak',(/'xi_u','eta_rho','time'/),(/xi_u,eta_rho,0/), NF_FTYPE)
        ierr=nf90_put_att(ncrst, varid, 'long_name',
     &       '<back-step mixed fast-time-averaged ubar(:,:,n-1)>')
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)

        varid = nccreate(ncrst,'DV_avg_bak',(/'xi_rho','eta_v','time'/),(/xi_rho,eta_v,0/), NF_FTYPE)
        ierr=nf90_put_att(ncrst, varid, 'long_name',
     &       '<back-step mixed fast-time-averaged vbar(:,:,n-1)>')
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)
#  elif defined PRED_COUPLED_MODE
        varid = nccreate(ncrst,'rufrc_bak',(/'xi_u','eta_rho','time'/),(/xi_u,eta_rho,0/), NF_FTYPE)
        ierr=nf90_put_att(ncrst, varid, 'long_name',
     &                        '3D-to-2D forcing, XI-component')
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)
        varid = nccreate(ncrst,'rvfrc_bak',(/'xi_rho','eta_v','time'/),(/xi_rho,eta_v,0/), NF_FTYPE)
        ierr=nf90_put_att(ncrst, varid, 'long_name',
     &                       '3D-to-2D forcing, ETA-component')
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)
#  endif
# endif  /* EXACT_RESTART */

! 3D momenta in XI- and ETA-directions.
        varid = nccreate(ncrst, vname(1,indxU),(/'xi_u','eta_rho','s_rho','time'/),(/xi_u,eta_rho,N,0/), NF_FTYPE)
        ierr=nf90_put_att (ncrst, varid, 'long_name',
     &                                     vname(2,indxU))
        ierr=nf90_put_att (ncrst, varid, 'units',
     &                                     vname(3,indxU))
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)

        varid = nccreate(ncrst,vname(1,indxV),(/'xi_rho','eta_v','s_rho','time'/),(/xi_rho,eta_v,N,0/), NF_FTYPE)
        ierr=nf90_put_att (ncrst, varid, 'long_name',
     &                                     vname(2,indxV))
        ierr=nf90_put_att (ncrst, varid, 'units',
     &                                     vname(3,indxV))
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)

! Tracer variables.
        do itrc=1,NT
          varid = nccreate(ncrst,t_vname(itrc),(/'xi_rho','eta_rho','s_rho','time'/),
     &                                         (/xi_rho,eta_rho,N,0/), NF_FTYPE)
          ierr=nf90_put_att (ncrst, varid, 'long_name',
     &                       t_lname(itrc))
          ierr=nf90_put_att (ncrst, varid, 'units', t_units(itrc))
          ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)
        enddo

# ifdef LMD_KPP
! Depth of surface planetary boundary layer.
!        ierr=nf_def_var (ncrst, vname(1,indxHbls)(1:lvar),
!     &                            NF_FTYPE, 3, r2dgrd, rstHbls)
        varid = nccreate(ncrst,vname(1,indxHbls),(/'xi_rho','eta_rho','time'/),(/xi_rho,eta_rho,0/), NF_FTYPE)
        ierr=nf90_put_att (ncrst, varid, 'long_name',
     &                                 vname(2,indxHbls))
        ierr=nf90_put_att (ncrst, varid, 'units',
     &                                 vname(3,indxHbls))
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)
# endif
# ifdef LMD_BKPP
! Thickness of bottom boundary layer.
!        ierr=nf_def_var (ncrst, vname(1,indxHbbl),
!     &                            NF_FTYPE, 3, r2dgrd, rstHbbl)
        varid = nccreate(ncrst,vname(1,indxHbbl),(/'xi_rho','eta_rho','time'/),(/xi_rho,eta_rho,0/), NF_FTYPE)
        ierr=nf90_put_att (ncrst, varid, 'long_name',
     &                                 vname(2,indxHbbl))
        ierr=nf90_put_att (ncrst, varid, 'units',
     &                                 vname(3,indxHbbl))
        ierr = nf90_put_att(ncrst,varid,'_FillValue', spv_set)
# endif
#endif /* SOLVE3D */

        ierr=nf90_enddef(ncrst)
#ifdef MPI_SILENT_MODE
        if (mynode==0) then
#endif
          write(*,'(7x,4A,1x,A,i4)') 'def_rst :: Created new ',
     &          'netCDF file ''', rstname(1:lrst), '''.' MYID
#ifdef MPI_SILENT_MODE
        endif
#endif

! Open an existing file and prepare for appending data.
! ==== == ======== ==== === ======= === ========= =====
! Check consistency of the dimensions of fields from the file with
! model dimensions. Determine the current size of unlimited dimension
! and set initial record [in the case of MPI serialized output, at
! this moment the last time record is assumed to be **partially**
! written by MPI processes with lower rank. Thus the next write is
! expected to be into the same record rather than next one (except
! MPI-master, who initializes the record).
!
! In the case when file is rejected (whether it cannot be opened, or
! something is wrong with its dimensions, create new file.

      elseif (ncrst == -1) then
        ierr=nf90_open(rstname, nf90_write, ncrst)
        if (ierr == nf90_noerr) then
          ierr=nf90_set_fill(ncrst, nf90_nofill, prev_fill_mode)
          if (ierr /= nf90_noerr) then
            write(*,'(1x,4A,I4)') '### ERROR: def_rst :: Cannot ',
     &        'switch to ''nf_nofill'' mode.', nf90_strerror(ierr)
     &         MYID
            goto 99
          endif

          ierr=checkdims(ncrst, rstname, rec)
          if (ierr == nf90_noerr) then
            if (nrpfrst==0) then
              ierr=rec+1 - nrecrst
            else
              ierr=rec+1 - (1+mod(nrecrst-1, abs(nrpfrst)))
            endif
            if (ierr > 0) then
              mpi_master_only write(*,
     &           '(/1x,A,I5,1x,3A/21x,A,I5,1x,A,I5,1x,A/21x,A/)')
     &           'WARNING: def_rst :: The actual number of records',
     &            rec, 'in netCDF file ''',  rstname(1:lrst), '''',
     &           'exceeds record', rec+1-ierr, '/', total_rec,
     &           'for restart initial conditions.',
     &           'All records beyond this will be overwriding.'
              rec=rec-ierr
            elseif (nrpfrst==0) then
              total_rec=rec+1           ! <-- set to the next record
#if defined MPI && !defined PARALLEL_FILES
              if (mynode > 0) total_rec=total_rec-1
#endif
            endif
            ierr=nf90_noerr
          endif
        endif

        if (ierr /= nf90_noerr) then
#if defined MPI && !defined PARALLEL_FILES
          if (mynode==0) then
            create_new_file=.true.
            goto 10
          else
            write(*,'(/1x,4A, 1x,A,I4/)')  '### ERROR: def_rst :: ',
     &        'Cannot open restart netCDF file ''', rstname(1:lrst),
     &        '''.' MYID
            goto 99                                     !--> ERROR
          endif
#else
          create_new_file=.true.
          goto 10
#endif
        endif

! Find netCDF variable IDs of evolving model variables:
! ---- ------ -------- --- -- -------- ----- ----------

! Time step indices and time:

        ierr=nf90_inq_varid(ncrst, 'time_step', rstTstep)
        if (ierr == nf90_noerr) then
          lvar=lenstr(vname(1,indxTime))
          ierr=nf90_inq_varid(ncrst, vname(1,indxTime)(1:lvar), rstTime)
          if (ierr /= nf90_noerr) then
            write(*,1) vname(1,indxTime)(1:lvar), rstname(1:lrst)
          endif
        else
          write(*,1) 'time_step', rstname(1:lrst)
        endif
        if (ierr /= nf90_noerr) goto 99                     !--> ERROR

! Free-surface and 2D momenta in XI- and ETA-directions

        lvar=lenstr(vname(1,indxZ))
        ierr=nf90_inq_varid(ncrst, vname(1,indxZ)(1:lvar), rstZ)
        if (ierr == nf90_noerr) then
          lvar=lenstr(vname(1,indxUb))
          ierr=nf90_inq_varid(ncrst, vname(1,indxUb)(1:lvar), rstUb)
          if (ierr == nf90_noerr) then
            lvar=lenstr(vname(1,indxVb))
            ierr=nf90_inq_varid(ncrst, vname(1,indxVb)(1:lvar), rstVb)
            if (ierr /= nf90_noerr) then
              write(*,1) vname(1,indxVb)(1:lvar), rstname(1:lrst)
            endif
          else
            write(*,1) vname(1,indxUb)(1:lvar), rstname(1:lrst)
          endif
        else
          write(*,1) vname(1,indxZ)(1:lvar),  rstname(1:lrst)
        endif
        if (ierr /= nf90_noerr) goto 99                     !--> ERROR

#ifdef SOLVE3D
# ifdef EXACT_RESTART
#  ifdef EXTRAP_BAR_FLUXES
        ierr=nf90_inq_varid(ncrst,   'DU_avg2',   rst_DU_avg2)
        if (ierr == nf90_noerr) then
          ierr=nf90_inq_varid(ncrst,  'DV_avg2',     rst_DV_avg2)
          if (ierr == nf90_noerr) then
            ierr=nf90_inq_varid(ncrst, 'DU_avg_bak', rst_DU_avg_bak)
            if (ierr == nf90_noerr) then
              ierr=nf90_inq_varid(ncrst, 'DV_avg_bak', rst_DV_avg_bak)
              if (ierr /= nf90_noerr) then
                write(*,1) 'DV_avg_bak', rstname(1:lrst)
              endif
            else
              write(*,1) 'DU_avg_bak', rstname(1:lrst)
            endif
          else
            write(*,1) 'DV_avg2', rstname(1:lrst)
          endif
        else
          write(*,1) 'DU_avg2', rstname(1:lrst)
        endif
#  elif defined PRED_COUPLED_MODE
        ierr=nf90_inq_varid(ncrst, 'rufrc_bak', rst_rufrc)
        if (ierr == nf90_noerr) then
          ierr=nf90_inq_varid(ncrst, 'rvfrc_bak', rst_rvfrc)
          if (ierr /= nf90_noerr) then
            write(*,1) 'rvfrc_bak', rstname(1:lrst)
          endif
        else
          write(*,1) 'rufrc_bak', rstname(1:lrst)
        endif
#  endif
        if (ierr /= nf90_noerr) goto 99                     !--> ERROR
# endif

! 3D momenta n XI- and ETA-directions.

        lvar=lenstr(vname(1,indxU))
        ierr=nf90_inq_varid (ncrst, vname(1,indxU)(1:lvar), rstU)
        if (ierr == nf90_noerr) then
          lvar=lenstr(vname(1,indxV))
          ierr=nf90_inq_varid (ncrst, vname(1,indxV)(1:lvar), rstV)
          if (ierr /= nf90_noerr) then
            write(*,1) vname(1,indxV)(1:lvar), rstname(1:lrst)
          endif
        else
          write(*,1) vname(1,indxU)(1:lvar), rstname(1:lrst)
        endif
        if (ierr /= nf90_noerr) goto 99                     !--> ERROR

! Tracer variables.

        do itrc=1,NT
          lvar=lenstr(vname(1,indxT+itrc-1))
          ierr=nf90_inq_varid(ncrst, vname(1,indxT+itrc-1)(1:lvar),
     &                                                  rstT(itrc))
          if (ierr /= nf90_noerr) then
            write(*,1) vname(1,indxT+itrc-1)(1:lvar), rstname(1:lrst)
            goto 99                                       !--> ERROR
          endif
        enddo

# ifdef LMD_KPP
! Depth of surface planetary boundary layer.

        lvar=lenstr(vname(1,indxHbls))
        ierr=nf90_inq_varid(ncrst, vname(1,indxHbls)(1:lvar), rstHbls)
        if (ierr /= nf90_noerr) then
          write(*,1) vname(1,indxHbls)(1:lvar), rstname(1:lrst)
          goto 99                                        !--> ERROR
        endif
# endif
# ifdef LMD_BKPP
! Thickness of bottom boundary layer.
        ierr=nf90_inq_varid(ncrst, vname(1,indxHbbl), rstHbbl)
        if (ierr /= nf90_noerr) then
        write(*,1) vname(1,indxHbbl), rstname
          goto 99                                       !--> ERROR
        endif
# endif
#endif
#ifdef MPI_SILENT_MODE
        if (mynode == 0) then
#endif
          write(*,'(7x,4A,I4,2x,A,I4)') 'def_rst :: open existing ',
     &            'file ''', rstname(1:lrst), ''' from rec =', rec
     &             MYID
#ifdef MPI_SILENT_MODE
        endif
#endif

#if defined MPI && !defined PARALLEL_FILES
      else
        ierr=nf90_open(rstname(1:lrst), nf_write, ncrst)
        if (ierr == nf90_noerr) then
          ierr=nf_set_fill(ncrst, nf_nofill, prev_fill_mode)
          if (ierr == nf90_noerr) then
# ifndef MPI_SILENT_MODE
            write(*,'(7x,4A,I4)') 'def_rst :: reopen existing file ''',
     &              rstname(1:lrst), '''.' MYID
# endif
          else
            write(*,'(1x,4A,I4)') '### ERROR: def_rst :: Cannot ',
     &        'switch to ''nf_nofill'' mode.', nf90_strerror(ierr)
     &         MYID
            goto 99
          endif
        else
          write(*,'(/1x,4A,1x,A,I4/)') '### ERROR: def_rst :: Cannot ',
     &            'reopen file ''', rstname(1:lrst),'''.'
     &             MYID
          goto 99                                         !--> ERROR
        endif
#endif
      endif              !<-- create_new_file

#ifdef PUT_GRID_INTO_RESTART

! Write grid variables.
! ----- ---- ----------

      if (total_rec == 1) call wrt_grid(ncrst, rstname, lrst)
#endif
  99  return                                              !--> ERROR
   1  format(/1x,'### ERROR: def_rst :: Cannot find variable ''', A,
     &                            ''' in netCDF file ''', A, '''.'/)

      end subroutine def_rst  !]

!----------------------------------------------------------------------
      subroutine wrt_his  ![
      ! writes requested model fields at requested levels into history netCDF file.
      !
      ! NOTE: wrt_his is called in middle of roms_step in main.F because
      !       even though u/ubar, etc, were computed at the end of the
      !       previous timestep, other variables like 'rho', 'Akv', etc
      !       have only just been computed.

      use param
      use tracers    !, only: wrt_his_trc
      use eos_vars
      use grid
      use mixing
      use ncvars
      use scalars
      use work_mod
      use mpi
      implicit none

      integer :: ierr, i, lstr, lvar, lenstr, var_id_tmp
     &        , start1d(2), count(2), ibuff(iaux), ncdf_write
      integer :: start(4)
      character(len=18) :: tstring
#ifdef SOLVE3D
      integer :: tile, itrc, j,k
#endif
#if defined MPI && !defined PARALLEL_FILES
      integer :: status(MPI_STATUS_SIZE), blank
      if (mynode > 0) call MPI_Recv (blank, 1, MPI_INTEGER,
     &            mynode-1, 2, ocean_grid_comm, status, ierr)
#endif

      call def_his(nrechis, ierr)         ! create/open history and set rechis index within the file.
      if (ierr /= nf90_noerr) goto 99
      lstr=lenstr(hisname)
      nrechis=max(nrechis,1)              !!! WARNING: Once time
      if (nrpfhis == 0) then              !!! stepping has been
        rechis=nrechis                    !!! started, it is assumed
      else                                !!! that the global history
        rechis=1+mod(nrechis-1, nrpfhis)  !!! rechis index "nrechis"
      endif                               !!! is advanced by main.

      ! write out evolving model variables:
      ibuff(1)=iic-1  ; ibuff(2)=nrecrst    ! Time step and
      ibuff(4:iaux)=0 ; ibuff(3)=nrechis    ! rechis numbers.
#ifdef AVERAGES
      ibuff(4)=nrecavg
#endif
      start1d(1)=1      ; count(1)=iaux
      start1d(2)=rechis ; count(2)=1
      start=1

      ierr=nf90_inq_varid(nchis, 'time_step', var_id_tmp)
      ierr=nf90_put_var(nchis, var_id_tmp, ibuff, start1d, count)
      if (ierr /= nf90_noerr) then
        write(*,'(/1x,3A,i6/11x,A,3x,A,i4/)') '### ERROR: wrt_his :: ',
     &       'Cannot write variable ''time_step'' into history file, ',
     &       'rec =', rechis, nf90_strerror(ierr) MYID
        goto 99                                           !--> ERROR
      endif


! Time: Note that (1) time is always double precision; and (2) which
! one to write into the file, time or tdays, is decided by attribute
! "units" stored as vname(3,indxTime).  When computng time is always
! in seconds, however it is more convenient to have it in days in all
! the files.

      ierr=nf90_inq_varid(nchis, 'ocean_time', var_id_tmp)
      ierr=nf90_put_var(nchis, var_id_tmp, time,(/rechis/))

! Barotropic mode variables: free-surface and 2D momentum
! components in XI-,ETA-directions.
! 'knew' index is normally 'n+1' timestep, but at this point
! its value hasn't evolved from 'n+1' of previous timestep,
! which makes it 'n' of this timestep still.

      start(3)=rechis
      if (wrthis(indxZ))  call ncwrite(nchis,vname(1,indxZ), zeta(i0:i1,j0:j1,knew),start)
      if (wrthis(indxUb)) call ncwrite(nchis,vname(1,indxUb),ubar( 1:i1,j0:j1,knew),start)
      if (wrthis(indxVb)) call ncwrite(nchis,vname(1,indxVb),vbar(i0:i1, 1:j1,knew),start)

#ifdef SOLVE3D
! 3D momentum components in XI- and ETA-directions:
! 'nstp' index is current timestep 'n', which was computed as the final 'nnew'
! in the previous timestep (same result).
! wrt_his called at the middle of next timestep as some variables only calculated there for t=n.

      start(3)=1; start(4)=rechis
      if (wrthis(indxU)) call ncwrite(nchis,vname(1,indxU),u( 1:i1,j0:j1,:,nstp),start)
      if (wrthis(indxV)) call ncwrite(nchis,vname(1,indxV),v(i0:i1, 1:j1,:,nstp),start)

      call wrt_his_trc(start)                              ! tracer variables

      if (wrthis(indxR)) then
# ifdef SPLIT_EOS
        call ncwrite(nchis,vname(1,indxR),rho1(i0:i1,j0:j1,:),start)
# else
        call ncwrite(nchis,vname(1,indxR), rho(i0:i1,j0:j1,:),start)
# endif
      endif
      if (wrthis(indxO)) then                              ! s-coordinate omega vertical velocity (m/s).
        do k=0,N
          do j=0,Mm+1
            do i=0,Lm+1
              work(i,j,k)=pm(i,j)*pn(i,j)*(We(i,j,k)+Wi(i,j,k))
            enddo
          enddo
        enddo
        call ncwrite(nchis,vname(1,indxO),work(i0:i1,j0:j1,:),start)
      endif
      if (wrthis(indxW)) then                              ! true vertical velocity (m/s).
        do tile=0,NSUB_X*NSUB_E-1
          call wvlcty (tile, work)
        enddo
        ! need 1:N as work is 1:N+1 in size
        ! put as 0:N-1 because z-dimension array is declared in correctly in wvlcty_tile.
        call ncwrite(nchis,vname(1,indxW),work(i0:i1,j0:j1,0:N-1),start)
      endif
      if (wrthis(indxAkv)) call ncwrite(nchis,vname(1,indxAkv),Akv(i0:i1,j0:j1,:),start)
      if (wrthis(indxAkt)) call ncwrite(nchis,vname(1,indxAkt),Akt(i0:i1,j0:j1,:,itemp),start)
# ifdef SALINITY
      if (wrthis(indxAks)) call ncwrite(nchis,vname(1,indxAks),Akt(i0:i1,j0:j1,:,isalt),start)
# endif
      start(3)=rechis                                      ! back to 2D vars
# ifdef LMD_KPP
      if (wrthis(indxHbls)) call ncwrite(nchis,vname(1,indxHbls),hbls(i0:i1,j0:j1,nstp),start)
# endif
# ifdef LMD_BKPP
      if (wrthis(indxHbbl)) call ncwrite(nchis,vname(1,indxHbbl),hbbl(i0:i1,j0:j1,nstp),start)
# endif
#endif
      start(3)=1; start(4)=rechis                          ! back to 3D vars
!#if defined BIOLOGY_BEC2 && defined BEC2_DIAG
!      call wrt_his_bgc_diag(rechis, ierr)
!#endif

  1   format(/1x, '### ERROR: wrt_his :: Cannot write variable ''',
     &              A, ''' into history file, rec =', i6, 3x,A,i4)
      goto 100
  99  if (may_day_flag == 0) may_day_flag=3
 100  continue

! Close or synchronize netCDF file to the disk to allow other
! processes to access data immediately after it is written.

#if defined MPI && !defined PARALLEL_FILES
      ierr=nf90_close (nchis)
      if (nrpfhis > 0 .and. rechis >= nrpfhis) nchis=-1
#else
      if (nrpfhis > 0 .and. rechis >= nrpfhis) then
        ierr=nf90_close (nchis)
        nchis=-1
      else
        ierr=nf90_sync(nchis)
      endif
#endif
      if (ierr == nf90_noerr) then
#ifdef MPI_SILENT_MODE
        if (mynode == 0) then
#endif
          write(tstring,'(F18.8)') tdays
          i=1
          do while (i < 18 .and. tstring(i:i) == ' ')
            i=i+1
          enddo
          write(*,'(7x,A,1x,A,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
     &      'wrt_his :: wrote history, tdays =', tstring(i:i+8),
     &      'step =', ibuff(1),  'rec =', rechis, '/',  nrechis
     &       MYID
#ifdef MPI_SILENT_MODE
        endif
#endif
      else
        write(*,'(/1x,2A/)')      '### ERROR: wrt_his :: Cannot ',
     &                       'synchronize or close history file.'
        if (may_day_flag == 0) may_day_flag=3
      endif
#if defined MPI && !defined PARALLEL_FILES
      if (mynode < NNODES-1) call MPI_Send (blank, 1, MPI_INTEGER,
     &                           mynode+1, 2, ocean_grid_comm, ierr)
#endif
      end subroutine wrt_his  !]

!----------------------------------------------------------------------
      subroutine wrt_rst  ![
      ! write model prognostic variables into restart netCDF file.

      use tracers, only: wrt_rst_trc
      use coupling
      use mixing
      use ncvars
      use scalars
      use mpi
      use read_write                                       ! temporary addtion until writing fixed

      implicit none

      integer :: record, start(2), count(2), ibuff(iaux),
     &           i, ierr, lstr, lvar, lenstr, ncdf_write, varid, itrc
      character(len=18) :: tstring

#if defined MPI && !defined PARALLEL_FILES
      integer :: status(MPI_STATUS_SIZE), blank
      if (mynode > 0) call MPI_Recv (blank, 1, MPI_INTEGER,
     &            mynode-1, 1, ocean_grid_comm, status, ierr)
#endif

! Create/open restart and set record index within the file.

      call def_rst(nrecrst, ierr)
      if (ierr /= nf90_noerr) goto 99
      lstr=lenstr(rstname)
      nrecrst=max(nrecrst,1)                  !!! WARNING: Here it
      if (nrpfrst == 0) then                  !!! is assumed that
        record=nrecrst                        !!! global restart
      else                                    !!! record index
        record=1+mod(nrecrst-1, abs(nrpfrst)) !!! nrecrst is already
      endif                                   !!! advanced by main.

! Write out evolving model variables:
! ----- --- -------- ----- ----------

      ibuff(1)=iic-1  ; ibuff(2)=nrecrst        ! Time step and
      ibuff(4:iaux)=0 ; ibuff(3)=nrechis        ! record numbers.
#ifdef AVERAGES
      ibuff(4)=nrecavg
#endif
      start(1)=1      ; count(1)=iaux
      start(2)=record ; count(2)=1

      ierr=nf90_inq_varid(ncrst,'time_step',varid)
      ierr=nf90_put_var(ncrst, varid, ibuff, start, count)
      if (ierr /= nf90_noerr) then
        write(*,'(/1x,3A,i6/11x,A,3x,A,i4/)') '### ERROR: wrt_rst :: ',
     &       'Cannot write variable ''time_step'' into restart file, ',
     &       'rec =', record, nf90_strerror(ierr) MYID
        goto 99                                           !--> ERROR
      endif

! Time: Note that (1) time is always double precision; and (2) which
! one to write into the file, time or tdays, is decided by attribute
! "units" stored as vname(3,indxTime).  When computng time is always
! in seconds, however it is more convenient to have it in days in all
! the files; (3) this code is designed for politeness - it reads time
! units from the initial condition file (see "get_init.F") and uses
! it throughout the entire run.

      ierr=nf90_inq_varid(ncrst, 'ocean_time', varid)
      if (vname(3,indxTime)(1:6) == 'second') then
        ierr=nf90_put_var(ncrst, varid, time,(/record/))
      else
        ierr=nf90_put_var(ncrst, varid, tdays,(/record/))
      endif
      if (ierr /= nf90_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(*,'(/1x,4A,i6/11x,A,3x,A,i4/)') '### ERROR: wrt_rst :: ',
     &        'Cannot write variable ''',    vname(1,indxTime)(1:lvar),
     &        ''' into restart file, rec =', record, nf90_strerror(ierr)
     &            MYID
        goto 99                                           !--> ERROR
      endif

! Free-surface and 2D momentum components in XI- and ETA-directions.
      call nc_write_var(ncrst, zeta(START_2D_ARRAY,knew), 1,
     &                    vname(1,indxZ), rp_var, record, ierr)
      if (ierr == nf90_noerr) then
        call nc_write_var(ncrst, ubar(START_2D_ARRAY,knew), 1,
     &                    vname(1,indxUb), up_var, record, ierr)
        if (ierr == nf90_noerr) then
          call nc_write_var(ncrst, vbar(START_2D_ARRAY,knew), 1,
     &                      vname(1,indxVb), vp_var, record, ierr)
          if (ierr /= nf90_noerr) then
            lvar=lenstr(vname(1,indxVb))
            write(*,1) vname(1,indxVb)(1:lvar), record MYID
          endif
        else
          lvar=lenstr(vname(1,indxUb))
          write(*,1) vname(1,indxUb)(1:lvar), record MYID
        endif
      else
        lvar=lenstr(vname(1,indxZ))
        write(*,1) vname(1,indxZ)(1:lvar), record MYID
      endif
      if (ierr /= nf90_noerr) goto 99                !--> ERROR

#ifdef SOLVE3D
# ifdef EXACT_RESTART
#  ifdef EXTRAP_BAR_FLUXES
      call nc_write_var(ncrst, DU_avg2, 1, 'DU_avg2', up_var, record, ierr)
      if (ierr == nf90_noerr) then
        call nc_write_var(ncrst, DV_avg2, 1, 'DV_avg2', vp_var, record, ierr)
        if (ierr == nf90_noerr) then
          call nc_write_var(ncrst, DU_avg_bak, 1, 'DU_avg_bak', up_var, record, ierr)
          if (ierr == nf90_noerr) then
            call nc_write_var(ncrst, DV_avg_bak, 1, 'DV_avg_bak', vp_var, record, ierr)
            if (ierr /= nf90_noerr) then
              write(*,1) 'DV_avg_bak', record MYID
            endif
          else
            write(*,1) 'DU_avg_bak', record MYID
          endif
        else
          write(*,1) 'DV_avg2', record MYID
        endif
      else
        write(*,1) 'DU_avg2', record MYID
      endif
#  elif defined PRED_COUPLED_MODE
      call nc_write_var(ncrst, rufrc_bak(START_2D_ARRAY,nstp), 1, 'rufrc', up_var, record, ierr)
      if (ierr == nf90_noerr) then
        call nc_write_var(ncrst, rvfrc_bak(START_2D_ARRAY,nstp), 1, 'rvfrc', vp_var, record, ierr)

        if (ierr /= nf90_noerr) then
          write(*,1) 'rvfrc_bak', record MYID
        endif
      else
        write(*,1) 'rufrc_bak', record MYID
      endif
#  endif
      if (ierr /= nf90_noerr) goto 99                !--> ERROR
# endif

! 3D momentum components in XI- and ETA-directions.
      call nc_write_var(ncrst, u(START_2D_ARRAY,1,nstp), N,
     &                    vname(1,indxU), up_var, record, ierr)
      if (ierr == nf90_noerr) then
        call nc_write_var(ncrst, v(START_2D_ARRAY,1,nstp), N,
     &                    vname(1,indxV), vp_var, record, ierr)
        if (ierr /= nf90_noerr) then
          lvar=lenstr(vname(1,indxV))
          write(*,1) vname(1,indxV)(1:lvar), record MYID
        endif
      else
        lvar=lenstr(vname(1,indxU))
        write(*,1) vname(1,indxU)(1:lvar), record MYID
      endif
      if (ierr /= nf90_noerr) goto 99                !--> ERROR

! Tracer variables.
      call wrt_rst_trc(record)

# ifdef LMD_KPP
! Depth of planetary boundary layer [m].
      call nc_write_var(ncrst, hbls(START_2D_ARRAY,nstp), 1,
     &                    vname(1,indxHbls), rp_var, record, ierr)
      if (ierr /= nf90_noerr) then
        lvar=lenstr(vname(1,indxHbls))
        write(*,1) vname(1,indxHbls)(1:lvar), record MYID
        goto 99                                         !--> ERROR
      endif
# endif
# ifdef LMD_BKPP
! Thickness of KPP bottom boundary layer [m].
      call nc_write_var(ncrst, hbbl(START_2D_ARRAY,nstp), 1,
     &                    vname(1,indxHbbl), rp_var, record, ierr)
      if (ierr /= nf90_noerr) then
        lvar=lenstr(vname(1,indxHbbl))
        write(*,1) vname(1,indxHbbl)(1:lvar), record MYID
        goto 99                                         !--> ERROR
      endif
# endif
#endif /* SOLVE3D */


  1   format(/1x, '### ERROR: wrt_rst :: Cannot write variable ''',
     &             A, ''' into restart file, rec =', i6, 3x,A,i4/)
      goto 100
  99  if (may_day_flag == 0) may_day_flag=3
 100  continue

! Close or synchronize netCDF file to disk to the allow other
! processes to access data immediately after it is written.

#if defined MPI && !defined PARALLEL_FILES
      ierr=nf90_close (ncrst)
      if (nrpfrst > 0 .and. record >= nrpfrst) ncrst=-1
#else
      if (nrpfrst > 0 .and. record >= nrpfrst) then
        ierr=nf90_close (ncrst)
        ncrst=-1
      else
        ierr=nf90_sync(ncrst)
      endif
#endif
      if (ierr == nf90_noerr) then
#ifdef MPI_SILENT_MODE
        if (mynode == 0) then
#endif
          write(tstring,'(F18.8)') tdays
          i=1
          do while (i<18 .and. tstring(i:i)==' ')
            i=i+1
          enddo
          write(*,'(7x,A,1x,A,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
     &      'wrt_rst :: wrote restart, tdays =', tstring(i:i+8),
     &      'step =', ibuff(1),  'rec =', record, '/',  nrecrst
     &       MYID
#ifdef MPI_SILENT_MODE
        endif
#endif
      else
        write(*,'(/1x,2A/)')    '### ERROR: wrt_rst :: Cannot ',
     &                          'synchronize/close restart file.'
        if (may_day_flag == 0) may_day_flag=3
      endif

#if defined MPI && !defined PARALLEL_FILES
      if (mynode < NNODES-1) call MPI_Send (blank, 1, MPI_INTEGER,
     &                        mynode+1, 1, ocean_grid_comm, ierr)
#endif

      end subroutine wrt_rst  !]

!----------------------------------------------------------------------

      end module ocean_vars
