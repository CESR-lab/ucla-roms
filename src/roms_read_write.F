      module roms_read_write
      ! ROMS specific read write functionality
      !
      ! output_root_name is used as prefix of all output files
      ! Length is limited to 26 characters, because output_file_names
      ! are appended with frame and node number and are limited to 32

      use netcdf
      use ncvars
      use scalars
      use nc_read_write
!     use param

      implicit none

#include "cppdefs.opt"

      private

      ! Type ncforce contains all required netcdf variables and input data array for frc variable.
      ! Set the name and time_name of variable to match what is in the input forcing files.
      ! Use same name as bulk_frc variable name, e.g. uwnd has nc_uwnd for netcdf vars.
      ! type (ncvs) :: nc_uwnd  = ncvs( name='uwnd',       time_name='wnd_time'  )
      type,public  :: ncforce
        character(len=20)                 :: vname     ! name of variable in input file
        character(len=20)                 :: tname     ! time variable name for variable
        real,dimension(:,:,:),allocatable :: vdata     ! currently must be exact size of sub-domain (no buffer)
        real,dimension(:,:),allocatable   :: cdata     ! Array to read in coarse surface forcing data
        integer                           :: coarse=0  ! Flag to indicate interpolation of the input data
        integer                           :: ifile=0   ! Index to file in list of filenames
        integer                           :: irec =0   ! Record number in file
        integer                           :: it1 = 1, it2 = 2              ! used to cycle between correct entries of 'data' above
        real,dimension(2)                 :: times = -99 ! [-99,-99]       ! stores 2 times that go with 'data' above
      end type ncforce

      integer, target :: test_var
      character(len=99), public  :: output_root_name

      integer,public :: i0,i1,j0,j1                        ! i0=0 (west edge), or i0=1 interior

      integer,parameter :: max_options_string=2000
      character(len=max_options_string),public :: surf_forcing_strings=' '

      public :: set_frc_data,create_filename
      public :: read_output_root_name
      public :: init_output_indices
      public :: find_new_record
      public :: nc_check_units
      public :: put_global_atts
      public :: store_string_att

      contains

! ----------------------------------------------------------------------
      subroutine init_output_indices  ![
      ! outputting variables requires buffer points
      ! on the domain boundary only.
      ! interior sub-domain is still  1:nx,1:ny
      ! but e.g. South-West corner is 0:nx,0:ny
      ! we write this as ncwrite( var(i0:i1,j0:j1) )

      use dimensions
      implicit none

      i0=1; i1=nx; j0=1; j1=ny
      if (.not.west_exchng)  i0=0
      if (.not.south_exchng) j0=0
      if (.not.east_exchng)  i1=nx+1
      if (.not.north_exchng) j1=ny+1

!      print *, mynode, 'i0=',i0,'j0=',j0,'i1=',i1,'j1=',j1
!      print *, mynode, 'nx=',ny,'j0=',j0,'i1=',i1,'j1=',j1

      end subroutine init_output_indices  !]

! ----------------------------------------------------------------------
      subroutine set_frc_data(nc,var1d,var2d)  ![
      implicit none

      ! input/outputs
      type(ncforce),               intent(inout) :: nc    ! derived type containing all neccesary supporting data
      real,dimension(:)  ,optional,intent(out)   :: var1d ! time interpolated forcing data
      real,dimension(:,:),optional,intent(out)   :: var2d ! time interpolated forcing data

      ! local
      real    :: cff1,cff2  ! for time interpolations
      real    :: modtime    ! model time
      logical :: d1

      integer           :: it1,it2
      real,dimension(2) :: vtimes

      if (present(var1d)) then
        d1 = .true.
      else
        d1 = .false.
      endif

      modtime = tdays + 0.5*dt*sec2day    ! input data times in days. Interpolated time is based on 1/2 step time

      it1 = nc%it1
      it2 = nc%it2

      if (nc%times(it2) < modtime) then  ! need to refresh data

        if (nc%ifile == 0) then   ! the first time, we need to fill both slots
          call fill_frc_slice(nc,modtime,it1,d1)
        else ! cycle it1,it2 between 1 and 2
          it1 = mod(it1,2) + 1
          it2 = mod(it2,2) + 1
        endif

        call fill_frc_slice(nc,modtime,it2,d1)
      endif 

      nc%it1 = it1
      nc%it2 = it2

      ! Temporal interpolation

      vtimes  = nc%times
      cff1=( vtimes(it2)-modtime )/( vtimes(it2)-vtimes(it1) )
      cff2=( modtime-vtimes(it1) )/( vtimes(it2)-vtimes(it1) )

      if (d1) then
        var1d = cff1*nc%vdata(:,1,it1) + cff2*nc%vdata(:,1,it2)
      else
        var2d = cff1*nc%vdata(:,:,it1) + cff2*nc%vdata(:,:,it2)
      endif


      if (vtimes(it1)>modtime.or.vtimes(it2)<modtime) then
        print * ,'set_frc_data :: ', nc%vname,frcfile(nc%ifile)
        print * ,'set_frc_data :: ', vtimes(it1), vtimes(it2), modtime
        error stop 'ERROR: set_frc_data :: time interpolation error'
      endif

      end subroutine set_frc_data  !]
      
! ----------------------------------------------------------------------
      subroutine find_new_record(vname,tname,time,ifile,irec,vtime) ![
      ! Finds the index for the filename and record number of the first record
      ! with a time that is larger than the model time. It will stores
      ! the corresponding time of the record in vtime.

      ! When called for the first time, it will find the last record with a 
      ! time smaller than the model time.
      implicit none

      ! input/output
      character(len=*)              ,intent(in)   :: vname   ! name of the variable
      character(len=*)              ,intent(in)   :: tname   ! name of the corresponding time
      real                          ,intent(in)   :: time    ! target time
      integer                       ,intent(inout):: ifile   ! index for file name in list
      integer                       ,intent(inout):: irec    ! variable record entry
      real                          ,intent(out)  :: vtime   ! read in variable time

      ! local
      integer, dimension(1) :: dimids     ! time dimension ID
      real    :: time_old
      logical :: found_var                ! is variable in file
      logical :: found_rec                ! is correct record found of var
      logical :: found_var_ever           ! if variable was found in any previous files
      integer :: ncid, ierr
      integer :: nfiles                   ! total number of forcing files
      integer :: irec_old, ifile_old
      integer :: max_recs                 ! total variable records in file
      integer :: vid,tid                  ! ID's of variable and its time
      real    :: cycle_length             ! cyclical data
      logical :: first                    ! first time call for this variable

      real,dimension(1)    :: ltime

      if (ifile==0) then                                             ! first time for this variable
        first = .true.
        ifile = 1
        irec_old  = 0
        ifile_old = 0
      else
        first = .false.
      endif
      
      found_rec      = .false.
      found_var_ever = .false.

      do while (.not.found_rec .and. ifile <= max_frc)               ! while record not found & still more files to check

        found_var = .false.                                          ! reset for new file
        ierr = nf90_open(frcfile(ifile), nf90_nowrite, ncid)
        if (ierr/=nf90_noerr) 
     &    call handle_ierr(ierr,'Find_new_record, opening: ',frcfile(ifile))

        ierr = nf90_inq_varid(ncid, vname,vid)                       ! check if this file contains the variable

        if (ierr==nf90_noerr) then                                   ! variable found, now check if file times contain model time

          found_var      = .true.
          found_var_ever = .true.
           
          ierr = nf90_inq_varid(ncid,tname,tid)                      ! get time ID
          if (ierr/=nf90_noerr) 
     &      call handle_ierr(ierr,'Find_new_record, searching: ',tname)

          ! Check if the time has a cycle length
          ierr=nf90_get_att(ncid, tid, 'cycle_length', cycle_length)
          if (ierr==nf90_noerr) then
            call find_new_record_cycle(ncid,tid,time,irec,vtime,cycle_length,first)
            found_rec = .true.
          elseif (ierr==nf90_enotatt) then                           ! non-cyclical var won't have attribute
            cycle_length = -1                                        ! nf90_enotatt = attribute not found
          else                                                       ! any other netcdf error is a problem.
            call handle_ierr(ierr,'find_new_record: cycle length att')
          endif

          ! Number of records in file
          ierr = nf90_inquire_variable(ncid,tid, dimids = dimids)
          ierr = nf90_inquire_dimension(ncid, dimids(1), len=max_recs)

          do while (irec<max_recs .and. .not. found_rec)             ! Search through the records until correct time is found

            irec=irec+1                                              ! why not use ncread????
            ierr=nf90_get_var(ncid,tid,ltime,(/irec/), (/1/))        ! Get vtime (nf90 needs array hence syntax)
            vtime = ltime(1)

            if (vtime>time) then                                     ! correct time
              found_rec = .true.
            else                                                     ! not yet
              irec_old = irec
              ifile_old = ifile                                      ! only needed initially. See label (A)
              time_old = vtime
            endif

          enddo                                                      ! while not found and irec<nrecs

        endif                                                        ! found var in file

        if (.not.found_rec) then                                     ! correct record not found, try next file
          ifile = ifile+1 
          irec = 0                                                   ! Reset irec for new file
        endif

        ierr=nf90_close (ncid)

      enddo                            ! search through files: while .not. found_rec .and. ifile<=nfiles.


      if (.not. found_var_ever) then
        write(*,'(/1x,4A,I3/)') 'ERROR: find_new_record:: ',
     &              'Could not find var: ', vname, ' mynode=', mynode
        error stop 'find_new_record: var not found in forcing files'
      endif
      if (.not. found_rec) then
        write(*,'(/1x,3A,I3/)')
     &  'ERROR: find_new_record: Found variable, but ran out of time records for ',
     &                                      vname, ' mynode=', mynode
        error stop 'find_new_record: found var, but ran out of time records'
      endif

      if (first) then  ! Handle initital timestep. Need vtime < model_time. As per (A) above.

        ! find_new_rec_cycle deals internally with this
        if(cycle_length<0) then     ! not cyclical data
          if (irec==1) then         ! model time between input files.
            ifile = ifile_old
          endif
          irec = irec_old
          vtime = time_old
        endif

      endif ! Handle initital timestep.

      if (irec==0) then
        write(*,'(/1x,4A,I3/)') 'ERROR: find_new_record:: ',
     &  'First available forcing record is past current time for var: ',
     &     vname, ' mynode=', mynode
        error stop 'First available forcing record is past current time'
      endif

      end subroutine find_new_record  !]

! ----------------------------------------------------------------------
      subroutine find_new_record_cycle(ncid,tid,time,irec,vtime,cycle_length,first) ![
      ! Find the first record that has vtime>time, or if it's
      ! the first time calling, the record preceeding that one.
      implicit none

      ! input/output
      integer                       ,intent(in)  :: ncid ! Handle of open netcdf file
      integer                       ,intent(in)  :: tid  ! ID of time variable
      real                          ,intent(in)  :: time ! target time
      integer                       ,intent(out) :: irec ! correct record number
      real                          ,intent(out) :: vtime! correct record time
      real                          ,intent(in)  :: cycle_length
      logical                       ,intent(in)  :: first

      ! local
      integer                       :: icycle,max_rec,ierr
      integer, dimension(1) :: dimids     ! time dimension ID
      real,allocatable,dimension(:) :: tim

      ! Number of records in file
      ierr = nf90_inquire_variable(ncid,tid, dimids = dimids)
      ierr = nf90_inquire_dimension(ncid, dimids(1), len=max_rec)
      
      allocate(tim(max_rec))
      ierr=nf90_get_var(ncid,tid,tim)

      icycle = floor((time - tim(max_rec))/cycle_length) + 1

      do irec = 1,max_rec
         vtime = tim(irec) + icycle*cycle_length
         if (vtime>time) then
           exit
         endif
      enddo

      if (first) then
        irec = irec-1;
        if (irec==0) then
          icycle = icycle-1
          irec = max_rec
        endif
        vtime = tim(irec) + icycle*cycle_length
      endif

      deallocate(tim)

      end subroutine find_new_record_cycle !]

! ----------------------------------------------------------------------
      subroutine fill_frc_slice(nc,modtime,it,d1) ![
      ! Fill a time slice of forcing data
      implicit none

      ! input/outputs
      type(ncforce)      ,intent(inout) :: nc      ! derived type containing all neccesary supporting data
      real               ,intent(in)    :: modtime ! current time of the model
      integer            ,intent(in)    :: it      ! index (1 or 2) of the forcing data
      logical            ,intent(in)    :: d1      ! 1 or 2 dimensional forcing data

      ! local
      integer :: ierr,ncid

      character(len=20)  :: vname  ! name of variable in input file
      character(len=20)  :: tname  ! time variable name for variable
      integer            :: irec,ifile

      irec  = nc%irec
      ifile = nc%ifile
      vname = nc%vname
      tname = nc%tname

      call find_new_record(vname,tname,modtime,ifile,irec,
     &                         nc%times(it) )

      ierr=nf90_open(frcfile(ifile),nf90_nowrite, ncid)
      if (nc%coarse) then
        print *,'coarse to fine not yet implemented here'
        stop
!       cdata => nc_var%cdata
!       ncread(ncid,vname,cdata,(/1,1,irec/))
!       call exchange_2d
!       call coarse2fine(cdata,vdata(:,:,it))
      else
        !! Allocate vdata to be the right size for the (DevinD - this should have been put somewhere clearer...
        !! subdomain (including whether we need a buffer or not
        if (d1) then
          call ncread(ncid,vname,nc%vdata(:,1,it),(/1,irec/))
        else
          call ncread(ncid,vname,nc%vdata(:,:,it),(/1,1,irec/))
        endif
      endif
      ierr = nf90_close(ncid)

      if (mynode == 0) 
     &      call display_read_time_to_log(vname,nc%times(it),irec )

      nc%irec  = irec
      nc%ifile = ifile

      end subroutine fill_frc_slice !]

! ----------------------------------------------------------------------
      subroutine display_read_time_to_log(var_name, var_time, var_irec)  ![
      implicit none   
        
      character(len=*), intent(in) :: var_name
      integer,          intent(in) :: var_irec
      real,             intent(in) :: var_time
          
      write(*,'(9x,A,A12,5x,A,G14.6,1x,A,I4)') 
     &  'set_frc :: ',         var_name,
     &  'input time (days) =', var_time,         ! G is number for time 
     &  'rec =',               var_irec

      end subroutine display_read_time_to_log  !]

! ----------------------------------------------------------------------
      subroutine read_output_root_name(keyword, kwlen)  ![
      ! read from .in file the root-name
      ! for all output netcdf result files
      ! this is called from read_inp.F
      ! The following needs to be in the roms.in input file to use this
      ! subroutine (file name needs to be 5 spaces from left margin):
      !output_root_name:
      !     rootname (<- insert desired rootname)

      implicit none

      ! inputs
      character(len=*) :: keyword
      integer          :: kwlen                            ! keyword length
      ! local
      integer           :: ierr = nf90_noerr
      integer,parameter :: input=15

      call cancel_kwd (keyword(1:kwlen), ierr)             ! remove keyword from keyword list for read_inp.F
                                                           ! to know it's accounted for
      read(input,'(A)',err=95) output_root_name            ! read & save filename root

      goto 100
                                                           ! error for read(input,'(A)',err=95)
  95  write(*,'(/1x,4A/)') '### ERROR: roms_read_write :: Cannot ',
     &                    'read entry ''', keyword(1:kwlen), '''.'
      error stop
 100  continue


      end subroutine read_output_root_name  !]

! ----------------------------------------------------------------------
      subroutine create_filename(fname,total_rec)  ![
      ! takes a filename for output and inserts the mpi node number 
      ! (if appropriate), and the label that matches the current record
      ! number
      implicit none

      ! input/output
      character(len=*),  intent(inout) :: fname              ! desired netcdf file name
      integer,           intent(in)    :: total_rec          ! total netcdf records so far for these variables

      ! local
      integer           :: lfnm, lvar, lenstr, lstr, ierr=0
      character(len=64) :: fname_tmp                         ! taken from read_inp.F

      fname_tmp=fname
      lstr=lenstr(fname_tmp)
#if defined MPI && defined PARALLEL_FILES
      call insert_node(fname_tmp, lstr, mynode, NNODES, ierr)  ! insert MPI node numbers to file name
#endif
      fname = fname_tmp(1:lstr)
      call insert_time_index(fname,lstr,total_rec,ierr)     ! insert time index to file name

      end subroutine create_filename !]

! ----------------------------------------------------------------------
      subroutine coarse2fine(cdata,fdata) ![
      ! interpolates from (possibly partitioned) coarse to fine surface forcing data
      implicit none


      ! input/output
      real,dimension(:,:),intent(in) :: cdata   ! coarse data
      real,dimension(:,:),intent(out):: fdata   ! fine data

      ! local
      integer :: nx,ny
      integer :: i,j
      integer :: ic,jc
      real  :: a,b

      a = 0.75
      b = 0.25


!     call exchange_2d(cdata)

      do j=1,ny,2
        do i=1,nx,2
          fdata(i  ,j  ) = a*(a*cdata(ic,jc  )+b*cdata(ic+1,jc  ))
     &                   + b*(a*cdata(ic,jc+1)+b*cdata(ic+1,jc+1))
          fdata(i+1,j  ) = a*(b*cdata(ic,jc  )+a*cdata(ic+1,jc  ))
     &                   + b*(b*cdata(ic,jc+1)+a*cdata(ic+1,jc+1))
          fdata(i  ,j+1) = b*(a*cdata(ic,jc  )+b*cdata(ic+1,jc  ))
     &                   + a*(a*cdata(ic,jc+1)+b*cdata(ic+1,jc+1))
          fdata(i+1,j+1) = b*(b*cdata(ic,jc  )+a*cdata(ic+1,jc  ))
     &                   + a*(b*cdata(ic,jc+1)+a*cdata(ic+1,jc+1))
        enddo
      enddo

      end subroutine coarse2fine !]

!-----------------------------------------------------------------------
      subroutine nc_check_units( ncid, var_name, var_units )  ![

      ! Ensure units in file are same as roms units

      implicit none

      ! input/output
      integer,          intent(in) :: ncid
      character(len=*), intent(in) :: var_name
      character(len=*), intent(in) :: var_units

      ! local
      integer           :: ierr !, ncid
      integer           :: var_id
      character(len=20) :: file_var_units ! variable units in forcing file

      ! file must be open already:

      ierr = nf90_inq_varid(ncid, var_name, var_id)                            ! get variable ID (varid)
      if (ierr/=0) call handle_ierr(ierr,'nc_check_units:: cannot find var')
      ierr = nf90_get_att(ncid, var_id, 'units', file_var_units)               ! get units
      if (ierr/=0) call handle_ierr(ierr,'nc_check_units:: cannot find var units')

      if (var_units /= file_var_units) then
        write(*,'(/3x,3A/5x,3A/)') 'ERROR: roms_read_write.F :: ',
     &  'nc_check_units: wrong units for var: ', var_name,
     &  trim(file_var_units) ,' -> input file should be -> ', var_units
        error stop                                                              ! stop simulation, fatal error
      end if

      end subroutine nc_check_units  !]

! ---------------------------------------------------------------------
      subroutine put_global_atts (ncid, ierr)  ![

      use param
!      use surf_flux !, only: dSSSdt
      use eos_vars
      use ncvars
      use scalars
      use scoord
      use strings
      use netcdf

      implicit none
      integer :: ncid, ierr

#if defined MPI && defined PARALLEL_FILES
      integer :: ibuff(4), npart

      npart=4                                    ! Put global attribute 'partition' which identifies subdomain
      ibuff(1)=mynode                            ! within the processor grid individually for each file.
      ibuff(2)=nnodes
      if (WESTERN_MPI_EDGE) then
        ibuff(3)=iSW_corn+iwest
      else
        ibuff(3)=iSW_corn+iwest+1
      endif
      if (SOUTHERN_MPI_EDGE) then
        ibuff(4)=jSW_corn+jsouth
      else
        ibuff(4)=jSW_corn+jsouth+1
      endif

      ierr=nf90_put_att(ncid, nf90_global, 'partition', ibuff)
#endif

      ierr=nf90_put_att(ncid, nf90_global, 'title',title)         ! Title of the run

#ifndef ANA_GRID
      ierr=nf90_put_att(ncid, nf90_global, 'grid_file',grdname)
#endif
#ifndef ANA_INITIAL
      ierr=nf90_put_att(ncid, nf90_global, 'init_file',ininame)
#endif

      ierr=nf90_put_att(ncid, nf90_global, 'ntimes', ntimes)         ! Time stepping parameters
      ierr=nf90_put_att(ncid, nf90_global, 'ndtfast',ndtfast)
      ierr=nf90_put_att(ncid, nf90_global, 'dt', dt)
      ierr=nf90_put_att(ncid, nf90_global, 'dtfast', dtfast)

! Vertical coordinate type, control parameters "theta_s", "theta_b",
! "hc", and stretching curves "Cs_w", "Cs_r" at W- and RHO-points.

#ifdef SOLVE3D
# ifdef VERT_COORD_TYPE_SM09
      ierr=nf90_put_att(ncid, nf90_global, 'VertCoordType', 'SM09')
# elif defined VERT_COORD_TYPE_KAU
      ierr=nf90_put_att(ncid, nf90_global, 'VertCoordType', 'KAU')
# elif defined VERT_COORD_TYPE_SH94
      ierr=nf90_put_att(ncid, nf90_global, 'VertCoordType', 'SH94')
# endif

      ierr=nf90_put_att(ncid, nf90_global, 'theta_s', theta_s)
      ierr=nf90_put_att(ncid, nf90_global, 'theta_b', theta_b)
      ierr=nf90_put_att(ncid, nf90_global, 'hc', hc)
      ierr=nf90_put_att(ncid, nf90_global, 'Cs_w', Cs_w)
      ierr=nf90_put_att(ncid, nf90_global, 'Cs_r', Cs_r)

! Boussinesq reference density and EOS coefficients for linear EOS.

      ierr=nf90_put_att(ncid, nf90_global, 'rho0', rho0)
      ierr=nf90_put_att (ncid, nf90_global, 'rho0_units', 'kg/m^3')
# ifndef NONLIN_EOS
      ierr=nf90_put_att(ncid, nf90_global,'Tcoef', Tcoef)
      ierr=nf90_put_att (ncid, nf90_global,'Tcoef_units', 'deg Celsius-1')

      ierr=nf90_put_att(ncid, nf90_global, 'T0', T0)
      ierr=nf90_put_att (ncid, nf90_global, 'T0_units', 'deg Celsius')

#  ifdef SALINITY
      ierr=nf90_put_att(ncid, nf90_global, 'Scoef', Scoef)
      ierr=nf90_put_att (ncid, nf90_global, 'Scoef_units', 'PSU-1')

      ierr=nf90_put_att(ncid, nf90_global, 'S0', S0)
      ierr=nf90_put_att (ncid, nf90_global, 'S0_units', 'PSU')
#  endif
# endif
#endif /* SOLVE3D */

#ifdef UV_VIS2
      ierr=nf90_put_att(ncid,nf90_global, 'visc2', visc2)            ! Horizontal viscosity and mixing coefficients.
      ierr=nf90_put_att (ncid,nf90_global, 'visc2_units', 'm^2/s')

      ierr=nf90_put_att(ncid, nf90_global, 'gamma2', gamma2)         ! Lateral slipperiness parameter
#endif
#ifdef SOLVE3D
# ifdef TS_DIF2
      ierr=nf90_put_att(ncid,nf90_global, 'tnu2', tnu2)
      ierr=nf90_put_att (ncid,nf90_global, 'tnu2_units', 'm^2/s')
# endif
# ifdef MY25_MIXING
#  ifdef Q2_DIF2
      ierr=nf90_put_att(ncid,nf90_global, 'q2nu2', q2nu2)
      ierr=nf90_put_att (ncid,nf90_global, 'q2nu2_units', 'm^2/s')
#  endif
#  ifdef Q2_DIF4
      ierr=nf90_put_att(ncid,nf90_global,'q2nu4', q2nu4)
      ierr=nf90_put_att (ncid,nf90_global,'q2nu4_units', 'm^4/s')
#  endif
# endif
# if !defined LMD_MIXING && !defined BVF_MIXING && !defined PP_MIXING\
                         && !defined MY2_MIXING && !defined MY25_MIXING

      ierr=nf90_put_att(ncid, nf90_global, 'Akv_bak', Akv_bak)       ! Background vertical viscosity tracer mixing
      ierr=nf90_put_att (ncid, nf90_global, 'Akv_bak_units', 'm^2/s')

      ierr=nf90_put_att(ncid, nf90_global, 'Akt_bak', Akt_bak)
      ierr=nf90_put_att (ncid, nf90_global, 'Akt_bak_units', 'm^2/s')
# endif
# ifdef MY25_MIXING
      ierr=nf90_put_att(ncid, nf90_global, 'Akq_bak',Akq_bak)
      ierr=nf90_put_att (ncid, nf90_global, 'Akq_bak_units', 'm^2/s')
# endif
#endif

#if  defined T_FRC_BRY || defined M2_FRC_BRY || defined TNUDGING \
  || defined Z_FRC_BRY || defined M3_FRC_BRY || defined M2NUDGING \
                       || defined M3NUDGING  || defined WKB_FRC_BRY

      ierr=nf90_put_att(ncid, nf90_global, 'ubind', ubind) ! Binding velocity for external data forcing at open boundaries
      ierr=nf90_put_att (ncid, nf90_global, 'ubind_units', 'm/s')
#endif

#ifdef SPONGE
      ierr=nf90_put_att(ncid, nf90_global, 'v_sponge', v_sponge)     ! Maximum viscosity in sponge layer
      ierr=nf90_put_att (ncid, nf90_global, 'v_sponge_units','m^2/s')
#endif

!#if defined SFLX_CORR && defined SALINITY
!      ierr=nf90_put_att(ncid, nf90_global, 'dSSSdt',
!     &                                        dSSSdt*(100.*day2sec)) ! Sea surface salinity correction "piston velocity"
!      ierr=nf90_put_att (ncid, nf90_global, 'dSSSdt_units','cm/day')
!#endif

      ierr=nf90_put_att(ncid,nf90_global, 'rdrg', rdrg)              ! Bottom drag parameters.
      ierr=nf90_put_att (ncid,nf90_global, 'rdrg_units', 'm/s')

      ierr=nf90_put_att(ncid,nf90_global, 'rdrg2',rdrg2)
      ierr=nf90_put_att (ncid,nf90_global, 'rdrg2_units',
     &                                      'nondimensional')
#ifdef SOLVE3D
      ierr=nf90_put_att(ncid,nf90_global, 'Zob', Zob)
      ierr=nf90_put_att (ncid,nf90_global, 'Zob_units', 'm')
#endif

      ierr=nf90_put_att (ncid, nf90_global, 'SRCS', srcs)  ! List of source-code file names and activated CPP-switches
      ierr=nf90_put_att (ncid, nf90_global, 'CPPS', cpps)

      ! string set added to get values from modules without struggling with circular references
      ! as those modules use roms_read_write so this module can't use them.
      ! rather the string is filled in the module's init and can be written here to all files.
      ierr=nf90_put_att (ncid, nf90_global, 'surf_forcing_strings', surf_forcing_strings)

      end subroutine put_global_atts  !]

! ---------------------------------------------------------------------
      subroutine store_string_att( char_string, string_add )  ![
      ! add settings that affect the solution to a string
      ! that will be added to global attributes of all nc files
      implicit none

      ! input
      character(len=*),intent(inout) :: char_string
      character(len=*),intent(in)    :: string_add

      ! local
      integer :: is, ie

      is=len_trim( char_string )
      is=is+2
      ie=len_trim( string_add )
      ie=ie+is-1

      if (ie > max_options_string) then
        print *, 'ERROR roms_read_write:: ',
     &    'store_string_att character full=', trim(string_add)
        stop 'ROMS stopped'
      endif

      char_string(is:ie) = string_add

      end subroutine store_string_att  !]

! ---------------------------------------------------------------------

      end module roms_read_write
