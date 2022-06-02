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
      use roms_read_write
      use dimensions
      implicit none

#ifndef ANA_BRY
      real(kind=8) :: bry_cycle, bry_time(2)
      integer :: itbry, ntbry, bry_ncycle, ibry, max_bry,
     &                  bry_rec,  bry_id, bry_time_id

# ifdef OBC_WEST
#  ifdef SOLVE3D
#   ifdef T_FRC_BRY
      integer t_west_id(NT)
#   endif
#  endif
# endif

# ifdef OBC_EAST
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      integer w_east_id
#   endif
#   ifdef T_FRC_BRY
      integer t_east_id(NT)
#   endif
#  endif
# endif

# ifdef OBC_SOUTH
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      integer w_south_id
#   endif
#   ifdef T_FRC_BRY
      integer t_south_id(NT)
#   endif
#  endif
# endif

# ifdef OBC_NORTH
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      integer w_north_id
#   endif
#   ifdef T_FRC_BRY
      integer t_north_id(NT)
#   endif
#  endif
# endif
#endif  /* ~ANA_BRY */


# ifdef OBC_WEST
#  ifdef Z_FRC_BRY
      real,allocatable,dimension(:)   :: zeta_west         ! all allocatabel arrays so could remove cppflag clutter...
#  endif
#  ifdef M2_FRC_BRY
      real,allocatable,dimension(:)   :: ubar_west, vbar_west
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      real,allocatable,dimension(:,:)   :: u_west, v_west, w_west
#   endif
#   ifdef T_FRC_BRY
      real,allocatable,dimension(:,:,:) :: t_west
      real,allocatable,dimension(:,:,:,:) :: t_west_dt
!      logical,dimension(nt)             :: t_w_default     ! whether boundary file contains tracer
#   endif
#  endif
# endif

# ifdef OBC_EAST
#  ifdef Z_FRC_BRY
      real,allocatable,dimension(:)   :: zeta_east
#  endif
#  ifdef M2_FRC_BRY
      real,allocatable,dimension(:)   :: ubar_east, vbar_east
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      real,allocatable,dimension(:,:)   :: u_east, v_east, w_east
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
#  endif
#  ifdef M2_FRC_BRY
      real,allocatable,dimension(:)   :: ubar_south, vbar_south
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      real,allocatable,dimension(:,:)   :: u_south, v_south, w_south
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
#  endif
#  ifdef M2_FRC_BRY
      real,allocatable,dimension(:)   :: ubar_north, vbar_north
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      real,allocatable,dimension(:,:)   :: u_north, v_north, w_north
#   endif
#   ifdef T_FRC_BRY
      real,allocatable,dimension(:,:,:)   :: t_north
      real,allocatable,dimension(:,:,:,:) :: t_north_dt
#   endif
#  endif
# endif

#if !defined ANA_BRY && ( defined OBC_WEST   || defined OBC_EAST  || \
                          defined OBC_SOUTH  || defined OBC_NORTH )  \
                     && ( defined M2_FRC_BRY || defined Z_FRC_BRY || \
                          defined M3_FRC_BRY || defined T_FRC_BRY )

      public :: get_bry_all, set_bry_all_tile

#endif

      type (ncforce) :: nc_z_w = ncforce( vname='zeta_west', tname='bry_time')
      type (ncforce) :: nc_ub_w = ncforce( vname='ubar_west', tname='bry_time')
      type (ncforce) :: nc_vb_w = ncforce( vname='vbar_west', tname='bry_time')
      type (ncforce) :: nc_u_w = ncforce( vname='u_west', tname='bry_time')
      type (ncforce) :: nc_v_w = ncforce( vname='v_west', tname='bry_time')
!      type (ncforce), dimension(nt) :: nc_t_w = ncforce( vname=' ', tname='bry_time')

      type (ncforce) :: nc_z_e = ncforce( vname='zeta_east', tname='bry_time')
      type (ncforce) :: nc_ub_e = ncforce( vname='ubar_east', tname='bry_time')
      type (ncforce) :: nc_vb_e = ncforce( vname='vbar_east', tname='bry_time')
      type (ncforce) :: nc_u_e = ncforce( vname='u_east', tname='bry_time')
      type (ncforce) :: nc_v_e = ncforce( vname='v_east', tname='bry_time')
!      type (ncforce), dimension(nt) :: nc_t_e = ncforce( vname=' ', tname='bry_time')

      type (ncforce) :: nc_z_s = ncforce( vname='zeta_south', tname='bry_time')
      type (ncforce) :: nc_ub_s = ncforce( vname='ubar_south', tname='bry_time')
      type (ncforce) :: nc_vb_s = ncforce( vname='vbar_south', tname='bry_time')
      type (ncforce) :: nc_u_s = ncforce( vname='u_south', tname='bry_time')
      type (ncforce) :: nc_v_s = ncforce( vname='v_south', tname='bry_time')
!      type (ncforce), dimension(nt) :: nc_t_s = ncforce( vname=' ', tname='bry_time')

      type (ncforce) :: nc_z_n = ncforce( vname='zeta_north', tname='bry_time')
      type (ncforce) :: nc_ub_n = ncforce( vname='ubar_north', tname='bry_time')
      type (ncforce) :: nc_vb_n = ncforce( vname='vbar_north', tname='bry_time')
      type (ncforce) :: nc_u_n = ncforce( vname='u_north', tname='bry_time')
      type (ncforce) :: nc_v_n = ncforce( vname='v_north', tname='bry_time')
!      type (ncforce), dimension(nt) :: nc_t_n = ncforce( vname=' ', tname='bry_time')

#ifdef NHMG
      type (ncforce) :: nc_w_w = ncforce( vname='w_west', tname='bry_time')
      type (ncforce) :: nc_w_e = ncforce( vname='w_east', tname='bry_time')
      type (ncforce) :: nc_w_s = ncforce( vname='w_south', tname='bry_time')
      type (ncforce) :: nc_w_n = ncforce( vname='w_north', tname='bry_time')
#endif

      contains

!----------------------------------------------------------------------
      subroutine init_arrays_boundary  ![
      implicit none

      ! local
      integer :: itrc, ifile, ncid, ierr, varid
      logical :: found_bry

# ifdef OBC_WEST
#  ifdef Z_FRC_BRY
      allocate( zeta_west(0:Mm+1) )
      allocate( nc_z_w%vdata(j0:j1,1,2) )
#  endif
#  ifdef M2_FRC_BRY
      allocate( ubar_west(0:Mm+1), vbar_west(0:Mm+1) )
      allocate( nc_ub_w%vdata(j0:j1,1,2) )
      allocate( nc_vb_w%vdata( 1:j1,1,2) )
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      allocate( u_west(0:Mm+1,N), v_west(0:Mm+1,N) )
      allocate( nc_u_w%vdata(j0:j1,nz,2) )
      allocate( nc_v_w%vdata( 1:j1,nz,2) )
#    ifdef NHMG
      allocate( w_west(0:Mm+1,N), nc_w_w%vdata(j0:j1,nz,2) )
#    endif
#   endif
#   ifdef T_FRC_BRY

      allocate( t_west_dt(0:Mm+1,N,2,NT) )
      allocate( t_west(0:Mm+1,N,NT) )
!      do itrc=1,nt
!        nc_t_w(itrc)%vname = trim(t_vname(itrc)) / / '_west'
!      enddo
!
!#    ifndef DEFAULT_BRY_VALUES
!      do itrc=1,nt
!        allocate( nc_t_w(itrc)%vdata(j0:j1,nz) )
!      enddo
!#    else
!      ! find first boundary file and check if tracer exists
!      ! assume all in one file so stop looking in additional files
!      ! if a boundary term is found in a file.
!      ifile=1; found_bry=.false.
!      do while (ifile <= max_frc .and. .not. found_bry)
!        ierr=nf90_open(frcfile(ifile), nf90_nowrite, ncid)
!        if (ierr/=0) call handle_ierr(ierr,'init_arrays_boundary, opening: ',frcfile(ifile))
!        do itrc=1,nt
!          ierr=nf90_inq_varid(ncid, nc_t_w(itrc)%vname, varid)
!          if (ierr==nf90_noerr) then
!            allocate( nc_t_w(itrc)%vdata(j0:j1,nz) )       ! SHOULD ONLY BE ALLOCATED FOR REQUIRED TRACERS, not all!
!            nc_t_w(itrc)%tbry = .true.
!            found_bry = .true.
!          endif
!
!          ! set the default here straight to the boundary value (not data)
!          ! then put if statement in set_frc_data and only do something if needed.
!
!          ! add tbry to ncforce, cleaner than extra array.
!
!        enddo
!        ifile=ifile+1
!      enddo
!      !
!      do itrc=1,nt
!        if (nc_t_w(itrc)%tbry) then
!          call set_defaults2?
!        endif
!      enddo
!#    endif /* DEFAULT_BRY_VALUES */

#   endif
#  endif
# endif

# ifdef OBC_EAST
#  ifdef Z_FRC_BRY
      allocate( zeta_east(0:Mm+1) )
      allocate( nc_z_e%vdata(j0:j1,1,2) )
#  endif
#  ifdef M2_FRC_BRY
      allocate( ubar_east(0:Mm+1), vbar_east(0:Mm+1) )
      allocate( nc_ub_e%vdata(j0:j1,1,2) )
      allocate( nc_vb_e%vdata( 1:j1,1,2) )
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      allocate( u_east(0:Mm+1,N), v_east(0:Mm+1,N) )
      allocate( nc_u_e%vdata(j0:j1,nz,2) )
      allocate( nc_v_e%vdata( 1:j1,nz,2) )
#    ifdef NHMG
      allocate( w_east(0:Mm+1,N), nc_w_e%vdata(j0:j1,nz,2) )
#    endif
#   endif
#   ifdef T_FRC_BRY
      allocate( t_east(0:Mm+1,N,NT), t_east_dt(0:Mm+1,N,2,NT) )
#   endif
#  endif
# endif

# ifdef OBC_SOUTH
#  ifdef Z_FRC_BRY
      allocate( zeta_south(0:Lm+1) )
      allocate( nc_z_s%vdata(i0:i1,1,2) )
#  endif
#  ifdef M2_FRC_BRY
      allocate( ubar_south(0:Lm+1), vbar_south(0:Lm+1) )
      allocate( nc_ub_s%vdata( 1:i1,1,2) )
      allocate( nc_vb_s%vdata(i0:i1,1,2) )
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      allocate( u_south(0:Lm+1,N), v_south(0:Lm+1,N) )
      allocate( nc_u_s%vdata( 1:i1,nz,2) )
      allocate( nc_v_s%vdata(i0:i1,nz,2) )
#    ifdef NHMG
      allocate( w_south(0:Lm+1,N), nc_w_s%vdata(i0:i1,nz,2) )
#    endif
#   endif
#   ifdef T_FRC_BRY
      allocate( t_south(0:Lm+1,N,NT), t_south_dt(0:Lm+1,N,2,NT) )
#   endif
#  endif
# endif

# ifdef OBC_NORTH
#  ifdef Z_FRC_BRY
      allocate( zeta_north(0:Lm+1) )
      allocate( nc_z_n%vdata(i0:i1,1,2) )
#  endif
#  ifdef M2_FRC_BRY
      allocate( ubar_north(0:Lm+1), vbar_north(0:Lm+1) )
      allocate( nc_ub_n%vdata( 1:i1,1,2) )
      allocate( nc_vb_n%vdata(i0:i1,1,2) )
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
      allocate( u_north(0:Lm+1,N), v_north(0:Lm+1,N) )
      allocate( nc_u_n%vdata( 1:i1,nz,2) )
      allocate( nc_v_n%vdata(i0:i1,nz,2) )
#    ifdef NHMG
      allocate( w_north(0:Lm+1,N), nc_w_n%vdata(i0:i1,nz,2) )
#    endif
#   endif
#   ifdef T_FRC_BRY
      allocate( t_north(0:Lm+1,N,NT), t_north_dt(0:Lm+1,N,2,NT) )
#   endif
#  endif
# endif

      end subroutine init_arrays_boundary !]

!----------------------------------------------------------------------
#if !defined ANA_BRY && ( defined OBC_WEST   || defined OBC_EAST  || \
                          defined OBC_SOUTH  || defined OBC_NORTH )  \
                     && ( defined M2_FRC_BRY || defined Z_FRC_BRY || \
                          defined M3_FRC_BRY || defined T_FRC_BRY )

      subroutine get_bry_all(ierr_ext)     ![ Read side boundary forcing

      use param
      use tracers, only: t_vname, set_bry_tracer_defaults
      use scalars
      use netcdf
      use roms_read_write

      implicit none                        ! fields from boundary file
      integer ierr_ext

# include "netcdf.inc"
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

      character(len=47) :: v_name ! Variable name after concatenating DevinD

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
      ierr=nf90_noerr ; ierr_all=nf90_noerr

      if (bry_rec < 0) then  !<-- initialization or file switch
        if (bry_id == -1) then
          ierr=nf90_open(bry_file(ibry), nf90_nowrite, bry_id)
          if (ierr /= nf90_noerr) write(*,'(/1x,4A/12x,A/)')
     &          '### ERROR: get_bry_all :: Cannot open file ''',
     &              trim(bry_file(ibry)),   ''' for reading.',
     &                                      nf90_strerror(ierr)

        endif
        if (ierr == nf90_noerr) then
          ierr=nf90_inq_varid(bry_id, 'bry_time',  bry_time_id)
          if (ierr == nf90_noerr) then
            call find_record(bry_id, bry_time_id,  bry_rec, ntbry,
     &              bry_ncycle, bry_cycle,  ibry, max_bry,  ierr)

            if (bry_rec < 0) then
              ierr=nf90_close(bry_id) ; bry_id=-1
              ibry=ibry+1 ;  goto 100  !--> proceed to the next file
            endif

          else
            write(*,1) 'bry_time', bry_file(ibry)(1:lsrc)
          endif
          if (ierr_all == nf90_noerr) ierr_all=ierr



# ifdef OBC_WEST
          if (west_bry_active) then
#  ifdef SOLVE3D
#   ifdef T_FRC_BRY

            do itrc=1,NT ! THIS CAN BE WRAPPER for all east/west/south...
              v_name = trim(t_vname(itrc)) / / '_west'
              ierr=nf90_inq_varid(bry_id,v_name, t_west_id(itrc))  ! look for tracer variable in bry file
              if (ierr /= nf90_noerr) then                         ! if not there then trigger default values
#    ifdef DEFAULT_BRY_VALUES
                if(mynode==0) write(*,11) t_vname(itrc),
     &                        bry_file(ibry)(1:lsrc) MYID
                t_west_id(itrc) = -100 ! Set to -100 to avoid read_bry_EW/NS
                ierr=nf90_noerr
#    else
                if (ierr_all == nf90_noerr) ierr_all=ierr
                write(*,1) t_vname(itrc), bry_file(ibry)(1:lsrc)
#    endif
              endif
            enddo

#   endif
#  endif
          endif !<-- west_bry_active
# endif


# ifdef OBC_EAST
          if (east_bry_active) then
#  ifdef SOLVE3D
#   ifdef T_FRC_BRY

            do itrc=1,NT ! THIS CAN BE WRAPPER for all east/west/south...
              v_name = trim(t_vname(itrc)) / / '_east'
              ierr=nf90_inq_varid(bry_id,v_name, t_east_id(itrc))
              if (ierr /= nf90_noerr) then
#    ifdef DEFAULT_BRY_VALUES
                if(mynode==0) write(*,11) t_vname(itrc),
     &                        bry_file(ibry)(1:lsrc) MYID
                t_east_id(itrc) = -100 ! Set to -100 to avoid read_bry_EW/NS
                ierr=nf90_noerr
#    else
                if (ierr_all == nf90_noerr) ierr_all=ierr
                write(*,1) t_vname(itrc), bry_file(ibry)(1:lsrc)
#    endif
              endif
            enddo
#   endif
#  endif
          endif  !<-- east_bry_active
# endif


# ifdef OBC_SOUTH
          if (south_bry_active) then
#  ifdef SOLVE3D
#   ifdef T_FRC_BRY

            do itrc=1,NT ! THIS CAN BE WRAPPER for all east/west/south...
              v_name = trim(t_vname(itrc)) / / '_south'
              ierr=nf90_inq_varid(bry_id,v_name, t_south_id(itrc))
              if (ierr /= nf90_noerr) then
#    ifdef DEFAULT_BRY_VALUES
                if(mynode==0) write(*,11) t_vname(itrc),
     &                        bry_file(ibry)(1:lsrc) MYID
                t_south_id(itrc) = -100 ! Set to -100 to avoid read_bry_EW/NS
                ierr=nf90_noerr
#    else
                if (ierr_all == nf90_noerr) ierr_all=ierr
                write(*,1) t_vname(itrc), bry_file(ibry)(1:lsrc)
#    endif
              endif
            enddo
#   endif
#  endif
          endif  !<-- south_bry_active
# endif


# ifdef OBC_NORTH
          if (north_bry_active) then
#  ifdef SOLVE3D
#   ifdef T_FRC_BRY

            do itrc=1,NT ! THIS CAN BE WRAPPER for all east/west/south...
              v_name = trim(t_vname(itrc)) / / '_north'
              ierr=nf90_inq_varid(bry_id,v_name, t_north_id(itrc))
              if (ierr /= nf90_noerr) then
#    ifdef DEFAULT_BRY_VALUES
                if(mynode==0) write(*,11) t_vname(itrc),
     &                      bry_file(ibry)(1:lsrc) MYID
                t_north_id(itrc) = -100 ! Set to -100 to avoid read_bry_EW/NS
                ierr=nf90_noerr
#    else
                if (ierr_all == nf90_noerr) ierr_all=ierr
                write(*,1) t_vname(itrc), bry_file(ibry)(1:lsrc)
#    endif
              endif
            enddo
#   endif
#  endif
          endif !<-- north_bry_active
# endif
        endif
      endif   !<-- bry_rec < 0

  1   format(1x,'### ERROR: get_bry_all :: Cannot find variable ''',
     &                                A, ''' in file ''', A, '''.')
  11  format(1x,'--- WARNING: get_bry_all :: Cannot find variable ''',
     &                                A, ''' in file ''', A, '''.',A,i5)

! Read data from the file:  Check if model time is bounded by past and
!----- ---- ---- --- -----  future data times: if not, increment record
! and cycling indices, flip time index and read a new portion of data.
! Repeat until model time falls between the two data times.

! IMPORTANT! Always read "bry_time" even if no actual data needs to be
! read: this is necessary to make sure that running out of data records
! condition is the same for all MPI nodes.

      do while(bry_time(itbry) < time+0.5*dt .and. ierr == nf90_noerr)
        call advance_record(bry_rec, ntbry,  bry_ncycle, bry_cycle,
     &                                       ibry, max_bry,  ierr)
        if (bry_rec < 0) then
          ierr=nf90_close(bry_id) ; bry_id=-1
          ibry=ibry+1 ;  goto 100  !--> proceed to the next file
        elseif (ierr /= nf90_noerr) then
          if (ierr_all == nf90_noerr) ierr_all=ierr
        endif

        if (ierr == nf90_noerr) then
          ierr=nf_get_var1_double(bry_id, bry_time_id, bry_rec, dtime)
          if (ierr == nf90_noerr) then
            itbry=min(3-itbry,ntbry)
            bry_time(itbry)=dtime*day2sec + bry_cycle*bry_ncycle

            if (need_bry) then
# ifdef OBC_WEST
              if (west_bry_active) then
#  ifdef SOLVE3D
#   ifdef T_FRC_BRY
                do itrc=1,NT
#    ifdef DEFAULT_BRY_VALUES
                  if(t_west_id(itrc)>=0) then ! read boundary values
                    ierr=read_bry_EW(bry_id, t_west_id(itrc),  bry_rec,
     &                           r_var, t_west_dt(0,1,itbry,itrc), N)
                  else ! Use default values
                    call set_bry_tracer_defaults(t_west_dt,itrc,itbry,Mm+1)
                    if(mynode==0) print 3, trim(t_vname(itrc))
!                    ierr=nf90_noerr
                  endif
#    else
                  ierr=read_bry_EW(bry_id, t_west_id(itrc),  bry_rec,
     &                           r_var, t_west_dt(0,1,itbry,itrc), N)
#    endif /* DEFAULT_BRY_VALUES */
                  if (ierr_all == nf90_noerr) ierr_all=ierr
                enddo
#   endif
#  endif
              endif
# endif

# ifdef OBC_EAST
              if (east_bry_active) then
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
#    ifdef NHMG
                ierr=read_bry_EW(bry_id, w_east_id, bry_rec, r_var,  ! why repeated?
     &                                   w_east_dt(0,1,itbry), N)
                if (ierr_all == nf90_noerr) ierr_all=ierr
#    endif
#   endif
#   ifdef T_FRC_BRY
                do itrc=1,NT
#    ifdef DEFAULT_BRY_VALUES
                  if(t_east_id(itrc)>=0) then ! Read boundary values
                    ierr=read_bry_EW(bry_id, t_east_id(itrc),  bry_rec,
     &                           r_var, t_east_dt(0,1,itbry,itrc), N)
                  else ! Use default values
                    call set_bry_tracer_defaults(t_east_dt,itrc,itbry,Mm+1)
                    if(mynode==0) print 3, trim(t_vname(itrc))
!                    ierr=nf90_noerr
                  endif
#    else
                  ierr=read_bry_EW(bry_id, t_east_id(itrc),  bry_rec,
     &                           r_var, t_east_dt(0,1,itbry,itrc), N)
#    endif /* DEFAULT_BRY_VALUES */
                  if (ierr_all == nf90_noerr) ierr_all=ierr
                enddo
#   endif
#  endif
              endif
# endif

# ifdef OBC_SOUTH
              if (south_bry_active) then
#  ifdef SOLVE3D
#   ifdef T_FRC_BRY
                do itrc=1,NT
#    ifdef DEFAULT_BRY_VALUES
                  if(t_south_id(itrc)>=0) then ! Read boundary values
                    ierr=read_bry_NS(bry_id, t_south_id(itrc),  bry_rec,
     &                             r_var, t_south_dt(0,1,itbry,itrc), N)
                  else ! Use default values
                    call set_bry_tracer_defaults(t_south_dt,itrc,itbry,Lm+1)  ! only need to do once
                    if(mynode==0) print 3, trim(t_vname(itrc))
!                    ierr=nf90_noerr
                  endif
#    else
                  ierr=read_bry_NS(bry_id, t_south_id(itrc),  bry_rec,
     &                           r_var, t_south_dt(0,1,itbry,itrc), N)
#    endif /* DEFAULT_BRY_VALUES */
                  if (ierr_all == nf90_noerr) ierr_all=ierr
                enddo
#   endif
#  endif
              endif
# endif

# ifdef OBC_NORTH
              if (north_bry_active) then
#  ifdef SOLVE3D
#   ifdef T_FRC_BRY
                do itrc=1,NT
#    ifdef DEFAULT_BRY_VALUES
                  if(t_north_id(itrc)>=0) then ! Read boundary values
                    ierr=read_bry_NS(bry_id, t_north_id(itrc),  bry_rec,
     &                           r_var, t_north_dt(0,1,itbry,itrc), N)
                  else ! Use default values
                    call set_bry_tracer_defaults(t_north_dt,itrc,itbry,Lm+1)
                    if(mynode==0) print 3, trim(t_vname(itrc))
!                    ierr=nf90_noerr
                  endif
#    else
                  ierr=read_bry_NS(bry_id, t_north_id(itrc),  bry_rec,
     &                           r_var, t_north_dt(0,1,itbry,itrc), N)
#    endif /* DEFAULT_BRY_VALUES */
                  if (ierr_all == nf90_noerr) ierr_all=ierr
                enddo
#   endif
#  endif
              endif
# endif
              ierr=ierr_all
              if (ierr == nf90_noerr) then
!# ifdef MPI_SILENT_MODE
!                if (mynode == 0) then
!# endif
!                  write(*,'(9x,A,3x,A,F12.4,2(1x,A,I5))')
!     &           'get_bry_all :: read all boundary fields',
!     &           'bry_time =', dtime, 'rec =', bry_rec MYID
!# ifdef MPI_SILENT_MODE
!                endif
!# endif
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

 3    format(1x,'--- WARNING get_bry_all: Using default value for ',A)
      if (ierr_ext == nf90_noerr) ierr_ext=ierr_all
# ifdef VERBOSE
      write(*,*) 'leaving get_bry_all' MYID
# endif

      end subroutine get_bry_all  !]

!----------------------------------------------------------------------
      subroutine set_bry_all_tile (istr,iend,jstr,jend, ierr)  ![
      use param
      use scalars
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
          call set_frc_data( nc_z_w, zeta_west(j0:j1), obry=2 )
!          call exch_bry_EW_tile(1, ny, zeta_west, 1)      ! DevinD - exchange did not make a difference
#  endif
#  ifdef M2_FRC_BRY
          call set_frc_data( nc_ub_w, ubar_west(j0:j1), obry=1 )
          call set_frc_data( nc_vb_w, vbar_west( 1:j1), obry=1 )
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
          call set_frc_data( nc_u_w, var2d=u_west(j0:j1,:), obry=1 )
          call set_frc_data( nc_v_w, var2d=v_west( 1:j1,:), obry=1 )
#     ifdef NHMG
          call set_frc_data( nc_w_w, var2d=w_west(j0:j1,:), obry=1 )
#     endif
#   endif
#   ifdef T_FRC_BRY
          do itrc=1,NT
            do k=1,N
              do j=jmin,jmax
                t_west(j,k,itrc)=cff1*t_west_dt(j,k,it1,itrc)
     &                          +cff2*t_west_dt(j,k,it2,itrc)
              enddo
            enddo
!            call set_frc_data( nc_t_w(itrc), var2d=t_west(j0:j1,:,nt) )
          enddo
#   endif
#  endif
        endif
# endif

# ifdef OBC_EAST
        if (EASTERN_EDGE) then
#  ifdef Z_FRC_BRY
          call set_frc_data( nc_z_e, zeta_east(j0:j1), obry=2 )
#  endif
#  ifdef M2_FRC_BRY
          call set_frc_data( nc_ub_e, ubar_east(j0:j1), obry=1 )
          call set_frc_data( nc_vb_e, vbar_east( 1:j1), obry=1 )
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
          call set_frc_data( nc_u_e, var2d=u_east(j0:j1,:), obry=1 )
          call set_frc_data( nc_v_e, var2d=v_east( 1:j1,:), obry=1 )
#     ifdef NHMG
          call set_frc_data( nc_w_e, var2d=w_east(j0:j1,:), obry=1 )
#     endif
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
          call set_frc_data( nc_z_s, zeta_south(i0:i1), obry=2 )
#  endif
#  ifdef M2_FRC_BRY
          call set_frc_data( nc_ub_s, ubar_south( 1:i1), obry=1 )
          call set_frc_data( nc_vb_s, vbar_south(i0:i1), obry=1 )
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
          call set_frc_data( nc_u_s, var2d=u_south( 1:i1,:), obry=1 )
          call set_frc_data( nc_v_s, var2d=v_south(i0:i1,:), obry=1 )
#     ifdef NHMG
          call set_frc_data( nc_w_s, var2d=w_south(i0:i1,:), obry=1 )
#     endif
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
          call set_frc_data( nc_z_n, zeta_north(i0:i1), obry=2 )
#  endif
#  ifdef M2_FRC_BRY
          call set_frc_data( nc_ub_n, ubar_north( 1:i1), obry=1 )
          call set_frc_data( nc_vb_n, vbar_north(i0:i1), obry=1 )
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
          call set_frc_data( nc_u_n, var2d=u_north( 1:i1,:), obry=1 )
          call set_frc_data( nc_v_n, var2d=v_north(i0:i1,:), obry=1 )
#     ifdef NHMG
          call set_frc_data( nc_w_n, var2d=w_north(i0:i1,:), obry=1 )
#     endif
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
           write(*,'(/3A/3(1x,A,F16.10)/)')            '### ERROR: ',
     &      'set_bry_all_tile :: current model time is out of bounds ',
     &      'of ''bry_time''.', 'bry_tstart =', bry_time(it1)*sec2day,
     &      'tdays =',  tdays,  'bry_tend =',   bry_time(it2)*sec2day
          ierr=ierr+1
          error stop
        endif
      endif  ! <-- need_bry

      end subroutine set_bry_all_tile  !]
#endif

!----------------------------------------------------------------------

      end module boundary
