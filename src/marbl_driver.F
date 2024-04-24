      module marbl_driver

!-----------------------------------------------------------------------
!     MODULE: marbl_driver
!      
!     DESCRIPTION:
!     This module interfaces between MARBL and ROMS.
!     It configures a MARBL "instance" with a particular number of tracers
!     and diagnostics, and passes this information back to ROMS for i/o.
!     It then (at each time step and gridpoint) receives tracer/forcing
!     information from ROMS and calculates interior tendencies of tracers
!     and surface fluxes. ROMS then updates its tracer and diagnostics
!     arrays based on these calculations.
!      
!     NOTES:
!     Configured to run with MARBL 0.45 on the "development" branch:
!     https://github.com/marbl-ecosys/MARBL/releases/tag/marbl0.45.0
!      
!     AUTHOR: Dafydd Stephenson
!     DATE: 2024-03-06
!-----------------------------------------------------------------------
            
#include "cppdefs.opt"
!     MARBL cpp def covering entire module, module empty if off
#ifdef MARBL

!     ----------------------------------------------------------------------

      use marbl_interface             , only: marbl_interface_class
#ifdef MARBL_DIAGS
      use MARBL_interface_public_types, only: marbl_diagnostics_type,
     &                                        marbl_saved_state_type
#endif            
      use param                       , only: mynode,nt,ntrc_bio,
     &                                        isalt,itemp,Lm,Mm
      use dimensions                  , only: nz
      use surf_flux                   , only: uwnd, vwnd
      use bgc_forces                  , only: pco2air,pco2air_alt,
     &                                        dust,iron
      use surf_flux                   , only: srflx
      use nc_read_write
      use scalars                     , only: nstp,nnew,dt,iic,ntstart
      use grid                        , only : rmask
#ifdef MPI
      use mpi
      use hidden_mpi_vars
#endif      
#ifdef NOX_FORCING
      use bgc_forces                  , only: nox
#endif
#ifdef NHY_FORCING
      use bgc_forces                  , only: nhy
#endif      
      use ocean_vars                  , only: Hz, z_r, z_w
      
!     ----------------------------------------------------------------------
      
      implicit none

      type(MARBL_interface_class) :: marbl_instance ! Object to interface between ROMS and MARBL
      integer                     :: nt_marbl       ! number of marbl tracers (according to MARBL)
      integer                     :: idx            ! Looping variable
       
!     Surface flux forcing variable indices (-1 means unused,default)
      integer,public ::    u10_sqr_ind=-1,       sss_ind=-1,      sst_ind=-1,
     &                       ifrac_ind=-1,  dust_dep_ind=-1,   fe_dep_ind=-1,
     &                    nox_flux_ind=-1,  nhy_flux_ind=-1, atmpress_ind=-1,
     &                        xco2_ind=-1,  xco2_alt_ind=-1

!     Interior tendency forcing variable indices (-1 means unused,default)
      integer,public ::       dustflux_ind=-1, PAR_col_frac_ind=-1,
     &                  surf_shortwave_ind=-1,       potemp_ind=-1,
     &                        salinity_ind=-1,     pressure_ind=-1,
     &                       fesedflux_ind=-1,    o2_scalef_ind=-1,
     &                    remin_scalef_ind=-1


#ifdef MARBL_DIAGS      
!     Keeping count of diagnostics:
!     These are the amount to allocate in ROMS...
      integer, parameter, public :: nr_marbl_diag_2d = 179
      integer, parameter, public :: nr_marbl_diag_3d = 164
!     ... and these are the amount supplied by MARBL (to be counted up)
      integer, public :: diag_cnt_sf_2d  ! Number of MARBL 2D surface flux diagnostic vars
      integer, public :: diag_cnt_sf_3d  ! Number of MARBL 3D surface flux diagnostic vars
      integer, public :: diag_cnt_it_2d  ! Number of MARBL 2D interior tendency diagnostic vars
      integer, public :: diag_cnt_it_3d  ! Number of MARBL 3D interior tendency diagnostic vars
      integer, public :: diag_cnt_2d_tot ! Total number of MARBL 2D diagnostic vars
      integer, public :: diag_cnt_3d_tot ! Total number of MARBL 3D diagnostic vars      
      integer, public :: diag_cnt_tot    ! Total number of MARBL diagnostic vars

!     Making local versions of these "bec2" arrays for bgc diagnostic output control
      logical, dimension(:), pointer :: wrt_bgc_diag_2d
      logical, dimension(:), pointer :: wrt_bgc_diag_3d
      integer, dimension(:), pointer :: idx_bgc_diag_2d
      integer, dimension(:), pointer :: idx_bgc_diag_3d      

      real,dimension(:,:,:,:),pointer :: diag_array_3d ! points to bec2_diag_3d
      real,dimension(:,:,:)  ,pointer :: diag_array_2d ! points to bec2_diag_2d
      
#endif      
!     Saved state variables: (mimicking structure of bgc_ecosys_vars for bec2_diag arrays)
      integer                             :: nr_marbl_ss_2d,nr_marbl_ss_3d       ! No. of saved_state vars
      parameter( nr_marbl_ss_2d=3, nr_marbl_ss_3d=2)                             ! Hardcode no. of vars (error if wrong)
      
      real,allocatable,dimension(:,:,:  ) :: marbl_saved_state_2d                ! Array to store them
      real,allocatable,dimension(:,:,:,:) :: marbl_saved_state_3d                ! " "
      character*72                        :: vname_marbl_ss_2d(4,nr_marbl_ss_2d) ! Metadata (shortname,longname,units)
      character*72                        :: vname_marbl_ss_3d(4,nr_marbl_ss_3d) ! " "


!     Integers to define where print statements should come from
      integer                             :: printi = 4    ! i location
      integer                             :: printj = 4    ! j location
      integer                             :: printnode = 7 ! cpu


! ----------------------------------------------------------------------
      
      contains
      
      subroutine marbldrv_configure_tracers(
     &     itot,t_vname,t_lname,t_units,t_tname,wrt_t,wrt_t_avg,t_ana_frc)
!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_configure_tracers
!     DESCRIPTION:
!     Set up MARBL instance and populate ROMS tracer metadata arrays
!      
!     METHOD:
!     1. Read a settings file (marbl_in) and pass the settings to the
!     uninitialised MARBL instance created above
!     2. Initialise the MARBL instance, pass metadata on expected tracers
!     from the MARBL instance back to ROMS
!     3. Check which surface flux and interior tendency forcing variables
!     are required by MARBL later in this module, keep/print this info
!                 
!     INPUT/OUTPUT:
!     - itot:
!          tracer/metadata array index variable
!     - t_vname,t_lname,t_units:
!          tracer short,long name,units metadata arrays
!     - wrt_t,wrt_t_avg:
!          logicals to determine if/how each tracer is written to output
!     - t_ana_frc:
!          integer to determine if each tracer's surface flux is
!          read in (0), or analytical (1,default)
!          
!     NOTES:
!     - MARBL requires domain geometry for initialisation, but the
!     associated ROMS variables have not yet been set at the time this
!     SR is called. A dummy array is instead passed to MARBL and the
!     geometry is updated with the correct arrays in a later SR before
!     any BGC calculations take place.
!       
!     
!-----------------------------------------------------------------------
      
      
      character*42, dimension(nt), intent(inout) :: t_vname, t_units
      character*60, dimension(nt), intent(inout) :: t_lname
      character*47, dimension(nt), intent(inout) :: t_tname
      character(len=256)                         :: namelist_line
      integer                    , intent(inout) :: itot
      integer , dimension(nt)    , intent(inout) :: t_ana_frc
      logical , dimension(nt)    , intent(inout) :: wrt_t, wrt_t_avg
      real    , dimension(nz)                    :: dummy_array
      integer :: marbl_settings_in, open_status, read_status, ierr
      
      dummy_array(:)=1.0        ! Temporary measure until we can access an allocated Hz,z_w,z_r
      marbl_settings_in = 111

!     1. Read in MARBL settings file (marbl_in), pass settings to MARBL:
!     ------------------------------------------------------------------
      

!     read the marbl_in into buffer
      if (mynode==0) then
         open(unit=marbl_settings_in, file='marbl_in', iostat=open_status)
         if (open_status .ne. 0) then ! can't find namelist file, skip
            write(*, '(A, I0, 3A)') "IO WARNING ", open_status,
     &           "could not open namelist file : ", "marbl_in",
     &           ". attempting to proceed with default MARBL settings"

         end if
      end if         
#ifdef MPI
         call MPI_Bcast(open_status,1,MPI_INTEGER,0,ocean_grid_comm,ierr)
#endif            
         if (open_status==0) then
            namelist_line = ''
            do                  !namelist line-by-line read
               if (mynode==0) then
                  read(marbl_settings_in,"(A)",iostat=read_status),
     &                 namelist_line
               end if
#ifdef MPI
               call MPI_Bcast(read_status,1,MPI_INTEGER,0,ocean_grid_comm,ierr)
#endif            
               
               if (read_status .ne. 0) then
                  exit          !loop
               else
#ifdef MPI
                  call MPI_Bcast(namelist_line,256,MPI_CHARACTER,0,ocean_grid_comm,ierr)
#endif                  
                  call marbl_instance%put_setting(namelist_line)
                  
                  if (mynode==printnode) print *,'marbl_in: ',namelist_line
               end if
               
            end do                     !/namelist line-by-line read

            if (is_iostat_end(read_status)) then
               if (mynode==printnode) then
                  write(*, '(3A)') "Successfully read ","marbl_in",
     &                 " to end of file"
               end if
            else
               if (mynode==printnode) then
                  write(*, '(A, I0, 2A)') "IO ERROR ", read_status,
     &                 " : ERROR reading MARBL namelist file  ",
     &                 "marbl_in"
               end if
               error stop
            end if
         
            close(marbl_settings_in)
            if (mynode==printnode) then
               print *, 'we successfully read and closed marbl_in'
            end if
         endif

!     2. Now generate a MARBL instance with these settings:
! ---------------------------------------------------------

      ! Don't yet have values for domain geometry, so init with dummy & update later:
      call marbl_instance%init(gcm_num_levels=nz,
     &     gcm_num_PAR_subcols = 1,
     &     gcm_num_elements_surface_flux = 1,
     &     gcm_delta_z = dummy_array(:), ! to update with Hz later
     &     gcm_zw = dummy_array(:),      ! to update with z_w later
     &     gcm_zt = dummy_array(:),      ! to update with z_r later
     &     lgcm_has_global_ops = .true.,
     &     unit_system_opt='mks')

!     Number of MARBL tracers according to MARBL
      nt_marbl=size(marbl_instance%tracer_metadata)
      
!     Check that the number of tracers in MARBL_instance agrees with ROMS' ntrc_bio:
      if ( nt_marbl .ne. ntrc_bio ) then
         write(*,'(7x,A,I4,A,I4,A,I4,A)')
     &        'ERROR: Allocated no. of MARBL tracers'
     &        ,ntrc_bio
     &        ,' does not match no. expected by MARBL: '
     &        ,nt_marbl
     &        ,' set ntrc_bio = ',nt_marbl
     &        ,' in param.opt and recompile'
         error stop 'ERROR reported from MARBL library'
      end if      
        
      call print_marbl_log(marbl_instance%StatusLog)
      call marbl_instance%StatusLog%erase()

!     Pass tracer metadata from new MARBL instance to ROMS for i/o:
      do idx=1,nt_marbl ! Loop over MARBL tracers...
         itot=itot+1    ! ...but use itot as index (to include other tracers)     
         t_vname(itot)=trim(marbl_instance%tracer_metadata(idx)%short_name)
         t_lname(itot)=trim(marbl_instance%tracer_metadata(idx)%long_name)
         t_units(itot)=trim(marbl_instance%tracer_metadata(idx)%units)
         
         t_tname(itot)=''
         wrt_t     (itot)=.false.
         wrt_t_avg (itot)=.true.
         t_ana_frc( itot)=1

      end do
      ! Read user file to determine which tracers to write out
      call marbldrv_read_tracer_output_list(t_vname,wrt_t)
      
!     3. Determine surf flux forcing variables MARBL needs
!     ---------------------------------------------------------
      
      if (mynode==printnode) then 
         print *, 'Here are the SF forcings MARBL requested'
      end if
      do idx=1,size(marbl_instance%surface_flux_forcings)
         if (mynode==printnode) then
            print *, 'var: ', (trim(MARBL_instance%surface_flux_forcings(idx)%metadata%varname))
            print *, 'units:',(trim(MARBL_instance%surface_flux_forcings(idx)%metadata%field_units))
         end if

!     Each SF variable has an index, set to -1 if not needed (default),
!     or set to the position in the MARBL instance if needed,
!     for later use in marbldrv_column_physics

         select case (trim(MARBL_instance%surface_flux_forcings(idx)%metadata%varname))
         case('u10_sqr')
            u10_sqr_ind = idx
         case('sss')
            sss_ind = idx
         case('sst')
            sst_ind = idx
         case('Ice Fraction')
            ifrac_ind = idx
         case('Dust Flux')
            dust_dep_ind = idx
         case('Iron Flux')
            fe_dep_ind = idx
         case('NOx Flux')
            nox_flux_ind = idx
         case('NHy Flux')
            nhy_flux_ind = idx
         case('Atmospheric Pressure')
            atmpress_ind = idx
         case('xco2')
            xco2_ind = idx
         case('xco2_alt_co2')
            xco2_alt_ind = idx
         case DEFAULT
            print *, 'Additional forcing requested but not indexed: ',
     &           trim(MARBL_instance%surface_flux_forcings(idx)
     &           %metadata%varname)
         end select
      enddo

!     4. Determine interior tendency forcing vars MARBL needs
!     ---------------------------------------------------------
      
      if (mynode==printnode) then
         print *, 'Here are the int. tend. forcings MARBL requested'
      end if

!     Each IT variable has an index, set to -1 if not needed (default),
!     or set to the position in the MARBL instance if needed,
!     for later use in marbldrv_column_physics
      
      do idx=1,size(MARBL_instance%interior_tendency_forcings)
         if (mynode==printnode) then
            print *,'var: '  ,(trim(MARBL_instance%interior_tendency_forcings(idx)%metadata%varname))
            print *,'units: ',(trim(MARBL_instance%interior_tendency_forcings(idx)%metadata%field_units))            
         end if
   
         select case (trim(MARBL_instance%interior_tendency_forcings(idx)%metadata%varname))
         case('Dust Flux')
            dustflux_ind = idx
         case('PAR Column Fraction')
            PAR_col_frac_ind = idx
         case('Surface Shortwave')
            surf_shortwave_ind = idx
         case('Potential Temperature')
            potemp_ind = idx
         case('Salinity')
            salinity_ind = idx
         case('Pressure')
            pressure_ind = idx
         case('Iron Sediment Flux')
            fesedflux_ind = idx
         case('O2 Consumption Scale Factor')
            o2_scalef_ind = idx
         case('Particulate Remin Scale Factor')
            remin_scalef_ind = idx
         case DEFAULT
            print *, 'Additional forcing requested but not indexed: ',
     &           trim(MARBL_instance%interior_tendency_forcings(idx)
     &           %metadata%varname)
         end select
      enddo

      
      end subroutine marbldrv_configure_tracers

      subroutine marbldrv_read_tracer_output_list(t_vname,wrt_t)

!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_read_tracer_output_list
!     DESCRIPTION: Read "marbl_tracer_output_list" file and set wrt_t
!                  logical array entries accordingly
!
!     METHOD:
!     1. Read the "marbl_tracer_output_list" file line-by-line
!     2. For each line not beginning with "!", compare that line with
!         every entry in "t_vname" (tracer shortname list)
!     3. If there is a match, set "wrt_t" to true for that tracer
!     4. If there are no matches, default to writing out  all tracers
!
!     NOTES:
!     - The possible options are listed in STDOUT or can be found with
!      $MARBL_ROOT/tests/regression_tests/requested_tracers/requested_tracers.py
!-----------------------------------------------------------------------

      
      character*42, dimension(nt), intent(inout) :: t_vname
      character(len=256)                         :: file_line
      logical , dimension(nt)    , intent(inout) :: wrt_t
      
      integer :: marbl_tracer_list, open_status, read_status,idx,ierr
      logical :: recognized_tracer, no_tracers_requested

      marbl_tracer_list = 1111
      no_tracers_requested = .true. ! If we find any requested tracers, flip to false

!     1. Read in MARBL tracer list file "marbl_tracer_output_list"
!     ------------------------------------------------------------------
      
      if (mynode==0) then
         open(unit=marbl_tracer_list, file='marbl_tracer_output_list',
     &        iostat=open_status)
         
         if (open_status .ne. 0) then ! can't find namelist file, skip and output all tracers
            write(*, '(A, I0, 3A)') "IO WARNING ", read_status,
     &           "could not open tracer list file : ",
     &           "marbl_tracer_output_list",
     &           ". All MARBL tracers will be written to output"
         end if
      end if
#ifdef MPI
      call MPI_Bcast(open_status,1,MPI_INTEGER,0,ocean_grid_comm,ierr)      
#endif      
      if (open_status==0) then
         file_line = ''
         do                     !line-by-line read
            if (mynode==0) then
               read(marbl_tracer_list,"(A)",iostat=read_status),
     &              file_line
            end if
#ifdef MPI
            call MPI_Bcast(read_status,1,MPI_INTEGER,0,ocean_grid_comm,ierr)            
#endif            
            if (read_status .ne. 0) then
               exit             !loop
            else
#ifdef MPI
               call MPI_Bcast(file_line,256,MPI_CHARACTER,0,ocean_grid_comm,ierr)
#endif                
!     Check if the line is a comment
               file_line=trim(adjustl(file_line))
               if (file_line(1:1) == '!') cycle ! Skip this iteration if it's a comment
!     Set the "wrt_t" variable
               recognized_tracer=.false.
               do idx=1,NT
                  if (trim(t_vname(idx)) == trim(file_line)) then
                     wrt_t(idx) =.true.
                     if (mynode==printnode) then
                        write(*,'(2A)'),
     &                       'REQUESTED MARBL TRACER OUTPUT: ' // trim(t_vname(idx))
                     end if 
                     recognized_tracer=.true.
                     no_tracers_requested=.false.
                     cycle
                  end if
               end do           ! idx
               if ( (.not. recognized_tracer) .and. (mynode==printnode))
     &              then
                  write(*,'(A)'), 'REQUESTED TRACER ' // trim(file_line) //
     &                 ' NOT RECOGNIZED. See tracer metadata section' //
     &                 ' of this output log for ' //
     &                 ' a list of valid entries.'
               end if
            end if              ! read_status=/=0
               
         end do                 !/namelist line-by-line read

         if (no_tracers_requested) then ! Write all of them to output
               
            wrt_t(:)=.true.
               
            if (mynode==0) then
               write(*,'(A)'), 'NO REQUESTED TRACERS FOUND IN ' //
     &              ' marbl_tracer_output_list. Defaulting to ' //
     &              ' writing all BGC tracers to output file.'
            end if
         end if
            
         if (is_iostat_end(read_status)) then
            if (mynode==printnode) then
               write(*, '(3A)') "Successfully read ",
     &              "marbl_tracer_output_list",
     &              " to end of file"
            end if
         else
            if (mynode==printnode) then
               write(*, '(A, I0, 2A)') "IO ERROR ", read_status,
     &              " : ERROR reading MARBL tracer list file  ",
     &              "marbl_tracer_output_list"
            end if
            error stop
         end if
         
         close(marbl_tracer_list)
         if (mynode==printnode) then
            write(*,'(2A)'), 'we successfully read and closed ',
     &           'marbl_tracer_output_list'
         end if
      endif

      end subroutine marbldrv_read_tracer_output_list

      
      subroutine marbldrv_configure_saved_state
      
!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_configure_saved_state
!     DESCRIPTION: Populate metadata arrays for MARBL saved state vars
!                  and allocate 2d and 3d saved state variable arrays
!
!     METHOD:
!     1. Verify ROMS' hardcoded number of saved state variables
!        (nr_marbl_ss_2d,nr_marbl_ss_3d in the module header)
!        corresponds to the number in the MARBL instance
!     2. Populate ROMS' metadata arrays for the MARBL saved state
!        based on information provided by the MARBL instance
!        (vname_marbl_ss_2d,vname_marbl_ss_3d in the module header)
!     3. Allocate arrays for storing MARBL saved state variables
!        (marbl_saved_state_2d, marbl_saved_state_3d)
!
!     NOTES:
!     - It is assumed that 3D surface flux saved state variables and
!       2D interior tendency saved state variables are possible for
!       purposes of generality, however, no such variables exist in
!       current versions of MARBL
!     
!-----------------------------------------------------------------------
      
      integer            :: nr_marbl_ss_2d_check,nr_marbl_ss_3d_check,
     &                      iss,i,j,k
      character(len=200) :: ss_varname

!     1. Check hardcoded array size matches MARBL expectations
!     ----------------------------------------------------------------------
      
!     Initialise "check" values to 0 and add 1 for each variable MARBL needs
      nr_marbl_ss_2d_check=0
      nr_marbl_ss_3d_check=0
      
!     Go through both MARBL saved state types and check if fields are
!     2d or 3d one by one (calculate nr_marbl_ss_2d,nr_marbl_ss_3d), then
!     increment the size check variables accordingly:

!     Add surface flux saved state variables to the counts:
            do iss=1,MARBL_instance%surface_flux_saved_state%saved_state_cnt
               if (MARBL_instance%surface_flux_saved_state%state(iss)%rank==2)
     &              then
                  nr_marbl_ss_2d_check=nr_marbl_ss_2d_check+1
               elseif (MARBL_instance%surface_flux_saved_state%state(iss)%rank==3)
     &                 then
                  nr_marbl_ss_3d_check=nr_marbl_ss_3d_check+1
               end if
            end do

!     Add interior tendency saved state variables to the counts:
            do iss=1,MARBL_instance%interior_tendency_saved_state%saved_state_cnt
               if (MARBL_instance%interior_tendency_saved_state%state(iss)%rank==2)
     &              then
                  nr_marbl_ss_2d_check=nr_marbl_ss_2d_check+1
               elseif (MARBL_instance%interior_tendency_saved_state%state(iss)%rank==3)
     &                 then
                  nr_marbl_ss_3d_check=nr_marbl_ss_3d_check+1
               end if
            end do

!     Throw error if counts don't match            
            if (nr_marbl_ss_2d/=nr_marbl_ss_2d_check) then
               write(*,'(7x,A,I4,A,I4,A)')
     &              'ERROR: Allocated no. of 2D marbl saved state vars '
     &              ,nr_marbl_ss_2d,
     &              ' does not match actual no. expected by MARBL: '
     &              ,nr_marbl_ss_2d_check
               error stop 'ERROR reported from MARBL library'
            elseif (nr_marbl_ss_3d/=nr_marbl_ss_3d_check) then
               write(*,'(7x,A,I4,A,I4,A)')
     &              'ERROR: Allocated no. of 3D marbl saved state vars '               
     &              ,nr_marbl_ss_3d,
     &              ' does not match actual no. expected by MARBL: '
     &              ,nr_marbl_ss_3d_check
               error stop 'ERROR reported from MARBL library'
            end if
               
!     2. Populate saved state metadata arrays in ROMS using MARBL instance
!     ----------------------------------------------------------------------
            
            ! surface flux saved state metadata (all 2D variables currently)
            do iss=1,MARBL_instance%surface_flux_saved_state%saved_state_cnt
               write(ss_varname, "(2A)")
     &              "MARBL_",
     &              trim(MARBL_instance%surface_flux_saved_state%state(iss)%short_name)
               vname_marbl_ss_2d(1,iss)=trim(ss_varname)
               vname_marbl_ss_2d(2,iss)=
     &              trim(MARBL_instance%surface_flux_saved_state%state(iss)%long_name)
               vname_marbl_ss_2d(3,iss)=
     &              trim(MARBL_instance%surface_flux_saved_state%state(iss)%units)
               vname_marbl_ss_2d(4,iss)='  '
               
               if (mynode==printnode) then
                  print *, 'MARBL SURFACE FLUX SAVED STATE VARIABLE ',
     &                 trim(vname_marbl_ss_2d(1,iss)),'  ',
     &                 trim(vname_marbl_ss_2d(2,iss)),'  ',
     &                 trim(vname_marbl_ss_2d(3,iss))
               end if
               
            end do
            
            ! interior tendency saved state metadata (all 3D variables currently)
            do iss=1,MARBL_instance%interior_tendency_saved_state%saved_state_cnt
               write(ss_varname, "(2A)")
     &              "MARBL_",
     &              trim(MARBL_instance%interior_tendency_saved_state%state(iss)%short_name)
               vname_marbl_ss_3d(1,iss)=trim(ss_varname)
               vname_marbl_ss_3d(2,iss)=
     &              trim(MARBL_instance%interior_tendency_saved_state%state(iss)%long_name)
               vname_marbl_ss_3d(3,iss)=
     &              trim(MARBL_instance%interior_tendency_saved_state%state(iss)%units)
               

               if (mynode==printnode) then
                  print *,
     &                 'MARBL INTERIOR TENDENCY SAVED STATE VARIABLE ',
     &                 trim(vname_marbl_ss_3d(1,iss)),'  ',
     &                 trim(vname_marbl_ss_3d(2,iss)),'  ',
     &                 trim(vname_marbl_ss_3d(3,iss))  
               end if
                  
            end do

!     3. Allocate saved state variable arrays in ROMS
!     ----------------------------------------------------------------------
            
            allocate( marbl_saved_state_2d(GLOBAL_2D_ARRAY,nr_marbl_ss_2d))
            allocate( marbl_saved_state_3d(GLOBAL_2D_ARRAY,nz,nr_marbl_ss_3d))
            if (mynode==printnode) then
               write(*,'(7x,A,I4,I4)') 'MARBL saved state allocation ::',
     &              nr_marbl_ss_2d,nr_marbl_ss_3d !3,2
            end if
            
            
      end subroutine marbldrv_configure_saved_state
!--------------------------------------------------------------------------------
      subroutine marbldrv_read_ss_vars_from_rst(
     &     ncid,rec_rst)

!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_read_saved_state_from_restart_file
!     DESCRIPTION:
!     Get MARBL saved state information from IC or restart file
!      
!     METHOD:
!     loop over 2D and 3D saved state metadata arrays populated in
!     marbldrv_configure_saved_state and get associated variables
!     using ncread. Called from inside `get_init` module.
!                 
!     INPUT:
!     - ncid: integer identifier of the restart netCDF file
!     - rec_rst: current file record number
      
!     NOTES:
!     - These code blocks adapted from get_init 
!     
!-----------------------------------------------------------------------
      
      use nc_read_write, only: ncread,handle_ierr
      use netcdf
      use dimensions, only : i0,i1,j0,j1
      use mpi_exchanges, only : exchange_xxx
!      use hidden_mpi_vars
      use param, only : N
      implicit none

      integer, intent(in) :: ncid,rec_rst
      integer :: ierr,varid,start(4)
      integer :: itrc,k
      
!     Read in MARBL saved_state variables here      
!     vname indices: 1=short name, 2=long name, 3=units
      start=1; start(3)=rec_rst ! Set "start" index to 1 for all dimensions except time (3) which starts at current record
!     Read in 2D saved state variables (for surface_flux_saved_state)
      do itrc=1,nr_marbl_ss_2d
         ierr=nf90_inq_varid (ncid, vname_marbl_ss_2d(1,itrc), varid) 
         if (ierr .eq. nf90_noerr) then
            call ncread(ncid,vname_marbl_ss_2d(1,itrc),
     &           marbl_saved_state_2d(i0:i1,j0:j1,itrc),start)
            
            marbl_saved_state_2d(i0:i1,j0:j1,itrc)=
     &         marbl_saved_state_2d(i0:i1,j0:j1,itrc)*rmask(i0:i1,j0:j1)
            call exchange_xxx(marbl_saved_state_2d(:,:,itrc))
         else
            marbl_saved_state_2d(GLOBAL_2D_ARRAY,itrc) = 0.0
            if(mynode==0) write(*,*) ' --- WARNING: '
     &         , trim(vname_marbl_ss_2d(1,itrc))            
     &         , ' not in initial file.  Initialized to 0.0'
            ierr=nf90_noerr
            
         endif
      end do
      start=1; start(4)=rec_rst                                       ! 3D vars

!     Read in 3d saved state variables (for interior_tendency_saved_state)
      do itrc=1,nr_marbl_ss_3d
         ierr=nf90_inq_varid (ncid, vname_marbl_ss_3d(1,itrc), varid)
         if (ierr == nf90_noerr) then
            call ncread(ncid, vname_marbl_ss_3d(1,itrc),
     &           marbl_saved_state_3d(i0:i1,j0:j1,:,itrc),start)
           do k=1,N
              marbl_saved_state_3d(i0:i1,j0:j1,k,itrc)=
     &             marbl_saved_state_3d(i0:i1,j0:j1,k,itrc)*rmask(i0:i1,j0:j1)
           enddo
           call exchange_xxx(marbl_saved_state_3d(:,:,:,itrc) )
         else
            marbl_saved_state_3d(GLOBAL_2D_ARRAY,:,itrc) = 0.0 ! Note ':' instead of '1' as in similar blocks
            if(mynode==0) write(*,*) ' --- WARNING: '
     &         , trim(vname_marbl_ss_3d(1,itrc))            
     &         , ' not in initial file.  Initialized to 0.0'
            ierr=nf90_noerr            
         endif
      end do
      
      
      end subroutine marbldrv_read_ss_vars_from_rst

      subroutine marbldrv_create_ss_vars_in_rst(ncid)
!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_create_ss_vars_in_rst
!     DESCRIPTION:
!     Create netcdf variables corresponding to MARBL saved state data
!     in ROMS' restart file
!      
!     METHOD:
!     loop over 2D and 3D saved state metadata arrays populated in
!     marbldrv_configure_saved_state create netcdf variables using nccreate
!     created in marbldrv_create_ss_vars_in_rst.
!     Called from inside `get_init` module.
!                 
!     INPUT:
!     - ncid: integer identifier of the restart netCDF file
!     - rec_rst: current file record number
      
!     NOTES:
!     - These code blocks adapted from get_init
!     - RELATED: marbldrv_read_ss_vars_from_rst,
!                marbldrv_write_ss_vars_to_rst
!     
!-----------------------------------------------------------------------

      use netcdf
      use nc_read_write
      use roms_read_write, only: dn_xr,dn_yr,dn_zr,dn_tm,xi_rho,eta_rho
      use param, only : N
      
      implicit none      

      integer, intent(in) :: ncid
      integer itrc,varid,ierr
      
      !     Add MARBL saved_state variables to restart file
      do itrc=1,nr_marbl_ss_3d
         varid = nccreate(ncid,vname_marbl_ss_3d(1,itrc),
     &                      (/dn_xr,dn_yr,dn_zr,dn_tm/),
     &                      (/xi_rho,eta_rho,N,0/), nf90_double)
        ierr=nf90_put_att (ncid, varid, 'long_name', vname_marbl_ss_3d(2,itrc))
        ierr=nf90_put_att (ncid, varid, 'units', vname_marbl_ss_3d(3,itrc))
      enddo
      do itrc=1,nr_marbl_ss_2d
         varid = nccreate(ncid,vname_marbl_ss_2d(1,itrc),
     &                      (/dn_xr,dn_yr,dn_tm/),
     &                      (/xi_rho,eta_rho,0/), nf90_double)
        ierr=nf90_put_att (ncid, varid, 'long_name', vname_marbl_ss_2d(2,itrc))
        ierr=nf90_put_att (ncid, varid, 'units', vname_marbl_ss_2d(3,itrc))
      enddo

      end subroutine
      
      subroutine marbldrv_write_ss_vars_to_rst(      
     &     ncid,rec_rst)

!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_write_ss_vars_to_rst
!     DESCRIPTION:
!     Add MARBL saved state data to ROMS' restart file
!      
!     METHOD:
!     loop over 2D and 3D saved state metadata arrays populated in
!     marbldrv_configure_saved_state and write to netcdf variables
!     created in marbldrv_create_ss_vars_in_rst.
!     Called from inside `get_init` module.
!                 
!     INPUT:
!     - ncid: integer identifier of the restart netCDF file
!     - rec_rst: current file record number
      
!     NOTES:
!     - These code blocks adapted from get_init
!     - RELATED: marbldrv_read_ss_vars_from_rst
!     
!-----------------------------------------------------------------------
      
      use nc_read_write, only: ncwrite
      use dimensions, only : i0,i1,j0,j1
      
      implicit none

      integer, intent(in) :: ncid,rec_rst
      integer :: start(4)
      integer :: itrc
      
        ! Add MARBL saved_state variables to restart file
        start=1; start(4)=rec_rst
        do itrc=1,nr_marbl_ss_3d
           call ncwrite(ncid,vname_marbl_ss_3d(1,itrc),
     &          marbl_saved_state_3d(i0:i1,j0:j1,:,itrc),start)
        enddo        
        start=1;start(3)=rec_rst        ! back to 2D vars
        do itrc=1,nr_marbl_ss_2d
           call ncwrite(ncid,vname_marbl_ss_2d(1,itrc),
     &          marbl_saved_state_2d(i0:i1,j0:j1,itrc),start)
        enddo                
      
      
      end subroutine marbldrv_write_ss_vars_to_rst

      
#ifdef MARBL_DIAGS
      subroutine marbldrv_configure_diagnostics(
     &     vname_bgc_diag_2d,vname_bgc_diag_3d
!     &     ,wrt_bgc_diag_2d,wrt_bgc_diag_3d
     &     )

!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_configure_diagnostics
!     DESCRIPTION:
!     Populate ROMS diagnostic metadata arrays using information on
!     MARBL diagnostics provided by the MARBL instance
!      
!     METHOD:
!     1. Verify ROMS' hardcoded number of diagnostics variables
!        (nr_bec2_diag_2d,nr_bec2_diag_3d in bgc_ecosys_vars.F)
!        corresponds to the number of diags in the MARBL instance
!     2. Populate ROMS BGC diagnostics metadata arrays
!        (vname_bec2_diag_2d,vname_bec2_diag_3d) using information
!        provided by the MARBL instance
!                 
!     INPUT/OUTPUT:
!     - vname_bgc_diag_2d,vname_bgc_diag_3d: 
!        diagnostics metadata arrays (short,long names and units)
!     - wrt_bgc_diag_2d,wrt_bgc_diag_3d:
!         logical arrays to determine which diagnostics are output
!         (currently hardcoded to all)
!          
!     NOTES:
!     - For generalisation purposes, it is assumed that 3D surface flux
!       diagnostics and 2D interior tendency diagnostics are possible.
!       However, all current versions of MARBL do not contain any.
!     
!-----------------------------------------------------------------------

      character*72, dimension(:,:), intent(inout) :: vname_bgc_diag_2d
      character*72, dimension(:,:), intent(inout) :: vname_bgc_diag_3d
!      logical,      dimension(:),   intent(inout) :: wrt_bgc_diag_2d
!      logical,      dimension(:),   intent(inout) :: wrt_bgc_diag_3d
      integer                                     :: m,diagidx2d,diagidx3d

!     1. Check hardcoded array size in ROMS matches MARBL expectations
!     ----------------------------------------------------------------------
      
!     Initialise "check" values to 0 and add 1 for each variable MARBL needs
      diag_cnt_sf_2d = 0 ! number of 2d surface flux diagnostics
      diag_cnt_sf_3d = 0 ! number of 3d surface flux diagnostics
      diag_cnt_it_2d = 0 ! number of 2d interior tendency diagnostics
      diag_cnt_it_3d = 0 ! number of 3d interior tendency diagnostics
            diag_cnt_tot = 0 ! total number of diagnostics
                diagidx2d = 1 ! index for incrementing 2d diagnostics
                diagidx3d = 1 ! index for incrementing 3d diagnostics

                
!     Go through both MARBL diagnostics types and check if fields are
!     2d or 3d one by one (calculate nr_marbl_ss_2d,nr_marbl_ss_3d), then
!     increment the size check variables accordingly:

!     Add surface flux diagnostics variables to the counts:      
      do m=1,size(MARBL_instance%surface_flux_diags%diags)
         diag_cnt_tot=diag_cnt_tot+1
         if (trim(MARBL_instance%surface_flux_diags%diags(m)
     &        %vertical_grid) .eq. "none") then ! 2D field
            diag_cnt_sf_2d=diag_cnt_sf_2d+1
         else
            diag_cnt_sf_3d=diag_cnt_sf_3d+1
         end if
      end do
      
!     Add interior tendency diagnostics variables to the counts:      
      do m=1,size(MARBL_instance%interior_tendency_diags%diags)
         diag_cnt_tot=diag_cnt_tot+1
         if (trim(MARBL_instance%interior_tendency_diags%diags(m)
     &        %vertical_grid) .eq. "none") then ! 2D field
            diag_cnt_it_2d=diag_cnt_it_2d+1
         else
            diag_cnt_it_3d=diag_cnt_it_3d+1
         end if
      end do

      diag_cnt_2d_tot = diag_cnt_sf_2d + diag_cnt_it_2d
      diag_cnt_3d_tot = diag_cnt_sf_3d + diag_cnt_it_3d
      
!     Throw error/advise if counts don't match                  
      if ( size(vname_bgc_diag_2d(1,:)) .ne.
     &     (diag_cnt_2d_tot) ) then
         
         write(*,'(7x,A,I4,A,I4,A,I4,A)')
     &        'ERROR: Allocated no. of 2D BGC diagnostics: '
     &        ,size(vname_bgc_diag_2d(1,:))
     &        ,' does not match no. expected by MARBL: '
     &        ,(diag_cnt_sf_2d+diag_cnt_it_2d)
     &        ,'. Set nr_bec2_diag_2d = '
     &        ,(diag_cnt_sf_2d+diag_cnt_it_2d)
     &        ,' in bgc_ecosys_vars.F and recompile'
         error stop 'ERROR reported from MARBL driver'
         
      elseif ( size(vname_bgc_diag_3d(1,:)) .ne.
     &        (diag_cnt_3d_tot) ) then
         
         write(*,'(7x,A,I4,A,I4,A,I4,A)')
     &        'ERROR: Allocated no. of 3D BGC diagnostics: '
     &        ,size(vname_bgc_diag_3d(1,:))
     &        ,' does not match no. expected by MARBL: '
     &        ,(diag_cnt_sf_3d+diag_cnt_it_3d)
     &        ,'. Set nr_bec2_diag_3d = '
     &        ,(diag_cnt_sf_3d+diag_cnt_it_3d)
     &        ,' in bgc_ecosys_vars.F and recompile'
         error stop 'ERROR reported from MARBL driver'
      end if

!     2. Fill BGC diagnostics metadata arrays in ROMS from MARBL instance
!     ----------------------------------------------------------------------

!     MARBL surface flux diagnostics:
      do m=1,size(MARBL_instance%surface_flux_diags%diags)
         ! Check if this diagnostic is 2D (true) or 3D (false/else)
         if (trim(MARBL_instance%surface_flux_diags%diags(m)%vertical_grid) .eq. "none") then ! 2D field
            vname_bgc_diag_2d(1,diagidx2d)=
     &           trim(MARBL_instance%surface_flux_diags%diags(m)%short_name)
            vname_bgc_diag_2d(2,diagidx2d)=
     &           trim(MARBL_instance%surface_flux_diags%diags(m)%long_name)
            vname_bgc_diag_2d(3,diagidx2d)=
     &           trim(MARBL_instance%surface_flux_diags%diags(m)%units)
            vname_bgc_diag_2d(4,diagidx2d)='  '
            wrt_bgc_diag_2d(diagidx2d)=.false.

            diagidx2d=diagidx2d+1
         else                   ! 3D field
            vname_bgc_diag_3d(1,diagidx3d)=
     &           trim(MARBL_instance%surface_flux_diags%diags(m)%short_name)
            vname_bgc_diag_3d(2,diagidx3d)=
     &           trim(MARBL_instance%surface_flux_diags%diags(m)%long_name)
            vname_bgc_diag_3d(3,diagidx3d)=
     &           trim(MARBL_instance%surface_flux_diags%diags(m)%units)
            vname_bgc_diag_3d(4,diagidx3d)='  '
            wrt_bgc_diag_3d(diagidx3d)=.false.
!
            diagidx3d=diagidx3d+1
         end if
      end do

!     MARBL interior tendency diagnostics
      do m=1,size(MARBL_instance%interior_tendency_diags%diags)
         if (trim(MARBL_instance%interior_tendency_diags%diags(m)%vertical_grid) .eq. "none") then ! 2D field
            vname_bgc_diag_2d(1,diagidx2d)=
     &           trim(MARBL_instance%interior_tendency_diags%diags(m)%short_name)
            vname_bgc_diag_2d(2,diagidx2d)=
     &           trim(MARBL_instance%interior_tendency_diags%diags(m)%long_name)
            vname_bgc_diag_2d(3,diagidx2d)=
     &           trim(MARBL_instance%interior_tendency_diags%diags(m)%units)
               vname_bgc_diag_2d(4,diagidx2d)='  '
               wrt_bgc_diag_2d(diagidx2d)=.false.

            diagidx2d=diagidx2d+1
         else                   ! 3D field
            vname_bgc_diag_3d(1,diagidx3d)=
     &           trim(MARBL_instance%interior_tendency_diags%diags(m)%short_name)
            vname_bgc_diag_3d(2,diagidx3d)=
     &           trim(MARBL_instance%interior_tendency_diags%diags(m)%long_name)
            vname_bgc_diag_3d(3,diagidx3d)=
     &            trim(MARBL_instance%interior_tendency_diags%diags(m)%units)
            vname_bgc_diag_3d(4,diagidx3d)='  '
            wrt_bgc_diag_3d(diagidx3d)=.false.

            diagidx3d=diagidx3d+1
         end if
      end do

!     Read user file to determine which diagnostics to write out
      call marbldrv_read_diagnostic_output_list(
     &     vname_bgc_diag_2d,vname_bgc_diag_3d
     &     )                    ! 
      
      if (mynode==printnode) then
         print *, 'number of 2D marbl diags', diag_cnt_2d_tot
         print *, 'number of 3D marbl diags', diag_cnt_3d_tot
      end if
      
      end subroutine marbldrv_configure_diagnostics


      subroutine marbldrv_read_diagnostic_output_list(
     &     vname_bgc_diag_2d,vname_bgc_diag_3d
     &     )                    ! 
!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_read_diagnostic_output_list
!     DESCRIPTION: Read "marbl_diagnostic_output_list" file and set 
!                  wrt_bec2_diag logical array entries accordingly
!
!     METHOD:
!     1. Read the "marbl_diagnostic_output_list" file line-by-line
!     2. For each line not beginning with "!", compare that line with
!         every entry in "vname_bec2_diag" (diagnostic shortname list)
!     3. If there is a match, set "wrt_bec2_diag" to true for that diag
!     4. If there are no matches, default to writing out  all diagnostics
!
!     NOTES:
!     - The possible options are listed in STDOUT or can be found with
!     $MARBL_ROOT/tests/regression_tests/requested_diagnostics/requested_diagnostics.py
!-----------------------------------------------------------------------

      character*72, dimension(:,:), intent(inout) :: vname_bgc_diag_2d
      character*72, dimension(:,:), intent(inout) :: vname_bgc_diag_3d

      character(len=256)                         :: file_line
      integer :: marbl_diag_list, read_status,open_status,idx,ierr
      logical :: recognized_diag, no_diags_requested

      marbl_diag_list = 11
      no_diags_requested = .true. ! If we find any requested tracers, flip to false

!     1. Read in MARBL diagnostics list file "marbl_diagnostic_output_list"
!     ------------------------------------------------------------------
      
      if (mynode==0) then
         open(unit=marbl_diag_list, file='marbl_diagnostic_output_list',
     &        iostat=open_status)
         
         if (open_status .ne. 0) then ! can't find file, skip and output all diags
            write(*, '(A, I0, 3A)') "IO WARNING ", open_status,
     &           "could not open diagnostics list file : ",
     &           "marbl_diagnostic_output_list",
     &           ". All MARBL diagnostics will be written to output"
         end if
      end if
#ifdef MPI
      call MPI_Bcast(open_status,1,MPI_INTEGER,0,ocean_grid_comm,ierr)
#endif
      if (open_status==0) then            
         file_line = ''
         do                     !line-by-line read
            if (mynode==0) then
               read(marbl_diag_list,"(A)",iostat=read_status),
     &              file_line
            end if
#ifdef MPI
            call MPI_Bcast(read_status,1,MPI_INTEGER,0,ocean_grid_comm,ierr)
#endif             
            if (read_status .ne. 0) then
               exit             !loop
            else
#ifdef MPI
               call MPI_Bcast(file_line,256,MPI_CHARACTER,0,ocean_grid_comm,ierr)
#endif                
!     Check if the line is a comment
               file_line=trim(adjustl(file_line))
               if (file_line(1:1) == '!') cycle ! Skip this iteration if it's a comment
!     Set the "wrt_t" variable
               recognized_diag=.false.
               do idx=1,size(vname_bgc_diag_2d(1,:))
                  if (trim(vname_bgc_diag_2d(1,idx)) == trim(file_line))
     &                 then
                     wrt_bgc_diag_2d(idx) =.true.
                     if (mynode==printnode) then
                        write(*,'(2A)'),
     &                       'REQUESTED MARBL DIAGNOSTIC OUTPUT: ' //
     &                       trim(vname_bgc_diag_2d(1,idx))
                     end if 
                     recognized_diag=.true.
                     no_diags_requested=.false.
                     cycle
                  end if
               end do           ! idx
            
               do idx=1,size(vname_bgc_diag_3d(1,:))
                  if (trim(vname_bgc_diag_3d(1,idx)) == trim(file_line))
     &                 then
                     wrt_bgc_diag_3d(idx) =.true.
                     if (mynode==printnode) then
                        write(*,'(2A)'),
     &                       'REQUESTED MARBL DIAGNOSTIC OUTPUT: ' //
     &                       trim(vname_bgc_diag_3d(1,idx))
                     end if 
                     recognized_diag=.true.
                     no_diags_requested=.false.
                     cycle
                  end if
               end do           ! idx

               if ( (.not. recognized_diag) .and. (mynode==printnode))
     &              then
               
                  write(*,'(A)'), 'REQUESTED DIAGNOSTIC ' //
     &                 trim(file_line) //
     &                 ' NOT RECOGNIZED. See ' //
     &                 ' $MARBL_ROOT/tests/regression_tests/' // 
     &                 'requested_diagnostics/requested_diagnostics.py ' //
     &                 'for a list of available diagnostics,' //
     &                 ' or run without the file ' //
     &                 'marbl_diagnostic_output_list to write all diags'
               end if
            end if              ! read_status=/=0
            
         end do                 !/namelist line-by-line read

         if (no_diags_requested) then ! Write all of them to output
               
            wrt_bgc_diag_2d(:)=.true.
            wrt_bgc_diag_3d(:)=.true.               
               
            if (mynode==0) then
               write(*,'(A)'), 'NO REQUESTED DIAGNOSTICS FOUND IN ' //
     &              ' marbl_diagnostic_output_list. Defaulting to ' //
     &              ' writing all BGC diagnostics to output file.'
            end if
         end if
         
         if (is_iostat_end(read_status)) then
            if (mynode==printnode) then
               write(*, '(3A)') "Successfully read ",
     &              "marbl_diagnostic_output_list",
     &              " to end of file"
            end if
         else
            if (mynode==printnode) then
               write(*, '(A, I0, 2A)') "IO ERROR ", read_status,
     &              " : ERROR reading MARBL diagnostic list file  ",
     &              "marbl_diagnostic_output_list"
            end if
            error stop
         end if
         
         close(marbl_diag_list)
         if (mynode==printnode) then
            write(*,'(2A)'), 'we successfully read and closed ',
     &           'marbl_diagnostic_output_list'
            
         end if
      endif         

      end subroutine marbldrv_read_diagnostic_output_list
      
#endif /* MARBL_DIAGS */
      
      subroutine marbldrv_column_physics(istr,iend,jstr,jend,tracer_array)
      
!-----------------------------------------------------------------------
!     SUBROUTINE: marbldrv_column_physics
      
!     DESCRIPTION:
!     Calculate surface fluxes, interior tracer tendencies,
!     saved state variables, and diagnostics using MARBL,
!     then update values in ROMS 
!      
!     METHOD:
!     1. Populate MARBL instance surface flux forcing values using ROMS
!        forcing variables
!     2. Populate MARBL surface flux saved state values from array in
!        ROMS (either from previous time step or read in restart file)
!     3. Compute surface fluxes using MARBL
!     4a. Update saved state array in ROMS using newly computed values
!     4b. Update diagnostics array in ROMS using newly computed values
!     5. Populate MARBL instance interior tendency forcing values
!        using ROMS ocean interior values
!     6. Populate MARBL interior tendency saved state values from
!        array in ROMS (either from previous time step or restart file)
!     7. Compute interior tendencies using MARBL
!     8a. Apply newly computed increments to ROMS tracer array
!     8b. Update saved state array in ROMS using newly computed values
!     8c. Update diagnostics array in ROMS using newly computed values      
!                 
!     INPUT/OUTPUT:
!     - istr,jstr,iend,jend: i,j indices on this node over which to loop
!     - tracer_array: the ROMS tracer array (t) to be updated by MARBL
!      
!     NOTES:
!     - All 3D variables are flipped when passing between MARBL and ROMS
!       as MARBL's indexing is 1 at surface and ROMS' is 1 at bottom.
!     - Ice fraction is everywhere set to 0 as no ice in ROMS currently
!     - Currently unable to use `press(:,:)` so atmpress is set to 1
!     - Currently unable to use `p(:,:,:)` so interior pressure is set
!       based on depth
!     - "PAR subcols" is unused as the dimension of this MARBL
!       variable is "ice categories" but we have no ice categories
!     - Bottom fluxes (and river inputs) are expected by MARBL but
!       currently are not set
!     - The interior tendency forcing variable "fesedflux" is unused.
!       It comes from the use of a 1/60* bathymetry dataset in MOM,
!       as  a sort of subgrid-scale way of putting iron where the bed
!       would be in reality even where the bed is unresolved
!       by suspending it at the equivalent point in the water column
!     - o2_scalef_ind and remin_scalef_ind are similarly unused.
!       They were implemented in MARBL to deal with OMZs in CESM2.

!     
!-----------------------------------------------------------------------
     
      integer, intent(in)                           :: istr,jstr,iend,jend

      real, dimension(GLOBAL_2D_ARRAY,nz,3,nt) ,
     &                            intent(inout)     :: tracer_array      
#ifdef MARBL_DIAGS
!      real, dimension(GLOBAL_2D_ARRAY,diag_cnt_2d_tot),
!     &                         intent(inout)        :: diag_array_2d
!      real, dimension(GLOBAL_2D_ARRAY,nz,diag_cnt_3d_tot),
!     &                         intent(inout)        :: diag_array_3d
      
      integer                                       :: diagidx2d,diagidx3d ! indices for stepping through diagnostics arrays
#endif
      integer                                       :: i,j,k,m
      logical                                       :: print_values = .false. ! Print example values to stdout during loop

      
!     1. Populate MARBL instance surface flux forcing values:
!     ----------------------------------------------------------------------
      do j=jstr,jend
         do i=istr,iend

            if (rmask(i,j)==0.) cycle
            
! sea surface salinity
            if (sss_ind > 0) then
               marbl_instance%surface_flux_forcings(sss_ind)%field_0d(1)
     &              =tracer_array(i,j,nz,nnew,isalt) ! (psu)
               if (print_values) call sf_test_print(sss_ind,i,j)
            end if
! sea surface temperature
            if (sst_ind > 0) then
               marbl_instance%surface_flux_forcings(sst_ind)%field_0d(1)
     &              =tracer_array(i,j,nz,nnew,itemp) ! (degC)
               if (print_values) call sf_test_print(sst_ind,i,j)
            end if
! ice fraction (unused)
            if (ifrac_ind > 0) then
               marbl_instance%surface_flux_forcings(ifrac_ind)%field_0d(1)=0
!     &              ifrac(i,j) ! (dimensionless)
               if (print_values) call sf_test_print(ifrac_ind,i,j)
            end if
! squared 10m windspeed            
            if (u10_sqr_ind > 0) then
               marbl_instance%surface_flux_forcings(u10_sqr_ind)%field_0d(1)=
     &               (uwnd(i,j)**2) + (vwnd(i,j)**2) ! (m2/s2)
               if (print_values) call sf_test_print(u10_sqr_ind,i,j)
            end if
! atmospheric pressure (constant)            
            if (atmpress_ind > 0) then
            marbl_instance%surface_flux_forcings(atmpress_ind)%field_0d(1)=1.
!     &              press(i,j) !(atm)
            if (print_values) call sf_test_print(atmpress_ind,i,j)            
         end if
! pco2         
            if (xco2_ind > 0) then
               marbl_instance%surface_flux_forcings(xco2_ind)%field_0d(1)=
     &              pco2air(i,j) !(ppm)
               if (print_values) call sf_test_print(xco2_ind,i,j)
            end if
! pco2 (alternative co2 forcing)
            if (xco2_alt_ind > 0) then
               marbl_instance%surface_flux_forcings(xco2_alt_ind)%field_0d(1)=
     &              pco2air_alt(i,j) !(ppm)
               if (print_values) call sf_test_print(xco2_alt_ind,i,j)
            end if
! dust
            if (dust_dep_ind > 0) then
               marbl_instance%surface_flux_forcings(dust_dep_ind)%field_0d(1)
     &              = dust(i,j) !(kg/m2/s)
               if (print_values) call sf_test_print(dust_dep_ind,i,j)
            end if
! iron deposition            
            if (fe_dep_ind > 0) then
               marbl_instance%surface_flux_forcings(fe_dep_ind)%field_0d(1)=
     &              iron(i,j)*0.01 ! nmol/cm2/s (ROMS) -> mmol/m2/s (MARBL MKS)
               if (print_values) call sf_test_print(fe_dep_ind,i,j)
            end if
! nox
#ifdef NOX_FORCING
            if (nox_flux_ind > 0) then
               marbl_instance%surface_flux_forcings(nox_flux_ind)%field_0d(1)=               
     &              nox(i,j)*71394.200220751 ! kg(N)/m2/s (ROMS) -> mmol/m2/s (MARBL MKS)
               if (print_values) call sf_test_print(nox_flux_ind,i,j)
            end if
#endif            
! nhy
#ifdef NHY_FORCING 
            if (nhy_flux_ind > 0) then
               marbl_instance%surface_flux_forcings(nhy_flux_ind)%field_0d(1)=
     &              nhy(i,j)*71394.200220751 ! kg(N)/m2/s (ROMS) -> mmol/m2/s (MARBL MKS)
               if (print_values) call sf_test_print(nhy_flux_ind,i,j)
            end if
#endif
! Surface tracer values            
            do m=1,nt_marbl
               marbl_instance%tracers_at_surface(1,m)=
     &              tracer_array(i,j,nz,nnew, m+(NT-nt_marbl))
               
!     Print surface values at desired location if requested
               if ( (i==printi) .and. (j==printj)
     &              .and. (mynode==printnode) .and. (print_values) ) then
                  print *, 'We just set surface ',
     &                 (trim(MARBL_instance%tracer_metadata(m)%short_name)),
     &                 ' as ',marbl_instance%tracers_at_surface(1,m),
     &                 (trim(MARBL_instance%tracer_metadata(m)%units))
               end if
!
            end do

!     2. Populate MARBL surface flux saved state values:
!     ----------------------------------------------------------------------
            
            do m=1,MARBL_instance%surface_flux_saved_state%saved_state_cnt
               MARBL_instance%surface_flux_saved_state%state(m)%field_2d(1)
     &              =marbl_saved_state_2d(i,j,m)               
               
!     Print surface flux saved state values at a chosen location if desired:
               if ( (i==printi) .and. (j==printj)
     &              .and. (mynode==printnode) .and. (print_values) ) then
                  print *, 'We just set surface ',
     &                 trim(MARBL_instance%surface_flux_saved_state%
     &                 state(m)%short_name),
     &                 ' as ',MARBL_instance%surface_flux_saved_state%
     &                 state(m)%field_2d(1),
     &                 trim(MARBL_instance%surface_flux_saved_state%
     &                 state(m)%units)
               end if
               
            end do
               
!     3. Compute surface fluxes using MARBL:
!     ----------------------------------------------------------------------
            
            call marbl_instance%surface_flux_compute()

            if (marbl_instance%StatusLog%labort_marbl) then
               print *, 'abort during surface flux compute'
               call marbl_instance%StatusLog%log_error_trace(
     &              "marbl_instance%surface_flux_compute()",
     &              "MARBL_tracers_column_physics")
               call print_marbl_log(marbl_instance%StatusLog,i,j)
            end if

            call marbl_instance%StatusLog%erase()

!     4a. Update saved state array in ROMS using newly computed values:
!     ----------------------------------------------------------------------            
            
            do m=1,MARBL_instance%surface_flux_saved_state%saved_state_cnt
               marbl_saved_state_2d(i,j,m)=
     &              MARBL_instance%surface_flux_saved_state%state(m)%field_2d(1)
            end do
            
!     4b. Update diagnostics array in ROMS using newly computed values:
!     ----------------------------------------------------------------------            
#ifdef MARBL_DIAGS
            diagidx2d=1 ! Initialise array indices to increment
            diagidx3d=1
            do m=1,size(MARBL_instance%surface_flux_diags%diags)
               if (trim(MARBL_instance%surface_flux_diags%diags(m
     &              )%vertical_grid) .eq. "none") then ! 2D field

                  if (wrt_bgc_diag_2d(diagidx2d)) then
                     diag_array_2d(i,j,idx_bgc_diag_2d(diagidx2d))=
     &                    real(MARBL_instance%surface_flux_diags%diags(m)%field_2d(1))
                  end if                     
                  diagidx2d=diagidx2d+1

               else             ! 3D field
                  if (wrt_bgc_diag_3d(diagidx3d)) then
                     diag_array_3d(i,j,:,idx_bgc_diag_2d(diagidx3d))=
     &                    real(MARBL_instance%surface_flux_diags%diags(m)%field_3d(1,:))
                  end if                     
                  diagidx3d=diagidx3d+1

               end if
            end do
#endif       
           
         end do   ! j=jstr,jend
      end do      ! i=istr,iend
      
!     5. Populate MARBL instance interior tendency forcing values:      
!     ----------------------------------------------------------------------

      do j=jstr,jend
         do i=istr,iend
            if (rmask(i,j)==0.) cycle
            
! First, update domain in MARBL instance to local geometry:
            
            MARBL_instance%domain%zw(:)      = -z_w(i,j,nz-1:0:-1) ! bottom interface depth
            MARBL_instance%domain%zt(:)      = -z_r(i,j,nz  :1:-1) ! centre depth
            MARBL_instance%domain%delta_z(:) = Hz(  i,j,nz  :1:-1) ! thickness
            MARBL_instance%domain%kmt        = nz                  ! number of active levels
            
! dust flux
            if (dustflux_ind > 0) then
               MARBL_instance%interior_tendency_forcings(dustflux_ind)%field_0d(1)=
     &              dust(i,j)
            end if
! PAR subcolumns (unused as the dimension of this field is the number of ice categories and we have none)
            if (PAR_col_frac_ind > 0) then
               continue
               !MARBL_instance%interior_tendency_forcings(PAR_col_frac_ind)%field_1d(1,:)=
            end if
! surface shortwave
            if (surf_shortwave_ind > 0) then               
               MARBL_instance%interior_tendency_forcings(surf_shortwave_ind)%field_1d(1,1)=
     &              srflx(i,j)
            end if
! potential temperature
            if (potemp_ind > 0) then
               MARBL_instance%interior_tendency_forcings(potemp_ind)%field_1d(1,:)
     &              =tracer_array(i,j,nz:1:-1,nnew,itemp)
               if (print_values) call if_test_print(potemp_ind,i,j)
            end if
! salinity
            if (salinity_ind > 0) then
               MARBL_instance%interior_tendency_forcings(salinity_ind)%field_1d(1,:)
     &              =tracer_array(i,j,nz:1:-1,nnew,isalt)
               if (print_values) call if_test_print(salinity_ind,i,j)               
            end if
! pressure (currently calculating from depth)
            if (pressure_ind > 0) then
               MARBL_instance%interior_tendency_forcings(pressure_ind)%field_1d(1,:)=
     &              -z_r(i,j,nz:1:-1)*0.1
!     &              =p(i,j,nz:1:-1)
               if (print_values) call if_test_print(pressure_ind,i,j)
            end if

! iron flux
            if (fesedflux_ind > 0) then
               MARBL_instance%interior_tendency_forcings(fesedflux_ind)%field_1d(1,:)=0.
               if (print_values) call if_test_print(fesedflux_ind,i,j)
            end if


! column tracers
            do m=1,nt_marbl
               marbl_instance%tracers(m,:)=                  
     &              tracer_array(i,j,nz:1:-1,nnew, m+(NT-nt_marbl) )
               
!     Print bottom value at a desired location if requested:
               if ( (i==printi) .and. (j==printj)
     &              .and. (mynode==printnode) .and. (print_values) ) then
                  print *, 'We just set bottom ',
     &                 (trim(MARBL_instance%tracer_metadata(m)%short_name)),
     &                 ' as ',marbl_instance%tracers(m,nz),
     &                 (trim(MARBL_instance%tracer_metadata(m)%units))                  
               end if

            end do

!     6. Populate MARBL interior tendency saved state values:
!     ----------------------------------------------------------------------
            
            
            do m=1,MARBL_instance%interior_tendency_saved_state%saved_state_cnt
               MARBL_instance%interior_tendency_saved_state%state(m)%field_3d(:,1)
     &              =marbl_saved_state_3d(i,j,nz:1:-1,m)
               
!     Print interior tendency saved state bottom values at a chosen location if desired:
               if ( (i==printi) .and. (j==printj)
     &              .and. (mynode==printnode) .and. (print_values) ) then
                  print *, 'We just set bottom ',
     &                 trim(MARBL_instance%
     &                 interior_tendency_saved_state%state(m)%long_name)
     &                 ,' as ',marbl_instance%
     &                 interior_tendency_saved_state%state(m)%field_3d(nz,1)
     &                 ,trim(MARBL_instance%interior_tendency_saved_state%state(m)%units) 
               end if
               
            end do    

!     Populate bottom fluxes (0 for now)
            marbl_instance%bot_flux_to_tend(:)=0.
            marbl_instance%bot_flux_to_tend(nz)=
     &           1./MARBL_instance%domain%delta_z(nz) 
            

!     7. Compute interior tendencies using MARBL:
!     ----------------------------------------------------------------------

            call marbl_instance%interior_tendency_compute()

            if (marbl_instance%StatusLog%labort_marbl) then
               print *, 'abort during interior tendency compute'
               call marbl_instance%StatusLog%log_error_trace(
     &             "marbl_instance%interior_tendency_compute()",
     &              "MARBL_tracers_column_physics")
               call print_marbl_log(marbl_instance%StatusLog,i,j)
               !error stop 'ERROR reported from MARBL library'                           
            end if
!            call print_marbl_log(marbl_instance%StatusLog,i,j)
            call marbl_instance%StatusLog%erase()
            
           
!     8a. Apply calculated increments to ROMS tracer array:
!     ----------------------------------------------------------------------

            do m=1,nt_marbl
               tracer_array(i,j,:,nnew, m+(NT-nt_marbl) ) =
     &              tracer_array(i,j,:,nnew, m+(NT-nt_marbl) ) +
     &              marbl_instance%interior_tendencies(m,nz:1:-1)*dt
            end do
            
!     8b. Update saved state array in ROMS using newly computed values:
!     ----------------------------------------------------------------------
            
            do m=1,MARBL_instance%interior_tendency_saved_state%saved_state_cnt
               marbl_saved_state_3d(i,j,:,m)=
     &             MARBL_instance%interior_tendency_saved_state%state(m)%field_3d(nz:1:-1,1)
            end do    
            
#ifdef MARBL_DIAGS

!     8c. Update diagnostics array in ROMS using newly computed values:
!     ----------------------------------------------------------------------
            
            diagidx2d=diag_cnt_sf_2d+1 ! Initialise array indices to increment
            diagidx3d=diag_cnt_sf_3d+1
            
            do m=1,size(MARBL_instance%interior_tendency_diags%diags)
               if (trim(MARBL_instance%interior_tendency_diags%diags(m
     &              )%vertical_grid) .eq. "none") then ! 2D field
                  if (wrt_bgc_diag_2d(diagidx2d)) then                  
                     diag_array_2d(i,j,idx_bgc_diag_2d(diagidx2d))=
     &                    real(MARBL_instance%interior_tendency_diags%diags(m)%field_2d(1))
                  end if                     
                  diagidx2d=diagidx2d+1
               else             ! 3D field
                  if (wrt_bgc_diag_3d(diagidx3d)) then                  
                     diag_array_3d(i,j,:,idx_bgc_diag_3d(diagidx3d))=
     &                    real(MARBL_instance%interior_tendency_diags%diags(m)%field_3d(:,1))
                  end if                      
                  diagidx3d=diagidx3d+1

               end if
            end do
#endif       
            
         end do                 !j=jstr,jend
      end do                    ! i=istr,iend


      
      
      end subroutine marbldrv_column_physics


      subroutine print_marbl_log(log_to_print,i,j)

!-----------------------------------------------------------------------
!     SUBROUTINE: print_marbl_log
!     DESCRIPTION:
!     Get MARBL status log info from MARBL instance and print to STDOUT.
!
!     METHOD:
!     Assemble the message from the i,j location (if known) the current
!     CPU (if not the main CPU), and the status log message from MARBL.
!     Print the message and abort the program if MARBL requires it
!      
!     INPUTS/OUTPUTS:
!     - log_to_print: the MARBL instance's StatusLog
!     - i,j (optional): the location of the message, if requested
!
!     NOTES:
!     - Adapted from the equivalent subroutine in MOM6
!     
!-----------------------------------------------------------------------
      
      use marbl_logging, only: marbl_status_log_entry_type, marbl_log_type
!     ----------------------------------------------------------------------
      integer              , optional, intent(in) :: i,j      
      class(marbl_log_type),           intent(in) :: log_to_print
      character(len=256)                          :: message_prefix,
     &                                               message_location,
     &                                               log_message
       
      type(marbl_status_log_entry_type), pointer  :: marbl_status_log_pointer
      
      marbl_status_log_pointer => log_to_print%FullLog

!     ----------------------------------------------------------------------

      
      do while (associated(marbl_status_log_pointer))
         write(message_prefix, "(A,I0,A)") "Task ",mynode,")"
      
! Check we are either on principal cpu OR all cpus can write:
         if ((.not. marbl_status_log_pointer%lonly_master_writes) .or.
     &        (mynode==printnode)) then
            
! Check if we know the location of this log message and write
            if (marbl_status_log_pointer%ElementInd .gt. 0) then
               if (present(i) .and. present(j)) then
                  write(message_location, "(A,I0,A,I0,A,I0)")
     &                 "Message from (i,j) (", i,",",j,") at level ",
     &                 marbl_status_log_pointer%ElementInd
               else
                  write(message_location, "(A)")
     &                 "grid cell responsible for message unknown"
               end if ! i,j present
               

! Also write the calling CPU unless on the principal CPU
               if (mynode==printnode) then
                  write(log_message,"(A)") trim(message_location)
               else
                  write(log_message,"(A,1X,A)")
     &                 trim(message_prefix), trim(message_location)
               end if           ! (mynode==printnode)
               print *, trim(log_message)
            end if              ! (ElementInd > 0)

! Now write the actual statuslog message from MARBL
            if (mynode==printnode) then
               write(log_message,"(A)")
     &              trim(marbl_status_log_pointer%LogMessage)
            else
               write(log_message,"(A,1X,A)")
     &              trim(message_prefix),
     &              trim(marbl_status_log_pointer%LogMessage)
            end if              ! (mynode==printnode)
            print *, trim(log_message)
         end if                 ! mynode 0 or not master task only
            
         marbl_status_log_pointer => marbl_status_log_pointer%next
      end do                    ! while(associated(marbl_status_log_pointer))

! If MARBL needs to abort, stop the program      
      if (log_to_print%labort_marbl) then
         error stop 'ERROR reported from MARBL library'            
      end if
         
         
      
      end subroutine print_marbl_log

      
      subroutine sf_test_print(ind,i,j)

!-----------------------------------------------------------------------
!     SUBROUTINE: sf_test_print
!     DESCRIPTION:
!     Print MARBL example surface flux forcing variable value to STDOUT
!
!     METHOD:
!     Check if we are at the example location (printi,printj) and on
!     the chosen CPU for printing (printnode), then print the value
!      
!     INPUTS/OUTPUTS:
!     - i,j : the current location, to be compared with printi,printj
!     - ind : the index of the surface flux forcing variable
      
!     NOTES:
!     
!-----------------------------------------------------------------------
      
      
      integer, intent(in) :: i,j,ind

      if (j==printj) then
         if (i==printi) then
            if (mynode==printnode) then
               print *, 'We just set ',
     &              (trim(MARBL_instance%surface_flux_forcings(ind)%metadata%varname)),
     &              ' as ',marbl_instance%surface_flux_forcings(ind)%field_0d(1),
     &              (trim(MARBL_instance%surface_flux_forcings(ind)%metadata%field_units))
            end if                     
         end if
      end if                              
      end subroutine sf_test_print

      subroutine if_test_print(ind,i,j)
!-----------------------------------------------------------------------
!     SUBROUTINE: if_test_print
!     DESCRIPTION:
!     Print MARBL example interior tendency forcing var value to STDOUT
!
!     METHOD:
!     Check if we are at the example location (printi,printj) and on
!     the chosen CPU for printing (printnode), then print bottom value
!      
!     INPUTS/OUTPUTS:
!     - i,j : the current location, to be compared with printi,printj
!     - ind : the index of the interior tendency forcing variable
      
!     NOTES:
!     
!-----------------------------------------------------------------------
      
      integer, intent(in) :: i,j,ind

      if (j==printj) then
         if (i==printi) then
            if (mynode==printnode) then
               print *, 'We just set bottom ',
     &              (trim(MARBL_instance%interior_tendency_forcings(ind)%metadata%varname)),
     &              ' as ',marbl_instance%interior_tendency_forcings(ind)%field_1d(1,nz)
            end if                     
         end if
      end if                              
      end subroutine if_test_print
#endif      
      end module marbl_driver

