      module pipe_frc

      ! Pipes forcing module
      ! --------------------

      ! for pipes:  Qbar(nx,y), the total pipe flux for that grid point ! (m3/s) function of time
      !             Qshape(npipes), the vertical shape of the pipe distribution ! fractions sum(Qshape) = 1

      ! Initial coding by Jeroen Molemaker & Devin Dollery (2020 Nov)

#include "cppdefs.opt"
! we have no PIPE_SOURCE flag here anymore, since we are moving away from cpp flags.
! Variables still need to be visible within code even when pipes are not used.


      use dimensions                                       ! has dimensions, a list of forcing files, and rmask
      use roms_read_write
      use nc_read_write
      use scalars

      implicit none

      private

#include "pipe_frc.opt"

      ! Variables used in the evolution equations
      integer,public  :: pidx ! Pipe index for looping through pipes
      integer,public,allocatable,dimension(:,:) :: pipe_idx      ! pipe indices at grid points
      real   ,public,allocatable,dimension(:,:) :: pipe_fraction ! pipe fractional flux at grid points
      real   ,public,allocatable,dimension(:,:) :: pipe_flx      ! pipe flux
      real   ,public, dimension(npip,N)         :: pipe_prf      ! Pipe vertical profile
      real,public,allocatable,dimension(:,:)    :: pipe_trc      ! Pipe tracer conc.
      real   ,public, dimension(npip)           :: pipe_vol      ! Pipe volume

     ! Netcdf names
      character(len=9)  :: pipe_flx_name = 'pipe_flux'  !! stored in the grid file
      character(len=11) :: pipe_vol_name = 'pipe_volume'!! stored in a forcing file
      character(len=11) :: pipe_trc_name = 'pipe_tracer'!! stored in a forcing file
      character(len=9)  :: pipe_tim_name = 'pipe_time'  !! stored in a forcing file
      character(len=5)  ::npipe_dim_name = 'npipe'      !! dimension name for number of pipes in file
      character(len=8)  :: ntrc_dim_name = 'ntracers'   !! dimension name for number of tracers in file


      logical :: init_pipe_done = .false.

      public set_pipe_frc
      public init_arrays_pipes

      contains

! ----------------------------------------------------------------------
      subroutine set_pipe_frc  ![
      ! set pipe forces (realistic and analytical)
      ! - read and interpolation all pipe forcing.
      ! - all pipe variables need time interpolation only
      !   here so can use same generic routine.
      ! - input data in days!
      implicit none

      ! local
      integer :: i,j

      if (.not. init_pipe_done) call init_pipe_frc

      if (p_analytical) then                                ! Set pipe flux volumes and tracer data

        call set_ana_pipe_frc

      else

        call set_frc_data(nc_pvol,pipe_vol)                ! set pipe volume for all pipes at current time

        call set_frc_data(nc_ptrc,var2d=pipe_trc)          ! set pipe tracers conc. for all pipes at current time

        call set_pipe_vert_prf                             ! set pipe vertical profiles

        do j = 1,ny                                        ! calculate fluxes with time interpolated volumes and conc.:
          do i = 1,nx
            if (pipe_idx(i,j) > 0.) then
              pipe_flx(i,j)=pipe_fraction(i,j)*pipe_vol( pipe_idx(i,j) )
            endif
          enddo
        enddo

      endif

      end subroutine set_pipe_frc !]
! ----------------------------------------------------------------------
      subroutine init_pipe_frc  ![

      ! Initialize pipe forcing:
      ! Read in a grid file with locations of pipes and flux contribution per cell.
      ! Done only once as pipe position does not change.
      ! Realistic case - stored as one value in NetCDF file where
      ! pipe grid point value = pidx + pipe_fraction

      use netcdf
      implicit none

      ! local
      integer :: ierr,ncid,v_id,i,j


      if (p_analytical) then

        ! pipe_flx is defined in ana_pipe_frc.h

      else

        ! Read 'pipe_flux' from grid file (Pipe_idx & pipe_fraction in one value)
        ierr=nf90_open(grdname, nf90_nowrite, ncid)
        ! Temporarily store as pipe_fraction to avoid extra array,
        ! but value still pipe_idx + pipe_fraction
        call ncread(ncid,pipe_flx_name, pipe_fraction(i0:i1,j0:j1))
        if(ierr/=0)
     &      call handle_ierr(ierr,'init_pipe_frc:: pipe not in file!')

        ! Separate pipe_idx & pipe_fraction:
        do j = 1,ny
          do i = 1,nx
            ! read in value = pidx + pipe_fraction = pipe_idx(i,j) + pipe_fraction(i,j)
            pipe_idx(i,j)=floor(pipe_fraction(i,j)-1e-5)
            pipe_fraction(i,j)=pipe_fraction(i,j)-pipe_idx(i,j)
          enddo
        enddo

      endif

      init_pipe_done = .true.
      if(mynode==0) write(*,'(/7x,A/)') 'pipe_frc: init pipe locations'

      end subroutine init_pipe_frc  !]
! ----------------------------------------------------------------------
      subroutine init_arrays_pipes  ![
      implicit none

      character(len=30) :: string
      allocate( pipe_idx(GLOBAL_2D_ARRAY) );      pipe_idx=0.
      allocate( pipe_fraction(GLOBAL_2D_ARRAY) ); pipe_fraction=0.
      allocate( pipe_flx(GLOBAL_2D_ARRAY) );      pipe_flx=0.
      allocate( pipe_trc(npip,nt) ) ;             pipe_trc=0.

      if (.not. p_analytical) then
        allocate(nc_pvol%vdata(npip,1 ,2))
        allocate(nc_ptrc%vdata(npip,nt,2))
      endif

      ! Print user options (pipe_frc.opt) to netcdf attributes
      pipe_frc_opt = ''
      write(string, '(A,I3)') '# of pipes: ', npip
      call store_string_att(pipe_frc_opt, string)
      if (p_analytical) then
        call store_string_att(pipe_frc_opt, ', analytical')
      else
        call store_string_att(pipe_frc_opt, ', real')
      endif
      pipe_frc_opt = trim(adjustl(pipe_frc_opt))

      end subroutine init_arrays_pipes  !]
! ----------------------------------------------------------------------
      subroutine set_pipe_vert_prf  ![
      ! set vertical discharge profile of pipes

      ! This is a time-independent equation for now, but
      ! it is a placeholder for a more sophisticated time-evolving
      ! profile in future as per the requirement.

      implicit none

      ! local
      integer :: i,j,ipipe

      ! Loop through all pipes and set the same profile:
      do ipipe=1,npip
        pipe_prf(ipipe,:)=0     ! Set all values to zero
        pipe_prf(ipipe,1)= 0.5  ! Dispersion profile bottom cell
        pipe_prf(ipipe,2)= 0.5  ! Dispersion profile 2nd from bottom cell
      enddo

      end subroutine set_pipe_vert_prf  !]
! ----------------------------------------------------------------------
      subroutine set_ana_pipe_frc  ![
      ! Analytical pipe forcing volume and tracer data
      ! Put here to avoid circular reference if in analytical.F

#include "ana_pipe_frc.h"

      end subroutine set_ana_pipe_frc  !]
! ----------------------------------------------------------------------

      end module pipe_frc
