      module pipe_frc

      ! Pipes forcing module
      ! --------------------

      ! for pipes:  Qbar(nx,y), the total pipe flux for that grid point ! (m3/s) function of time
      !             Qshape(npipes), the vertical shape of the pipe distribution ! fractions sum(Qshape) = 1

      ! Initial coding by Jeroen Molemaker & Devin Dollery (2020 Nov)

#include "cppdefs.opt"
#if defined PIPE_SOURCE
  
      use dimensions !! has dimensions, a list of forcing files, and rmask
      use read_write
      use scalars

      implicit none

      private

#include "pipe_frc.opt"

      ! Variables used in the evolution equations
      integer,public  :: pidx ! Pipe index for looping through pipes
      integer,public, dimension(GLOBAL_2D_ARRAY) :: pipe_idx  ! pipe indices at grid points
      real   ,public, dimension(GLOBAL_2D_ARRAY) :: pipe_fraction ! pipe fractional flux at grid points

      real   ,public, dimension(GLOBAL_2D_ARRAY) :: pipe_flx  ! pipe flux
      real   ,public, dimension(npip,N)          :: pipe_prf  ! Pipe vertical profile
      real   ,public, dimension(npip,nt)         :: pipe_trc  ! Pipe tracer conc.

      real   ,public, dimension(npip)            :: pipe_vol  ! Pipe volume

     ! Netcdf names
      character(len=9)  :: pipe_flx_name = 'pipe_flux'  !! stored in the grid file
      character(len=11) :: pipe_vol_name = 'pipe_volume'!! stored in a forcing file
      character(len=11) :: pipe_trc_name = 'pipe_tracer'!! stored in a forcing file
      character(len=9)  :: pipe_tim_name = 'pipe_time'  !! stored in a forcing file
      character(len=5)  ::npipe_dim_name = 'npipe'      !! dimension name for number of pipes in file
      character(len=8)  :: ntrc_dim_name = 'ntracers'   !! dimension name for number of tracers in file
      ! Netcdf variables used for both pipe_vol & pipe_trc
      real, dimension(npip,2)    :: pipe_vol_data ! pvol - pipe data at input times
      real, dimension(npip,nt,2) :: pipe_trc_data ! ptrc - pipe tracer data at input times
      integer, dimension(2) :: p_file_indx = 0
      integer, dimension(2) :: p_irec      = 0
      integer, dimension(2) :: p_it1 = 1, p_it2 = 2
      real, dimension(2,2)  :: p_times = -99 ! [-99,-99]


      public set_pipe_frc_tile

      contains

! ----------------------------------------------------------------------
      subroutine set_pipe_frc_tile(istr,iend,jstr,jend)

      ! SET PIPE FORCES (REALISTIC OR ANALYTICAL FORCING):
      ! - read and interpolation all pipe forcing.
      ! - All pipe variables need time interpolation only
      !   here so can use same generic routine.
      ! - Input data in days!

      implicit none

      ! input/outputs
      integer,intent(in) :: istr,iend,jstr,jend
      ! local
      integer :: i, j

      ! Initialize pipe flux locations (One time only)
      if(FIRST_TIME_STEP) call init_pipe_frc(istr,iend,jstr,jend)

      ! Set pipe flux volumes and tracer data
      if(p_analytical) then
!       if(mynode==0) print*, 'pipes: going to set_ana_pipe_frc'

        call set_ana_pipe_frc(istr,iend,jstr,jend)

        ! pipe vertical profile (pipe_prf) is currently set in analytical.F
        
      else

        ! Set pipe volume for all pipes at current time
        call set_small_arrays(
     &          pipe_vol_name,  pipe_tim_name, ! Text names
     &          pipe_vol_data,  pipe_vol,      ! Variable arrays
     &          p_times(:,1),                  ! Input times
     &          npip,           1,             ! dim1 and dim2=1 if 1D
     &          p_file_indx(1), p_irec(1),     ! File indx & Current input record
     &          p_it1(1),       p_it2(1) )     ! Time index placeholders                   ! Online spatial interpolation flag

!        if(mynode==0) print*, 'pipe_vol= ', pipe_vol ! debug

        ! Set pipe tracers conc. for all pipes at current time
        call set_small_arrays(
     &          pipe_trc_name,  pipe_tim_name, ! Text names
     &          pipe_trc_data,  pipe_trc,      ! Variable arrays
     &          p_times(:,2),                  ! Input times
     &          npip,           nt,            ! dim1 and dim2
     &          p_file_indx(2), p_irec(2),     ! File indx & Current input record
     &          p_it1(2),       p_it2(2) )     ! Time index placeholders

!        if(mynode==0) print*, 'pipe_trc(1,1)= ', pipe_trc(1,1) ! debug
!        if(mynode==0) print*, 'pipe_trc(1,2)= ', pipe_trc(1,2) ! debug

        ! Set pipe vertical profiles
        call set_pipe_vert_prf

        ! Now calculated fluxes with time interpolated volumes and conc.:
        do j = jstr,jend
          do i = istr,iend
            pipe_flx(i,j)=pipe_fraction(i,j)*pipe_vol( pipe_idx(i,j) )
          enddo
        enddo

      endif

      end subroutine set_pipe_frc_tile
! ----------------------------------------------------------------------

      subroutine init_pipe_frc(istr,iend,jstr,jend)

      ! Initialize pipe forcing:
      ! Read in a grid file with locations of pipes and flux contribution per cell.
      ! Done only once as pipe position does not change.
      ! Realistic case - stored as one value in NetCDF file where
      ! pipe grid point value = pidx + pipe_fraction

      use netcdf
      implicit none

      ! inputs
      integer,intent(in) :: istr,iend,jstr,jend
      ! local
      integer :: ierr,ncid,v_id,i,j

      ! Should put this in wrapper in read_write.F as 'nc_read_once'
      if (p_analytical) then

        ! pipe_flx is defined in ana_grid()

      else

        ! Read 'pipe_flux' from grid file (Pipe_idx & pipe_fraction in one value)
        ierr=nf90_open(grdname, nf90_nowrite, ncid)
        ierr=nf90_inq_varid(ncid, pipe_flx_name, v_id)
        ! Temporarily store as pipe_fraction to avoid extra array,
        ! but value still pipe_idx + pipe_fraction
        ierr=ncdf_read_mod(ncid, v_id, 0, rp_var, pipe_fraction(START_2D_ARRAY), 1)
        if(ierr/=0) call handle_ierr(ierr,'init_pipe_frc:: pipe not in file!')

        ! Separate pipe_idx & pipe_fraction:
        do j = jstr,jend
          do i = istr,iend
            ! read in value = pidx + pipe_fraction = pipe_idx(i,j) + pipe_fraction(i,j)
            pipe_idx(i,j)=floor(pipe_fraction(i,j)-1e-5)
            pipe_fraction(i,j)=pipe_fraction(i,j)-pipe_idx(i,j)
          enddo
        enddo

      endif

      if(mynode==0) print*, 'pipes: init pipe locations'
          
      end subroutine init_pipe_frc

! ----------------------------------------------------------------------

      subroutine set_pipe_vert_prf
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

      end subroutine set_pipe_vert_prf

! ----------------------------------------------------------------------
# if defined PIPE_SOURCE
      subroutine set_ana_pipe_frc(istr,iend,jstr,jend)
      ! Analytical pipe forcing volume and tracer data
      ! Put here to avoid circular reference if in analytical.F

#include "pipe_frc_ana.opt"

      end subroutine set_ana_pipe_frc
#endif /* PIPE_SOURCE */
! ----------------------------------------------------------------------

# endif /* ifdef PIPE_SOURCE */
      end module pipe_frc
