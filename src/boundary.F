      module boundary

      ! boundary forcing module

#include "cppdefs.opt"

      use param
      use roms_read_write
      use dimensions
      use roms_mpi
      use tracers
      implicit none

      ! Set module to private and make vars explicitly public

      real,allocatable,dimension(:)     :: zeta_west
      real,allocatable,dimension(:)     :: ubar_west, vbar_west
      real,allocatable,dimension(:,:)   :: u_west, v_west, w_west
      real,allocatable,dimension(:,:,:) :: t_west

      real,allocatable,dimension(:)     :: zeta_east
      real,allocatable,dimension(:)     :: ubar_east, vbar_east
      real,allocatable,dimension(:,:)   :: u_east, v_east, w_east
      real,allocatable,dimension(:,:,:) :: t_east

      real,allocatable,dimension(:)     :: zeta_south
      real,allocatable,dimension(:)     :: ubar_south, vbar_south
      real,allocatable,dimension(:,:)   :: u_south, v_south, w_south
      real,allocatable,dimension(:,:,:) :: t_south

      real,allocatable,dimension(:)     :: zeta_north
      real,allocatable,dimension(:)     :: ubar_north, vbar_north
      real,allocatable,dimension(:,:)   :: u_north, v_north, w_north
      real,allocatable,dimension(:,:,:) :: t_north

      type (ncforce) :: nc_z_w = ncforce( vname='zeta_west', tname='bry_time')
      type (ncforce) :: nc_ub_w = ncforce( vname='ubar_west', tname='bry_time')
      type (ncforce) :: nc_vb_w = ncforce( vname='vbar_west', tname='bry_time')
      type (ncforce) :: nc_u_w = ncforce( vname='u_west', tname='bry_time')
      type (ncforce) :: nc_v_w = ncforce( vname='v_west', tname='bry_time')
      type (ncforce),dimension(nt) :: nc_t_w=ncforce(vname=' ',tname='bry_time')

      type (ncforce) :: nc_z_e = ncforce( vname='zeta_east', tname='bry_time')
      type (ncforce) :: nc_ub_e = ncforce( vname='ubar_east', tname='bry_time')
      type (ncforce) :: nc_vb_e = ncforce( vname='vbar_east', tname='bry_time')
      type (ncforce) :: nc_u_e = ncforce( vname='u_east', tname='bry_time')
      type (ncforce) :: nc_v_e = ncforce( vname='v_east', tname='bry_time')
      type (ncforce),dimension(nt) :: nc_t_e=ncforce(vname=' ',tname='bry_time')

      type (ncforce) :: nc_z_s = ncforce( vname='zeta_south', tname='bry_time')
      type (ncforce) :: nc_ub_s = ncforce( vname='ubar_south', tname='bry_time')
      type (ncforce) :: nc_vb_s = ncforce( vname='vbar_south', tname='bry_time')
      type (ncforce) :: nc_u_s = ncforce( vname='u_south', tname='bry_time')
      type (ncforce) :: nc_v_s = ncforce( vname='v_south', tname='bry_time')
      type (ncforce),dimension(nt) :: nc_t_s=ncforce(vname=' ',tname='bry_time')

      type (ncforce) :: nc_z_n = ncforce( vname='zeta_north', tname='bry_time')
      type (ncforce) :: nc_ub_n = ncforce( vname='ubar_north', tname='bry_time')
      type (ncforce) :: nc_vb_n = ncforce( vname='vbar_north', tname='bry_time')
      type (ncforce) :: nc_u_n = ncforce( vname='u_north', tname='bry_time')
      type (ncforce) :: nc_v_n = ncforce( vname='v_north', tname='bry_time')
      type (ncforce),dimension(nt) :: nc_t_n=ncforce(vname=' ',tname='bry_time')

#ifdef NHMG
      type (ncforce) :: nc_w_w = ncforce( vname='w_west', tname='bry_time')
      type (ncforce) :: nc_w_e = ncforce( vname='w_east', tname='bry_time')
      type (ncforce) :: nc_w_s = ncforce( vname='w_south', tname='bry_time')
      type (ncforce) :: nc_w_n = ncforce( vname='w_north', tname='bry_time')
      logical,dimension(4) :: w_bry_dat = .false.  ! w_west,east,south,north found in bry files
#endif

      public :: set_bry_all

      contains

!----------------------------------------------------------------------
      subroutine init_arrays_boundary  ![
      implicit none

      ! local
      integer :: itrc
      if (obc_west.or.obc_east.or.obc_north.or.obc_south) then
        bc_options = ''
      endif
      if (obc_west)  call store_string_att(bc_options,'OBC_WEST, ')
      if (obc_east)  call store_string_att(bc_options,'OBC_EAST, ')
      if (obc_north) call store_string_att(bc_options,'OBC_NORTH, ')
      if (obc_south) call store_string_att(bc_options,'OBC_SOUTH, ')
      bc_options = trim(adjustl(bc_options))
# ifdef OBC_M3ORLANSKI
      call store_string_att(bc_options,'OBC_M3ORLANSKI, ')
#  endif
# ifdef OBC_M3SPECIFIED
      call store_string_att(bc_options,'OBC_M3SPECIFIED, ')
#  endif
# ifdef OBC_TORLANSKI
      call store_string_att(bc_options,'OBC_TORLANSKI, ')
#  endif
# ifdef OBC_TSPECIFIED
      call store_string_att(bc_options,'OBC_TSPECIFIED, ')
#  endif

      ! To read input data from a file, we need 2 time-slices
      ! (..%vdata) and a data array to which we interpolate in
      ! time.

      if (obc_west) then
        allocate(zeta_west(0:ny+1),nc_z_w%vdata(j0:j1,1,2) )
        allocate(ubar_west(0:ny+1),nc_ub_w%vdata(j0:j1,1,2) )
        allocate(vbar_west(0:ny+1),nc_vb_w%vdata( 1:j1,1,2) )
        allocate(u_west(0:ny+1,nz),nc_u_w%vdata(j0:j1,nz,2) )
        allocate(v_west(0:ny+1,nz),nc_v_w%vdata( 1:j1,nz,2) )
# ifdef NHMG
        allocate( w_west(0:ny+1,nz), nc_w_w%vdata(j0:j1,nz,2) )
        if (inode==0) then
          w_bry_dat(1) = in_frc_files(nc_w_w%vname)
          if (.not.w_bry_dat(1) .and. jnode==0)
     &      write(*,*) 'boundary: w_west not in bry files, set to 0'
        endif
# endif
        allocate( t_west(0:ny+1,nz,nt) )
        do itrc=1,nt
          if (t_bry_dat(itrc)) then
            nc_t_w(itrc)%vname = trim(t_vname(itrc)) // '_west'
            allocate( nc_t_w(itrc)%vdata(j0:j1,nz,2) )
          endif
        enddo
      endif

      if (obc_east) then
        allocate( zeta_east(0:ny+1) )
        allocate( nc_z_e%vdata(j0:j1,1,2) )
        allocate( ubar_east(0:ny+1), vbar_east(0:ny+1) )
        allocate( nc_ub_e%vdata(j0:j1,1,2) )
        allocate( nc_vb_e%vdata( 1:j1,1,2) )
        allocate( u_east(0:ny+1,nz), v_east(0:ny+1,nz) )
        allocate( nc_u_e%vdata(j0:j1,nz,2) )
        allocate( nc_v_e%vdata( 1:j1,nz,2) )
#    ifdef NHMG
        allocate( w_east(0:ny+1,nz), nc_w_e%vdata(j0:j1,nz,2) )
        if (inode==npx-1) then
          w_bry_dat(2) = in_frc_files(nc_w_e%vname)
          if (.not.w_bry_dat(2) .and. jnode==0)
     &      write(*,*) 'boundary: w_east not in bry files, set to 0'
        endif
#    endif
        allocate( t_east(0:ny+1,nz,nt) )
        do itrc=1,nt
          if (t_bry_dat(itrc)) then
            nc_t_e(itrc)%vname = trim(t_vname(itrc)) / / '_east'
            allocate( nc_t_e(itrc)%vdata(j0:j1,nz,2) )
          endif
        enddo
      endif

      if (obc_south) then
        allocate( zeta_south(0:nx+1) )
        allocate( nc_z_s%vdata(i0:i1,1,2) )
        allocate( ubar_south(0:nx+1), vbar_south(0:nx+1) )
        allocate( nc_ub_s%vdata( 1:i1,1,2) )
        allocate( nc_vb_s%vdata(i0:i1,1,2) )
        allocate( u_south(0:nx+1,nz), v_south(0:nx+1,nz) )
        allocate( nc_u_s%vdata( 1:i1,nz,2) )
        allocate( nc_v_s%vdata(i0:i1,nz,2) )
#    ifdef NHMG
        allocate( w_south(0:nx+1,nz), nc_w_s%vdata(i0:i1,nz,2) )
        if (jnode==0) then
          w_bry_dat(3) = in_frc_files(nc_w_s%vname)
          if (.not.w_bry_dat(3) .and. inode==0)
     &      write(*,*) 'boundary: w_south not in bry files, set to 0'
        endif
#    endif
        allocate( t_south(0:nx+1,nz,nt) )
        do itrc=1,nt
          if (t_bry_dat(itrc)) then
            nc_t_s(itrc)%vname = trim(t_vname(itrc)) / / '_south'
            allocate( nc_t_s(itrc)%vdata(i0:i1,nz,2) )
          endif
        enddo
      endif

      if (obc_north) then
        allocate( zeta_north(0:nx+1) )
        allocate( nc_z_n%vdata(i0:i1,1,2) )
        allocate( ubar_north(0:nx+1), vbar_north(0:nx+1) )
        allocate( nc_ub_n%vdata( 1:i1,1,2) )
        allocate( nc_vb_n%vdata(i0:i1,1,2) )
        allocate( u_north(0:nx+1,nz), v_north(0:nx+1,nz) )
        allocate( nc_u_n%vdata( 1:i1,nz,2) )
        allocate( nc_v_n%vdata(i0:i1,nz,2) )
#    ifdef NHMG
        allocate( w_north(0:nx+1,nz), nc_w_n%vdata(i0:i1,nz,2) )
        if (jnode==npy-1) then
          w_bry_dat(4) = in_frc_files(nc_w_n%vname)
          if (.not.w_bry_dat(4) .and. inode==0)
     &      write(*,*) 'boundary: w_north not in bry files, set to 0'
        endif
#    endif
        allocate( t_north(0:nx+1,nz,nt) )
        do itrc=1,nt
          if (t_bry_dat(itrc)) then
            nc_t_n(itrc)%vname = trim(t_vname(itrc)) / / '_north'
            allocate( nc_t_n(itrc)%vdata(i0:i1,nz,2) )
          endif
        enddo
      endif

      end subroutine init_arrays_boundary !]
!----------------------------------------------------------------------
      subroutine set_bry_all  ![
      use param
      use scalars
      implicit none

      integer :: itrc

#  ifdef ANA_BRY
      return
#  endif
      if (obc_west.and.(inode==0)) then
        call set_frc_data( nc_z_w, zeta_west(j0:j1), obry=2 )
        call set_frc_data( nc_ub_w, ubar_west(j0:j1), obry=1 )
        call set_frc_data( nc_vb_w, vbar_west( 1:j1), obry=1 )
        call set_frc_data( nc_u_w, var2d=u_west(j0:j1,:), obry=1 )
        call set_frc_data( nc_v_w, var2d=v_west( 1:j1,:), obry=1 )
#     ifdef NHMG
        if (w_bry_dat(1)) then
          call set_frc_data( nc_w_w, var2d=w_west(j0:j1,:), obry=1 )
        else
          w_west = 0
        endif
#     endif
        do itrc=1,nt
          if (t_bry_dat(itrc)) then
            call set_frc_data( nc_t_w(itrc), var2d=t_west(j0:j1,:,itrc), obry=1 )
          else
            t_west(j0:j1,:,itrc) = t_bry_val(itrc)
          endif
        enddo
      endif

      if (obc_east.and.inode.eq.npx-1) then
        call set_frc_data( nc_z_e, zeta_east(j0:j1), obry=2 )
        call set_frc_data( nc_ub_e, ubar_east(j0:j1), obry=1 )
        call set_frc_data( nc_vb_e, vbar_east( 1:j1), obry=1 )
        call set_frc_data( nc_u_e, var2d=u_east(j0:j1,:), obry=1 )
        call set_frc_data( nc_v_e, var2d=v_east( 1:j1,:), obry=1 )
#     ifdef NHMG
        if (w_bry_dat(2)) then
          call set_frc_data( nc_w_e, var2d=w_east(j0:j1,:), obry=1 )
        else
          w_east = 0
        endif
#     endif
        do itrc=1,nt
          if (t_bry_dat(itrc)) then
            call set_frc_data( nc_t_e(itrc), var2d=t_east(j0:j1,:,itrc), obry=1 )
          else
            t_east(j0:j1,:,itrc) = t_bry_val(itrc)
          endif
        enddo
      endif

      if (obc_south.and.jnode.eq.0) then
        call set_frc_data( nc_z_s, zeta_south(i0:i1), obry=2 )
        call set_frc_data( nc_ub_s, ubar_south( 1:i1), obry=1 )
        call set_frc_data( nc_vb_s, vbar_south(i0:i1), obry=1 )
        call set_frc_data( nc_u_s, var2d=u_south( 1:i1,:), obry=1 )
        call set_frc_data( nc_v_s, var2d=v_south(i0:i1,:), obry=1 )
#     ifdef NHMG
        if (w_bry_dat(3)) then
          call set_frc_data( nc_w_s, var2d=w_south(i0:i1,:), obry=1 )
        else
          w_south = 0
        endif
#     endif
        do itrc=1,nt
          if (t_bry_dat(itrc)) then
            call set_frc_data( nc_t_s(itrc), var2d=t_south(i0:i1,:,itrc), obry=1 )
          else
            t_south(i0:i1,:,itrc) = t_bry_val(itrc)
          endif
        enddo
      endif

      if (obc_north.and.jnode.eq.npy-1) then
        call set_frc_data( nc_z_n, zeta_north(i0:i1), obry=2 )
        call set_frc_data( nc_ub_n, ubar_north( 1:i1), obry=1 )
        call set_frc_data( nc_vb_n, vbar_north(i0:i1), obry=1 )
        call set_frc_data( nc_u_n, var2d=u_north( 1:i1,:), obry=1 )
        call set_frc_data( nc_v_n, var2d=v_north(i0:i1,:), obry=1 )
#     ifdef NHMG
        if (w_bry_dat(4)) then
          call set_frc_data( nc_w_n, var2d=w_north(i0:i1,:), obry=1 )
        else
          w_north = 0
        endif
#     endif
        do itrc=1,nt
          if (t_bry_dat(itrc)) then
            call set_frc_data( nc_t_n(itrc), var2d=t_north(i0:i1,:,itrc), obry=1 )
          else
            t_north(i0:i1,:,itrc) = t_bry_val(itrc)
          endif
        enddo
      endif

      end subroutine set_bry_all  !]
!----------------------------------------------------------------------
      logical function in_frc_files(vname) ![
      ! True if vname is found in any of the forcing files
      use netcdf
      implicit none
      character(len=*),intent(in) :: vname
      integer :: ifile,ncid,vid,ierr

      in_frc_files = .false.
      do ifile=1,max_frc
        ierr = nf90_open(frcfile(ifile),nf90_nowrite,ncid)
        if (ierr/=nf90_noerr) cycle
        ierr = nf90_inq_varid(ncid,vname,vid)
        in_frc_files = (ierr==nf90_noerr)
        ierr = nf90_close(ncid)
        if (in_frc_files) return
      enddo

      end function in_frc_files !]

!----------------------------------------------------------------------

      end module boundary
