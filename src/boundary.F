      module boundary

      ! boundary forcing module

#include "cppdefs.opt"

      use param
      use roms_read_write
      use dimensions
      use hidden_mpi_vars
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

      public :: set_bry_all

      type (ncforce) :: nc_z_w = ncforce( vname='zeta_west', tname='bry_time')
      type (ncforce) :: nc_ub_w = ncforce( vname='ubar_west', tname='bry_time')
      type (ncforce) :: nc_vb_w = ncforce( vname='vbar_west', tname='bry_time')
      type (ncforce) :: nc_u_w = ncforce( vname='u_west', tname='bry_time')
      type (ncforce) :: nc_v_w = ncforce( vname='v_west', tname='bry_time')
      type (ncforce), dimension(nt) :: nc_t_w = ncforce( vname=' ', tname='bry_time')

      type (ncforce) :: nc_z_e = ncforce( vname='zeta_east', tname='bry_time')
      type (ncforce) :: nc_ub_e = ncforce( vname='ubar_east', tname='bry_time')
      type (ncforce) :: nc_vb_e = ncforce( vname='vbar_east', tname='bry_time')
      type (ncforce) :: nc_u_e = ncforce( vname='u_east', tname='bry_time')
      type (ncforce) :: nc_v_e = ncforce( vname='v_east', tname='bry_time')
      type (ncforce), dimension(nt) :: nc_t_e = ncforce( vname=' ', tname='bry_time')

      type (ncforce) :: nc_z_s = ncforce( vname='zeta_south', tname='bry_time')
      type (ncforce) :: nc_ub_s = ncforce( vname='ubar_south', tname='bry_time')
      type (ncforce) :: nc_vb_s = ncforce( vname='vbar_south', tname='bry_time')
      type (ncforce) :: nc_u_s = ncforce( vname='u_south', tname='bry_time')
      type (ncforce) :: nc_v_s = ncforce( vname='v_south', tname='bry_time')
      type (ncforce), dimension(nt) :: nc_t_s = ncforce( vname=' ', tname='bry_time')

      type (ncforce) :: nc_z_n = ncforce( vname='zeta_north', tname='bry_time')
      type (ncforce) :: nc_ub_n = ncforce( vname='ubar_north', tname='bry_time')
      type (ncforce) :: nc_vb_n = ncforce( vname='vbar_north', tname='bry_time')
      type (ncforce) :: nc_u_n = ncforce( vname='u_north', tname='bry_time')
      type (ncforce) :: nc_v_n = ncforce( vname='v_north', tname='bry_time')
      type (ncforce), dimension(nt) :: nc_t_n = ncforce( vname=' ', tname='bry_time')

#ifdef NHMG
      type (ncforce) :: nc_w_w = ncforce( vname='w_west', tname='bry_time')
      type (ncforce) :: nc_w_e = ncforce( vname='w_east', tname='bry_time')
      type (ncforce) :: nc_w_s = ncforce( vname='w_south', tname='bry_time')
      type (ncforce) :: nc_w_n = ncforce( vname='w_north', tname='bry_time')
#endif

      contains

!----------------------------------------------------------------------
      subroutine init_arrays_boundary  ![
      use tracers, only: t_vname
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
      allocate( t_west(0:Mm+1,N,NT) )
      do itrc=1,nt
        nc_t_w(itrc)%vname = trim(t_vname(itrc)) / / '_west'
        allocate( nc_t_w(itrc)%vdata(j0:j1,nz,2) )
      enddo
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
      allocate( t_east(0:Mm+1,N,NT) )
      do itrc=1,nt
        nc_t_e(itrc)%vname = trim(t_vname(itrc)) / / '_east'
        allocate( nc_t_e(itrc)%vdata(j0:j1,nz,2) )
      enddo
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
      allocate( t_south(0:Lm+1,N,NT) )
      do itrc=1,nt
        nc_t_s(itrc)%vname = trim(t_vname(itrc)) / / '_south'
        allocate( nc_t_s(itrc)%vdata(i0:i1,nz,2) )
      enddo
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
      allocate( t_north(0:Lm+1,N,NT) )
      do itrc=1,nt
        nc_t_n(itrc)%vname = trim(t_vname(itrc)) / / '_north'
        allocate( nc_t_n(itrc)%vdata(i0:i1,nz,2) )
      enddo
#   endif
#  endif
# endif

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
#  ifdef Z_FRC_BRY
        call set_frc_data( nc_z_w, zeta_west(j0:j1), obry=2 )
!       call exch_bry_EW_tile(1, ny, zeta_west, 1)      ! DevinD - exchange did not make a difference
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
!       call set_frc_data( nc_w_w, var2d=w_west(j0:j1,:), obry=1 )
        w_west = 0
#     endif
#   endif
#   ifdef T_FRC_BRY
        do itrc=1,nt
          call set_frc_data( nc_t_w(itrc), var2d=t_west(j0:j1,:,itrc), obry=1 )
        enddo
#   endif
#  endif
      endif

      if (obc_east.and.(.not.east_exchng)) then
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
!       call set_frc_data( nc_w_e, var2d=w_east(j0:j1,:), obry=1 )
        w_east = 0
#     endif
#   endif
#   ifdef T_FRC_BRY
        do itrc=1,NT
          call set_frc_data( nc_t_e(itrc), var2d=t_east(j0:j1,:,itrc), obry=1 )
        enddo
#   endif
#  endif
      endif

      if (obc_south.and.(.not.south_exchng)) then
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
!       call set_frc_data( nc_w_s, var2d=w_south(i0:i1,:), obry=1 )
        w_south = 0
#     endif
#   endif
#   ifdef T_FRC_BRY
        do itrc=1,nt
          call set_frc_data( nc_t_s(itrc), var2d=t_south(i0:i1,:,itrc), obry=1 )
        enddo
#   endif
#  endif
      endif

      if (obc_north.and.(.not.north_exchng)) then
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
!       call set_frc_data( nc_w_n, var2d=w_north(i0:i1,:), obry=1 )
        w_north = 0
#     endif
#   endif
#   ifdef T_FRC_BRY
        do itrc=1,NT
          call set_frc_data( nc_t_n(itrc), var2d=t_north(i0:i1,:,itrc), obry=1 )
        enddo
#   endif
#  endif
      endif

      end subroutine set_bry_all  !]

!----------------------------------------------------------------------

      end module boundary
