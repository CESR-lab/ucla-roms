! This module created to avoid circular dependencies with the MARBL BGC driver
      module get_bulk_wnd
#include "cppdefs.opt"      
#if defined BULK_FRC
      use param
      use read_write, only: set_frc_var_tile, ncvs
      
      implicit none
!      private
      
! next line borrowed straight from bulk_frc.opt      
#include "bulk_wnd.opt"  

      real,allocatable,dimension(:,:) ::  uwnd             ! time interpolated u-wind
      real,allocatable,dimension(:,:) ::  vwnd             ! time interpolated v-wind

      public init_arrays_bulk_wnd
      public set_bulk_wnd_tile
      
      contains

! --------------------------------------------------------------------------------      
      subroutine init_arrays_bulk_wnd
      implicit none

      allocate(  uwnd(GLOBAL_2D_ARRAY) )
      allocate(  vwnd(GLOBAL_2D_ARRAY) )

      allocate( nc_uwnd%data( GLOBAL_2D_ARRAY,2) )
      allocate( nc_vwnd%data( GLOBAL_2D_ARRAY,2) )

      end subroutine init_arrays_bulk_wnd
! --------------------------------------------------------------------------------
      subroutine set_bulk_wnd_tile(istr,iend,jstr,jend)

      implicit none
      integer, intent(in) :: istr,iend,jstr,jend
      
      ! Wind (get_bulk_wnd.F & set_bulk_wnd() in old code)
      call set_frc_var_tile( istr, iend, jstr, jend, nc_uwnd,  uwnd,  interp_frc )  ! u-wind
      call set_frc_var_tile( istr, iend, jstr, jend, nc_vwnd,  vwnd,  interp_frc )  ! v-wind
     
      end subroutine set_bulk_wnd_tile
#endif /* defined BULK_FRC */      
      end module get_bulk_wnd

