! This module created to avoid circular dependencies with the MARBL BGC driver
      module get_bulk_wnd
#include "cppdefs.opt"      

      use param
!      use read_write, only: set_frc_var_tile, ncvs
      use roms_read_write
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

      allocate( nc_uwnd%vdata( GLOBAL_2D_ARRAY,2) )
      allocate( nc_vwnd%vdata( GLOBAL_2D_ARRAY,2) )

      if (interp_frc==1) then 
        nc_uwnd%coarse = 1
        nc_vwnd%coarse = 1
      endif

      
      end subroutine init_arrays_bulk_wnd
! --------------------------------------------------------------------------------
      subroutine set_bulk_wnd_tile(istr,iend,jstr,jend)

      implicit none
      integer, intent(in) :: istr,iend,jstr,jend
      
      ! Wind (get_bulk_wnd.F & set_bulk_wnd() in old code)
      !call set_frc_var_tile( istr, iend, jstr, jend, nc_uwnd,  uwnd,  interp_frc ) ! u-wind
!call set_frc_var_tile( istr, iend, jstr, jend, nc_vwnd,  vwnd,  interp_frc ) ! v-wind
      call set_frc_data(nc_uwnd,uwnd,'r')
      call set_frc_data(nc_vwnd,vwnd,'r')

     
      end subroutine set_bulk_wnd_tile

      end module get_bulk_wnd

