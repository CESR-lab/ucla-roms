      module bgc_shared_vars
      ! formerly bgc_ecosys_bec2.h

c # TODO marbldrve_configure_diagnostics for pointers example

#include "cppdefs.opt"
#if defined(BIOLOGY_BEC2) || defined(MARBL)

      use param
      use tracers               ! for iPO4, etc, indices of bgc tracers (formerly in param.h)
      use roms_read_write
#if defined(BIOLOGY_BEC2)
      use bec2_vars
#elif defined(MARBL)
      use marbl_driver,only: marbldrv_configure_saved_state
# ifdef MARBL_DIAGS
     &     ,marbldrv_configure_diagnostics
     &     ,nr_marbl_diag_2d,nr_marbl_diag_3d
     &     ,marbl_diag_2d,marbl_diag_3d
     &     ,idx_marbl_diag_2d,idx_marbl_diag_3d
     &     ,wrt_marbl_diag_2d,wrt_marbl_diag_3d
# endif
#endif

      implicit none

#include "bgc.opt"


#if defined(BEC2_DIAG) || defined(MARBL_DIAGS)
!
! Diagnostic variables appearing in average and history files:
!
      integer :: nr_bgc_wrdiag_2d, nr_bgc_wrdiag_3d, nr_bgc_wrdiag
#if defined(BEC2_DIAG)
      integer, parameter :: nr_bgc_diag_2d = nr_bec2_diag_2d
      integer, parameter :: nr_bgc_diag_3d = nr_bec2_diag_3d
#elif defined(MARBL_DIAGS)
      ! Just hardcoding this for now but would be nice to obtain it dynamically
      integer, parameter ::  nr_bgc_diag_3d = nr_marbl_diag_3d
      integer, parameter ::  nr_bgc_diag_2d = nr_marbl_diag_2d
#endif
      integer, parameter :: nr_bgc_diag = nr_bgc_diag_2d + nr_bgc_diag_3d

      real,allocatable,dimension(:,:,:,:),target   :: bgc_diag_3d
      real,allocatable,dimension(:,:,:)  ,target   :: bgc_diag_2d

! Control BEC2_DIAG Output vars
      logical, target :: wrt_bgc_diag_2d(nr_bgc_diag_2d)
      logical, target :: wrt_bgc_diag_3d(nr_bgc_diag_3d)
      integer, target :: idx_bgc_diag_2d(nr_bgc_diag_2d)
      integer, target :: idx_bgc_diag_3d(nr_bgc_diag_3d)
      ! Arrays storing information (name, unit, fill value) about each diagnostic variable:
      character*72, target ::  vname_bgc_diag_2d(4,nr_bgc_diag_2d)
      character*72, target ::  vname_bgc_diag_3d(4,nr_bgc_diag_3d)

      real,allocatable,dimension(:,:,:,:),target :: bgc_diag_3d_avg
      real,allocatable,dimension(:,:,:)  ,target :: bgc_diag_2d_avg

#endif /* BEC2_DIAG  || MARBL_DIAGS */

      contains

#if defined(BEC2_DIAG) || defined(MARBL_DIAGS)
      subroutine find_write_bgc_diag
      implicit none

      integer idiag
#if defined(MARBL_DIAGS)

      wrt_marbl_diag_2d => wrt_bgc_diag_2d
      wrt_marbl_diag_3d => wrt_bgc_diag_3d
      idx_marbl_diag_2d => idx_bgc_diag_2d
      idx_marbl_diag_3d => idx_bgc_diag_3d


      call marbldrv_configure_diagnostics(
     &     vname_bgc_diag_2d,vname_bgc_diag_3d
     & )

#elif defined(BEC2_DIAG)

      wrt_bec2_diag_2d => wrt_bgc_diag_2d
      wrt_bec2_diag_3d => wrt_bgc_diag_3d
      idx_bec2_diag_2d => idx_bgc_diag_2d
      idx_bec2_diag_3d => idx_bgc_diag_3d
      vname_bec2_diag_2d => vname_bgc_diag_2d
      vname_bec2_diag_3d => vname_bgc_diag_3d

#include "bgc_2Ddiagnostics.opt"
#include "bgc_3Ddiagnostics.opt"

#endif /* MARBL_DIAGS */


      nr_bgc_wrdiag_2d=0
      do idiag=1,nr_bgc_diag_2d
        if (wrt_bgc_diag_2d(idiag)) then
           nr_bgc_wrdiag_2d=nr_bgc_wrdiag_2d+1
           idx_bgc_diag_2d(idiag)=nr_bgc_wrdiag_2d
           if (mynode == 0) then
              write(*,'(7x,A,I3,I3,A,A)') '2D BGC diags :: ',
     &        nr_bgc_wrdiag_2d, idiag, ' ', vname_bgc_diag_2d(1,idiag)
           endif
        else
           idx_bgc_diag_2d(idiag)=0
        endif
      enddo

      nr_bgc_wrdiag_3d=0
      do idiag=1,nr_bgc_diag_3d
        if (wrt_bgc_diag_3d(idiag)) then
           nr_bgc_wrdiag_3d=nr_bgc_wrdiag_3d+1
           idx_bgc_diag_3d(idiag)=nr_bgc_wrdiag_3d
           if (mynode == 0) then
              write(*,'(7x,A,I3,I3,A,A)') '3D BGC diags :: ',
     &        nr_bgc_wrdiag_3d, idiag, ' ',vname_bgc_diag_3d(1,idiag)
           endif
        else
           idx_bgc_diag_3d(idiag)=0
        endif
      enddo


      end subroutine find_write_bgc_diag  !]
#endif /* BEC2_DIAG || MARBL_DIAGS */


      subroutine init_arrays_bgc_shared_vars
      implicit none

#if defined(BEC2_DIAG) || defined (MARBL_DIAGS)

      call find_write_bgc_diag

      allocate( bgc_diag_3d(GLOBAL_2D_ARRAY,N,nr_bgc_wrdiag_3d) )
      allocate( bgc_diag_2d(GLOBAL_2D_ARRAY,nr_bgc_wrdiag_2d) )

      bgc_diag_2d=0.0
      bgc_diag_3d=0.0
#if defined(MARBL_DIAGS)
      marbl_diag_2d => bgc_diag_2d
      marbl_diag_3d => bgc_diag_3d
#elif defined(BEC2_DIAG)
      bec2_diag_2d => bgc_diag_2d
      bec2_diag_3d => bgc_diag_3d
#endif

      if (wrt_avg_dia) then
        allocate( bgc_diag_3d_avg(GLOBAL_2D_ARRAY,N,nr_bgc_wrdiag_3d) )
        allocate( bgc_diag_2d_avg(GLOBAL_2D_ARRAY,nr_bgc_wrdiag_2d) )
        bgc_diag_2d_avg=0.0
        bgc_diag_3d_avg=0.0
      endif

      if (mynode == 0) then
              write(*,'(7x,A,I3,I3)') 'BGC diags allocation :: ',
     &        nr_bgc_wrdiag_2d, nr_bgc_wrdiag_3d
      endif

#endif /* BEC2_DIAG || MARBL_DIAGS */

#if defined(MARBL)
!     Allocate arrays for MARBL saved state variables
      call marbldrv_configure_saved_state
#elif defined(BIOLOGY_BEC2)
      call init_arrays_bec2_vars ! bec2 version
#endif

      end subroutine init_arrays_bgc_shared_vars

#endif  /* BIOLOGY_BEC2 || MARBL  */
      end module bgc_shared_vars
