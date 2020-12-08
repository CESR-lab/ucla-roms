#include "cppdefs.h"
      module diagnostics
#if defined DIAGNOSTICS
  
       use dimensions !! has dimensions, a list of forcing files, and rmask
       implicit none

       private
       ! Udiag(1): du/dt
       ! Mdiag(2): hydrostatic pressure gradient
       ! Udiag(3): cor(u) 
       ! Udiag(4): adv(u) (non-dissipative part)
       ! Udiag(5): horizontal mixing (including numerical dissipation)
       ! Udiag(6): vertical mixing   (including numerical dissipation)
       ! Udiag(7): non-hydrostatic pressure gradient

       ! Tdiag(1): dt/dt
       ! Tdiag(2): adv(t) (non-dissipative part)
       ! Tdiag(3): adv(t) (dissipative part)
       ! Tdiag(4): vertical mixing
       ! Tdiag(5): sources/sinks

       logical,public :: diag_time = .true.
       logical,public :: do_diags


       integer,public,parameter :: ndm = 6  !! number of momentum diag terms
       integer,public,parameter :: ndt = 5  !! number of tracer diag terms

       real,public, dimension(nx,ny,nz,3,ndm)  :: Mdiag  !! u,v,w

       real,public, dimension(nx,ny,nz,ndt,nt) :: Tdiag
       

       ! These are private to this module
       ! Names of diagnostic variables

       character(len=5) :: udiag_name = 'udiag' !! 
       character(len=5) :: vdiag_name = 'vdiag' !! 
       character(len=5) :: wdiag_name = 'wdiag' !! 
       character(len=5) :: tdiag_name = 'tdiag' !! 

       public write_diagnostics

      contains
! ----------------------------------------------------------------------
        subroutine write_diagnostics
         use netcdf
         implicit none

        end subroutine write_diagnostics
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
        subroutine create_diagnostics_file
         use netcdf
         implicit none



        end subroutine create_diagnostics_file
       ! --------------------------------------------------------------------
#endif
      end module diagnostics
