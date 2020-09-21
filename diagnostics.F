#include "cppdefs.h"
      module diagnostics
#if defined DIAGNOSTICS
  
       use dimensions !! has dimensions, a list of forcing files, and rmask
       implicit none

       private
       ! Udiag(1): du/dt
       ! Udiag(2): hydrostatic pressure gradient (includes WEC terms)
       ! Udiag(3): cor(u)  + curvilinear
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

# ifdef NHMG
       real,public, dimension(nx,ny,nz,3,7) :: Mdiag
# else
       real,public, dimension(nx,ny,nz,2,6) :: Mdiag
# endif
       integer,public,parameter :: ndt = 5  !! number of tracer diag terms
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

         write(20,999) Mdiag(10,10,10,1,:)
 999  format(7e12.4)

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
