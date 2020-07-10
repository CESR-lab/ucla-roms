#include "cppdefs.h"
      module diagnostics
#if defined DIAGNOSTICS
  
       use dimensions !! has dimensions, a list of forcing files, and rmask
       implicit none

       private
       ! Udiag(1): du/dt
       ! Udiag(2): adv(u) (non-dissipative part)
       ! Udiag(3): adv(u) (dissipative part)
       ! Udiag(4): cor(u) 
       ! Udiag(5): dp_hydro
       ! Udiag(6): vertical mixing
       ! Udiag(7): dp_nonhydro

       ! Tdiag(1): dt/dt
       ! Tdiag(2): adv(t) (non-dissipative part)
       ! Tdiag(3): adv(t) (dissipative part)
       ! Tdiag(4): vertical mixing
       ! Tdiag(5): sources/sinks

       integer,public,parameter, ndt = 5  !! number of tracer diag terms
# ifdef NHMG
       integer,public, parameter, ndu = 7 !! number of momentum diag terms
       real,public, dimension(nx,ny,nz+1,ndu)  :: Wdiag
# else
       integer,public, parameter, ndu = 6
# endif
       real,public, dimension(nx+1,ny,nz,ndu)  :: Udiag
       real,public, dimension(nx,ny+1,nz,ndu)  :: Vdiag
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
