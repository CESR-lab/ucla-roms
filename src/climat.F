      module climat

#include "cppdefs.opt"

      use param
      implicit none

! This is include file "climat"
!----- -- ------- ---- ----------
! Free surface climatology:
! ==== ======= ============
!   ssh        sea surface height climatology at current time-step.
!   Znudgcof   inverse relaxation time [1/sec] for nudging toward
!                               free surface climatological fields.
!   sshg       two-time-level array to hold climatological data for
!                                                     free surface.
!   tssh       time of read in sea surface height climatology.

#if defined M2NUDGING && !defined M2_FRC_BRY
      real ssh(GLOBAL_2D_ARRAY)
# ifndef ANA_SSH
#  if defined SSH_DATA || defined ALL_DATA
      real sshg(GLOBAL_2D_ARRAY,2)
 
      real(kind=8) ssh_cycle, ssh_time(2)
      integer itssh, ntssh, ssh_ncycle, ssh_rec, ssh_id, ssh_tid
 
#   undef SSH_DATA
#  endif /* SSH_DATA */
# endif /* !ANA_SSH */
#endif
!
! Temperature, salinity climatology:
!------------- -------- ------------
!   tclm       climatology for tracer variables at current time-step.
!   Tnudgcof   inverse relaxation time [1/sec] for nudging toward
!                                       tracer climatological fields.
!   tclima     two-time-level array to hold climatological data for
!                                               tracer variables.
!   ttclm      time of read in climatology for tracer type variables.
!
! WARNING: CPP-switch TNUDGING introduces forcing at the boundary
! only, while TCLIMATOLOGY activates spatially nonuniform nudging
! inside the domain using Tnudgcof specified in set_nudgcof.
!
#ifdef SOLVE3D
# if defined TCLIMATOLOGY || defined TNUDGING 
      real tclm(GLOBAL_2D_ARRAY,N,NT)
#  ifdef TCLIMATOLOGY 
      real Tnudgcof(GLOBAL_2D_ARRAY,NT)
#  endif
#  ifndef ANA_TCLIMA
#   if defined TCLIMA_DATA || defined ALL_DATA
      real tclima(GLOBAL_2D_ARRAY,N,2,NT)
 
      real(kind=8) tclm_cycle(NT), tclm_time(2,NT)
      integer tclm_ncycle(NT), tclm_rec(NT), ittclm(NT),
     &        nttclm(NT),      tclm_tid(NT), tclm_id(NT)
#    undef TCLIMA_DATA
#   endif
#  endif
# endif
#endif
!
! barotropic and baroclinic velocity climatology:
! ========== === ========== ======== ===========
!   ubclm     climatology for bar. u-velocity at current time-step.
!   vbclm     climatology for bar. v-velocity at current time-step.
!   uclm      climatology for u-velocity at current time-step.
!   vclm      climatology for v-velocity at current time-step.
! 
!   ubclima   two-time-level array to hold climatological data
!   vbclima
!   uclima
!   vclima
!
#ifdef UCLIMATOLOGY
      real ubclm(GLOBAL_2D_ARRAY)
      real vbclm(GLOBAL_2D_ARRAY)
      real ubclima(GLOBAL_2D_ARRAY,2)
      real vbclima(GLOBAL_2D_ARRAY,2)
 
# ifdef SOLVE3D
      real uclm(GLOBAL_2D_ARRAY,N)
      real vclm(GLOBAL_2D_ARRAY,N)
      real uclima(GLOBAL_2D_ARRAY,N,2)
      real vclima(GLOBAL_2D_ARRAY,N,2)
# endif
 
      real(kind=8) uclm_cycle,  uclm_time(2)
      integer uclm_ncycle,   uclm_rec,  ituclm,  ntuclm,
     &   uclm_tid, ubclm_id, vbclm_id,  uclm_id, vclm_id
 
#endif /* UCLIMATOLOGY */

      end module climat
