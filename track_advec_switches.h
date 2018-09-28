! This part tracks the status of CPP-switches related to advection.
! Because some of these switches may be configured inside the include
! files where they are used, while others may be configured outside
! (usually in the file where the include file is included), the same
! include file may yield different numerical algorithms depending on
! where it is used and the CPP-settings defined there. Therefore, this
! part of the tracking algorithm is used in three places: pre_step,
! step3d_uv1, step3d_t in somewhat redundant way to detect all the
! possible algorithmic variations.

! Switches controlling advection terms in momentum equation and the
! associated top and bottom vertical boundary conditions appearing in
! "compute_horiz_rhs_uv_terms.h" and "compute_vert_rhs_uv_terms.h"

# ifdef UPSTREAM_UV
      is=ie+2 ; ie=is+10
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='UPSTREAM_UV'
# endif
# ifdef SPLINE_UV
      is=ie+2 ; ie=is+8
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='SPLINE_UV'
# endif
# ifdef NEUMANN_UV
      is=ie+2 ; ie=is+9
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='NEUMANN_UV'
# endif

! Switches controlling tracer advection and the associated vertical
! top/bottom boundary conditions inside "compute_horiz_tracer_fluxes.h"
! and "compute_vert_tracer_fluxes.h"

# ifdef AKIMA
      is=ie+2 ; ie=is+4
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='AKIMA'
# endif
# ifdef UPSTREAM_TS
      is=ie+2 ; ie=is+10
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='UPSTREAM_TS'
# endif
# ifdef AKIMA_V
      is=ie+2 ; ie=is+6
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='AKIMA_V'
# endif
# ifdef SPLINE_TS
      is=ie+2 ; ie=is+8
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='SPLINE_TS'
# endif
# ifdef NEUMANN_TS
      is=ie+2 ; ie=is+9
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='NEUMANN_TS'
# endif
