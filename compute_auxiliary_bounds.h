! Auxiliary module "compute_auxiliary_bounds.h":
!---------- ------ -----------------------------
! Compute derived bounds for the loop indices over a subdomain
! "tile". The extended bounds [labelled by suffix R] are designed to
! cover also the outer ghost points, if the subdomain "tile" is
! adjacent to a PHYSICAL boundary. (NOTE: istrR,iendR,jstrR,jendR
! computed by this module DO NOT COVER ghost points associated with
! periodic boundaries (if any) or with 2-point computational marhins
! of MPI subdomains.
! 
! This module also computes loop-bounds for U- and V-type variables
! which belong to the interior of the computational domain. These are
! labelled by suffixes U,V and they step one grid point inward from
! the side of the subdomain adjacent to the physical boundary.
! Conversely, for an internal subdomain [which does not have segments
! of physical boundary] all variables with suffixes R,U,V are set to
! the same values are the corresponding non-suffixed variables.
! 
! Because this module also contains type declarations for these
! bounds, it must be included just after the last type declaration
! inside a subroutine, but before the first executable statement.
 
#ifdef EW_PERIODIC
# undef istrU
# define istrU istr
# undef istrR
# define istrR istr
# undef iendR
# define iendR iend
#else
      integer istrU, istrR, iendR
#endif
 
#ifdef NS_PERIODIC
# undef jstrV
# define jstrV jstr
# undef jstrR
# define jstrR jstr
# undef jendR
# define jendR jend
#else
      integer jstrV, jstrR, jendR
#endif
 
#ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        istrR=istr-1
        istrU=istr+1
      else
        istrR=istr
        istrU=istr
      endif
      if (EASTERN_EDGE) then
        iendR=iend+1
      else
        iendR=iend
      endif
#endif
 
#ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        jstrR=jstr-1
        jstrV=jstr+1
      else
        jstrR=jstr
        jstrV=jstr
      endif
      if (NORTHERN_EDGE) then
        jendR=jend+1
      else
        jendR=jend
      endif
#endif
 
