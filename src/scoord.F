      module scoord

#include "cppdefs.opt"

      use param
      implicit none

! Control parameters for vertical coordinate transform: "theta_s" and
! "theta_b" are surface and bottom refinement coefficients for Cs=Cs(s)
! curves; "hc" is critical depth [meters] above which the vertical grid
! spacing remains approximately uniform while becoming stretched below.
! Once Cs(s) is specified, the unperturbed (zeta=0) vertical coordinate
! transformation z=z(s) is defined as either,
!
!               z = hc*s + (h-hc)*Cs(s)                  SH94
! or
!                        hc*s + h*Cs(s)
!               z = h * ----------------                 SM09
!                            hc + h
!
! depending on setting VERT_COORD_TYPE_XXXX. The upper transform has
! the limitation of
!
!                       hc < hmin
!
! (otherwise z=z(s) looses monotonicity resulting in negative grid
! spacings) which severely restricts the choice of hc; in the second
! case "hc" can be selected independently from minimum depth "hmin".
!
! Setting theta_s=theta_b=0 in the case of SH94 yields vertically
! uniform spacing (plain sigma); For both cases increase of "theta_s"
! leads to a more concentrated resolution toward the surface; typically
! 0=<theta_s<10 for both types; Optimally chosen "hc" should be
! comparable to the expected thermocline depth.
! Arrays "Cs_w" and "Cs_r" are stretching curves, Cs=Cs(s), where "s"
! is sigma-coordinate, -1 < s < 0, while "_w" and "_r" stand for
! vertical W- and RHO-point locations.

#ifdef SOLVE3D
      real theta_s,theta_b, hc, Cs_w(0:N), Cs_r(N)
#endif

      end module scoord
