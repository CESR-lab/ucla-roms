#ifndef LINEAR_DRAG_ONLY
        if (Zob > 0.) then
          do j=jstrV-1,jend
            do i=istrU-1,iend
# define nrdg nstp
              cff=sqrt( 0.333333333333*(
     &              u(i,j,1,nrdg)**2 +u(i+1,j,1,nrdg)**2
     &                     +u(i,j,1,nrdg)*u(i+1,j,1,nrdg)
     &               +v(i,j,1,nrdg)**2+v(i,j+1,1,nrdg)**2
     &                     +v(i,j,1,nrdg)*v(i,j+1,1,nrdg)
     &                                               ) )

# undef nrdg

! Canonical log-layer drag law formula valid for Zob << Hz only.   It
! should NEVER be used as it may cause instability in extremely shallow
! areas due to the singularity as Zob approaches Hz  (e.g., Yusuke PV4
! configuration in the extreme low-tide phase).

c**           rd(i,j)=rdrg + cff*(vonKar/log(Hz(i,j,1)/Zob))**2

! The following formula comes from the finite-volume interpretation of
! log-layer in the bottom-most grid box derived WITHOUT using Zob << Hz
! assumption. It covers the entire range of from Zob << Hz to Zob >> Hz
! and correctly asymptotes to the viscous limit in the latter extreme.
! The next formula after it is a simplified version which nevertheless
! yields correct asymptotic limits in both extremes, safely passing
! through Zob=Hz transition, and closely matching the upper formula
! there: setting Zob=Hz makes [vonKar/(2*ln 2-1)]^2 = [vonKar/0.386]^2
! upper vs. [vonKar/ln(3/2)]^2 = [vonKar/0.405]^2 lower.

c**           rd(i,j)=rdrg + cff*(  vonKar/(   (1.+Zob/Hz(i,j,1))
c**                                  *log(1.+Hz(i,j,1)/Zob) -1. ) )**2

!  Removed the addition of a linear 'background' drag
c             rd(i,j)=rdrg+cff*( vonKar/log(1.+0.5*Hz(i,j,1)/Zob) )**2

              rd(i,j)= cff*( vonKar/log(1.+0.5*Hz(i,j,1)/Zob) )**2
# ifdef WEC

!! JM look at this. Not working atm. Intended to unify rd between wec and nowec
! Combined wave-current bottom drag coefficient rd at rho-points
! Based on formulation in Soulsby, 1995 (See Uchiyama, 2008, eq. 8)
!  ubmag = magnitude of orbital bottom velocity
!!       Uorb1 = 0.5*fr(i,j)*wh(i,j)/max(sinh(min(kD(i,j),khmax)),eps)
!! we need both orbital velocity and frq of wave field

!          tau_c = rd(i,j)*cff
!          tau_w = 0.695*(ubmag**1.48)*((Zob*fr(i,j))**0.52)

!          rd(i,j) = rd(i,j)*(1.0 + 1.2*((tau_w/(tau_w+tau_c))**3.2))

# endif


# if !defined IMPLICIT_BOTTOM_DRAG
#  if !defined IMPLCT_NO_SLIP_BTTM_BC
                                                     ! must have
              rd(i,j)=min(rd(i,j), 0.8*Hz(i,j,1)/dt) ! restriction
                                                     ! for stability
#  endif
# endif
            enddo
          enddo
        else  !<-- Zob > 0.
#endif /* ! LINEAR_DRAG_ONLY */

          do j=jstrV-1,jend
            do i=istrU-1,iend
              rd(i,j)=rdrg

              rd(i,j)=min(rd(i,j), 0.8*Hz(i,j,1)/dt)

            enddo
          enddo

#ifndef LINEAR_DRAG_ONLY
        endif  !<-- Zob > 0
#endif

! Save "rd" into shared array "r_D" for subsequent use in barotropic
! mode. Note that "rd" is computed with an extra row of points on the
! western and southern sides because this is needed for interpolation
! to U- and V-points.  Array "r_D" saved here inherits these extra
! rows into its periodic and/or computational margins, so there is no
! need for exchange call.

        call ext_copy_prv2shr_2d_tile(istr,iend,jstr,jend, rd,r_D)

