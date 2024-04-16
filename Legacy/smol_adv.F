#include "cppdefs.h"
#ifdef SMOLARKIEWICZ
      subroutine smol_adv (itrc,trc,advtrc)
!
!--------------------------------------------------------------------
!  This subroutine evaluates horizontal and vertical advection terms
!  for tracers using the  Smolarkiewicz (1983) advection scheme.  It
!  uses  an  upstream  advection  scheme  with a  second  corrective
!  upstream step to reduce the implicit diffusion. An anti-diffusion
!  velocity is  computed  and  used in the  second  pass through the
!  advection operator.
!
!  On Input:
!
!     itrc      Tracer type array index.
!     trc       Tracer (scalar) field to advect.
!
!  On Output:
!
!     advtrc    Advected tracer field.
!
!  Calls:   smol_adiff, smol_ups
!
!
!  References:
!
!    Smolarkiewicz, P.K., 1983: A simple positive definite advection
!     scheme with small implicit diffusion. Monthly Weather Review,
!     111, 479-486.
!
!    Smolarkiewicz, P.K., 1984:  A fully  multidimensional positive
!     definite  advection  transport algorithm  with small implicit
!     diffusion. Journal of computational Physics, 54, 325-362.
!
!    Smolarkiewicz, P.K. and W.W. Grabowski, 1990:  The multidimen-
!     sional positive definite advection transport algorithm:  Non-
!     oscillatory  option.  Journal  of  Computational Physics, 86,
!     355-375.
!--------------------------------------------------------------------
!
      implicit none
# include "param.h"
# include "pconst.h"
# include "grid.h"
# include "ocean.h"
# include "scalars.h"
# include "work.h"

      integer i, iter, itrc, j, k
      real advtrc(0:L,0:M,N), trc(0:L,0:M,N),
     &         ua(L,0:M,N), va(0:L,M,N), wa(0:L,0:M,0:N)
      equivalence (ua,b3d), (va,c3d), (wa,d3d)
!
!  Load tracer to advect into array "advtrc" which share the same
!  memory space than utility work array "a3d".  See calling routine
!  "step3d".
!
      do k=1,N
        do j=0,M
          do i=0,L
            advtrc(i,j,k)=trc(i,j,k)
          enddo
        enddo
      enddo
c
c  Load advective horizontal velocities into "ua" and "va" arrays.
c
      do k=1,N
        do j=0,M
          do i=1,L
            ua(i,j,k)=u(i,j,k)
          enddo
        enddo
        do j=1,M
          do i=0,L
            va(i,j,k)=v(i,j,k)
          enddo
        enddo
      enddo
c
c  Load advective vertical velocity into "wa". Convert to true velocity.
c  Recall that w=omega*Hz/(m*n) and omega has units of 1/second.
c
      do k=1,Nm
        do j=0,M
          do i=0,L
            wa(i,j,k)=w(i,j,k)*pm(i,j)*pn(i,j)
          enddo
        enddo
      enddo
      do j=0,M
        do i=0,L
          wa(i,j,0)=0.
          wa(i,j,N)=0.
        enddo
      enddo
!
!  Advection scheme: The number of iterations ADV_ORD is optional,
!  the minimum is two.  Each additional iteration increases the
!  solution accuracy. In the first iteration an "upstream" scheme
!  is used in its classical sense, while each following corrective
!  iteration re-applies the scheme but with a specially defined
!  "anti-diffusion" velocity.
!
      do iter=1,adv_ord
!
!  Advect using the upstream advection operator. The array "advtrc"
!  is overwritten with the new advected tracer values.
!
        call smol_ups (itrc,advtrc,ua,va,wa)
!
!  Calculate "anti-diffusion velocities".  Previous advected
!  velocities (ua,va,wa) are overwritten with the new anti-diffusion
!  velocities.
!
        if (iter<adv_ord) then
          call smol_adiff (trc,advtrc,ua,va,wa)
        endif
      enddo
!
!  Compute tracer advection contribution.
!
      do k=1,N
        do j=0,M
          do i=0,L
            advtrc(i,j,k)=advtrc(i,j,k)-trc(i,j,k)
          enddo
        enddo
      enddo
#else
      subroutine smol_adv
#endif /* SMOLARKIEWICZ */
      return
      end

