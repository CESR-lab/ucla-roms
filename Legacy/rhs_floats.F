!
! $Id: rhs_floats.F,v 2.0 2004/03/29 20:11:03 capet Exp $
!
#include "cppdefs.h"
#ifdef FLOATS
#ifdef SOLVE3D


      subroutine rhs_floats (uf,vf,Wf, nfltmax,indx)
      implicit none
#  include "param.h"
#  include "grid.h"
#  include "floats.h"
#  include "scalars.h"
#  include "ocean3d.h"
      integer nfltmax,indx(nfltmax)
      real uf(GLOBAL_2D_ARRAY,  N)
CSDISTRIBUTE_RESHAPE A(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real vf(GLOBAL_2D_ARRAY,  N)
CSDISTRIBUTE_RESHAPE A(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real Wf(GLOBAL_2D_ARRAY,0:N)
CSDISTRIBUTE_RESHAPE A(BLOCK_PATTERN,*) BLOCK_CLAUSE
      integer id,iflt, i1,i2, j1,j2, k1,k2, kh1,kh2
      real p1,p2, q1,q2, r1,r2, cff1,cff2

      do id=1,nfltmax
        iflt=indx(id)

        k1=int(track(izgrd,nfp1,iflt))
        r2=track(izgrd,nfp1,iflt) - float(k1)
        r1=1.-r2
        k1=max(k1,   0)
        k2=min(k1+1, N)

        i1=int(track(ixgrd,nfp1,iflt))
        i2=i1+1
        p2=track(ixgrd,nfp1,iflt)-float(i1)
        p1=1.-p2

        j1=int(track(iygrd,nfp1,iflt))
        j2=j1+1
        q2=track(iygrd,nfp1,iflt)-float(j1)
        q1=1.0-q2

        kh1=max(k1,1)
        kh2=min(k2+1,N)
#  ifdef MASKING
        cff1=2.*q1*( p1*pm(i1,j1)*pn(i1,j1)*rmask(i1,j1)*(
     &                 r1*Wf(i1,j1,k1)/(Hz(i1,j1,kh1 )+Hz(i1,j1,k1 ))
     &                +r2*Wf(i1,j1,k2)/(Hz(i1,j1,k2-1)+Hz(i1,j1,kh2))
     &                                                              )

     &             +p2*pm(i2,j1)*pn(i2,j1)*rmask(i2,j1)*(
     &                 r1*Wf(i2,j1,k1)/(Hz(i2,j1,kh1 )+Hz(i2,j1,k1 ))
     &                +r2*Wf(i2,j1,k2)/(Hz(i2,j1,k2-1)+Hz(i2,j1,kh2))
     &                                                             ))

     &     +2.*q2*( p1*pm(i1,j2)*pn(i1,j2)*rmask(i1,j2)*(
     &                 r1*Wf(i1,j2,k1)/(Hz(i1,j2,kh1 )+Hz(i1,j2,k1 ))
     &                +r2*Wf(i1,j2,k2)/(Hz(i1,j2,k2-1)+Hz(i1,j2,kh2))
     &                                                              )

     &             +p2*pm(i2,j2)*pn(i2,j2)*rmask(i2,j2)*(
     &                 r1*Wf(i2,j2,k1)/(Hz(i2,j2,kh1 )+Hz(i2,j2,k1 ))
     &                +r2*Wf(i2,j2,k2)/(Hz(i2,j2,k2-1)+Hz(i2,j2,kh2))
     &                                                             ))

        cff2=q1*(p1*rmask(i1,j1) + p2*rmask(i2,j1))
     &      +q2*(p1*rmask(i1,j2) + p2*rmask(i2,j2))

        if (cff2>0.) then
          track(izrhs,nfp1,iflt)=cff1/cff2
        else
          track(izrhs,nfp1,iflt)=0.0
        endif
#  else
        track(izrhs,nfp1,iflt)=2.*q1*( p1*pm(i1,j1)*pn(i1,j1)*(
     &                 r1*Wf(i1,j1,k1)/(Hz(i1,j1,kh1 )+Hz(i1,j1,k1 ))
     &                +r2*Wf(i1,j1,k2)/(Hz(i1,j1,k2-1)+Hz(i1,j1,kh2))
     &                                                              )

     &             +p2*pm(i2,j1)*pn(i2,j1)*(
     &                 r1*Wf(i2,j1,k1)/(Hz(i2,j1,kh1 )+Hz(i2,j1,k1 ))
     &                +r2*Wf(i2,j1,k2)/(Hz(i2,j1,k2-1)+Hz(i2,j1,kh2))
     &                                                             ))

     &     +2.*q2*( p1*pm(i1,j2)*pn(i1,j2)*(
     &                 r1*Wf(i1,j2,k1)/(Hz(i1,j2,kh1 )+Hz(i1,j2,k1 ))
     &                +r2*Wf(i1,j2,k2)/(Hz(i1,j2,k2-1)+Hz(i1,j2,kh2))
     &                                                              )

     &             +p2*pm(i2,j2)*pn(i2,j2)*(
     &                 r1*Wf(i2,j2,k1)/(Hz(i2,j2,kh1 )+Hz(i2,j2,k1 ))
     &                +r2*Wf(i2,j2,k2)/(Hz(i2,j2,k2-1)+Hz(i2,j2,kh2))
     &                                                             ))
#  endif

        k1=int(track(izgrd,nfp1,iflt)+0.5)
        r2=track(izgrd,nfp1,iflt)+0.5 - float(k1)
        r1=1.-r2
        k1=max(k1,   1)
        k2=min(k1+1, N)

        i1=int(track(ixgrd,nfp1,iflt)+0.5)
        i2=i1+1
        p2=track(ixgrd,nfp1,iflt)+0.5 - float(i1)
        p1=1.-p2

        j1=int(track(iygrd,nfp1,iflt))
        j2=j1+1
#  ifdef MASKING
!
! Note that the following code segment yields
!
!           rmask(i1,j1) rmask(i1,j2)       q1   q2
!           ----------------------------------------------
!               1            1        0< q1,q2 <1, q1+q2=1
!               1            0               1    0
!               0            1               0    1
!               0            0               0    0
!
        q1=rmask(i1,j1)*( 1. - rmask(i1,j2)*(
     &     track(iygrd,nfp1,iflt)-float(j1) ))
        q2=rmask(i1,j2)*(1.-q1)
#  else
        q2=track(iygrd,nfp1,iflt)-float(j1)
        q1=1.-q2
#  endif
        track(ixrhs,nfp1,iflt)=0.5*q1*(
     &                   p1*(r1*uf(i1,j1,k1)+r2*uf(i1,j1,k2))
     &                               *(pm(i1-1,j1)+pm(i1,j1))

     &                  +p2*(r1*uf(i2,j1,k1)+r2*uf(i2,j1,k2))
     &                               *(pm(i2-1,j1)+pm(i2,j1))

     &      ) + 0.5*q2*( p1*(r2*uf(i1,j2,k1)+r2*uf(i1,j2,k2))
     &                               *(pm(i1-1,j2)+pm(i1,j2))

     &                  +p2*(r1*uf(i2,j2,k1)+r2*uf(i2,j2,k2))
     &                               *(pm(i2-1,j2)+pm(i2,j2))
     &                                                      )


        i1=int(track(ixgrd,nfp1,iflt))
        i2=i1+1
#  ifdef MASKING
        p1=rmask(i1,j1)*( 1. - rmask(i2,j1)*(
     &     track(ixgrd,nfp1,iflt)-float(i1) ))
        p2=rmask(i2,j1)*(1.-p1)

#  else
        p2=track(ixgrd,nfp1,iflt)-float(i1)
        p1=1.-p2
#  endif
        j1=int(track(iygrd,nfp1,iflt)+0.5)
        j2=j1+1
        q2=track(iygrd,nfp1,iflt)+0.5 - float(j1)
        q1=1.-q2

        track(iyrhs,nfp1,iflt)=0.5*p1*(
     &                   q1*(r1*vf(i1,j1,k1)+r2*vf(i1,j1,k2))
     &                               *(pn(i1,j1-1)+pn(i1,j1))

     &                  +q2*(r1*vf(i1,j2,k1)+r2*vf(i1,j2,k2))
     &                               *(pn(i1,j2-1)+pn(i1,j2))

     &      ) + 0.5*p2*( q1*(r1*vf(i2,j1,k1)+r2*vf(i2,j1,k2))
     &                               *(pn(i2,j1-1)+pn(i2,j1))

     &                  +q2*(r1*vf(i2,j2,k1)+r2*vf(i2,j2,k2))
     &                               *(pn(i2,j2-1)+pn(i2,j2))
     &                                                      )
      enddo
      return
      end

# else



      subroutine rhs_floats (uf,vf, nfltmax,indx)
      implicit none
#  include "param.h"
#  include "grid.h"
#  include "floats.h"
#  include "scalars.h"
      integer nfltmax,indx(nfltmax)
      real uf(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE A(BLOCK_PATTERN) BLOCK_CLAUSE
      real vf(GLOBAL_2D_ARRAY)
CSDISTRIBUTE_RESHAPE A(BLOCK_PATTERN) BLOCK_CLAUSE
      integer id,iflt, i1,i2, j1,j2
      real p1,p2, q1,q2

      do id=1,nfltmax
        iflt=indx(id)

        i1=int(track(ixgrd,nfp1,iflt)+0.5)
        i2=i1+1
        p2=track(ixgrd,nfp1,iflt)+0.5 - float(i1)
        p1=1.-p2

        j1=int(track(iygrd,nfp1,iflt))
        j2=j1+1
#  ifdef MASKING
        q1=rmask(i1,j1)*( 1. - rmask(i1,j2)*(
     &     track(iygrd,nfp1,iflt)-float(j1) ))
        q2=rmask(i1,j2)*(1.-q1)
#  else
        q2=track(iygrd,nfp1,iflt)-float(j1)
        q1=1.-q2
#  endif
        track(ixrhs,nfp1,iflt)=0.5*q1*(
     &                   p1*uf(i1,j1)*(pm(i1-1,j1)+pm(i1,j1))
     &                  +p2*uf(i2,j1)*(pm(i2-1,j1)+pm(i2,j1))
     &      ) + 0.5*q2*( p1*uf(i1,j2)*(pm(i1-1,j2)+pm(i1,j2))
     &                  +p2*uf(i2,j2)*(pm(i2-1,j2)+pm(i2,j2))
     &                                                      )


        i1=int(track(ixgrd,nfp1,iflt))
        i2=i1+1
#  ifdef MASKING
        p1=rmask(i1,j1)*( 1. - rmask(i2,j1)*(
     &     track(ixgrd,nfp1,iflt)-float(i1) ))
        p2=rmask(i2,j1)*(1.-p1)

#  else
        p2=track(ixgrd,nfp1,iflt)-float(i1)
        p1=1.-p2
#  endif
        j1=int(track(iygrd,nfp1,iflt)+0.5)
        j2=j1+1
        q2=track(iygrd,nfp1,iflt)+0.5 - float(j1)
        q1=1.-q2

        track(iyrhs,nfp1,iflt)=0.5*p1*(
     &                   q1*vf(i1,j1)*(pn(i1,j1-1)+pn(i1,j1))
     &                  +q2*vf(i1,j2)*(pn(i1,j2-1)+pn(i1,j2))
     &      ) + 0.5*p2*( q1*vf(i2,j1)*(pn(i2,j1-1)+pn(i2,j1))
     &                  +q2*vf(i2,j2)*(pn(i2,j2-1)+pn(i2,j2))
     &                                                      )
      enddo
      return
      end
# endif  /* SOLVE3D */
#else
      subroutine rhs_floats_empty
      end
#endif    /* FLOATS */
