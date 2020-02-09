!
! $Id: interp_rho.F,v 2.1 2004/04/09 21:23:33 capet Exp $
!
#include "cppdefs.h"
#ifdef FLOATS
      subroutine interp_r3d_type (A, ifld, nfltmax, indx)
      implicit none
# include "param.h"
# include "grid.h"
# include "floats.h"
# include "scalars.h"
      integer ifld, nfltmax, indx(nfltmax)
      real A(GLOBAL_2D_ARRAY,N)
      integer id,iflt, i1,i2, j1,j2, k1,k2
      real p1,p2, q1,q2, r1,r2, cff1, cff2

      do id=1,nfltmax
        iflt=indx(id)

        k1=int(track(izgrd,nfp1,iflt)+0.5)
        r2=track(izgrd,nfp1,iflt)+0.5 - float(k1)
        r1=1.-r2
        k1=max(k1,   1)
        k2=min(k1+1, N)

        i1=int(track(ixgrd,nfp1,iflt))
        i2=i1+1
        p2=track(ixgrd,nfp1,iflt)-float(i1)
        p1=1.-p2

        j1=int(track(iygrd,nfp1,iflt))
        j2=j1+1
        q2=track(iygrd,nfp1,iflt)-float(j1)
        q1=1.0-q2

#ifdef MASKING
        cff1=p1*q1*rmask(i1,j1)*(r1*A(i1,j1,k1)+r2*A(i1,j1,k2))
     &      +p2*q1*rmask(i2,j1)*(r1*A(i2,j1,k1)+r2*A(i2,j1,k2))
     &      +p1*q2*rmask(i1,j2)*(r1*A(i1,j2,k1)+r2*A(i1,j2,k2))
     &      +p2*q2*rmask(i2,j2)*(r1*A(i2,j2,k1)+r2*A(i2,j2,k2))

        cff2=q1*(p1*rmask(i1,j1) + p2*rmask(i2,j1))
     &      +q2*(p1*rmask(i1,j2) + p2*rmask(i2,j2))

        if (cff2>0.) then
          trackaux(ifld,iflt)=cff1/cff2
        else
          trackaux(ifld,iflt)=0.0
        endif
#else
        trackaux(ifld,iflt)=p1*q1*(r1*A(i1,j1,k1)+r2*A(i1,j1,k2))
     &                       +p2*q1*(r1*A(i2,j1,k1)+r2*A(i2,j1,k2))
     &                       +p1*q2*(r1*A(i1,j2,k1)+r2*A(i1,j2,k2))
     &                       +p2*q2*(r1*A(i2,j2,k1)+r2*A(i2,j2,k2))
#endif
      enddo
      return
      end


      subroutine interp_w3d_type (A, ifld, nfltmax,indx)
      implicit none
# include "param.h"
# include "floats.h"
# include "scalars.h"
      integer ifld, nfltmax, indx(nfltmax)
      real A(GLOBAL_2D_ARRAY,0:N)
      integer id,iflt, i1,i2, j1,j2, k1,k2
      real p1,p2, q1,q2, r1,r2

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

        trackaux(ifld,iflt)=p1*q1*(r1*A(i1,j1,k1)+r2*A(i1,j1,k2))
     &                     +p2*q1*(r1*A(i2,j1,k1)+r2*A(i2,j1,k2))
     &                     +p1*q2*(r1*A(i1,j2,k1)+r2*A(i1,j2,k2))
     &                     +p2*q2*(r1*A(i2,j2,k1)+r2*A(i2,j2,k2))
      enddo
      return
      end


      subroutine interp_r2d_type (A, ifld, nfltmax, indx)
      implicit none
# include "param.h"
# include "floats.h"
# include "scalars.h"
      integer ifld, nfltmax, indx(nfltmax)
      real A(GLOBAL_2D_ARRAY)
      integer id,iflt, i1,i2, j1,j2
      real p1,p2, q1,q2

      do id=1,nfltmax
        iflt=indx(id)

        i1=int(track(ixgrd,nfp1,iflt))
        i2=i1+1
        p2=track(ixgrd,nfp1,iflt)-float(i1)
        p1=1.-p2

        j1=int(track(iygrd,nfp1,iflt))
        j2=j1+1
        q2=track(iygrd,nfp1,iflt)-float(j1)
        q1=1.0-q2

        trackaux(ifld,iflt)=q1*(p1*A(i1,j1) + p2*A(i2,j1))
     &                     +q2*(p1*A(i1,j2) + p2*A(i2,j2))
      enddo
      return
      end

      subroutine fill_ini
      implicit none
# include "param.h"
# include "floats.h"
# include "scalars.h"
# include "ncvars_floats.h"
      integer iflt
      do iflt=1,nfloats
         trackaux(ifdpt,iflt)=Tinfo(izgrd,iflt) ! Either the depth or the

         if (wrtflt(indxfltGrd)) then
           trackaux(izgrd,iflt)=Tinfo(izgrd,iflt)! sigma level position is
                                                 ! unavailaible
           trackaux(ixgrd,iflt)=Tinfo(ixgrd,iflt)! for initial storage
           trackaux(iygrd,iflt)=Tinfo(iygrd,iflt)! for initial storage
         endif
      enddo
      return
      end

      subroutine interp_r2d_type_ini (A, ifld)
      implicit none
# include "param.h"
# include "floats.h"
# include "scalars.h"
      integer ifld
      real A(GLOBAL_2D_ARRAY)
      integer iflt, i1,i2, j1,j2
      real p1,p2, q1,q2

      do iflt=1,nfloats

        i1=int(Tinfo(ixgrd,iflt))
        i2=i1+1
        p2=Tinfo(ixgrd,iflt)-float(i1)
        p1=1.-p2

        j1=int(Tinfo(iygrd,iflt))
        j2=j1+1
        q2=Tinfo(iygrd,iflt)-float(j1)
        q1=1.0-q2

        trackaux(ifld,iflt)=q1*(p1*A(i1,j1) + p2*A(i2,j1))
     &                       +q2*(p1*A(i1,j2) + p2*A(i2,j2))
      enddo
      return
      end

#else
      subroutine interp_floats_empty
      end
#endif
