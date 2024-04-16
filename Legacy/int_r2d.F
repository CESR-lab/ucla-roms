#define GRID_LEVEL 1
#define MAX_GRID_LEVEL 2
#include "cppdefs.h"
#if GRID_LEVEL < MAX_GRID_LEVEL

      subroutine int_r2d_west (istr,iend,jstr,jend, q,buff)
      implicit none
# include "param.h"
# include "grid.h"
      integer istr,iend,jstr,jend, jmin,jmax, i,j,jc
      real q(GLOBAL_2D_ARRAY)
      real buff(PRIVATE_2D_SCRATCH_ARRAY)
      real qC,qxR,qxL, qyR,qyL, qxyRR,qxyLR,qxyRL,qxyLL,
     &     cff, OneSixth, OneFourth, OneThird, TwoThird
      parameter (OneSixth=1./6., OneFourth=0.25,
     &           OneThird=1./3., TwoThird=2./3.)

      jmin=jmin_child+(jstr-1)/3
      jmax=jmin_child+(jend-1)/3
      i=imin_child-1

      do j=jmin,jmax
        if (rmask(i,j)>0.) then
          qC=q(i,j)
        else
          cff=rmask(i+1,j)+rmask(i+1,j)+rmask(i,j+1)+rmask(i,j-1)
     &                     +0.7071*( rmask(i+1,j+1)+rmask(i-1,j+1)
     &                             +rmask(i+1,j-1)+rmask(i-1,j-1))
          if (cff>0.) then
            qC=(  rmask(i+1,j)*q(i+1,j) +rmask(i+1,j)*q(i+1,j)
     &           +rmask(i,j+1)*q(i,j+1) +rmask(i,j-1)*q(i,j-1)
     &     +0.7071*(
     &        rmask(i+1,j+1)*q(i+1,j+1) +rmask(i-1,j+1)*q(i-1,j+1)
     &       +rmask(i+1,j-1)*q(i+1,j-1) +rmask(i-1,j-1)*q(i-1,j-1)
     &                                                     ))/cff
          else
            qC=0.
          endif
        endif

        qxR=(q(i+1,j)-qC)*umask(i+1,j)
        qxL=(qC-q(i-1,j))*umask(i  ,j)
        qyR=(q(i,j+1)-qC)*vmask(i,j+1)
        qyL=(qC-q(i,j-1))*vmask(i,j  )
        qxyRR=(q(i+1,j+1)-q(i,j+1))*umask(i+1,j+1)-qxR
     &       +(q(i+1,j+1)-q(i+1,j))*vmask(i+1,j+1)-qyR
        qxyLR=(q(i,j+1)-q(i-1,j+1))*umask(i,j+1)-qxL
     &       +qyR-(q(i-1,j+1)-q(i-1,j))*vmask(i-1,j+1)
        qxyRL=qxR-(q(i+1,j-1)-q(i,j-1))*umask(i+1,j-1)
     &       +(q(i+1,j)-q(i+1,j-1))*vmask(i+1,j)-qyL
        qxyLL=qxL-(q(i,j-1)-q(i-1,j-1))*umask(i,j-1)
     &       +qyL-(q(i-1,j)-q(i-1,j-1))*vmask(i-1,j)

        jc=2+3*(j-jmin_child)
        buff(istr-2,jc  )=qC
        buff(istr-2,jc+1)=qC +OneSixth*(qyR+qyL+OneThird*(qyR-qyL))
        buff(istr-2,jc-1)=qC -OneSixth*(qyR+qyL-OneThird*(qyR-qyL))

        buff(istr-1,jc  )=qC +OneSixth*(qxR+qxL+OneThird*(qxR-qxL))
        buff(istr-1,jc+1)=qC +OneSixth*( qxR+qxL +qyR+qyL
     &                        +OneThird*( qxR-qxL +qyR-qyL
     &                        +OneFourth*( qxyRR+qxyLR+qxyRL+qxyLL
     &                          +TwoThird*( qxyRR     -qxyLL
     &                           +OneSixth*( qxyRR-qxyLR-qxyRL+qxyLL
     &                                                         )))))
        buff(istr-1,jc-1)=qC +OneSixth*( qxR+qxL -qyR-qyL
     &                        +OneThird*( qxR-qxL +qyR-qyL
     &                        -OneFourth*( qxyRR+qxyLR+qxyRL+qxyLL
     &                          +TwoThird*(      -qxyLR+qxyRL
     &                           -OneSixth*( qxyRR-qxyLR-qxyRL+qxyLL
     &                                                         )))))
      enddo
      return
      end



      subroutine int_r2d_east (istr,iend,jstr,jend, q,buff)
      implicit none
# include "param.h"
# include "grid.h"
      integer istr,iend,jstr,jend, jmin,jmax, i,j,jc
      real q(GLOBAL_2D_ARRAY)
      real buff(PRIVATE_2D_SCRATCH_ARRAY)
      real qC,qxR,qxL, qyR,qyL, qxyRR,qxyLR,qxyRL,qxyLL,
     &     cff, OneSixth, OneFourth, OneThird, TwoThird
      parameter (OneSixth=1./6., OneFourth=0.25,
     &           OneThird=1./3., TwoThird=2./3.)

      jmin=jmin_child+(jstr-1)/3
      jmax=jmin_child+(jend-1)/3
      i=imax_child

      do j=jmin,jmax
        if (rmask(i,j)>0.) then
          qC=q(i,j)
        else
          cff=rmask(i+1,j)+rmask(i+1,j)+rmask(i,j+1)+rmask(i,j-1)
     &                     +0.7071*( rmask(i+1,j+1)+rmask(i-1,j+1)
     &                             +rmask(i+1,j-1)+rmask(i-1,j-1))
          if (cff>0.) then
            qC=(  rmask(i+1,j)*q(i+1,j) +rmask(i+1,j)*q(i+1,j)
     &           +rmask(i,j+1)*q(i,j+1) +rmask(i,j-1)*q(i,j-1)
     &     +0.7071*(
     &        rmask(i+1,j+1)*q(i+1,j+1) +rmask(i-1,j+1)*q(i-1,j+1)
     &       +rmask(i+1,j-1)*q(i+1,j-1) +rmask(i-1,j-1)*q(i-1,j-1)
     &                                                     ))/cff
          else
            qC=0.
          endif
        endif

        qxR=(q(i+1,j)-qC)*umask(i+1,j)
        qxL=(qC-q(i-1,j))*umask(i  ,j)
        qyR=(q(i,j+1)-qC)*vmask(i,j+1)
        qyL=(qC-q(i,j-1))*vmask(i,j  )
        qxyRR=(q(i+1,j+1)-q(i,j+1))*umask(i+1,j+1)-qxR
     &       +(q(i+1,j+1)-q(i+1,j))*vmask(i+1,j+1)-qyR
        qxyLR=(q(i,j+1)-q(i-1,j+1))*umask(i,j+1)-qxL
     &       +qyR-(q(i-1,j+1)-q(i-1,j))*vmask(i-1,j+1)
        qxyRL=qxR-(q(i+1,j-1)-q(i,j-1))*umask(i+1,j-1)
     &       +(q(i+1,j)-q(i+1,j-1))*vmask(i+1,j)-qyL
        qxyLL=qxL-(q(i,j-1)-q(i-1,j-1))*umask(i,j-1)
     &       +qyL-(q(i-1,j)-q(i-1,j-1))*vmask(i-1,j)

        jc=2+3*(j-jmin_child)
        buff(iend+2,jc  )=qC
        buff(iend+2,jc+1)=qC +OneSixth*(qyR+qyL+OneThird*(qyR-qyL))
        buff(iend+2,jc-1)=qC -OneSixth*(qyR+qyL-OneThird*(qyR-qyL))

        buff(iend+1,jc  )=qC -OneSixth*(qxR+qxL-OneThird*(qxR-qxL))
        buff(iend+1,jc+1)=qC +OneSixth*( -qxR-qxL +qyR+qyL
     &                        +OneThird*( qxR-qxL  +qyR-qyL
     &                        -OneFourth*( qxyRR+qxyLR+qxyRL+qxyLL
     &                          +TwoThird*(      qxyLR-qxyRL
     &                           -OneSixth*( qxyRR-qxyLR-qxyRL+qxyLL
     &                                                         )))))
        buff(iend+1,jc-1)=qC -OneSixth*( qxR+qxL +qyR+qyL
     &                        -OneThird*( qxR-qxL +qyR-qyL
     &                        +OneFourth*( qxyRR+qxyLR+qxyRL+qxyLL
     &                          -TwoThird*( qxyRR-qxyLL
     &                           -OneSixth*( qxyRR-qxyLR-qxyRL+qxyLL
     &                                                         )))))
      enddo
      return
      end



      subroutine int_r2d_south (istr,iend,jstr,jend, q,buff)
      implicit none
# include "param.h"
# include "grid.h"
      integer istr,iend,jstr,jend, imin,imax, i,j, ic
      real q(GLOBAL_2D_ARRAY)
      real buff(PRIVATE_2D_SCRATCH_ARRAY)
      real qC,qxR,qxL, qyR,qyL, qxyRR,qxyLR,qxyRL,qxyLL,
     &     cff, OneSixth, OneFourth, OneThird, TwoThird
      parameter (OneSixth=1./6., OneFourth=0.25,
     &           OneThird=1./3., TwoThird=2./3.)

      imin=imin_child+(jstr-1)/3
      imax=imin_child+(jend-1)/3
      j=jmin_child-1

      do i=imin,imax
        if (rmask(i,j)>0.) then
          qC=q(i,j)
        else
          cff=rmask(i+1,j)+rmask(i+1,j)+rmask(i,j+1)+rmask(i,j-1)
     &                     +0.7071*( rmask(i+1,j+1)+rmask(i-1,j+1)
     &                             +rmask(i+1,j-1)+rmask(i-1,j-1))
          if (cff>0.) then
            qC=(  rmask(i+1,j)*q(i+1,j) +rmask(i+1,j)*q(i+1,j)
     &           +rmask(i,j+1)*q(i,j+1) +rmask(i,j-1)*q(i,j-1)
     &     +0.7071*(
     &        rmask(i+1,j+1)*q(i+1,j+1) +rmask(i-1,j+1)*q(i-1,j+1)
     &       +rmask(i+1,j-1)*q(i+1,j-1) +rmask(i-1,j-1)*q(i-1,j-1)
     &                                                     ))/cff
          else
            qC=0.
          endif
        endif

        qxR=(q(i+1,j)-qC)*umask(i+1,j)
        qxL=(qC-q(i-1,j))*umask(i  ,j)
        qyR=(q(i,j+1)-qC)*vmask(i,j+1)
        qyL=(qC-q(i,j-1))*vmask(i,j  )
        qxyRR=(q(i+1,j+1)-q(i,j+1))*umask(i+1,j+1)-qxR
     &       +(q(i+1,j+1)-q(i+1,j))*vmask(i+1,j+1)-qyR
        qxyLR=(q(i,j+1)-q(i-1,j+1))*umask(i,j+1)-qxL
     &       +qyR-(q(i-1,j+1)-q(i-1,j))*vmask(i-1,j+1)
        qxyRL=qxR-(q(i+1,j-1)-q(i,j-1))*umask(i+1,j-1)
     &       +(q(i+1,j)-q(i+1,j-1))*vmask(i+1,j)-qyL
        qxyLL=qxL-(q(i,j-1)-q(i-1,j-1))*umask(i,j-1)
     &       +qyL-(q(i-1,j)-q(i-1,j-1))*vmask(i-1,j)

        ic=2+3*(i-imin_child)
        buff(ic  ,jstr-2)=qC
        buff(ic+1,jstr-2)=qC +OneSixth*(qxR+qxL+OneThird*(qxR-qxL))
        buff(ic-1,jstr-2)=qC -OneSixth*(qxR+qxL-OneThird*(qxR-qxL))

        buff(ic  ,jstr-1)=qC +OneSixth*(qyR+qyL+OneThird*(qyR-qyL))
        buff(ic+1,jstr-1)=qC +OneSixth*( qxR+qxL +qyR+qyL
     &                        +OneThird*( qxR-qxL +qyR-qyL
     &                        +OneFourth*( qxyRR+qxyLR+qxyRL+qxyLL
     &                          +TwoThird*( qxyRR     -qxyLL
     &                           +OneSixth*( qxyRR-qxyLR-qxyRL+qxyLL
     &                                                         )))))
        buff(ic-1,jstr-1)=qC +OneSixth*( -qxR-qxL +qyR+qyL
     &                        +OneThird*( qxR-qxL  +qyR-qyL
     &                        -OneFourth*( qxyRR+qxyLR+qxyRL+qxyLL
     &                          +TwoThird*(      qxyLR-qxyRL
     &                           -OneSixth*( qxyRR-qxyLR-qxyRL+qxyLL
     &                                                         )))))
      enddo
      return
      end



      subroutine int_r2d_north (istr,iend,jstr,jend, q,buff)
      implicit none
# include "param.h"
# include "grid.h"
      integer istr,iend,jstr,jend, imin,imax, i,j, ic
      real q(GLOBAL_2D_ARRAY)
      real buff(PRIVATE_2D_SCRATCH_ARRAY)
      real qC,qxR,qxL, qyR,qyL, qxyRR,qxyLR,qxyRL,qxyLL,
     &     cff, OneSixth, OneFourth, OneThird, TwoThird
      parameter (OneSixth=1./6., OneFourth=0.25,
     &           OneThird=1./3., TwoThird=2./3.)

      imin=imin_child+(jstr-1)/3
      imax=imin_child+(jend-1)/3
      j=jmax_child

      do i=imin,imax
        if (rmask(i,j)>0.) then
          qC=q(i,j)
        else
          cff=rmask(i+1,j)+rmask(i+1,j)+rmask(i,j+1)+rmask(i,j-1)
     &                     +0.7071*( rmask(i+1,j+1)+rmask(i-1,j+1)
     &                             +rmask(i+1,j-1)+rmask(i-1,j-1))
          if (cff>0.) then
            qC=(  rmask(i+1,j)*q(i+1,j) +rmask(i+1,j)*q(i+1,j)
     &           +rmask(i,j+1)*q(i,j+1) +rmask(i,j-1)*q(i,j-1)
     &     +0.7071*(
     &        rmask(i+1,j+1)*q(i+1,j+1) +rmask(i-1,j+1)*q(i-1,j+1)
     &       +rmask(i+1,j-1)*q(i+1,j-1) +rmask(i-1,j-1)*q(i-1,j-1)
     &                                                     ))/cff
          else
            qC=0.
          endif
        endif

        qxR=(q(i+1,j)-qC)*umask(i+1,j)
        qxL=(qC-q(i-1,j))*umask(i  ,j)
        qyR=(q(i,j+1)-qC)*vmask(i,j+1)
        qyL=(qC-q(i,j-1))*vmask(i,j  )
        qxyRR=(q(i+1,j+1)-q(i,j+1))*umask(i+1,j+1)-qxR
     &       +(q(i+1,j+1)-q(i+1,j))*vmask(i+1,j+1)-qyR
        qxyLR=(q(i,j+1)-q(i-1,j+1))*umask(i,j+1)-qxL
     &       +qyR-(q(i-1,j+1)-q(i-1,j))*vmask(i-1,j+1)
        qxyRL=qxR-(q(i+1,j-1)-q(i,j-1))*umask(i+1,j-1)
     &       +(q(i+1,j)-q(i+1,j-1))*vmask(i+1,j)-qyL
        qxyLL=qxL-(q(i,j-1)-q(i-1,j-1))*umask(i,j-1)
     &       +qyL-(q(i-1,j)-q(i-1,j-1))*vmask(i-1,j)

        ic=2+3*(i-imin_child)
        buff(ic  ,jend+2)=qC
        buff(ic+1,jend+2)=qC +OneSixth*(qxR+qxL+OneThird*(qxR-qxL))
        buff(ic-1,jend+2)=qC -OneSixth*(qxR+qxL-OneThird*(qxR-qxL))

        buff(ic  ,jend+1)=qC -OneSixth*(qyR+qyL-OneThird*(qyR-qyL))
        buff(ic+1,jend+1)=qC +OneSixth*( qxR+qxL -qyR-qyL
     &                        +OneThird*( qxR-qxL +qyR-qyL
     &                        -OneFourth*( qxyRR+qxyLR+qxyRL+qxyLL
     &                          +TwoThird*(      -qxyLR+qxyRL
     &                           -OneSixth*( qxyRR-qxyLR-qxyRL+qxyLL
     &                                                         )))))
        buff(ic-1,jend+1)=qC -OneSixth*( qxR+qxL +qyR+qyL
     &                        -OneThird*( qxR-qxL +qyR-qyL
     &                        +OneFourth*( qxyRR+qxyLR+qxyRL+qxyLL
     &                          -TwoThird*( qxyRR-qxyLL
     &                           -OneSixth*( qxyRR-qxyLR-qxyRL+qxyLL
     &                                                         )))))
      enddo
      return
      end
#endif




