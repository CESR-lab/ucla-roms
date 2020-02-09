!
! $Id: init_floats.F,v 2.2 2004/04/15 00:26:08 capet Exp $
!
#include "cppdefs.h"
       subroutine init_floats
#ifdef FLOATS
!
!================================================== John M. Klinck ===
!  Copyright (c) 2000 Rutgers/UCLA                                   !
!================================================ Hernan G. Arango ===
!                                                                    !
!  This routine reads in and process initial float locations from    !
!  input floats script.                                              !
!                                                                    !
!=====================================================================

      implicit none
# include "param.h"
# include "grid.h"
# include "ocean3d.h"
# include "ncvars_floats.h"
# include "scalars.h"
# include "floats.h"
# include "init_floats.h"

      integer Ncount, i, icard, iunit, j, k, lstr, mc, nc
      integer index(Mfloats)
      integer lenstr
      real xfloat, yfloat, zfloat
      real Ip0(Mfloats), Jp0(Mfloats), lat(Mfloats), lon(Mfloats)
      character*35 frmt
      parameter (iunit=50)

      k=0

!---------------------------------------------------------------------
!  Read in initial float locations.
!---------------------------------------------------------------------

      lstr=lenstr(fposnam)
      open(iunit, file=fposnam(1:lstr), form='formatted',
     &                             status='old', err=195)
c     write(*,10)
c 10  format(/,' FLOATS processing parameters:',/)
!
!  Read input parameters according to their input card number.
!
      icard=0
      do while (icard<99)
!
!  Read in floats identification title.
!
        if (icard==1) then
          read(iunit,'(a)',err=70) Ftitle
          lstr=lenstr(Ftitle)
          write(*,'(1x,3A)')   '(', Ftitle(1:lstr), ')'

!  Read in initial floats location.

        elseif (icard==2) then
          i=0
          nfloats=0
c         write(*,30)
c30       format(1x,'Fcoor',2x,'Ftype',2x,'Fcount',2x,'Ft0',
c    &        2x,'Fx0',2x,'Fy0',2x,'Fz0',2x,'Fgrd',2x,'Fdt',
c    &                           2x,'Fdx',2x,'Fdy',2x,'Fdz',/)
          do while (.true.)
            i=i+1
            read(iunit,*,err=40) Ft0(i), Fx0(i), Fy0(i), Fz0(i),
     &                    Fgrd(i),Fcoor(i), Ftype(i), Fcount(i),
     &                           Fdt(i), Fdx(i), Fdy(i), Fdz(i)
            nfloats=nfloats+Fcount(i)
c# if !defined SPHERICAL
c           if (Fcoor(i)==1) then
c             frmt='(i1,i2,i5,f10.4,2f8.0,f8.2,2x,i1,4f9.3)'
c           else
c            frmt='(i1,i2,i5,f10.4,3f8.2,2x,i1,4f9.3)'
c           endif
c# else
c             frmt='(i1,i2,i5,f10.4,3f8.2,2x,i1,4f9.3)'
c# endif
c           write(*,frmt) Fcoor(i), Ftype(i), Fcount(i),
c    &                         Ft0(i), Fx0(i), Fy0(i), Fz0(i),
c    &                   Fgrd(i),Fdt(i), Fdx(i), Fdy(i), Fdz(i)
          enddo
  40      if (Ft0(i)/=99.0) then
            write(*,50) icard, i, fposnam
  50        format(/,' INIT_FLOATS - error while reading input card: ',
     &             i2, ', floats location entry: ',i3,/,15x,
     &             'in input script: ',a)
            stop
          else
            Ncount=i-1
            i_floats=i-1
            goto 90
          endif
          if (Mfloats<nfloats) then
            write(*,60) Mfloats, nfloats
  60        format(/,' INIT_FLOATS - too small dimension parameter,',
     &             ' Mfloats',2i6,/,15x,'change file  param.h  and',
     &             ' recompile.')
            stop
          endif
        endif

!  Read last input card ID.

        read(iunit,*,err=70) icard
      enddo
      goto 90
!
!  Error while reading input parameters.
!
  70  write(*,80) icard, fposnam
  80  format(/,' INIT_FLOATS - error while reading input card: ',
     &       i2,15x,'in input script: ',a)

  90  close(iunit)
      write(*,100) nfloats
 100  format(/,2x,i6,4x,'nfloats',t26,
     &       'Number of float trajectories to compute.',/)

!---------------------------------------------------------------------
! Set initial float location.
!---------------------------------------------------------------------
!
! Set time of float release (seconds after model initialization) and
! initial float horizontal positions (grid units).  Fill the initial
! vertical level or depth position.

      mc=0
      nc=0
      do i=1,Ncount
        if (Fcount(i)==1) then
          nc=nc+1
          Tinfo(itstr,nc)=(tdays+Ft0(i))*day2sec
          Tinfo(izgrd,nc)=Fz0(i)
          Tinfo(igrd,i)=FLOAT(Fgrd(i)) ! each float belongs to this grid initially
          if (Fcoor(i)==0) then
             Tinfo(ixgrd,nc)=MIN(MAX(0.5,Fx0(i)),FLOAT(LLm)+0.5)
             Tinfo(iygrd,nc)=MIN(MAX(0.5,Fy0(i)),FLOAT(MMm)+0.5)
          else
            mc=mc+1
            lon(mc)=Fx0(i)
            lat(mc)=Fy0(i)
            index(mc)=nc
          endif
        elseif (Fcount(i)>1) then
          do j=1,Fcount(i)
            nc=nc+1
            Tinfo(igrd,nc)=FLOAT(Fgrd(i)) ! each float belongs to this grid initially
            if (Fdt(i)>0.0) then
              Tinfo(itstr,nc)=(tdays+Ft0(i)+FLOAT(j-1)*Fdt(i))*
     &                        day2sec
              Tinfo(izgrd,nc)=Fz0(i)
              if (Fcoor(i)==0) then
                Tinfo(ixgrd,nc)=MIN(MAX(0.5,Fx0(i)),FLOAT(LLm)+0.5)
                Tinfo(iygrd,nc)=MIN(MAX(0.5,Fy0(i)),FLOAT(MMm)+0.5)
              else
                mc=mc+1
                lon(mc)=Fx0(i)
                lat(mc)=Fy0(i)
                index(mc)=nc
              endif
            else
              Tinfo(itstr,nc)=(tdays+Ft0(i))*day2sec
              if (Fdz(i)==0.0) then
                Tinfo(izgrd,nc)=Fz0(i)
              else
                if (Fz0(i)>0.0) then
                  zfloat=Fz0(i)+FLOAT(j-1)*Fdz(i)
                  Tinfo(izgrd,nc)=MIN(MAX(0.0,zfloat),FLOAT(N))
                else
                  Tinfo(izgrd,nc)=Fz0(i)+FLOAT(j-1)*Fdz(i)
                endif
              endif
              if (Fcoor(i)==0) then
                  xfloat=Fx0(i)+FLOAT(j-1)*Fdx(i)
                  Tinfo(ixgrd,nc)=MIN(MAX(0.5,xfloat),
     &                                FLOAT(LLm)+0.5)
                  yfloat=Fy0(i)+FLOAT(j-1)*Fdy(i)
                  Tinfo(iygrd,nc)=MIN(MAX(0.5,yfloat),
     &                                FLOAT(MMm)+0.5)
              else
                mc=mc+1
                index(mc)=nc
                lon(mc)=Fx0(i)+FLOAT(j-1)*Fdx(i)
                lat(mc)=Fy0(i)+FLOAT(j-1)*Fdy(i)
              endif
            endif
          enddo
        endif
      enddo
!
!  Set number of floats trajectories to compute.
!
      nfloats=nc
!
!  If applicable, convert floats initial (lon,lat) positions to grid
!  units.
!
      if (mc>0) then
# ifdef SPHERICAL
          call hindices (Ip0,Jp0,lon,lat,mc,lonr,latr)
# else
          call hindices (Ip0,Jp0,lon,lat,mc,xr,yr)
# endif

        do i=1,mc
          nc=index(i)
          Tinfo(ixgrd,nc)=MIN(MAX(0.5,Ip0(i)),FLOAT(LLm)+0.5)
          Tinfo(iygrd,nc)=MIN(MAX(0.5,Jp0(i)),FLOAT(MMm)+0.5)
        enddo
      endif
!
      return

 195  write(*,205) fposnam
 205  format(/,'float file ',A,/,' not found => float initialization',
     &   ' using restart file')

          call get_init_flts
          ldefflt= .false.

C! Test de l'initialisation des flotteurs xa
C      do nc=1,nfloats
C      write(*,53) nc, Tinfo(ixgrd,nc), Tinfo(iygrd,nc),
C     &    Tinfo(izgrd,nc), NINT(Tinfo(igrd,nc)), Tinfo(itstr,nc)
C      enddo
C! Fin test
#endif /* FLOATS */
      return
      end
#ifdef FLOATS
      subroutine hindices (Ipos,Jpos,Xpos,Ypos,Npos,Xgrd,Ygrd)
!
!================================================ Hernan G. Arango ===
!  Copyright (c) 2000 Rutgers/UCLA                                   !
!======================================== Alexander F. Shchepetkin ===
!                                                                    !
!  Given position vectors Xpos and Ypos of size Npos,  this routine  !
!  finds the corresponding indices Ipos and Jpos of the  model grid  !
!  (Xgrd,Ygrd) cell containing each requested position.              !
!                                                                    !
!  Calls:    Try_Range                                               !
!                                                                    !
!=====================================================================
!
#ifdef AGRIF
        USE Agrif_Util
#endif
      implicit none
# include "param.h"
!
      logical found
      logical Try_Range
      integer Imax, Imin, Jmax, Jmin, Npos, i0, j0, k
      real Ipos(Npos), Jpos(Npos), Xpos(Npos), Ypos(Npos),
     &  deltax, deltay, dx1, dy1, dx2, dy2, c1, c2
      real Xgrd(GLOBAL_2D_ARRAY),
     &     Ygrd(GLOBAL_2D_ARRAY)
!
!---------------------------------------------------------------------
!  Determine grid cell indices containing requested position points.
!  Then, interpolate to fractional cell position.
!---------------------------------------------------------------------
!
!  Initialize all indices.
!
      do k=1,Npos
        Ipos(k)=0.0
        Jpos(k)=0.0
      enddo
!
!  Check each position to find if it falls inside the whole domain.
!  Once it is stablished that it inside, find the exact cell to which
!  it belongs by successively dividing the domain by a half (binary
!  search).
!
      do k=1,Npos

      found=Try_Range(0,LLm+1,0,MMm+1,Xpos(k),Ypos(k),Xgrd,Ygrd)

      if (found) then
          Imin=0
          Imax=LLm+1
          Jmin=0
          Jmax=MMm+1
          do while (((Imax-Imin)>1).or.((Jmax-Jmin)>1))
            if ((Imax-Imin)>1) then
              i0=(Imin+Imax)/2
              found=Try_Range(Imin,i0,Jmin,Jmax,Xpos(k),Ypos(k),
     &                        Xgrd,Ygrd)
              if (found) then
                Imax=i0
              else
                Imin=i0
              endif
            endif
            if ((Jmax-Jmin)>1) then
              j0=(Jmin+Jmax)/2
              found=Try_Range(Imin,Imax,Jmin,j0,Xpos(k),Ypos(k),
     &                        Xgrd,Ygrd)
              if (found) then
                Jmax=j0
              else
                Jmin=j0
              endif
            endif
          enddo

! Improved interpolation block (02/06/03). Not totally
! accurate in case where the grid is strongly deformed.
! The nonalignment of xsi (eta) axis with longitude
! (latitude) direction is taken into account.

           dy1=Ygrd(Imin,Jmin+1)-Ygrd(Imin,Jmin)
           dx1=Xgrd(Imin,Jmin+1)-Xgrd(Imin,Jmin)
           dy2=Ygrd(Imin+1,Jmin)-Ygrd(Imin,Jmin)
           dx2=Xgrd(Imin+1,Jmin)-Xgrd(Imin,Jmin)

           c1=Xpos(k)        *dy1-Ypos(k)        *dx1
           c2=Xgrd(Imin,Jmin)*dy2-Ygrd(Imin,Jmin)*dx2
           deltax=(c1*dx2-c2*dx1)/(dx2*dy1-dy2*dx1)
           deltax=(deltax-Xgrd(Imin,Jmin))/dx2
           Ipos(k)=FLOAT(Imin)+MIN(MAX(0.0,deltax),1.0)

           c1=Xgrd(Imin,Jmin)*dy1-Ygrd(Imin,Jmin)*dx1
           c2=Xpos(k)        *dy2-Ypos(k)        *dx2
           deltay=(c1*dy2-c2*dy1)/(dx2*dy1-dy2*dx1)
           deltay=(deltay-Ygrd(Imin,Jmin))/dy1
           Jpos(k)=FLOAT(Jmin)+MIN(MAX(0.0,deltay),1.0)

        endif

      enddo
      return
      end
      function Try_Range (Imin,Imax,Jmin,Jmax,Xo,Yo,Xgrd,Ygrd)
!
!================================================ Hernan G. Arango ===
!  Copyright (c) 2000 Rutgers/UCLA                                   !
!======================================== Alexander F. Shchepetkin ===
!                                                                    !
!  Given a grided domain with matrix coordinates Xgrd and Ygrd, this !
!  function finds if the point (Xo,Yo)  is inside the box defined by !
!  the requested corners (Imin,Jmin) and (Imax,Jmax). It will return !
!  logical switch  Try_Range=.true.  if (Xo,Yo) is inside, otherwise !
!  it will return false.                                             !
!                                                                    !
!  Calls:   inside                                                   !
!                                                                    !
!=====================================================================
!
      implicit none
#include "param.h"
!
      logical Try_Range, inside
      integer Imax, Imin, Jmax, Jmin, Nb, NX, i, j, shft
      parameter(NX=2*(LLm+2)+2*(MMm+2)+1)
      real Xb(NX), Yb(NX), Xo, Yo
      real Xgrd(GLOBAL_2D_ARRAY),
     &     Ygrd(GLOBAL_2D_ARRAY)
!
!---------------------------------------------------------------------
!  Define closed polygon.
!---------------------------------------------------------------------
!
!  Note that the last point (Xb(Nb),Yb(Nb)) does not repeat first
!  point (Xb(1),Yb(1)).  Instead, in function inside, it is implied
!  that the closing segment is (Xb(Nb),Yb(Nb))-->(Xb(1),Yb(1)). In
!  fact, function inside sets Xb(Nb+1)=Xb(1) and Yb(Nb+1)=Yb(1).
!
      Nb=2*(Jmax-Jmin+Imax-Imin)
      shft=1-Imin
      do i=Imin,Imax-1
        Xb(i+shft)=Xgrd(i,Jmin)
        Yb(i+shft)=Ygrd(i,Jmin)
      enddo
      shft=1-Jmin+Imax-Imin
      do j=Jmin,Jmax-1
        Xb(j+shft)=Xgrd(Imax,j)
        Yb(j+shft)=Ygrd(Imax,j)
      enddo
      shft=1+Jmax-Jmin+2*Imax-Imin
      do i=Imax,Imin+1,-1
        Xb(shft-i)=Xgrd(i,Jmax)
        Yb(shft-i)=Ygrd(i,Jmax)
      enddo
      shft=1+2*Jmax-Jmin+2*(Imax-Imin)
      do j=Jmax,Jmin+1,-1
        Xb(shft-j)=Xgrd(Imin,j)
        Yb(shft-j)=Ygrd(Imin,j)
      enddo
!
!---------------------------------------------------------------------
!  Check if point (Xo,Yo) is inside of the defined polygon.
!---------------------------------------------------------------------
!
      Try_Range=inside(Xo,Yo,Xb,Yb,Nb)
      return
      end
      function inside (Xo,Yo,Xb,Yb,Nb)
!
!================================================ Hernan G. Arango ===
!  Copyright (c) 2000 Rutgers/UCLA                                   !
!======================================== Alexander F. Shchepetkin ===
!                                                                    !
!  Given the vectors Xb and Yb of size Nb, defining the coordinates  !
!  of a closed polygon,  this function find if the point (Xo,Yo) is  !
!  inside the polygon.  If the point  (Xo,Yo)  falls exactly on the  !
!  boundary of the polygon, it still considered inside.              !
!                                                                    !
!  This algorithm does not rely on the setting of  Xb(Nb)=Xb(1) and  !
!  Yb(Nb)=Yb(1).  Instead, it assumes that the last closing segment  !
!  is (Xb(Nb),Yb(Nb)) --> (Xb(1),Yb(1)).                             !
!                                                                    !
!  Reference:                                                        !
!                                                                    !
!    Reid, C., 1969: A long way from Euclid. Oceanography EMR,       !
!      page 174.                                                     !
!                                                                    !
!  Algorithm:                                                        !
!                                                                    !
!  The decision whether the point is  inside or outside the polygon  !
!  is done by counting the number of crossings from the ray (Xo,Yo)  !
!  to (Xo,-infinity), hereafter called meridian, by the boundary of  !
!  the polygon.  In this counting procedure,  a crossing is counted  !
!  as +2 if the crossing happens from "left to right" or -2 if from  !
!  "right to left". If the counting adds up to zero, then the point  !
!  is outside.  Otherwise,  it is either inside or on the boundary.  !
!                                                                    !
!  This routine is a modified version of the Reid (1969) algorithm,  !
!  where all crossings were counted as positive and the decision is  !
!  made  based on  whether the  number of crossings is even or odd.  !
!  This new algorithm may produce different results  in cases where  !
!  Xo accidentally coinsides with one of the (Xb(k),k=1:Nb) points.  !
!  In this case, the crossing is counted here as +1 or -1 depending  !
!  of the sign of (Xb(k+1)-Xb(k)).  Crossings  are  not  counted if  !
!  Xo=Xb(k)=Xb(k+1).  Therefore, if Xo=Xb(k0) and Yo>Yb(k0), and if  !
!  Xb(k0-1) < Xb(k0) < Xb(k0+1),  the crossing is counted twice but  !
!  with weight +1 (for segments with k=k0-1 and k=k0). Similarly if  !
!  Xb(k0-1) > Xb(k0) > Xb(k0+1), the crossing is counted twice with  !
!  weight -1 each time.  If,  on the other hand,  the meridian only  !
!  touches the boundary, that is, for example, Xb(k0-1) < Xb(k0)=Xo  !
!  and Xb(k0+1) < Xb(k0)=Xo, then the crossing is counted as +1 for  !
!  segment k=k0-1 and -1 for segment k=k0, resulting in no crossing. !
!                                                                    !
!  Note 1: (Explanation of the logical condition)                    !
!                                                                    !
!  Suppose  that there exist two points  (x1,y1)=(Xb(k),Yb(k))  and  !
!  (x2,y2)=(Xb(k+1),Yb(k+1)),  such that,  either (x1 < Xo < x2) or  !
!  (x1 > Xo > x2).  Therefore, meridian x=Xo intersects the segment  !
!  (x1,y1) -> (x2,x2) and the ordinate of the point of intersection  !
!  is:                                                               !
!                                                                    !
!                 y1*(x2-Xo) + y2*(Xo-x1)                            !
!             y = -----------------------                            !
!                          x2-x1                                     !
!                                                                    !
!  The mathematical statement that point  (Xo,Yo)  either coinsides  !
!  with the point of intersection or lies to the north (Yo>=y) from  !
!  it is, therefore, equivalent to the statement:                    !
!                                                                    !
!         Yo*(x2-x1) >= y1*(x2-Xo) + y2*(Xo-x1),   if   x2-x1 > 0    !
!  or                                                                !
!         Yo*(x2-x1) <= y1*(x2-Xo) + y2*(Xo-x1),   if   x2-x1 < 0    !
!                                                                    !
!  which, after noting that  Yo*(x2-x1) = Yo*(x2-Xo + Xo-x1) may be  !
!  rewritten as:                                                     !
!                                                                    !
!        (Yo-y1)*(x2-Xo) + (Yo-y2)*(Xo-x1) >= 0,   if   x2-x1 > 0    !
!  or                                                                !
!        (Yo-y1)*(x2-Xo) + (Yo-y2)*(Xo-x1) <= 0,   if   x2-x1 < 0    !
!                                                                    !
!  and both versions can be merged into  essentially  the condition  !
!  that (Yo-y1)*(x2-Xo)+(Yo-y2)*(Xo-x1) has the same sign as x2-x1.  !
!  That is, the product of these two must be positive or zero.       !
!                                                                    !
!=====================================================================
!
      implicit none
!
      logical inside
      integer Nb, Nstep, crossings, i, inc, k, kk, nc
      parameter (Nstep=128)
      integer index(Nstep)
      real Xb(Nb+1), Yb(Nb+1), Xo, Yo, dx1, dx2, dxy
!
!---------------------------------------------------------------------
!  Find intersections.
!---------------------------------------------------------------------
!
!  Set crossings counter and close the contour of the polygon.
!
      crossings=0
      Xb(Nb+1)=Xb(1)
      Yb(Nb+1)=Yb(1)
!
!  The search is optimized.  First select the indices of segments
!  where Xb(k) is different from Xb(k+1) and Xo falls between them.
!  Then, further investigate these segments in a separate loop.
!  Doing it in two stages takes less time because the first loop is
!  pipelined.
!
      do kk=0,Nb-1,Nstep
        nc=0
        do k=kk+1,MIN(kk+Nstep,Nb)
          if (((Xb(k+1)-Xo)*(Xo-Xb(k))>=0.0).and.
     &        (Xb(k)/=Xb(k+1))) then
            nc=nc+1
            index(nc)=k
          endif
        enddo
        do i=1,nc
          k=index(i)
          if (Xb(k)/=Xb(k+1)) then
            dx1=Xo-Xb(k)
            dx2=Xb(k+1)-Xo
            dxy=dx2*(Yo-Yb(k))-dx1*(Yb(k+1)-Yo)
            inc=0
            if ((Xb(k)==Xo).and.(Yb(k)==Yo)) then
              crossings=1
              goto 10
            elseif (((dx1==0.0).and.(Yo>=Yb(k  ))).or.
     &              ((dx2==0.0).and.(Yo>=Yb(k+1)))) then
              inc=1
            elseif ((dx1*dx2>0.0).and.             ! See Note 1
     &              ((Xb(k+1)-Xb(k))*dxy>=0.0)) then
              inc=2
            endif
            if (Xb(k+1)>Xb(k)) then
              crossings=crossings+inc
            else
              crossings=crossings-inc
            endif
          endif
        enddo
      enddo
!
!  Determine if point (Xo,Yo) is inside of closed polygon.
!
  10  if (crossings==0) then
        inside=.false.
      else
        inside=.true.
      endif
      return
      end
#endif /* FLOATS */
