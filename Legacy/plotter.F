! The purpose of this code is to visualize ROMS standard output
! integral characteristics, such as time histories of kinetic
! energy, etc. This code is linked to NCAR graphics librart and
! therefore it must be a single precision code.
!

#include "cppdefs.h"

      implicit none
      integer*4, parameter :: iin=11, ndots_max=200000,
     &           solid_line=1, dashed_line=2, dotted_line=3
      real*4, parameter :: angle=0., left=-1., ctr=0., right=+1.,
     &           vL=0.15, vR=0.975, vB=0.15, vT=0.975,
     &        thick_line=2.2,   medium_line=1.1,   thin_line=0.6,
     &        large_font=0.025, medium_font=0.011, small_font=0.009

      character filename*128, label*6, string*64
      integer*4 last_arg, ne, is,ie, i,j,n,
     &        ndots, iic, icolor, iocheck,ipa, lab_int, ilab
      real*4 xx(ndots_max), tmin,tmax, t,dt, tdays,  xp,
     &       yy(ndots_max), emin,emax, y,de, PA(16), yp,
     &         level
#ifdef INT_IARGC
      integer iargc
#endif

      data icolor/0/  ipa/1/
      logical start_range, inside_range, log_scale

c--#define POSTSCRIPT
#ifdef POSTSCRIPT
      integer ngpswk
      call gopks(6,i)
      call ngsetc('ME','energy.ps')
      call gopwk(1,2,ngpswk('PS','PORTRAIT','COLOR')) !!'MONOCHROME'))
      call gacwk(1)
#else
      call opngks
#endif
      last_arg=iargc()
      if (last_arg==0) stop '### ERROR: NO FILE TO PROCESS'

#if defined SEAMOUNT
      tmin=0.
      tmax=250.

      ipa=1

      if (ipa==1 .or. ipa==2) then
        log_scale=.true.
        emin=3.E-6
        emax=1.1E-3
         emin=3.E-10

c        emin=3.E-10
c        emax=0.12E-3


c        emin=3.E-7    ! Comparison
c        emax=0.3E-3   ! among 30.XXXh


      elseif (ipa>=3 .and. ipa<=5) then
        log_scale=.false.
        emin=0.
        emax=0.35
        emax=0.175   ! 30.XXXh
      elseif (ipa==6) then
        log_scale=.false.
        emin=0.
        emax=1.05
      endif

#elif defined DAMEE_B
      tmin=0.
      tmax=tmin + 360.*20.
# ifdef REST_STATE_TEST
      emin=0.2E-5
      emax=2.0E-4
      log_scale=.true.
# else
      emin=0.5E-3
      emax=2.25E-3
      log_scale=.false.
# endif
#elif defined PACIFIC
      tmin=15.
      tmax=10800.
      ipa=1
      if (ipa==1) then
        emax=0.0024
        emin=0.0010
      elseif (ipa==2) then
        emax=0.0007
        emin=0.00035
      endif
#else
#error No application specified: Check main switch in "cppdefs.h".
#endif

      ne=7
      de=(emax-emin)/float(ne)

      write(*,*) vL,vR, vB,vT, tmin,tmax, emin,emax

      call set (vL,vR, vB,vT, tmin,tmax, emin,emax, 1)
      call gsclip (0)
      call gslwsc (thin_line)
      call line(tmin,emin, tmax,emin)
      call line(tmin,emin, tmin,emax)

      if (tmax-tmin>719.) then             ! Select time interval
        dt=360.                               ! [days] for horizontal
      elseif (tmax-tmin>89.) then          ! scale avoiding over
        dt=30.                                ! loading of the plot.
      elseif (tmax-tmin>29.) then
        dt=10.
      else
        dt=1.
      endif

      lab_int=1+int(medium_font*(tmax-tmin)/(0.1*dt))
      ilab=0
      t=tmin                                   ! Draw horizontal
      do while (t<=tmax)                     ! axis and label it.
        call gsln (solid_line)                 ! Automatically select
        if (mod(ilab, lab_int)==0) then      ! labeling interval
          if (tmax-tmin > 3600.) then       ! so the neighboring
            write(label,'(I6)') int(t-tmin)/360
          else
            write(label,'(I6)') int(t-tmin)
          endif
          is=1
          do while (label(is:is)==' ' .and. is<6)
            is=is+1
          enddo                                ! labels do not touch
          yp=emin-0.03*(emax-emin)             ! each other.
          call gslwsc (thick_line)
          call plchhq(t,yp, label(is:6), medium_font,angle,ctr)
          call gslwsc (thin_line)
          call gslwsc (medium_line)
        else
          call gslwsc (thin_line)
        endif
        yp=emin-0.01*(emax-emin)
        call line(t,emin, t,yp)
        call gsln (dotted_line)
        call line (t,emin, t,emax)
        t=t+dt
        ilab=ilab+1
      enddo

      if (ipa==1) then
        string=':F9:E:BF7:kin:N:,m:S:2:N:s:S:-2'
      elseif (ipa==2) then
        string=':F9:E:BF7:bar:N:,m:S:2:N:s:S:-2'
      elseif (ipa==3) then
        string=':F9:V:BF7:max:N:,ms:S:-1'
      elseif (ipa==4) then
        string=':F9:<V>:BF7:max:N:,ms:S:-1'
      elseif (ipa==5) then
        string=':F9:(V-<V>):BF7:max:N:,ms:S:-1'
      elseif (ipa==6) then
        string=':F9:E:BF7:bar:N:/:F9:E:BF7:kin:N:'

      elseif (ipa==3) then
        string=':F9:E:BF7:bcn:N:,m:S:2:N:s:S:-2'
      endif

      call gsln (solid_line)
      call gslwsc (thick_line)
      xp=tmin+0.025*(tmax-tmin)
      call plchhq (xp,emax, string, large_font,angle,left)
      if (log_scale) then
        call set (vL,vR, vB,vT, tmin,tmax, emin,emax, 2)
        do j=-10,+10
          level=10.D0**j
          if (emin<=level .and. level<=emax) then
            write(label,'(PE6.0)') level
            call gslwsc (thick_line)
            xp=tmin-0.025*(tmax-tmin)
            call plchhq (xp,level, label, medium_font,angle,right)
            call gslwsc (thin_line)
            xp=tmin-0.015*(tmax-tmin)
            call line(tmin,level, xp, level)
            call gsln (dotted_line)
            call line(tmin,level, tmax,level)
            call gsln (solid_line)
            do i=2,9
              y=0.1*level*float(i)
              if (emin<=y .and. y<=emax) then
                xp=tmin-0.005*(tmax-tmin)
                call line(tmin,y, xp,y)
              endif
            enddo
          endif
        enddo
      else
        do j=1,ne
          level=emin+j*de
          write(label,'(F6.4)') level
          call gslwsc (thick_line)
          xp=tmin-0.025*(tmax-tmin)
          call plchhq(xp,level, label, medium_font,angle,right)
          call gslwsc (thin_line)
          xp= tmin-0.015*(tmax-tmin)
          call line(tmin,level, xp, level)
        enddo
      endif
      call gsclip (1)
c*      goto 100

      do n=1,last_arg
        filename=' '
        call getarg(n,filename)
        write(*,'(1x,2A)') 'Processing ',  filename
        open(unit=iin,file=filename,form='formatted')
        ndots=0
        start_range=.false.
        inside_range=.true.
        do while (inside_range)
          if (start_range) then
            read(iin,*,iostat=iocheck,end=98) iic,
#ifdef SEAMOUNT
     &                        tdays, (PA(i), i=1,6)
#else
     &                        tdays, (PA(i), i=1,4)
#endif
            if (iocheck==0) then
              if (tdays>tmax) inside_range=.false.
              if (ndots < ndots_max) then
                ndots=ndots+1
                xx(ndots)=tdays
                yy(ndots)=PA(ipa)
              else
                write(*,'(/1x,2A/16x,A,I9,1x,A/)') '### WARNING: ',
     &                'Exceeding maximum number of points allowed',
     &                'ndots_max =', ndots_max, 'record truncated.'
                goto 99
              endif
            endif
          else
            read(iin,'(A)',iostat=iocheck) string
            if (string(2:5)=='STEP') start_range=.true.
          endif
        enddo
  98    close (iin)
        write(*,*) 'ndots=',ndots
  99    if (ndots>1) then
c          if (icolor==0 .or. xx(1)-tmin < 0.01*(tmax-tmin)) then
            call set_new_color (icolor, last_arg)
c          endif
          call gslwsc (medium_line)
          call curved (xx,yy,ndots)
c**       do i=1,ndots
c**         write(*,*) i, xx(i),yy(i)
c**       enddo
        endif

        is=1
        do while (filename(is:is)==' ')
          is=is+1
        enddo
        i=is+1
        do while (filename(i:i)/=' ')
          ie=i
          i=i+1
        enddo
        if (filename(is:is+4)=='trace') is=is+5
        if (filename(is:is)=='.'.or.filename(is:is)=='_') is=is+1
        xp=xx(ndots)+0.01*(tmax-tmin)
        call plchhq (xp,yy(ndots), filename(is:ie),
     &                       small_font,angle,left)
      enddo

 100  call frame
#ifdef POSTSCRIPT
      call gdawk(1)
      call gclwk(1)
      call gclks
#else
      call clsgks
      call system ('idt gmeta &')
#endif
      stop
      end

c--#define  POSTSCRIPT

      subroutine set_new_color (icolor, maxcolor)
      implicit none
      integer*4 icolor, maxcolor
      real*4 frac, R,G,B
#ifdef POSTSCRIPT
      integer*4 ifrac
      real*4 RGB(3,30)
      data RGB/0.15,0.00,0.37,  0.07,0.05,0.47,  0.00,0.10,0.56,
     &         0.07,0.20,0.64,  0.20,0.30,0.71,  0.25,0.40,0.77,
     &         0.25,0.50,0.80,  0.35,0.60,0.85,  0.50,0.70,0.90,
     &         0.50,0.80,0.75,  0.50,0.80,0.60,  0.50,0.80,0.55,
     &         0.50,0.80,0.30,  0.50,0.85,0.30,  0.50,0.90,0.30,
     &         0.65,0.90,0.15,  0.80,0.90,0.00,  0.80,0.85,0.00,
     &         0.80,0.80,0.00,  0.85,0.75,0.00,  0.90,0.70,0.00,
     &         0.93,0.60,0.00,  0.95,0.50,0.00,  0.95,0.40,0.00,
     &         0.93,0.30,0.00,  0.88,0.15,0.00,  0.82,0.00,0.04,
     &         0.73,0.00,0.13,  0.66,0.00,0.20,  0.60,0.00,0.30/
#else
      real*4 Hue,Sat,Val
      parameter (Sat=1., Val=1.)
#endif
      icolor=max(1,icolor)    ! Check, whether color index has been
      icolor=icolor+1         ! initialized, advance it, then define
      if (maxcolor>1) then
        frac=float(icolor-2)/float(maxcolor-1)-0.5
      else
        frac=0.               ! fractional variable "frac" within the
      endif                   ! range -1/2 < frac < 1/2
#ifdef POSTSCRIPT
      frac=15.5 + 28.999*frac
      ifrac=int(frac)
      frac=frac-float(ifrac)
      write(*,*) 'set_new_color, ifrac,frac =', ifrac,frac

      R=(1.-frac)*RGB(1,ifrac)+frac*RGB(1,ifrac+1)
      G=(1.-frac)*RGB(2,ifrac)+frac*RGB(2,ifrac+1)
      B=(1.-frac)*RGB(3,ifrac)+frac*RGB(3,ifrac+1)
#else
      Hue=135.-300.*frac            ! This algorithm sets bright,
      if (Hue<0.) then           ! the most distinctive colors,
        Hue=Hue+360.                ! centered around grass green
      elseif (Hue>360.) then     ! (Hue=135) using color-wheel
        Hue=Hue-360.                ! method.
      endif
      call hsvrgb (Hue,Sat,Val, R,G,B)
      write(*,*) 'set_new_color, Hue,Sat,Val =', Hue,Sat,Val
#endif
      write(*,*) 'set_new_color, icolor=', icolor, 'R,G,B =', R,G,B
      call gscr (1, icolor, R,G,B)
      call gsplci (icolor)
      end

