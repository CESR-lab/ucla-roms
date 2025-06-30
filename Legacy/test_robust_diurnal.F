      program test_robust_diurnal

! A program for testing ROBUST_DIURNAL_SRFLUX diurnal cycle modulation
! algorithm in "set_srflux.F". This program not a part of ROMS code per
! se.  The idealized diurnal cycle consists of 12-hour, cos^2-shaped
! pulses with 12-hour pauses in between,
!
!          tdays <= 1/4   Ampl = 0
!   1/4  < tdays  < 3/4   Ampl = 2 + 2*cos[2*pi*(2*tdays-1)]
!   3/4 <= tdays          Ampl = 0
!
! within each day. The above Amp=Amp(tdays) has the property that
! its integral over one day is exactly equal to 1.  The purpose of
! the algorithm below is not to sample Amp(tdays) at the discrete
! moments corresponding to each discrete moment t=t(it)+0.5*dt, but
! integrate it analytically within each time step [for each interval
! from time(it) to time(it)+dt] so as the result the discrete sum
! of Ampl(it) over each pulse EXACTLY adds up to 1 without relying
! on the smallness of time step relative to 12-hour pulse duration.
! The algorithm below utilized
!
! int_t1^t2 Amp(t) dt = 2*(t2-t2) + 1/(2*pi)*{ sin[2*pi*(2*t2-1)]
!                                             -sin[2*pi*(2*t1-1)] }
!
! where a special care needs to be taken to restrict t1,t2 to be
! within the range where Ampl >= 0, and furthermore, algebraically
! transform it to avoid to avoid calculation of the difference of
! two sin-function of close arguments in order to avoid roundoff
! errors, sin(x2)-sin(x1) = 2*cos((x1+x2)/2)*sin((x2-x2)/2).

      implicit none
      real*8, parameter :: pi=3.14159265358979323D0,
     &                         sec2day=1.D0/86400.D0
      real*8 dt,tdays, sum, dtau, cff, tau1, tau2, Ampl
      integer it

      dt = 5432.1D0

      sum=0.D0
      do it=1,60
        tdays=dt*sec2day * dble(it-1)

#define COS_CYCLE
#ifdef COS_CYCLE

/*
        cff=tdays-dble(int(tdays)) +0.5D0*dt*sec2day
        if (0.25D0 < cff .and. cff < 0.75D0) then
          Ampl=pi*cos(pi*(2.D0*cff-1.D0))
        else
          Ampl=0.D0
        endif
*/

      dtau=2.D0*pi * dt*sec2day
      cff=4.D0*pi*(tdays-dble(int(tdays))-0.5D0)
      tau1=max(-pi, min(pi, cff))
      tau2=max(-pi, min(pi, cff+2.D0*dtau))
c**   Ampl=(tau2-tau1 +sin(tau2)-sin(tau1))/dtau
      cff=tau2-tau1
      Ampl=( cff + 2.D0*cos(0.5D0*(tau1+tau2))*sin(0.5D0*cff) )/dtau


#elif defined COS2_CYCLE

! Straightforward computation of the above formula sampled at
! instantaneous moments of time.

/*
        cff=tdays-dble(int(tdays)) +0.5D0*dt*sec2day
        if (0.25D0 < cff .and. cff < 0.75D0) then
          Ampl=2.D0+2.D0*cos(2.D0*pi*(2.D0*cff-1.D0))
        else
          Ampl=0.D0
        endif
*/


!>>>>>>>>>>>>> code segment being tested >>>>>>>>>>>>>>>
      dtau=2.D0*pi * dt*sec2day
      cff=4.D0*pi*(tdays-dble(int(tdays))-0.5D0)
      tau1=max(-pi, min(pi, cff))
      tau2=max(-pi, min(pi, cff+2.D0*dtau))
c**   Ampl=(tau2-tau1 +sin(tau2)-sin(tau1))/dtau
      cff=tau2-tau1
      Ampl=( cff + 2.D0*cos(0.5D0*(tau1+tau2))*sin(0.5D0*cff) )/dtau
!<<<<<<<<<<<<< end of tested code segment <<<<<<<<<<

#endif
        write(*,*) tdays, Ampl
        sum=sum + Ampl

        if (int(tdays+dt*sec2day) > int(tdays)) then
          write(*,*) 'sum = ',sum*dt*sec2day
          sum=0.D0
        endif
      enddo
      end
