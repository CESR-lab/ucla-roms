#include "cppdefs.h"
#ifdef SOLVE3D
                                         ! Setup S-coordinate system:
      subroutine set_scoord              ! input: theta_s, theta_b, hc
      implicit none                      ! output: Cs_w(0:N), Cs_r(1:N)
# include "param.h"
# include "scoord.h"
# include "scalars.h"
# include "comm_vars.h"
      integer k                          ! Compute vertical stretching
      real ds,sc, z1,zhc,z2,z3, CSF      ! curves at W- and RHO-points,

      ds=1.D0/dble(N)                    !      -1 < Cs_r,Cs_w < 0,
      Cs_w(N)=0.D0
      do k=N-1,1,-1                      ! then print out z-coordinate
        sc=ds*dble(k-N)                  ! of S-levels for minimum, 
        Cs_w(k)=CSF(sc, theta_s,theta_b) ! h=hc, half way, and maximum 
      enddo                              ! depth.
      Cs_w(0)=-1.D0
      do k=1,N
        sc=ds*(dble(k-N)-0.5D0)
        Cs_r(k)=CSF(sc, theta_s,theta_b)
      enddo


      mpi_master_only write(*,'(/1x,A/,/2x,A,7x,A/)')
     &        'Vertical S-coordinate system (z at W-points):',
     &             'level   S-coord    Cs-curve    Z at hmin',
     &                        'at hc    half way     at hmax'
      do k=N,0,-1
        sc=ds*dble(k-N)
# if defined VERT_COORD_TYPE_SM09 || defined VERT_COORD_TYPE_KAU
        z1=hmin*(hc*sc + hmin*Cs_w(k))/(hc+hmin)
        zhc=0.5*hc*(sc + Cs_w(k))
        z2=0.5*hmax*(hc*sc + 0.5*hmax*Cs_w(k))/(hc+0.5*hmax)
        z3=hmax*(hc*sc + hmax*Cs_w(k))/(hc+hmax)
# else
        zhc=hc*sc
        z1=zhc +     (hmin-hc)*Cs_w(k)
        z2=zhc + (0.5*hmax-hc)*Cs_w(k)
        z3=zhc +     (hmax-hc)*Cs_w(k)
# endif
        if (hc < 1.E+4) then
          mpi_master_only write(*,'(I7,F11.6,F12.7,4F12.3)')
     &                              k, ds*(k-N),Cs_w(k), z1,zhc,z2,z3
        else
         mpi_master_only write(*,'(I7,F11.6,F12.7,F12.3,12x,2F12.3)')
     &                              k, ds*(k-N),Cs_w(k), z1,    z2,z3
        endif
      enddo
# if ! defined VERT_COORD_TYPE_SM09 && ! defined VERT_COORD_TYPE_KAU
      if (hc > hmin) then
        write(*,'(/1x,2A,F7.2/26x,A,F7.2/)') '### ERROR: set_scoord ',
     &         ':: Specified S-coordinate critical depth  hc   =', hc,
     &          'exceeds minimum unmasked topography.  hmin =', hmin
        may_day_flag=9
      endif
# endif
      end

! Vertical stretching functions: In principle stretching curves can be
! selected independently from the vertical transformation type, which
! would require a separate name for CPP-switch here.  This however has
! little practical value (merely to reproduce my own Pacific solutions
! of 2005 KPP presentation where the SM09 coordinate transformation
! formula was used in combination with legacy SH94 stretching curves),
! so for simplicity selection of SM09 stretching curve is tied with
! SM09 coordinate transform.



# if defined VERT_COORD_TYPE_SM09 || defined VERT_COORD_TYPE_KAU
                                              ! Note that mathematical
      function CSF(sc, theta_s,theta_b)       ! limits of csrf,CSF for
      implicit none                           ! theta_s, theta_b --> 0
      real*8 CSF, sc, theta_s,theta_b,csrf    ! match that under "else"
                                              ! logical branches.
      if (theta_s > 0.D0) then
        csrf=(1.D0-cosh(theta_s*sc))/(cosh(theta_s)-1.D0)
      else
        csrf=-sc**2                           ! Reference: This form of
      endif                                   ! CSF corresponds exactly
      if (theta_b > 0.D0) then
        CSF=(exp(theta_b*csrf)-1.D0)/(1.D0-exp(-theta_b))
      else
        CSF=csrf                              ! to Eq.(2.4) from SM2009 
      endif                                   ! article published in
      end                                     ! J. Comput. Phys.


# elif defined THIS_SHOULD_BE_NEVER_DEFINED
                                              ! Alternative variants
      function CSF (sc, theta_s,theta_b)      ! for reference only
      implicit none
      real*8 CSF, sc, theta_s,theta_b, sc1, csrf,cbot, wgth
      real*8, parameter :: alfa=1.D0, beta=1.D0, albt=alfa/beta

      if (theta_s > 0.D0) then
        csrf=(1.D0-cosh(theta_s*sc))/(cosh(theta_s)-1.D0)
      else
        csrf=-sc**2
      endif

      sc1=sc+1.D0
      if (theta_b > 0.D0) then
        cbot=sinh(theta_b*sc1)/sinh(theta_b)-1.D0
      else
        cbot=sc
      endif

      wgth=sc1**alfa * ( 1.D0+ albt*(1.D0-sc1**beta) )
      CSF=wgth*csrf + (1.D0-wgth)*cbot
      end


      function CSF (sc, theta_s,theta_b)
      implicit none
      real*8 CSF, sc, theta_s,theta_b, sc1, csrf,cbot, wgth
      real*8, parameter :: alfa=1.D0, beta=1.D0, albt=alfa/beta

      if (theta_s > 0.D0) then
        csrf=(1.D0-cosh(theta_s*sc))/(cosh(theta_s)-1.D0)
        if (theta_b > 0.D0) then
          sc1=sc+1.D0
          cbot=sinh(theta_b*sc1)/sinh(theta_b)-1.D0
          wgth=sc1**alfa *( 1.D0 +  albt*( 1.D0-sc1**beta ))
          CSF=wgth*csrf + (1.D0-wgth)*cbot
        else
          CSF=csrf
        endif
      else
        CSF=sc
      endif
      end

# else
      function CSF (sc, theta_s,theta_b)      ! Legacy S-coordinate
      real*8 CSF, sc, theta_s,theta_b         ! stretching from SH94
      if (theta_s > 0.D0) then
        CSF=(1.D0-theta_b)*sinh(theta_s*sc)/sinh(theta_s)
     &           +theta_b*( 0.5D0*tanh(theta_s*(sc+0.5D0))
     &                       /tanh(0.5*theta_s)   -0.5D0 )
      else
        CSF=sc
      endif
      end
# endif

      subroutine check_scoord_switches(ierr)
      implicit none
      integer ierr, is,ie, lenstr
# include "param.h"
# include "strings.h"
      ie=lenstr(cpps)
      is=ie+2 ; ie=is+10
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='<scoord.h>'
# ifdef VERT_COORD_TYPE_SM09
      is=ie+2 ; ie=is+19
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='VERT_COORD_TYPE_SM09'
# endif
# ifdef VERT_COORD_TYPE_KAU
      is=ie+2 ; ie=is+18
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='VERT_COORD_TYPE_KAU'
# endif
      return
  99  mpi_master_only write(*,'(/1x,2A/12x,A/)')      '### ERROR: ',
     &  'Insufficient length of string "cpps" in file "strings.h".',
     &        'Increase parameter "max_opt_size" it and recompile.'
      ierr=ierr+1
      end
#else
      subroutine set_scoord_empty
      end
#endif /* SOLVE3D */
