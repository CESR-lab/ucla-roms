#include "cppdefs.opt"
#ifdef SOLVE3D
                                         ! Setup S-coordinate system:
      subroutine set_scoord              ! input: theta_s, theta_b, hc

      use param
      use comm_vars
      use scalars
      use scoord

      implicit none                      ! output: Cs_w(0:N), Cs_r(1:N)

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


      mpi_master_only write(*,'(/1x,A/)')
     &        'Vertical S-coordinate system (z at W-points):'

      mpi_master_only write(*,
     &     '(A7,2x,A11,2x,A14,2x,A16,2x,A14,2x,A14,2x,A14)')
     &     'level', 'S-coord', 'Cs-curve', 'Z at hmin',
     &     'at hc', 'half way', 'at hmax'


      do k=N,0,-1
        sc=ds*dble(k-N)
        z1=hmin*(hc*sc + hmin*Cs_w(k))/(hc+hmin)
        zhc=0.5*hc*(sc + Cs_w(k))
        z2=0.5*hmax*(hc*sc + 0.5*hmax*Cs_w(k))/(hc+0.5*hmax)
        z3=hmax*(hc*sc + hmax*Cs_w(k))/(hc+hmax)

        mpi_master_only write(*,
     &'(I7,2X,ES11.4,2X,ES14.7,2X,ES16.7,2X,ES14.7,2X,ES14.7,2X,ES14.7)')
     &       k, ds*(k-N), Cs_w(k), z1, zhc, z2, z3
      enddo
      end


      function CSF(sc, theta_s,theta_b)
      implicit none
      real*8 CSF, sc, theta_s,theta_b,csrf

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

#else
      subroutine set_scoord_empty
      end
#endif /* SOLVE3D */
