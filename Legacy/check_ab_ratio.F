      implicit none
      real Tt, Ts, Tp, ratio
      real, parameter ::
     &       A0=+0.665157E-01,  A1=+0.170907E-01,  A2=-0.203814E-03,
     &       A3=+0.298357E-05,  A4=-0.255019E-07,  B0=+0.378110E-02,
     &       B1=-0.846960E-04,  CO=-0.678662E-05,  D0=+0.380374E-04,
     &       D1=-0.933746E-06,  D2=+0.791325E-08,  E0=-0.164759E-06,
     &       F0=-0.251520E-11,  G0=+0.512857E-12,  H0=-0.302285E-13,
     &                                             Smean=35.0

      Ts=40.0 - Smean    !  Check Value: ab_ratio=0.34763 (Celsius/PSU)
      Tt=10.0            !  ab_ratio=0.34763 (Celsius/PSU)
      Tp=4000.0          !  at s=40.0, t=10.0, p=4000.0
 

               ratio=A0+Tt*(A1+Tt*(A2+Tt*(A3+Tt*A4)))
     &                 +Ts*(B0+Tt*B1+Ts*CO)
     &                 +Tp*(D0+Tt*(D1+Tt*D2)+Ts*E0
     &                        +Tp*(Ts*F0+Tt*Tt*G0+Tp*H0))

      write(*,*) ratio

      ratio = 0.34763
      write(*,*) ratio

      stop
      end
 
 
 
 
