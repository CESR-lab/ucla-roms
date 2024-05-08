      implicit none
      real beta,  A0, A1, A2, A3, A4, B0, B1, CO, D0, D1, D2, E0,
     &     alpha, F0, G0, H0, Q0, Q1, Q3, Q4, R1, R2, U0, U1, U2,
     &            V0, V1, W0, X0, Y0, Z0, Smean,  Tt, Ts, Tp, ratio
      parameter(A0=+0.665157E-01, A1=+0.170907E-01, A2=-0.203814E-03,
     &          A3=+0.298357E-05, A4=-0.255019E-07, B0=+0.378110E-02,
     &          B1=-0.846960E-04, CO=-0.678662E-05, D0=+0.380374E-04,
     &          D1=-0.933746E-06, D2=+0.791325E-08, E0=-0.164759E-06,
     &          F0=-0.251520E-11, G0=+0.512857E-12, H0=-0.302285E-13,
!-
     &          Q0=+0.785567e+01, Q1=-0.301985e-01, Q3=+0.555579e-03,
     &          Q4=-0.415613e-05, R1=-0.356603e-02, R2=+0.788212e-04,
     &          U0=-0.121555e-03, U1=+0.192867e-05, U2=-0.213127e-07,
     &          V0=+0.176621e-08, V1=-0.175379e-10, W0=+0.408195e-06,
     &          X0=+0.515032e-04, Y0=-0.602281e-11, Z0=+0.121551e-13,
     &          Smean=35.0)
!-
      Tt=10.0
      Ts=40.0 - Smean
      Tp=4000.0
!-
               beta=Q0+Tt*(Q1+Tt*(Q3+Tt*Q4))
     &                +Ts*(R1+Tt*R2+Ts*X0)
     &                +Tp*( U0+Tt*(U1+Tt*U2)+Ts*W0
     &                        +Tp*(V0+Tt*V1+Ts*Y0+Tp*Z0))
 
 
          ratio=A0+Tt*(A1+Tt*( A2+Tt*(A3+Tt*A4)))
     &            +Ts*(B0+Tt*B1+Ts*CO)
     &            +Tp*( D0+Tt*(D1+Tt*D2)+Ts*E0
     &                    +Tp*(Ts*F0+Tt*Tt*G0+Tp*H0))
 
          alpha=ratio*beta
      write(*,*) beta, ratio
 
c?? THE ORIGINAL VERSION
 
          beta=Q0+Tt*(Q1+Tt*(Q3+Tt*Q4))+
     &              Ts*(R1+Tt*R2)+
     &              Tp*(U0+Tt*(U1+Tt*U2))+
     &              Tp*Tp*(V0+Tt*V1)+
     &              Ts*Tp*W0+Ts*Ts*X0+
     &              Tp*Tp*Ts*Y0+
     &              Tp*Tp*Tp*Z0
 
          ratio=(A0+Tt*(A1+Tt*(A2+Tt*(A3+
     &                            Tt*A4)))+
     &                Ts*(B0+Tt*B1)+Ts*Ts*CO)+
     &               Tp*(D0+Tt*(D1+Tt*D2))+
     &               Tp*Ts*E0+
     &               Tp*Tp*(Ts*F0+Tt*Tt*G0+
     &                                Tp*H0)
 
 
 
          alpha=ratio*beta
 
      write(*,*) beta, ratio
      stop
      end
 
c**     7.246469      0.4151567
 
 
