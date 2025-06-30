      implicit none
      real A00, A01, A02, A03, A04, A10, A11, A12, A13, Tt,
     &          AS0, AS1, AS2, B00, B01, B02, B03, B10, B11,  Ts,
     &          B12, B13, BS1, E00, E01, E02, E10, E11, E12, sqrtTs
      parameter (A00=+19092.56 ,  A01=+209.8925   , A02=-3.041638 ,
     &          A03=-1.852732e-3, A04=-1.361629e-5, A10=104.4077   ,
     &          A11=-6.500517   , A12=+0.1553190  , A13=2.326469e-4 ,
     &          AS0=-5.587545   , AS1=+0.7390729  , AS2=-1.909078e-2,
     &          B00=+4.721788e-1, B01=+1.028859e-2, B02=-2.512549e-4,
     &          B03=-5.939910e-7, B10=-1.571896e-2, B11=-2.598241e-4,
     &          B12=+7.267926e-6, BS1=+2.042967e-3,
     &          E00=+1.045941e-5, E01=-5.782165e-10,E02=+1.296821e-7,
     &          E10=-2.595994e-7, E11=-1.248266e-9, E12=-3.508914e-9)
      real QR,  Q01, Q02, Q03, Q04, Q05, Q10, Q11,
     &          Q12, Q13, Q14, QS0, QS1, QS2, Q20
      parameter (QR=+999.842594 , Q01=+6.793952e-2, Q02=-9.095290e-3,
     &          Q03=+1.001685e-4, Q04=-1.120083e-6, Q05=+6.536332e-9,
     &          Q10=+0.824493   , Q11=-4.08990e-3 , Q12=+7.64380e-5 ,
     &          Q13=-8.24670e-7 , Q14=+5.38750e-9 , QS0=-5.72466e-3 ,
     &          QS1=+1.02270e-4 , QS2=-1.65460e-6 , Q20=+4.8314e-4 )
      real rho0
      parameter (rho0=1000.)
!-----------------------
      integer i,j,k
      parameter (i=1,j=1,k=1)
      real rho(i,j,k), z_r(i,j,k), rho1(i,k), K0(i,k),K1(i,k),K2(i,k)

      Ts=40.0
      Tt=10.0
      z_r(i,j,k)=-4000.0
      call original (Ts,Tt,z_r)

              sqrtTs=sqrt(Ts)
!----------------------
              K0(i,k)=A00+Tt*(A01+Tt*(A02+Tt*(A03+Tt*A04)))
     &                   +Ts*(A10+Tt*(A11+Tt*(A12+Tt*A13))
     &                       +sqrtTs*(AS0+Tt*(AS1+Tt*AS2)))
              K1(i,k)=B00+Tt*(B01+Tt*(B02+Tt*B03))
     &                   +Ts*(B10+Tt*(B11+Tt*B12)+sqrtTs*BS1)
              K2(i,k)=E00+Tt*(E01+Tt*E02)
     &                   +Ts*(E10+Tt*(E11+Tt*E12))
              rho1(i,k)=QR+Tt*(Q01+Tt*(Q02+Tt*(Q03+Tt*(Q04+Tt*Q05))))
     &                    +Ts*(Q10+Tt*(Q11+Tt*(Q12+Tt*(Q13+Tt*Q14)))
     &                          +sqrtTs*(QS0+Tt*(QS1+Tt*QS2))+Ts*Q20)
!-
              rho(i,j,k)=rho1(i,k)/( 1.0+0.1*z_r(i,j,k)/( K0(i,k)
     &                       -z_r(i,j,k)*(K1(i,k)-z_r(i,j,k)*K2(i,k))
     &                                                       ))-rho0
!---------------------

      write(*,*) rho(i,j,k)
      stop
      end
!----------------
      subroutine original (Ts,Tt,z_r)
      implicit none
      real A0, A1, A2, A3, A4, B0, B1, B2, B3, D0, D1, D2, E0, E1,
     &     E2, E3, F0, F1, F2, G0, G1, G2, G3, H0, H1, H2, Q0, Q1,
     &     Q3, Q4, Q5, Q6, U0, U1, U2, U3, U4, V0, V1, V2, W0

      parameter (A0=+19092.56   , A1=+209.8925   , A2=+3.041638  ,
     &           A3=-1.852732e-3, A4=+1.361629e-5, B0=+104.4077   ,
     &           B1=+6.500517   , B2=+0.1553190  , B3=-2.326469e-4,
     &           D0=-5.587545   , D1=+0.7390729  , D2=+1.909078e-2,
     &           E0=+4.721788e-1, E1=+1.028859e-2, E2=-2.512549e-4,
     &           E3=+5.939910e-7, F0=-1.571896e-2, F1=+2.598241e-4,
     &           F2=-7.267926e-6, G0=+2.042967e-3, G1=+1.045941e-5,
     &           G2=5.782165e-10, G3=+1.296821e-7, H0=-2.595994e-7,
     &           H1=-1.248266e-9, H2=-3.508914e-9)
      parameter (Q0=+999.842594 , Q1=+6.793952e-2, Q3=-9.095290e-3,
     &           Q4=+1.001685e-4, Q5=-1.120083e-6, Q6=+6.536332e-9,
     &           U0=+0.824493   , U1=-4.08990e-3 , U2=+7.64380e-5 ,
     &           U3=-8.24670e-7 , U4=+5.38750e-9 , V0=-5.72466e-3 ,
     &           V1=+1.02270e-4 , V2=-1.65460e-6 , W0=+4.8314e-4 )

      real rho0
      parameter (rho0=1000.)
      integer i,j,k
      parameter (i=1,j=1,k=1)
      real dena(i,j,k), a3d(i,j,k),z_r(i,j,k),Ts(i,j),Tt(i,j),Tz(i,j)


            Tz(i,j)=z_r(i,j,k)

            a3d(i,j,k)=A0+
     &                 Tt(i,j)*(A1-Tt(i,j)*(A2-Tt(i,j)*(A3-
     &                          Tt(i,j)*A4)))+
     &                 Ts(i,j)*(B0-Tt(i,j)*(B1-Tt(i,j)*(B2-
     &                          Tt(i,j)*B3)))+
     &                 sqrt(Ts(i,j)*Ts(i,j)*Ts(i,j))*(D0+Tt(i,j)*(D1-
     &                                                Tt(i,j)*D2))-
     &                 Tz(i,j)*(E0+Tt(i,j)*(E1+Tt(i,j)*(E2-
     &                          Tt(i,j)*E3)))-
     &                 Tz(i,j)*Ts(i,j)*(F0-Tt(i,j)*(F1+Tt(i,j)*F2))-
     &                 Tz(i,j)*sqrt(Ts(i,j)*Ts(i,j)*Ts(i,j))*G0+
     &                 Tz(i,j)*Tz(i,j)*(G1-Tt(i,j)*(G2-Tt(i,j)*G3))+
     &                 Tz(i,j)*Tz(i,j)*Ts(i,j)*(H0+Tt(i,j)*(H1+
     &                                          Tt(i,j)*H2))


            dena(i,j,k)=(Q0+
     &                   Tt(i,j)*(Q1+Tt(i,j)*(Q3+Tt(i,j)*(Q4+
     &                            Tt(i,j)*(Q5+Tt(i,j)*Q6))))+
     &                   Ts(i,j)*(U0+Tt(i,j)*(U1+Tt(i,j)*(U2+
     &                            Tt(i,j)*(U3+Tt(i,j)*U4))))+
     &                   sqrt(Ts(i,j)*Ts(i,j)*Ts(i,j))*(V0+Tt(i,j)*
     &                                            (V1+Tt(i,j)*V2))+
     &                   W0*Ts(i,j)*Ts(i,j)) /
     &                  (1.0+0.1*Tz(i,j)/a3d(i,j,k)) - rho0

      write(*,*) dena(i,j,k)
      end






