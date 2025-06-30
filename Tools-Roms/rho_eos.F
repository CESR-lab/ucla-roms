      subroutine rho_eos (Lm,Mm,N, T,S, z_r,z_w, rp_r,rp_w)


#define SPLIT_EOS
#define NONLIN_EOS
#define DUKO_2001

      implicit none
      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k
      real*4   T(0:Lm+1,0:Mm+1,N),    S(0:Lm+1,0:Mm+1,N)
      real*8 z_r(0:Lm+1,0:Mm+1,N),  z_w(0:Lm+1,0:Mm+1,0:N),
     &      rp_r(0:Lm+1,0:Mm+1,N), rp_w(0:Lm+1,0:Mm+1,0:N)

      real*8 rho1(0:Lm+1,N), qp1(0:Lm+1,N), bvf(0:Lm+1,0:N)

      real*8, parameter ::     r00=999.842594D0,    r01=6.793952D-2,
     &      r02=-9.095290D-3,   r03=1.001685D-4,    r04=-1.120083D-6,
     &                                              r05=6.536332D-9,
     &      r10=0.824493D0,     r11=-4.08990D-3,    r12=7.64380D-5,
     &                          r13=-8.24670D-7,    r14=5.38750D-9,
     &      rS0=-5.72466D-3,    rS1=1.02270D-4,     rS2=-1.65460D-6,
     &      r20=4.8314D-4,

     &      K00=19092.56D0,     K01=209.8925D0,     K02=-3.041638D0,
     &                          K03=-1.852732D-3,   K04=-1.361629D-5,
     &      K10=104.4077D0,     K11=-6.500517D0,    K12=0.1553190D0,
     &                                              K13=2.326469D-4,
     &      KS0=-5.587545D0,    KS1=+0.7390729D0,   KS2=-1.909078D-2

#  ifndef SPLIT_EOS
      real*8, parameter ::      B00=0.4721788D0,
     &      B01=0.01028859D0,   B02=-2.512549D-4,   B03=-5.939910D-7,
     &      B10=-0.01571896D0,  B11=-2.598241D-4,   B12=7.267926D-6,
     &                          BS1=2.042967D-3,

     &      E00=+1.045941D-5,   E01=-5.782165D-10,  E02=+1.296821D-7,
     &      E10=-2.595994D-7,   E11=-1.248266D-9,   E12=-3.508914D-9
#  endif
      real*8, parameter :: qp2=0.0000172
      real*8 dpth,cff, Tt,Ts,sqrtTs, rho0, K0, dr00

      integer numthreads, trd, chunk_size, margin, jstr,jend
C$    integer omp_get_num_threads, omp_get_thread_num



#  ifdef DUKO_2001
      real*8 K0_Duk, rho1_0

      Tt=3.8D0
      Ts=34.5D0
      sqrtTs=sqrt(Ts)

      cff=r00-1000.D0
      rho1_0=cff +Tt*( r01+Tt*( r02+Tt*( r03+Tt*( r04+Tt*r05 ))))
     &                            +Ts*( R10+Tt*( r11+Tt*( r12+Tt*(
     &                                              r13+Tt*r14 )))
     &                   +sqrtTs*( rS0+Tt*( rS1+Tt*rS2 ))+Ts*r20 )

      K0_Duk= Tt*( K01+Tt*( K02+Tt*( K03+Tt*K04 )))
     &       +Ts*( K10+Tt*( K11+Tt*( K12+Tt*K13 ))
     &            +sqrtTs*( KS0+Tt*( KS1+Tt*KS2 )))

c**   write(*,'(2(4x,A,F16.8))') 'rho1_0=',rho1_0, 'K0_Duk=',K0_Duk
      rho0=rho1_0
#  endif
      dr00=r00-1000.D0-rho0



      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1


      numthreads=1
C$    numthreads=omp_get_num_threads()
      trd=0
C$    trd=omp_get_thread_num()
      chunk_size=(jmax-jmin + numthreads)/numthreads
      margin=(chunk_size*numthreads -jmax+jmin-1)/2
      jstr=max(  trd   *chunk_size  -margin, jmin )
      jend=min( (trd+1)*chunk_size-1-margin, jmax )

c**    write(*,*) 'trd,jstr,jend =', trd,jstr,jend



      do j=jstr,jend
        do k=1,N
          do i=imin,imax
            Tt=T(i,j,k)
            Ts=max(S(i,j,k),0.)
            sqrtTs=sqrt(Ts)

            rho1(i,k)=( dr00 +Tt*( r01+Tt*( r02+Tt*( r03+Tt*(
     &                                         r04+Tt*r05 ))))

     &                       +Ts*( r10+Tt*( r11+Tt*( r12+Tt*(
     &                                          r13+Tt*r14 )))
     &                            +sqrtTs*(rS0+Tt*(
     &                                 rS1+Tt*rS2 ))+Ts*r20 ))
#  ifdef MASKING
     &                                            *rmask(i,j)
#  endif


            K0= Tt*( K01+Tt*( K02+Tt*( K03+Tt*K04 )))
     &         +Ts*( K10+Tt*( K11+Tt*( K12+Tt*K13 ))
     &              +sqrtTs*( KS0+Tt*( KS1+Tt*KS2 )))

            qp1(i,k)= 0.1D0*(rho0+rho1(i,k))*(K0_Duk-K0)
     &                          /((K00+K0)*(K00+K0_Duk))
#   ifdef MASKING
     &                   *rmask(i,j)
#   endif
          enddo
        enddo

        cff=1.D0
        do k=1,N-1
          do i=imin,imax
            dpth=z_w(i,j,N)-0.5*(z_r(i,j,k+1)+z_r(i,j,k))

            bvf(i,k)=-cff*( rho1(i,k+1)-rho1(i,k)   ! Elementary
     &                      +(qp1(i,k+1)-qp1(i,k))  ! adiabatic
     &                      *dpth*(1.-2.*qp2*dpth)  ! difference
     &                 )/(z_r(i,j,k+1)-z_r(i,j,k))
          enddo
        enddo
        do i=imin,imax
          bvf(i,N)=bvf(i,N-1)
          bvf(i,0)=bvf(i,1)
          rp_w(i,j,N)=rho1(i,N) + (rho1(i,N)-rho1(i,N-1))
     &                           *(z_w(i,j,N)-z_r(i,j,N))
     &                         /(z_r(i,j,N)-z_r(i,j,N-1))

          rp_r(i,j,N)=rp_w(i,j,N) + bvf(i,N)*(z_w(i,j,N)-z_r(i,j,N))
        enddo

        do k=N-1,1,-1
          do i=imin,imax
            rp_r(i,j,k)=rp_r(i,j,k+1) +bvf(i,k)*( z_r(i,j,k+1)
     &                                             -z_r(i,j,k))

            rp_w(i,j,k)=( rp_r(i,j,k)*(z_r(i,j,k+1)-z_w(i,j,k))
     &                   +rp_r(i,j,k+1)*(z_w(i,j,k)-z_r(i,j,k))
     &                              )/(z_r(i,j,k+1)-z_r(i,j,k))
          enddo
        enddo

        do i=imin,imax
          rp_w(i,j,0)=rp_r(i,j,1) +bvf(i,1)*(z_r(i,j,1)-z_w(i,j,0))
        enddo
      enddo ! <-- j
      return
      end
