
        do i=istr,iend
          DC(i,0) =dt*pm(i,j)*pn(i,j)
          FC(i,N-1)= 0.5*dt*(Akv(i,j,N)+Akv(i,j,N-1))
     &                      /Hz(i,j,N)

          WC(i,N-1)= DC(i,0)*0.5*(Wi(i,j,N)+Wi(i,j,N-1))

          cff=1./( 0.5*(Hz(i,j,N)                )
     &                    +FC(i,N-1)-min(WC(i,N-1),0.) )

          CF(i,N-1)=cff*(  FC(i,N-1)+max(WC(i,N-1),0.) )

          DC(i,N)=cff*w(i,j,N,nnew)
        enddo

        do k=N-1,2,-1      !--> forward elimination for w
          do i=istr,iend
            FC(i,k-1)= 0.5*dt*(Akv(i,j,k)+Akv(i,j,k-1))!! flux coef between w(k) and w(k-1)
     &                   /Hz(i,j,k)

            WC(i,k-1)= DC(i,0)*0.5*(Wi(i,j,k)+Wi(i,j,k-1))

            cff=1./( 0.5*(Hz(i,j,k)+Hz(i,j,k-1))
     &                            +FC(i,k-1)-min(WC(i,k-1),0.)
     &                            +FC(i,k  )+max(WC(i,k),0.)
     &                     -CF(i,k)*(FC(i,k)-min(WC(i,k),0.))
     &                                                       )
            CF(i,k-1)=cff*( FC(i,k-1)+max(WC(i,k-1),0.) )

            DC(i,k)=cff*( w(i,j,k,nnew)
     &                 +DC(i,k+1)*(FC(i,k)-min(WC(i,k),0.)))
          enddo
        enddo
        !! change this you feel energetic. Instead of no-flux at the bottom we
        !! have w(0) = 0 (different BC for tri-diag system)
        do i=istr,iend
          w(i,j,1,nnew)=( w(i,j,k,nnew)  +DC(i,2)*(FC(i,1)-min(WC(i,1),0.))
     &                         )/( 0.5*(Hz(i,j,2)+Hz(i,j,1))
     &                                       +FC(i,1)+max(WC(i,1),0.)
     &                              -CF(i,1)*(FC(i,1)-min(WC(i,1),0.))
     &                                                               )
        enddo
        do k=2,N,+1          !--> backsubstitution for w
          do i=istr,iend
            w(i,j,k,nnew)=DC(i,k) +CF(i,k-1)*w(i,j,k-1,nnew)
          enddo
        enddo
        !------- end computing of w(:,:,:,nnew) -----
