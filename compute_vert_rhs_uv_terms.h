#ifdef UV_ADV

! Compute and add in vertical advection terms:

# define SPLINE_UV
c--# define NEUMANN_UV

# ifdef SPLINE_UV
        do i=istrU,iend
          DC(i,1)=0.5625*(Hz(i  ,j,1)+Hz(i-1,j,1))
     &           -0.0625*(Hz(i+1,j,1)+Hz(i-2,j,1))
#  if defined NEUMANN_UV
          CF(i,1)=0.5 ;  FC(i,0)=1.5*u(i,j,1,nrhs)
#  else
          CF(i,1)=1.  ;  FC(i,0)=2.0*u(i,j,1,nrhs)
#  endif
        enddo
        do k=1,N-1,+1    !--> recursive
          do i=istrU,iend
            DC(i,k+1)=0.5625*(Hz(i  ,j,k+1)+Hz(i-1,j,k+1))
     &               -0.0625*(Hz(i+1,j,k+1)+Hz(i-2,j,k+1))

            cff=1./(2.*DC(i,k)+DC(i,k+1)*(2.-CF(i,k)))
            CF(i,k+1)=cff*DC(i,k)
            FC(i,k)=cff*( 3.*( DC(i,k  )*u(i,j,k+1,nrhs)
     &                        +DC(i,k+1)*u(i,j,k  ,nrhs))
     &                              -DC(i,k+1)*FC(i,k-1))
          enddo
        enddo               !--> discard DC, keep CF,FC
        do i=istrU,iend
#  if defined NEUMANN_UV
          FC(i,N)=(3.*u(i,j,N,nrhs)-FC(i,N-1))/(2.-CF(i,N))
#  else
          FC(i,N)=(2.*u(i,j,N,nrhs)-FC(i,N-1))/(1.-CF(i,N))
#  endif
          DC(i,N)=0.        !<-- uppermost W*U flux
        enddo
        do k=N-1,1,-1       !--> recursive
          do i=istrU,iend
            FC(i,k)=FC(i,k)-CF(i,k+1)*FC(i,k+1)

#  ifdef MASKING
            DC(i,k)=FC(i,k) * 0.5*( We(i,j,k)+We(i-1,j,k) -0.125*(
     &                        (We(i+1,j,k)-We(i  ,j,k))*umask(i+1,j)
     &                       -(We(i-1,j,k)-We(i-2,j,k))*umask(i-1,j)
     &                                                          ))
#  else
            DC(i,k)=FC(i,k)*( 0.5625*(We(i  ,j,k)+We(i-1,j,k))
     &                       -0.0625*(We(i+1,j,k)+We(i-2,j,k)))
#  endif

            ru(i,j,k+1)=ru(i,j,k+1) -DC(i,k+1)+DC(i,k)
          enddo
        enddo                       !--> discard CF,FC
        do i=istrU,iend
          ru(i,j,1)=ru(i,j,1) -DC(i,1)
        enddo                          !--> discard DC
# else
        do k=2,N-2
          do i=istrU,iend
            FC(i,k)=( 0.5625*(u(i,j,k  ,nrhs)+u(i,j,k+1,nrhs))
     &               -0.0625*(u(i,j,k-1,nrhs)+u(i,j,k+2,nrhs)))
     &                      *( 0.5625*(We(i  ,j,k)+We(i-1,j,k))
     &                        -0.0625*(We(i+1,j,k)+We(i-2,j,k)))
          enddo
        enddo
        do i=istrU,iend
          FC(i,N)=0.
          FC(i,N-1)=( 0.5625*(u(i,j,N-1,nrhs)+u(i,j,N,nrhs))
     &                 -0.0625*(u(i,j,N-2,nrhs)+u(i,j,N,nrhs)))
     &                 *( 0.5625*(We(i  ,j,N-1)+We(i-1,j,N-1))
     &                   -0.0625*(We(i+1,j,N-1)+We(i-2,j,N-1)))

          FC(i,  1)=( 0.5625*(u(i,j,  1,nrhs)+u(i,j,2,nrhs))
     &                 -0.0625*(u(i,j,  1,nrhs)+u(i,j,3,nrhs)))
     &                     *( 0.5625*(We(i  ,j,1)+We(i-1,j,1))
     &                       -0.0625*(We(i+1,j,1)+We(i-2,j,1)))
          FC(i,0)=0.
        enddo
c*      do k=1,N-1
c*        do i=istrU,iend
c*          FC(i,k)=0.25*(u(i,j,k,nrhs)+u(i,j,k+1,nrhs))
c*     &                        *(We(i,j,k)+We(i-1,j,k))
c*        enddo
c*      enddo
c*      do i=istrU,iend
c*        FC(i,0)=0.
c*        FC(i,N)=0.
c*      enddo
        do k=1,N
          do i=istrU,iend
            ru(i,j,k)=ru(i,j,k)-FC(i,k)+FC(i,k-1)
          enddo
        enddo               !--> discard FC
# endif

        if (j >= jstrV) then
# ifdef SPLINE_UV
          do i=istr,iend
            DC(i,1)=0.5625*(Hz(i  ,j,1)+Hz(i,j-1,1))
     &             -0.0625*(Hz(i,j+1,1)+Hz(i,j-2,1))
#  if defined NEUMANN_UV
            CF(i,1)=0.5 ;  FC(i,0)=1.5*v(i,j,1,nrhs)
#  else
            CF(i,1)=1.  ;  FC(i,0)=2.0*v(i,j,1,nrhs)
#  endif
          enddo
          do k=1,N-1,+1       !--> recursive
            do i=istr,iend
              DC(i,k+1)=0.5625*(Hz(i  ,j,k+1)+Hz(i,j-1,k+1))
     &                 -0.0625*(Hz(i,j+1,k+1)+Hz(i,j-2,k+1))

              cff=1./(2.*DC(i,k)+DC(i,k+1)*(2.-CF(i,k)))
              CF(i,k+1)=cff*DC(i,k)
              FC(i,k)=cff*( 3.*( DC(i,k  )*v(i,j,k+1,nrhs)
     &                          +DC(i,k+1)*v(i,j,k  ,nrhs))
     &                                -DC(i,k+1)*FC(i,k-1))
            enddo
          enddo               !--> discard DC, keep CF,FC
          do i=istr,iend
#  if defined NEUMANN_UV
            FC(i,N)=(3.*v(i,j,N,nrhs)-FC(i,N-1))/(2.-CF(i,N))
#  else
            FC(i,N)=(2.*v(i,j,N,nrhs)-FC(i,N-1))/(1.-CF(i,N))
#  endif
            DC(i,N)=0.        !<-- uppermost W*V flux
          enddo
          do k=N-1,1,-1       !--> recursive
            do i=istr,iend
              FC(i,k)=FC(i,k)-CF(i,k+1)*FC(i,k+1)

#  ifdef MASKING
              DC(i,k)=FC(i,k) * 0.5*( We(i,j,k)+We(i,j-1,k) -0.125*(
     &                         (We(i,j+1,k)-We(i,j  ,k))*vmask(i,j+1)
     &                        -(We(i,j-1,k)-We(i,j-2,k))*vmask(i,j-1)
     &                                                           ))
#  else
              DC(i,k)=FC(i,k)*( 0.5625*(We(i,j  ,k)+We(i,j-1,k))
     &                         -0.0625*(We(i,j+1,k)+We(i,j-2,k)))
#  endif

              rv(i,j,k+1)=rv(i,j,k+1) -DC(i,k+1)+DC(i,k)
            enddo
          enddo               !--> discard CF,FC
          do i=istr,iend
            rv(i,j,1)=rv(i,j,1) -DC(i,1)
          enddo                         !--> discard DC

# else
          do k=2,N-2
            do i=istr,iend
              FC(i,k)=( 0.5625*(v(i,j,k ,nrhs)+v(i,j,k+1,nrhs))
     &                 -0.0625*(v(i,j,k-1,nrhs)+v(i,j,k+2,nrhs)))
     &                       *( 0.5625*(We(i,j  ,k)+We(i,j-1,k))
     &                         -0.0625*(We(i,j+1,k)+We(i,j-2,k)))
            enddo
          enddo
          do i=istr,iend
            FC(i,N)=0.
            FC(i,N-1)=(  0.5625*(v(i,j,N-1,nrhs)+v(i,j,N,nrhs))
     &                  -0.0625*(v(i,j,N-2,nrhs)+v(i,j,N,nrhs)))
     &                  *( 0.5625*(We(i,j  ,N-1)+We(i,j-1,N-1))
     &                    -0.0625*(We(i,j+1,N-1)+We(i,j-2,N-1)))

            FC(i,  1)=(  0.5625*(v(i,j,  1,nrhs)+v(i,j,2,nrhs))
     &                  -0.0625*(v(i,j,  1,nrhs)+v(i,j,3,nrhs)))
     &                      *( 0.5625*(We(i,j  ,1)+We(i,j-1,1))
     &                        -0.0625*(We(i,j+1,1)+We(i,j-2,1)))
            FC(i,0)=0.
          enddo
c*        do k=1,N-1
c*          do i=istr,iend
c*            FC(i,k)=0.25*(v(i,j,k,nrhs)+v(i,j,k+1,nrhs))
c*     &                          *(We(i,j,k)+We(i,j-1,k))
c*          enddo
c*        enddo
c*        do i=istr,iend
c*          FC(i,0)=0.
c*          FC(i,N)=0.
c*        enddo
          do k=1,N
            do i=istr,iend
              rv(i,j,k)=rv(i,j,k)-FC(i,k)+FC(i,k-1)
            enddo
          enddo
# endif
        endif !<-- j >= jstrV
#endif /* UV_ADV */
