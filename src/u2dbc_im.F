#include "cppdefs.opt"

      subroutine u2dbc(grad)

! Set lateral boundary conditions for the barotropic (i.e. vertically
! integrated) x-component velocity ubar(:,:,knew).

      use param
#ifdef WEC
      use wec_frc
#endif
      use boundary
      use grid
      use ocean_vars
      use scalars

      implicit none

      ! import/export
      real,dimension(GLOBAL_2D_ARRAY) :: grad

      ! Local
      integer i,j
      real :: cx,cy, cext, cff, dft,dfx,dfy, zx,hx
      real, parameter :: eps=1.D-33

! For open bc's: m2flather,m2specified,no-gradient

# if defined OBC_M2FLATHER
      logical,parameter ::obc_m2flather   = .true.
      logical,parameter ::obc_m2specified = .false.
# elif defined OBC_M2SPECIFIED
      logical,parameter ::obc_m2flather   = .false.
      logical,parameter ::obc_m2specified = .true.
# else   /* no-gradient at obc */
      logical,parameter ::obc_m2flather   = .false.
      logical,parameter ::obc_m2specified = .false.
# endif

      if (west_bnd) then
        if (obc_west) then
          if (obc_m2flather) then
            do j=1,ny
              cff=0.5*(h(0,j)+h(1,j))
              hx=sqrt(g/cff)
              cx=dtfast*cff*hx*0.5*(pm(0,j)+pm(1,j))

              zx=(0.5+cx)*zeta(1,j,kstp)+(0.5-cx)*zeta(0,j,kstp)
              if (cx > 0.292893218813452) then
                zx=zx + ( zeta(1,j,knew) +cx*zeta(0,j,kstp)
     &                              -(1.+cx)*zeta(1,j,kstp)
     &                      )*(1.-0.292893218813452/cx)**2
              endif

              ubar(1,j,knew)=0.5*( (1.-cx)*ubar(1,j,kstp)
     &                                 +cx*ubar(2,j,kstp)
     &                   +ubar_west(j) -hx*(zx-zeta_west(j)) )
     &                                     *umask(1,j)
            enddo

          elseif (obc_m2specified) then
            do j=1,ny
              ubar(1,j,knew)=ubar_west(j)*umask(1,j)
            enddo
          else  ! no gradient
            do j=1,ny
              ubar(1,j,knew)=ubar(2,j,knew)*umask(1,j)
            enddo
          endif
        else ! closed boundary
          do j=1,ny
            ubar(1,j,knew)=0.
          enddo
        endif   ! obc_west
      endif     ! west_bnd

      if (east_bnd) then
        if (obc_east) then
          if (obc_m2flather) then
            do j=1,ny
              cff=0.5*(h(nx,j)+h(nx+1,j))
              hx=sqrt(g/cff)
              cx=dtfast*cff*hx*0.5*(pm(nx,j)+pm(nx+1,j))

              zx=(0.5+cx)*zeta(nx,j,kstp)+(0.5-cx)*zeta(nx+1,j,kstp)
              if (cx > 0.292893218813452) then
                zx=zx + ( zeta(nx,j,knew) +cx*zeta(nx+1,j,kstp)
     &                               -(1.+cx)*zeta(nx,j,kstp)
     &                           )*(1.-0.292893218813452/cx)**2
              endif

              ubar(nx+1,j,knew)=0.5*( (1.-cx)*ubar(nx+1,j,kstp)
     &                                     +cx*ubar(nx,j,kstp)
     &                       +ubar_east(j) +hx*(zx-zeta_east(j)) )
     &                                         *umask(nx+1,j)
            enddo
          elseif (obc_m2specified) then

            ubar(nx+1,1:ny,knew)=ubar_east(1:ny)
     &                        *umask(nx+1,1:ny)
          else  ! No gradient
            ubar(nx+1,1:ny,knew)=ubar(nx,1:ny,knew)
     &                         *umask(nx+1,1:ny)
          endif
        else ! Closed boundary
          ubar(nx+1,1:ny,knew)=0.
        endif  ! obc_east
      endif  ! east_bnd


      if (south_bnd) then
        if (obc_south) then
          if (obc_m2flather) then
            ! even though this is part of the m2flather logical path,
            ! technically, we're use orlanski conditions for the boundaries
            ! where ubar is tangetial rather than normal
            do i=0,nx
              grad(i,0)=ubar(i+1,0,kstp)-ubar(i,0,kstp)
              grad(i,1)=ubar(i+1,1,kstp)-ubar(i,1,kstp)
            enddo
            do i=1,nx
              cx=-0.125*dtfast*(vbar(i,1,kstp)+vbar(i-1,1,kstp))
     &          *(pn(i,0)+pn(i-1,0)+pn(i,1)+pn(i-1,1))

              cy= 0.125*dtfast*(ubar(i,0,kstp)+ubar(i,1,kstp))
     &          *(pm(i,0)+pm(i-1,0)+pm(i,1)+pm(i-1,1))

              if (cx > 0.) then
                cext=0.
              else
                cext=-cx ; cx=0.
              endif

              ubar(i,0,knew)=(1.-cx)*( ubar(i,0,kstp)
     &                          -max(cy,0.)*grad(i-1,0)
     &                          -min(cy,0.)*grad(i  ,0)
     &                                                     )
     &                       +cx*(        ubar(i,1,kstp)
     &                            -max(cy,0.)*grad(i-1,1)
     &                            -min(cy,0.)*grad(i  ,1)
     &                                                     )

              ubar(i,0,knew)=(1.-cext)*ubar(i,0,knew)
     &                                   +cext*ubar_south(i)
              ubar(i,0,knew)=ubar(i,0,knew)*umask(i,0)
            enddo
          elseif (obc_m2specified) then
            do i=1,nx
              ubar(i,0,knew)=ubar_south(i)*umask(i,0)
            enddo
          else ! no-gradient
            do i=1,nx
              ubar(i,0,knew)=ubar(i,1,knew)*umask(i,0)
            enddo
          endif  ! obc_m2flather, obc_m2specified, or no-gradient
        else ! slip or no-slip at wall
          do i=i0,i1  
            ubar(i,0,knew)=gamma2*ubar(i,1,knew)
     &                                *umask(i,0)
          enddo
        endif   ! obc_south
      endif   ! south_bnd



      if (north_bnd) then
        if (obc_north) then
          if (obc_m2flather) then
            do i=0,nx
              grad(i,ny  )=ubar(i+1,ny  ,kstp)-ubar(i,ny,kstp  )
              grad(i,ny+1)=ubar(i+1,ny+1,kstp)-ubar(i,ny+1,kstp)
            enddo
            do i=1,nx

              cx=0.125*dtfast*(vbar(i,ny+1,kstp)+vbar(i-1,ny+1,kstp))
     &             *(pn(i,ny)+pn(i-1,ny)+pn(i,ny+1)+pn(i-1,ny+1))
  
              cy=0.125*dtfast*(ubar(i,ny,kstp)+ubar(i,ny+1,kstp))
     &             *(pm(i,ny)+pm(i-1,ny)+pm(i,ny+1)+pm(i-1,ny+1))

              if (cx > 0.) then
                cext=0.
              else
                cext=-cx ; cx=0.
              endif

              ubar(i,ny+1,knew)=(1.-cx)*( ubar(i,ny+1,kstp)
     &                            -max(cy,0.)*grad(i-1,ny+1)
     &                            -min(cy,0.)*grad(i  ,ny+1)
     &                                                     )
     &                         +cx*(        ubar(i,ny,kstp)
     &                              -max(cy,0.)*grad(i-1,ny)
     &                            -min(cy,0.)*grad(i  ,ny)
     &                                                     )

              ubar(i,ny+1,knew)=(1.-cext)*ubar(i,ny+1,knew)
     &                                  +cext*ubar_north(i)
              ubar(i,ny+1,knew)=ubar(i,ny+1,knew)*umask(i,ny+1)
            enddo
          elseif (obc_m2specified) then
            do i=1,nx
              ubar(i,ny+1,knew)=ubar_north(i)*umask(i,ny+1)
            enddo
          else  ! no-gradient
            do i=1,nx
              ubar(i,ny+1,knew)=ubar(i,ny,knew)*umask(i,ny+1)
            enddo
          endif
        else  ! slip or no-slip at closed boundary
          do i=i0,i1
            ubar(i,ny+1,knew)=gamma2*ubar(i,ny,knew)
     &                                *umask(i,ny+1)
          enddo
        endif   ! obc_north
      endif  ! north_bnd

      ! Corners between adjacent open boundaries
      if (west_bnd.and.south_bnd.and.obc_west.and.obc_south) then
        ubar(1,0,knew)=0.5*( ubar(2,0,knew)
     &                                  +ubar(1,1,knew))
      endif
      if (east_bnd.and.south_bnd.and.obc_east.and.obc_south) then
        ubar(nx+1,0,knew)=0.5*( ubar(nx,0,knew)
     &                                +ubar(nx+1,1,knew))
      endif
      if (west_bnd.and.north_bnd.and.obc_west.and.obc_north) then
        ubar(1,ny+1,knew)=0.5*( ubar(2,ny+1,knew)
     &                                  +ubar(1,ny,knew))
      endif
      if (east_bnd.and.north_bnd.and.obc_east.and.obc_north) then
        ubar(nx+1,ny+1,knew)=0.5*( ubar(nx,ny+1,knew)
     &                                +ubar(nx+1,ny,knew))
      endif

      end

