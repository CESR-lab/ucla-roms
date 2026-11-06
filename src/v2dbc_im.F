#include "cppdefs.opt"

      subroutine v2dbc(grad)

! Set lateral boundary conditions for the barotropic (i.e. vertically
! integrated) y-component velocity vbar(:,:,knew).

      use param
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

! Normal faces: south & north -- v is the normal component

      if (south_bnd) then
        if (obc_south) then
          if (obc_m2flather) then
            do i=1,nx
              cff=0.5*(h(i,0)+h(i,1))
              hx=sqrt(g/cff)
              cx=dtfast*cff*hx*0.5*(pn(i,0)+pn(i,1))

              zx=(0.5+cx)*zeta(i,1,kstp)+(0.5-cx)*zeta(i,0,kstp)
              if (cx > 0.292893218813452) then
                zx=zx + ( zeta(i,1,knew) +cx*zeta(i,0,kstp)
     &                              -(1.+cx)*zeta(i,1,kstp)
     &                      )*(1.-0.292893218813452/cx)**2
              endif

              vbar(i,1,knew)=0.5*( (1.-cx)*vbar(i,1,kstp)
     &                                 +cx*vbar(i,2,kstp)
     &                  +vbar_south(i) -hx*(zx-zeta_south(i)) )
     &                                     *vmask(i,1)
            enddo

          elseif (obc_m2specified) then
            do i=1,nx
              vbar(i,1,knew)=vbar_south(i)*vmask(i,1)
            enddo
          else  ! no gradient
            do i=1,nx
              vbar(i,1,knew)=vbar(i,2,knew)*vmask(i,1)
            enddo
          endif
        else ! closed boundary
          do i=1,nx
            vbar(i,1,knew)=0.
          enddo
        endif   ! obc_south
      endif     ! south_bnd

      if (north_bnd) then
        if (obc_north) then
          if (obc_m2flather) then
            do i=1,nx
              cff=0.5*(h(i,ny)+h(i,ny+1))
              hx=sqrt(g/cff)
              cx=dtfast*cff*hx*0.5*(pn(i,ny)+pn(i,ny+1))

              zx=(0.5+cx)*zeta(i,ny,kstp)+(0.5-cx)*zeta(i,ny+1,kstp)
              if (cx > 0.292893218813452) then
                zx=zx + ( zeta(i,ny,knew) +cx*zeta(i,ny+1,kstp)
     &                               -(1.+cx)*zeta(i,ny,kstp)
     &                           )*(1.-0.292893218813452/cx)**2
              endif

              vbar(i,ny+1,knew)=0.5*( (1.-cx)*vbar(i,ny+1,kstp)
     &                                     +cx*vbar(i,ny,kstp)
     &                      +vbar_north(i) +hx*(zx-zeta_north(i)) )
     &                                         *vmask(i,ny+1)
            enddo
          elseif (obc_m2specified) then
            vbar(1:nx,ny+1,knew)=vbar_north(1:nx)
     &                        *vmask(1:nx,ny+1)
          else  ! No gradient
            vbar(1:nx,ny+1,knew)=vbar(1:nx,ny,knew)
     &                         *vmask(1:nx,ny+1)
          endif
        else ! Closed boundary
          vbar(1:nx,ny+1,knew)=0.
        endif  ! obc_north
      endif    ! north_bnd


! Tangential faces: west & east -- v is the tangential component

      if (west_bnd) then
        if (obc_west) then
          if (obc_m2flather) then
            ! even though this is part of the m2flather logical path,
            ! technically, we're using orlanski conditions for the boundaries
            ! where vbar is tangential rather than normal
            do j=0,ny
              grad(0,j)=vbar(0,j+1,kstp)-vbar(0,j,kstp)
              grad(1,j)=vbar(1,j+1,kstp)-vbar(1,j,kstp)
            enddo
            do j=1,ny
              cx=-0.125*dtfast*(ubar(1,j,kstp)+ubar(1,j-1,kstp))
     &          *(pm(0,j)+pm(0,j-1)+pm(1,j)+pm(1,j-1))

              cy= 0.125*dtfast*(vbar(0,j,kstp)+vbar(1,j,kstp))
     &          *(pn(0,j)+pn(0,j-1)+pn(1,j)+pn(1,j-1))

              if (cx > 0.) then
                cext=0.
              else
                cext=-cx ; cx=0.
              endif

              vbar(0,j,knew)=(1.-cx)*( vbar(0,j,kstp)
     &                          -max(cy,0.)*grad(0,j-1)
     &                          -min(cy,0.)*grad(0,j  )
     &                                                     )
     &                       +cx*(        vbar(1,j,kstp)
     &                            -max(cy,0.)*grad(1,j-1)
     &                            -min(cy,0.)*grad(1,j  )
     &                                                     )

              vbar(0,j,knew)=(1.-cext)*vbar(0,j,knew)
     &                                   +cext*vbar_west(j)
              vbar(0,j,knew)=vbar(0,j,knew)*vmask(0,j)
            enddo
          elseif (obc_m2specified) then
            do j=1,ny
              vbar(0,j,knew)=vbar_west(j)*vmask(0,j)
            enddo
          else ! no-gradient
            do j=1,ny
              vbar(0,j,knew)=vbar(1,j,knew)*vmask(0,j)
            enddo
          endif  ! obc_m2flather, obc_m2specified, or no-gradient
        else ! slip or no-slip at wall
          do j=1,j1
            vbar(0,j,knew)=gamma2*vbar(1,j,knew)
     &                                *vmask(0,j)
          enddo
        endif   ! obc_west
      endif     ! west_bnd

      if (east_bnd) then
        if (obc_east) then
          if (obc_m2flather) then
            do j=0,ny
              grad(nx  ,j)=vbar(nx  ,j+1,kstp)-vbar(nx  ,j,kstp)
              grad(nx+1,j)=vbar(nx+1,j+1,kstp)-vbar(nx+1,j,kstp)
            enddo
            do j=1,ny

              cx=0.125*dtfast*(ubar(nx+1,j,kstp)+ubar(nx+1,j-1,kstp))
     &           *(pm(nx,j)+pm(nx,j-1)+pm(nx+1,j)+pm(nx+1,j-1))

              cy=0.125*dtfast*(vbar(nx,j,kstp)+vbar(nx+1,j,kstp))
     &           *(pn(nx,j)+pn(nx,j-1)+pn(nx+1,j)+pn(nx+1,j-1))

              if (cx > 0.) then
                cext=0.
              else
                cext=-cx ; cx=0.
              endif

              vbar(nx+1,j,knew)=(1.-cx)*( vbar(nx+1,j,kstp)
     &                            -max(cy,0.)*grad(nx+1,j-1)
     &                            -min(cy,0.)*grad(nx+1,j  )
     &                                                     )
     &                         +cx*(        vbar(nx,j,kstp)
     &                              -max(cy,0.)*grad(nx,j-1)
     &                              -min(cy,0.)*grad(nx,j  )
     &                                                     )

              vbar(nx+1,j,knew)=(1.-cext)*vbar(nx+1,j,knew)
     &                                   +cext*vbar_east(j)
              vbar(nx+1,j,knew)=vbar(nx+1,j,knew)*vmask(nx+1,j)
            enddo
          elseif (obc_m2specified) then
            do j=1,ny
              vbar(nx+1,j,knew)=vbar_east(j)*vmask(nx+1,j)
            enddo
          else ! no-gradient
            do j=1,ny
              vbar(nx+1,j,knew)=vbar(nx,j,knew)*vmask(nx+1,j)
            enddo
          endif
        else ! slip or no-slip at closed boundary
          do j=1,j1
            vbar(nx+1,j,knew)=gamma2*vbar(nx,j,knew)
     &                                *vmask(nx+1,j)
          enddo
        endif   ! obc_east
      endif     ! east_bnd

      ! Corners between adjacent open boundaries
      if (west_bnd.and.south_bnd.and.obc_west.and.obc_south) then
        vbar(0,1,knew)=0.5*( vbar(0,2,knew)
     &                                  +vbar(1,1,knew))
      endif
      if (east_bnd.and.south_bnd.and.obc_east.and.obc_south) then
        vbar(nx+1,1,knew)=0.5*( vbar(nx+1,2,knew)
     &                                +vbar(nx,1,knew))
      endif
      if (west_bnd.and.north_bnd.and.obc_west.and.obc_north) then
        vbar(0,ny+1,knew)=0.5*( vbar(0,ny,knew)
     &                                  +vbar(1,ny+1,knew))
      endif
      if (east_bnd.and.north_bnd.and.obc_east.and.obc_north) then
        vbar(nx+1,ny+1,knew)=0.5*( vbar(nx+1,ny,knew)
     &                                +vbar(nx,ny+1,knew))
      endif

      end subroutine v2dbc
