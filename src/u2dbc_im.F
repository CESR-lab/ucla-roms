#include "cppdefs.opt"

      subroutine u2dbc_tile (grad)

! Set lateral boundary conditions for the barotropic (i.e. vertically
! integrated) XI-component velocity ubar(:,:,knew).

      use param
#ifdef WEC
      use wec_frc
#endif
      use boundary
      use grid
      use ocean_vars
      use scalars
      use hidden_mpi_vars

      implicit none

      ! in 
      real, dimension(1-bf:nx +bf,1-bf:ny +bf) :: grad

      ! local
      integer i,j
      real :: cx,cy,cext,cff,dft,dfx,dfy,zx,hx
      real, parameter :: eps=1.D-33

      if (inode==0.and.(.not.ew_periodic)) then            ! Western edge
        if (obc_west.and.obc_m2flather) then               ! open boundary conditions
                                                           ! (m2flather)
          do j=1,ny
            cff=0.5*(h(0,j)+h(1,j))
            hx=sqrt(g/cff)
            cx=dtfast*cff*hx*0.5*(pm(0,j)+pm(1,j))
            zx=(0.5+cx)*zeta(1,j,kstp)+(0.5-cx)*zeta(0,j,kstp)
            if (cx > 0.292893218813452) then
              zx=zx + ( zeta(1,j,knew) +cx*zeta(0,j,kstp)
     &                            -(1.+cx)*zeta(1,j,kstp)
     &                        )*(1.-0.292893218813452/cx)**2
            endif
            ubar(1,j,knew)=0.5*( (1.-cx)*ubar(1,j,kstp)
     &                               +cx*ubar(2,j,kstp)
     &             +ubar_west(j) -hx*(zx-zeta_west(j)))
     &                                      *umask(1,j)
          enddo

        elseif (obc_west.and.(.not.obc_m2flather)) then    ! zero gradient OBC
          do j=1,nx
            ubar(1,j,knew)=ubar(2,j,knew)
          enddo

        else                                               ! closed boundary conditions
          do j=1,ny
            ubar(1,j,knew)=0.                              ! (no-flux, default)
          enddo                                   
        endif                                              ! <-- end of boundary conditions
      endif                                                ! <-- end of Western edge


      if (inode==np_xi-1.and.(.not.ew_periodic)) then      ! Eastern edge    
        if (obc_east.and.obc_m2flather) then               ! open boundary conditions
                                                           ! (m2flather)
          do j=1,ny
            cff=0.5*(h(nx,j)+h(nx+1,j))
            hx=sqrt(g/cff)
            cx=dtfast*cff*hx*0.5*(pm(nx,j)+pm(nx+1,j))
            zx=(0.5+cx)*zeta(nx,j,kstp)+(0.5-cx)*zeta(nx+1,j,kstp)
            if (cx > 0.292893218813452) then
              zx=zx + ( zeta(nx,j,knew) +cx*zeta(nx+1,j,kstp)
     &                               -(1.+cx)*zeta(nx,j,kstp)
     &                            )*(1.-0.292893218813452/cx)**2
            endif
            ubar(nx+1,j,knew)=0.5*( (1.-cx)*ubar(nx+1,j,kstp)
     &                                    +cx*ubar(nx,j,kstp)
     &                   +ubar_east(j) +hx*(zx-zeta_east(j)))
     &                                         *umask(nx+1,j)
          enddo

        elseif (obc_east.and.(.not.obc_m2flather)) then    ! zero gradient OBC
          do j=1,ny
            ubar(ny+1,j,knew)=ubar(ny,j,knew)
          enddo

        else                                               ! closed boundary conditions 
          do j=1,ny
            ubar(nx+1,j,knew)=0.
          enddo
        endif                                              ! <-- end of boundary conditions
      endif                                                ! <-- end of Eastern edge



      if (jnode==0.and.(.not.ns_periodic)) then            ! Southern edge
        if (obc_south.and.obc_m2flather) then              ! open boundary conditions 
          ! NOTE m2orlanski is used instead of m2flather for tangential bounds
          do i=0,nx
            grad(i,0)=ubar(i+1,0,kstp)-ubar(i,0,kstp)
            grad(i,1)=ubar(i+1,1,kstp)-ubar(i,1,kstp)
          enddo
          do i=1,nx
            cx=-0.125*dtfast*(vbar(i,1,kstp)+vbar(i-1,1,kstp))
     &      *(pn(i,0)+pn(i-1,0)+pn(i,1)+pn(i-1,1))

            cy= 0.125*dtfast*(ubar(i,0,kstp)+ubar(i,1,kstp))
     &      *(pm(i,0)+pm(i-1,0)+pm(i,1)+pm(i-1,1))
            if (cx > 0.) then
              cext=0.
            else
              cext=-cx ; cx=0.
            endif
            ubar(i,0,knew)=(1.-cx)*( ubar(i,0,kstp)
     &                      -max(cy,0.)*grad(i-1,0)
     &                      -min(cy,0.)*grad(i  ,0)
     &                                            )
     &                  +cx*(        ubar(i,1,kstp)
     &                      -max(cy,0.)*grad(i-1,1)
     &                      -min(cy,0.)*grad(i  ,1)
     &                                            )
            ubar(i,0,knew)=(1.-cext)*ubar(i,0,knew)
     &                                  +cext*ubar_south(i)
            ubar(i,0,knew)=ubar(i,0,knew)*umask(i,0)
          enddo

        elseif (obc_south.and.(.not.obc_m2flather)) then   ! other OBC
          do i=1,nx
            ubar(i,0,knew)=ubar(i,1,knew)*umask(i,0)
          enddo

        else                                               ! closed boundary conditions
          do i=1,i1
            ubar(i,0,knew)=gamma2*ubar(i,1,knew)
     &                               *umask(i,0)
          enddo
        endif                                              ! <-- end of boundary conditions
      endif                                                ! <-- end of Southern egde

      if (jnode==np_eta-1.and.(.not.ns_periodic)) then     ! Northern edge
        if (obc_north.and.obc_m2flather) then              ! open boundary conditions 
          ! NOTE m2orlanski is used instead of m2flather for tangential bounds
          do i=0,nx
            grad(i,ny)=ubar(i+1,ny,kstp)-ubar(i,ny,kstp)
            grad(i,ny+1)=ubar(i+1,ny+1,kstp)-ubar(i,ny+1,kstp)
          enddo
          do i=1,nx
            cx=0.125*dtfast*(vbar(i,ny+1,kstp)+vbar(i-1,ny+1,kstp))
     &         *(pn(i,ny)+pn(i-1,ny)+pn(i,ny+1)+pn(i-1,ny+1))

            cy=0.125*dtfast*(ubar(i,ny,kstp)+ubar(i,ny+1,kstp))
     &         *(pm(i,ny)+pm(i-1,ny)+pm(i,ny+1)+pm(i-1,ny+1))
            if (cx > 0.) then
              cext=0.
            else
              cext=-cx ; cx=0.
            endif
            ubar(i,ny+1,knew)=(1.-cx)*( ubar(i,ny+1,kstp)
     &                         -max(cy,0.)*grad(i-1,ny+1)
     &                         -min(cy,0.)*grad(i  ,ny+1)
     &                                                  )
     &                       +cx*(        ubar(i,ny,kstp)
     &                           -max(cy,0.)*grad(i-1,ny)
     &                           -min(cy,0.)*grad(i  ,ny)
     &                                                  )
            ubar(i,ny+1,knew)=(1.-cext)*ubar(i,ny+1,knew)
     &                                  +cext*ubar_north(i)
            ubar(i,ny+1,knew)=ubar(i,ny+1,knew)*umask(i,ny+1)
          enddo

        elseif (obc_north.and.(.not.obc_m2flather)) then   ! other OBC
          do i=1,nx
            ubar(i,ny+1,knew)=ubar(i,ny+1,knew)*umask(i,ny+1)
          enddo

        else                                               ! closed boundary conditions 
          do i=1,i1
            ubar(i,ny+1,knew)=gamma2*ubar(i,ny,knew)
     &                                *umask(i,ny+1)
          enddo
        endif                                              ! <-- end of boundary conditions
      endif                                                ! <-- end of Northern edge


      ! Periodic boundaries
      ! ======== ==========

      if (ew_periodic) then
        if (jnode==0) then                                 ! Southern "edge"
          do i=1,nx
            ubar(i,0,knew)=gamma2*ubar(i,1,knew)
     &                               *umask(i,0)
          enddo          
        elseif (jnode==np_eta-1) then                      ! Northern "edge"
          do i=1,nx                                        
            ubar(i,ny+1,knew)=gamma2*ubar(i,ny,knew)
     &                                *umask(i,ny+1)
          enddo
        endif
      endif

      ! Corners between adjacent open boundaries
      ! ======= ======= ======== ==== ==========

      ! South West
      if (obc_south.and.obc_west) then
        if (inode==0.and.jnode==0) then
          ubar(1,0,knew)=0.5*( ubar(2,0,knew)
     &                      +ubar(1,1,knew) )
        endif
      endif

      ! South East
      if (obc_south.and.obc_east) then
        if (inode==np_xi-1.and.jnode==0) then
          ubar(nx+1,0,knew)=0.5*( ubar(nx,0,knew)
     &                       +ubar(nx+1,1,knew) )
        endif
      endif

      ! North West
      if (obc_north.and.obc_west) then
        if (inode==0.and.jnode==np_eta-1) then
          ubar(1,ny+1,knew)=0.5*( ubar(2,ny+1,knew)
     &                           +ubar(1,ny,knew) )
        endif
      endif

      ! North East
      if (obc_north.and.obc_east) then
        if (inode==np_xi-1.and.jnode==np_eta-1) then
          ubar(nx+1,ny+1,knew)=0.5*( ubar(nx,ny+1,knew)
     &                            +ubar(nx+1,ny,knew) )
        endif
      endif

! Daniel add
!====================================================================
!       NO LAGRANGIAN FLUX CONDITION AT LAND-SEA INTERFACES
!====================================================================
#if defined WEC && defined MASKING
      do i=1,j1
        do j=j0,j1
          ubar(i,j,knew)=ubar(i,j,knew)*umask(i,j)
     &                  +ust2d(i,j)*(umask(i,j)-1.0)
        enddo
      enddo
#endif

      end

