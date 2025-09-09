#include "cppdefs.opt"

      subroutine v2dbc_tile (grad)

! Set lateral boundary conditions for the barotropic (i.e. vertically
! integrated) ETA-component velocity vbar(:,:,knew).

      use param
#ifdef WEC
      use wec_frc
#endif
      use boundary
      use grid
      use ocean_vars
      use scalars

      implicit none

      ! in 
      real, dimension(1-bf:nx +bf,1-bf:ny +bf) :: grad

      ! local
      integer i,j
      real :: cx,cy,cext,cff,dft,dfx,dfy,zx,hx
      real, parameter :: eps=1.D-33
     
      if (jnode==0.and.(.not.ns_periodic)) then            ! Southern edge 
        if (obc_south.and.obc_m2flather) then              ! open boundary conditions 
                                                           ! (m2flather)
          do i=1,nx
            cff=0.5*(h(i,0)+h(i,1))
            hx=sqrt(g/cff)
            cx=dtfast*cff*hx*0.5*(pn(i,0)+pn(i,1))
            zx=(0.5+cx)*zeta(i,1,kstp)+(0.5-cx)*zeta(i,0,kstp)
            if (cx > 0.292893218813452) then
              zx=zx + ( zeta(i,1,knew) +cx*zeta(i,0,kstp)
     &                            -(1.+cx)*zeta(i,1,kstp)
     &                        )*(1.-0.292893218813452/cx)**2
            endif
            vbar(i,1,knew)=0.5*( (1.-cx)*vbar(i,1,kstp)
     &                               +cx*vbar(i,2,kstp)
     &           +vbar_south(i) -hx*(zx-zeta_south(i)))
     &                                      *vmask(i,1)
          enddo

        elseif (obc_south.and.(.not.obc_m2flather)) then   ! other OBC
          do i=1,nx
            vbar(i,1,knew)=vbar(i,2,knew)*vmask(i,1)
          enddo

        else                                               ! closed boundary conditions
          do i=1,nx
#  ifdef WEC
            vbar(i,1,knew)=-vst2d(i,1)                     ! no Lagrangian flux
     &                     *vmask(i,1)
#  else
            vbar(i,1,knew)=0.                              ! (no-flux, default)
#  endif 
          enddo                                
        endif                                              ! <-- end of boundary conditions
      endif                                                ! <-- end of Southern edge



      if (jnode==np_eta-1.and.(.not.ns_periodic)) then     ! Northern edge  
        if (obc_north.and.obc_m2flather) then              ! open boundary conditions
                                                           ! (m2flather)
          do i=1,nx
            cff=0.5*(h(i,ny)+h(i,ny+1))
            hx=sqrt(g/cff)
            cx=dtfast*cff*hx*0.5*(pn(i,ny)+pn(i,ny+1))
            zx=(0.5+cx)*zeta(i,ny,kstp)+(0.5-cx)*zeta(i,ny+1,kstp)
            if (cx > 0.292893218813452) then
              zx=zx + ( zeta(i,ny,knew) +cx*zeta(i,ny+1,kstp)
     &                               -(1.+cx)*zeta(i,ny,kstp)
     &                            )*(1.-0.292893218813452/cx)**2
            endif
            vbar(i,ny+1,knew)=0.5*( (1.-cx)*vbar(i,ny+1,kstp)
     &                                    +cx*vbar(i,ny,kstp)
     &                 +vbar_north(i) +hx*(zx-zeta_north(i)))
     &                                         *vmask(i,ny+1)
          enddo

        elseif (obc_south.and.(.not.obc_m2flather)) then   ! other OBC
          do i=1,nx
            vbar(i,ny+1,knew)=vbar(i,ny,knew)*vmask(i,ny+1)
          enddo

        else                                               ! closed boundary conditions
          do i=1,nx
#  ifdef WEC
            vbar(i,ny+1,knew)=-vst2d(i,ny+1)               ! no Lagrangian flux
     &                        *vmask(i,ny+1)
#  else
            vbar(i,ny+1,knew)=0.                         ! (no-flux, default)
#  endif                        
          enddo                               
        endif                                              ! <-- end of boundary conditions 
      endif                                                ! <-- end of Northern edge


 
      if (inode==0.and.(.not.ew_periodic)) then            ! Western edge    
        if (obc_west.and.obc_m2flather) then               ! open boundary condition
          ! NOTE m2orlanski is used instead of m2flather for tangential bounds
          do j=0,ny
            grad(0,j)=vbar(0,j+1,kstp)-vbar(0,j,kstp)
            grad(1,j)=vbar(1,j+1,kstp)-vbar(1,j,kstp)
          enddo
          do j=1,ny
            cx=-0.125*dtfast*(ubar(1,j,kstp)+ubar(1,j-1,kstp))
     &      *(pm(0,j)+pm(0,j-1)+pm(1,j)+pm(1,j-1))
            cy= 0.125*dtfast*(vbar(0,j,kstp)+vbar(1,j,kstp))
     &      *(pn(0,j)+pn(0,j-1)+pn(1,j)+pn(1,j-1))
            if (cx > 0.) then
              cext=0.
            else
              cext=-cx ; cx=0.
            endif
            vbar(0,j,knew)=(1.-cx)*( vbar(0,j,kstp)
     &                      -max(cy,0.)*grad(0,j-1)
     &                      -min(cy,0.)*grad(0,j  )
     &                                            )
     &                  +cx*(        vbar(1,j,kstp)
     &                      -max(cy,0.)*grad(1,j-1)
     &                      -min(cy,0.)*grad(1,j  )
     &                                            )
            vbar(0,j,knew)=(1.-cext)*vbar(0,j,knew)
     &                           +cext*vbar_west(j)
            vbar(0,j,knew)=vbar(0,j,knew)*vmask(0,j)
          enddo

        elseif (obc_west.and.(.not.obc_m2flather)) then   ! other OBC
          do j=1,ny
            vbar(0,j,knew)=vbar(0,j,knew)*vmask(0,j)
          enddo

        else                                               ! closed boundary conditions   
          do j=1,j1
            vbar(0,j,knew)=gamma2*vbar(1,j,knew)
     &                               *vmask(0,j)
          enddo
        endif                                              ! <-- end of boundary conditions
      endif                                                ! <-- end of Western edge



      if (inode==np_xi-1.and.(.not.ew_periodic)) then      ! Eastern edge
        if (obc_east.and.obc_m2flather) then               ! open boundary conditions
          ! NOTE m2orlanski is used instead of m2flather for tangential bounds
          do j=0,ny
            grad(nx,j)=vbar(nx,j+1,kstp)-vbar(nx,j,kstp)
            grad(nx+1,j)=vbar(nx+1,j+1,kstp)-vbar(nx+1,j,kstp)
          enddo
          do j=1,ny
            cx=0.125*dtfast*(ubar(nx+1,j,kstp)+ubar(nx+1,j-1,kstp))
     &         *(pm(nx,j)+pm(nx,j-1)+pm(nx+1,j)+pm(nx+1,j-1))
            cy=0.125*dtfast*(vbar(nx,j,kstp)+vbar(nx+1,j,kstp))
     &         *(pn(nx,j)+pn(nx,j-1)+pn(nx+1,j)+pn(nx+1,j-1))
            if (cx > 0.) then
              cext=0.
            else
              cext=-cx ; cx=0.
            endif
            vbar(nx+1,j,knew)=(1.-cx)*( vbar(nx+1,j,kstp)
     &                         -max(cy,0.)*grad(nx+1,j-1)
     &                         -min(cy,0.)*grad(nx+1,j  )
     &                                                  )
     &                       +cx*(        vbar(nx,j,kstp)
     &                           -max(cy,0.)*grad(nx,j-1)
     &                           -min(cy,0.)*grad(nx,j  )
     &                                                  )
            vbar(nx+1,j,knew)=(1.-cext)*vbar(nx+1,j,knew)
     &                                 +cext*vbar_east(j)
            vbar(nx+1,j,knew)=vbar(nx+1,j,knew)*vmask(nx+1,j)
          enddo

        elseif (obc_east.and.(.not.obc_m2flather)) then    ! other OBC
          do j=1,ny
            vbar(nx+1,j,knew)=vbar(nx+1,j,knew)*vmask(nx+1,j)
          enddo

        else                                               ! closed boundary conditions
          do j=1,j1
            vbar(nx+1,j,knew)=gamma2*vbar(nx,j,knew)
     &                                *vmask(nx+1,j)
          enddo
        endif                                              ! <-- end of boundary conditions
      endif                                                ! <-- end of Eastern Edge


      ! Periodic boundaries
      ! ======== ==========

      if (ns_periodic) then                                ! periodic boundary conditions
        if (inode==0) then                                 ! Western "edge"
          do j=1,ny
            vbar(0,j,knew)=gamma2*vbar(1,j,knew)
     &                               *vmask(0,j)
          enddo
        elseif (inode==np_xi-1) then                       ! Eastern "edge"
          do j=1,ny
            vbar(nx+1,j,knew)=gamma2*vbar(nx,j,knew)
     &                                *vmask(nx+1,j)
          enddo
        endif
      endif


      ! Corners between adjacent open boundaries
      ! ======= ======= ======== ==== ==========

      ! South West
      if (obc_south.and.obc_west) then
        if (inode==0.and.jnode==0) then
          vbar(0,1,knew)=0.5*(vbar(0,2,knew)
     &                      +vbar(1,1,knew))
        endif
      endif

      ! South East
      if (obc_south.and.obc_east) then
        if (inode==np_xi-1.and.jnode==0) then
          vbar(nx+1,1,knew)=0.5*( vbar(nx+1,2,knew)
     &                            +vbar(nx,1,knew))
        endif
      endif

      ! North West
      if (obc_north.and.obc_west) then
        if (inode==0.and.jnode==np_eta-1) then
          vbar(0,ny+1,knew)=0.5*( vbar(0,ny,knew)
     &                        +vbar(1,ny+1,knew))
        endif
      endif

      ! North East
      if (obc_north.and.obc_east) then
        if (inode==np_xi-1.and.jnode==np_eta-1) then
          vbar(nx+1,ny+1,knew)=0.5*( vbar(nx+1,ny,knew)
     &                             +vbar(nx,ny+1,knew))
        endif
      endif

! Daniel add
!====================================================================
!       NO LAGRANGIAN FLUX CONDITION AT LAND-SEA INTERFACES
!====================================================================
#if defined WEC && defined MASKING
      do i=i0,i1
        do j=1,j1
          vbar(i,j,knew)=vbar(i,j,knew)*vmask(i,j)
     &                  +vst2d(i,j)*(vmask(i,j)-1.0)
        enddo
      enddo
#endif

      end

