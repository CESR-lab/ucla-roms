      subroutine zetabc (zeta_new)

      ! Set boundary conditions for free surface zeta_new

      use param
      use boundary
      use grid
      use ocean_vars
      use scalars

      implicit none
      
      ! input/output
      real,dimension(1-bf:nx +bf,1-bf:ny +bf),intent(inout) :: zeta_new

      ! local
      real :: cx, cff,dft,dfx
      integer i,j
      real, parameter :: eps=1.E-33

      if (inode==0.and.(.not.ew_periodic)) then
        if (obc_west) then
          do j=0,ny                                        ! West Flather obc
            cx=dtfast*pm(1,j)*sqrt(g*h(1,j))
            zeta_new(0,j)=( (1.-cx)*zeta(0,j,kstp)
     &                        +cx*zeta(1,j,kstp) )
     &                                 *rmask(0,j)
          enddo
        else
          do j=0,ny                                        ! Western edge
            zeta_new(0,j)=zeta_new(1,j)                    ! zero gradient BC
     &                      *rmask(0,j)
          enddo
        endif
      endif

      if (inode==np_xi-1.and.(.not.ew_periodic)) then
        if (obc_east) then
          do j=0,ny                                        ! East Flather obc
            cx=dtfast*pm(nx,j)*sqrt(g*h(nx,j))
            zeta_new(nx+1,j)=( (1.-cx)*zeta(nx+1,j,kstp)
     &                               +cx*zeta(nx,j,kstp) )
     &                                      *rmask(nx+1,j)
          enddo
        else
          do j=0,ny                                        ! Eastern edge
            zeta_new(nx+1,j)=zeta_new(nx,j)                ! zero gradient BC
     &                       *rmask(nx+1,j)
          enddo
        endif
      endif


      if (jnode==0.and.(.not.ns_periodic)) then
        if (obc_south) then
          do i=0,nx                                        ! South Flather obc
            cx=dtfast*pn(i,1)*sqrt(g*h(i,1))
            zeta_new(i,0)=( (1.-cx)*zeta(i,0,kstp)
     &                        +cx*zeta(i,1,kstp) )
     &                                 *rmask(i,0)
          enddo
        else
          do i=0,nx                                        ! Southern edge
            zeta_new(i,0)=zeta_new(i,1)                    ! zero gradient BC
     &                      *rmask(i,0)
          enddo
        endif
      endif


      if (jnode==np_eta-1.and.(.not.ns_periodic)) then
        if (obc_north) then
          do i=0,nx                                        ! North Flather obc
            cx=dtfast*pn(i,ny)*sqrt(g*h(i,ny))
            zeta_new(i,ny+1)=( (1.-cx)*zeta(i,ny+1,kstp)
     &                             +cx*zeta(i,ny,kstp) )
     &                                    *rmask(i,ny+1)
          enddo
        else
          do i=0,nx                                        ! Northern edge
            zeta_new(i,ny+1)=zeta_new(i,ny)                ! zero gradient BC
     &                       *rmask(i,ny+1)
          enddo
        endif
      endif

      ! Corners between adjacent open boundaries
      if ((.not.ew_periodic).and.(.not.ns_periodic)) then
        if (jnode==0.and.inode==0) then
          zeta_new(0,0)=0.5*( zeta_new(1,0)
     &                       +zeta_new(0,1) )
        endif
        if (jnode==0.and.inode==np_xi-1) then
          zeta_new(nx+1,0)=0.5*( zeta_new(nx,0)
     &                        +zeta_new(nx+1,1) )
        endif
        if (jnode==np_eta-1.and.inode==0) then
          zeta_new(0,ny+1)=0.5*( zeta_new(1,ny+1)
     &                            +zeta_new(0,ny) )
        endif
        if (jnode==np_eta-1.and.inode==np_xi-1) then
          zeta_new(nx+1,ny+1)=0.5*( zeta_new(nx,ny+1)
     &                             +zeta_new(nx+1,ny) )
        endif
      endif

      end
