#include "cppdefs.opt"
                                            ! Setting up curvilinear
      subroutine setup_grid1                ! grid: Compute various
      use param
      use roms_mpi
      implicit none                         ! combinations of metric

      call setup_Corls_tile
      call setup_grid1_tile

      end

      subroutine setup_Corls_tile
      use param
      use grid
      use scalars

      implicit none

      integer :: i,j

#ifdef SPHERICAL
      logical,parameter :: latlon  = .true.
#else
      logical,parameter :: latlon  = .false.
#endif
      logical,parameter :: nontrad = .false.

      real :: cff1, cff2

! Using extended bounds amounts to filling 
! the full buffers of the array. Strictly 
! speaking, we can't do that at the physical
! boundaries since there is no information beyond
! a single buffer point. Ignoring that fact,
! the loop ranges become much simpler at the costs
! of doing wasted computation at the physical boundaries
      ! Filling the full buffers negates the need for mpi exchanges afterwards
      ! We can extend the computation for the full buffer, which would negate
      ! the need for an mpi_exchange afterwards

      ! Coriolis frequency at North Pole is 4*pi/T where T is sidereal day
      cff1=4.*pi/(23.9344699D0*3600.D0) ! We're pedantic about day length

      ! Also, we're recomputing f and ignoring f from the grid file
      ! If we want to avoid mpi_exchanges we should explicitly add in the 
      ! buffer size bf
      do j= 1-bf,ny+bf
        do i= 1-bf,nx+bf
          if (latlon) then ! Grid has lat/lon coordinates
            f(i,j)=cff1*sin(deg2rad*latr(i,j))
            if (nontrad) then ! include horizontal components of coriolis
              cff2=cff1*cos(deg2rad*latr(i,j))
              f_XI(i,j)= cff2*sin(angler(i,j))
              f_ETA(i,j)=cff2*cos(angler(i,j))
            endif
          endif
          fomn(i,j)=f(i,j)/(pm(i,j)*pn(i,j))
        enddo
      enddo
 
      end subroutine setup_Corls_tile


      subroutine setup_grid1_tile
      use param
      use grid
      use scalars
      use roms_mpi
      implicit none
      integer :: i,j

#ifdef MASKING
      real cff1, cff2
#endif

      ! Rely on mpi exchange to fill buffer values
      ! need to compute the physical boundary buffers
      ! since we can't get those from mpi
      ! using i0,i1 maps onto istrR,iendR in this case.

      !  Compute dx (1/pm) and dy (1/pn) at rho-points.
      do j=j0,j1
        do i=i0,i1
          dm_r(i,j)=1./pm(i,j)
          dn_r(i,j)=1./pn(i,j)
        enddo
      enddo

#if (defined CURVGRID && defined UV_ADV)

!  Compute d(1/n)/d(xi) and d(1/m)/d(eta) tems, both at RHO-points.

      do j=j0,j1
        do i=1,nx
          dndx(i,j)=0.5/pn(i+1,j)-0.5/pn(i-1,j)
        enddo
      enddo
      do j=1,ny
        do i=i0,i1
          dmde(i,j)=0.5/pm(i,j+1)-0.5/pm(i,j-1)
        enddo
      enddo

#endif

!  Compute m/n at horizontal U-points.

      do j=j0,j1
        do i=i0,i1
           pmon_u(i,j)=(pm(i,j)+pm(i-1,j))/(pn(i,j)+pn(i-1,j))
           dm_u(i,j)=2.0/(pm(i,j)+pm(i-1,j))
           dn_u(i,j)=2.0/(pn(i,j)+pn(i-1,j))
           pn_u(i,j)=0.5*(pn(i,j)+pn(i-1,j))
           umask(i,j)=rmask(i,j)*rmask(i-1,j)
           iA_u(i,j)=0.25*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
        enddo
      enddo

!  Compute n/m at horizontal V-points.

      do j=j0,j1
        do i=i0,i1
          pnom_v(i,j)=(pn(i,j)+pn(i,j-1))/(pm(i,j)+pm(i,j-1))
          dm_v(i,j)=2.0/(pm(i,j)+pm(i,j-1))
          dn_v(i,j)=2.0/(pn(i,j)+pn(i,j-1))
          pm_v(i,j)=0.5*(pm(i,j)+pm(i,j-1))
          vmask(i,j)=rmask(i,j)*rmask(i,j-1)
          iA_v(i,j)=0.25*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
        enddo
      enddo

! Compute n/m and m/n at horizontal PSI-points.
! Set mask according to slipperness parameter gamma.

      do j=1,j1
        do i=1,i1
          dm_p(i,j)=4./(pm(i,j)+pm(i,j-1)+pm(i-1,j)+pm(i-1,j-1))
          dn_p(i,j)=4./(pn(i,j)+pn(i,j-1)+pn(i-1,j)+pn(i-1,j-1))

#ifdef MASKING
c          pmask(i,j)=rmask(i,j)*rmask(i-1,j)*rmask(i,j-1)
c     &                                      *rmask(i-1,j-1)
C***      if (gamma2 < 0.) pmask(i,j)=2.-pmask(i,j)

! Set no-slip boundary conditions on land-mask boundaries
! regardless of supplied value of gamma2.


          cff1=1.  !<-- computation of off-diagonal nonlinear terms
          cff2=2.

          if (rmask(i-1,j  )>0.5 .and. rmask(i,j  )>0.5 .and.
     &        rmask(i-1,j-1)>0.5 .and. rmask(i,j-1)>0.5) then
            pmask(i,j)=1.

          elseif(rmask(i-1,j  )<0.5 .and. rmask(i,j  )>0.5 .and.
     &           rmask(i-1,j-1)>0.5 .and. rmask(i,j-1)>0.5) then
            pmask(i,j)=cff1
          elseif(rmask(i-1,j  )>0.5 .and. rmask(i,j  )<0.5 .and.
     &           rmask(i-1,j-1)>0.5 .and. rmask(i,j-1)>0.5) then
            pmask(i,j)=cff1
          elseif(rmask(i-1,j  )>0.5 .and. rmask(i,j  )>0.5 .and.
     &           rmask(i-1,j-1)<0.5 .and. rmask(i,j-1)>0.5) then
            pmask(i,j)=cff1
          elseif(rmask(i-1,j  )>0.5 .and. rmask(i,j  )>0.5 .and.
     &           rmask(i-1,j-1)>0.5 .and. rmask(i,j-1)<0.5) then
            pmask(i,j)=cff1


          elseif(rmask(i-1,j  )>0.5 .and. rmask(i,j  )<0.5 .and.
     &           rmask(i-1,j-1)>0.5 .and. rmask(i,j-1)<0.5) then
            pmask(i,j)=cff2
          elseif(rmask(i-1,j  )<0.5 .and. rmask(i,j  )>0.5 .and.
     &           rmask(i-1,j-1)<0.5 .and. rmask(i,j-1)>0.5) then
            pmask(i,j)=cff2
          elseif(rmask(i-1,j  )>0.5 .and. rmask(i,j  )>0.5 .and.
     &           rmask(i-1,j-1)<0.5 .and. rmask(i,j-1)<0.5) then
            pmask(i,j)=cff2
          elseif(rmask(i-1,j  )<0.5 .and. rmask(i,j  )<0.5 .and.
     &           rmask(i-1,j-1)>0.5 .and. rmask(i,j-1)>0.5) then
            pmask(i,j)=cff2

          else
            pmask(i,j)=0.
          endif
#endif
        enddo
      enddo

#ifdef EXCHANGE
      call exchange_xxx(dm_r,dn_r,dm_p)
      call exchange_xxx(dn_p)
      call exchange_xxx(dm_u, dn_u, iA_u)
      call exchange_xxx(dm_v, dn_v, iA_v)
      call exchange_xxx(pmon_u, pnom_v)

# if defined CURVGRID && defined UV_ADV
      call exchange_xxx(dndx, dmde)
# endif

# ifdef MASKING
      call exchange_xxx(rmask,umask)
      call exchange_xxx(vmask,pmask)
# endif
#endif /* EXCHANGE */
      end
