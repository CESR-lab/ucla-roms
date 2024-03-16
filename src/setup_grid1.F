#include "cppdefs.opt"
                                            ! Setting up curvilinear
      subroutine setup_grid1(tile)          ! grid: Compute various
      use param
      use hidden_mpi_vars
      implicit none                         ! combinations of metric
      integer tile                          ! terms.
#include "compute_tile_bounds.h"
      call setup_Corls_tile(istr,iend,jstr,jend)
      call setup_grid1_tile(istr,iend,jstr,jend)
      end

      subroutine setup_Corls_tile(istr,iend,jstr,jend)
      use param
      use grid
      use scalars
      implicit none
      integer istr,iend,jstr,jend, i,j

#ifdef SPHERICAL
      real cff1, cff2
#endif

#include "compute_extended_bounds.h"

! Set f and f/mn at horizontal RHO-points.  Note that Coriolis arrays
! "f" and "fomn" must computed over the extended ranges of horizontal
! indices and must NOT be exchanged/communicated the standard way
! because of potential conflict with periodicity (if selected); same
! is true for the coordinate variables).

#ifdef SPHERICAL
! Coriolis frequency at North Pole is 4*pi/T where T is sidereal day

      cff1=4.*pi/(23.9344699D0*3600.D0)
      do j=jstrR,jendR
        do i=istrR,iendR
          f(i,j)=cff1*sin(deg2rad*latr(i,j))
          fomn(i,j)=f(i,j)/(pm(i,j)*pn(i,j))
# ifdef NON_TRADITIONAL
          cff2=cff1*cos(deg2rad*latr(i,j))
          f_XI(i,j)= cff2*sin(angler(i,j))
          f_ETA(i,j)=cff2*cos(angler(i,j))
# endif
        enddo
      enddo
#else
      do j=jstrR,jendR
        do i=istrR,iendR
          fomn(i,j)=f(i,j)/(pm(i,j)*pn(i,j))
        enddo
      enddo
#endif
      end


      subroutine setup_grid1_tile (istr,iend,jstr,jend)
      use param
      use grid
      use scalars
      implicit none
      integer :: istr,iend,jstr,jend, i,j

#ifdef MASKING
      real cff1, cff2
#endif

#include "compute_auxiliary_bounds.h"

!  Compute n/m and m/n; all at horizontal RHO-points.

      do j=jstrR,jendR
        do i=istrR,iendR
          dm_r(i,j)=1./pm(i,j)
          dn_r(i,j)=1./pn(i,j)
        enddo
      enddo

#if (defined CURVGRID && defined UV_ADV)

!  Compute d(1/n)/d(xi) and d(1/m)/d(eta) tems, both at RHO-points.

      do j=jstrR,jendR
        do i=istr,iend
          dndx(i,j)=0.5/pn(i+1,j)-0.5/pn(i-1,j)
        enddo
      enddo
      do j=jstr,jend
        do i=istrR,iendR
          dmde(i,j)=0.5/pm(i,j+1)-0.5/pm(i,j-1)
        enddo
      enddo

#endif

!  Compute m/n at horizontal U-points.

      do j=jstrR,jendR
        do i=istr,iendR
           pmon_u(i,j)=(pm(i,j)+pm(i-1,j))/(pn(i,j)+pn(i-1,j))
           dm_u(i,j)=2.0/(pm(i,j)+pm(i-1,j))
           dn_u(i,j)=2.0/(pn(i,j)+pn(i-1,j))
           pn_u(i,j)=0.5*(pn(i,j)+pn(i-1,j))
#ifdef MASKING
           umask(i,j)=rmask(i,j)*rmask(i-1,j)
#endif
           iA_u(i,j)=0.25*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
        enddo
      enddo

!  Compute n/m at horizontal V-points.

      do j=jstr,jendR
        do i=istrR,iendR
          pnom_v(i,j)=(pn(i,j)+pn(i,j-1))/(pm(i,j)+pm(i,j-1))
          dm_v(i,j)=2.0/(pm(i,j)+pm(i,j-1))
          dn_v(i,j)=2.0/(pn(i,j)+pn(i,j-1))
          pm_v(i,j)=0.5*(pm(i,j)+pm(i,j-1))

#ifdef MASKING
          vmask(i,j)=rmask(i,j)*rmask(i,j-1)
#endif
          iA_v(i,j)=0.25*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
        enddo
      enddo

! Compute n/m and m/n at horizontal PSI-points.
! Set mask according to slipperness parameter gamma.

      do j=jstr,jendR
        do i=istr,iendR
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
