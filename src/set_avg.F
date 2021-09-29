#include "cppdefs.opt"
#ifdef AVERAGES

      subroutine set_avg (tile)

      use param
      use wec_frc          ! , only: wstokes_wec_tile, indxWST
      use hidden_mpi_vars
      use ncvars
      use private_scratch
      use work_mod

      implicit none
      integer tile
# include "compute_tile_bounds.h"
# ifdef SOLVE3D
      if (wrtavg(indxW)) then
        call wvlcty_tile(istr,iend,jstr,jend, work, A2d(1,1),
     &                                    A2d(1,2), A2d(1,3))
      endif

#  ifdef WEC
      if (wrtavg(indxWST)) then
        call wstokes_wec_tile(istr,iend,jstr,jend)
      endif
#  endif
# endif
      call set_avg_tile (istr,iend,jstr,jend)
      return
      end

      subroutine set_avg_tile(istr,iend,jstr,jend)

! Compute time-averaged fields within a tile.
! ------- ------------- ------ ------ - -----
! Because of the syncronization issues, a delayed-mode procedure is
! used for averaging.  This implies that all fields to be averaged are
! sampled during the next time step, rather than at the end of the time
! step when they are computed.
!
! Although this algorithm results in somewhat awkward controlling logic
! it has the advantage that that all fields to be sampled correspond to
! exactly the same time, which is time step "n".  Particularly, this is
! done this way because vertical velocity corresponding to the newly
! computed horizontal velocities becomes available only during the
! following time step. The same applies to the density field.
!
! The algorithm consists of three logical blocks: (1) initialization
! of the average arrays: when  mod(iic-1,navg)==1  the target arrays
! are set to the first contribution; (2) accumulation of averaged data,
! when mod(iic-1,navg)>1; and (3) adding the last contribution and
! scaling.

      use param
# ifdef WEC
      use wec_frc, only: set_avg_wec_tile
# endif

      use tracers, only: set_avg_trc
      use bgc ! set_avg_bgc_diag
      use surf_flux ! set_avg_surf_flux
      use averages
      use eos_vars
      use grid
      use mixing
      use ncvars
      use ocean2d
      use ocean3d
      use scalars
      use work_mod

      implicit none
      integer istr,iend,jstr,jend, i,j
# ifdef SOLVE3D
     &                       , itrc, k
# endif
      real cff,cff1

# include "compute_auxiliary_bounds.h"

      if (iic > ntsavg) then
        cff=1./dble(navg)    !<-- scaling
        if (mod(iic-ntsavg,navg) == 1) then
          if (ZEROTH_TILE) then
            mpi_master_only write(*,'(7x,A,I8,2x,A,I8,2x,A,I8)')
     &              'set_avg :: started averaging at iic =', iic,
     &                          'ntsavg=', ntsavg, 'navg=', navg
          endif
          cff1=0.   !--> initialize
        else
          cff1=1.   !--> keep adding
        endif

        if (ZEROTH_TILE) time_avg=cff1*time_avg + cff*time

        if (wrtavg(indxZ)) then
          do j=jstrR,jendR
            do i=istrR,iendR
              zeta_avg(i,j)=cff1*zeta_avg(i,j) + cff*zeta(i,j,knew)
            enddo
          enddo
        endif
        if (wrtavg(indxUb)) then
          do j=jstrR,jendR
            do i=istrR,iendR
              ubar_avg(i,j)=cff1*ubar_avg(i,j) + cff*ubar(i,j,knew)
            enddo
          enddo
        endif
        if (wrtavg(indxVb)) then
          do j=jstrR,jendR
            do i=istrR,iendR
              vbar_avg(i,j)=cff1*vbar_avg(i,j) +cff*vbar(i,j,knew)
            enddo
          enddo
        endif
# ifdef SOLVE3D
        if (wrtavg(indxU)) then
          do k=1,N
            do j=jstrR,jendR
              do i=istrR,iendR
                u_avg(i,j,k)=cff1*u_avg(i,j,k) +cff*u(i,j,k,nstp)
              enddo
            enddo
          enddo
        endif
        if (wrtavg(indxV)) then
          do k=1,N
            do j=jstrR,jendR
              do i=istrR,iendR
                v_avg(i,j,k)=cff1*v_avg(i,j,k) + cff*v(i,j,k,nstp)
              enddo
            enddo
          enddo
        endif
        ! Set tracer averages
        call set_avg_trc(istrR,iendR,jstrR,jendR,cff,cff1)

        if (wrtavg(indxR)) then
          do k=1,N
            do j=jstrR,jendR
              do i=istrR,iendR
#  ifdef SPLIT_EOS
                rho_avg(i,j,k)=cff1*rho_avg(i,j,k) +cff*( rho1(i,j,k)
     &                                        -qp1(i,j,k)*z_r(i,j,k) )
#  else
                rho_avg(i,j,k)=cff1*rho_avg(i,j,k) +cff*rho(i,j,k)
#  endif
              enddo
            enddo
          enddo
        endif
        if (wrtavg(indxO)) then
          if (mod(iic-ntsavg,navg) == 1) then
            do k=0,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  w_avg(i,j,k)=We(i,j,k)+Wi(i,j,k)
                enddo
              enddo
            enddo
          elseif (mod(iic-ntsavg,navg) > 1) then
            do k=0,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  w_avg(i,j,k)=w_avg(i,j,k) +We(i,j,k)+Wi(i,j,k)
                enddo
              enddo
            enddo
          elseif (mod(iic-ntsavg,navg) == 0) then
            do k=0,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  w_avg(i,j,k)=cff*( w_avg(i,j,k)+We(i,j,k)+Wi(i,j,k)
     &                                              )*pm(i,j)*pn(i,j)
                enddo
              enddo
            enddo
          endif
        endif
        if (wrtavg(indxW)) then
          do k=1,N
            do j=jstrR,jendR
              do i=istrR,iendR
                wvl_avg(i,j,k)=cff1*wvl_avg(i,j,k) +cff*work(i,j,k)
              enddo
            enddo
          enddo
        endif
        if (wrtavg(indxAkv)) then
          do k=0,N
            do j=jstrR,jendR
              do i=istrR,iendR
                akv_avg(i,j,k)=cff1*akv_avg(i,j,k) +cff*Akv(i,j,k)
              enddo
            enddo
          enddo
        endif
        if (wrtavg(indxAkt)) then
          do k=0,N
            do j=jstrR,jendR
              do i=istrR,iendR
                akt_avg(i,j,k)=cff1*akt_avg(i,j,k)+cff*Akt(i,j,k,itemp)
              enddo
            enddo
          enddo
        endif
#  ifdef SALINITY
        if (wrtavg(indxAks)) then
          do k=0,N
            do j=jstrR,jendR
              do i=istrR,iendR
                aks_avg(i,j,k)=cff1*aks_avg(i,j,k)+cff*Akt(i,j,k,isalt)
              enddo
            enddo
          enddo
        endif
#  endif
#  ifdef LMD_KPP
        if (wrtavg(indxHbls)) then
          do j=jstrR,jendR
            do i=istrR,iendR
              hbl_avg(i,j)=cff1*hbl_avg(i,j) +cff*hbls(i,j,nstp)
            enddo
          enddo
        endif
#  endif
#  ifdef LMD_BKPP
        if (wrtavg(indxHbbl)) then
          do j=jstrR,jendR
            do i=istrR,iendR
              hbbl_avg(i,j)=cff1*hbbl_avg(i,j) +cff*hbbl(i,j,nstp)
            enddo
          enddo
        endif
#  endif
#  if defined BIOLOGY_BEC2 && defined BEC2_DIAG
        call set_avg_bgc_diag(istrR,iendR,jstrR,jendR,cff,cff1)
#  endif
# endif /* SOLVE3D */

# ifdef WEC 
        call set_avg_wec_tile(istr,iend,jstr,jend,cff,cff1
     &                                 ,jstrR,jendR,istrR,iendR)
# endif /* WEC */

      ! Set surface flux averages of sustr, etc
#ifdef SURF_FLUX_OUTPUT_AVG
      call set_avg_surf_flux(istrR,iendR,jstrR,jendR,cff,cff1)
#endif

        if (mod(iic-ntsavg,navg) == 0) then
          if (ZEROTH_TILE) then
            mpi_master_only write(*,'(7x,A,I8,2x,A,I8,2x,A,I5)')
     &              'set_avg :: finished averaging at iic=', iic,
     &                         'ntsavg=', ntsavg,  'navg=', navg
          endif
        endif
      endif       !<-- iic > ntsavg
      end
#else
      subroutine set_avg_empty
      end
#endif /* AVERAGES */
