#include "cppdefs.h"
#ifdef AVERAGES

      subroutine set_avg (tile)

#ifdef WEC
      use wec_frc, only: wstokes_wec_tile
#endif

      implicit none
      integer tile
# include "param.h"
# ifdef SOLVE3D
#  include "ncvars.h"
#  include "work.h"
#  include "private_scratch.h"
# endif
# include "compute_tile_bounds.h"
# ifdef SOLVE3D
      if (wrtavg(indxW)) then
        call wvlcty_tile(istr,iend,jstr,jend, work, A2d(1,1),
     &                                    A2d(1,2), A2d(1,3))
      endif

! Daniel add
#  ifdef WEC
      if (wrtavg(indxWST)) then
        ! ** wstokes_tile needs to be added to wec module **
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

! DevinD add
# ifdef WEC
      use wec_frc
# endif

      implicit none
      integer istr,iend,jstr,jend, i,j
# ifdef SOLVE3D
     &                       , itrc, k
# endif
      real cff,cff1
# include "param.h"
# include "scalars.h"
# include "ncvars.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "eos_vars.h"
# include "mixing.h"
# include "averages.h"
# ifdef SOLVE3D
#  include "work.h"
# endif

# include "compute_auxiliary_bounds.h"

      if (iic > ntsavg) then
        cff=1./dble(navg)    !<-- scaling
        if (mod(iic-ntsavg,navg) == 1) then
          if (ZEROTH_TILE) then
            mpi_master_only write(*,'(7x,A,I8,2x,A,I8,2x,A,I5)')
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
        do itrc=1,NT
          if (wrtavg(indxT+itrc-1)) then
            do k=1,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  t_avg(i,j,k,itrc)=cff1*t_avg(i,j,k,itrc)
     &                              +cff*t(i,j,k,nstp,itrc)
                enddo
              enddo
            enddo
          endif
        enddo
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
# endif /* SOLVE3D */
! Daniel add
# ifdef WEC 
        if (wrtavg(indxSUP)) then
          do j=jstrR,jendR
            do i=istrR,iendR
              sup_avg(i,j)= cff1*sup_avg(i,j) +cff*sup(i,j)
            enddo
          enddo
        endif
        if (wrtavg(indxUST2D)) then
          do j=jstrR,jendR
            do i=istrR,iendR
              ust2d_avg(i,j)= cff1*ust2d_avg(i,j) +cff*ust2d(i,j)
            enddo
          enddo
        endif
        if (wrtavg(indxVST2D)) then
          do j=jstrR,jendR
            do i=istrR,iendR
              vst2d_avg(i,j)= cff1*vst2d_avg(i,j) +cff*vst2d(i,j)
            enddo
          enddo
        endif

#  ifdef SOLVE3D
        if (wrtavg(indxUST)) then
          do k=1,N
            do j=jstrR,jendR
              do i=istrR,iendR
                ust_avg(i,j,k)= cff1*ust_avg(i,j,k) +cff*ust(i,j,k)
              enddo
            enddo
          enddo
        endif
        if (wrtavg(indxVST)) then
          do k=1,N
            do j=jstrR,jendR
              do i=istrR,iendR
                vst_avg(i,j,k)= cff1*vst_avg(i,j,k) +cff*vst(i,j,k)
              enddo
            enddo
          enddo
        endif
        if (wrtavg(indxWST)) then
          do k=1,N
            do j=jstrR,jendR
              do i=istrR,iendR
                wst_avg(i,j,k)= cff1*wst_avg(i,j,k) +cff*wst(i,j,k)
              enddo
            enddo
          enddo
        endif
        if (wrtavg(indxAkb)) then
          do k=0,N
            do j=jstrR,jendR
              do i=istrR,iendR
                akb_avg(i,j,k)= cff1*akb_avg(i,j,k) +cff*Akb(i,j,k)
              enddo
            enddo
          enddo
        endif
        if (wrtavg(indxAkw)) then
          do k=0,N
            do j=jstrR,jendR
              do i=istrR,iendR
                akw_avg(i,j,k)= cff1*akw_avg(i,j,k) +cff*Akw(i,j,k)
              enddo
            enddo
          enddo
        endif
        if (wrtavg(indxKVF)) then
          do k=1,N
            do j=jstrR,jendR
              do i=istrR,iendR
                kvf_avg(i,j,k)= cff1*kvf_avg(i,j,k) +cff*kvf(i,j,k)
              enddo
            enddo
          enddo
        endif
        if (wrtavg(indxCALP)) then
          do j=jstrR,jendR
            do i=istrR,iendR
              calp_avg(i,j)= cff1*calp_avg(i,j) +cff*calP(i,j)
            enddo
          enddo
        endif
        if (wrtavg(indxKAPS)) then
          do j=jstrR,jendR
            do i=istrR,iendR
              kaps_avg(i,j)= cff1*kaps_avg(i,j) +cff*Kapsrf(i,j)
            enddo
          enddo
        endif
#  endif  /* SOLVE3D */
# endif   /* WEC */



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
