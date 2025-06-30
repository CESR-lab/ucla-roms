#include "cppdefs.opt"
#if defined SOLVE3D && ( defined TS_DIF2 || defined TS_DIF4 )

      subroutine t3dmix
      use param
      use hidden_mpi_vars
      use private_scratch
      implicit none
      integer,save:: tile=0
      integer itrc

# include "compute_tile_bounds.h"
      do itrc=1,nt
        call t3dmix_tile (istr,iend,jstr,jend, itrc, A2d(1,1),
     &                                     A2d(1,2), A2d(1,3))
      enddo
      end

      subroutine t3dmix_tile (istr,iend,jstr,jend, itrc, FX,FE,LapT)

! Computes horizontal mixing of tracers along S-surfaces.

      use param
      use grid
      use mixing
      use ocean_vars
      use scalars
      use tracers
      use diagnostics

      implicit none
      integer istr,iend,jstr,jend, itrc, i,j,k
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FX,FE, LapT
# ifdef TS_DIF4
      integer imin,imax,jmin,jmax
# endif

# include "compute_auxiliary_bounds.h"

      do k=1,nz

# ifdef TS_DIF2

! Add in horizontal Laplacian diffusion along constant S-surfaces.
!--------------------------------------------------------------------
! Compute XI- and ETA-components of diffusive tracer flux.

        do j=jstr,jend
          do i=istr,iend+1
            FX(i,j)=0.25*(diff2(i,j,itrc)+diff2(i-1,j,itrc))
     &                   *pmon_u(i,j)*(Hz(i,j,k)+Hz(i-1,j,k))*(
     &                     t(i,j,k,nrhs,itrc)-t(i-1,j,k,nrhs,itrc)
     &                                                              )
#  ifdef MASKING
     &                                                   *umask(i,j)
#  endif
          enddo
        enddo
        do j=jstr,jend+1
          do i=istr,iend
            FE(i,j)=0.25*(diff2(i,j,itrc)+diff2(i,j-1,itrc))
     &                     *pnom_v(i,j)*(Hz(i,j,k)+Hz(i,j-1,k))*(
     &                      t(i,j,k,nrhs,itrc)-t(i,j-1,k,nrhs,itrc)
     &                                                              )
#  ifdef MASKING
     &                                                   *vmask(i,j)
#  endif
          enddo
        enddo
        if (diag_trc .and. wrt_t_dia(itrc) .and. calc_diag) then
          do j=jstr,jend
            do i=istr,iend+1
              Tdiag(i,j,k,tmixx,td(itrc)) =
     &           Tdiag(i,j,k,tmixx,td(itrc)) - FX(i,j)
            enddo
          enddo
          do j=jstr,jend+1
            do i=istr,iend
              Tdiag(i,j,k,tmixy,td(itrc)) =
     &           Tdiag(i,j,k,tmixy,td(itrc)) - FE(i,j)
            enddo
          enddo
        endif

! Add in horizontal diffusion of tracer [T m^3/s].

        do j=jstr,jend
          do i=istr,iend
            t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+dt*pm(i,j)*pn(i,j)
     &                         *(FX(i+1,j)-FX(i,j)+FE(i,j+1)-FE(i,j))
     &                                                     /Hz(i,j,k)
          enddo
        enddo
# endif /* TS_DIF2 */

# ifdef TS_DIF4
#  ifdef EW_PERIODIC
      imin=istr-1
      imax=iend+1
#  else
      if (WESTERN_EDGE) then
        imin=istr
      else
        imin=istr-1
      endif
      if (EASTERN_EDGE) then
        imax=iend
      else
        imax=iend+1
      endif
#  endif
#  ifdef NS_PERIODIC
      jmin=jstr-1
      jmax=jend+1
#  else
      if (SOUTHERN_EDGE) then
        jmin=jstr
      else
        jmin=jstr-1
      endif
      if (NORTHERN_EDGE) then
        jmax=jend
      else
        jmax=jend+1
      endif
#  endif


!  Compute/Add in horizontal biharmonic diffusion along constant
!  S-surfaces.
!-----------------------------------------------------------------
!  The biharmonic operator is computed by applying the Laplacian
!  operator twice.
!
!  Compute horizontal tracer flux in the XI-direction at U-points.

        do j=jmin,jmax
          do i=imin,imax+1
            FX(i,j)=pmon_u(i,j)*0.5*(Hz(i,j,k)+Hz(i-1,j,k))*(
     &                      t(i,j,k,nrhs,itrc)-t(i-1,j,k,nrhs,itrc)
     &                                                              )
#  ifdef MASKING
     &                                                   *umask(i,j)
#  endif
          enddo
        enddo

!  Compute horizontal tracer flux in the ETA-direction at V-points.

        do j=jmin,jmax+1
          do i=imin,imax
            FE(i,j)=pnom_v(i,j)*0.5*(Hz(i,j,k)+Hz(i,j-1,k))*(
     &                      t(i,j,k,nrhs,itrc)-t(i,j-1,k,nrhs,itrc)
     &                                                            )
#  ifdef MASKING
     &                                                 *vmask(i,j)
#  endif
          enddo
        enddo

!  Compute the first Laplacian, without mixing coefficient.
!  Multiply by the metrics of the second Laplacian.
!  Save into work array "LapT".

        do j=jmin,jmax
          do i=imin,imax
            LapT(i,j)=( FX(i+1,j)-FX(i,j) +FE(i,j+1)-FE(i,j) )
     &                               *pm(i,j)*pn(i,j)/Hz(i,j,k)
          enddo
        enddo

!  Apply boundary conditions (except periodic; closed or gradient)
!  to the first Laplacian.

#  ifndef EW_PERIODIC
        if (WESTERN_EDGE) then
          do j=jstr-1,jend+1
#   ifdef WESTERN_WALL
            LapT(istr-1,j)=0.
#   else
            LapT(istr-1,j)=LapT(istr,j)
#   endif
          enddo
        endif
        if (EASTERN_EDGE) then
          do j=jstr-1,jend+1
#   ifdef EASTERN_WALL
            LapT(iend+1,j)=0.
#   else
            LapT(iend+1,j)=LapT(iend,j)
#   endif
          enddo
        endif
# endif
# ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do i=istr-1,iend+1
#  ifdef SOUTHERN_WALL
            LapT(i,jstr-1)=0.
#  else
            LapT(i,jstr-1)=LapT(i,jstr)
#  endif
          enddo
        endif
        if (NORTHERN_EDGE) then
          do i=istr-1,iend+1
#  ifdef NORTHERN_WALL
            LapT(i,jend+1)=0.
#  else
            LapT(i,jend+1)=LapT(i,jend)
#  endif
          enddo
        endif
# endif /* !NS_PERIODIC */

! Compute FX=d(LapT)/d(xi) and FE=d(LapT)/d(eta) terms.  After that
! compute and add in biharmonic mixing [T m^3/s]. Multiply by mixing
! coefficient.

        do j=jstr,jend
          do i=istr,iend+1
            FX(i,j)=0.25*(diff4(i,j,itrc)+diff4(i-1,j,itrc))
     &                   *pmon_u(i,j)*(Hz(i,j,k)+Hz(i-1,j,k))
     &                               *(LapT(i,j)-LapT(i-1,j))
# ifdef MASKING
     &                                           *umask(i,j)
# endif
          enddo
        enddo
        do j=jstr,jend+1
          do i=istr,iend
            FE(i,j)=0.25*(diff4(i,j,itrc)+diff4(i,j-1,itrc))
     &                   *pnom_v(i,j)*(Hz(i,j,k)+Hz(i,j-1,k))
     &                               *(LapT(i,j)-LapT(i,j-1))
# ifdef MASKING
     &                                           *vmask(i,j)
# endif
          enddo
        enddo
        do j=jstr,jend
          do i=istr,iend
            t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+dt*pm(i,j)*pn(i,j)
     &                         *(FX(i,j)-FX(i+1,j)+FE(i,j)-FE(i,j+1))
     &                                                     /Hz(i,j,k)
          enddo
        enddo
#endif /* TS_DIF4 */
      enddo

# ifdef EXCHANGE
        call exchange_xxx(t(:,:,:,nnew,itrc) )
# endif
      end
#else
      subroutine t3dmix_empty
      end
#endif
