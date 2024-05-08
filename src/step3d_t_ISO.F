#include "cppdefs.opt"
#ifdef SOLVE3D

#ifndef ADV_ISONEUTRAL
# define UPSTREAM_TS
#endif
c--# define AKIMA
c--# define CONST_TRACERS

! The following two CPPs activate switching triad and stabilizing
! algorithms.  They must always be defined whenever the isoneutral
! hyperdiffusion is selected.


# ifdef ADV_ISONEUTRAL
#  define SW_TRIADS
#  define STABILIZE
# endif
                                  ! Note that arrays A2d(1,1:4) within
      subroutine step3d_t(tile)   ! the list of agruments in call below

      use param
      use hidden_mpi_vars
      use private_scratch

      implicit none               ! are repeated to use the same memory
      integer tile                ! for vertical and horizontal scratch
                                  ! arrays.

# include "compute_tile_bounds.h"

      call step3d_t_iso_tile( istr,iend,jstr,jend,
     &                       A2d(1,1), A2d(1,2), A2d(1,3), A2d(1,4),
     &                       A2d(1,1), A2d(1,2), A2d(1,3)
# ifdef ADV_ISONEUTRAL
     &                     , A2d(1,4), A2d(1,6), A2d(1,8), A2d(1,10)
     &                                                    , A3d(1,1)
#  ifdef STABILIZE
     &                                                    , A3d(1,2)
#  endif
# endif
     &                                                             )
      end

      subroutine step3d_t_iso_tile( istr,iend,jstr,jend, WC,FC,CF,DC
     &                                                 ,  FX,FE, wrk1
# ifdef ADV_ISONEUTRAL
     &                                           , FSC,dTdz,dTdx,dTde
     &                                                         , LapT
#  ifdef STABILIZE
     &                                                         ,  Akz
#  endif
# endif
     &                                                             )

      use param
      use pipe_frc
      use river_frc
      use surf_flux, only: stflx, srflx
      use tracers
#ifdef MARBL      
      use marbl_driver, only: marbldrv_column_physics
#endif      
      use diagnostics
      use eos_vars
      use grid
      use mixing
      use ocean_vars
      use scalars

      implicit none
      integer istr,iend,jstr,jend, imin,imax,jmin,jmax, i,j,k!, td
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N) :: WC,FC,CF,DC
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY)     :: FX,FE, wrk1
# ifdef ADV_ISONEUTRAL
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,2) :: FSC,  dTdz,
     &                                               dTdx, dTde
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,N) :: LapT
#  ifdef STABILIZE
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,0:N) :: Akz
#  endif
# endif
      real cff
      real, parameter :: epsil=1.E-33
      integer itrc, iAkt
# ifdef ADV_ISONEUTRAL
      integer k1,k2
      real sumX,sumE
#  ifdef SW_TRIADS
      integer idx,ide
      real, dimension(0:4), parameter :: wgt=(/ 0., 1., 0.5,
     &                           0.3333333333333333, 0.25 /)
#  endif
#  ifdef  STABILIZE
      real cff2, cffX, s2_XLU,s2_XRU, s2_ELU,s2_ERU,
     &           cffE, s2_XLL,s2_XRL, s2_ELL,s2_ERL
#  endif
# endif

# include "compute_auxiliary_bounds.h"

# ifndef EW_PERIODIC
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
# else
      imin=istr-1
      imax=iend+1
# endif
# ifndef NS_PERIODIC
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
# else
      jmin=jstr-1
      jmax=jend+1
# endif
#  define ntdf nstp
c--#  define ntdf nrhs

! Note that k-loop is recursive in the case of ADV_ISONEUTRAL switch
! active, while and "k1" and "k2" are alternating indices, 1-2-1-2-...
! etc, designed to mimic full 3D arrays, while using only two vertical
! slices. This applies to arrays  dTdx,dTde (both are at vertical RHO-
! points) and  FSC,dTdz (at vertical W-points). The index mapping
! rules for k1,k2 are as follows:
!   (i) k2 is always above k1;
!  (ii) arrays of W-type are computed ahead at k+1 --> k2
! (iii) arrays of RHO-type are computed at k -

!     if (mynode==1) then
!       write(60,*) 't:          ',t(10,10,20,nnew,1)
!     endif

      do itrc=1,nt

# ifdef ADV_ISONEUTRAL
        k2=1
        do k=0,N,+1       !--> recursive
          k1=k2
          k2=3-k1
          if (k > 0) then
# else
        do k=1,N
# endif

! Compute and add in contribution due to horizontal advection fluxes.
! Note that at this moment the content of array t(:,:,:,nnew,itrc) is
! concentration of tracer "itrc" taken at time corresponding to "nstp"
! and already MULTIPLIED by Hz(:,:,:) at "nstp" in "pre_step3d_tile".

# include "compute_horiz_tracer_fluxes.h"

# ifdef DIAGNOSTICS
#  ifndef ADV_ISONEUTRAL
          if (diag_trc .and. wrt_t_dia(itrc) .and. calc_diag) then
            call diag_t_adv_hc4(k,itrc)
          endif
#  endif
# endif

            do j=jstr,jend
              do i=istr,iend
                t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc) -dt*pm(i,j)
     &                                *pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                          +FE(i,j+1)-FE(i,j)
     &                                                           )
              enddo
            enddo                 !-->  discard FX,FE

# ifdef DIAGNOSTICS
#  ifdef ADV_ISONEUTRAL
            if (diag_trc .and. wrt_t_dia(itrc) .and. calc_diag) then
              do j= 1,ny
                do i=1,nx+1
                  Tdiag(i,j,k,tadvx,td(itrc)) = FX(i,j)
                enddo
              enddo
              do j= 1,ny+1
                do i=1,nx
                  Tdiag(i,j,k,tadvy,td(itrc)) = FE(i,j)
                enddo
              enddo
            endif
#  else
            if (diag_trc .and. wrt_t_dia(itrc) .and. calc_diag) then
              do j= 1,ny
                do i=1,nx+1
                  Tdiag(i,j,k,tmixx,td(itrc)) = FX(i,j)
     &             - Tdiag(i,j,k,tadvx,td(itrc))
                enddo
              enddo
              do j= 1,ny+1
                do i=1,nx
                  Tdiag(i,j,k,tmixy,td(itrc)) = FE(i,j)
     &             - Tdiag(i,j,k,tadvy,td(itrc))
                enddo
              enddo
            endif
#  endif
# endif

# ifdef ADV_ISONEUTRAL
          endif !<-- k > 0

! Computation of the rotated (isoneutral) biharmonic operator expects
! the following five fields be defined at entry:
!
!         d rho |(ad)           d rho |(ad)           [ d rho|(ad)]^-1
!  dRdx = ------|        dRde = ------|        idRz = [ -----|    ]
!          d xi | s             d eta | s             [  d z |    ]
!
!                    diff3u               diff3v
!
!  WARNING: it is assumed here that diff3u and diff3v are the square
! roots of the hyperdiffusivities
!
!
!                     /  |       Schematic of variable placement on
!                   /    |       the vertical grid and explanation of
!                 /      |       meaning of alternating indices "k1"
!               /        |       and "k2" in the code below:
!             /          |
!           /            |
!         /        dTdx(i+1,j,k2)       ------------------
!       /                |             |                  |
!      |                 |             |                  |
!      |                 |             |                  |
!      |                 |             |                  |
!      |                /|       dTdx(i,j,k2)     dTdx(i+1,j,k2)  k+1
!      |              /  |             |                  |
! dTdx(i,j,k2)      /    |             |                  |
!      |          /      |             |                  |
!      |   dTdz(i,j,k2)  |             |                  |
!      |      /          |              -- dTdz(i,j,k2) --        k+1/2
!      |    /      dTdx(i+1,j,k1)      |                  |
!      |  /              |             |                  |
!      |/                |             |                  |
!      |                 |             |                  |
!      |                 |       dTdx(i,j,k1)     dTdx(i+1,j,k1)  k
!      |                 |             |                  |
!      |                /              |                  |
!      |              /                |                  |
! dTdx(i,j,k1)      /                  |                  |
!      |          /                     -- dTdz(i,j,k1) --        k-1/2
!      |    dTdz(i,j,k1)
!      |      /
!      |    /
!      |  /   Physical grid,                logical grid
!      |/   X-Z cross-section


!  Add in horizontal biharmonic diffusion along rotated surfaces.
!  The biharmonic operator is computed by applying the rotated
!  Laplacian operator twice.
!

! THE FIRST ROTATED LAPLACIAN OPERATOR:
!---- ----- ------- --------- ----------


! Initialize computation of FSC flux appearing due to projection of
! along-isopycnal fluxes onto vertical direction by setting first to
! just metric term FSC=1/[drho/dz]. This metric term will later be
! converted later to the actual flux.

          if (k == 0) then
            do j=jmin-1,jmax+1
              do i=imin-1,imax+1
                FSC(i,j,k2)=0.
                dTdz(i,j,k2)=idRz(i,j,1)*( t(i,j,2,ntdf,itrc)
     &                                   - t(i,j,1,ntdf,itrc) )
              enddo
            enddo
          elseif (k < N) then
            do j=jmin-1,jmax+1
              do i=imin-1,imax+1
                FSC(i,j,k2)=idRz(i,j,k)*(z_r (i,j,k+1)-z_r (i,j,k))
                dTdz(i,j,k2)=idRz(i,j,k)*( t(i,j,k+1,ntdf,itrc)
     &                                       - t(i,j,k,ntdf,itrc) )
              enddo
            enddo
          elseif (k == N) then
            do j=jmin-1,jmax+1
              do i=imin-1,imax+1
                FSC(i,j,k2)=0.
                dTdz(i,j,k2)=dTdz(i,j,k1)
              enddo
            enddo
          endif

! Note that dTdx,dTde are computed one-step-ahead along k-loop before
! they are used to compute FX,FE below.  This migh seem to be off, but
! it should be kept in mind that dTdx,dTde are needed for dual purpose:
! (i) to compute FX,FE, and eventually, to compute slope-contributions
! to FSC further below, where both slices of dTdx,dTde will be used
! simultaneously.

          if (k < N) then
            do j=jmin,jmax
              do i=imin,imax+1
                dTdx(i,j,k2)=0.5*(pm(i,j)+pm(i-1,j))*(
     &            t(i,j,k+1,ntdf,itrc)-t(i-1,j,k+1,ntdf,itrc) )
#  ifdef MASKING
     &                                              *umask(i,j)
#  endif
              enddo
            enddo

            do j=jmin,jmax+1
              do i=imin,imax
                dTde(i,j,k2)=0.5*(pn(i,j)+pn(i,j-1))*(
     &            t(i,j,k+1,ntdf,itrc)-t(i,j-1,k+1,ntdf,itrc) )
#  ifdef MASKING
     &                                              *vmask(i,j)
#  endif
              enddo
            enddo
          endif  !<-- k < N

! Compute horizontal components of isopycnic diffusive flux

          if (k > 0) then
            do j=jmin,jmax
              do i=imin,imax+1
                FX(i,j)=diff3u(i,j,k)*0.5*(Hz(i,j,k)+Hz(i-1,j,k))
     &                                                  *dn_u(i,j)*(
     &                       dTdx(i,j,k1)
#  ifdef SW_TRIADS
     &       -0.5*( min(dRdx(i,j,k),0.)*(dTdz(i-1,j,k1)+dTdz(i,j,k2))
     &             +max(dRdx(i,j,k),0.)*(dTdz(i-1,j,k2)+dTdz(i,j,k1))
#  else
     &                -0.25*dRdx(i,j,k)*( dTdz(i-1,j,k1)+dTdz(i,j,k2)
     &                                   +dTdz(i-1,j,k2)+dTdz(i,j,k1)
#  endif
     &                                                            ) )
              enddo
            enddo

            do j=jmin,jmax+1
              do i=imin,imax
                FE(i,j)=diff3v(i,j,k)*0.5*(Hz(i,j,k)+Hz(i,j-1,k))
     &                                                  *dm_v(i,j)*(
     &                       dTde(i,j,k1)
#  ifdef SW_TRIADS
     &       -0.5*(min(dRde(i,j,k),0.)*(dTdz(i,j-1,k1)+dTdz(i,j,k2))
     &            +max(dRde(i,j,k),0.)*(dTdz(i,j-1,k2)+dTdz(i,j,k1)))
#  else
     &                -0.25*dRde(i,j,k)*(dTdz(i,j-1,k1)+dTdz(i,j,k2)
     &                                  +dTdz(i,j-1,k2)+dTdz(i,j,k1))
#  endif
     &                                                             )
              enddo
            enddo  !--> keep dTdx,dTde

! Finalize computation of the vertical component FSC

            if (k < N) then
              do j=jmin,jmax
                do i=imin,imax
#  ifdef SW_TRIADS
                  sumX=0. ; idx=0

                  if (dRdx(i,j,k) < 0.) then
                   sumX=diff3u(i,j,k)*dRdx(i,j,k)*(
     &                    dRdx(i,j,k)*dTdz(i,j,k2)-dTdx(i,j,k1) )
                   idx=1
                  endif

                  if (dRdx(i,j,k+1) > 0.) then
                   sumX=sumX + diff3u(i,j,k+1)*dRdx(i,j,k+1)*(
     &                           dRdx(i,j,k+1)*dTdz(i,j,k2)
     &                                        -dTdx(i,j,k2) )
                   idx=idx+1
                  endif

                  if (dRdx(i+1,j,k+1) < 0.) then
                   sumX=sumX + diff3u(i+1,j,k+1)*dRdx(i+1,j,k+1)*(
     &                           dRdx(i+1,j,k+1)*dTdz(i  ,j,k2)
     &                                          -dTdx(i+1,j,k2) )
                   idx=idx+1
                  endif

                  if (dRdx(i+1,j,k) > 0.) then
                   sumX=sumX + diff3u(i+1,j,k)*dRdx(i+1,j,k)*(
     &                           dRdx(i+1,j,k)*dTdz(i  ,j,k2)
     &                                        -dTdx(i+1,j,k1) )
                   idx=idx+1
                  endif
#  else
                  sumX=diff3u(i,j,k)*dRdx(i,j,k)*(
     &              dRdx(i,j,k)*dTdz(i,j,k2)-dTdx(i,j,k1) )

     &              +diff3u(i,j,k+1)*dRdx(i,j,k+1)*(
     &                   dRdx(i,j,k+1)*dTdz(i,j,k2)-dTdx(i,j,k2) )

     &              +diff3u(i+1,j,k+1)*dRdx(i+1,j,k+1)*(
     &               dRdx(i+1,j,k+1)*dTdz(i,j,k2)-dTdx(i+1,j,k2) )

     &              +diff3u(i+1,j,k)*dRdx(i+1,j,k)*(
     &               dRdx(i+1,j,k)*dTdz(i,j,k2)-dTdx(i+1,j,k1) )
#  endif

#  ifdef SW_TRIADS
                  sumE=0. ; ide=0

                  if (dRde(i,j,k) < 0.) then
                    sumE=diff3v(i,j,k)*dRde(i,j,k)*(
     &                     dRde(i,j,k)*dTdz(i,j,k2)
     &                                -dTde(i,j,k1) )
                   ide=1
                  endif

                  if (dRde(i,j,k+1) > 0.) then
                    sumE=sumE + diff3v(i,j,k+1)*dRde(i,j,k+1)*(
     &                            dRde(i,j,k+1)*dTdz(i,j,k2)
     &                                         -dTde(i,j,k2) )
                   ide=ide+1
                  endif

                  if (dRde(i,j+1,k+1) <  0.) then
                    sumE=sumE + diff3v(i,j+1,k+1)*dRde(i,j+1,k+1)*(
     &                            dRde(i,j+1,k+1)*dTdz(i,j  ,k2)
     &                                           -dTde(i,j+1,k2) )
                   ide=ide+1
                  endif

                  if (dRde(i,j+1,k) > 0.) then
                    sumE=sumE + diff3v(i,j+1,k)*dRde(i,j+1,k)*(
     &                            dRde(i,j+1,k)*dTdz(i,j  ,k2)
     &                                         -dTde(i,j+1,k1) )
                    ide=ide+1
                  endif

                  FSC(i,j,k2)=FSC(i,j,k2)*(sumX*wgt(idx)+sumE*wgt(ide))
#  else
                  sumE=diff3v(i,j,k)*dRde(i,j,k)*(
     &                  dRde(i,j,k)*dTdz(i,j,k2)-dTde(i,j,k1) )

     &                +diff3v(i,j,k+1)*dRde(i,j,k+1)*(
     &                  dRde(i,j,k+1)*dTdz(i,j,k2)-dTde(i,j,k2) )

     &                +diff3v(i,j+1,k+1)*dRde(i,j+1,k+1)*(
     &                dRde(i,j+1,k+1)*dTdz(i,j,k2)-dTde(i,j+1,k2) )

     &              +diff3v(i,j+1,k)*dRde(i,j+1,k)*(
     &                  dRde(i,j+1,k)*dTdz(i,j,k2)-dTde(i,j+1,k1) )

                  FSC(i,j,k2)=FSC(i,j,k2)*0.25*(sumX+sumE)
#  endif
                enddo
              enddo
            endif     !<--   k < N
            do j=jmin,jmax
              do i=imin,imax
                LapT(i,j,k)=( pm(i,j)*pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                         +FE(i,j+1)-FE(i,j)
     &                      )+FSC(i,j,k2)-FSC(i,j,k1) )/Hz(i,j,k)
              enddo
            enddo
          endif   ! <-- k > 0
        enddo   ! --> k


! Apply lateral boundary conditions to the Laplacian (for the cases
! other than periodic: closed or gradient.

#  ifndef EW_PERIODIC
        if (WESTERN_EDGE) then
          do k=1,N
            do j=jmin,jmax
#   ifndef OBC_WEST
              LapT(istr-1,j,k)=0.
#   else
              LapT(istr-1,j,k)=LapT(istr,j,k)
#   endif
            enddo
          enddo
        endif
        if (EASTERN_EDGE) then
          do k=1,N
            do j=jmin,jmax
#   ifndef OBC_EAST
              LapT(iend+1,j,k)=0.
#   else
              LapT(iend+1,j,k)=LapT(iend,j,k)
#   endif
            enddo
          enddo
        endif
#  endif
#  ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do k=1,N
            do i=imin,imax
#   ifndef OBC_SOUTH
              LapT(i,jstr-1,k)=0.
#   else
              LapT(i,jstr-1,k)=LapT(i,jstr,k)
#   endif
            enddo
          enddo
        endif
        if (NORTHERN_EDGE) then
          do k=1,N
            do i=imin,imax
#   ifndef OBC_NORTH
              LapT(i,jend+1,k)=0.
#   else
              LapT(i,jend+1,k)=LapT(i,jend,k)
#   endif
            enddo
          enddo
        endif
#  endif /* !NS_PERIODIC */

! The second rotated Laplacian operator:
!---- ------ ------- --------- ----------

        k2=1
        do k=0,N,+1
          k1=k2
          k2=3-k1

          if (k == 0) then
            do j=jstr-1,jend+1
              do i=istr-1,iend+1
                FSC(i,j,k2)=0.
                dTdz(i,j,k2)=idRz(i,j,1)*(LapT(i,j,2)-LapT(i,j,1))
              enddo
            enddo
          elseif (k < N) then
            do j=jstr-1,jend+1
              do i=istr-1,iend+1
                FSC(i,j,k2)= idRz(i,j,k)*(z_r (i,j,k+1)-z_r (i,j,k))
                dTdz(i,j,k2)=idRz(i,j,k)*(LapT(i,j,k+1)-LapT(i,j,k))
              enddo
            enddo
          elseif (k == N) then
            do j=jmin-1,jmax+1
              do i=imin-1,imax+1
                FSC(i,j,k2)=0.
                dTdz(i,j,k2)=dTdz(i,j,k1)
              enddo
            enddo
          endif

          if (k < N) then
            do j=jstr,jend
              do i=istr,iend+1
                dTdx(i,j,k2)=0.5*(pm(i,j)+pm(i-1,j))*( LapT(i,j,k+1)
     &                                              -LapT(i-1,j,k+1))
#  ifdef MASKING
     &                                                    *umask(i,j)
#  endif
              enddo
            enddo
            do j=jstr,jend+1
              do i=istr,iend
                dTde(i,j,k2)=0.5*(pn(i,j)+pn(i,j-1))*( LapT(i,j,k+1)
     &                                              -LapT(i,j-1,k+1))
#  ifdef MASKING
     &                                                    *vmask(i,j)
#  endif
              enddo
            enddo
          endif

          if (k > 0) then
            do j=jstr,jend
              do i=istr,iend+1
                FX(i,j)=-diff3u(i,j,k)*0.5*(Hz(i,j,k)+Hz(i-1,j,k))
     &                                                 *dn_u(i,j)*(
     &                      dTdx(i,j,k1)
#  ifdef SW_TRIADS
     &     -0.5*(min(dRdx(i,j,k),0.)*(dTdz(i-1,j,k1)+dTdz(i,j,k2))
     &          +max(dRdx(i,j,k),0.)*(dTdz(i-1,j,k2)+dTdz(i,j,k1)))
#  else
     &              -0.25*dRdx(i,j,k)*(dTdz(i-1,j,k1)+dTdz(i,j,k2)
     &                                +dTdz(i-1,j,k2)+dTdz(i,j,k1))
#  endif
     &                                                           )
              enddo
            enddo

            do j=jstr,jend+1
              do i=istr,iend
                FE(i,j)=-diff3v(i,j,k)*0.5*(Hz(i,j,k)+Hz(i,j-1,k))
     &                                                 *dm_v(i,j)*(
     &                      dTde(i,j,k1)
#  ifdef SW_TRIADS
     &     -0.5*(min(dRde(i,j,k),0.)*(dTdz(i,j-1,k1)+dTdz(i,j,k2))
     &          +max(dRde(i,j,k),0.)*(dTdz(i,j-1,k2)+dTdz(i,j,k1)))
#  else
     &              -0.25*dRde(i,j,k)*(dTdz(i,j-1,k1)+dTdz(i,j,k2)
     &                                +dTdz(i,j-1,k2)+dTdz(i,j,k1))
#  endif
     &                                                           )
              enddo
            enddo

            if (k < N) then    ! Compute the vertical component...
              do j=jstr,jend
                do i=istr,iend

#  ifdef STABILIZE

                  s2_XLL=(FSC(i,j,k2)*dRdx(i,j,  k  ))**2
                  s2_XLU=(FSC(i,j,k2)*dRdx(i,j,  k+1))**2
                  s2_XRU=(FSC(i,j,k2)*dRdx(i+1,j,k+1))**2
                  s2_XRL=(FSC(i,j,k2)*dRdx(i+1,j,k  ))**2

                  s2_ELL=(FSC(i,j,k2)*dRde(i,j  ,k  ))**2
                  s2_ELU=(FSC(i,j,k2)*dRde(i,j  ,k+1))**2
                  s2_ERU=(FSC(i,j,k2)*dRde(i,j+1,k+1))**2
                  s2_ERL=(FSC(i,j,k2)*dRde(i,j+1,k  ))**2

c**               cff=1./(z_r(i,j,k+1)-z_r(i,j,k))
                  cff=2./(Hz(i,j,k+1)+Hz(i,j,k))

                  cff2=cff**2
                  cffX=pm(i,j)**2
                  cffE=pn(i,j)**2

#   ifdef SW_TRIADS
                  Akz(i,j,k) = 15.*(  max( diff3u(i  ,j,k  )*s2_XLL,
     &                                     diff3u(i+1,j,k  )*s2_XRL,
     &                                     diff3u(i  ,j,k+1)*s2_XLU,
     &                                     diff3u(i+1,j,k+1)*s2_XRU)

     &                               +max( diff3v(i,j  ,k  )*s2_ELL,
     &                                     diff3v(i,j+1,k  )*s2_ERL,
     &                                     diff3v(i,j  ,k+1)*s2_ELU,
     &                                     diff3v(i,j+1,k+1)*s2_ERU)
     &                 )*(
     &                     max( diff3u(i  ,j,k  )*(cffX + cff2*s2_XLL),
     &                          diff3u(i  ,j,k+1)*(cffX + cff2*s2_XLU),
     &                          diff3u(i+1,j,k+1)*(cffX + cff2*s2_XRU),
     &                          diff3u(i+1,j,k  )*(cffX + cff2*s2_XRL))

     &                    +max( diff3v(i,j  ,k  )*(cffE + cff2*s2_ELL),
     &                          diff3v(i,j  ,k+1)*(cffE + cff2*s2_ELU),
     &                          diff3v(i,j+1,k+1)*(cffE + cff2*s2_ERU),
     &                          diff3v(i,j+1,k  )*(cffE + cff2*s2_ERL))
     &                   )
#   else
                  Akz(i,j,k) = 0.5*(  diff3u(i  ,j,k  )*s2_XLL
     &                               +diff3u(i+1,j,k  )*s2_XRL
     &                               +diff3u(i  ,j,k+1)*s2_XLU
     &                               +diff3u(i+1,j,k+1)*s2_XRU

     &                               +diff3v(i,j,k  )*s2_ELL
     &                               +diff3v(i,j+1,k)*s2_ERL
     &                               +diff3v(i,j,k+1)*s2_ELU
     &                               +diff3v(i,j+1,k+1)*s2_ERU
     &                      )*(
     &                          diff3u(i  ,j,k  )*(cffX + cff2*s2_XLL)
     &                         +diff3u(i  ,j,k+1)*(cffX + cff2*s2_XLU)
     &                         +diff3u(i+1,j,k+1)*(cffX + cff2*s2_XRU)
     &                         +diff3u(i+1,j,k  )*(cffX + cff2*s2_XRL)

     &                         +diff3v(i,j  ,k  )*(cffE + cff2*s2_ELL)
     &                         +diff3v(i,j  ,k+1)*(cffE + cff2*s2_ELU)
     &                         +diff3v(i,j+1,k+1)*(cffE + cff2*s2_ERU)
     &                         +diff3v(i,j+1,k  )*(cffE + cff2*s2_ERL)
     &                       )
#   endif
#  endif  /* STABILIZE */



#  ifdef SW_TRIADS
                  sumX=0. ; idx=0

                  if (dRdx(i,j,k) < 0.) then
                    sumX=diff3u(i,j,k)*dRdx(i,j,k)*(
     &                                 dRdx(i,j,k)*dTdz(i,j,k2)
     &                                            -dTdx(i,j,k1) )
                    idx=1
                  endif

                  if (dRdx(i,j,k+1) > 0.) then
                    sumX=sumX + diff3u(i,j,k+1)*dRdx(i,j,k+1)*(
     &                            dRdx(i,j,k+1)*dTdz(i,j,k2)
     &                                         -dTdx(i,j,k2) )
                    idx=idx+1
                  endif

                  if (dRdx(i+1,j,k+1) < 0.) then
                    sumX=sumX + diff3u(i+1,j,k+1)*dRdx(i+1,j,k+1)*(
     &                            dRdx(i+1,j,k+1)*dTdz(i  ,j,k2)
     &                                           -dTdx(i+1,j,k2) )
                    idx=idx+1
                  endif

                  if (dRdx(i+1,j,k  ) > 0.) then
                    sumX=sumX + diff3u(i+1,j,k)*dRdx(i+1,j,k)*(
     &                            dRdx(i+1,j,k)*dTdz(i  ,j,k2)
     &                                         -dTdx(i+1,j,k1) )
                    idx=idx+1
                  endif
#  else
                  sumX=diff3u(i  ,j,k  )*dRdx(i  ,j,k  )
     &                 *(dRdx(i  ,j,k  )*dTdz(i  ,j,k2 )-dTdx(i  ,j,k1))
     &                +diff3u(i  ,j,k+1)*dRdx(i  ,j,k+1)
     &                 *(dRdx(i  ,j,k+1)*dTdz(i  ,j,k2 )-dTdx(i  ,j,k2))
     &                +diff3u(i+1,j,k+1)*dRdx(i+1,j,k+1)
     &                 *(dRdx(i+1,j,k+1)*dTdz(i  ,j,k2 )-dTdx(i+1,j,k2))
     &                +diff3u(i+1,j,k  )*dRdx(i+1,j,k  )
     &                 *(dRdx(i+1,j,k  )*dTdz(i  ,j,k2 )-dTdx(i+1,j,k1))
#  endif
#  ifdef SW_TRIADS
                  sumE=0. ; ide=0

                  if (dRde(i,j,k) < 0.) then
                   sumE=diff3v(i,j,k)*dRde(i,j,k)*(
     &                                dRde(i,j,k)*dTdz(i,j,k2)
     &                                           -dTde(i,j,k1) )
                   ide=1
                  endif

                  if (dRde(i,j,k+1) > 0.) then
                   sumE=sumE + diff3v(i,j,k+1)*dRde(i,j,k+1)*(
     &                           dRde(i,j,k+1)*dTdz(i,j,k2)
     &                                        -dTde(i,j,k2) )
                   ide=ide+1
                  endif

                  if (dRde(i,j+1,k+1) < 0.) then
                   sumE=sumE + diff3v(i,j+1,k+1)*dRde(i,j+1,k+1)*(
     &                           dRde(i,j+1,k+1)*dTdz(i,j  ,k2)
     &                                          -dTde(i,j+1,k2) )
                   ide=ide+1
                  endif

                  if (dRde(i,j+1,k  ) > 0.) then
                   sumE=sumE + diff3v(i,j+1,k)*dRde(i,j+1,k)*(
     &                           dRde(i,j+1,k)*dTdz(i,j  ,k2)
     &                                        -dTde(i,j+1,k1) )
                   ide=ide+1
                  endif
#  else
                  sumE=diff3v(i,j  ,k  )*dRde(i,j  ,k  )
     &                 *(dRde(i,j  ,k  )*dTdz(i,j  ,k2 )-dTde(i,j ,k1))
     &                +diff3v(i,j  ,k+1)*dRde(i,j  ,k+1)
     &                 *(dRde(i,j  ,k+1)*dTdz(i,j  ,k2 )-dTde(i,j ,k2))
     &                +diff3v(i,j+1,k+1)*dRde(i,j+1,k+1)
     &                 *(dRde(i,j+1,k+1)*dTdz(i,j  ,k2)-dTde(i,j+1,k2))
     &                +diff3v(i,j+1,k  )*dRde(i,j+1,k )
     &                 *(dRde(i,j+1,k  )*dTdz(i,j  ,k2)-dTde(i,j+1,k1))
#  endif


#  ifdef SW_TRIADS
                 FSC(i,j,k2)=-FSC(i,j,k2)*(sumX*wgt(idx)+sumE*wgt(ide))
#  else
                 FSC(i,j,k2)=-FSC(i,j,k2)*0.25*(sumX+sumE)
#  endif
#  ifdef STABILIZE
     &                         -cff*Akz(i,j,k)*( t(i,j,k+1,ntdf,itrc)
     &                                            -t(i,j,k,ntdf,itrc))
#  endif

                enddo
              enddo
            endif    !<-- k < N
            do j=jstr,jend
              do i=istr,iend
                t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc) +dt*( pm(i,j)
     &               *pn(i,j)*(FX(i+1,j)-FX(i,j) +FE(i,j+1)-FE(i,j))
     &                                    +FSC(i,j,k2)-FSC(i,j,k1) )
              enddo
            enddo

# ifdef DIAGNOSTICS
            if (diag_trc .and. wrt_t_dia(itrc) .and. calc_diag) then
              do j=1,ny
                do i=1,nx+1
                  Tdiag(i,j,k,tmixx,td(itrc)) = -FX(i,j)
                enddo
              enddo
              do j=1,ny+1
                do i=1,nx
                  Tdiag(i,j,k,tmixy,td(itrc)) = -FE(i,j)
                enddo
              enddo
              do j=1,ny
                do i=1,nx
                  Tdiag(i,j,k,tmixz,td(itrc)) = -FSC(i,j,k2)/(pm(i,j)*pn(i,j))
                enddo
              enddo
            endif
# endif /* DIAGNOSTICS */

          endif   !<-- k > 0
# endif /* ADV_ISONEUTRAL */
        enddo      !<-- k



! Compute and apply vertical advective fluxes.

# ifdef DIAGNOSTICS
        if (diag_trc .and. wrt_t_dia(itrc) .and. calc_diag) then
          call diag_t_adv_vc4(itrc )
        endif
# endif

        do j=jstr,jend

# include "compute_vert_tracer_fluxes.h"

# ifdef DIAGNOSTICS
          if (diag_trc .and. wrt_t_dia(itrc) .and. calc_diag) then
            do k=1,nz
              do i=1,nx
                Tdiag(i,j,k,tadvz,td(itrc)) = FC(i,k)
                VFlxD(i,k) = t(i,j,k,nnew,itrc)
              enddo
            enddo
          endif

# endif


          do k=1,nz
            do i=istr,iend
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc) -dt*pm(i,j)*pn(i,j)
     &                                              *(FC(i,k)-FC(i,k-1))

              if (pipe_source) then
                if (pipe_idx(i,j)>0.) then                           ! pipe_source parameter will
                pidx = pipe_idx(i,j)                                 ! remove branching if false
                t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)
     &            +dt*pm(i,j)*pn(i,j)*
     &              pipe_flx(i,j)*pipe_prf(pidx,k)*pipe_trc(pidx,itrc)
                endif
              endif

            enddo
          enddo


! Start computation of vertical diffusion: add surface forcing fluxes

          !! JM add line for bottom tracer flux
          do i=istr,iend
            t(i,j,N,nnew,itrc)=t(i,j,N,nnew,itrc)+dt*stflx(i,j,itrc)
!           t(i,j,1,nnew,itrc)=t(i,j,1,nnew,itrc)+dt*btflx(i,j,itrc)
          enddo

# ifdef LMD_KPP
! Add the solar radiation flux in temperature equation. Also compute
! the nonlocal transport flux for unstable (convective) forcing
! conditions into matrix DC when using the Large et al. 1994 KPP
! scheme.

          if (itrc == itemp) then
            do k=N-1,1,-1
              do i=istr,iend
                cff=srflx(i,j)*swr_frac(i,j,k)
#  ifdef LMD_NONLOCAL
     &                 -ghat(i,j,k)*(stflx(i,j,itemp)-srflx(i,j))
#  endif
                t(i,j,k+1,nnew,itemp)=t(i,j,k+1,nnew,itemp) -dt*cff
                t(i,j,k  ,nnew,itemp)=t(i,j,k  ,nnew,itemp) +dt*cff
              enddo
            enddo


c?? WARNING: the above implies that light (if any) reaching all the
c?? way to the bottom is entirely absorbed within the botom-most grid
c?? box, rather that reflected/scattered back to the water column. It
c?? is not clear, however, how to deal with this situation correctly
c??
c??         do i=istr,iend
c??           cff=srflx(i,j)*swr_frac(i,j,0)
c??#  ifdef LMD_NONLOCAL
c??  &               -ghat(i,j,0)*(stflx(i,j,itemp)-srflx(i,j))
c??#  endif
c??           t(i,j,1,nnew,itemp)=t(i,j,1,nnew,itemp) -dt*cff
c??         enddo

#  if defined LMD_NONLOCAL && defined SALINITY
          elseif (itrc == isalt) then
            do k=N-1,1,-1
              do i=istr,iend
                cff=-dt*ghat(i,j,k)*stflx(i,j,isalt)
                t(i,j,k+1,nnew,isalt)=t(i,j,k+1,nnew,isalt) -cff
                t(i,j,k  ,nnew,isalt)=t(i,j,k  ,nnew,isalt) +cff
              enddo
            enddo
#  endif
          endif

# endif

! Perform an implicit time step for vertical diffusion
!
!   dq(k)     1     [         q(k+1)-q(k)             q(k)-q(k-1) ]
!  ------ = ----- * [ Akt(k)* ----------- - Akt(k-1)* ----------- ]
!    dt     Hz(k)   [            dz(k)                   dz(k-1)  ]
!
! where q(k) represents tracer field t(:,:,k,:,itrc). Doing so
! implies solution of a tri-diagonal system
!
!     -FC(k-1)*q_new(k-1) +[Hz(k)+FC(k-1)+FC(k)]*q_new(k)
!                       -FC(k)*q_new(k+1) = Hz(k)*q_old(k)
!
!                dt*Akt(k)
! where FC(k) = ----------- is normalized diffusivity coefficient
!                  dz(k)
!
! defined at W-points; q_new(k) is the new-time-step (unknown) tracer
! field; q_old(k) is old-time-step tracer (known).   As long as Akt(k)
! is nonnegative, the tri-diagonal matrix is diagonally dominant which
! guarantees stability of the Gaussian elimination algorithm.  Top and
! bottom boundary conditions are assumed to be no-flux, effectively
! Akt(N)=Akt(0)=0, hence FC(N)=FC(1)=0. This leads to equations for
! top and bottom grid boxes;
!
!   -FC(N-1)*q_new(N-1) +[Hz(N)+FC(N-1)]*q_new(N) = Hz(N)*q_old(N)
!
!          [Hz(1)+FC(1)]*q_new(1) -FC(1)*q_new(2) = Hz(1)*q_old(1)
!
! The FC(N)=FC(0)=0 boundary conditions does not mean that physical
! boundary conditions are no flux: the forcing fluxes have been
! applied explicitly above.  Instead, the no-flux condition should
! be interpreted as that the implicit step merely redistributes the
! tracer concentration throughout the water column. At this moment
! the content of array t(:,:,:,nnew,itrc) has meaning of Hz*tracer.
! After the implicit step it becomes just tracer.


          iAkt=min(itrc,iTandS)

          do i=istr,iend
            DC(i,0)=dt*pm(i,j)*pn(i,j) !<--horizontal metric

# if defined ADV_ISONEUTRAL && defined STABILIZE
            FC(i,1)=2.*dt*(Akt(i,j,1,iAkt)+Akz(i,j,1))
# else
            FC(i,1)=2.*dt* Akt(i,j,1,iAkt)
# endif
     &                           /(Hz(i,j,1)+Hz(i,j,2))
            WC(i,1)=DC(i,0)*Wi(i,j,1)

            cff=1./(Hz(i,j,1) +FC(i,1)+max(WC(i,1),0.))
            CF(i,1)=cff*(      FC(i,1)-min(WC(i,1),0.))
            DC(i,1)=cff*t(i,j,1,nnew,itrc)
          enddo
          do k=2,N-1,+1
            do i=istr,iend
# if defined ADV_ISONEUTRAL && defined STABILIZE
              FC(i,k)=2.*dt*(Akt(i,j,k,iAkt)+Akz(i,j,k))
# else
              FC(i,k)=2.*dt* Akt(i,j,k,iAkt)
# endif
     &                            /(Hz(i,j,k)+Hz(i,j,k+1))

              WC(i,k)=Wi(i,j,k)*DC(i,0)

              cff=1./( Hz(i,j,k) +FC(i,k)+max(WC(i,k),0.)
     &                         +FC(i,k-1)-min(WC(i,k-1),0.)
     &              -CF(i,k-1)*(FC(i,k-1)+max(WC(i,k-1),0.))
     &                                                     )

              CF(i,k)=cff*(       FC(i,k)-min(WC(i,k),0.) )

              DC(i,k)=cff*( t(i,j,k,nnew,itrc) +DC(i,k-1)*(
     &                          FC(i,k-1)+max(WC(i,k-1),0.) ))
            enddo
          enddo          !--> discard DC(:,0)
          do i=istr,iend
            t(i,j,N,nnew,itrc)=( t(i,j,N,nnew,itrc) +DC(i,N-1)*(
     &                               FC(i,N-1)+max(WC(i,N-1),0.) )
     &                   )/( Hz(i,j,N) +FC(i,N-1)-min(WC(i,N-1),0.)
     &                      -CF(i,N-1)*(FC(i,N-1)+max(WC(i,N-1),0.))
     &                                                             )
# ifdef MASKING
     &                                                  *rmask(i,j)
# endif
          enddo
          do k=N-1,1,-1
            do i=istr,iend
              t(i,j,k,nnew,itrc)=(DC(i,k)+CF(i,k)*t(i,j,k+1,nnew,itrc))
# ifdef MASKING
     &                                                      *rmask(i,j)
# endif
            enddo
          enddo           !--> discard FC,WC,CF,DC

# ifdef DIAGNOSTICS
          if (diag_trc .and. wrt_t_dia(itrc) .and. calc_diag) then 
            do k=1,nz
              do i=istr,iend
                VFlxD(i,k) =
     &            (VFlxD(i,k) - Hz(i,j,k)*t(i,j,k,nnew,itrc))/dt
              enddo
            enddo
            ! Convert divergence to fluxes on w-points
            do i=istr,iend
              ZFlx(i,1)= VFlxD(i,1)/(pm(i,j)*pn(i,j)) 
              do k=2,nz
                ZFlx(i,k)= VFlxD(i,k)/(pm(i,j)*pn(i,j)) + ZFlx(i,k-1)
              enddo
            enddo
            if (j==20.and.mynode==0) then
!              print *,maxval(Tdiag(1:nx,j,:,tmixz,td(itrc)))
            endif
            do i=istr,iend
              do k=1,nz
!               Tdiag(i,j,k,vmixd,td(itrc)) = 
!    &               ZFlx(i,k) - Tdiag(i,j,k,tadvz,td(itrc)) 
!               Tdiag(i,j,k,tmixz,td(itrc)) = 
!    &               ZFlx(i,k) - Tdiag(i,j,k,tadvz,td(itrc)) 

#  ifdef ADV_ISONEUTRAL
                Tdiag(i,j,k,tmixz,td(itrc)) = Tdiag(i,j,k,tmixz,td(itrc)) +
     &               ZFlx(i,k) - Tdiag(i,j,k,tadvz,td(itrc)) 
#  else
                Tdiag(i,j,k,tmixz,td(itrc)) = 
     &               ZFlx(i,k) - Tdiag(i,j,k,tadvz,td(itrc)) 
#  endif
              enddo
            enddo

          endif
# endif

#ifdef CONST_TRACERS
          do k=1,N
            do i=istr,iend
              t(i,j,k,nnew,itrc)=t(i,j,k,nstp,itrc)
            enddo
          enddo
#endif
        enddo      !<-- j

        call t3dbc_tile (istr,iend,jstr,jend, itrc, wrk1)
# ifdef DIAGNOSTICS
!       if (diag_trc .and. wrt_t_dia(itrc) .and. calc_diag)
!    &    call set_diags_t_bc(itrc,td)
# endif

      enddo     !<-- itrc

# if defined(MARBL)

      if (mynode==0) then
      print *, 'doing BGC with MARBL'
      end if
      
!     Pass indices, tracer and mask arrays, diagnostics arrays, and saved state arrays to MARBL
!      if (mynode==0) then
!         print *, shape(t), ' is the shape of t'
!         print *, lbound(t), ' is the lower bound of t'
!      end if
      call marbldrv_column_physics(istr,iend,jstr,jend,t)
# elif defined(BIOLOGY_BEC2)
      if (mynode==0) then
      print *, 'doing BGC with BEC'
      end if
      
      call ecosys_bec2_tile(istr,iend,jstr,jend) ! BEC, 2014
# endif

# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      ! for BGC, we should be able to pack more tracer arrays
      ! in a single mpi_exchange
      do itrc=1,nt-1,2
        call exchange_xxx( t(:,:,:,nnew,itrc),
     &                     t(:,:,:,nnew,itrc+1) )
      enddo
      if (mod(nt,2) == 1) then
        call exchange_xxx( t(:,:,:,nnew,nt) )
      endif
# endif
      end

      subroutine check_step_t_switches(ierr)

! The following routine keeps track of the status of CPP settings in
! "compute_horiz_tracer_fluxes.h" and "compute_horiz_vert_fluxes.h".
! This code must be placed here rather than in a separate file so it
! is exposed to the relevant CPP-settings.   It does not affect any
! model results, other than signature in global attribute "CPPS" in
! the output netCDF files.

      use param
      use strings

      implicit none
      integer ierr, is,ie, lenstr

      ie=lenstr(cpps)
      is=ie+2 ; ie=is+15
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='<step3d_t_ISO.F>'

# include "track_advec_switches.h"

      return
  99  mpi_master_only write(*,'(/1x,2A/12x,A/)')      '### ERROR: ',
     &  'Insufficient length of string "cpps" in file "strings".',
     &        'Increase parameter "max_opt_size" it and recompile.'
      ierr=ierr+1
      end
#else
      subroutine step3d_t_empty
      end
#endif  /* SOLVE3D */
