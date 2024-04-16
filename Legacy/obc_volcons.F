#include "cppdefs.h"
#ifdef OBC_VOLCONS

      subroutine obc_flux_tile (istr,iend,jstr,jend, dFlx)

! Compute integral crossection of open all boundaries "bc_crss" and
! integral mass flux "bc_flux" incoming into model domain (hence
! POSITIVE "bc_flux" means GAIN in mass) in order to enforce global
! mass conservation constraint. This is done by first computing the
! integral crossection of- and fluxes across-  the segments of open
! boundaries within subdomains [tiles], then making global summation
! over the subdomains. At the end, global integrated flux is divided
! by the global cross-section to compute the correction velocity
! "ubar_xs" needed to cancel exactly global mass gain.

      implicit none
      integer istr,iend,jstr,jend,  i,j, inc,jnc, nsubs
      real*QUAD dFlx(PRIVATE_2D_SCRATCH_ARRAY), my_crss, my_flux
# include "param.h"
# include "grid.h"
# include "ocean2d.h"
# include "scalars.h"
# ifdef MPI
#  include "mpif.h"
      integer size, step, status(MPI_STATUS_SIZE), ierr
      real*QUAD buff(2)
# endif
                          ! Reset partial sums over the segments
      my_crss=QuadZero    ! of open boundaries, then compute partial
      my_flux=QuadZero    ! cross-sections and mass fluxes over the
                          ! subdomains
# ifdef OBC_WEST
      if (WESTERN_EDGE) then
        do j=jstr,jend
          dFlx(istr,j)=0.5*( zeta(istr-1,j,knew)+zeta(istr,j,knew)
     &                        +h(istr-1,j)+h(istr,j) )*dn_u(istr,j)
#  ifdef MASKING
     &                                               *umask(istr,j)
#  endif
          dFlx(istr-1,j)=dFlx(istr,j)*ubar(istr,j,knew)
c**       my_crss=my_crss + dFlx(istr  ,j)
c**       my_flux=my_flux + dFlx(istr-1,j)
        enddo
        jnc=1
        do while(jstr<=jend-jnc)
          do j=jstr,jend-jnc,2*jnc
            dFlx(istr  ,j) = dFlx(istr  ,j) + dFlx(istr  ,j+jnc)
            dFlx(istr-1,j) = dFlx(istr-1,j) + dFlx(istr-1,j+jnc)
          enddo
          jnc=2*jnc
        enddo
        my_crss=my_crss + dFlx(istr  ,jstr)
        my_flux=my_flux + dFlx(istr-1,jstr)
      endif
# endif
# ifdef OBC_EAST
      if (EASTERN_EDGE) then
        do j=jstr,jend
          dFlx(iend,j)=0.5*( zeta(iend,j,knew)+zeta(iend+1,j,knew)
     &                      +h(iend,j)+h(iend+1,j) )*dn_u(iend+1,j)
#  ifdef MASKING
     &                                             *umask(iend+1,j)
#  endif
          dFlx(iend+1,j)=dFlx(iend,j)*ubar(iend+1,j,knew)
c**       my_crss=my_crss + dFlx(iend  ,j)
c**       my_flux=my_flux - dFlx(iend+1,j)
        enddo
        jnc=1
        do while(jstr<=jend-jnc)
          do j=jstr,jend-jnc,2*jnc
            dFlx(iend  ,j) = dFlx(iend  ,j) + dFlx(iend  ,j+jnc)
            dFlx(iend+1,j) = dFlx(iend+1,j) + dFlx(iend+1,j+jnc)
          enddo
          jnc=2*jnc
        enddo
        my_crss=my_crss + dFlx(iend  ,jstr)
        my_flux=my_flux - dFlx(iend+1,jstr)
      endif
# endif
# ifdef OBC_SOUTH
      if (SOUTHERN_EDGE) then
        do i=istr,iend
          dFlx(i,jstr)=0.5*( zeta(i,jstr,knew)+zeta(i,jstr-1,knew)
     &                        +h(i,jstr)+h(i,jstr-1) )*dm_v(i,jstr)
#  ifdef MASKING
     &                                               *vmask(i,jstr)
#  endif
          dFlx(i,jstr-1)=dFlx(i,jstr)*vbar(i,jstr,knew)
c**       my_crss=my_crss + dFlx(i,jstr  )
c**       my_flux=my_flux + dFlx(i,jstr-1)
        enddo
        inc=1
        do while(istr<=iend-inc)
          do i=istr,iend-inc,2*inc
            dFlx(i,jstr  ) = dFlx(i,jstr  ) + dFlx(i+inc,jstr  )
            dFlx(i,jstr-1) = dFlx(i,jstr-1) + dFlx(i+inc,jstr-1)
          enddo
          inc=2*inc
        enddo
        my_crss=my_crss + dFlx(istr,jstr  )
        my_flux=my_flux + dFlx(istr,jstr-1)
      endif
# endif
# ifdef OBC_NORTH
      if (NORTHERN_EDGE) then
        do i=istr,iend
          dFlx(i,jend)=0.5*( zeta(i,jend,knew)+zeta(i,jend+1,knew)
     &                      +h(i,jend)+h(i,jend+1) )*dm_v(i,jend+1)
#  ifdef MASKING
     &                                             *vmask(i,jend+1)
#  endif
          dFlx(i,jend+1)=dFlx(i,jend)*vbar(i,jend+1,knew)

c**       my_crss=my_crss + dFlx(i,jend  )
c**       my_flux=my_flux - dFlx(i,jend+1)
        enddo
        inc=1
        do while(istr<=iend-inc)
          do i=istr,iend-inc,2*inc
             dFlx(i,jend  ) = dFlx(i,jend  ) + dFlx(i+inc,jend  )
             dFlx(i,jend+1) = dFlx(i,jend+1) + dFlx(i+inc,jend+1)
          enddo
          inc=2*inc
        enddo
        my_crss=my_crss + dFlx(istr,jend  )
        my_flux=my_flux - dFlx(istr,jend+1)
      endif
# endif

      if (SINGLE_TILE_MODE) then    ! Global summation is performed
        nsubs=1                     ! in two stages, first within
      else                          ! each shared memory (MPI-nodes),
        nsubs=NSUB_X*NSUB_E         ! then across  MPI nodes (Reduce
      endif                         ! -- Broadcast sequence).

C$OMP CRITICAL (obc_flx_cr_rgn)
        if (bc_count==0) then     ! Counter tile_count identifies
          bc_flux=QuadZero          ! the last thread (within each
          bc_crss=QuadZero          ! MPI-process) leaving critical
        endif                       ! region. This thread (which is
        bc_count=bc_count+1         ! not necessarily master thread
        bc_crss=bc_crss+my_crss     ! within its MPI-process) is
        bc_flux=bc_flux+my_flux     ! responsible for communication
        if (bc_count==nsubs) then ! exchange between MPI-nodes.
          bc_count=0
# ifdef MPI
          size=NNODES               ! Reduce-broadcast sequence...
  1        step=(size+1)/2
            if (mynode>=step .and. mynode<size) then
              buff(1)=bc_crss
              buff(2)=bc_flux
              call MPI_Send (buff,  4, MPI_DOUBLE_PRECISION,
     &             mynode-step, 17, ocean_grid_comm,      ierr)
            elseif (mynode < size-step) then
              call MPI_Recv (buff,  4, MPI_DOUBLE_PRECISION,
     &             mynode+step, 17, ocean_grid_comm, status, ierr)

              bc_crss=bc_crss+buff(1)
              bc_flux=bc_flux+buff(2)
            endif
           size=step
          if (size>1) goto 1

          buff(1)=bc_crss
          buff(2)=bc_flux
          call MPI_Bcast(buff, 4, MPI_DOUBLE_PRECISION,
     &                         0, ocean_grid_comm, ierr)
          bc_crss=buff(1)
          bc_flux=buff(2)
# endif
          ubar_xs=bc_flux/bc_crss
        endif
C$OMP END CRITICAL (obc_flx_cr_rgn)
      end
#endif
#if defined OBC_FLUX_CORR || defined OBC_VOLCONS

      subroutine set_DUV_bc_tile (istr,iend,jstr,jend, Drhs, DUon,
     &                                                       DVom)

! Correct vertically integrated mass fluxes DUon and DVom along the
! open boundaries in such a way that the integral volume is conserved.
! This is done by applying ubar_xs correction to the velocities.

      implicit none
# include "param.h"
      integer istr,iend,jstr,jend, i,j
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: Drhs,DUon,DVom
# include "grid.h"
# include "ocean2d.h"
# include "scalars.h"

# include "compute_auxiliary_bounds.h"

# ifdef OBC_WEST
      if (WESTERN_EDGE) then
        do j=jstr-1,jendR
          DUon(istr,j)=DUon(istr,j) - ubar_xs*dn_u(istr,j)
     &                         *0.5*(Drhs(istr,j)+Drhs(istr-1,j))
#  ifdef MASKING
     &                                             *umask(istr,j)
#  endif
        enddo
      endif
# endif
# ifdef OBC_EAST
      if (EASTERN_EDGE) then
        do j=jstr-1,jendR
          DUon(iend+1,j)=DUon(iend+1,j) + ubar_xs*dn_u(iend+1,j)
     &                         *0.5*(Drhs(iend+1,j)+Drhs(iend,j))
#  ifdef MASKING
     &                                           *umask(iend+1,j)
#  endif
        enddo
      endif
# endif
# ifdef OBC_SOUTH
      if (SOUTHERN_EDGE) then
        do i=istr-1,iendR
          DVom(i,jstr)=DVom(i,jstr) - ubar_xs*dm_v(i,jstr)
     &                         *0.5*(Drhs(i,jstr)+Drhs(i,jstr-1))
#  ifdef MASKING
     &                                             *vmask(i,jstr)
#  endif
        enddo
      endif
# endif
# ifdef OBC_NORTH
      if (NORTHERN_EDGE) then
        do i=istr-1,iendR
          DVom(i,jend+1)=DVom(i,jend+1) + ubar_xs*dm_v(i,jend+1)
     &                         *0.5*(Drhs(i,jend+1)+Drhs(i,jend))
#  ifdef MASKING
     &                                           *vmask(i,jend+1)
#  endif
        enddo
      endif
# endif
      end
#else
      subroutine obc_volcons_empty
      end
#endif
