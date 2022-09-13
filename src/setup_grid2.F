#include "cppdefs.opt"
#define SUM_BY_PAIRS

      subroutine setup_grid2 (tile)
      use param
      use hidden_mpi_vars
      use private_scratch
      implicit none
      integer tile

#include "compute_tile_bounds.h"
      call setup_grid2_tile (istr,iend,jstr,jend, A2d(1,1),A2d(1,3))
      end

      subroutine setup_grid2_tile (istr,iend,jstr,jend, dA,dV)

! Find minimum and maximum depth of model topography; minimum and
! maximum horizontal XI- and ETA-grid intervals around h-points;
! Check minimum and maximum values of Courant number for barotropic
! mode (based on surface gravity wave speed, barotropic time step
! and local grid intervals); combined crossection of all open
! boundaries.

      use param
      use comm_vars
      use grid
      use scalars
      use mpi

      implicit none
      integer istr,iend,jstr,jend, i,j, nsub
      real(kind=QUAD), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: dA,dV
      real(kind=QUAD)  my_area,  my_volume,  my_crss
      real my_hmax, my_grdmax, my_Cu_max, cff,
     &     my_hmin, my_grdmin, my_Cu_min, my_Cu_Cor
#ifdef SUM_BY_PAIRS
      integer is,isize,itg, js,jsize,jtg
# if defined OBC_WEST || defined OBC_EAST
     &                    , jnc
# endif
# if defined OBC_SOUTH || defined OBC_NORTH
     &                    , inc
# endif
#endif
#ifdef MPI
# if QUAD==16
#  define MPI_BUFF_TYPE MPI_REAL16
# else
#  define MPI_BUFF_TYPE MPI_REAL8
# endif
      integer size, step, itag, status(MPI_STATUS_SIZE), ierr
      real*QUAD buff(16)
#endif

#include "compute_auxiliary_bounds.h"

      my_hmin=+1.E+20              ! Here two-dimensional Courant
      my_hmax=-1.E+20              ! number is defined as
      my_grdmin=+1.E+20            !
      my_grdmax=-1.E+20            ! Cu = c*dt * sqrt(1/dx^2+1/dy^2)
      my_Cu_min=+1.E+20            !
      my_Cu_max=-1.E+20            ! where c=sqrt(g*h) is phase
      my_Cu_Cor=-1.E+20            ! speed for barotropic mode,
                                   ! and dx,dy are grid spacings
      do j=jstrR,jendR             ! in each direction.
        do i=istrR,iendR
#ifdef MASKING
          if (rmask(i,j) > 0.5) then
#endif
            my_hmin=min(my_hmin, h(i,j))
            my_hmax=max(my_hmax, h(i,j))

            cff=1./sqrt(pm(i,j)*pn(i,j))
            my_grdmin=min(my_grdmin, cff)
            my_grdmax=max(my_grdmax, cff)

            cff=dtfast*sqrt(g*h(i,j)*(pm(i,j)*pm(i,j)+pn(i,j)*pn(i,j)))
            my_Cu_min=min(my_Cu_min, cff)
            my_Cu_max=max(my_Cu_max, cff)

            cff=dt*abs(f(i,j))
            my_Cu_Cor=max(my_Cu_Cor, cff)
#ifdef MASKING
          endif
#endif
        enddo
      enddo                        ! Compute unmasked total area and

                                   ! unperturbed volume of model grid
      do j=jstr,jend               ! (assuming zero free-surface).
        do i=istr,iend
#ifdef MASKING
          dA(i,j)=rmask(i,j)/(pm(i,j)*pn(i,j))
#else
          dA(i,j)=1./(pm(i,j)*pn(i,j))
#endif
          dV(i,j)=dA(i,j)*h(i,j)
        enddo
      enddo

#ifdef SUM_BY_PAIRS
! Perform summation of all dA within the tile using split-directional
! reduction-by-pairs algorithm to avoid accumulation of roundoff errors.

      isize=iend-istr  ; jsize=jend-jstr

      do while (isize>0 .or. jsize>0)
        if (jsize > 0) then
          js=(jsize+1)/2-1
          do j=0,js
            jtg=jstr+j
            do i=istr,istr+isize
              dA(i,jtg)=dA(i,jtg+j)+dA(i,jtg+j+1)
              dV(i,jtg)=dV(i,jtg+j)+dV(i,jtg+j+1)
            enddo
          enddo
          if (2*js+1 < jsize) then
            js=js+1
            jtg=jstr+js
            do i=istr,istr+isize
              dA(i,jtg)=dA(i,jtg+js)
              dV(i,jtg)=dV(i,jtg+js)
            enddo
          endif
          jsize=js
        endif

        if (isize > 0) then
          is=(isize+1)/2-1
          do j=jstr,jstr+jsize
            do i=0,is
              itg=istr+i
              dA(itg,j)=dA(itg+i,j)+dA(itg+i+1,j)
              dV(itg,j)=dV(itg+i,j)+dV(itg+i+1,j)
            enddo
          enddo
          if (2*is+1 < isize) then
            is=is+1
            itg=istr+is
            do j=jstr,jstr+jsize
              dA(itg,j)=dA(itg+is,j)
              dV(itg,j)=dV(itg+is,j)
            enddo
          endif
          isize=is
        endif
      enddo

      my_area  =dA(istr,jstr)
      my_volume=dV(istr,jstr)

#else
      my_area=QuadZero                  ! Straightforward summation:
      my_volume=QuadZero                ! this should produce the same
      do j=jstr,jend                    ! result as the above with the
        do i=istr,iend                  ! only difference is due to
          my_area  =my_area  +dA(i,j)   ! roundoff errors.
          my_volume=my_volume+dV(i,j)
        enddo
      enddo
#endif
                                        ! Compute integral crossection
      my_crss=QuadZero                  ! of all open boundaries...
#ifdef OBC_WEST
      if (WESTERN_EDGE) then
        do j=jstr,jend
          dA(istr,j)=0.5*(h(istr-1,j)+h(istr,j))*dn_u(istr,j)
# ifdef MASKING
     &                                          *umask(istr,j)
# endif
# ifndef SUM_BY_PAIRS
          my_crss=my_crss + dA(istr,j)
# endif
        enddo
# ifdef SUM_BY_PAIRS
        jnc=1
        do while(jstr<=jend-jnc)
          js=2*jnc
          do j=jstr,jend-jnc,js
            dA(istr,j) = dA(istr,j) + dA(istr,j+jnc)
          enddo
          jnc=js
        enddo
        my_crss=my_crss + dA(istr,jstr)
# endif
      endif
#endif
#ifdef OBC_EAST
      if (EASTERN_EDGE) then
        do j=jstr,jend
          dA(iend,j)=0.5*(h(iend,j)+h(iend+1,j))*dn_u(iend+1,j)
# ifdef MASKING
     &                                         *umask(iend+1,j)
# endif
# ifndef SUM_BY_PAIRS
          my_crss=my_crss + dA(iend,j)
# endif
        enddo
# ifdef SUM_BY_PAIRS
        jnc=1
        do while(jstr<=jend-jnc)
          js=2*jnc
          do j=jstr,jend-jnc,js
            dA(iend,j) = dA(iend,j) + dA(iend,j+jnc)
          enddo
          jnc=js
        enddo
        my_crss=my_crss + dA(iend,jstr)
# endif
      endif
#endif
#ifdef OBC_SOUTH
      if (SOUTHERN_EDGE) then
        do i=istr,iend
          dA(i,jstr)=0.5*(h(i,jstr)+h(i,jstr-1))*dm_v(i,jstr)
# ifdef MASKING
     &                                         *vmask(i,jstr)
# endif
# ifndef SUM_BY_PAIRS
          my_crss=my_crss + dA(i,jstr)
# endif
        enddo
# ifdef SUM_BY_PAIRS
        inc=1
        do while(istr<=iend-inc)
          is=2*inc
          do i=istr,iend-inc,is
            dA(i,jstr) = dA(i,jstr) + dA(i+inc,jstr)
          enddo
          inc=is
        enddo
        my_crss=my_crss + dA(istr,jstr)
# endif
      endif
#endif
#ifdef OBC_NORTH
      if (NORTHERN_EDGE) then
        do i=istr,iend
          dA(i,jend)=0.5*(h(i,jend)+h(i,jend+1))*dm_v(i,jend+1)
# ifdef MASKING
     &                                         *vmask(i,jend+1)
# endif
# ifndef SUM_BY_PAIRS
          my_crss=my_crss + dA(i,jend)
# endif
        enddo
# ifdef SUM_BY_PAIRS
        inc=1
        do while(istr<=iend-inc)
          is=2*inc
          do i=istr,iend-inc,is
            dA(i,jend) = dA(i,jend) + dA(i+inc,jend)
          enddo
          inc=is
        enddo
        my_crss=my_crss + dA(istr,jend)
# endif
      endif
#endif

! Thus far everything has ben computed within the tile. Next: combine
! the results to get global min/max and summations.  This needs to be
! done in two stages, first among all tiles which belong to the same
! shared memory group (MPI-node); within each MPI-node, then across
! MPI nodes (Reduce--Broadcast sequence).

      if (SINGLE_TILE_MODE) then
        nsub=1
      else
        nsub=NSUB_X*NSUB_E
      endif

C$OMP CRITICAL (grd2_cr_rgn)
        if (tile_count==0) then       ! Global min, max, and sum
          hmin=my_hmin                  ! operations within each shared
          hmax=my_hmax                  ! memory group (an MPI-node).
          grdmin=my_grdmin              ! Counter "tile_count" is to
          grdmax=my_grdmax              ! identify the first thread
          Cg_min=my_Cu_min              ! entering the critical region
          Cg_max=my_Cu_max              ! and the last one to leave
          Cu_Cor=my_Cu_Cor              ! (neither one is necessarily
          area=my_area                  ! the master thread within its
          volume=my_volume              ! MPI-process). In the case of
          bc_crss=my_crss               ! MPI code the last thread is
        else                            ! responsible for communication
          hmin=min(hmin, my_hmin)       ! between MPI nodes to finalize
          hmax=max(hmax, my_hmax)       ! the global operations.
          grdmin=min(grdmin, my_grdmin)
          grdmax=max(grdmax, my_grdmax)
          Cg_min=min(Cg_min, my_Cu_min)
          Cg_max=max(Cg_max, my_Cu_max)
          Cu_Cor=max(Cu_Cor, my_Cu_Cor)
          area=area+my_area
          volume=volume+my_volume
          bc_crss=bc_crss+my_crss
        endif

        tile_count=tile_count+1

        if (tile_count==nsub) then
          tile_count=0
#ifdef MPI
          size=NNODES
          do while (size>1)
           step=(size+1)/2
            if (mynode>=step .and. mynode<size) then
              buff(1)=hmin
              buff(2)=hmax               ! This is similar MPI_Reduce
              buff(3)=grdmin             ! MIN or MAX operation, except
              buff(4)=grdmax             ! that MIN/MAX and QUAD
              buff(5)=Cg_min             ! summations come together.
              buff(6)=Cg_max
              buff(7)=Cu_Cor
              buff(8)=area
              buff(9)=volume
              buff(10)=bc_crss

              itag=mynode+300
              call MPI_Send (buff, 10, MPI_BUFF_TYPE, mynode-step,
     &                       itag, ocean_grid_comm,          ierr)
            elseif (mynode < size-step) then
              itag=mynode+step+300
              call MPI_Recv (buff, 10, MPI_BUFF_TYPE, mynode+step,
     &                       itag, ocean_grid_comm,  status, ierr)
              cff=buff(1)
              hmin=  min(hmin,   cff)
              cff=buff(2)                 ! Copying buff() into
              hmax=  max(hmax,   cff)     ! "cff" is needed here
              cff=buff(3)                 ! to strip down QUAD
              grdmin=min(grdmin, cff)     ! precision in case if
              cff=buff(4)                 ! intrinsic min/max
              grdmax=max(grdmax, cff)     ! operations in some
              cff=buff(5)                 ! compilers cannot
              Cg_min=min(Cg_min, cff)     ! handle QUADs.

              cff=buff(6)
              Cg_max=max(Cg_max, cff)
              cff= buff(7)
              Cu_Cor=max(Cu_Cor, cff)
              area=area + buff(8)
              volume=volume + buff(9)
              bc_crss=bc_crss + buff(10)
            endif
           size=step
          enddo      !<--  while(size>1)

          buff(1)=hmin
          buff(2)=hmax
          buff(3)=grdmin
          buff(4)=grdmax
          buff(5)=Cg_min
          buff(6)=Cg_max
          buff(7)=Cu_Cor
          buff(8)=area
          buff(9)=volume
          buff(10)=bc_crss

          call MPI_Bcast(buff, 10, MPI_BUFF_TYPE, 0,
     &                       ocean_grid_comm, ierr)
          hmin=  buff(1)
          hmax=  buff(2)
          grdmin=buff(3)
          grdmax=buff(4)
          Cg_min=buff(5)
          Cg_max=buff(6)
          Cu_Cor=buff(7)
          area=buff(8)
          volume=buff(9)
          bc_crss=buff(10)

          if (mynode==0) then
#endif
            write(*,'(1x,A,F12.6,3x,A,ES14.7,5x,A,ES23.16)')
     &       'hmin =' ,hmin,  'grdmin =', grdmin, 'area =',   area
            write(*,'(1x,A,F12.6,3x,A,ES14.7,3x,A,ES23.16)')
     &       'hmax =' ,hmax,  'grdmax =', grdmax, 'volume =', volume
            write(*,'(43x,A,ES23.16)')        'open_cross =', bc_crss
            write(*,'(1x,A,F10.7,3x,A,F10.7,3x,A,F10.7)')
     &       'Cg_max =',Cg_max, 'Cg_min =',Cg_min, 'Cu_Cor =',Cu_Cor
#ifdef MPI
          endif
#endif
        endif
C$OMP END CRITICAL (grd2_cr_rgn)

      end
