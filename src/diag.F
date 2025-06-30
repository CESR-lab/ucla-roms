#include "cppdefs.opt"

! Compute and report various diagnostics: volume-averaged kinetic and
! potential energy, horizontally averaged free-surface perturbation;
! maximum velocity, maximum Courant number and it location, etc. Since
! this operation involves computation of  global sums, it is done in
! three stages: at first, summation within the current tile [subdomain
! of indices (istr:iend,jstr:jend)] done independently by an individual
! thread (in 3D cases this stage also performs vertical integration);
! then summation across the threads of the same MPI process (if any);
! and, finally, MPI reduction to compute global integrals.


#define SUM_BY_PAIRS
#define MAX_ADV_CFL
c--#define MAX_BY_VERT_CFL_ONLY

#define REPORT_MAX_UV
#define REPORT_POTEN_ENRG


      subroutine diag
      use param
      use hidden_mpi_vars
      use private_scratch
      use scalars
      implicit none


      integer   i,j
      integer,save:: tile=0

#include "compute_tile_bounds.h"

      j=min(iic,2*ninfo)-1             ! This introduces "logarithmic"
      i=1                              ! scaling of diag call frequency
      do while(i<j)                    ! at the beginning of model run.
        i=2*i
      enddo
      if (mod(iic,min(i,ninfo)) == 0) then
        call diag_tile (istr,iend,jstr,jend, A2d(1,1),A2d(1,2),A2d(1,3)
#ifdef SOLVE3D
     &               , A2d(1,4), A2d(1,5), A2d(1,6), A2d(1,7), A2d(1,8)
#endif
     &                                                                )
      endif
      end

      subroutine diag_tile (istr,iend,jstr,jend, dVol, ke,pe
#ifdef SOLVE3D
     &                            , ke2b, ke3bc, kesrf, ub,vb
#endif
     &                                                      )

      use param
      use comm_vars
      use diag_vars
      use eos_vars
      use grid
      use ocean_vars
      use scalars
      use mpi

      implicit none

#include "diag.opt"

      integer istr,iend,jstr,jend
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: dVol, ke, pe
#ifdef SOLVE3D
     &                             , ke2b, ke3bc, kesrf, ub, vb
#endif
#ifdef MPI
# if QUAD==16
#  define MPI_BUFF_TYPE MPI_REAL16
# else
#  define MPI_BUFF_TYPE MPI_REAL8
# endif
      integer size, step, itag, status(MPI_STATUS_SIZE)
      real*QUAD buff(16)
#endif

      real cff, dA, v2, my_v2d_max
      real*QUAD my_avzeta, my_ke,   my_pe
#ifdef SOLVE3D
     &        , my_ke3bc,  my_ke2b, my_kesrf
#endif
      integer i,j,k, nsubs, ierr, is,ie

#ifdef SOLVE3D
      real v2bc
# ifdef MAX_ADV_CFL
      real my_Cu_Adv, my_Cu_W, ciV, cx, cw
      integer my_i_cmax, my_j_cmax, my_k_cmax
# else
      real my_v3d_max, my_v3bc_max
# endif
#endif
#ifdef SUM_BY_PAIRS
      integer isize,itg, js,jsize,jtg
#endif
      integer, parameter :: max_check_line=128
      character(len=max_check_line) check_line
      character(len=18) tstring

      ierr=0
      my_v2d_max=0.
#ifdef SOLVE3D
# ifdef MAX_ADV_CFL
      my_Cu_Adv=0.  ;  my_Cu_W=0.
      my_i_cmax=0   ;  my_j_cmax=0  ;  my_k_cmax=0
# else
      my_v3d_max=0. ;  my_v3bc_max=0.
# endif
      do j=jstr,jend+1
        do i=istr,iend+1
          ub(i,j)=(Hz(i,j,N)+Hz(i-1,j,N))*u(i,j,N,nstp)
          vb(i,j)=(Hz(i,j,N)+Hz(i,j-1,N))*v(i,j,N,nstp)
        enddo
        do k=N-1,2,-1
          do i=istr,iend+1
            ub(i,j)=ub(i,j)+(Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k,nstp)
            vb(i,j)=vb(i,j)+(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,nstp)
          enddo
        enddo
        do i=istr,iend+1
          ub(i,j)=(ub(i,j)+(Hz(i,j,1)+Hz(i-1,j,1))*u(i,j,1,nstp))
     &         /(z_w(i,j,N)+z_w(i-1,j,N)-z_w(i,j,0)-z_w(i-1,j,0))
          vb(i,j)=(vb(i,j)+(Hz(i,j,1)+Hz(i,j-1,1))*v(i,j,1,nstp))
     &         /(z_w(i,j,N)+z_w(i,j-1,N)-z_w(i,j,0)-z_w(i,j-1,0))
        enddo
      enddo

      cff=g/rho0
      do j=jstr,jend
        do i=istr,iend
          v2=0.5*(ub(i,j)**2+ub(i+1,j)**2 +vb(i,j)**2+vb(i,j+1)**2)
          my_v2d_max=max(my_v2d_max, v2)

          ke(i,j)=0.
          pe(i,j)=0.5*g*z_w(i,j,N)*z_w(i,j,N)

          ke2b(i,j)=0.5*(z_w(i,j,N)-z_w(i,j,0))*v2
          ke3bc(i,j)=0.
          kesrf(i,j)=0.25*( u(i,j,N,nstp)**2 + u(i+1,j,N,nstp)**2
     &                     +v(i,j,N,nstp)**2 + v(i,j+1,N,nstp)**2)
        enddo
        do k=N,1,-1
          do i=istr,iend
            v2=0.5*( u(i,j,k,nstp)**2 + u(i+1,j,k,nstp)**2
     &              +v(i,j,k,nstp)**2 + v(i,j+1,k,nstp)**2)

            v2bc=0.5*( (u(i  ,j,k,nstp)-ub(i  ,j))**2
     &                +(u(i+1,j,k,nstp)-ub(i+1,j))**2
     &                +(v(i,j  ,k,nstp)-vb(i,j  ))**2
     &                +(v(i,j+1,k,nstp)-vb(i,j+1))**2)

! Find maximum Advective Courant number: in the code segment below it
! is defined in the finite-volume sense as the sum of fluxes directed
! OUTWARD from grid box Hz divided by its control volume. Essentially
! it is the fraction of water replaced within the grid box during one
! time step.   Note that under this definition, Cu=1 is the condition
! when flux-split first-order upstream advection scheme looses its
! positive definiteness property (hence stability).  Reported values
! are: "Cu_Adv" is full (tri-dimensional) Courant number; "i,j,k" are
! coordinates where its maximum occurs;  "Cu_W" is contribution into
! "Cu_Adv" due to vertical velocity.

# ifdef MAX_ADV_CFL
#  ifdef MASKING
            ciV=dt*rmask(i,j)*pm(i,j)*pn(i,j)/Hz(i,j,k)
#  else
            ciV=dt*pm(i,j)*pn(i,j)/Hz(i,j,k)
#  endif
#  ifdef VERT_EXPL_ONLY
            cw=ciV*( max(We(i,j,k), 0.) -min(We(i,j,k-1), 0.) )
#  else
            cw=ciV*( max(We(i,j,k)+Wi(i,j,k), 0.)
     &                      -min(We(i,j,k-1)+Wi(i,j,k-1), 0.) )
#  endif

            cx=cw+ciV*( max(FlxU(i+1,j,k), 0.) -min(FlxU(i,j,k), 0.)
     &                 +max(FlxV(i,j+1,k), 0.) -min(FlxV(i,j,k), 0.))

#  ifdef MAX_BY_VERT_CFL_ONLY
            if (cw > my_Cu_W) then
#  else
            if (cx > my_Cu_Adv) then
#  endif
              my_Cu_Adv=cx ; my_Cu_W=cw
              my_i_cmax=i  ; my_j_cmax=j ; my_k_cmax=k
            endif
# else
            my_v3d_max  = max(my_v3d_max,  v2  )
            my_v3bc_max = max(my_v3bc_max, v2bc)
# endif
            ke(i,j)=ke(i,j) + 0.5*v2*Hz(i,j,k)

            pe(i,j)=pe(i,j) + cff*Hz(i,j,k)
# ifdef SPLIT_EOS
     &      *(rho1(i,j,k)+qp1(i,j,k)*(z_w(i,j,N)-z_r(i,j,k)))
# else
     &                                            *rho(i,j,k)
# endif
     &                               *(z_r(i,j,k)-z_w(i,j,0))

            ke3bc(i,j)=ke3bc(i,j) + 0.5*v2bc*Hz(i,j,k)
          enddo
        enddo
        do i=istr,iend
# ifdef MASKING
          dA=rmask(i,j)/(pm(i,j)*pn(i,j))
# else
          dA=1./(pm(i,j)*pn(i,j))
# endif
          dVol(i,j) = dA * z_w(i,j,N)
          ke(i,j)   = dA * ke(i,j)
          pe(i,j)   = dA * pe(i,j)
          ke2b(i,j) = dA * ke2b(i,j)
          ke3bc(i,j)= dA * ke3bc(i,j)
          kesrf(i,j)= dA * kesrf(i,j)
        enddo
      enddo
#else
      cff=0.5*g
      do j=jstr,jend
        do i=istr,iend
          v2=0.5*( ubar(i,j,kstp)**2+ubar(i+1,j,kstp)**2
     &            +vbar(i,j,kstp)**2+vbar(i,j+1,kstp)**2)
          my_v2d_max=max(my_v2d_max, v2)
# ifdef MASKING
          dA=rmask(i,j)/(pm(i,j)*pn(i,j))
# else
          dA=1./(pm(i,j)*pn(i,j))
# endif
          dVol(i,j)=dA*zeta(i,j,kstp)
          ke(i,j)=dA*0.5*(zeta(i,j,kstp)+h(i,j))*v2
          pe(i,j)=cff*dVol(i,j)*zeta(i,j,kstp)
        enddo
      enddo
#endif /* SOLVE3D */
#ifdef SUM_BY_PAIRS

      isize=iend-istr                  ! Horizontal summation within
      jsize=jend-jstr                  ! the tile: compute partial sums
                                       ! using a split-directional
      do while (isize>0 .or. jsize>0)  ! reduction-by-pairs algorithm
        if (jsize>0) then              ! in order to avoid accumulation
          js=(jsize+1)/2-1             ! of roundoff errors.
          do j=0,js
            jtg=jstr+j
            do i=istr,istr+isize
               dVol(i,jtg) =  dVol(i,jtg+j) +  dVol(i,jtg+j+1)
                 ke(i,jtg) =    ke(i,jtg+j) +    ke(i,jtg+j+1)
                 pe(i,jtg) =    pe(i,jtg+j) +    pe(i,jtg+j+1)
# ifdef SOLVE3D
               ke2b(i,jtg) =  ke2b(i,jtg+j) +  ke2b(i,jtg+j+1)
              ke3bc(i,jtg) = ke3bc(i,jtg+j) + ke3bc(i,jtg+j+1)
              kesrf(i,jtg) = kesrf(i,jtg+j) + kesrf(i,jtg+j+1)
# endif
            enddo
          enddo
          if (2*js+1 < jsize) then     ! carry-over points which do not
            js=js+1                    ! have pair to be added together
            jtg=jstr+js
            do i=istr,istr+isize
               dVol(i,jtg) =  dVol(i,jtg+js)
                 ke(i,jtg) =    ke(i,jtg+js)
                 pe(i,jtg) =    pe(i,jtg+js)
# ifdef SOLVE3D
               ke2b(i,jtg) =  ke2b(i,jtg+js)
              ke3bc(i,jtg) = ke3bc(i,jtg+js)
              kesrf(i,jtg) = kesrf(i,jtg+js)
# endif
            enddo
          endif
          jsize=js
        endif

        if (isize>0) then
          is=(isize+1)/2-1
          do j=jstr,jstr+jsize
            do i=0,is
              itg=istr+i
               dVol(itg,j) =  dVol(itg+i,j) +  dVol(itg+i+1,j)
                 ke(itg,j) =    ke(itg+i,j) +    ke(itg+i+1,j)
                 pe(itg,j) =    pe(itg+i,j) +    pe(itg+i+1,j)
# ifdef SOLVE3D
               ke2b(itg,j) =  ke2b(itg+i,j) +  ke2b(itg+i+1,j)
              ke3bc(itg,j) = ke3bc(itg+i,j) + ke3bc(itg+i+1,j)
              kesrf(itg,j) = kesrf(itg+i,j) + kesrf(itg+i+1,j)
# endif
            enddo
          enddo
          if (2*is+1 < isize) then
            is=is+1
            itg=istr+is
            do j=jstr,jstr+jsize
               dVol(itg,j) =  dVol(itg+is,j)
                 ke(itg,j) =    ke(itg+is,j)
                 pe(itg,j) =    pe(itg+is,j)
# ifdef SOLVE3D
               ke2b(itg,j) =  ke2b(itg+is,j)
              ke3bc(itg,j) = ke3bc(itg+is,j)
              kesrf(itg,j) = kesrf(itg+is,j)
# endif
            enddo
          endif
          isize=is
        endif
      enddo  !<-- while()
#else
! Alternative: a two-stage (one for each direction) straightforward
! summation algorithm as an effort to mitigate roundoff errors.

      do j=jstr+1,jend
        do i=istr,iend
          dVol(i,jstr)  = dVol(i,jstr)  + dVol(i,j)
          ke(i,jstr)    = ke(i,jstr)    + ke(i,j)
          pe(i,jstr)    = pe(i,jstr)    + pe(i,j)
# ifdef SOLVE3D
          ke2b(i,jstr)  = ke2b(i,jstr)  + ke2b(i,j)
          ke3bc(i,jstr) = ke3bc(i,jstr) + ke3bc(i,j)
          kesrf(i,jstr) = kesrf(i,jstr) + kesrf(i,j)
# endif
        enddo
      enddo
      do i=istr+1,iend
        dVol(istr,jstr)  = dVol(istr,jstr)  + dVol(i,jstr)
        ke(istr,jstr)    = ke(istr,jstr)    + ke(i,jstr)
        pe(istr,jstr)    = pe(istr,jstr)    + pe(i,jstr)
# ifdef SOLVE3D
        ke2b(istr,jstr)  = ke2b(istr,jstr)  + ke2b(i,jstr)
        ke3bc(istr,jstr) = ke3bc(istr,jstr) + ke3bc(i,jstr)
        kesrf(istr,jstr) = kesrf(istr,jstr) + kesrf(i,jstr)
# endif
      enddo
#endif /* SUM_BY_PAIRS */

      my_avzeta=dVol(istr,jstr)
      my_ke=ke(istr,jstr)
      my_pe=pe(istr,jstr)
#ifdef SOLVE3D
      my_ke2b=ke2b(istr,jstr)
      my_ke3bc=ke3bc(istr,jstr)
      my_kesrf=kesrf(istr,jstr)
#endif
      if (SINGLE_TILE_MODE) then
        nsubs=1
      else
        nsubs=NSUB_X*NSUB_E
      endif

! Perform global summation: whoever gets first to the critical region
! resets global sums before global summation starts; after the global
! summation is completed, thread, which is the last one to enter the
! critical region, finalizes the computation of diagnostics and
! prints them out.

C$OMP CRITICAL (diag_cr_rgn)
      if (tile_count == 0) then              ! Initialize global sums
        avzeta=my_avzeta                     ! for the summation among
        avke=my_ke                           ! the threads within the
        avpe=my_pe                           ! shared-memory group
        v2d_max=my_v2d_max                   ! (within MPI process).
#ifdef SOLVE3D
        avke2b=my_ke2b
        avke3bc=my_ke3bc
        avkesrf=my_kesrf
# ifdef MAX_ADV_CFL
        Cu_Adv=my_Cu_Adv
        Cu_W=my_Cu_W
        i_cx_max=my_i_cmax
        j_cx_max=my_j_cmax
        k_cx_max=my_k_cmax
# else
        v3d_max=my_v3d_max
        v3bc_max=my_v3bc_max
# endif
#endif
      else                                   ! Perform global summation
        avzeta=avzeta+my_avzeta              ! among the threads within
        avke=avke+my_ke                      ! MPI process.
        avpe=avpe+my_pe
        v2d_max=max(v2d_max, my_v2d_max)
#ifdef SOLVE3D
        avke2b=avke2b+my_ke2b
        avke3bc=avke3bc+my_ke3bc
        avkesrf=avkesrf+my_kesrf
# ifdef MAX_ADV_CFL
#  ifdef MAX_BY_VERT_CFL_ONLY
        if (my_Cu_W > Cu_W) then
#  else
        if (my_Cu_Adv > Cu_Adv) then
#  endif
          Cu_Adv=my_Cu_Adv
          Cu_W=my_Cu_W
          i_cx_max=my_i_cmax
          j_cx_max=my_j_cmax
          k_cx_max=my_k_cmax
        endif
# else
        v3d_max=max(v3d_max, my_v3d_max)
        v3bc_max=max(v3bc_max,my_v3bc_max)
# endif
#endif
      endif
                                             ! This counter identifies
      tile_count=tile_count+1                ! the last thread entering
                                             ! critical region, which
      if (tile_count == nsubs) then          ! ever it happens to be,
        tile_count=0                         ! not necessarily master.
#ifdef MPI
# ifdef SOLVE3D
#  ifdef MAX_ADV_CFL
        i_cx_max=i_cx_max + iSW_corn         ! convert to absolute
        j_cx_max=j_cx_max + jSW_corn         ! indices
#  endif
# endif
        if (NNODES>1) then                   ! Perform global summation
          size=NNODES                        ! among all MPI processes
          do while (size>1)
            step=(size+1)/2
            if (mynode>=step .and. mynode<size) then
              buff(1)=may_day_flag
              buff(2)=avzeta
              buff(3)=avke                   ! This is essentially an
              buff(4)=avpe                   ! MPI_Reduce operation,
              buff(5)=avke2b                 ! although with a rather
              buff(6)=avke3bc                ! complicated function.
              buff(7)=avkesrf
              buff(8)=v2d_max
# ifdef SOLVE3D
#  ifdef MAX_ADV_CFL
              buff(9)=Cu_Adv
              buff(10)=Cu_W
              buff(11)=i_cx_max
              buff(12)=j_cx_max
              buff(13)=k_cx_max
#  else
              buff(9)=v3d_max
              buff(10)=v3bc_max
#  endif
# endif

c**     write(*,*) 'sending  ',mynode,'->',mynode-step,' size=',size
c**     write(*,*) 'sending  avke =', buff(3)/volume

              itag=mynode+400
              call MPI_Send (buff, 14, MPI_BUFF_TYPE, mynode-step,
     &                       itag, ocean_grid_comm,          ierr)
            elseif (mynode < size-step) then
              itag=mynode+400+step
              call MPI_Recv (buff, 14, MPI_BUFF_TYPE, mynode+step,
     &                       itag, ocean_grid_comm,  status, ierr)

c**     write(*,*) 'received ',mynode,'<-',mynode+step,' size=',size
c**     write(*,*) 'received avke =', buff(3)/volume

              i=buff(1)+0.0001
              if (i>0 .and. may_day_flag==0) may_day_flag=i
              avzeta=avzeta+buff(2)
              avke=avke+buff(3)              ! Here "v2" is used as
              avpe=avpe+buff(4)              ! a scratch variable to
              avke2b=avke2b+buff(5)          ! avoid type conflict
              avke3bc=avke3bc+buff(6)        ! in the case if "buff"
              avkesrf=avkesrf+buff(7)        ! is 16-Byte, while
              v2=buff(8)                     ! "v2", "v2d_max", etc
              v2d_max=max(v2d_max, v2)       ! are only 8-byte long.
# ifdef SOLVE3D
#  ifdef MAX_ADV_CFL
#   ifdef MAX_BY_VERT_CFL_ONLY
              if (buff(10) > Cu_W) then
#   else
              if (buff(9) > Cu_Adv) then
#   endif
                Cu_Adv=buff(9)
                Cu_W=buff(10)
                i_cx_max=buff(11) +0.0001
                j_cx_max=buff(12) +0.0001
                k_cx_max=buff(13) +0.0001
              endif
#  else
              v2=buff(9)
              v3d_max=max(v3d_max, v2)
              v2=buff(19)
              v3bc_max=max(v3bc_max, v2)
#  endif
# endif
            endif
            size=step
          enddo  !<-- while (size>1)
        endif !<-- NNODES>1

        if (mynode==0) then
#endif
          avke=avke/(volume+avzeta)            ! Compute and print
          avpe=avpe/(volume+avzeta)            ! global diagnostics
          avke2b=avke2b/(volume+avzeta)        ! (the last thread of
          avke3bc=avke3bc/(volume+avzeta)      ! master MPI process
          avkesrf=avkesrf/area                 ! only).
          v2d_max=sqrt(v2d_max)
          avzeta=avzeta/area
#ifdef SOLVE3D
# ifndef MAX_ADV_CFL
          v3d_max=sqrt(v3d_max)
          v3bc_max=sqrt(v3bc_max)
# endif
#endif


#ifdef SOLVE3D
# if defined MAX_ADV_CFL
          if (first_time==0) then
            first_time=1
            if (code_check) then
              write(*,3) 'STEP',  'KINETIC_ENRG',                  ! Special for code checks
     &                'BAROTR_KE', 'MAX_ADV_CFL', 'MAX_VERT_CFL'
            else
              write(*,1) 'STEP', 'time[DAYS]',  'KINETIC_ENRG',    !  Normal log file output
     &                'BAROTR_KE', 'MAX_ADV_CFL', 'MAX_VERT_CFL',
     &                                     'i_cx', 'j_cx', 'k_c'
            endif
          endif
          if (code_check) then
            write(*,4) iic, avke,  avke2b, Cu_Adv, Cu_W
          else
            write(check_line,2,iostat=ierr)  avke,  avke2b,  Cu_Adv,
     &                        Cu_W, i_cx_max,  j_cx_max,  k_cx_max
          endif

  1   format(/1x,A,2x,A,1x,A,5x,A,8x,A,2(5x,A),3x,A,2(1x,A))
  2   format(ES18.11, ES17.10, 2F16.12, 2I7,I5, I3)
  3   format(/2x,A,1x,A,5x,A,8x,A,2(5x,A),3x,A,2(1x,A))
  4   format(I2,1x,4ES23.16)

# elif defined REPORT_MAX_UV && defined REPORT_POTEN_ENRG

          if (first_time==0) then
            first_time=1
            write(*,1) 'STEP', 'time[DAYS]', 'KINETIC_ENRG',
     &                   'BAROTR_KE', 'MAX_U3D', 'MAX_UBAR',
     &                   'MAX_UBCL', 'POTEN_ENRG'
C$   &                                                 , 'trd'
          endif
          write(check_line,2,iostat=ierr)  avke,    avke2b,
     &                       v3d_max, v2d_max, v3bc_max, avpe
C$   &                                            , proc(2)
  1       format(/1x,A,2x,A,1x,A,5x,A,6x,A,7x,A,6x,A,7x,A,3x,A)
  2       format(ES18.11, ES16.9, 3(ES15.8), ES16.8,I3)

# elif defined REPORT_MAX_UV
          if (first_time==0) then
            first_time=1
            write(*,1) 'STEP', 'time[DAYS]', 'KINETIC_ENRG',
     &                 'BAROTR_ENRG', 'MAX_U3D', 'MAX_UBAR',
     &                 'MAX_UBCL', 'FREE_SURFACE'
C$   &                                                 , 'trd'
          endif
          write(check_line,2,iostat=ierr)  avke, avke2b, v3d_max,
     &                                  v2d_max, v3bc_max, avzeta
C$   &                                                 , proc(2)
  1       format(/1x,A,2x,A,1x,A,5x,A,4x,A,7x,A,6x,A,7x,A,1x,A)
  2       format(ES18.11,ES16.9,3(ES15.8),ES16.8,I3)
# else
          if (first_time == 0) then
            first_time=1
            write(*,'(1x,A,2x,A,1x,A,5x,A,4x,A,4x,A,1x,A)')
     &              'STEP',       'time[DAYS]', 'KINETIC_ENRG',
     &              'SURFACE_KE', 'POTEN_ENRG', 'FREE_SURFACE'
C$   &                                                  ,'trd'
          endif
          write(check_line,'(ES18.11,3(ES15.8),I3)',iostat=ierr)
     &                               avke, avkesrf, avpe, avzeta
C$   &                                                 , proc(2)
# endif
#else
          if (first_time == 0) then
            first_time=1
            write(*,'(1x,A,2x,A,1x,A,5x,A,4x,A,9x,A,1x,A)')
     &                   'STEP',    'time[DAYS]', 'KINETIC_ENRG',
     &                   'POTEN_ENRG', 'KE+PE',   'FREE_SURFACE'
C$   &                                                 , 'trd'
          endif
          write(check_line,'(ES18.11,3(ES15.8),I3)',iostat=ierr)
     &                             avke, avpe, avke+avpe, avzeta
C$   &                                                 , proc(2)
#endif /* SOLVE3D */
          if (.not.code_check) then
          ie=max_check_line
          do while (check_line(ie:ie)==' ' .and. ie>0)
            ie=ie-1
          enddo
          is=1
          do while (check_line(is:is)==' ' .and. is<ie)
            is=is+1
          enddo                            ! Suppress floating-point
          i=is-1                           ! "E"s in order to shorten
          do while (i<ie)                  ! printed diagnostic line
            i=i+1
            if (check_line(i:i)=='E' .or. check_line(i:i)=='e') then
              check_line(i:ie-1)=check_line(i+1:ie)
              check_line(ie:ie)=' '
              ie=ie-1

! Appearance of any non-digit character other than blank,"+", "-", or
! "." indicates floating point exception condition (usually caused by
! blowup).    Force writing history file at every time step and set
! may_day_flag to signal termination.  The code behavior in this
! situation is to write the last state of all variables into history
! file and quit.

            elseif (check_line(i:i)<'0' .or. check_line(i:i)>'9') then
              if ( check_line(i:i) /= ' ' .and.
     &             check_line(i:i) /= '+'  .and.
     &             check_line(i:i) /= '-'  .and.
     &             check_line(i:i) /= '.') then
                if (may_day_flag == 0) then
                  write(*,*) '---> setting may day flag'
                  may_day_flag=1    !--> BLOW UP
                endif
              endif
            endif
          enddo  ! <-- while(..)

          write(tstring,'(F18.8)') tdays
          i=1
          do while (i<18 .and. tstring(i:i)==' ')
            i=i+1
          enddo
          write(*,'(I7,2(1x,A))') iic, tstring(i:i+8),
     &                                check_line(is:ie)
          endif ! .not. code_check

          call flush(6)


#ifdef MPI
        endif    ! <-- mynode==0

        buff(1)=may_day_flag

        call MPI_Bcast(buff,2, MPI_BUFF_TYPE, 0,ocean_grid_comm,ierr)

        may_day_flag=buff(1) +0.0001
#endif
      endif  ! <-- tile_count==nsubs
C$OMP END CRITICAL (diag_cr_rgn)
      diag_sync=.true.     !<-- to signal that "diag" has been called
      end

