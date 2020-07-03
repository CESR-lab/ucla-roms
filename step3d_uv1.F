#include "cppdefs.h"
#ifdef SOLVE3D
# define UPSTREAM_UV

      subroutine step3d_uv1(tile)
      implicit none
      integer tile
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call step3d_uv1_tile( istr,iend,jstr,jend,  A3d(1,1), A3d(1,2),
     &                            A2d(1,1),A2d(1,2),A2d(1,3),A2d(1,4),
     &                            A2d(1,1),A2d(1,2),A2d(1,3),A2d(1,4),
     &                                             A2d(1,5),A2d(1,6)
# ifdef NHMG
     &                                            ,A3d(1,5)
# endif
     &                                                              )
      end

      subroutine step3d_uv1_tile( istr,iend,jstr,jend, ru,rv,
     &                                            WC,FC,CF,DC,
     &                             UFx,UFe,VFx,VFe, wrk1,wrk2 
# ifdef NHMG
     &                                            ,rw
# endif
     &                                                              )
#ifdef NHMG
      use nhmg, only : nhmg_solve,halo
      use mg_grids
      use mg_tictoc, only : tic, toc
#endif
      implicit none
# include "param.h"
      integer istr,iend,jstr,jend, imin,imax,jmin,jmax, i,j,k
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,N) :: ru,rv
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N) :: WC,FC,CF,DC
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: UFx,UFe,VFx,VFe,
     &                                                   wrk1,wrk2
      real cff
      real, parameter ::  delta=0.1666666666666667 ! delta=0.125
# ifdef UPSTREAM_UV
     &                  , gamma=0.3333333333333333 ! gamma=0.25
# endif
# ifdef NHMG
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,0:N) :: rw
      real  Flxw,Uflxw,Vflxw
      real  ub,ut
!     real  div,dmax
!     integer :: ierr !,imx,jmx,kmx
#endif

# if defined NHMG
#  include "ocean2d.h"
# endif
# include "grid.h"
# include "ocean3d.h"
# include "coupling.h"
# include "mixing.h"
# include "forces.h"
# include "scalars.h"
# if defined ADV_ISONEUTRAL && defined EXCHANGE
#  include "eos_vars.h"
# endif

# include "compute_auxiliary_bounds.h"


# ifdef NHMG

# define FX UFx
# define FE VFe

      do k=1,N
# include "compute_horiz_rhs_w_terms.h"
      enddo

# undef FE
# undef FX


      do j=jstr,jend

# include "compute_vert_rhs_w_terms.h"
!  what is the dimension of w at this point?
!  Since rw is m4/s2, this would make w in dz*w, just like u and v
        do i=istr,iend
          DC(i,0)=dt*pm(i,j)*pn(i,j)
        enddo
        do k=1,N
          do i=istr,iend
            w(i,j,k,nnew)=w(i,j,k,nnew) +DC(i,0)*rw(i,j,k)
          enddo
        enddo

      enddo
# endif  /* NHMG */

      do k=1,N
# include "compute_horiz_rhs_uv_terms.h"
      enddo

      do j=jstr,jend

# include "compute_vert_rhs_uv_terms.h"

        do i=istrU,iend
          DC(i,0)=dt*0.25*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))

          u(i,j,N,nnew)=u(i,j,N,nnew) +DC(i,0)*ru(i,j,N)
     &                                     +dt*sustr(i,j)
          u(i,j,1,nnew)=u(i,j,1,nnew)+DC(i,0)*ru(i,j,1)

          rufrc(i,j)=ru(i,j,1) +ru(i,j,N) +dn_u(i,j)*dm_u(i,j)*(
     &                                                 sustr(i,j)
     &                 -0.5*(r_D(i-1,j)+r_D(i,j))*u(i,j,1,nstp) )
        enddo
        do k=2,N-1
          do i=istrU,iend
            u(i,j,k,nnew)=u(i,j,k,nnew)+DC(i,0)*ru(i,j,k)

            rufrc(i,j)=rufrc(i,j) +ru(i,j,k)
          enddo
        enddo
        if (j >= jstrV) then
          do i=istr,iend
            DC(i,0)=dt*0.25*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))

            v(i,j,N,nnew)=v(i,j,N,nnew) +DC(i,0)*rv(i,j,N)
     &                                       +dt*svstr(i,j)
            v(i,j,1,nnew)=v(i,j,1,nnew) +DC(i,0)*rv(i,j,1)
     &                    -0.5*(r_D(i,j-1)+r_D(i,j))*v(i,j,1,nstp)
            rvfrc(i,j)=rv(i,j,1) +rv(i,j,N) +dm_v(i,j)*dn_v(i,j)*(
     &                                                   svstr(i,j)
     &                   -0.5*(r_D(i,j-1)+r_D(i,j))*v(i,j,1,nstp) )
          enddo
          do k=2,N-1
            do i=istr,iend
              v(i,j,k,nnew)=v(i,j,k,nnew)+DC(i,0)*rv(i,j,k)

              rvfrc(i,j)=rvfrc(i,j) +rv(i,j,k)
            enddo
          enddo
        endif

      enddo  !<-- j
      !! Looks like we have u*dz, v*dz, w*dz at the moment

# if defined NHMG
      ! Translate to fluxes, different than in pre_step, we have dz*u
      ! instead of u

      do k=1,N
        do j=jstr,jend
          do i=istr,iend
            u(i,j,k,nnew) = u(i,j,k,nnew)*dn_u(i,j)
          enddo
        enddo
        do j=jstr,jend
          do i=istr,iend
            v(i,j,k,nnew) = v(i,j,k,nnew)*dm_v(i,j)
          enddo
        enddo
        if (k==N) then
          do j=Jstr,Jend
            do i=istr,iend
              w(i,j,k,nnew) = w(i,j,k,nnew)*dm_r(i,j)*dn_r(i,j)
     &            /(0.5*(            Hz(i,j,k)))
            enddo
          enddo
        else
          do j=Jstr,Jend
            do i=istr,iend
              w(i,j,k,nnew) = w(i,j,k,nnew)*dm_r(i,j)*dn_r(i,j)
     &          /(0.5*(Hz(i,j,k+1)+Hz(i,j,k)))
            enddo
          enddo
        endif
      enddo

#  ifdef EXCHANGE
      !! Todo: we should be able to avoid this. We only need u and v at the boundaries
      !! Todo: computational ranges for u and v need to be extended to Iend+1, Jend+1
      call exchange_2_tile (istr,iend,jstr,jend,
     &                         u(START_2D_ARRAY,1,nnew),N,
     &                         v(START_2D_ARRAY,1,nnew),N)
#  endif

      !! only do this for neumann nhmg conditions
      ! Compute a 'barotropic' correction to w such that it matches
      ! the current ubar,vbar divergence (see Molemaker et al., 2018)
      if (surface_neumann) then

        nh_ubar = 0.
        nh_vbar = 0.
        do k=1,N
          do j=Jstr,Jend
            do i=Istr,Iend+1
              nh_ubar(i,j) = nh_ubar(i,j) + u(i,j,k,nnew)
            enddo
          enddo
          do j=Jstr,Jend+1
            do i=Istr,Iend
              nh_vbar(i,j) = nh_vbar(i,j) + v(i,j,k,nnew)
            enddo
          enddo
        enddo
        do j=Jstr,Jend
          do i=Istr,Iend
            nh_wcor(i,j) = w(i,j,N,nnew) +
     &         (nh_ubar(i+1,j)-nh_ubar(i,j)+nh_vbar(i,j+1)-nh_vbar(i,j))
          enddo
        enddo
        do k=1,N
          do j=Jstr,Jend
            do i=Istr,Iend
              w(i,j,k,nnew) = w(i,j,k,nnew) - nh_wcor(i,j)
     &                       * (z_w(i,j,k)-z_w(i,j,0))
     &                       / (z_w(i,j,N)-z_w(i,j,0))
            enddo
          enddo
        enddo
      endif

!     if (mynode==0) then
!       print *,'----------------------'
!       ub = u(31,10,1,nnew)*pn_u(31,10)/(0.5*(Hz(31,10,1)+Hz(30,10,1)))
!       ut = u(31,10,n,nnew)*pn_u(31,10)/(0.5*(Hz(31,10,n)+Hz(30,10,n)))
!       print 99,ub,ut,ut/ub
!       ub = u(31,10,n,nnew)*pn_u(31,10)/(0.5*(Hz(31,10,n)+Hz(30,10,n)))
!       ut = w(1,10,n,nnew)*pn(10,1)*pm(10,1)
!       print 98,ub,ut,ut/ub
!     endif
! 99  format( 'uv1 befor ub,ut,ratio: ',3e13.4)
! 98  format( 'uv1 befor ut,wt,ratio: ',3e13.4)

      call nhmg_solve(Lm,Mm,N,halo,padd_X,padd_E,
     &           u(:,:,:,nnew),v(:,:,:,nnew),w(:,:,:,nnew) )

      ! At this point u,v,w are fluxes with the dimension of m^3/s
      do k=1,N
        do j=jstr,jend
          do i=Istr,Iend
            u(i,j,k,nnew) = (u(i,j,k,nnew) + grid(1)%du(k,j,i) )*pn_u(i,j)

            rufrc(i,j) = rufrc(i,j) + dn_u(i,j)*grid(1)%du(k,j,i)/dt 
            nhdu(i,j,k,iprec1) = grid(1)%du(k,j,i)/dt !stored as a trend
          enddo
        enddo

        do j=jstr,jend
          do i=Istr,Iend
            v(i,j,k,nnew) = (v(i,j,k,nnew) + grid(1)%dv(k,j,i) )*pm_v(i,j)

            rvfrc(i,j) = rvfrc(i,j) + dm_v(i,j)*grid(1)%dv(k,j,i)/dt 
            nhdv(i,j,k,iprec1) = grid(1)%dv(k,j,i)/dt !stored as a trend
          enddo
        enddo

        if (k==N) then
          do j=jstr,jend
            do i=Istr,Iend
              w(i,j,k,nnew) = (w(i,j,k,nnew) + grid(1)%dw(k+1,j,i) )
     &              * (pm(i,j)*pn(i,j))*(0.5*(            Hz(i,j,k)))

              nhdw(i,j,k,iprec1) = grid(1)%dw(k+1,j,i)/dt !stored as a trend
            enddo
          enddo
        else
          do j=jstr,jend
            do i=Istr,Iend
              w(i,j,k,nnew) = (w(i,j,k,nnew) + grid(1)%dw(k+1,j,i) )
     &              * (pm(i,j)*pn(i,j))*(0.5*(Hz(i,j,k+1)+Hz(i,j,k)))

              nhdw(i,j,k,iprec1) = grid(1)%dw(k+1,j,i)/dt !stored as a trend
            enddo
          enddo
        endif
      enddo

# endif
      !! velocities are now in m2/s (dz*u)

!     if (mynode==0) then
!       print 93, ru(31,10,1)*dt*0.5*(pm(31,10)+pm(30,10))
!    &              *pn_u(31,10)/(0.5*(Hz(31,10,1)+Hz(30,10,1)))
!    &           ,grid(1)%du(1,10,31)
!    &              *pn_u(31,10)/(0.5*(Hz(31,10,1)+Hz(30,10,1)))
!       ub = u(31,10,1,nnew)/(0.5*(Hz(31,10,1)+Hz(30,10,1)))
!       ut = u(31,10,n,nnew)/(0.5*(Hz(31,10,n)+Hz(30,10,n)))
!       print 97, ub,ut,ut/ub
!       ub = u(31,10,n,nnew)/(0.5*(Hz(31,10,n)+Hz(30,10,n)))
!       ut = w(1,10,n,nnew)/(0.5*Hz(1,10,n))
!       print 96,ub,ut,ut/ub
!       print 92, grid(1)%dw(n,10,1)/dt,dt
!       write(20,91), grid(1)%dw(n,10,1)/dt

!     endif
! 93  format( 'uv1 after ru,du: ',2e13.4)
! 97  format( 'uv1 after ub,ut,ratio: ',3e13.4)
! 96  format( 'uv1 after ut,wt,ratio: ',3e13.4)
! 92  format( 'uv1 after rw, dt: ',2e15.6)
! 91  format(2e15.6)

# undef r_D
# if defined ADV_ISONEUTRAL && defined EXCHANGE
      call exchange_2_tile (istr,iend,jstr,jend, dRdx,N, dRde,N)
# endif
      end

      subroutine check_step_uv1_switches(ierr)

! This code must be placed here rather than in a separate file in order
! to be exposed to the relevant CPP-macros defined locally in this file
! above.  It does not affect any model results, other than signature in
! global attribute "CPPS" in output netCDF files.

      implicit none
      integer ierr, is,ie, lenstr
# include "param.h"
# include "strings.h"
      ie=lenstr(cpps)
      is=ie+2 ; ie=is+13
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='<step3d_uv1.F>'

# include "track_advec_switches.h"

      return
  99  mpi_master_only write(*,'(/1x,2A/12x,A/)')      '### ERROR: ',
     &  'Insufficient length of string "cpps" in file "strings.h".',
     &        'Increase parameter "max_opt_size" it and recompile.'
      ierr=ierr+1
      end

#else
      subroutine step3d_uv1_empty
      end
#endif
