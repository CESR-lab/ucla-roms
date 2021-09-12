#include "cppdefs.h"
#ifdef SOLVE3D
# define UPSTREAM_UV

      subroutine step3d_uv1(tile)
      use param
      implicit none
      integer tile

# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call step3d_uv1_tile( istr,iend,jstr,jend,    A3d(1,1),A3d(1,2),
     &                            A2d(1,1),A2d(1,2),A2d(1,3),A2d(1,4),
     &                            A2d(1,1),A2d(1,2),A2d(1,3),A2d(1,4),
     &                                              A2d(1,5),A2d(1,6)
# ifdef NHMG
     &                                             ,A3d(1,5)
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

      use param
# ifdef NHMG
      use nhmg
      use mg_grids
      use mg_tictoc, only : tic, toc
# endif
# ifdef WEC
      use wec_frc
# endif
      use diagnostics
      use surf_flux, only: sustr, svstr

      implicit none
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
      real  div,dmax
!     integer :: ierr !,imx,jmx,kmx
#endif

#  include "ocean2d.h"
# include "grid.h"
# include "ocean3d.h"
# include "coupling.h"
# include "mixing.h"
# include "scalars.h"
# if defined ADV_ISONEUTRAL && defined EXCHANGE
#  include "eos_vars.h"
# endif

! Taken from rhs3d35S.F for WEC, used in compute_horiz_rhs_uv_terms.h
! but needs to be declared here in declarations section not in
! the compute_horiz_rhs_uv_terms.h file else wont compile.
# ifdef WEC
      real cff1
# endif

# include "compute_auxiliary_bounds.h"

# ifdef DIAGNOSTICS
      if (diag_uv .and. calc_diag) then                                                      ! called here not in prsgrd as prsgrd is called twice from main.
        Mdiag(istrU:iend,jstr :jend,:,1,iprsgr)=ru(istrU:iend,jstr :jend,:)    ! ru is scratch size not global size
        Mdiag(istr :iend,jstrV:jend,:,2,iprsgr)=rv(istr :iend,jstrV:jend,:)    ! interior points not bry
        diag_icori = .true.  ! triggers diag calculation in compute_horiz_uv... as also called from pre_step3d, but only want this time.
      endif
      if (prs_flx) pg_cor = .false. ! ensure prs_flx not calculated in predictor step of prsgrd.
#  ifdef NHMG
      if (calc_diag) then
        Wdiag(istr:iend,jstr:jend,1:N,iwprsgr) =rw(istr :iend,jstr :jend,1:N)  ! Note 1:N
      endif
#  endif
# endif

# ifdef NHMG

#  define FX UFx
#  define FE VFe

      do k=1,N
#  include "compute_horiz_rhs_w_terms.h"
      enddo

#  undef FE
#  undef FX

      do j=jstr,jend

#  include "compute_vert_rhs_w_terms.h"
        do i=istr,iend
          DC(i,0)=dt*pm(i,j)*pn(i,j)
        enddo
        do k=1,N           !! w is dz*w after this
          do i=istr,iend
            w(i,j,k,nnew)=w(i,j,k,nnew) +DC(i,0)*rw(i,j,k)
          enddo
        enddo

      enddo
#  ifdef DIAGNOSTICS
      if (calc_diag) call set_diags_w_at_uv1( istr, iend, jstr, jend, DC, rw )
#  endif
# endif  /* NHMG */

      do k=1,N
# include "compute_horiz_rhs_uv_terms.h"
      enddo

# ifdef WEC
      call compute_wec_rhs_uv_terms(istr,iend,jstr,jend,ru,rv,wrk1,wrk2)
# endif

# ifdef DIAGNOSTICS
        if (diag_uv .and. calc_diag) then
          call set_diags_u_4th_adv( istrU, istr, iend, jstr, jend )
          call set_diags_v_4th_adv( jstrV, istr, iend, jstr, jend )
        endif
# endif

      do j=jstr,jend

# include "compute_vert_rhs_uv_terms.h"

! Apply wave bottom streaming stress as a body force in the bottom-most
! grid cells only (if !defined BODY_FRICTION)
!! JM move this bit to inside compute_wec_rhs_uv_terms...
# ifdef WEC
#  if defined WAVE_FRICTION && !defined BODY_FRICTION
        !! JM add this somehow like a bustr in a manner like below
        !! or move into compute_wec_rhs_uv.
        do i=istrU,iend
          ru(i,j,1) =ru(i,j,1) + dm_u(i,j)*dn_u(i,j)*frc2dx(i,j)
        enddo
        do i=istr,iend
          rv(i,j,1) =rv(i,j,1) + dm_v(i,j)*dn_v(i,j)*frc2de(i,j)
        enddo
#  endif
# endif 

        do i=istrU,iend
          DC(i,0)=dt*0.25*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))

          u(i,j,N,nnew)=u(i,j,N,nnew) +DC(i,0)*ru(i,j,N)
     &                                     +dt*sustr(i,j)
#  if defined WEC && defined SURFACE_BREAK
        !! JM we will add brk2dx to sustr
     &                                    + dt * brk2dx(i,j)
#  endif
          u(i,j,1,nnew)=u(i,j,1,nnew)+DC(i,0)*ru(i,j,1)

          rufrc(i,j)=ru(i,j,1) +ru(i,j,N) +dn_u(i,j)*dm_u(i,j)*(
     &                                                 sustr(i,j)
#   if defined WEC && defined SURFACE_BREAK
     &                                      + brk2dx(i,j)
#   endif
     &                 -0.5*(r_D(i-1,j)+r_D(i,j))*u(i,j,1,nstp) )
        enddo
        do k=2,N-1
          do i=istrU,iend
            u(i,j,k,nnew)=u(i,j,k,nnew)+DC(i,0)*ru(i,j,k)

            rufrc(i,j)=rufrc(i,j) +ru(i,j,k)
          enddo
        enddo
# ifdef DIAGNOSTICS
        if (diag_uv .and. calc_diag)
     &    call set_diags_uv_dissipation(istr, iend, jstr, jend, istrU, j, 1, DC, ru, sustr, u)
# endif
        if (j >= jstrV) then
          do i=istr,iend
            DC(i,0)=dt*0.25*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))

            v(i,j,N,nnew)=v(i,j,N,nnew) +DC(i,0)*rv(i,j,N)
     &                                       +dt*svstr(i,j)
#  ifdef WEC
     &                                       +dt*brk2de(i,j)
#  endif
            v(i,j,1,nnew)=v(i,j,1,nnew) +DC(i,0)*rv(i,j,1)

            rvfrc(i,j)=rv(i,j,1) +rv(i,j,N) +dm_v(i,j)*dn_v(i,j)*(
     &                                                   svstr(i,j)
# ifdef WEC
     &                                               +brk2de(i,j)
# endif
     &                   -0.5*(r_D(i,j-1)+r_D(i,j))*v(i,j,1,nstp) )
          enddo
          do k=2,N-1
            do i=istr,iend
              v(i,j,k,nnew)=v(i,j,k,nnew)+DC(i,0)*rv(i,j,k)

              rvfrc(i,j)=rvfrc(i,j) +rv(i,j,k)
            enddo
          enddo
# ifdef DIAGNOSTICS
          if (diag_uv .and. calc_diag)
     &      call set_diags_uv_dissipation(istr, iend, jstr, jend, istr, j, 2, DC, rv, svstr, v)
# endif
        endif  !<-- if j>=jstrV

      enddo  !<-- j
      !! velocities are in the form:  u*dz, v*dz, w*dz

# if defined NHMG

#  ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do k=1,N
          do j=jstr,jend
            u(istr,j,k,nnew) = u(istr,j,k,nnew)*
     &         0.5*(Hz(istr-1,j,k)+Hz(istr,j,k))
          enddo
        enddo
      endif
      if (EASTERN_EDGE) then
        do k=1,N
          do j=jstr,jend
            u(iend+1,j,k,nnew) = u(iend+1,j,k,nnew)*
     &         0.5*(Hz(iend+1,j,k)+Hz(iend,j,k))
          enddo
        enddo
      endif
#  endif
#  ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do k=1,N
          do i=istr,iend
            v(i,jstr,k,nnew) = v(i,jstr,k,nnew)*
     &         0.5*(Hz(i,jstr-1,k)+Hz(i,jstr,k))
          enddo
        enddo
      endif
      if (NORTHERN_EDGE) then
        do k=1,N
          do i=istr,iend
            v(i,jend+1,k,nnew) = v(i,jend+1,k,nnew)*
     &         0.5*(Hz(i,jend+1,k)+Hz(i,jend,k))
          enddo
        enddo
      endif
#  endif


      ! Translate to fluxes, different than in pre_step, we have dz*u
      ! instead of u

      do k=1,N
        do j=jstr,jend
          do i=istr,iend+1
            u(i,j,k,nnew) = u(i,j,k,nnew)*dn_u(i,j)
          enddo
        enddo
        do j=jstr,jend+1
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

      iprec1 = 1+mod(iprec1,2)
      iprec2 = 1+mod(iprec2,2)

# endif /* NHMG */


# undef r_D
# if defined ADV_ISONEUTRAL && defined EXCHANGE
       !!!! JM not clear what this does !!!
      call exchange_2_tile (istr,iend,jstr,jend, dRdx,N, dRde,N)
# endif

      end

      subroutine check_step_uv1_switches(ierr)

! This code must be placed here rather than in a separate file in order
! to be exposed to the relevant CPP-macros defined locally in this file
! above.  It does not affect any model results, other than signature in
! global attribute "CPPS" in output netCDF files.

      use param
      implicit none
      integer ierr, is,ie, lenstr
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
