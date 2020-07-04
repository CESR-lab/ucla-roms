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
     &                                             A2d(1,5),A2d(1,6))
      end

      subroutine step3d_uv1_tile( istr,iend,jstr,jend, ru,rv,
     &                                            WC,FC,CF,DC,
     &                             UFx,UFe,VFx,VFe, wrk1,wrk2)

# ifdef WEC
      use wec_frc
# endif

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
# include "grid.h"
# include "ocean3d.h"
# include "coupling.h"
# include "mixing.h"
# include "forces.h"
# include "scalars.h"
# if defined ADV_ISONEUTRAL && defined EXCHANGE
#  include "eos_vars.h"
# endif

# ifdef WEC
#  include "ocean2d.h"
#  ifdef VF_ROBUST
      real gam1,gam2
#  endif
      real cff1
# endif

# include "compute_auxiliary_bounds.h"

       do k=1,N
# include "compute_horiz_rhs_uv_terms.h"
       enddo

! Compute dynamic bottom drag coefficient.  Note that "rd" must be
! aliased to a horizontal array beyond #4 to avoid overlap with one
! of the vertical arrays used below. Because "rd" computed inside the
! include file has extra rows of points on the western and southern
! sides, it can interpolated to U- and V-points without being
! synchronized(exchanged). It is also copied into shared array "r_D".
! If recomputed here, it is the private array which must be used in
! the code below, so "r_D" needs to be redefined accordingly.

# ifdef RECOMPUTE_RD
#  define dtau dt
#  define rd wrk1
#  include "compute_rd_bott_drag.h"
#  undef rd
#  define r_D wrk1
# endif

      do j=jstr,jend

# include "compute_vert_rhs_uv_terms.h"
      
! Daniel add
# ifdef WEC
      enddo !---> end j-loop

      call compute_wec_rhs_uv_terms(istr,iend,jstr,jend,ru,rv,wrk1,wrk2)

      do j=jstr,jend ! ---> restart j-loop
# endif

! Daniel add
! Apply wave bottom streaming stress as a body force in the bottom-most
! grid cells only (if !defined BODY_FRICTION)
# ifdef WEC
#  if defined WAVE_FRICTION && !defined BODY_FRICTION
        do i=istrU,iend
          ru(i,j,1) =ru(i,j,1) + dm_u(i,j)*dn_u(i,j)*frc2dx(i,j)
        enddo
        do i=istr,iend
          rv(i,j,1) =rv(i,j,1) + dm_v(i,j)*dn_v(i,j)*frc2de(i,j)
        enddo
#  endif
# endif 


# ifdef IMPLICIT_BOTTOM_DRAG
        do i=istrU,iend
          DC(i,0)=dt*0.25*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))

          FC(i,N-1)= 2.*dt*(Akv(i,j,N-1)+Akv(i-1,j,N-1))
     &                      /(  Hz(i,j,N  )+Hz(i-1,j,N  )
     &                         +Hz(i,j,N-1)+Hz(i-1,j,N-1))

          WC(i,N-1)= DC(i,0)*0.5*(Wi(i,j,N-1)+Wi(i-1,j,N-1))

          cff=1./( 0.5*(Hz(i,j,N)+Hz(i-1,j,N))
     &                   +FC(i,N-1)-min(WC(i,N-1),0.) )

          CF(i,N-1)=cff*( FC(i,N-1)+max(WC(i,N-1),0.) )

          DC(i,N)=cff*( u(i,j,N,nnew) +DC(i,0)*ru(i,j,N)
!Daniel add
#  ifdef WIND_AT_RHO_POINTS

     &                    +dt*( 0.5 * (sustr(i,j)+sustr(i-1,j))
#    if defined WEC && defined SURFACE_BREAK
     &                            + brk2dx(i,j) )
#    endif 
#  else
     &                                       +dt*(sustr(i,j)
#    if defined WEC && defined SURFACE_BREAK
     &                            + brk2dx(i,j)
#    endif
     &                                                   )
#  endif
     &                                                  )
        enddo
        do k=N-1,2,-1      !--> forward elimination
          do i=istrU,iend
            FC(i,k-1)= 2.*dt*(Akv(i,j,k-1)+Akv(i-1,j,k-1))
     &                        /(  Hz(i,j,k  )+Hz(i-1,j,k  )
     &                           +Hz(i,j,k-1)+Hz(i-1,j,k-1))

            WC(i,k-1)= DC(i,0)*0.5*(Wi(i,j,k-1)+Wi(i-1,j,k-1))

            cff=1./( 0.5*(Hz(i,j,k)+Hz(i-1,j,k))
     &                           +FC(i,k-1)-min(WC(i,k-1),0.)
     &                             +FC(i,k)+max(WC(i,k),0.)
     &                    -CF(i,k)*(FC(i,k)-min(WC(i,k),0.))
     &                                                      )
            CF(i,k-1)=cff*(   FC(i,k-1)+max(WC(i,k-1),0.) )

            DC(i,k)=cff*( u(i,j,k,nnew) +DC(i,0)*ru(i,j,k)
     &                 +DC(i,k+1)*(FC(i,k)-min(WC(i,k),0.)) )
          enddo
        enddo
        do i=istrU,iend
          DC(i,1)=( u(i,j,1,nnew) +DC(i,0)*ru(i,j,1)
     &                      +DC(i,2)*(FC(i,1)-min(WC(i,1),0.))
     &                          )/( 0.5*(Hz(i,j,1)+Hz(i-1,j,1))
#  ifdef IMPLCT_NO_SLIP_BTTM_BC
     &                            +0.5*dt*(r_D(i,j)+r_D(i-1,j))
#  endif
     &                                +FC(i,1)+max(WC(i,1),0.)
     &                       -CF(i,1)*(FC(i,1)-min(WC(i,1),0.))
     &                                                       )

          u(i,j,1,nnew)=DC(i,1) * 0.5*(Hz(i,j,1)+Hz(i-1,j,1))

#  ifdef CORR_COUPLED_MODE
          rufrc(i,j)=ru(i,j,1) +dm_u(i,j)*dn_u(i,j)*(
#   ifdef WIND_AT_RHO_POINTS
     &                   0.5*(sustr(i-1,j)+sustr(i,j))
#   else
     &                                      sustr(i,j)
#   endif
! Daniel add
#   ifdef WEC
     &                                     + brk2dx(i,j)
#   endif     
     &            -0.5*(r_D(i-1,j)+r_D(i,j))*DC(i,1) )
#  endif
        enddo
        do k=2,N,+1          !--> backsubstitution
          do i=istrU,iend
            DC(i,k)=DC(i,k) +CF(i,k-1)*DC(i,k-1)

            u(i,j,k,nnew)=DC(i,k) * 0.5*(Hz(i,j,k)+Hz(i-1,j,k))

#  ifdef CORR_COUPLED_MODE
            rufrc(i,j)=rufrc(i,j) +ru(i,j,k)
#  endif
          enddo
        enddo

        if (j >= jstrV) then
          do i=istr,iend
            DC(i,0)=dt*0.25*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))

            FC(i,N-1)= 2.*dt*(Akv(i,j,N-1)+Akv(i,j-1,N-1))
     &                        /(  Hz(i,j,N  )+Hz(i,j-1,N  )
     &                           +Hz(i,j,N-1)+Hz(i,j-1,N-1))

            WC(i,N-1)= DC(i,0)*0.5*(Wi(i,j,N-1)+Wi(i,j-1,N-1))

            cff=1./( 0.5*(Hz(i,j,N)+Hz(i,j-1,N))
     &                     +FC(i,N-1)-min(WC(i,N-1),0.) )

            CF(i,N-1)=cff*( FC(i,N-1)+max(WC(i,N-1),0.) )

            DC(i,N)=cff*( v(i,j,N,nnew) +DC(i,0)*rv(i,j,N)
!Daniel add
#  ifdef WIND_AT_RHO_POINTS

     &                    +dt*( 0.5 * (svstr(i,j-1)+svstr(i,j))
#    if defined WEC && defined SURFACE_BREAK
     &                            + brk2de(i,j) )
#    endif 
#  else
     &                                       +dt*(svstr(i,j)
#    if defined WEC && defined SURFACE_BREAK
     &                            + brk2de(i,j)
#    endif
     &                                                   )
#  endif
     &                                                    )
          enddo
          do k=N-1,2,-1      !--> forward elimination
            do i=istr,iend
              FC(i,k-1)= 2.*dt*(Akv(i,j,k-1)+Akv(i,j-1,k-1))
     &                          /(  Hz(i,j,k  )+Hz(i,j-1,k  )
     &                             +Hz(i,j,k-1)+Hz(i,j-1,k-1))

              WC(i,k-1)= DC(i,0)*0.5*(Wi(i,j,k-1)+Wi(i,j-1,k-1))

              cff=1./( 0.5*(Hz(i,j,k)+Hz(i,j-1,k))
     &                              +FC(i,k-1)-min(WC(i,k-1),0.)
     &                                +FC(i,k)+max(WC(i,k),0.)
     &                       -CF(i,k)*(FC(i,k)-min(WC(i,k),0.))
     &                                                        )
              CF(i,k-1)=cff*( FC(i,k-1)+max(WC(i,k-1),0.) )

              DC(i,k)=cff*( v(i,j,k,nnew) +DC(i,0)*rv(i,j,k)
     &                   +DC(i,k+1)*(FC(i,k)-min(WC(i,k),0.)) )
            enddo
          enddo
          do i=istr,iend
            DC(i,1)=( v(i,j,1,nnew) +DC(i,0)*rv(i,j,1)
     &                       +DC(i,2)*(FC(i,1)-min(WC(i,1),0.))
     &                            )/( 0.5*(Hz(i,j,1)+Hz(i,j-1,1))
#  ifdef IMPLCT_NO_SLIP_BTTM_BC
     &                              +0.5*dt*(r_D(i,j)+r_D(i,j-1))
#  endif
     &                                  +FC(i,1)+max(WC(i,1),0.)
     &                         -CF(i,1)*(FC(i,1)-min(WC(i,1),0.))
     &                                                          )
            v(i,j,1,nnew)=DC(i,1) * 0.5*(Hz(i,j,1)+Hz(i,j-1,1))
#  ifdef CORR_COUPLED_MODE
            rvfrc(i,j)=rv(i,j,1) +dm_v(i,j)*dn_v(i,j)*(
#   ifdef WIND_AT_RHO_POINTS
     &                     0.5*(svstr(i,j-1)+svstr(i,j))
#   else
     &                                        svstr(i,j)
#   endif
! Daniel add
# ifdef WEC
     &                                      + brk2de(i,j)
# endif
     &              -0.5*(r_D(i,j-1)+r_D(i,j))*DC(i,1) )
#  endif
          enddo
          do k=2,N,+1          !--> backsubstitution
            do i=istr,iend
              DC(i,k)=DC(i,k) +CF(i,k-1)*DC(i,k-1)

              v(i,j,k,nnew)=DC(i,k) * 0.5*(Hz(i,j,k)+Hz(i,j-1,k))
#  ifdef CORR_COUPLED_MODE
              rvfrc(i,j)=rvfrc(i,j) + rv(i,j,k)
#  endif
            enddo
          enddo
        endif   !<-- j >= jstrV

# else  /* --> !IMPLICIT_BOTTOM_DRAG */

        do i=istrU,iend
          DC(i,0)=dt*0.25*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))

          u(i,j,N,nnew)=u(i,j,N,nnew) +DC(i,0)*ru(i,j,N)
#  ifdef WIND_AT_RHO_POINTS
     &                  +0.5*dt*(sustr(i-1,j)+sustr(i,j))
#  else
     &                                     +dt*sustr(i,j)
#  endif
!Daniel add
#  ifdef WEC
     &                                    + dt * brk2dx(i,j)
#  endif
          u(i,j,1,nnew)=u(i,j,1,nnew)+DC(i,0)*ru(i,j,1)
#  ifndef IMPLCT_NO_SLIP_BTTM_BC
     &                   -0.5*(r_D(i-1,j)+r_D(i,j))*u(i,j,1,nstp)
#  endif
#  ifdef CORR_COUPLED_MODE
          rufrc(i,j)=ru(i,j,1) +ru(i,j,N) +dn_u(i,j)*dm_u(i,j)*(
#   ifdef WIND_AT_RHO_POINTS
     &                              0.5*(sustr(i-1,j)+sustr(i,j))
#   else
     &                                                 sustr(i,j)
#   endif
! Daniel add
#   ifdef WEC
     &                                      + brk2dx(i,j)
#   endif
     &                 -0.5*(r_D(i-1,j)+r_D(i,j))*u(i,j,1,nstp) )
#  endif
        enddo
        do k=2,N-1
          do i=istrU,iend
            u(i,j,k,nnew)=u(i,j,k,nnew)+DC(i,0)*ru(i,j,k)
#  ifdef CORR_COUPLED_MODE
            rufrc(i,j)=rufrc(i,j) +ru(i,j,k)
#  endif
          enddo
        enddo
        if (j >= jstrV) then
          do i=istr,iend
            DC(i,0)=dt*0.25*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))

            v(i,j,N,nnew)=v(i,j,N,nnew) +DC(i,0)*rv(i,j,N)
#  ifdef WIND_AT_RHO_POINTS
     &                    +0.5*dt*(svstr(i,j-1)+svstr(i,j))
#  else
     &                                       +dt*svstr(i,j)
#  endif
! Daniel add
#  ifdef WEC
     &                                      +dt * brk2de(i,j)
#  endif
            v(i,j,1,nnew)=v(i,j,1,nnew) +DC(i,0)*rv(i,j,1)
#  ifndef IMPLCT_NO_SLIP_BTTM_BC
     &                    -0.5*(r_D(i,j-1)+r_D(i,j))*v(i,j,1,nstp)
#  endif
#  ifdef CORR_COUPLED_MODE
            rvfrc(i,j)=rv(i,j,1) +rv(i,j,N) +dm_v(i,j)*dn_v(i,j)*(
#   ifdef WIND_AT_RHO_POINTS
     &                                0.5*(svstr(i,j-1)+svstr(i,j))
#   else
     &                                                   svstr(i,j)
#   endif
!Daniel add
# ifdef WEC
     &                                               +brk2de(i,j)
# endif
     &                   -0.5*(r_D(i,j-1)+r_D(i,j))*v(i,j,1,nstp) )
#  endif
          enddo
          do k=2,N-1
            do i=istr,iend
              v(i,j,k,nnew)=v(i,j,k,nnew)+DC(i,0)*rv(i,j,k)
#  ifdef CORR_COUPLED_MODE
              rvfrc(i,j)=rvfrc(i,j) +rv(i,j,k)
#  endif
            enddo
          enddo
        endif

# endif /* IMPLICIT_BOTTOM_DRAG */
      enddo  !<-- j
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
