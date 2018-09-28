#include "cppdefs.h"
#ifdef SOLVE3D

      subroutine step3d_uv2(tile)
      implicit none
      integer tile
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call step3d_uv2_tile( istr,iend,jstr,jend, A2d(1,1), A2d(1,2),
     &                                           A2d(1,3), A2d(1,4))
      end

      subroutine step3d_uv2_tile(istr,iend,jstr,jend, WC,FC,CF,DC)
      implicit none
# include "param.h"
      integer istr,iend,jstr,jend, i,j,k
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N) :: WC,FC,CF,DC
      real cff

# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "coupling.h"
# include "scalars.h"
# if defined ADV_ISONEUTRAL || !defined IMPLICIT_BOTTOM_DRAG
#  include "mixing.h"
# endif
# ifdef ADV_ISONEUTRAL
#  include "eos_vars.h"
      real dRz,dRx_max, dpth, r0g, cfs,cfb
      real, parameter :: epsil=1.E-33
#  define  ALPHA_MAX 4.
#  define GAMMA 0.0833333333333
# endif
# ifdef PSOURCE
#  include "sources.h"
      integer is
# endif

# include "compute_auxiliary_bounds.h"

      do j=jstr,jend

# ifndef IMPLICIT_BOTTOM_DRAG

! Apply implicit time step for vertical viscosty term to horizontal
! velocity components,
!
!   du(k)     1     [         u(k+1)-u(k)             u(k)-u(k-1) ]
!  ------ = ----- * [ Akv(k)* ----------- - Akv(k-1)* ----------- ]
!    dt     Hz(k)   [           dz(k)                   dz(k-1)   ]
!
! where u(k) represents u,v(:,:,:,nnew) velocity components; Hz(k) and
! dz(k) are grid spacing between vertical W and RHO-points respectively
! interpolated to horizontal U,V locations. Backward Euler implicit
! time step requires solution of a tri-diagonal system,
!
!     -FC(k-1)*u_new(k-1) +[Hz(k)+FC(k-1)+FC(k)]*u_new(k)
!                       -FC(k)*u_new(k+1) = Hz(k)*u_old(k)
!
!                dt*Akv(k)
! where FC(k) = ----------- is normalized vertical viscosity defined
!                  dz(k)
!
! at W-points; u_new(k) is the new-time-step (unknown) velocity;
! u_old(k) is old-time-step tracer (known).  As long as vertical
! viscosity Ak(k) is nonnegative, the tri-diagonal matrix is
! diagonally dominant which guarantees stability of Gaussian
! elimination procedure, (e.g., Richtmeyer annd  Morton,  1967).
! Since top (wind) and bottom (drag) boundary stresses are applied
! explicitly, for the purpose of solving the tri-diagonal system
! boundary conditions are effectively no-stress Akv(N)=Akt(0)=0, so
! that the corresponding FC(N)=FC(1)=0. This leads to equations for
! top and bottom grid boxes,
!
!   -FC(N-1)*u_new(N-1) +[Hz(N)+FC(N-1)]*u_new(N) = Hz(N)*u_old(N)
!
!          [Hz(1)+FC(1)]*u_new(1) -FC(1)*u_new(2) = Hz(1)*u_old(1)
!
! Basically FC(N)=FC(1)=0 conditions should be interpreted as the
! implicit step merely redistributes the horizontal momentum
! components u,v throughout the vertical column, while conserving
! its vertical integral. At this moment the content of arrays
! u,v(:,:,:,nnew) has meaning of Hz*velocity. After the implicit
! step they becomes just velocity, which is temporarily loaded into
! array DC(:,1:N)
!>
! Also perform coupling between the barotropic and baroclinic modes.
! For the purpose of the computational efficiency (cache utilization)
! the implicit time step and coupling are partially overlapped. To
! distinguish the two operations in the code below we place separator
! (!>): lines above the separator (!>) belong to the tri-diagonal
! solver, while operations below it belong to the coupling.
! The coupling procedure makes vertical integration of the newly
! computed velocities [which ends up being stored in the scratch
! array DC(:,1:N) after the tri-diagonal problem is resolved] to
! compute its vertical mean [stored as DC(:,0)], and replace it with
! that from the 2D barotropic sub-model ubar=DU_avg1/(D*dn_u), where
! D=CF(:,0) is the total depth of water column.
!
!
! Time step momentum equation in the XI-direction
!----- ---- -------- -------- -- --- ------------

        do i=istrU,iend
          DC(i,0)=dt*0.25*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))

          FC(i,N-1)= 2.*dt*(Akv(i,j,N-1)+Akv(i-1,j,N-1))
     &                      /(  Hz(i,j,N  )+Hz(i-1,j,N  )
     &                         +Hz(i,j,N-1)+Hz(i-1,j,N-1))

          WC(i,N-1)= DC(i,0)*0.5*(Wi(i,j,N-1)+Wi(i-1,j,N-1))

          cff=1./( 0.5*(Hz(i,j,N)+Hz(i-1,j,N))
     &                   +FC(i,N-1)-min(WC(i,N-1),0.) )

          CF(i,N-1)=cff*( FC(i,N-1)+max(WC(i,N-1),0.) )

          DC(i,N)=cff*u(i,j,N,nnew)
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

            DC(i,k)=cff*( u(i,j,k,nnew)
     &                 +DC(i,k+1)*(FC(i,k)-min(WC(i,k),0.)) )
          enddo
        enddo       !--> discard WC,DC(i,0)
        do i=istrU,iend
          DC(i,1)=( u(i,j,1,nnew) +DC(i,2)*(FC(i,1)-min(WC(i,1),0.))
     &                                )/( 0.5*(Hz(i,j,1)+Hz(i-1,j,1))
#  ifdef IMPLCT_NO_SLIP_BTTM_BC
     &                                  +0.5*dt*(r_D(i,j)+r_D(i-1,j))
#  endif
     &                                      +FC(i,1)+max(WC(i,1),0.)
     &                             -CF(i,1)*(FC(i,1)-min(WC(i,1),0.))
     &                                                       )
!>
          CF(i,0)=0.5*(Hz(i,j,1)+Hz(i-1,j,1))    !--> coupling
          DC(i,0)=CF(i,0)*DC(i,1)
        enddo
        do k=2,N,+1          !--> backsubstitution
          do i=istrU,iend
            DC(i,k)=DC(i,k) +CF(i,k-1)*DC(i,k-1) !<-- backsubstitution
!>
            cff=0.5*(Hz(i,j,k)+Hz(i-1,j,k))      !--> coupling
            CF(i,0)=CF(i,0)+cff                  !    [ vertical    ]
            DC(i,0)=DC(i,0)+cff*DC(i,k)          !    [ integration ]
          enddo
        enddo
        do i=istrU,iend
          DC(i,0)=( DC(i,0)*dn_u(i,j)-DU_avg1(i,j)  !<-- mismatch
     &                       )/(CF(i,0)*dn_u(i,j))
        enddo
        do k=N,1,-1
          do i=istrU,iend
            u(i,j,k,nnew)=(DC(i,k)-DC(i,0))      !<-- correction
#  ifdef MASKING
     &                          *umask(i,j)
#  endif
          enddo
        enddo        !--> discard DC,CF

# else /* IMPLICIT_BOTTOM_DRAG */

        do i=istrU,iend                          !--> coupling only
          CF(i,0)=0.5*(Hz(i,j,N)+Hz(i-1,j,N))
          DC(i,0)=u(i,j,N,nnew)
          u(i,j,N,nnew)=u(i,j,N,nnew)/CF(i,0)
        enddo
        do k=N-1,1,-1
          do i=istrU,iend
            cff=0.5*(Hz(i,j,k)+Hz(i-1,j,k))
            CF(i,0)=CF(i,0) +cff
            DC(i,0)=DC(i,0) +u(i,j,k,nnew)
            u(i,j,k,nnew)=u(i,j,k,nnew)/cff
          enddo
        enddo
        do i=istrU,iend
          DC(i,0)=( DC(i,0)*dn_u(i,j)-DU_avg1(i,j)  !<-- mismatch
     &                       )/(CF(i,0)*dn_u(i,j))
        enddo
        do k=1,N
          do i=istrU,iend
            u(i,j,k,nnew)=(u(i,j,k,nnew)-DC(i,0))
#  ifdef MASKING
     &                                *umask(i,j)
#  endif
          enddo
        enddo

# endif

        if (j>=jstrV) then

# ifndef IMPLICIT_BOTTOM_DRAG

! Time step momentum equation in the ETA-direction.
!----- ---- -------- -------- -- --- --------------
! (All steps are similar to XI-component above.)

          do i=istr,iend
            DC(i,0)=dt*0.25*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))

            FC(i,N-1)= 2.*dt*(Akv(i,j,N-1)+Akv(i,j-1,N-1))
     &                        /(  Hz(i,j,N  )+Hz(i,j-1,N  )
     &                           +Hz(i,j,N-1)+Hz(i,j-1,N-1))

            WC(i,N-1)= DC(i,0)*0.5*(Wi(i,j,N-1)+Wi(i,j-1,N-1))

            cff=1./( 0.5*(Hz(i,j,N)+Hz(i,j-1,N))
     &                     +FC(i,N-1)-min(WC(i,N-1),0.) )

            CF(i,N-1)=cff*( FC(i,N-1)+max(WC(i,N-1),0.) )

            DC(i,N)=cff*v(i,j,N,nnew)
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

              DC(i,k)=cff*( v(i,j,k,nnew)
     &                   +DC(i,k+1)*(FC(i,k)-min(WC(i,k),0.)) )
            enddo
          enddo
          do i=istr,iend
            DC(i,1)=( v(i,j,1,nnew) +DC(i,2)*(FC(i,1)-min(WC(i,1),0.))
     &                                  )/( 0.5*(Hz(i,j,1)+Hz(i,j-1,1))
#  ifdef IMPLCT_NO_SLIP_BTTM_BC
     &                                    +0.5*dt*(r_D(i,j)+r_D(i,j-1))
#  endif
     &                                        +FC(i,1)+max(WC(i,1),0.)
     &                               -CF(i,1)*(FC(i,1)-min(WC(i,1),0.))
     &                                                              )
!>
            CF(i,0)=0.5*(Hz(i,j,1)+Hz(i,j-1,1))   !--> coupling
            DC(i,0)=CF(i,0)*DC(i,1)
          enddo
          do k=2,N,+1          !--> backsubstitution
            do i=istr,iend
              DC(i,k)=DC(i,k)+CF(i,k-1)*DC(i,k-1) !<-- backsubstitution
!>
              cff=0.5*(Hz(i,j,k)+Hz(i,j-1,k))     !--> coupling
              CF(i,0)=CF(i,0)+cff                 !    [ vertical    ]
              DC(i,0)=DC(i,0)+cff*DC(i,k)         !    [ integration ]
            enddo
          enddo
          do i=istr,iend
            DC(i,0)=( DC(i,0)*dm_v(i,j)-DV_avg1(i,j)  !<-- mismatch
     &                         )/(CF(i,0)*dm_v(i,j))
          enddo
          do k=N,1,-1
            do i=istr,iend
              v(i,j,k,nnew)=(DC(i,k)-DC(i,0))     !<-- correction
#  ifdef MASKING
     &                            *vmask(i,j)
#  endif
            enddo
          enddo       !--> discard DC,CF(:,0)

# else  /* IMPLICIT_BOTTOM_DRAG */

          do i=istr,iend
            CF(i,0)=0.5*(Hz(i,j,N)+Hz(i,j-1,N))
            DC(i,0)=v(i,j,N,nnew)
            v(i,j,N,nnew)=v(i,j,N,nnew)/CF(i,0)
          enddo
          do k=N-1,1,-1
            do i=istr,iend
              cff=0.5*(Hz(i,j,k)+Hz(i,j-1,k))
              CF(i,0)=CF(i,0) +cff
              DC(i,0)=DC(i,0) +v(i,j,k,nnew)
              v(i,j,k,nnew)=v(i,j,k,nnew)/cff
            enddo
          enddo      !--> discard FC, keep DC,CF(:,0)
          do i=istr,iend
            DC(i,0)=( DC(i,0)*dm_v(i,j)-DV_avg1(i,j)
     &                         )/(CF(i,0)*dm_v(i,j))
          enddo
          do k=1,N
            do i=istr,iend
              v(i,j,k,nnew)=(v(i,j,k,nnew)-DC(i,0))
#  ifdef MASKING
     &                                  *vmask(i,j)
#  endif
            enddo
          enddo       !--> discard DC,CF(:,0)
# endif
        endif  !<-- j>=jstrV
      enddo      ! <-- j

! Set PHYSICAL lateral boundary conditions.

      call u3dbc_tile (istr,iend,jstr,jend, DC)
      call v3dbc_tile (istr,iend,jstr,jend, DC)

# ifdef PSOURCE
      do is=1,Nsrc           ! Apply point mass sources associated
        i=Isrc(is)           ! with river runoff simulations.
        j=Jsrc(is)
        if (istrR<=i .and. i<=iendR .and.
     &      jstrR<=j .and. j<=jendR) then
          if (Dsrc(is)==0) then
            do k=1,N
              u(i,j,k,nnew)=2.*Qsrc(is,k)/( dn_u(i,j)*(
     &                       z_w(i-1,j,k)-z_w(i-1,j,k-1)
     &                      +z_w(i  ,j,k)-z_w(i  ,j,k-1)
     &                                                ))
            enddo
          else
            do k=1,N
              v(i,j,k,nnew)=2.*Qsrc(is,k)/( dm_v(i,j)*(
     &                       z_w(i,j-1,k)-z_w(i,j-1,k-1)
     &                      +z_w(i,j  ,k)-z_w(i,j  ,k-1)
     &                                                ))
            enddo
          endif
        endif
      enddo
# endif

! Coupling 2D and 3D momentum equations:
!--------- -- --- -- -------- ----------
! Compute inaccurate vertical integral of 3D velocity field, and
! replace it with the one computed by barotropic mode. After that
! compute mass fluxes through grid box faces.
!
! Meaning of scratch variables below:
!
! DC(i,k) [where k=1,N]  height of grid box for U- or V-cell.
! DC(i,0) total depth of water column at horizontal U- or V-points.
! FC(i,0) vertically integrated mass flux, error, and correction for
!                       the new time step velocity u,v(:,:,:,nnew)
! CF(i,0) vertically integrated mass flux/error/correction for the
!             intermediate time step [n+1/2] mass fluxes FlxU,FlxV.
!
! This procedure also replaces instantaneous (with respect to the
! fast time step) vertically integrated (barotropic) velocities with
! their values based on fast-time-averaged mass fluxes. These are to
! be used as initial conditions for the barotropic mode at the new
! time step.
!
! Explanation of horizontal loop indices: in the case of periodic
! boundaries (in either direction) the coupling and computation of
! mass fluxes FlxU,FlxV is performed within the internal range of
! indices (excluding ghost zones), after that the ghost points for
! the newly computed arrays are exchanged; in the case of nonperiodic
! boundaries the coupling is done over the extended range of indices
! (that is including boundary points).

# ifdef EW_PERIODIC
#  define IU_RANGE istr,iend
#  define IV_RANGE istr,iend
# else
#  define IU_RANGE istr,iendR
#  define IV_RANGE istrR,iendR
# endif

# ifdef NS_PERIODIC
#  define J_RANGE jstr,jend
# else
#  define J_RANGE jstrR,jendR
# endif

      do j=J_RANGE                  !--> XI-component
        do i=IU_RANGE
          DC(i,N)=0.5*(Hz(i,j,N)+Hz(i-1,j,N))*dn_u(i,j)
          DC(i,0)=DC(i,N)
          FC(i,0)=DC(i,N)*u(i,j,N,nnew)
        enddo
        do k=N-1,1,-1
          do i=IU_RANGE
            DC(i,k)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)
            DC(i,0)=DC(i,0)+DC(i,k)
            FC(i,0)=FC(i,0)+DC(i,k)*u(i,j,k,nnew)
          enddo
        enddo
        do i=IU_RANGE
          DC(i,0)=1./DC(i,0)
          ubar(i,j,knew)=DC(i,0)*DU_avg1(i,j)
          FC(i,0)=DC(i,0)*(FC(i,0)-DU_avg1(i,j))
!>
          CF(i,0)=0.   !<-- initialize vertical integral
        enddo          !--> keep DC(:,0) for later
        do k=N,1,-1
          do i=IU_RANGE
            u(i,j,k,nnew)=(u(i,j,k,nnew)-FC(i,0))
# ifdef MASKING
     &                                  *umask(i,j)
# endif
c--# define TR
# ifdef TR
            CF(i,k)=0.5*DC(i,k)*(u(i,j,k,nstp)+u(i,j,k,nnew))
# else
c--#define EPSIL 0.125
c--#define DELTA 0.75

#  define EPSIL 0.36
#  define DELTA 0.28

c>># define EPSIL 0.2105
c>># define DELTA 0.5790

            CF(i,k)=DELTA*FlxU(i,j,k) + EPSIL*DC(i,k)
     &                         *(u(i,j,k,nstp)+u(i,j,k,nnew))
# endif
            CF(i,0)=CF(i,0)+CF(i,k)    !<-- vertical integration
          enddo
        enddo
        do i=IU_RANGE
          CF(i,0)=DC(i,0)*(CF(i,0)-DU_avg2(i,j))  !<-- mismatch
        enddo
        do k=1,N,+1
          do i=IU_RANGE
            FlxU(i,j,k)=CF(i,k)-DC(i,k)*CF(i,0)   !<-- correction

# ifdef ADV_ISONEUTRAL
            diff3u(i,j,k)=sqrt(GAMMA*abs(u(i,j,k,nnew))*dm_u(i,j))
     &                                                  *dm_u(i,j)
# endif
          enddo
        enddo

        if (j >= jstr) then         !--> ETA-component
          do i=IV_RANGE
            DC(i,N)=0.5*(Hz(i,j,N)+Hz(i,j-1,N))*dm_v(i,j)
            DC(i,0)=DC(i,N)
            FC(i,0)=DC(i,N)*v(i,j,N,nnew)
          enddo
          do k=N-1,1,-1
            do i=IV_RANGE
              DC(i,k)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)
              DC(i,0)=DC(i,0)+DC(i,k)
              FC(i,0)=FC(i,0)+DC(i,k)*v(i,j,k,nnew)
            enddo
          enddo
          do i=IV_RANGE
            DC(i,0)=1./DC(i,0)
            vbar(i,j,knew)=DC(i,0)*DV_avg1(i,j)
            FC(i,0)=DC(i,0)*(FC(i,0)-DV_avg1(i,j))
!>
            CF(i,0)=0.
          enddo        !--> keep DC(:,0) for later
          do k=N,1,-1
            do i=IV_RANGE
              v(i,j,k,nnew)=(v(i,j,k,nnew)-FC(i,0))
# ifdef MASKING
     &                                    *vmask(i,j)
# endif
# ifdef TR
              CF(i,k)=0.5*DC(i,k)*(v(i,j,k,nstp)+v(i,j,k,nnew))
# else
              CF(i,k)=DELTA*FlxV(i,j,k) + EPSIL*DC(i,k)
     &                           *(v(i,j,k,nstp)+v(i,j,k,nnew))
# endif
              CF(i,0)=CF(i,0)+CF(i,k)
            enddo
          enddo
          do i=IV_RANGE
            CF(i,0)=DC(i,0)*(CF(i,0)-DV_avg2(i,j))
          enddo
          do k=1,N,+1
            do i=IV_RANGE
              FlxV(i,j,k)=CF(i,k)-DC(i,k)*CF(i,0)

# ifdef ADV_ISONEUTRAL
              diff3v(i,j,k)=sqrt(GAMMA*abs(v(i,j,k,nnew))*dn_v(i,j))
     &                                                   *dn_v(i,j)
# endif
            enddo
          enddo
# ifdef ADV_ISONEUTRAL
          if (j <= jend) then
            r0g=rho0/g
            do k=1,N-1
              do i=istr,iend
#  define ISO
#  ifdef ISO
#    ifdef SPLIT_EOS
                dpth=-0.5*(z_r(i,j,k+1)+z_r(i,j,k))
                dRz=rho1(i,j,k)-rho1(i,j,k+1)+(qp1(i,j,k)-qp1(i,j,k+1))
     &                                      *dpth*(1.-2.*qp2*dpth)
#    else
                dRz=rho(i,j,k)-rho(i,j,k+1)
#    endif
                dRz=max(dRz,0.) +r0g*f(i,j)**2
     &                        *(z_r(i,j,k+1)-z_r(i,j,k))
#  else
                dRz=z_r(i,j,k+1)-z_r(i,j,k)
#  endif


c*              dRx_max=(z_r(i,j,k+1)-z_r(i,j,k))*max(
c*   &                        abs(dRdx(i,j,k  )), abs(dRdx(i+1,j,k  )),
c*   &                        abs(dRdx(i,j,k+1)), abs(dRdx(i+1,j,k+1)),
c*
c*   &                        abs(dRde(i,j,k  )), abs(dRde(i,j+1,k  )),
c*   &                        abs(dRde(i,j,k+1)), abs(dRde(i,j+1,k+1))
c*   &                                                               )

               dRx_max=max(

     &       dm_u(i  ,j)*max(abs(dRdx(i  ,j,k)), abs(dRdx(i  ,j,k+1))),
     &       dm_u(i+1,j)*max(abs(dRdx(i+1,j,k)), abs(dRdx(i+1,j,k+1))),

     &       dn_v(i,j  )*max(abs(dRde(i,j  ,k)), abs(dRde(i,j  ,k+1))),
     &       dn_v(i,j+1)*max(abs(dRde(i,j+1,k)), abs(dRde(i,j+1,k+1)))
     &                                                               )


/*
                cff=(z_w(i,j,N)-z_w(i,j,k))/max( hbls(i,j,3-nstp),
c**  &                                     z_w(i,j,N)-z_w(i,j,N-5) )
c??     &                                     z_w(i,j,N)-z_w(i,j,N-1) )
c??     &                                      20. ) !<-- 20 meters
     &                                      50. ) !<-- 50 meters
                if (cff < 1.) then
                  cff=ALPHA_MAX*cff*cff*(3.-2.*cff)
                else
                  cff=ALPHA_MAX
                endif

                cff=ALPHA_MAX*cff*cff*(3.-2.*cff)
*/


                cfs=min(1., (z_w(i,j,N)-z_w(i,j,k))/max( 50.,
     &                                    hbls(i,j,3-nstp) ))
                cfb=min(1., (z_w(i,j,k)-z_w(i,j,0))/max( 50.,
     &                                    hbbl(i,j,3-nstp) )) 

                cff=ALPHA_MAX*cfs*(2.-cfs) * cfb*(2.-cfb)
 

c               idRz(i,j,k)=1./dRz

               idRz(i,j,k)=cff/max(cff*dRz, dRx_max, epsil)

c             idRz(i,j,k)=ALPHA_MAX/max(ALPHA_MAX*dRz,dRx_max,epsil)

              enddo
            enddo
          endif      !<-- j <= jend
# endif
        endif    !<-- j >= jstr
      enddo    ! <-- j

# undef IU_RANGE
# undef IV_RANGE
# undef J_RANGE

! Exchange periodic boundaries and computational margins.

# ifdef EXCHANGE
      call exchange_3_tile (istr,iend,jstr,jend,   FlxU,N,
     &                         u(START_2D_ARRAY,1,nnew), N,
     &                         ubar(START_2D_ARRAY,knew),1)

      call exchange_3_tile (istr,iend,jstr,jend,   FlxV,N,
     &                         v(START_2D_ARRAY,1,nnew), N,
     &                         vbar(START_2D_ARRAY,knew),1)

#  ifdef ADV_ISONEUTRAL
      call exchange_3_tile (istr,iend,jstr,jend, idRz,N+1,
     &                                 diff3u,N, diff3v,N)
#  endif
# endif
      end


      subroutine check_step_uv2_switches(ierr)

! This code must be placed here rather than in a separate file in order
! to be exposed to the relevant CPP-macros defined locally in this file
! above.  It does not affect any model results, other than signature in
! global attribute "CPPS" in output netCDF files.

      implicit none
      integer ierr, is,ie, lenstr
# include "param.h"
# include "strings.h"
      integer lstr
      character(len=32) tmpstr
      ie=lenstr(cpps)
      is=ie+2 ; ie=is+13
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='<step3d_uv2.F>'
# ifdef DELTA
      write(tmpstr,*) DELTA ;  lstr=lenstr(tmpstr)
      do while(tmpstr(lstr:lstr)=='0' .and. lstr>3)
        lstr=lstr-1
      enddo
      is=ie+2 ; ie=is+5+lstr
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='DELTA='/ /tmpstr(1:lstr)
# endif
# ifdef EPSIL
      write(tmpstr,*) EPSIL  ;  lstr=lenstr(tmpstr)
      do while(tmpstr(lstr:lstr)=='0' .and. lstr>3)
        lstr=lstr-1
      enddo
      is=ie+2 ; ie=is+5+lstr
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='EPSIL='/ /tmpstr(1:lstr)
# endif
# ifdef GAMMA
      write(tmpstr,*) GAMMA ; lstr=lenstr(tmpstr)
      do while(tmpstr(lstr:lstr)=='0' .and. lstr>3)
        lstr=lstr-1
      enddo
      is=ie+2 ; ie=is+5+lstr
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='GAMMA='/ /tmpstr(1:lstr)
# endif
# ifdef ALPHA_MAX
      write(tmpstr,'(F16.12)') ALPHA_MAX
      lstr=lenstr(tmpstr)
      do while(tmpstr(lstr:lstr) == '0' .and. lstr > 3)
        lstr=lstr-1
      enddo
      is=ie+2 ; ie=is+9+lstr
      if (ie > max_opt_size) goto 99
      cpps(is:ie)='ALPHA_MAX='/ /tmpstr(1:lstr)
# endif
      return
  99  mpi_master_only write(*,'(/1x,2A/12x,A/)')      '### ERROR: ',
     &  'Insufficient length of string "cpps" in file "strings.h".',
     &        'Increase parameter "max_opt_size" it and recompile.'
      ierr=ierr+1
      end
#else
      subroutine step3d_uv_empty
      end
#endif /* SOLVE3D */
