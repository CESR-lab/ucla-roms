      module advection

      ! Contains various advection algorithms for momentum and tracers
      !
      ! Pay attention!
      ! The uv_hadv/uv_vadv routines add a flux divergence
      ! to the incoming rhs terms (ru,rv). The t_hadv/t_vadv routines
      ! overwrite advective fluxes. (Fx,Fy,Fz)

      use dimensions
      use diagnostics
      use ocean_vars
      use grid
      use tracers

      implicit none

      private

#include "advection.opt"

      ! These control whether we use 'real' or 'modified' up3 and c4
      real,parameter :: coef_up3 = 1./6. ! coef_up3 = 1./8. for 'real' up3
      real,parameter :: coef_c4 = 1./12. ! coef_c4 = 1./16. for 'real' c4
      real,parameter :: delta = 0.1666666666666667 !! delta=0.125
      real,parameter :: gamma = 0.3333333333333333 !! gamma=0.25

      real,allocatable,dimension(:,:) :: wrk1
      real,allocatable,dimension(:,:) :: wrk2
      real,allocatable,dimension(:,:) :: UFx
      real,allocatable,dimension(:,:) :: UFe
      real,allocatable,dimension(:,:) :: VFx
      real,allocatable,dimension(:,:) :: VFe
      real,allocatable,dimension(:,:) :: FC
      real,allocatable,dimension(:,:) :: crv
      real,allocatable,dimension(:,:) :: CF
      real,allocatable,dimension(:,:) :: DC

      public init_advection
      public uv_cor_curv
      public uv_hadv_pre
      public uv_hadv_cor
      public uv_vadv_pre
      public uv_vadv_cor
      public t_hadv_pre
      public t_hadv_cor
      public t_vadv_pre
      public t_vadv_cor
      public uv_hadv_c4
      public uv_vadv_c4

      contains

!------------------------------------------------------------------------------
      subroutine init_advection ![
      ! allocate work space and check options

      if (.not.uv_adv.and.curvgrid) then
        print *, 'Curvgrid cant be true when uv_adv is false'
        stop
      endif

      allocate(wrk1(GLOBAL_2D_ARRAY))
      allocate(wrk2(GLOBAL_2D_ARRAY))
      allocate(UFx(GLOBAL_2D_ARRAY))
      allocate(UFe(GLOBAL_2D_ARRAY))
      allocate(VFx(GLOBAL_2D_ARRAY))
      allocate(VFe(GLOBAL_2D_ARRAY))
      allocate(crv(nx,nz))
      allocate(FC(nx,0:nz))
      allocate(CF(nx,0:nz))
      allocate(DC(nx,0:nz))

      end subroutine init_advection !]
!------------------------------------------------------------------------------
      subroutine uv_cor_curv(ru,rv,k) ![
      ! Add Coriolis terms and the contribution to advection that is
      ! associated with curvilinear horizontal coordinates.

!DIR$ ATTRIBUTES FORCEINLINE :: uv_cor_curv
      ! Import
      integer, intent(in) :: k
      real,dimension(GLOBAL_2D_ARRAY,nz) :: ru,rv
      ! Local
      integer :: i,j
      real    :: cff,cff1

      if ((.not.uv_cor).and.(.not.curvgrid)) return

      do j=0,ny
        do i=0,nx
          cff=0.5*Hz(i,j,k)*(
# ifdef UV_COR
     &              fomn(i,j)
# endif
# if (defined CURVGRID && defined UV_ADV)
     &             +0.5*( (v(i,j,k,nrhs)+v(i,j+1,k,nrhs))*dndx(i,j)
     &                   -(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))*dmde(i,j))
# endif
     &                                                             )
# ifdef WEC
#  if (defined CURVGRID && defined UV_ADV)
     ! Define a cff1 that is added to UFx, VFe for stokes terms
          cff1=0.5*Hz(i,j,k)*(
     &          0.5*( dndx(i,j)*(vst(i,j,k)+vst(i,j+1,k))
     &               -dmde(i,j)*(ust(i,j,k)+ust(i+1,j,k)) ))
#  else
          cff1 = 0.0
#  endif
          UFx(i,j)=(cff+cff1)*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))
          UFe(i,j)=cff*(vst(i,j,k)+vst(i,j+1,k))
          VFe(i,j)=(cff+cff1)*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))
          VFx(i,j)=cff*(ust(i,j,k)+ust(i+1,j,k))
# else /* no WEC */

            UFx(i,j)=cff*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))
            VFe(i,j)=cff*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))
# endif
        enddo
      enddo
      do j=1,ny
        do i=1,nx
          ru(i,j,k)=ru(i,j,k)+0.5*(UFx(i,j)+UFx(i-1,j))
# ifdef WEC
     &               + 0.5*(UFe(i,j)+UFe(i-1,j))
# endif
        enddo
      enddo
      do j=1,ny
        do i=1,nx
            rv(i,j,k)=rv(i,j,k)-0.5*(VFe(i,j)+VFe(i,j-1))
# ifdef WEC
     &               -0.5*(VFx(i,j)+VFx(i,j-1))
# endif
        enddo
      enddo

      if (CORR_STAGE) then
        if (diag_uv.and.calc_diag) then
          Udiag(1:nx,1:ny,k,icori) = 0.5*(UFx(1:nx,1:ny)+UFx(0:nx-1,1:ny))*dxdyi_u
          Vdiag(1:nx,1:ny,k,icori) =-0.5*(VFe(1:nx,1:ny)+VFe(1:nx,0:ny-1))*dxdyi_v
        endif
      endif

      end subroutine uv_cor_curv !]
!------------------------------------------------------------------------------
      subroutine uv_hadv_pre(ru,rv,k) ![
      ! Import
      integer, intent(in) :: k
      real,dimension(GLOBAL_2D_ARRAY,nz) :: ru,rv

      if (.not.uv_adv) return

      select case (hadv_scheme)
      case (adv_up3)
        ! Leave the diffusive part of up3 for the corrector step
        call uv_hadv_c4(ru,rv,k)
      case (adv_c4)
        call uv_hadv_c4(ru,rv,k)
      case default
        print *, 'advection scheme not implemented yet'
        stop
      end select

      end subroutine  uv_hadv_pre !]
!------------------------------------------------------------------------------
      subroutine uv_hadv_cor(ru,rv,k) ![
      ! Momentum advection during the corrector step

      ! Import/export
      integer, intent(in) :: k
      real,dimension(GLOBAL_2D_ARRAY,nz) :: ru,rv

      if (.not.uv_adv) return

      select case (hadv_scheme)
      case (adv_up3)
        call uv_hadv_up3(ru,rv,k)
      case (adv_c4)
        call uv_hadv_c4(ru,rv,k)
      case default
        print *, 'advection scheme not implemented yet'
        stop
      end select

      end subroutine  uv_hadv_cor !]
!------------------------------------------------------------------------------
      subroutine uv_vadv_pre(ru,rv,j) ![
      ! Vertical advective fluxes for horizontal momentum equations
      ! during the predictor step

      ! Import
      integer, intent(in) :: j
      real,dimension(GLOBAL_2D_ARRAY,nz) :: ru,rv

      if (.not.uv_adv) return

      select case (vadv_scheme)
      case (adv_up3)
        ! Leave the diffusive part of up3 for the corrector step
        call uv_vadv_c4(ru,rv,j)
      case (adv_c4)
        call uv_vadv_c4(ru,rv,j)
      case (adv_spline)
        call uv_vadv_spline(ru,rv,j)
      case default
        print *, 'advection scheme not implemented yet'
        stop
      end select

      end subroutine  uv_vadv_pre !]
!------------------------------------------------------------------------------
      subroutine uv_vadv_cor(ru,rv,j) ![
      ! Vertical advective fluxes for horizontal momentum equations
      ! during the corrector step

      ! Import
      integer, intent(in) :: j
      real,dimension(GLOBAL_2D_ARRAY,nz) :: ru,rv

      if (.not.uv_adv) return

      select case (vadv_scheme)
      case (adv_up3)
        call uv_vadv_up3(ru,rv,j)
      case (adv_c4)
        call uv_vadv_c4(ru,rv,j)
      case (adv_spline)
        call uv_vadv_spline(ru,rv,j)
      case default
        print *, 'advection scheme not implemented yet'
        stop
      end select

      end subroutine  uv_vadv_cor !]
!------------------------------------------------------------------------------
      subroutine t_hadv_pre(Fx,Fy,k,itrc) ![
!DIR$ ATTRIBUTES FORCEINLINE :: t_hadv_pre

      ! Horizontal tracer advection fluxes during predictor step

      ! Import
      integer, intent(in) :: k,itrc
      real,dimension(GLOBAL_2D_ARRAY) :: Fx,Fy
      ! local

      select case (hadv_scheme)
      case (adv_up3)
        ! Leave the diffusive part of up3 for the corrector step
        call t_hadv_c4(Fx,Fy,k,itrc)
      case (adv_c4)
        call t_hadv_c4(Fx,Fy,k,itrc)
      case default
        print *, 'advection scheme not implemented yet'
        stop
      end select

      end subroutine  t_hadv_pre !]
!----------------------------------------------------------------
      subroutine t_hadv_cor(Fx,Fy,k,itrc) ![
!DIR$ ATTRIBUTES FORCEINLINE :: t_hadv_cor

      ! Horizontal tracer advection fluxes during corrector step

      ! Import
      integer, intent(in) :: k,itrc
      real,dimension(GLOBAL_2D_ARRAY) :: Fx,Fy
      ! local

      select case (hadv_scheme)
      case (adv_up3)
        if (adv_isoneutral) then
          call t_hadv_c4(Fx,Fy,k,itrc)
        else
          call t_hadv_up3(Fx,Fy,k,itrc)
        endif
      case (adv_c4)
        call t_hadv_c4(Fx,Fy,k,itrc)
      case default
        print *, 'advection scheme not implemented yet'
        stop
      end select

      end subroutine  t_hadv_cor !]
!----------------------------------------------------------------
      subroutine t_vadv_pre(Fz,j,itrc) ![
!DIR$ ATTRIBUTES FORCEINLINE :: t_vadv_pre

      ! Vertical tracer advection fluxes during predictor step

      ! Import
      integer, intent(in) :: j,itrc
      real,dimension(1-bf:nx+bf,0:nz) :: Fz
      ! local

      select case (vadv_scheme)
      case (adv_up3)
        ! Leave the diffusive part of up3 for the corrector step
        call t_vadv_c4(Fz,j,itrc)
      case (adv_c4)
        call t_vadv_c4(Fz,j,itrc)
      case (adv_spline)
        call t_vadv_spline(Fz,j,itrc)
      case default
        print *, 'advection scheme not implemented yet'
        stop
      end select

      end subroutine  t_vadv_pre !]
!----------------------------------------------------------------
      subroutine t_vadv_cor(Fz,j,itrc) ![
!DIR$ ATTRIBUTES FORCEINLINE :: t_vadv_cor

      ! Vertical tracer advection fluxes during corrector step

      ! Import
      integer, intent(in) :: j,itrc
      real,dimension(1-bf:nx+bf,0:nz) :: Fz
      ! local

      select case (vadv_scheme)
      case (adv_up3)
        call t_vadv_up3(Fz,j,itrc)
      case (adv_c4)
        call t_vadv_c4(Fz,j,itrc)
      case (adv_spline)
        call t_vadv_spline(Fz,j,itrc)
      case default
        print *, 'advection scheme not implemented yet'
        stop
      end select

      end subroutine  t_vadv_cor !]
!----------------------------------------------------------------
      subroutine uv_hadv_up3(ru,rv,k) ![
      ! computes 1 level of up3 advection

      ! Upwind  3rd order:         :  ( 3*f(k+1) + 6*f(k) -f(k-1))/8 (we>0)
      ! 'Modified' upwind 3rd order:  ( 2*f(k+1) + 5*f(k) -f(k-1))/6 (we>0)
      ! we get the modified version if coef_up3 is set to 1/6 in advection.opt


      ! Import
      integer, intent(in) :: k
      real,dimension(GLOBAL_2D_ARRAY,nz) :: ru,rv
      ! local
      integer :: i,j
      real    :: cff

      do j=1,ny
        do i=0,nx+1
          wrk1(i,j)=u(i-1,j,k,nrhs)-2.*u(i,j,k,nrhs)+u(i+1,j,k,nrhs)
          wrk2(i,j)=FlxU(i-1,j,k) -2.*FlxU(i,j,k) +FlxU(i+1,j,k)
        enddo
      enddo

      if (west_bnd) then
        do j=1,ny
          wrk1(1,j) = wrk1(2,j)
          wrk2(1,j) = wrk2(2,j)
        enddo
      endif
      if (east_bnd) then
        do j=1,ny
          wrk1(nx+1,j) = wrk1(nx,j)
          wrk2(nx+1,j) = wrk2(nx,j)
        enddo
      endif
      do j=1,ny
        do i=0,nx
          cff=FlxU(i,j,k)+FlxU(i+1,j,k)-delta*( wrk2(i  ,j)
     &                                            +wrk2(i+1,j))
          UFx(i,j)=0.25*( cff*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))
     &                          -gamma*( max(cff,0.)*wrk1(i  ,j)
     &                                  +min(cff,0.)*wrk1(i+1,j)
     &                                                      ))
        enddo
      enddo

# define vee wrk1
# define Hvee wrk2
      do j=0,ny+1
        do i=1,nx
          vee(i,j)=v(i,j-1,k,nrhs)-2.*v(i,j,k,nrhs)+v(i,j+1,k,nrhs)
          Hvee(i,j)=FlxV(i,j-1,k) -2.*FlxV(i,j,k)  +FlxV(i,j+1,k)
        enddo
      enddo

      if (south_bnd) then
        do i=1,nx
          vee(i,1) =vee(i,2)
          Hvee(i,1)=Hvee(i,2)
        enddo
      endif
      if (north_bnd) then
        do i=1,nx
          vee(i,ny+1)=vee(i,ny)
          Hvee(i,ny+1)=Hvee(i,ny)
        enddo
      endif

      do j=0,ny
        do i=1,nx
          cff=FlxV(i,j,k)+FlxV(i,j+1,k)-delta*( Hvee(i,j  )
     &                                           +Hvee(i,j+1))
          VFe(i,j)=0.25*( cff*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))
     &                          -gamma*( max(cff,0.)*vee(i,j  )
     &                                  +min(cff,0.)*vee(i,j+1)
     &                                                      ))
        enddo
      enddo
# undef Hvee
# undef vee

# define uee wrk1
# define Hvxx wrk2

      do j=0,ny+1
        do i=1,nx
          uee(i,j)=u(i,j-1,k,nrhs)-2.*u(i,j,k,nrhs)+u(i,j+1,k,nrhs)
        enddo
      enddo
      if (south_bnd) then
        do i=1,nx
          uee(i,0)=uee(i,1)
        enddo
      endif
      if (north_bnd) then
        do i=1,nx
          uee(i,ny+1)=uee(i,ny)
        enddo
      endif

      do j=1,ny+1
        do i=0,nx
         Hvxx(i,j)=FlxV(i-1,j,k)-2.*FlxV(i,j,k)+FlxV(i+1,j,k)
        enddo
      enddo
      do j=1,ny+1
        do i=1,nx
          cff=FlxV(i,j,k)+FlxV(i-1,j,k)-delta*( Hvxx(i  ,j)
     &                                           +Hvxx(i-1,j))
          UFe(i,j)=0.25*( cff*(u(i,j,k,nrhs)+u(i,j-1,k,nrhs))
     &                          -gamma*( max(cff,0.)*uee(i,j-1)
     &                                  +min(cff,0.)*uee(i,j  )
     &                                                      ))
        enddo
      enddo
# undef Hvxx
# undef uee

# define vxx wrk1
# define Huee wrk2

      do j=1,ny
        do i=0,nx+1
          vxx(i,j)=v(i-1,j,k,nrhs)-2.*v(i,j,k,nrhs)+v(i+1,j,k,nrhs)
        enddo
      enddo
      if (west_bnd) then
        do j=1,ny
          vxx(0,j)=vxx(1,j)
        enddo
      endif
      if (east_bnd) then
        do j=1,ny
          vxx(nx+1,j)=vxx(nx,j)
        enddo
      endif
      do j=0,ny
        do i=1,nx+1
         Huee(i,j)=FlxU(i,j-1,k)-2.*FlxU(i,j,k)+FlxU(i,j+1,k)
        enddo
      enddo
      do j=1,ny
        do i=1,nx+1
          cff=FlxU(i,j,k)+FlxU(i,j-1,k)-delta*( Huee(i,j  )
     &                                           +Huee(i,j-1))
          VFx(i,j)=0.25*( cff*(v(i,j,k,nrhs)+v(i-1,j,k,nrhs))
     &                          -gamma*( max(cff,0.)*vxx(i-1,j)
     &                                  +min(cff,0.)*vxx(i  ,j)
     &                                                      ))
        enddo
      enddo
# undef Huee
# undef vxx
      do j=1,ny
        do i=1,nx
          ru(i,j,k)=ru(i,j,k)-UFx(i,j  )+UFx(i-1,j)
     &                         -UFe(i,j+1)+UFe(i  ,j)
        enddo
      enddo
      do j=1,ny
        do i=1,nx
          rv(i,j,k)=rv(i,j,k)-VFx(i+1,j)+VFx(i,j  )
     &                         -VFe(i  ,j)+VFe(i,j-1)
        enddo
      enddo

      end subroutine  uv_hadv_up3 !]
!----------------------------------------------------------------
      subroutine uv_hadv_c4(ru,rv,k) ![
! Add horizontal advection of momentum: compute diagonal [UFx,VFe]
! and off-diagonal [UFe,VFx] components of momentum flux tensor due
! to horizontal advection; after that add their divergence to r.h.s.
      !       coef_c4 = 1./12.  ! 1./16.

      ! Centered  4th order:           (-f(k+2) +9*f(k+1) + 9*f(k) -f(k-1))/16
      ! Centered 'modified' 4th order  (-f(k+2) +7*f(k+1) + 7*f(k) -f(k-1))/12
      ! modified when coef_c4 = 1/12, normal for coef_c4=1./16

      ! Import
      integer, intent(in) :: k
      real,dimension(GLOBAL_2D_ARRAY,nz) :: ru,rv
      ! local
      integer :: i,j
      real    :: cff

      do j=1,ny
        do i=0,nx+1
          wrk1(i,j)=u(i-1,j,k,nrhs)-2.*u(i,j,k,nrhs)+u(i+1,j,k,nrhs)
          wrk2(i,j)=FlxU(i-1,j,k) -2.*FlxU(i,j,k) +FlxU(i+1,j,k)
        enddo
      enddo

      if (west_bnd) then
        do j=1,ny
          wrk1(1,j) = wrk1(2,j)
          wrk2(1,j) = wrk2(2,j)
        enddo
      endif
      if (east_bnd) then
        do j=1,ny
          wrk1(nx+1,j) = wrk1(nx,j)
          wrk2(nx+1,j) = wrk2(nx,j)
        enddo
      endif
      do j=1,ny
        do i=0,nx
          UFx(i,j)=0.25*( u(i,j,k,nrhs)+u(i+1,j,k,nrhs)
     &                         -delta*(wrk1(i,j)+wrk1(i+1,j))
     &                  )*( FlxU(i,j,k)+FlxU(i+1,j,k)
     &                      -delta*(wrk2(i,j)+wrk2(i+1,j)))
        enddo
      enddo

# define vee wrk1
# define Hvee wrk2
      do j=0,ny+1
        do i=1,nx
          vee(i,j)=v(i,j-1,k,nrhs)-2.*v(i,j,k,nrhs)+v(i,j+1,k,nrhs)
          Hvee(i,j)=FlxV(i,j-1,k) -2.*FlxV(i,j,k)  +FlxV(i,j+1,k)
        enddo
      enddo

      if (south_bnd) then
        do i=1,nx
          vee(i,1) =vee(i,2)
          Hvee(i,1)=Hvee(i,2)
        enddo
      endif
      if (north_bnd) then
        do i=1,nx
          vee(i,ny+1)=vee(i,ny)
          Hvee(i,ny+1)=Hvee(i,ny)
        enddo
      endif

      do j=0,ny
        do i=1,nx
          VFe(i,j)=0.25*( v(i,j,k,nrhs)+v(i,j+1,k,nrhs)
     &                        -delta*(vee(i,j)+vee(i,j+1))
     &                   )*( FlxV(i,j,k)+FlxV(i,j+1,k)
     &                      -delta*(Hvee(i,j)+Hvee(i,j+1)))
        enddo
      enddo
# undef Hvee
# undef vee

# define uee wrk1
# define Hvxx wrk2
      do j=0,ny+1
        do i=1,nx
          uee(i,j)=u(i,j-1,k,nrhs)-2.*u(i,j,k,nrhs)+u(i,j+1,k,nrhs)
        enddo
      enddo
      if (south_bnd) then
        do i=1,nx
          uee(i,0)=uee(i,1)
        enddo
      endif
      if (north_bnd) then
        do i=1,nx
          uee(i,ny+1)=uee(i,ny)
        enddo
      endif

      do j=1,ny+1
        do i=0,nx
         Hvxx(i,j)=FlxV(i-1,j,k)-2.*FlxV(i,j,k)+FlxV(i+1,j,k)
        enddo
      enddo
      do j=1,ny+1
        do i=1,nx
          UFe(i,j)=0.25*( u(i,j,k,nrhs)+u(i,j-1,k,nrhs)
     &                        -delta*(uee(i,j)+uee(i,j-1))
     &                  )*( FlxV(i,j,k)+FlxV(i-1,j,k)
     &                     -delta*(Hvxx(i,j)+Hvxx(i-1,j)))
        enddo
      enddo
# undef Hvxx
# undef uee

# define vxx wrk1
# define Huee wrk2

      do j=1,ny
        do i=0,nx+1
          vxx(i,j)=v(i-1,j,k,nrhs)-2.*v(i,j,k,nrhs)+v(i+1,j,k,nrhs)
        enddo
      enddo
      if (west_bnd) then
        do j=1,ny
          vxx(0,j)=vxx(1,j)
        enddo
      endif
      if (east_bnd) then
        do j=1,ny
          vxx(nx+1,j)=vxx(nx,j)
        enddo
      endif
      do j=0,ny
        do i=1,nx+1
         Huee(i,j)=FlxU(i,j-1,k)-2.*FlxU(i,j,k)+FlxU(i,j+1,k)
        enddo
      enddo
      do j=1,ny
        do i=1,nx+1
          VFx(i,j)=0.25*( v(i,j,k,nrhs)+v(i-1,j,k,nrhs)
     &                        -delta*(vxx(i,j)+vxx(i-1,j))
     &                  )*( FlxU(i,j,k)+FlxU(i,j-1,k)
     &                     -delta*(Huee(i,j)+Huee(i,j-1)))
        enddo
      enddo
# undef Huee
# undef vxx
      do j=1,ny
        do i=1,nx
          ru(i,j,k)=ru(i,j,k)-UFx(i,j  )+UFx(i-1,j)
     &                         -UFe(i,j+1)+UFe(i  ,j)
        enddo
      enddo
      do j=1,ny
        do i=1,nx
          rv(i,j,k)=rv(i,j,k)-VFx(i+1,j)+VFx(i,j  )
     &                         -VFe(i  ,j)+VFe(i,j-1)
        enddo
      enddo

      end subroutine  uv_hadv_c4 !]
!----------------------------------------------------------------
      subroutine uv_vadv_up3(ru,rv,j) ![
!DIR$ ATTRIBUTES FORCEINLINE :: uv_vadv_up3

! This version contains up3 advection in the vertical direction
!       coef_up3 = 1./6.  ! 1./8.

       ! Upwind  3rd order:         :  ( 3*f(k+1) + 6*f(k) -f(k-1))/8 (we>0)
       ! 'Modified' upwind 3rd order:  ( 2*f(k+1) + 5*f(k) -f(k-1))/6 (we>0)
       ! we get the modified version if coef_up3 is set to 1/6 in advection.opt

      ! Import
      integer, intent(in) :: j
      real,dimension(GLOBAL_2D_ARRAY,nz) :: ru,rv
      ! local
      integer :: i,k
      real :: Flxw

!      u-component -------------------------------------------------------

        do k=1,nz-1
          do i = 1,nx
            FC(i,k)=(u(i,j,k+1,nrhs)-u(i,j,k,nrhs))
          enddo
        enddo
        do i=1,nx
          FC(i,nz) = 2*FC(i,nz-1)-FC(i,nz-2)
          FC(i,0 ) = 2*FC(i,1) -FC(i,2)
        enddo

        do k=1,nz
          do i=1,nx
            DC(i,k) = FC(i,k)-FC(i,k-1)
          enddo
        enddo

        do k=1,nz-1
          do i=1,nx
            Flxw = 0.5*(We(i,j,k)+We(i-1,j,k))
            FC(i,k)=0.5*(u(i,j,k,nrhs)+u(i,j,k+1,nrhs))
     &                                                *Flxw
     &                  - coef_up3*( DC(i,k  )*max(Flxw,0.)
     &                              +DC(i,k+1)*min(Flxw,0.))
          enddo
        enddo
!       flxw>0 :   3/6*(u(k+1) + u(k) )   1/6* ( u(k+1) - 2*u(k) - u(k-1) )

        FC(1:nx,0)  = 0
        FC(1:nx,nz) = 0

        do k=1,nz
          do i=1,nx
            ru(i,j,k)=ru(i,j,k)-FC(i,k)+FC(i,k-1)
          enddo
        enddo

!      v-component -------------------------------------------------------

        do k=1,nz-1
          do i = 1,nx
            FC(i,k)=(v(i,j,k+1,nrhs)-v(i,j,k,nrhs))
          enddo
        enddo
        do i=1,nx
          FC(i,nz) = 2*FC(i,nz-1)-FC(i,nz-2)
          FC(i,0 ) = 2*FC(i,1) -FC(i,2)
        enddo

        do k=1,nz
          do i=1,nx
            crv(i,k) = FC(i,k)-FC(i,k-1)
             !FC(i,k) = FC(i,k+1)-FC(i,k)
          enddo
        enddo

        do k=1,nz-1
          do i=1,nx
            Flxw = 0.5*(We(i,j,k)+We(i,j-1,k))
            FC(i,k)=0.5*(v(i,j,k,nrhs)+v(i,j,k+1,nrhs))
     &                                                *Flxw
     &                 - coef_up3*( crv(i,k  )*max(Flxw,0.)
     &                             +crv(i,k+1)*min(Flxw,0.))
          enddo
        enddo
        FC(1:nx,0)  = 0
        FC(1:nx,nz) = 0

        do k=1,nz
          do i=1,nx
            rv(i,j,k)=rv(i,j,k)-FC(i,k)+FC(i,k-1)
          enddo
        enddo
      end subroutine uv_vadv_up3 !]
!----------------------------------------------------------------
      subroutine uv_vadv_c4(ru,rv,j) ![
      ! This version contains c4 advection in the vertical direction
      !       coef_c4 = 1./12.  ! 1./16.

      ! Centered  4th order:           (-f(k+2) +9*f(k+1) + 9*f(k) -f(k-1))/16
      ! Centered 'modified' 4th order  (-f(k+2) +7*f(k+1) + 7*f(k) -f(k-1))/12
      ! modified when coef_c4 = 1/12, normal for coef_c4=1./16

      ! Import
      integer, intent(in) :: j
      real,dimension(GLOBAL_2D_ARRAY,nz) :: ru,rv
      ! local
      integer :: i,k

!     u-component -------------------------------------------------------

      do k= 2,nz-1
        do i=1,nx
          FC(i,k)= u(i,j,k+1,nrhs)-u(i,j,k-1,nrhs)
        enddo
      enddo
      do i=1,nx
        FC(i,1 ) = u(i,j,2 ,nrhs)-u(i,j,1 ,nrhs)
        FC(i,nz) = u(i,j,nz,nrhs)-u(i,j,nz-1,nrhs)
      enddo

      do k=1,nz-1
        do i=1,nx
          FC(i,k)= 0.5*(We(i,j,k)+We(i-1,j,k))*(
     &               0.5*(u(i,j,k,nrhs)+u(i,j,k+1,nrhs))
     &              -coef_c4*(FC(i,k+1)-FC(i,k)) )
        enddo
      enddo

      FC(1:nx,0)  = 0
      FC(1:nx,nz) = 0

      do k=1,nz
        do i=1,nx
          ru(i,j,k)=ru(i,j,k)-FC(i,k)+FC(i,k-1)
        enddo
      enddo

!     v-component -------------------------------------------------------

      do k= 2,nz-1
        do i=1,nx
          FC(i,k)= v(i,j,k+1,nrhs)-v(i,j,k-1,nrhs)
        enddo
      enddo
      do i=1,nx
        FC(i,1 ) = v(i,j,2 ,nrhs)-v(i,j,1 ,nrhs)
        FC(i,nz) = v(i,j,nz,nrhs)-v(i,j,nz-1,nrhs)
      enddo
      do k=1,nz-2
        do i=1,nx
          FC(i,k)= 0.5*(We(i,j,k)+We(i-1,j,k))*(
     &               0.5*(v(i,j,k,nrhs)+v(i,j,k+1,nrhs))
     &              -coef_c4*(FC(i,k+1)-FC(i,k)) )
        enddo
      enddo

      FC(1:nx,0)  = 0
      FC(1:nx,nz) = 0

      do k=1,nz
        do i=1,nx
          rv(i,j,k)=rv(i,j,k)-FC(i,k)+FC(i,k-1)
        enddo
      enddo

      end subroutine uv_vadv_c4 !]
!----------------------------------------------------------------
      subroutine uv_vadv_spline(ru,rv,j) ![
!DIR$ ATTRIBUTES FORCEINLINE :: uv_vadv_spline

      ! Vertical advection using a spline algorithm

      ! Import
      integer, intent(in) :: j
      real,dimension(GLOBAL_2D_ARRAY,nz) :: ru,rv
      ! local
      integer :: i,k
      real :: cff

        do i=1,nx
          DC(i,1)=0.5625*(Hz(i  ,j,1)+Hz(i-1,j,1))
     &           -0.0625*(Hz(i+1,j,1)+Hz(i-2,j,1))
          CF(i,1)=1.  ;  FC(i,0)=2.0*u(i,j,1,nrhs)
        enddo
        do k=1,nz-1,+1    !--> recursive
          do i=1,nx
            DC(i,k+1)=0.5625*(Hz(i  ,j,k+1)+Hz(i-1,j,k+1))
     &               -0.0625*(Hz(i+1,j,k+1)+Hz(i-2,j,k+1))

            cff=1./(2.*DC(i,k)+DC(i,k+1)*(2.-CF(i,k)))
            CF(i,k+1)=cff*DC(i,k)
            FC(i,k)=cff*( 3.*( DC(i,k  )*u(i,j,k+1,nrhs)
     &                        +DC(i,k+1)*u(i,j,k  ,nrhs))
     &                              -DC(i,k+1)*FC(i,k-1))
          enddo
        enddo               !--> discard DC, keep CF,FC
        do i=1,nx
          FC(i,nz)=(2.*u(i,j,nz,nrhs)-FC(i,nz-1))/(1.-CF(i,nz))
          DC(i,nz)=0.        !<-- uppermost W*U flux
        enddo
        do k=nz-1,1,-1       !--> recursive
          do i=1,nx
            FC(i,k)=FC(i,k)-CF(i,k+1)*FC(i,k+1)

            DC(i,k)=FC(i,k) * 0.5*( We(i,j,k)+We(i-1,j,k) -0.125*(
     &                        (We(i+1,j,k)-We(i  ,j,k))*umask(i+1,j)
     &                       -(We(i-1,j,k)-We(i-2,j,k))*umask(i-1,j)
     &                                                          ))
            ru(i,j,k+1)=ru(i,j,k+1) -DC(i,k+1)+DC(i,k)
          enddo
        enddo                       !--> discard CF,FC
        do i=1,nx
          ru(i,j,1)=ru(i,j,1) -DC(i,1)
        enddo                          !--> discard DC

          do i=1,nx
            DC(i,1)=0.5625*(Hz(i  ,j,1)+Hz(i,j-1,1))
     &             -0.0625*(Hz(i,j+1,1)+Hz(i,j-2,1))
            CF(i,1)=1.  ;  FC(i,0)=2.0*v(i,j,1,nrhs)
          enddo
          do k=1,nz-1,+1       !--> recursive
            do i=1,nx
              DC(i,k+1)=0.5625*(Hz(i  ,j,k+1)+Hz(i,j-1,k+1))
     &                 -0.0625*(Hz(i,j+1,k+1)+Hz(i,j-2,k+1))

              cff=1./(2.*DC(i,k)+DC(i,k+1)*(2.-CF(i,k)))
              CF(i,k+1)=cff*DC(i,k)
              FC(i,k)=cff*( 3.*( DC(i,k  )*v(i,j,k+1,nrhs)
     &                          +DC(i,k+1)*v(i,j,k  ,nrhs))
     &                                -DC(i,k+1)*FC(i,k-1))
            enddo
          enddo               !--> discard DC, keep CF,FC
          do i=1,nx
            FC(i,nz)=(2.*v(i,j,nz,nrhs)-FC(i,nz-1))/(1.-CF(i,nz))
            DC(i,nz)=0.        !<-- uppermost W*V flux
          enddo
          do k=nz-1,1,-1       !--> recursive
            do i=1,nx
              FC(i,k)=FC(i,k)-CF(i,k+1)*FC(i,k+1)

              DC(i,k)=FC(i,k) * 0.5*( We(i,j,k)+We(i,j-1,k) -0.125*(
     &                         (We(i,j+1,k)-We(i,j  ,k))*vmask(i,j+1)
     &                        -(We(i,j-1,k)-We(i,j-2,k))*vmask(i,j-1)
     &                                                           ))
              rv(i,j,k+1)=rv(i,j,k+1) -DC(i,k+1)+DC(i,k)
            enddo
          enddo               !--> discard CF,FC
          do i=1,nx
            rv(i,j,1)=rv(i,j,1) -DC(i,1)
          enddo                         !--> discard DC

      end subroutine uv_vadv_spline !]
!----------------------------------------------------------------
      subroutine t_hadv_up3(Fx,Fy,k,itrc) ![
!DIR$ ATTRIBUTES FORCEINLINE :: t_hadv_up3

      ! Horizontal tracer advection fluxes using a upwind 3rd order algorithm

      ! Import
      integer, intent(in) :: k,itrc
      real,dimension(GLOBAL_2D_ARRAY) :: Fx,Fy
      ! local
      integer :: i,j

      do j=1,ny
        do i=0,nx+2
          Fx(i,j)=(t(i,j,k,nrhs,itrc)-t(i-1,j,k,nrhs,itrc))*umask(i,j)
        enddo
      enddo
      if (west_bnd) then
        do j=1,ny
          Fx(0,j)=Fx(1,j)
        enddo
      endif
      if (east_bnd) then
        do j=1,ny
          Fx(nx+2,j)=Fx(nx+1,j)
        enddo
      endif
      do j=1,ny
        do i=0,nx+1
          wrk1(i,j)=Fx(i+1,j)-Fx(i,j)
        enddo
      enddo
      do j=1,ny
        do i=1,nx+1
          Fx(i,j)=0.5*(t(i,j,k,nrhs,itrc)+t(i-1,j,k,nrhs,itrc))*FlxU(i,j,k)
     &                    -coef_up3*( wrk1(i-1,j)*max(FlxU(i,j,k),0.)
     &                               +wrk1(i  ,j)*min(FlxU(i,j,k),0.))
        enddo
      enddo

      do j=0,ny+2
        do i=1,nx
          Fy(i,j)=(t(i,j,k,nrhs,itrc)-t(i,j-1,k,nrhs,itrc))*vmask(i,j)
        enddo
      enddo
      if (south_bnd) then
        do i=1,nx
          Fy(i,0)=Fy(i,1)
        enddo
      endif
      if (north_bnd) then
        do i=1,nx
          Fy(i,ny+2)=Fy(i,ny+1)
        enddo
      endif
      do j=0,ny+1
        do i=1,nx
          wrk1(i,j)=Fy(i,j+1)-Fy(i,j)
        enddo
      enddo
      do j=1,ny+1
        do i=1,nx
          Fy(i,j)=0.5*(t(i,j,k,nrhs,itrc)+t(i,j-1,k,nrhs,itrc))*FlxV(i,j,k)
     &                    -coef_up3*( wrk1(i,j-1)*max(FlxV(i,j,k),0.)
     &                               +wrk1(i,j  )*min(FlxV(i,j,k),0.))
        enddo
      enddo

      end subroutine t_hadv_up3 !]
!----------------------------------------------------------------
      subroutine t_hadv_c4(Fx,Fy,k,itrc) ![
!DIR$ ATTRIBUTES FORCEINLINE :: t_hadv_c4

      ! Horizontal tracer advection fluxes using a centered 4th order algorithm

      ! Import
      integer, intent(in) :: k,itrc
      real,dimension(GLOBAL_2D_ARRAY) :: Fx,Fy
      ! local
      integer :: i,j

      do j=1,ny
        do i=0,nx+2
          Fx(i,j)=(t(i,j,k,nrhs,itrc)-t(i-1,j,k,nrhs,itrc))*umask(i,j)
        enddo
      enddo
      if (west_bnd) then
        do j=1,ny
          Fx(0,j)=Fx(1,j)
        enddo
      endif
      if (east_bnd) then
        do j=1,ny
          Fx(nx+2,j)=Fx(nx+1,j)
        enddo
      endif

      do j=1,ny
        do i=0,nx+1
          wrk1(i,j)=0.5*(Fx(i+1,j)+Fx(i,j))
        enddo
      enddo
      do j=1,ny
        do i=1,nx+1
          Fx(i,j)=0.5*( t(i,j,k,nrhs,itrc)+t(i-1,j,k,nrhs,itrc)
     &                -0.3333333333333333*(wrk1(i,j)-wrk1(i-1,j))
     &                                              )*FlxU(i,j,k)
        enddo
      enddo
      do j=0,ny+2
        do i=1,nx
          Fy(i,j)=(t(i,j,k,nrhs,itrc)-t(i,j-1,k,nrhs,itrc))*vmask(i,j)
        enddo
      enddo
      if (south_bnd) then
        do i=1,nx
          Fy(i,0)=Fy(i,1)
        enddo
      endif
      if (north_bnd) then
        do i=1,nx
          Fy(i,ny+2)=Fy(i,ny+1)
        enddo
      endif

      do j=0,ny+1
        do i=1,nx
          wrk1(i,j)=0.5*(Fy(i,j+1)+Fy(i,j))
        enddo
      enddo
      do j=1,ny+1
        do i=1,nx
          Fy(i,j)=0.5*( t(i,j,k,nrhs,itrc)+t(i,j-1,k,nrhs,itrc)
     &               -0.3333333333333333*(wrk1(i,j)-wrk1(i,j-1))
     &                                             )*FlxV(i,j,k)
        enddo
      enddo

      end subroutine t_hadv_c4 !]
!----------------------------------------------------------------
      subroutine t_vadv_spline(Fz,j,itrc) ![
!DIR$ ATTRIBUTES FORCEINLINE :: t_vadv_spline

      ! Vertical tracer advection fluxes using a spline algorithm


      ! Import
      integer, intent(in) :: j,itrc
      real,dimension(1-bf:nx+bf,0:nz) :: Fz
      ! local
      integer :: i,k
      real :: cff

      do i=1,nx
        CF(i,1)=1.
        Fz(i,0)=2.0*t(i,j,1,nrhs,itrc)
      enddo
      do k=1,nz-1,+1    !--> recursive
        do i=1,nx
          cff=1./(2.*Hz(i,j,k)+Hz(i,j,k+1)*(2.-CF(i,k)))
          CF(i,k+1)=cff*Hz(i,j,k)
          Fz(i,k)=cff*( 3.*( Hz(i,j,k  )*t(i,j,k+1,nrhs,itrc)
     &                       +Hz(i,j,k+1)*t(i,j,k  ,nrhs,itrc))
     &                                  -Hz(i,j,k+1)*Fz(i,k-1))
        enddo
      enddo
      do i=1,nx
        Fz(i,nz)=(2.*t(i,j,nz,nrhs,itrc)-Fz(i,nz-1))/(1.-CF(i,nz))
      enddo
      do k=nz-1,0,-1    !<-- recursive
        do i=1,nx
          Fz(i,k)=Fz(i,k)-CF(i,k+1)*Fz(i,k+1)

          Fz(i,k+1)=Fz(i,k+1)*We(i,j,k+1)  ! Convert interface value
        enddo                               ! into vertical flux
      enddo
      do i=1,nx
        Fz(i,nz)=0.                         ! Set top and bottom
        Fz(i,0)=0.                          ! boundary conditions.
      enddo

      end subroutine t_vadv_spline !]
!----------------------------------------------------------------
      subroutine t_vadv_up3(Fz,j,itrc) ![
!DIR$ ATTRIBUTES FORCEINLINE :: t_vadv_up3

      ! Vertical tracer advection fluxes using a upwind 3rd order algorithm

      ! Import
      integer, intent(in) :: j,itrc
      real,dimension(1-bf:nx+bf,0:nz) :: Fz
      ! local
      integer :: i,k

      do k=1,nz-1
        do i = 1,nx
          FC(i,k)=(t(i,j,k+1,nrhs,itrc)-t(i,j,k,nrhs,itrc))
        enddo
      enddo
      do i=1,nx
        FC(i,nz)= 2*FC(i,nz-1)-FC(i,nz-2)
        FC(i,0) = 0
      enddo

      do k=1,nz
        do i=1,nx
          DC(i,k) = FC(i,k)-FC(i,k-1)
        enddo
      enddo

      do k= 1,nz-1
        do i=1,nx
          Fz(i,k)=0.5*(t(i,j,k,nrhs,itrc)+t(i,j,k+1,nrhs,itrc))
     &                                                  *We(i,j,k)
     &                    -coef_up3*( DC(i,k  )*max(We(i,j,k),0.)
     &                               +DC(i,k+1)*min(We(i,j,k),0.))
        enddo
      enddo

      Fz(1:nx,0 ) = 0
      Fz(1:nx,nz) = 0

      end subroutine t_vadv_up3 !]
!----------------------------------------------------------------
      subroutine t_vadv_c4(Fz,j,itrc) ![
!DIR$ ATTRIBUTES FORCEINLINE :: t_vadv_c4

      ! Vertical tracer advection fluxes using a centered 4th order algorithm

      ! Import
      integer, intent(in) :: j,itrc
      real,dimension(1-bf:nx+bf,0:nz) :: Fz
      ! local
      integer :: i,k

      do k=2,nz-1
        do i = 1,nx
          FC(i,k)= t(i,j,k+1,nrhs,itrc)-t(i,j,k-1,nrhs,itrc)
        enddo
      enddo
      do i=1,nx
        FC(i,1) = t(i,j,2 ,nrhs,itrc)-t(i,j,1,nrhs,itrc)
        FC(i,nz)= t(i,j,nz,nrhs,itrc)-t(i,j,nz-1,nrhs,itrc)
      enddo

      do k=1,nz-1
        do i=1,nx
          Fz(i,k)=  We(i,j,k)*(
     &         0.5*(t(i,j,k,nrhs,itrc)+t(i,j,k+1,nrhs,itrc))
     &                  -coef_c4*(FC(i,k+1)-FC(i,k)) )
        enddo
      enddo

      Fz(1:nx,0)  = 0
      Fz(1:nx,nz) = 0

      end subroutine t_vadv_c4 !]
!----------------------------------------------------------------
      end module advection
