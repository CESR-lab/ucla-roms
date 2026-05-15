      module advection

      ! contains various advection algorithms 


      use dimensions
      use diagnostics
      use ocean_vars
      use grid

      implicit none

      private

#include "advection.opt"

      real,parameter :: delta = 0.1666666666666667 !! delta=0.125
      real,parameter :: gamma = 0.3333333333333333 !! gamma=0.25

      real,allocatable,dimension(:,:) :: wrk1
      real,allocatable,dimension(:,:) :: wrk2
      real,allocatable,dimension(:,:) :: UFx
      real,allocatable,dimension(:,:) :: UFe
      real,allocatable,dimension(:,:) :: VFx
      real,allocatable,dimension(:,:) :: VFe

      public init_advection
      public uv_cor_curv
      public uv_hor_adv

!DIR$ ATTRIBUTES INLINE :: uv_cor_curv, uv_hor_adv

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

      end subroutine init_advection !]
!------------------------------------------------------------------------------
      subroutine uv_cor_curv(ru,rv,k) ![
      ! Add Coriolis terms and contribution to advection that is
      ! associated with curvilinear horizontal coordinates.

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
      subroutine uv_hor_adv(ru,rv,k) ![
! Add horizontal advection of momentum: compute diagonal [UFx,VFe]
! and off-diagonal [UFe,VFx] components of momentum flux tensor due
! to horizontal advection; after that add their divergence to r.h.s.

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
# ifdef UPSTREAM_UV
          cff=FlxU(i,j,k)+FlxU(i+1,j,k)-delta*( wrk2(i  ,j)
     &                                            +wrk2(i+1,j))
          UFx(i,j)=0.25*( cff*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))
     &                          -gamma*( max(cff,0.)*wrk1(i  ,j)
     &                                  +min(cff,0.)*wrk1(i+1,j)
     &                                                      ))
# else
          UFx(i,j)=0.25*( u(i,j,k,nrhs)+u(i+1,j,k,nrhs)
     &                         -delta*(wrk1(i,j)+wrk1(i+1,j))
     &                  )*( FlxU(i,j,k)+FlxU(i+1,j,k)
     &                      -delta*(wrk2(i,j)+wrk2(i+1,j)))
# endif
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
# ifdef UPSTREAM_UV
          cff=FlxV(i,j,k)+FlxV(i,j+1,k)-delta*( Hvee(i,j  )
     &                                           +Hvee(i,j+1))
          VFe(i,j)=0.25*( cff*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))
     &                          -gamma*( max(cff,0.)*vee(i,j  )
     &                                  +min(cff,0.)*vee(i,j+1)
     &                                                      ))
# else
          VFe(i,j)=0.25*( v(i,j,k,nrhs)+v(i,j+1,k,nrhs)
     &                        -delta*(vee(i,j)+vee(i,j+1))
     &                   )*( FlxV(i,j,k)+FlxV(i,j+1,k)
     &                      -delta*(Hvee(i,j)+Hvee(i,j+1)))
# endif
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
# ifdef UPSTREAM_UV
          cff=FlxV(i,j,k)+FlxV(i-1,j,k)-delta*( Hvxx(i  ,j)
     &                                           +Hvxx(i-1,j))
          UFe(i,j)=0.25*( cff*(u(i,j,k,nrhs)+u(i,j-1,k,nrhs))
     &                          -gamma*( max(cff,0.)*uee(i,j-1)
     &                                  +min(cff,0.)*uee(i,j  )
     &                                                      ))
# else
          UFe(i,j)=0.25*( u(i,j,k,nrhs)+u(i,j-1,k,nrhs)
     &                        -delta*(uee(i,j)+uee(i,j-1))
     &                  )*( FlxV(i,j,k)+FlxV(i-1,j,k)
     &                     -delta*(Hvxx(i,j)+Hvxx(i-1,j)))
# endif
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
# ifdef UPSTREAM_UV
          cff=FlxU(i,j,k)+FlxU(i,j-1,k)-delta*( Huee(i,j  )
     &                                           +Huee(i,j-1))
          VFx(i,j)=0.25*( cff*(v(i,j,k,nrhs)+v(i-1,j,k,nrhs))
     &                          -gamma*( max(cff,0.)*vxx(i-1,j)
     &                                  +min(cff,0.)*vxx(i  ,j)
     &                                                      ))
# else
          VFx(i,j)=0.25*( v(i,j,k,nrhs)+v(i-1,j,k,nrhs)
     &                        -delta*(vxx(i,j)+vxx(i-1,j))
     &                  )*( FlxU(i,j,k)+FlxU(i,j-1,k)
     &                     -delta*(Huee(i,j)+Huee(i,j-1)))
# endif
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

      end subroutine  uv_hor_adv !]
      !----------------------------------------------------------------

      end module advection
