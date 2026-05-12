#include "cppdefs.opt"
#if defined SOLVE3D && defined LMD_MIXING
# define SMOOTH_RIG  /* local Ri number smoothing */

      subroutine lmd_vmix(tind)
      use param
      use roms_mpi
      use private_scratch
      implicit none
      integer,save :: tile=0
      integer      :: tind

# include "compute_tile_bounds.h"
      call lmd_vmix_tile (istr,iend,jstr,jend,  A3d(1,1), A3d(1,2),
     &                                           A3d(1,3), A3d(1,4),
     &                                 A2d(1,1), A2d(1,2), A2d(1,3),
     &                                                         tind)
# if defined LMD_KPP || defined LMD_BKPP
      call lmd_kpp_tile  (istr,iend,jstr,jend,  A3d(1,1),  A3d(1, 2),
     &                                                     A3d(1, 3),

     &                    A2d(1, 1), A2d(1, 2), A2d(1, 3), A2d(1, 4),
     &                    A2d(1, 5), A2d(1, 6), A2d(1, 7), A2d(1, 8),
     &                    A2d(1, 9), A2d(1,10), A2d(1,11), A2d(1,12),
     &                    A2d(1,13), A2d(1,14), A2d(1,15), A2d(1,16),
     &                    A2d(1,17), A2d(1,18), iA2d(1,1), iA2d(1,2),
     &                                                          tind)
# endif
      end

      subroutine lmd_vmix_tile (istr,iend,jstr,jend, Kv,Kt,Ks, Rig,
     &                                              FX,FE,FE1,tind)

! Compute vertical mixing coefficients for momentum and tracers at the
! ocean interior using the Large, McWilliams, and Doney (1994) mixing
! parameterization.
!
! output:   Kv       vertical viscosity coefficient [m^2/s]
!           Kt       vertical heat conductivity coefficient
!                          for potential temperature [m^2/s]
!           Ks       vertical diffusion for salinity [m^2/s]
!
!  Reference: Large, W. G., J. C. McWilliams, and S. C. Doney, 1994:
!             A Review and model with a nonlocal boundary layer
!             parameterization.  Reviews of Geophysics, 32,363-403.

      use param
      use tracers, only: t
      use grid
      use mixing
      use ocean_vars
      use scalars
      use roms_mpi

      implicit none

      integer  :: tind
      integer istr,iend,jstr,jend, i,j,k
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,0:N) :: Kv,Kt,Ks, Rig
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: FX,FE,FE1

      real nu_sx, cff,cff1, dudz,dvdz
      real, parameter ::     ! Critical gradient Richardson number
     &      Ri0=0.7,         ! below which shear instabilty occurs.

     &      nu0m= 1.e-2,     ! Maximum viscosity and diffusivity
     &      nu0s= 1.e-2,     ! due to shear instability [m^2/s];

     &      nuwm=1.0e-4,     ! Interior viscosity and diffusivity
     &      nuws=0.1e-4,     ! due to wave breaking, [m^2/s]

     &      nu0c=0.1,        ! Convective adjustment [m^2/s];

     &      lmd_nu=1.5e-6,   ! Molecular viscosity [m^2/s];

                             ! Value of double-diffusive density
     &      lmd_Rrho0=1.9,   ! ratio where diffusivity goes to zero
                             !                     in salt fingering.
     &      lmd_nuf=10.0e-4, ! Scaling factors for double diffusion
     &      lmd_fdd=0.7,     ! coefficient in salt fingering.

     &      lmd_tdd1=0.909,
     &      lmd_tdd2=4.6,    ! Double diffusion constants for
     &      lmd_tdd3=0.54,   ! temperature (Marmorino and Caldwell,
     &      lmd_sdd1=0.15,   ! 1976) and salinity (Fedorov, 1988)
     &      lmd_sdd2=1.85,
     &      lmd_sdd3=0.85,
     &      eps=1.E-14,
     &      Lturb = 10.      ! Length scale of bottom boundary layer 
                             ! where turbulence is supressed

      real dist,mult

# ifdef LMD_RIMIX
c----#  define tind nstp
#  ifdef SMOOTH_RIG
#   ifdef EW_PERIODIC
#    define I_EXT_RANGE istr-1,iend+1
#   else
#    define I_EXT_RANGE imin,imax
      integer imin,imax
#   endif
#   ifdef NS_PERIODIC
#    define J_EXT_RANGE jstr-1,jend+1
#   else
#    define J_EXT_RANGE jmin,jmax
      integer jmin,jmax
#   endif
#   ifndef EW_PERIODIC
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
#   endif
#   ifndef NS_PERIODIC
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
#   endif
#  else
#   define I_EXT_RANGE istr,iend
#   define J_EXT_RANGE jstr,jend
#  endif  /* SMOOTH_RIG */
# endif  /* LMD_RIMIX */

!                                    ! Compute horizontal velocity
! Compute local Richardson number:   ! shear  (du/dz)^2+(dv/dz)^2  at
!-------- ----- ---------- -------   ! horizontal RHO- and vertical
!                                    ! W-points. Then compute gradient
      do k=1,N-1                     ! Richardson number (already
# ifdef LMD_RIMIX
        do j=J_EXT_RANGE             ! divided by its critical value.
          do i=I_EXT_RANGE
            cff=0.5/(z_r(i,j,k+1)-z_r(i,j,k))
            dudz=cff*( u(i  ,j,k+1,tind)-u(i  ,j,k,tind)
     &                +u(i+1,j,k+1,tind)-u(i+1,j,k,tind))
            dvdz=cff*( v(i,j  ,k+1,tind)-v(i,j  ,k,tind)
     &                +v(i,j+1,k+1,tind)-v(i,j+1,k,tind))

            Rig(i,j,k)=bvf(i,j,k)/( Ri0*max(
     &                    dudz*dudz+dvdz*dvdz, 1.E-10 ))
          enddo
        enddo

#  ifdef SMOOTH_RIG
#   ifndef EW_PERIODIC
        if (WESTERN_EDGE) then
          do j=J_EXT_RANGE
            Rig(istr-1,j,k)=Rig(istr,j,k)
          enddo
        endif
        if (EASTERN_EDGE) then
          do j=J_EXT_RANGE
            Rig(iend+1,j,k)=Rig(iend,j,k)
          enddo
        endif
#   endif
#   ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do i=I_EXT_RANGE
            Rig(i,jstr-1,k)=Rig(i,jstr,k)
          enddo
        endif
        if (NORTHERN_EDGE) then
          do i=I_EXT_RANGE
            Rig(i,jend+1,k)=Rig(i,jend,k)
          enddo
        endif
#    ifndef EW_PERIODIC
        if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
          Rig(istr-1,jstr-1,k)=Rig(istr,jstr,k)
        endif
        if (WESTERN_EDGE .and. NORTHERN_EDGE) then
          Rig(istr-1,jend+1,k)=Rig(istr,jend,k)
        endif
        if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
          Rig(iend+1,jstr-1,k)=Rig(iend,jstr,k)
        endif
        if (EASTERN_EDGE .and. NORTHERN_EDGE) then
          Rig(iend+1,jend+1,k)=Rig(iend,jend,k)
        endif
#    endif
#   endif

       cff=1./12.
       cff1=3./16.
                                               ! Isotropic smoothing by
        do j=jstr-1,jend+1                     ! first computing masked
          do i=istr,iend+1                     ! elementary differences
            FX(i,j)=(Rig(i,j,k)-Rig(i-1,j,k))  ! in each direction,
#  ifdef MASKING
     &                            *umask(i,j)  ! after which expressing
#  endif
          enddo                                ! the transversal terms
        enddo                                  ! via the very same
        do j=jstr,jend+1                       ! masked differences.
          do i=istr-1,iend+1
            FE1(i,j)=(Rig(i,j,k)-Rig(i,j-1,k))
#  ifdef MASKING
     &                             *vmask(i,j)
#  endif
          enddo
          do i=istr,iend
            FE(i,j)=FE1(i,j) + cff*( FX(i+1,j)+FX(i  ,j-1)
     &                              -FX(i  ,j)-FX(i+1,j-1))
          enddo
        enddo
        do j=jstr,jend
          do i=istr,iend+1
            FX(i,j)=FX(i,j) + cff*( FE1(i,j+1)+FE1(i-1,j  )
     &                             -FE1(i,j  )-FE1(i-1,j+1))
          enddo
          do i=istr,iend
            Rig(i,j,k)=Rig(i,j,k) + cff1*( FX(i+1,j)-FX(i,j)
     &                                    +FE(i,j+1)-FE(i,j))
          enddo
        enddo              !--> discard FX,FE,FE1

#  endif /* SMOOTH_RIG */
#  undef I_EXT_RANGE
#  undef J_EXT_RANGE
# endif /* LMD_RIMIX */

! Compute "interior" viscosities and diffusivities everywhere as
! superposition of three processes: local Richardson number instability
! due to resolved vertical shear, internal wave breaking,
! and double diffusion.

        do j=jstr,jend
          do i=istr,iend
# ifdef LMD_RIMIX
            cff=min(1., max(0., Rig(i,j,k)))  ! Compute mixing due
            nu_sx=1. - cff*cff                ! to shear instability
            nu_sx=nu_sx*nu_sx*nu_sx           ! and internal wave
                                              ! breaking.
            Kv(i,j,k)=nuwm + nu0m*nu_sx
            Kt(i,j,k)=nuws + nu0s*nu_sx
# else
            Kv(i,j,k)=nuwm                    ! internal waves
            Kt(i,j,k)=nuws                    ! only
# endif
            ! Convective adjustment will now be added everywhere
            ! but 'tailed' near the bottom
            ! It is not generally recommended to use...
# if defined LMD_CONVEC 
            if (Rig(i,j,k) < 0.) then
              Kv(i,j,k)=Kv(i,j,k) + nu0c
              Kt(i,j,k)=Kt(i,j,k) + nu0c
            endif
# endif
            Ks(i,j,k)=Kt(i,j,k)
          enddo
        enddo

      enddo     ! <-- k

! Supress mixing away near the bottom. 
      do k=1,nz-1
        do j=1,ny
          do i=1,nx
            dist = z_w(i,j,k)-z_w(i,j,0)
            if (dist<Lturb) then
              mult = sin(0.5*pi*(z_w(i,j,k) - z_w(i,j,0))/Lturb)
              Kv(i,j,k) = Kv(i,j,k)*mult
              Kt(i,j,k) = Kt(i,j,k)*mult
            endif
          enddo
        enddo
      enddo

      do j=1,ny
        do i=1,nx
          Kv(i,j,nz)=Kv(i,j,nz-1)
          Kt(i,j,nz)=Kt(i,j,nz-1)
          Kv(i,j, 0)=Kv(i,j,   1)
          Kt(i,j, 0)=Kt(i,j,   1)
        enddo
      enddo

# if defined LMD_KPP || defined LMD_BKPP
      ! vertical smoothing of interior mixing
      do k=1,nz-1
        do j=jstr,jend
          do i=istr,iend
            Kv(i,j,k)=0.5*Kv(i,j,k)+0.25*Kv(i,j,k-1)+0.25*Kv(i,j,k+1)
            Kt(i,j,k)=0.5*Kt(i,j,k)+0.25*Kt(i,j,k-1)+0.25*Kt(i,j,k+1)
          enddo
        enddo
      enddo
# else /* no kpp at all */

! Finalize: Copy everything into shared arrays:

      do k=1,N-1
        do j=jstr,jend
          do i=istr,iend
            AKv(i,j,k) = 0.5*Kv(i,j,k)+0.25*Kv(i,j,k-1)+0.25*Kv(i,j,k+1)
            AKt(i,j,k) = 0.5*Kt(i,j,k)+0.25*Kt(i,j,k-1)+0.25*Kt(i,j,k+1)
          enddo
        enddo
      enddo

      !! Do we really need to exchange the mixing coefficients?
      call exchange_xxx(Akv,Akt)

# endif /* LMD_KPP || LMD_BKPP */
      end
#else
      subroutine lmd_vmix_empty
      end
#endif /* LMD_MIXING */
