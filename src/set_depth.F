#include "cppdefs.opt"

      module zstar_geom  ![
      ! zeta* support (design D, minimal): zflt/zflt_bak hold the two
      ! most recent fast-filtered zetas (recorded by set_depth at the
      ! end of each fast window); zeta* = 1.5*zflt - 0.5*zflt_bak is
      ! centered at n+1/2 of the upcoming fast window.  The corrector
      ! prsgrd shifts its surface pressure by g*(zeta*-z_w(:,:,nz));
      ! step3d_uv1 removes the exact depth integral of that shift
      ! from rufrc/rvfrc, so the barotropic forcing is unchanged
      ! identically while the divergence seen by the NH projection
      ! carries the n+1/2 free surface.
      ! (see Documentation/readme-specifics/README.time_levels)
      use param
      use dimensions
      implicit none

      real,allocatable,dimension(:,:) :: zflt, zflt_bak, zstar
      logical :: zstar_init = .false.

      end module zstar_geom  !]



      subroutine set_depth ![

! Compute evolving z-coordinates "z_r" and "z_w" and the heights
! of vertical grid boxes "Hz" from fast-time-averaged free surface,
! topography h(i,j), and vertical coordinate transformation formula.

      use param
      use wec_frc
      use coupling
      use grid
      use ocean_vars
      use scalars
      use scoord
      use roms_mpi
      use zstar_geom

      implicit none
      integer ::  i,j,k
      real cff_r,cff1_r, cff_w,cff1_w, ds

                                       ! During initialization or
      if (iic == 0) then               ! restart compute barotropic
        do j=j0,j1                     ! fluxes and place them into
          do i=i0,i1                   ! arrays for fast-time-averaged
            hinv(i,j)=1./(h(i,j)+hc)   ! fluxes. Later during the first
          enddo                        ! procedure uses these DU_avg1,
        enddo                          ! DV_avg1 as input to compute
                                       ! ubar,vbar as it does so during
        do j=j0,j1                    ! all subsequent steps.
          do i=1,i1
            DU_avg1(i,j)=0.5*( h(i-1,j)+h(i,j)+zeta(i,j,1)
     &            +zeta(i-1,j,1) )*dn_u(i,j)*( ubar(i,j,1)
# ifdef WEC
     &                                        + ust2d(i,j)
# endif
     &                                                   )
          enddo
        enddo

        do j=1,j1
          do i=i0,i1
            DV_avg1(i,j)=0.5*( h(i,j)+h(i,j-1)+zeta(i,j,1)
     &            +zeta(i,j-1,1) )*dm_v(i,j)*( vbar(i,j,1)
# ifdef WEC
     &                                        + vst2d(i,j)
# endif
     &                                                   )
          enddo
        enddo

      endif  ! first time

      ds=1.D0/dble(nz)
      do j=j0,j1
        do i=i0,i1
          z_w(i,j,0)=-h(i,j)
        enddo

        do k=1,nz,+1   !--> irreversible because of recursion in Hz

          cff_w=hc*ds* dble(k-nz)
          cff_r=hc*ds*(dble(k-nz)-0.5)

          cff1_w=Cs_w(k)
          cff1_r=Cs_r(k)

          do i=i0,i1

            z_w(i,j,k)=zeta(i,j,knew) +(zeta(i,j,knew)+h(i,j))
     &                        *(cff_w+cff1_w*h(i,j))*hinv(i,j)

            z_r(i,j,k)=zeta(i,j,knew) +(zeta(i,j,knew)+h(i,j))
     &                        *(cff_r+cff1_r*h(i,j))*hinv(i,j)

            Hz(i,j,k)=z_w(i,j,k)-z_w(i,j,k-1)
          enddo
        enddo
      enddo
# ifdef EXCHANGE
      if (iic == 0) then
        call exchange_xxx(hinv)
      endif
      call exchange_xxx(z_w,z_r,Hz)
# endif

#if defined NHMG || defined NONTRAD_COR
!     Compute slopes
      do k = 1,nz
         do j = 0,ny+1
            do i = 0,nx+1
               dzdxi (i,j,k)=0.5*(z_r(i+1,j,k)-z_r(i-1,j,k))*pm(i,j)
               dzdeta(i,j,k)=0.5*(z_r(i,j+1,k)-z_r(i,j-1,k))*pn(i,j)
            enddo
         enddo
      enddo

      ! Physical boundaries: the boundary-ring (i=1/nx, j=1/ny) slopes
      ! are valid as computed -- the centered difference only reaches the
      ! ghost ring, whose z_r comes from real data (no-gradient zeta/h
      ! fill).  Only the ghost ring itself lacks a centered difference;
      ! mirror the boundary ring into it, consistent with the dz halo
      ! treatment in libnhmg (fill_outer_halos).  Zeroing the boundary
      ! ring here makes the NHMG operator inconsistent along the walls
      ! and drives spurious circulation, O(amp^2) in surface amplitude.
      if (inode.eq.0.and.(.not.ew_periodic)) then
         do k = 1,nz
            do j = 0,ny+1
               dzdxi (0,j,k) = dzdxi (1,j,k)
               dzdeta(0,j,k) = dzdeta(1,j,k)
            enddo
         enddo
       endif
      if (inode.eq.npx-1.and.(.not.ew_periodic)) then
         do k = 1,nz
            do j = 0,ny+1
               dzdxi (nx+1,j,k) = dzdxi (nx,j,k)
               dzdeta(nx+1,j,k) = dzdeta(nx,j,k)
            enddo
         enddo
       endif
      if (jnode.eq.0.and.(.not.ns_periodic)) then
         do k = 1,nz
            do i = 0,nx+1
               dzdxi (i,0,k) = dzdxi (i,1,k)
               dzdeta(i,0,k) = dzdeta(i,1,k)
            enddo
         enddo
       endif
      if (jnode.eq.npy-1.and.(.not.ns_periodic)) then
         do k = 1,nz
            do i = 0,nx+1
               dzdxi (i,ny+1,k) = dzdxi (i,ny,k)
               dzdeta(i,ny+1,k) = dzdeta(i,ny,k)
            enddo
         enddo
       endif

#endif /* NHMG or NONTRAD_COR */


      ! zstar prototype: record fast-filtered zeta history (set_depth
      ! runs at the end of each fast window and at init, so
      ! zeta(:,:,knew) here is the filtered free surface).
      if (.not.zstar_init) then
        allocate( zflt    (GLOBAL_2D_ARRAY) )
        allocate( zflt_bak(GLOBAL_2D_ARRAY) )
        allocate( zstar   (GLOBAL_2D_ARRAY) )
        zflt = zeta(:,:,knew)
        zstar_init = .true.
      endif
      zflt_bak = zflt
      zflt = zeta(:,:,knew)

      end subroutine set_depth !]
! ----------------------------------------------------------------------------
      subroutine set_zstar ![
      ! compute zeta* = 1.5*zflt - 0.5*zflt_bak (centered n+1/2)
      use dimensions
      use roms_mpi
      use zstar_geom
      implicit none
      zstar = 1.5*zflt - 0.5*zflt_bak
      call exchange_xxx(zstar)
      end subroutine set_zstar !]
! ----------------------------------------------------------------------------
      subroutine set_HUV ![

      use param
      use wec_frc
      use grid
      use ocean_vars
      use scalars
      use roms_mpi

      implicit none
      integer :: i,j,k


      do k=1,nz
        do j=j0,j1
          do i=1,i1
            FlxU(i,j,k)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)

     &                                       *( u(i,j,k,nrhs)
# ifdef WEC
     &                                            +ust(i,j,k)
# endif
     &                                                      )
! JM make sure that ust,vst are properly masked in wec_frc and remove
# if defined WEC && defined MASKING
            FlxU(i,j,k)=FlxU(i,j,k)*umask(i,j)
# endif
          enddo
        enddo
        do j=1,j1
          do i=i0,i1
            FlxV(i,j,k)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)
! Daniel add/edit for WEC
!     &                                         *v(i,j,k,nrhs)
     &                                       *( v(i,j,k,nrhs)
# ifdef WEC
     &                                           + vst(i,j,k)
# endif
     &                                                      )
# if defined WEC && defined MASKING
            FlxV(i,j,k)=FlxV(i,j,k)*vmask(i,j)
# endif
          enddo
        enddo
      enddo
# ifdef EXCHANGE
      call exchange_xxx(FlxU,FlxV)
# endif
      end subroutine set_HUV !]
! ----------------------------------------------------------------------------
      subroutine set_HUV1 ![

      use param
      use wec_frc
      use coupling
      use grid
      use ocean_vars
      use scalars
      use dimensions
      use diagnostics
      use roms_mpi

      implicit none
      integer :: i,j,k
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:nz) :: DC,FC
      real cff
      integer kbak

# if defined EXTRAP_BAR_FLUXES && defined KEEP_CORIOLIS
      kbak=knew-1
      if (kbak < 1) kbak=4
# endif

      do j=j0,j1      ! big j-loop

        do i=1,i1
          DC(i,nz)=0.5*(Hz(i,j,nz)+Hz(i-1,j,nz))*dn_u(i,j)
          DC(i,0)=DC(i,nz)
          FC(i,0)=DC(i,nz)*u(i,j,nz,nnew)
        enddo
        do k=nz-1,1,-1
          do i=1,i1
            DC(i,k)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)
            DC(i,0)=DC(i,0)+DC(i,k)
            FC(i,0)=FC(i,0)+DC(i,k)*u(i,j,k,nnew)
          enddo
        enddo


# ifdef CORR_COUPLED_MODE
#  ifdef EXTRAP_BAR_FLUXES
        if (FIRST_TIME_STEP) then
#  endif
          do i=1,i1
            FC(i,0)=(FC(i,0)-DU_avg1(i,j))/DC(i,0)
#  ifdef WEC
     &                                 +ust2d(i,j)
#  endif
          enddo
          !! JM: Very first time step only

#  ifdef EXTRAP_BAR_FLUXES
! PAC23 verified setting with theoretical alpha_max=1.0877.
! This should match the simlar-labelled setting in step2d_FB.F.
! All the other settings below  are to be considered as
! experimental.

#   define NOW 3.63D0
#   define MID 4.47D0
#   define BAK 2.05D0


!  DC(i,0) is the barotropic surface H*dy = dy*integral(Hz)
!  FC(i,0) is the barotropic flux dy*integral(Hz*u)

        else ! not first time step
          do i=1,i1
#   ifdef KEEP_CORIOLIS
            cff=1./DC(i,0)
            ubar(i,j,kbak)=cff*( NOW*DU_avg1(i,j)
     &                          -MID*DU_avg2(i,j)
     &                       +BAK*DU_avg_bak(i,j)
     &                                          )
            FC(i,0)=cff*FC(i,0)-ubar(i,j,kbak)  ! now FC(i,0) is mismatch in[m/s]
#   else
            FC(i,0)=( FC(i,0) -NOW*DU_avg1(i,j)
     &                        +MID*DU_avg2(i,j)
     &                     -BAK*DU_avg_bak(i,j)
     &               )/DC(i,0)
#   endif
          enddo
        endif
#  endif
     !! now FC(i,0) => mismatch between  int(U*dz)/int(dz) and ubar
# else
        do i=1,i1
          FC(i,0)=(FC(i,0)-DU_avg2(i,j))/DC(i,0) ! FC(i,0) is mismatch in [m/s]
#  ifdef WEC
     &                               +ust2d(i,j)
#  endif
        enddo
# endif

        do k=1,nz
          do i=1,i1
!!JM make sure that ust, vst are masked in wec_frc and fix this!!
            u(i,j,k,nnew)=(u(i,j,k,nnew)-FC(i,0)) !! removing ubar mismatch from u
# ifdef MASKING
     &                                *umask(i,j)
#  ifdef WEC
     &               +ust(i,j,k)*(umask(i,j)-1.0)
#  endif
# endif
            FlxU(i,j,k)=DC(i,k)*( u(i,j,k,nnew)
# ifdef WEC
     &                                +ust(i,j,k)
# endif
     &                                          )
          enddo
        enddo

        if (j >= 1) then
          do i=i0,i1
            DC(i,nz)=0.5*(Hz(i,j,nz)+Hz(i,j-1,nz))*dm_v(i,j)
            DC(i,0)=DC(i,nz)
            FC(i,0)=DC(i,nz)*v(i,j,nz,nnew)
          enddo
          do k=nz-1,1,-1
            do i=i0,i1
              DC(i,k)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)
              DC(i,0)=DC(i,0)+DC(i,k)
              FC(i,0)=FC(i,0)+DC(i,k)*v(i,j,k,nnew)
            enddo
          enddo
# ifdef CORR_COUPLED_MODE
#  ifdef EXTRAP_BAR_FLUXES
          if (FIRST_TIME_STEP) then
#  endif
            do i=i0,i1
              FC(i,0)=(FC(i,0)-DV_avg1(i,j))/DC(i,0)
#  ifdef WEC
     &                                   +vst2d(i,j)
#  endif
            enddo

#  ifdef EXTRAP_BAR_FLUXES
          else
            do i=i0,i1
#   ifdef KEEP_CORIOLIS
              cff=1./DC(i,0)
              vbar(i,j,kbak)=cff*( NOW*DV_avg1(i,j)
     &                            -MID*DV_avg2(i,j)
     &                         +BAK*DV_avg_bak(i,j)
     &                                            )
              FC(i,0)=cff*FC(i,0)-vbar(i,j,kbak)
#   else
              FC(i,0)=( FC(i,0) -NOW*DV_avg1(i,j)
     &                          +MID*DV_avg2(i,j)
     &                       -BAK*DV_avg_bak(i,j)
     &                 )/DC(i,0)
! Daniel add --> not sure if WEC should be in here. DevinD also not sure but suspect so.
#ifdef WEC
     &                   +vst2d(i,j)
#endif
#   endif
            enddo
          endif
#  endif
# else
          do i=i0,i1
            FC(i,0)=(FC(i,0)-DV_avg2(i,j))/DC(i,0)
#  ifdef WEC
     &                                 +vst2d(i,j)
#  endif
          enddo
# endif
          do k=1,nz
            do i=i0,i1
!! JM same thing as above for ust. (masking in wec_frc)
              v(i,j,k,nnew)=(v(i,j,k,nnew)-FC(i,0))
# ifdef MASKING
     &                                  *vmask(i,j)
#  ifdef WEC
     &                 +vst(i,j,k)*(vmask(i,j)-1.0)
#  endif
# endif
              FlxV(i,j,k)=DC(i,k)*( v(i,j,k,nnew)
# ifdef WEC
     &                                +vst(i,j,k)
# endif
     &                                          )

            enddo
          enddo
        endif
      enddo  ! <-- j-loop

# ifdef EXCHANGE
#  if defined EXTRAP_BAR_FLUXES && defined KEEP_CORIOLIS
      call exchange_xxx(FlxU,u(:,:,:,nnew),ubar)
      call exchange_xxx(FlxV,v(:,:,:,nnew),vbar)
#  else
      call exchange_xxx(FlxU,FlxV,u(:,:,:,nnew),v(:,:,:,nnew))
#  endif
# endif
      end subroutine set_HUV1 !]

