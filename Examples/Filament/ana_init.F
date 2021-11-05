#include "cppdefs.opt"
#ifdef ANA_INITIAL
                                  ! Set initial conditions for momentum
      subroutine ana_init (tile)  ! free surface, and tracer variables

      use param
      use hidden_mpi_vars
      use private_scratch

      implicit none               ! using analytical expressions.
      integer tile

# include "compute_tile_bounds.h"
      call ana_init_tile( istr,iend,jstr,jend,
# ifdef SOLVE3D
     &                                 A3d(1,1)
# else
     &                                 A2d(1,1)
# endif
     &                                        )
      end

      subroutine ana_init_tile (istr,iend,jstr,jend, AA)

      use param
      use tracers, only: t
      use boundary
      use climat
      use eos_vars
      use grid
      use mixing
      use ocean2d
      use ocean3d
      use scalars
      use read_write
      use netcdf

      implicit none

      integer :: istr,iend,jstr,jend, i,j,k, itrc, i0,j0
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,0:N) :: AA

      real :: Ampl, fac, x,y, x0,y0, cff0,cff1,cff2,cff3

      ! Devin added:
      real :: cff_avg
      real, dimension(GLOBAL_2D_ARRAY)    :: h_sbl
      real, dimension(GLOBAL_2D_ARRAY, N) :: dbdx_xr    ! db/dx at xi_rho, s_rho. v sits at xi_rho, and constant in y
      real, allocatable, dimension(:,:,:) :: b, a1, a2


      ! Daniel values:
      ! - SMC Density Initial Condition:
      real :: b0          = 6.4e-3
      real :: B_cff       = 0.025    ! gamma in JCM
      real :: lambda_inv  = 10.0      ! m. scale of transition from surface to pycnocline
      real :: Nb          = 1.0e-7   ! surface s^-2
      real :: N0          = 3.4e-5   ! background
      real :: h0          = 60.0     ! m. boundary layer beyond filament
      real :: dh0         = 15.      ! m
      real :: dh2         = 10.0     ! m. Buffer between pycnocline and deeper stratification.
      real :: N2          = 1.5e-4   ! Daniel said try 0
      real :: h2          = 70.0     !
      real :: L           = 2000.0     ! m
      real :: bs0         = 8.51e-3  ! ms-2

# include "compute_auxiliary_bounds.h"

      allocate(  b(GLOBAL_2D_ARRAY, N) )
      allocate( a1(GLOBAL_2D_ARRAY, N) )
      allocate( a2(GLOBAL_2D_ARRAY, N) )

      ! Can also just read this in...
      b=0.0
      do j=jstrR-1,jendR+1     ! j is length of filament. Go beyond buffer for averaging to v.
        do i=istrR-1,iendR+1   ! go beyond buffer to allow for averaging later. Since analytical, doesn't matter.

          h_sbl(i,j) = h0  + (dh0 * exp(-(( xr(i,j) /L )**2)))

          do k=1,N

            a1(i,j,k) = Nb * (z_r(i,j,k) + HD)
            a2(i,j,k) = (0.5 * N0) * ( (1+ B_cff) * z_r(i,j,k) - ( 1- B_cff)
     &          *( h_sbl(i,j) + lambda_inv * log(cosh((1./lambda_inv) *(z_r(i,j,k) + h_sbl(i,j))))))

            b(i,j,k) = bs0 + a1(i,j,k) + a2(i,j,k)

            ! buoyancy do as temperature:
            t(i,j,k,1,itemp) = b(i,j,k) / (2.0e-4 * 9.81)

          enddo
        enddo
      enddo


      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
!            t(i,j,k,1,itemp)=20.0
            t(i,j,k,1,isalt)=36.
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
          enddo
        enddo
      enddo

      ! Set geostrophic flow in u-direction based on temperature as proxy for buoyancy:
      dbdx_xr=0.0
      v(:,:,:,1)=0.0
      vbar(:,:,1)=0.0
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            dbdx_xr(i,j,k) = ( b(i+1,j,k) - b(i-1,j,k) ) * pm(i,j) * 0.5      ! pm * 0.5 = / 2 * dx
            if(k==1) then
              v(i,j,k,1)= dbdx_xr(i,j,k) * Hz(i,j,k) / f(i,j)
            else
              v(i,j,k,1)= v(i,j,k-1,1) + dbdx_xr(i,j,k) * Hz(i,j,k) / f(i,j)
            endif
            vbar(i,j,1)=vbar(i,j,1)+v(i,j,k,1) * Hz(i,j,k)! remove as wrong
          enddo
        enddo
      enddo
      vbar(:,:,1) = vbar(:,:,1) / HD  ! average vertical integral over depth
                                      ! excludes zeta since 'v' was calculated without it, but should be there

      ! zeta needed for vbar calculation:
      zeta(:,:,1)=0.
      do j=jstrR,jendR
        do i=istrR,iendR                                              ! N.B.!! actually need v at surface!

          !zeta(i,j,1)=zeta(i-1,j,1) + f(i,j) * v(i,j,N,1) * dm_r(i,j) / g

          ! rather use trapezoidal rule for integration.
          ! but since v and z both at xi_rho need to use 0.25 (1+2+1)
          zeta(i,j,1)=zeta(i-1,j,1) + 0.25 * f(i,j) * dm_r(i,j) / g   ! so should be extrapolating
     &               * ( v(i-1,j,N,1) + 2*v(i,j,N,1) + v(i+1,j,N,1) ) ! missing i-1 beyond buffer
                                                                      ! only e-6 difference in zeta with trapezoidal rule.
          ubar(i,j,1)=0.
          !vbar(i,j,2)=vbar(i,j,1)  ! DevinD experiment, made no difference
        enddo
      enddo

      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            u(i,j,k,1)=0.
!            v(i,j,k,1)=0.
            u(i,j,k,2)=u(i,j,k,1)  ! copy from reading above
            v(i,j,k,2)=v(i,j,k,1)
!            zeta(i,j,2)=zeta(i,j,1)  ! DevinD added, made no difference
          enddo
        enddo
      enddo

      ! DevinD: need to exchange any extra variables I might configure here!

! The rest of the code is independent of the specific configuration.

! Analytical problems with open boundaries may require external forcing
! data to be applied at inflow side.  For the problems of with time
! independent inflow conditions it is sufficient to save initial
! condition values at the boundary point into boundary forcing arrays,
! which will remain constant thereafter and be used as the source of
! forcing data.  Note that this automatically guarantees consistency
! between the initial state and side boundary conditions.


# if defined OBC_WEST  || defined OBC_EAST ||\
     defined OBC_SOUTH || defined OBC_NORTH
#  ifdef ANA_BRY
#   ifdef OBC_WEST
      if (WESTERN_EDGE) then
#    ifdef M2_FRC_BRY
        do j=jstrR,jendR                     ! Save perimeter values
          zeta_west(j)=zeta(istr-1,j,1)      ! of the initial fields
          ubar_west(j)=ubar(istr  ,j,1)      ! to be used as external
          vbar_west(j)=vbar(istr-1,j,1)      ! forcing data for side
        enddo                                ! boundaries.
#    endif
#    ifdef SOLVE3D
        do k=1,N
          do j=jstrR,jendR
#     ifdef M3_FRC_BRY
            u_west(j,k)=u(istr  ,j,k,1)
            v_west(j,k)=v(istr-1,j,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_west(j,k,itemp)=t(istr-1,j,k,1,itemp)
#     endif
          enddo
        enddo
#    endif
      endif
#   endif

#   ifdef OBC_EAST
      if (EASTERN_EDGE) then
#    ifdef M2_FRC_BRY
        do j=jstrR,jendR
          zeta_east(j)=zeta(iend+1,j,1)
          ubar_east(j)=ubar(iend+1,j,1)
          vbar_east(j)=vbar(iend+1,j,1)
        enddo
#    endif
#    ifdef SOLVE3D
        do k=1,N
          do j=jstrR,jendR
#     ifdef M3_FRC_BRY
            u_east(j,k)=u(iend+1,j,k,1)
            v_east(j,k)=v(iend+1,j,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_east(j,k,itemp)=t(iend+1,j,k,1,itemp)
#     endif
          enddo
        enddo
#    endif
      endif
#   endif


#   ifdef OBC_SOUTH
      if (SOUTHERN_EDGE) then
#    ifdef M2_FRC_BRY
        do i=istrR,iendR
          zeta_south(i)=zeta(i,jstr-1,1)
          ubar_south(i)=ubar(i,jstr-1,1)
          vbar_south(i)=vbar(i,jstr  ,1)
        enddo
#    endif
#    ifdef SOLVE3D
        do k=1,N
          do i=istrR,iendR
#     ifdef M3_FRC_BRY
            u_south(i,k)=u(i,jstr-1,k,1)
            v_south(i,k)=v(i,jstr  ,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_south(i,k,itemp)=t(i,jstr-1,k,1,itemp)
#     endif
          enddo
        enddo
#    endif
      endif
#   endif

#   ifdef OBC_NORTH
      if (NORTHERN_EDGE) then
#    ifdef M2_FRC_BRY
        do i=istrR,iendR
          zeta_north(i)=zeta(i,jend+1,1)
          ubar_north(i)=ubar(i,jend+1,1)
          vbar_north(i)=vbar(i,jend+1,1)
        enddo
#    endif
#    ifdef SOLVE3D
        do k=1,N
          do i=istrR,iendR
#     ifdef M3_FRC_BRY
            u_north(i,k)=u(i,jend+1,k,1)
            v_north(i,k)=v(i,jend+1,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_north(i,k,itemp)=t(i,jend+1,k,1,itemp)
#     endif
          enddo
        enddo
#    endif
      endif
#   endif
#  else /* alternative to ANA_BRY */
#   ifdef UCLIMATOLOGY
      do j=jstrR,jendR                   ! Save initial data into
        do i=istrR,iendR                 ! climatology for the
          ubclm(i,j)=ubar(i,j,1)         ! subsequent use as inflow
          vbclm(i,j)=vbar(i,j,1)         ! boundary conditions
        enddo
      enddo
#   endif
#   ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
#    ifdef UCLIMATOLOGY
            uclm(i,j,k)=u(i,j,k,1)
            vclm(i,j,k)=v(i,j,k,1)
#    endif
#    ifdef TCLIMATOLOGY
            tclm(i,j,k,itemp)=t(i,j,k,1,itemp)
#    endif
          enddo
        enddo
      enddo
#   endif /* SOLVE3D */
#  endif /* ANA_BRY */
# endif /* at least one of OBC_XXXXX */


# ifdef EXCHANGE
      call exchange2d_3_tile(istr,iend,jstr,jend,
     &                       zeta(START_2D_ARRAY,1),
     &                       ubar(START_2D_ARRAY,1),
     &                       vbar(START_2D_ARRAY,1))
#  ifdef SOLVE3D
      call exchange_2_tile(istr,iend,jstr,jend,
     &                     u(START_2D_ARRAY,1,1), N,
     &                     v(START_2D_ARRAY,1,1), N)
      do itrc=1,NT
        call exchange_tile(istr,iend,jstr,jend,
     &                t(START_2D_ARRAY,1,1,itrc), N)
      enddo
#  endif
#  ifdef ISWAKE
#   ifdef UCLIMATOLOGY
      call exchange2d_2_tile(istr,iend,jstr,jend, ubclm,vbclm)
#   endif
#   ifdef SOLVE3D
#    ifdef UCLIMATOLOGY
      call exchange_2_tile(istr,iend,jstr,jend, uclm,N, vclm,N)
#    endif
#    ifdef TCLIMATOLOGY
      call exchange_tile(istr,iend,jstr,jend,
     &                    tclm(START_2D_ARRAY,1,itemp), N)
#    endif
#   endif
#  elif defined COLD_FILAMENT
      call exchange_tile(istr,iend,jstr,jend, Akv,N+1)

#  endif
# endif  /* EXCHANGE */


# if defined OBC_WEST && (defined NS_PERIODIC || defined MPI)
      if (WESTERN_EDGE) then
#  ifdef M2_FRC_BRY
        call exch_bry_EW_tile (jstr,jend, zeta_west,1)
        call exch_bry_EW_tile (jstr,jend, ubar_west,1)
        call exch_bry_EW_tile (jstr,jend, vbar_west,1)
#  endif
#  ifdef M3_FRC_BRY
        call exch_bry_EW_tile (jstr,jend, u_west,N)
        call exch_bry_EW_tile (jstr,jend, v_west,N)
#  endif
#  ifdef T_FRC_BRY
        do itrc=1,NT
          call exch_bry_EW_tile (jstr,jend, t_west(0,1,itrc),N)
        enddo
#  endif
      endif
# endif

# if defined OBC_EAST && (defined NS_PERIODIC || defined MPI)
      if (EASTERN_EDGE) then
#  ifdef M2_FRC_BRY
        call exch_bry_EW_tile (jstr,jend, zeta_east,1)
        call exch_bry_EW_tile (jstr,jend, ubar_east,1)
        call exch_bry_EW_tile (jstr,jend, vbar_east,1)
#  endif
#  ifdef M3_FRC_BRY
        call exch_bry_EW_tile (jstr,jend, u_east,N)
        call exch_bry_EW_tile (jstr,jend, v_east,N)
#  endif
#  ifdef T_FRC_BRY
        do itrc=1,NT
          call exch_bry_EW_tile (jstr,jend, t_east(0,1,itrc),N)
        enddo
#  endif
      endif
# endif

# if defined OBC_SOUTH && (defined EW_PERIODIC || defined MPI)
      if (SOUTHERN_EDGE) then
#  ifdef M2_FRC_BRY
        call exch_bry_NS_tile (istr,iend, zeta_south,1)
        call exch_bry_NS_tile (istr,iend, ubar_south,1)
        call exch_bry_NS_tile (istr,iend, vbar_south,1)
#  endif
#  ifdef M3_FRC_BRY
        call exch_bry_NS_tile (istr,iend, u_south,N)
        call exch_bry_NS_tile (istr,iend, v_south,N)
#  endif
#  ifdef T_FRC_BRY
        do itrc=1,NT
          call exch_bry_NS_tile (istr,iend, t_south(0,1,itrc),N)
        enddo
#  endif
      endif
# endif

# if defined OBC_NORTH && (defined EW_PERIODIC || defined MPI)
      if (NORTHERN_EDGE) then
#  ifdef M2_FRC_BRY
        call exch_bry_NS_tile (istr,iend, zeta_north,1)
        call exch_bry_NS_tile (istr,iend, ubar_north,1)
        call exch_bry_NS_tile (istr,iend, vbar_north,1)
#  endif
#  ifdef M3_FRC_BRY
        call exch_bry_NS_tile (istr,iend, u_north,N)
        call exch_bry_NS_tile (istr,iend, v_north,N)
#  endif
#  ifdef T_FRC_BRY
        do itrc=1,NT
          call exch_bry_NS_tile (istr,iend, t_north(0,1,itrc),N)
        enddo
#  endif
      endif
# endif
      end
#else
      subroutine ana_init_empty
      end
#endif /* ANA_INITIAL */
