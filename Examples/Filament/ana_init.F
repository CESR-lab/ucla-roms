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

      real :: Ampl, fac, x,y, x0,y0, cff, cff0,cff1,cff2,cff3

      ! Devin added:
!      real, dimension(GLOBAL_1DX_ARRAY)   :: bs, h2d                 ! h -> h2d
      real               :: imp, ksi, FF, cff_avg   ! F -> FF
      integer            :: ncid, varid, ierr
      logical, parameter :: read_inputs = .false.
      real, dimension(GLOBAL_2D_ARRAY)    :: h_sbl
      real, dimension(GLOBAL_2D_ARRAY, N) :: dbdx_xr
      real, allocatable, dimension(:,:,:) :: b, a1, a2, a3  ! db/dx at xi_rho, s_rho. v sits at xi_rho, and constant in y


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
      allocate( a3(GLOBAL_2D_ARRAY, N) )

      if (read_inputs) then  ! read in data from input file
        ! Read in initial conditions from python TTW script:
        ierr=nf90_open('python/input_ideal.nc', nf90_nowrite, ncid)
        call nc_read_var(ncid, ubar(:,:,1), 1, 'ubar', up_var, 0, ierr)
        if(ierr/=0) call handle_ierr(ierr,'ana_init:: error reading!')
        call nc_read_var(ncid, vbar(:,:,1), 1, 'vbar', vp_var, 0, ierr)
        if(ierr/=0) call handle_ierr(ierr,'ana_init:: error reading!')
        call nc_read_var(ncid, u(:,:,:,1),  N, 'u',    up_var, 0, ierr)
        if(ierr/=0) call handle_ierr(ierr,'ana_init:: error reading!')
        call nc_read_var(ncid, v(:,:,:,1),  N, 'v',    vp_var, 0, ierr)
        if(ierr/=0) call handle_ierr(ierr,'ana_init:: error reading!')
      endif


      ! Can also just read this in...
      b=0.0
      do j=jstrR,jendR         ! j is length of filament
        do i=istrR-1,iendR+1   ! go beyond buffer to allow for averaging later. Since analytical, doesn't matter.

!          bs(i)  = bs0 - (db0 * exp(-(( xr(i,j)/1000 /L)**2)))  ! DH uses KM
!          h2d(i) = h0  + (dh0 * exp(-(( xr(i,j)/1000 /L)**2)))
          h_sbl(i,j) = h0  + (dh0 * exp(-(( xr(i,j) /L )**2)))
!          h_sbl(i,j) = h0  - (dh0 * exp(-(( xr(i,j) /L )**2)))



          ! finding a
!          imp=N02*(h3-h2d(i))/(bs(i)-b3)
!          aa_cff=x1
!          do while ( mu*aa_cff / ( exp( aa_cff )-1.0 ) > imp )
!            aa_cff = aa_cff + 0.01 ! precision (maybe use 10 (1000x 0.01)
!          enddo
          !if(mynode==0 .and. j==1) print *, 'Final aa_cff=', aa_cff

          do k=1,N

            a1(i,j,k) = Nb * (z_r(i,j,k) + HD)
            a2(i,j,k) = (0.5 * N0) * ( (1+ B_cff) * z_r(i,j,k) - ( 1- B_cff)
     &          *( h_sbl(i,j) + lambda_inv * log(cosh((1./lambda_inv) *(z_r(i,j,k) + h_sbl(i,j))))))
!            a3(i,j,k) = (0.5 * N2) * ( (1+ B_cff) * z_r(i,j,k) - ( 1- B_cff)
!     &          *( h_sbl(i,j) + dh2 + lambda_inv2 * log(cosh((1./lambda_inv2) *(z_r(i,j,k) + h_sbl(i,j) + dh2)))))
            b(i,j,k) = bs0 + a1(i,j,k) + a2(i,j,k) ! - a3(i,j,k)

            ! buoyancy do as temperature:
            t(i,j,k,1,itemp) = b(i,j,k) / (2.0e-4 * 9.81)

!            if     ( z_r(i,j,k) > -h2d(i) ) then                   ! first layer by surface
!              t(i,j,k,1,itemp)=bs(i)*1000.0
!            elseif ( z_r(i,j,k) < -h3     ) then                   ! background / bottom depth
!              t(i,j,k,1,itemp)=N02*( HD + z_r(i,j,k) ) *1000.0
!            else                                                   ! middle layer
!              ksi=(-z_r(i,j,k)-h2d(i)) / (h3-h2d(i))
!              FF=( exp(aa_cff*(1.-(ksi**mu)))-1.) / ( exp(aa_cff)-1.)
!              t(i,j,k,1,itemp)= ( ((bs(i)-b3)*FF) + b3 ) *1000.0
!              !if(j==1 .and. k==100) print *, 'FF=', FF
!            endif

          enddo
        enddo
      enddo

!      error stop 'STOP HERE!!!'

      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
!            t(i,j,k,1,itemp)=20.0
!            t(i,j,k,1,itemp)=4.+10.*exp(z_r(i,j,k)/50.)
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
      cff_avg = 1.0/N
      do k=1,N
        do j=jstrR,jendR              ! pn * 0.5 = / 2 * dx
          do i=istrR,iendR            ! loop ranges too long because of averaging but leave for now
            dbdx_xr(i,j,k) = ( t(i+1,j,k,1,itemp) - t(i-1,j,k,1,itemp) ) * pm(i,j) *0.5 * 0.001 ! divide by 1000 for scale
            if(k==1) then
              v(i,j,k,1)= dbdx_xr(i,j,k) * Hz(i,j,k) / f(i,j)
            else
              v(i,j,k,1)= v(i,j,k-1,1) + dbdx_xr(i,j,k) * Hz(i,j,k) / f(i,j)
            endif
            vbar(i,j,1)=vbar(i,j,1)+v(i,j,k,1) * cff_avg
          enddo
        enddo
      enddo

      zeta(:,:,1)=0.
      do j=jstrR,jendR
        do i=istrR,iendR
          zeta(i,j,1)=zeta(i-1,j,1) + f(i,j) * v(i,j,N,1) * dm_r(i,j) / g
!          zeta(i,j,1)=zeta(i-1,j,1) + f(i,j) * vbar(i,j,1) * dm_r(i,j) / g
          ubar(i,j,1)=0.           ! non-trivial initialization required.
          !vbar(i,j,1)=0.          ! Note: A code to initialize T [and S]
        enddo                      ! must always be supplied for 3D
      enddo                        ! applications.

      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            u(i,j,k,1)=0.
!            v(i,j,k,1)=0.
            u(i,j,k,2)=u(i,j,k,1)  ! copy from reading above
            v(i,j,k,2)=v(i,j,k,1)
          enddo
        enddo
      enddo


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
