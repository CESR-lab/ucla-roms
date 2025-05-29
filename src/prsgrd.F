#include "cppdefs.opt"

      subroutine prsgrd
      use param
      use hidden_mpi_vars
      use private_scratch
      use dimensions
      use tides
      use diagnostics
      use coupling
      use ocean_vars
      use eos_vars
      use grid
      use scalars
#ifdef WEC
      use wec_frc, cff1 => cff_ext, cff2 => cff_ext2
#endif

      implicit none

      ! local
      integer         :: i,j,k
      real            :: grho, HalfGRho, cff, cfr
      real, parameter :: OneFifth=0.2, OneTwelfth=1./12., epsil=0.
      real,dimension(:,:),pointer   :: dR,dZ,
     &                                 FC,dZx,rx,dRx
      real,dimension(:,:,:),pointer :: ru,rv,P
# ifdef SPLIT_EOS
      real,dimension(:,:,:),pointer :: rho
      real dpth
# endif
# ifdef NHMG
      real,dimension(:,:,:),pointer :: rw
# endif
# ifdef ADV_ISONEUTRAL
      real :: r0g
# endif

      ! set pointers/scratch space
      dR(1-bf:nx+bf,0:N)         => scratch_2d(:,1)
      dZ(1-bf:nx+bf,0:N)         => scratch_2d(:,2)
      FC(1-bf:nx+bf,1-bf:ny+bf)  => scratch_2d(:,3)
      rx(1-bf:nx+bf,1-bf:ny+bf)  => scratch_2d(:,4)
      dRx(1-bf:nx+bf,1-bf:ny+bf) => scratch_2d(:,5)
      ru(1-bf:nx+bf,1-bf:ny+bf,0:nz)  => scratch_3d(:,1)
      rv(1-bf:nx+bf,1-bf:ny+bf,0:nz)  => scratch_3d(:,2)
      P(1-bf:nx+bf,1-bf:ny+bf,1:nz)   => scratch_3d(:,3)
# ifdef SPLIT_EOS
      rho(1-bf:nx+bf,1-bf:ny+bf,1:nz) => scratch_3d(:,4)
# endif
# ifdef NHMG
      rw(1-bf:nx+bf,1-bf:ny+bf,0:nz)  => scratch_3d(:,5)
# endif

! A non-conservative Density-Jacobian scheme using cubic polynomial
! fits for rho and z_r as functions of nondimensional coordinates xi,
! eta, and s (basically their respective fortran indices). The cubic
! polynomials are constructed by specifying first derivatives of
! interpolated fields on co-located (non-staggered) grid locations.
! The derivatives are computed using harmonic (rather that algebraic)
! averaging of elementary differences, which guarantees monotonicity
! of the resultant interpolant.

! In the code below, if CPP-switch SPLIT_EOS is defined, the Equation
! of State (EOS) is assumed to have form
!
!       rho(T,S,z) = rho1(T,S) + qp1(T,S)*dpth*[1.-qp2*dpth]
!
! where "rho1" is density at computed 1 atm pressure; "qp1" is
! compressibility coefficient, which does not depend on pressure;
! "qp2" is just a constant; and  dpth=-z is absolute depth playing
! the role of Boussinesq proxy for EOS pressure [hence "dpth" is
! positive, increasing downward. It should be noted that dpth=-z,
! counting from the unperturbed state of free surface, instead of
! dpth=(zeta-z) to compy with Boussinesq energetic consistencty,
! c.f., Young, 2010]. In this case
!
!   d rho    d rho1   d qp1                                    d z
!  ------- = ------ + ----- *dpth*[..] - qp1*[1.-2.*qp2*dpth]*------
!   d s,x     d s,x   d s,x                                    d s,x
!
!           |<--- adiabatic part --->|  |<--- compressible part --->|
!
! where the first two terms constitute "adiabatic derivative" of
! density, which is subject to harmonic averaging, while the last
! term is added in later. This approach guarantees that density
! profile reconstructed by cubic polynomial maintains its positive
! stratification in physical sense as long as discrete values of
! density are positively stratified.
!
! The scheme retains exact antisymmetry J(rho,z_r)=-J(z_r,rho) [with
! the exception of harmonic averaging algorithm in the case when CPP-
! switch SPLIT_EOS is defined, see above]. Setting parameter OneFifth
! (see above) to zero makes scheme become identical to the standard
! POM density Jacobian.

! Preliminary step (same for XI- and ETA-components):
!------------ ---- ----- --- --- --- ----------------

      grho=g/rho0
      HalfGRho=0.5*grho

      do j=0,ny
        do k=1,nz-1
          do i=0,nx
            dZ(i,k)=z_r(i,j,k+1)-z_r(i,j,k)
# ifdef SPLIT_EOS
c**         dpth=z_w(i,j,N)-0.5*(z_r(i,j,k+1)+z_r(i,j,k))
            dpth=          -0.5*(z_r(i,j,k+1)+z_r(i,j,k))

            dR(i,k)=rho1(i,j,k+1)-rho1(i,j,k)            ! elementary
     &              +(qp1(i,j,k+1)-qp1(i,j,k))           ! adiabatic
     &                     *dpth*(1.-qp2*dpth)           ! difference
# else
            dR(i,k)=rho(i,j,k+1)-rho(i,j,k)
# endif
          enddo
        enddo
        ! Extend derivatives top and bottom
        ! dR and dZ are at 'w-points'
        do i=0,nx
          dR(i,nz)=dR(i,nz-1)
          dR(i,0 )=dR(i,1)
          dZ(i,nz)=dZ(i,nz-1)
          dZ(i,0 )=dZ(i,1)
        enddo
        ! Irreversible due to in-place overwriting of dZ and dR
        do k=nz,1,-1 
          do i=0,nx
            cff=2.*dZ(i,k)*dZ(i,k-1)
            dZ(i,k)=cff/(dZ(i,k)+dZ(i,k-1))

            cfr=2.*dR(i,k)*dR(i,k-1)
            if (cfr > epsil) then
              dR(i,k)=cfr/(dR(i,k)+dR(i,k-1))
            else
              dR(i,k)=0.
            endif
# ifdef SPLIT_EOS
c**         dpth=z_w(i,j,N)-z_r(i,j,k)
            dpth=          -z_r(i,j,k)
            dR(i,k)=dR(i,k)  -qp1(i,j,k)*dZ(i,k)*(1.-2.*qp2*dpth)
            rho(i,j,k)=rho1(i,j,k) +qp1(i,j,k)*dpth*(1.-qp2*dpth)
# endif
          enddo
        enddo

        do i=0,nx
          P(i,j,N)=g*z_w(i,j,N) + grho*( rho(i,j,N)
     &       +0.5*(rho(i,j,N)-rho(i,j,N-1))*(z_w(i,j,N)-z_r(i,j,N))
     &          /(z_r(i,j,N)-z_r(i,j,N-1)) )*(z_w(i,j,N)-z_r(i,j,N))
          if (pot_tides) then
            P(i,j,N) = P(i,j,N) - g*ptide(i,j)
          endif
        enddo
        do k=nz-1,1,-1
          do i=0,nx
            P(i,j,k)=P(i,j,k+1)+HalfGRho*( (rho(i,j,k+1)+rho(i,j,k))
     &                                     *(z_r(i,j,k+1)-z_r(i,j,k))

     &     -OneFifth*( (dR(i,k+1)-dR(i,k))*( z_r(i,j,k+1)-z_r(i,j,k)
     &                              -OneTwelfth*(dZ(i,k+1)+dZ(i,k)) )

     &                -(dZ(i,k+1)-dZ(i,k))*( rho(i,j,k+1)-rho(i,j,k)
     &                              -OneTwelfth*(dR(i,k+1)+dR(i,k)) )
     &                                                             ))
          enddo
        enddo
      enddo   !<-- j

      nullify(dZ)
      dZx(1-bf:nx+bf,1-bf:ny+bf) => scratch_2d(:,2)

#ifdef DIAGNOSTICS
      if (CORR_STAGE) then
        if (allocated(rho_d)) then    
          rho_d = rho(1:nx,1:ny,:)
        endif
        if (allocated(p_int)) then    
          p_int = 0
          do k=1,nz
           p_int = p_int + p(1:nx,1:ny,k)*Hz(1:nx,1:ny,k)
          enddo
        endif
      endif
#endif

! Compute XI-component of pressure gradient term:
!-------- ------------ -- -------- -------- -----

      do k=N,1,-1
        do j=1,ny
          do i=0,nx+1
            FC(i,j)=(z_r(i,j,k)-z_r(i-1,j,k))
# ifdef MASKING
     &                              *umask(i,j)
# endif
# ifdef SPLIT_EOS
c**         dpth=0.5*( z_w(i,j,N)+z_w(i-1,j,N)
c**   &               -z_r(i,j,k)-z_r(i-1,j,k))

            dpth=-0.5*(z_r(i,j,k)+z_r(i-1,j,k))

            rx(i,j)=( rho1(i,j,k)-rho1(i-1,j,k)          ! elementary
     &                +(qp1(i,j,k)-qp1(i-1,j,k))         ! adiabatic
     &                     *dpth*(1.-qp2*dpth) )         ! difference
# else
            rx(i,j)=(rho(i,j,k)-rho(i-1,j,k))
# endif
# ifdef MASKING
     &                              *umask(i,j)
# endif
          enddo
        enddo

        if (inode==0.and.(.not.ew_periodic)) then         ! Extrapolate elementary
          do j=1,ny               ! differences near physical
            FC(1,j)=FC(2,j)            ! boundaries to compensate.
            rx(1,j)=rx(2,j)            ! for reduced loop ranges.
          enddo
        endif
        if (inode==np_xi-1) then
          do j=1,ny
            FC(nx+1,j)=FC(nx,j)
            rx(nx+1,j)=rx(nx,j)
          enddo
        endif

        do j=1,ny
          do i=0,nx
            cff=2.*FC(i,j)*FC(i+1,j)
            if (cff > epsil) then
              dZx(i,j)=cff/(FC(i,j)+FC(i+1,j))
            else
              dZx(i,j)=0.
            endif

            cfr=2.*rx(i,j)*rx(i+1,j)
            if (cfr > epsil) then
              dRx(i,j)=cfr/(rx(i,j)+rx(i+1,j))
            else
              dRx(i,j)=0.
            endif
# ifdef SPLIT_EOS
            dRx(i,j)=dRx(i,j) -qp1(i,j,k)*dZx(i,j)
     &                      *(1.+2.*qp2*z_r(i,j,k))
# endif
          enddo               !--> discard FC, rx

          do i=1,nx
            ru(i,j,k)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*(
     &                              P(i-1,j,k)-P(i,j,k)-HalfGRho*(

     &            (rho(i,j,k)+rho(i-1,j,k))*(z_r(i,j,k)-z_r(i-1,j,k))

     &   -OneFifth*( (dRx(i,j)-dRx(i-1,j))*( z_r(i,j,k)-z_r(i-1,j,k)
     &                            -OneTwelfth*(dZx(i,j)+dZx(i-1,j)) )

     &              -(dZx(i,j)-dZx(i-1,j))*( rho(i,j,k)-rho(i-1,j,k)
     &                            -OneTwelfth*(dRx(i,j)+dRx(i-1,j)) )
     &                                                            )))
          enddo
        enddo

# ifdef ADV_ISONEUTRAL
        if (CORR_STAGE) then   !<-- corrector stage only
          r0g=rho0/g
          do j=1,ny
            do i=1,i1
#define ISO
#  ifdef ISO
c              dRdx(i,j,k)=-0.5*(pm(i,j)+pm(i-1,j))*rx(i,j)

              dRdx(i,j,k)=0.5*(pm(i,j)+pm(i-1,j))*(
     &                     r0g*0.25*(f(i,j)+f(i-1,j))**2
     &                        *(z_r(i,j,k)-z_r(i-1,j,k))

c    &                  -rx(i,j)

c    &          -0.5*rx(i,j)-0.25*(dRx(i,j)+dRx(i-1,j))

c    &                       -0.5*(dRx(i,j)+dRx(i-1,j))

     &          -0.5*rx(i,j)-0.25*(rx(i-1,j)+rx(i+1,j))

c    &        -0.75*rx(i,j)-0.125*(rx(i-1,j)+rx(i+1,j))

     &                                                )
#  else
              dRdx(i,j,k)=0.5*(pm(i,j)+pm(i-1,j))
     &                 *(z_r(i,j,k)-z_r(i-1,j,k))
#  endif
            enddo
          enddo
        endif
# endif



! ETA-component of pressure gradient term:
!-------------- -- -------- -------- -----

        do j=0,ny+1
          do i=1,nx
            FC(i,j)=(z_r(i,j,k)-z_r(i,j-1,k))
# ifdef MASKING
     &                              *vmask(i,j)
# endif
# ifdef SPLIT_EOS
c**         dpth=0.5*( z_w(i,j,N)+z_w(i,j-1,N)
c**  &                -z_r(i,j,k)-z_r(i,j-1,k))

            dpth=-0.5*(z_r(i,j,k)+z_r(i,j-1,k))

            rx(i,j)=( rho1(i,j,k)-rho1(i,j-1,k)          ! elementary
     &                +(qp1(i,j,k)-qp1(i,j-1,k))         ! adiabatic
     &                     *dpth*(1.-qp2*dpth) )         ! difference
# else
            rx(i,j)=(rho(i,j,k)-rho(i,j-1,k))
# endif
# ifdef MASKING
     &                              *vmask(i,j)
# endif
          enddo
        enddo

        if (jnode==0.and.(.not.ns_periodic)) then
          do i=1,nx
            FC(i,1)=FC(i,2)
            rx(i,1)=rx(i,2)
          enddo
        endif
        if (jnode==np_eta-1) then
          do i=1,nx
            FC(i,ny+1)=FC(i,ny)
            rx(i,ny+1)=rx(i,ny)
          enddo
        endif

        do j=0,ny
          do i=1,nx
            cff=2.*FC(i,j)*FC(i,j+1)
            if (cff > epsil) then
              dZx(i,j)=cff/(FC(i,j)+FC(i,j+1))
            else
              dZx(i,j)=0.
            endif

            cfr=2.*rx(i,j)*rx(i,j+1)
            if (cfr > epsil) then
              dRx(i,j)=cfr/(rx(i,j)+rx(i,j+1))
            else
              dRx(i,j)=0.
            endif
# ifdef SPLIT_EOS
            dRx(i,j)=dRx(i,j) -qp1(i,j,k)*dZx(i,j)
     &                         *(1.+2.*qp2*z_r(i,j,k))
# endif
          enddo               !--> discard FC, rx

          if (j >= 1) then
            do i=1,nx
              rv(i,j,k)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)*(
     &                             P(i,j-1,k)-P(i,j,k) -HalfGRho*(

     &            (rho(i,j,k)+rho(i,j-1,k))*(z_r(i,j,k)-z_r(i,j-1,k))

     &   -OneFifth*( (dRx(i,j)-dRx(i,j-1))*( z_r(i,j,k)-z_r(i,j-1,k)
     &                            -OneTwelfth*(dZx(i,j)+dZx(i,j-1)) )

     &              -(dZx(i,j)-dZx(i,j-1))*( rho(i,j,k)-rho(i,j-1,k)
     &                            -OneTwelfth*(dRx(i,j)+dRx(i,j-1)) )
     &                                                            )))

            enddo
          endif
        enddo   ! <-- j

# ifdef ADV_ISONEUTRAL
        if (CORR_STAGE) then   !<-- corrector stage only
          r0g=rho0/g
          do j=1,j1
            do i=1,nx
#ifdef ISO
c              dRde(i,j,k)=-0.5*(pn(i,j)+pn(i,j-1))*rx(i,j)

              dRde(i,j,k)=0.5*(pn(i,j)+pn(i,j-1))*(
     &                      r0g*0.25*(f(i,j)+f(i,j-1))**2
     &                        *(z_r(i,j,k)-z_r(i,j-1,k))

c    &                -rx(i,j)

c    &            -0.5*rx(i,j)-0.25*(dRx(i,j)+dRx(i,j-1))

c    &                         -0.5*(dRx(i,j)+dRx(i,j-1))

     &            -0.5*rx(i,j)-0.25*(rx(i,j-1)+rx(i,j+1))

c    &            -0.75*rx(i,j)-0.125*(rx(i,j-1)+rx(i,j+1))

     &                                                    )
#else
              dRde(i,j,k)=0.5*(pn(i,j)+pn(i,j-1))
     &                 *(z_r(i,j,k)-z_r(i,j-1,k))
#endif
            enddo
          enddo
        endif
# endif


      enddo   !<-- k

# ifdef DIAGNOSTICS
      if (CORR_STAGE) then
        if ((calc_diag.and.diag_uv)) then
          do k=1,nz
            Udiag(:,:,k,ipgr) = ru(1:nx,1:ny,k)*dxdyi_u*umask(1:nx,1:ny)
            Vdiag(:,:,k,ipgr) = rv(1:nx,1:ny,k)*dxdyi_v*vmask(1:nx,1:ny)
          enddo
        endif
        if (diag_pflx) then
          call calc_diag_pflux(p)
        endif
      endif
# endif

# ifdef NHMG
!
!--------------------------------------------------------------------
! hydrostatic pressure gradient for the vertical flux
!       dpdz = -zx*dpdx - zy*dpdy (See Molemaker et al., 2018)
!--------------------------------------------------------------------
!
!     We are re-using dR as a work array 

      do j=Jstr,Jend
         do k = 1,N
            do i=Istr,Iend
              dR(i,k) = dzdxi(i,j,k) *(                                
     $              ru(i  ,j,k  )/(Hz(i,j,k)+Hz(i-1,j,k)) 
     $            + ru(i+1,j,k  )/(Hz(i,j,k)+Hz(i+1,j,k)))
     $            +     dzdeta(i,j,k) *(                               
     $              rv(i,j  ,k  )/(Hz(i,j,k)+Hz(i,j-1,k)) 
     $            + rv(i,j+1,k  )/(Hz(i,j,k)+Hz(i,j+1,k)))
            enddo
         enddo

         !! rw is a V*dwdt units object (m4/s2), just like ru, rv
         do k = 1,N-1 
            do i=Istr,Iend
              rw(i,j,k)=-0.25*(dR(i,k+1)+dR(i,k))*(Hz(i,j,k+1)+Hz(i,j,k))
            enddo
         enddo
         do i=Istr,Iend
            rw(i,j,N) = -0.5*dR(i,N)*Hz(i,j,N)
         enddo
      enddo
#  ifdef DIAGNOSTICS_NHMG
      if (diag_uv.and.calc_diag) then
        Wdiag(:,:,:,iwprsgr) =rw(1:nx,1:ny,:)*dxdyi
      endif
#  endif
# endif /* NHMG */


#ifdef WEC
! ------------------------------------
!add in wave effect using POM-Jacobian
! ------------------------------------
#  define DC rx
      !! FC contains the vertical vortex force K (eq(5) from uchiyama2010)
      !! which behaves like a hydrostatic pressure
      !! JM We should be able to add these term directly to the
      !! hydrostatic pressure above
      !! In addition, instead of K, we should compute  ust.u_z in place
      do k=N,1,-1           !--> irreversible
        if (k.eq.N) then
          do j=0,ny
            do i=0,nx
              FC(i,j)=-(z_w(i,j,N)-z_r(i,j,N))
     &                          *(1.5*kvf(i,j,N)-0.5*kvf(i,j,N-1))
     &          *(z_w(i,j,N)-z_r(i,j,N))/(z_r(i,j,N)-z_r(i,j,N-1))
              DC(i,j)=-g*sup(i,j)-calP(i,j)+Kapsrf(i,j)  !! additional surface pressure gradient
            enddo
          enddo
        else
          do j=0,ny
            do i=0,nx
              FC(i,j)=FC(i,j)-0.5*(kvf(i,j,k+1)+kvf(i,j,k))
     &                           *(z_r(i,j,k+1)-z_r(i,j,k))
            enddo
          enddo
        endif

        do j=1,ny
          do i=1,nx
            cff  = 0.5*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)
            cff1 = cff*( FC(i-1,j)-FC(i,j) +0.5*( kvf(i,j,k)
     &           +kvf(i-1,j,k) )*(z_r(i,j,k)-z_r(i-1,j,k)) )
            cff2 = cff*(DC(i-1,j)-DC(i,j))
            ru(i,j,k)=ru(i,j,k) + cff1 + cff2
          enddo
        enddo
        do j=1,ny
          do i=1,nx
            cff  = 0.5*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)
            cff1 = cff*( FC(i,j-1)-FC(i,j) +0.5*( kvf(i,j,k)
     &           +kvf(i,j-1,k) )*(z_r(i,j,k)-z_r(i,j-1,k)) )
            cff2 = cff*(DC(i,j-1)-DC(i,j))
            rv(i,j,k)=rv(i,j,k) + cff1 + cff2
          enddo
        enddo
      enddo   !<-- k
#  undef DC
#endif /* WEC */

      nullify(dR,dZx,FC,rx,ru,rv,P)
# ifdef SPLIT_EOS
      nullify(rho)
# endif
# ifdef NHMG
      nullify(rw)
# endif

      end
