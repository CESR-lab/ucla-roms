      module wec_frc
      ! Wave and Current module which provides all the terms that enter
      ! in the momentum equations

       use dimensions !! has dimensions, a list of forcing files and rmask

       implicit none

       private

      ! Variables usable by code that uses this module
       real, public, dimension(nx,ny) :: sup            !! Zeta hat  (many place)
       real, public, dimension(nx,ny) :: brk2dx, brk2de !! show up like wind stress
       real, public, dimension(nx,ny) :: frc2dx,frc2de  !! wave drag, bottom only
       real, public, dimension(nx,ny) :: calP           !! pressure correction (prsgrd)
       real, public, dimension(nx,ny) :: Kapsrf         !! bernouille head?  not sure (prsgrd)
!      real, public, dimension(nx,ny) :: lm             !! wavelength for diffusivity (not used)
       
       real, public, dimension(nx+1,ny)    :: ust2d
       real, public, dimension(nx,ny+1)    :: vst2d
       real, public, dimension(nx-1,ny,nz) :: ust
       real, public, dimension(nx,ny-1,nz) :: vst
       real, public, dimension(nx,ny,nz)   :: brk3dx, brk3de, wst, Kappa,
     &                                     kvf, Akw, frc3dx, frc3de
       real, public, dimension(nx,ny,nz+1) :: Akb

       real, public :: vstu,ustv,dudx,dvdx,dude,dvde


       !logical :: analytical = .false.  ! analytical forcing
       
       ! Can be called from outside this module
       !public ana_wec_tile
       !public wstokes_tile

      contains 
        ! -------------------------------------
        subroutine ana_wec_tile(istr,jstr,iend,jend)
         !analytical wec forcing
         ! Set everything to zero for testing
         implicit none
         real(kind=8), intent(in) :: istr,jstr,iend,jend
         integer :: i, j, k

          ! rho-point variables
          do i=istr,iend
            do j=jstr,jend
               sup(i,j)    = 0.D0
               brk2dx(i,j) = 0.D0
               brk2de(i,j) = 0.D0
               frc2dx(i,j) = 0.D0
               frc2de(i,j) = 0.D0
               calP(i,j)   = 0.D0
               Kapsrf(i,j) = 0.D0
!               lm(i,j)     = 0.D0
               do k=1,nz
                  brk3dx(i,j,k) = 0.D0
                  brk3de(i,j,k) = 0.D0
                  wst(i,j,k)    = 0.D0
                  Kappa(i,j,k)  = 0.D0 !may not need this variable
                  kvf(i,j,k)    = 0.D0
                  Akw(i,j,k)    = 0.D0
                  frc3dx(i,j,k) = 0.D0
                  frc3de(i,j,k) = 0.D0
                enddo
                do k=0,nz
                   Akb(i,j,k) = 0.D0
                enddo
             enddo
          enddo

          !u-point variables
          do i=istr,iend-1
           do j=jstr,jend
              ust2d(i,j) = 0.D0
              do k=1,nz
                 ust(i,j,k) = 0.D0
              enddo
           enddo
          enddo
 
          !v-point variables
          do i=istr,iend
           do j=jstr,jend-1
              vst2d(i,j) = 0.D0
              do k=1,nz
                 vst(i,j,k) =0.D0
              enddo
            enddo
          enddo
 
        end subroutine ana_wec_tile
        ! ----------------------------------------
                 
!        subroutine wstokes_tile(istr,jstr,iend,jend, pm, pn)
         ! Compute vertical stokes velocity         
!         implicit none
!         real(kind=8), intent(in) :: istr,jstr,iend,jend
!         integer :: i, j, k
!         real Wrk(0:nz)

!         do j=jstr,jend
!           do i=istr,iend
!              Wrk(0) = 0.D0
!              do k=1,nz,+1 ! --> recursive
!                  Wrk(k) = Wrk(k-1) -Hz(i,j,k)
!     &       *( pm(i,j)*(ust(i+1,j,k)-ust(i,j,k))
!     &         +pn(i,j)*(vst(i,j+1,k)-vst(i,j,k)) )
!             enddo 
!           wst(i,j,nz)=+0.375*Wrk(nz) +0.75*Wrk(nz-1) -0.125*Wrk(nz-2)
!           do k=nz-1,2,-1
!              wst(i,j,k)=+0.5625*(Wrk(k  )+Wrk(k-1))
!     &                 -0.0625*(Wrk(k+1)+Wrk(k-2))
!           enddo
!           wst(i,j,1)= -0.125*Wrk(2) +0.75*Wrk(1) +0.375*Wrk(0)
!          enddo
!         enddo
!# ifndef EW_PERIODIC
!         if (WESTERN_EDGE) then                       ! Set lateral
!           do k=1,nz                                   ! boundary
!             do j=jstr,jend                           ! conditions
!               wst(istr-1,j,k)=wst(istr,j,k)
!           enddo
!          enddo
!         endif
!        if (EASTERN_EDGE) then
!          do k=1,nz
!            do j=jstr,jend
!              wst(iend+1,j,k)=wst(iend,j,k)
!            enddo
!          enddo
!        endif
!# endif
!# ifndef NS_PERIODIC
!         if (SOUTHERN_EDGE) then
!           do k=1,nz
!             do i=istr,iend
!               wst(i,jstr-1,k)=wst(i,jstr,k)
!             enddo
!          enddo
!         endif
!         if (NORTHERN_EDGE) then
!           do k=1,nz
!             do i=istr,iend
!               wst(i,jend+1,k)=wst(i,jend,k)
!             enddo
!           enddo
!         endif
!#  ifndef EW_PERIODIC
!         if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
!           do k=1,nz
!             wst(istr-1,jstr-1,k)=wst(istr,jstr,k)
!           enddo
!         endif
!         if (WESTERN_EDGE .and. NORTHERN_EDGE) then
!           do k=1,nz
!             wst(istr-1,jend+1,k)=wst(istr,jend,k)
!          enddo
!         endif
!         if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
!           do k=1,nz
!             wst(iend+1, jstr-1,k)=wst(iend,jstr,k)
!           enddo
!         endif
!         if (EASTERN_EDGE .and. NORTHERN_EDGE) then
!           do k=1,nz
!             wst(iend+1,jend+1,k)=wst(iend,jend,k)
!           enddo
!         endif
!#  endif
!# endif
   
!        end subroutine wstokes_tile
        ! ----------------------------------------
        
      end module wec_frc
