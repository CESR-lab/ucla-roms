      module scoord

      use dimensions
      use param
      use netcdf
      use roms_mpi
      use comm_vars

      implicit none

      private

! Control parameters for vertical coordinate transform: "theta_s" and
! "theta_b" are surface and bottom refinement coefficients for Cs=Cs(s)
! curves; "hc" is critical depth [meters] above which the vertical grid
! spacing remains approximately uniform while becoming stretched below.
! Once Cs(s) is specified, the unperturbed (zeta=0) vertical coordinate
! transformation z=z(s) is defined as either,
!
!               z = hc*s + (h-hc)*Cs(s)                  SH94
! or
!                        hc*s + h*Cs(s)
!               z = h * ----------------                 SM09
!                            hc + h
!
! depending on setting VERT_COORD_TYPE_XXXX. The upper transform has
! the limitation of
!
!                       hc < hmin
!
! (otherwise z=z(s) looses monotonicity resulting in negative grid
! spacings) which severely restricts the choice of hc; in the second
! case "hc" can be selected independently from minimum depth "hmin".
!
! Setting theta_s=theta_b=0 in the case of SH94 yields vertically
! uniform spacing (plain sigma); For both cases increase of "theta_s"
! leads to a more concentrated resolution toward the surface; typically
! 0=<theta_s<10 for both types; Optimally chosen "hc" should be
! comparable to the expected thermocline depth.
! Arrays "Cs_w" and "Cs_r" are stretching curves, Cs=Cs(s), where "s"
! is sigma-coordinate, -1 < s < 0, while "_w" and "_r" stand for
! vertical W- and RHO-point locations.

      real,public :: theta_s,theta_b,hc
      real,public,dimension(nz)  :: Cs_r
      real,public,dimension(0:nz):: Cs_w

      public set_scoord
      public check_scoord

      contains

! ---------------------------------------------------------------
      subroutine set_scoord ![
      ! Takes theta_s,theta_b, hc which are initialized
      ! in read_inp and produces Cs_r and Cs_w

      implicit none                      ! output: Cs_w(0:N), Cs_r(1:N)

      ! Local
      integer k                          ! Compute vertical stretching
      real ds,sc, z1,zhc,z2,z3           ! curves at W- and RHO-points,

      ds=1.D0/dble(nz)                   !      -1 < Cs_r,Cs_w < 0,
      Cs_w(nz)=0.D0
      do k=nz-1,1,-1                     ! then print out z-coordinate
        sc=ds*dble(k-nz)                 ! of S-levels for minimum,
        Cs_w(k)=CSF(sc, theta_s,theta_b) ! h=hc, half way, and maximum
      enddo                              ! depth.
      Cs_w(0)=-1.D0
      do k=1,nz
        sc=ds*(dble(k-nz)-0.5D0)
        Cs_r(k)=CSF(sc, theta_s,theta_b)
      enddo


      if (mynode==0) write(*,'(/1x,A/)')
     &        'Vertical S-coordinate system (z at W-points):'

      if (mynode==0) write(*,
     &     '(A7,2x,A11,2x,A14,2x,A16,2x,A14,2x,A14,2x,A14)')
     &     'level', 'S-coord', 'Cs-curve', 'Z at hmin',
     &     'at hc', 'half way', 'at hmax'


      do k=nz,0,-1
        sc=ds*dble(k-nz)
        z1=hmin*(hc*sc + hmin*Cs_w(k))/(hc+hmin)
        zhc=0.5*hc*(sc + Cs_w(k))
        z2=0.5*hmax*(hc*sc + 0.5*hmax*Cs_w(k))/(hc+0.5*hmax)
        z3=hmax*(hc*sc + hmax*Cs_w(k))/(hc+hmax)

        if (mynode==0)  write(*,
     &'(I7,2X,ES11.4,2X,ES14.7,2X,ES16.7,2X,ES14.7,2X,ES14.7,2X,ES14.7)')
     &       k, ds*(k-nz), Cs_w(k), z1, zhc, z2, z3
      enddo

      end subroutine set_scoord !]
! ---------------------------------------------------------------
      function CSF(sc, theta_s,theta_b) ![

      ! This form of CSF corresponds to eq. (2.4) from
      ! Shchepetkin and McWilliams (2009) J. Comput. Phys.

      implicit none

      ! import/export
      real,intent(in) :: sc,theta_s,theta_b
      real :: CSF

      ! local
      real :: csrf


      if (theta_s > 0.D0) then
        csrf=(1.D0-cosh(theta_s*sc))/(cosh(theta_s)-1.D0)
      else
        csrf=-sc**2
      endif
      if (theta_b > 0.D0) then
        CSF=(exp(theta_b*csrf)-1.D0)/(1.D0-exp(-theta_b))
      else
        CSF=csrf
      endif
      end function CSF !]
! ---------------------------------------------------------------
      subroutine check_scoord(fname)  ![
      ! Verify that the Cs_w/Cs_r coefficients are matching the
      ! stretching parameters that are used for this run
      ! Use for bry and ini data files
      ! must run set_scoord first

      implicit none

      ! import
      character(len=*) fname

      ! local

      integer        :: ncid,ierr,varid
      real,parameter :: eps =1.0e-6
      real,dimension(:),allocatable :: dat

      ierr = nf90_open(trim(fname), nf90_nowrite, ncid)
      if (ierr /= 0) then
        if (mynode == 0) print *, 'ERROR: cannot open ', trim(fname)
        return    ! or ROMS_ABORT per the error-handling design
      endif

      ierr = nf90_inq_varid(ncid,'Cs_w',varid)
      if (ierr/=0) then
        if (mynode.eq.0) print *,'WARNING: Cs_w not present in ',trim(fname)
      else
        allocate(dat(0:nz))
        ierr = nf90_get_var(ncid,varid,dat)
        if (ierr /= 0) then
          if (mynode==0) print *, 'ERROR: failed to read Cs_w from ', trim(fname)
          stop
        endif
        if (maxval(abs(Cs_w-dat))>eps) then
          if (mynode==0)  write(*,'(/1x,4A/)') '### ERROR: ',
     &      'Mismatch in ''Cs_w'' ', 'values in file ''', trim(fname), '''.'
          stop
        endif
        deallocate(dat)
      endif

      ierr = nf90_inq_varid(ncid,'Cs_r',varid)
      if (ierr/=0) then
        if (mynode.eq.0) print *,'WARNING: Cs_r not present in ',trim(fname)
      else
        allocate(dat(nz))
        ierr = nf90_get_var(ncid,varid,dat)
        if (ierr /= 0) then
          if (mynode==0) print *, 'ERROR: failed to read Cs_r from ', trim(fname)
          stop
        endif
        if (maxval(abs(Cs_r-dat))>eps) then
          if (mynode==0)  write(*,'(/1x,4A/)') '### ERROR: ',
     &      'Mismatch in ''Cs_r'' ', 'values in file ''', trim(fname), '''.'
          stop
        endif
        deallocate(dat)
      endif

      ierr = nf90_close(ncid)

      end  subroutine check_scoord !]
! ---------------------------------------------------------------

      end module scoord
