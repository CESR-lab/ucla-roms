! The content of this file is a collection of routines associated
! with vertical coordinate transformation, and sould matches exactly
! the actual transform in ROMS code.

c--#define VERBOSE
#define VERT_COORD_TYPE_SM09

      subroutine set_scoord(theta_s,theta_b, N, Cs_r,Cs_w)
      implicit none
      integer N, k                     ! Setup S-coordinate system:
      real(kind=8) theta_s,   theta_b, ! compute vertical stretching
     &             Cs_w(0:N), Cs_r(N), ! curves Cs_w,Cs_r at W- and
     &                ds, sc, CSF      ! RHO-points.
      ds=1.D0/dble(N)                  ! input:  theta_s, theta_b, N
      Cs_w(N)=0.D0                     ! output: Cs_w, Cs_r
      do k=N-1,1,-1                    !    -1 < Cs_r,Cs_w < 0
        sc=ds*dble(k-N)
        Cs_w(k)=CSF(sc, theta_s,theta_b)
      enddo
      Cs_w(0)=-1.D0
      do k=1,N
        sc=ds*(dble(k-N)-0.5D0)
        Cs_r(k)=CSF(sc, theta_s,theta_b)
      enddo
#ifdef VERBOSE
      write(*,*) 'theta_s=', theta_s, ' theta_b=', theta_b
      write(*,'(12x,A,9x,A,10x,A)') 'Cs_r(k)', 'dCs/ds', 'd^2Cs/ds^2'

      do k=N,1,-1
        write(*,'(4x,I3,3F16.8)')  k, Cs_r(k), (Cs_w(k)-Cs_w(k-1))/ds,
     &                 -4.D0*(Cs_w(k)-2.D0*Cs_r(k)+Cs_w(k-1))/(ds*ds)
      enddo
#endif
      end

      subroutine print_scoord(N, Cs_w, hc,hmin,hmax)
      implicit none
      integer N, k
      real(kind=8) Cs_w(0:N), hc,hmin,hmax,  ds, sc, z1,zhc,z2,z3
      ds=1.D0/dble(N)
      write(*,'(/1x,A/,/2x,A,7x,A/)')
     &        'Vertical S-coordinate system (z at W-points):',
     &             'level   S-coord    Cs-curve    Z at hmin',
     &                        'at hc    half way     at hmax'
      do k=N,0,-1
        sc=ds*dble(k-N)
# ifdef VERT_COORD_TYPE_SM09
        z1=hmin*(hc*sc + hmin*Cs_w(k))/(hc+hmin)
        zhc=0.5*hc*(sc + Cs_w(k))
        z2=0.5*hmax*(hc*sc + 0.5*hmax*Cs_w(k))/(hc+0.5*hmax)
        z3=hmax*(hc*sc + hmax*Cs_w(k))/(hc+hmax)
# else
        zhc=hc*sc
        z1=zhc +     (hmin-hc)*Cs_w(k)
        z2=zhc + (0.5*hmax-hc)*Cs_w(k)
        z3=zhc +     (hmax-hc)*Cs_w(k)
# endif
        if (hc < 1.E+4) then
          write(*,'(I7,F11.6,F12.7,4F12.3)')
     &               k, ds*(k-N),Cs_w(k), z1,zhc,z2,z3
        else
          write(*,'(I7,F11.6,F12.7,F12.3,12x,2F12.3)')
     &               k, ds*(k-N),Cs_w(k), z1,    z2,z3
        endif
      enddo
      end



! Vertical stretching functions: In principle stretching curves can be
! selected independently from the vertical transformation type, which
! would require a separate name for CPP-switch here.  This however has
! little practical value (merely to reproduce my own Pacific solutions
! of 2005 KPP presentation where the SM09 coordinate transformation
! formula was used in combination with legacy SH94 stretching curves),
! so for simplicity selection of SM09 stretching curve is tied with
! SM09 coordinate transform.



#ifdef VERT_COORD_TYPE_SM09
                                               ! Note that mathematical
      function CSF(sc, theta_s,theta_b)        ! limits of csrf,CSF for
      implicit none                            ! theta_s, theta_b --> 0
      real(kind=8) CSF,sc,theta_s,theta_b,csrf ! match that under else
                                               ! logical branches.
      if (theta_s > 0.D0) then
        csrf=(1.D0-cosh(theta_s*sc))/(cosh(theta_s)-1.D0)
      else
        csrf=-sc**2
      endif
      if (theta_b > 0.D0) then
        CSF=(exp(theta_b*csrf)-1.D0)/(1.D0-exp(-theta_b))
      else
        CSF=csrf                               ! Reference: This form of
      endif                                    ! CSF exactly corresponds
      end                                      ! to Eq.(2.4) from SM2009


#elif defined THIS_SHOULD_BE_NEVER_DEFINED
                                               ! Alternative variants
      function CSF (sc, theta_s,theta_b)       ! for reference only
      implicit none
      real(kind=8) CSF, sc, theta_s,theta_b, sc1, csrf,cbot, wgth
      real(kind=8), parameter :: alfa=1.D0, beta=1.D0, albt=alfa/beta

      if (theta_s > 0.D0) then
        csrf=(1.D0-cosh(theta_s*sc))/(cosh(theta_s)-1.D0)
      else
        csrf=-sc**2
      endif

      sc1=sc+1.D0
      if (theta_b > 0.D0) then
        cbot=sinh(theta_b*sc1)/sinh(theta_b)-1.D0
      else
        cbot=sc
      endif

      wgth=sc1**alfa * ( 1.D0+ albt*(1.D0-sc1**beta) )
      CSF=wgth*csrf + (1.D0-wgth)*cbot
      end


      function CSF (sc, theta_s,theta_b)
      implicit none
      real(kind=8) CSF, sc, theta_s,theta_b, sc1, csrf,cbot, wgth
      real(kind=8), parameter :: alfa=1.D0, beta=1.D0, albt=alfa/beta

      if (theta_s > 0.D0) then
        csrf=(1.D0-cosh(theta_s*sc))/(cosh(theta_s)-1.D0)
        if (theta_b > 0.D0) then
          sc1=sc+1.D0
          cbot=sinh(theta_b*sc1)/sinh(theta_b)-1.D0
          wgth=sc1**alfa *( 1.D0 +  albt*( 1.D0-sc1**beta ))
          CSF=wgth*csrf + (1.D0-wgth)*cbot
        else
          CSF=csrf
        endif
      else
        CSF=sc
      endif
      end

# else
                                                    ! Legacy SH94
      function CSF(sc, theta_s,theta_b)             ! S-coordinate
      real(kind=8) CSF, sc, theta_s,theta_b         ! stretching
                                                    ! formula
      if (theta_s > 0.D0) then
        CSF=(1.D0-theta_b)*sinh(theta_s*sc)/sinh(theta_s)
     &           +theta_b*( 0.5D0*tanh(theta_s*(sc+0.5D0))
     &                       /tanh(0.5*theta_s)   -0.5D0 )
      else
        CSF=sc
      endif
      end

# endif





! The following routine just retrieves Cs_r,Cs_w, hc, VertCoordType
! from a previously opened necCDF file. Note that it is preferable to
! read Cs-curves from the file (here saved as attributes) rather than
! try to recompute them from transform parameters "theta_s, theta_b"
! because the specific details of the transformation may be changed
! at any time.


      subroutine read_scoord(ncid, N, Cs_r, Cs_w, hc, VertCoordType)
      implicit none
      integer ncid, N, ierr
      real(kind=8) Cs_r(N), Cs_w(0:N), hc
      character(len=*) VertCoordType
#include "netcdf.inc"

      VertCoordType=' '
!     ierr=nf_get_att_text (ncid, nf_global, 'VertCoordType',
!    &                                        VertCoordType)
      VertCoordType='NEW'
      ierr=nf_noerr
      if (ierr == nf_noerr) then
        write(*,*) 'VertCoordType=', VertCoordType
        if (VertCoordType=='NEW') VertCoordType='SM09'
      else
        write(*,*) 'Global attribute VertCoordType not found.'
        ierr=nf_get_att_text (ncid, nf_global, 'sc_r', Cs_r)
        if (ierr == nf_noerr) then
          write(*,*) 'Found  ''sc_r'':==>  old s-coord.'
        else
          write(*,*) '''sc_r'' not found ==> VertCoordType=''SM09''.'
          VertCoordType='SM09'
        endif
      endif

      ierr=nf_get_att_double (ncid, nf_global, 'Cs_r', Cs_r)
      if (ierr == nf_noerr) then
        ierr=nf_get_att_double (ncid, nf_global, 'Cs_w', Cs_w)
        if (ierr == nf_noerr) then
          ierr=nf_get_att_double (ncid, nf_global, 'hc',  hc)
          if (ierr == nf_noerr) then

            write(*,*) 'read Cs_r,Cs_w,hc from netCDF file.'
            return  !--> successful return

          else
            write(*,'(/1x,2A/12x,A/)')      '### ERROR: Cannot get ',
     &                  'global attribute ''hc''', nf_strerror(ierr)
          endif
        else
          write(*,'(/1x,2A/12x,A/)') '### ERROR: Cannot get global ',
     &                       'attribute ''Cs_w''', nf_strerror(ierr)
        endif
      else
        write(*,'(/1x,2A/12x,A/)')   '### ERROR: Cannot get global ',
     &                       'attribute ''Cs_r''', nf_strerror(ierr)
      endif
      stop
      end


      subroutine set_depth(Lm,Mm,N, Cs_r,Cs_w,hc, VertCoordType,
     &                                               h, z_r,z_w)
      implicit none
      integer Lm,Mm,N, i,j,k, ierr
      real(kind=8) h(0:Lm+1,0:Mm+1),       hinv(0:Lm+1), hc,     ds,
     &             z_r(0:Lm+1,0:Mm+1,N),   Cs_r(N),      Cs_r_k, cff_r,
     &             z_w(0:Lm+1,0:Mm+1,0:N), Cs_w(0:N),    Cs_w_k, cff_w
      character(len=*) VertCoordType

#ifdef VERBOSE
      integer lenstr
      write(*,'(1x,4A,F12.4,1x,A,3I5)')           'enter set_depth: ',
     &   'VertCoordType = ''', VertCoordType(1:lenstr(VertCoordType)),
     &                           ''' hc =', hc, 'Lm,Mm,N =', Lm,Mm,N
c**   write(*,'(1x,A,80F10.6)') 'Cs_r =', (Cs_r(k), k=1,N)
#endif
      if (VertCoordType == 'SM09') then
        write(*,*) 'new VertCoordType = SM09'
      else
        write(*,*) 'defaulting to SH94 S-coordinate'
      endif

      ds=1.D0/dble(N)
      do j=0,Mm+1
        if (VertCoordType == 'SM09') then
          Cs_r_k=Cs_r(N) ;  cff_r=-0.5D0*hc*ds
          do i=0,Lm+1
            z_w(i,j,N)=0.D0
            hinv(i)=h(i,j)/(hc+h(i,j))
            z_r(i,j,N)=hinv(i)*( cff_r + Cs_r_k*h(i,j) )
            z_w(i,j,0)=-h(i,j)
          enddo
          do k=N-1,1,-1
            Cs_w_k=Cs_w(k)        ; Cs_r_k=Cs_r(k)
            cff_w=hc*ds*dble(k-N) ; cff_r=hc*ds*(dble(k-N)-0.5D0)
            do i=0,Lm+1
              z_w(i,j,k)=hinv(i)*( cff_w + Cs_w_k*h(i,j) )
              z_r(i,j,k)=hinv(i)*( cff_r + Cs_r_k*h(i,j) )
            enddo
          enddo
        else
          ierr=0
          Cs_r_k=Cs_r(N)  ; cff_r=hc*(-0.5D0*ds -Cs_r_k)
          do i=0,Lm+1
            if (hc > h(i,j)) ierr=1
            z_w(i,j,N)=0.D0
            z_r(i,j,N)=cff_r + Cs_r_k*h(i,j)
            z_w(i,j,0)=-h(i,j)
          enddo
          if (ierr /= 0) then
            write(*,'(/1x,2A/)')  '### ERROR: set_depth: hc > hmin ',
     &             'occurrence while selecting legacy s-coordinate.'
            stop
          endif
          do k=N-1,1,-1
            Cs_w_k=Cs_w(k) ; cff_w=hc*(ds* dble(k-N)        -Cs_w_k)
            Cs_r_k=Cs_r(k) ; cff_r=hc*(ds*(dble(k-N)-0.5D0) -Cs_r_k)
            do i=0,Lm+1
              z_w(i,j,k)=cff_w + Cs_w_k*h(i,j)
              z_r(i,j,k)=cff_r + Cs_r_k*h(i,j)
            enddo
          enddo
        endif
      enddo !<-- j
#ifdef VERBOSE
      write(*,'(1x,A)')      '...leaving set_depth'
#endif
      end



      subroutine compute_ubar(Lm,Mm,N, z_w, u,ubar)
      implicit none
      integer Lm,Mm,N, i,j,k
      real(kind=8) z_w(0:Lm+1,0:Mm+1,0:N)
      real(kind=4) u(1:Lm+1,0:Mm+1,N), ubar(1:Lm+1,0:Mm+1)
      do j=0,Mm+1
        do i=1,Lm+1
          ubar(i,j)=0.0
        enddo
        do k=N,1,-1
          do i=1,Lm+1
            ubar(i,j)=ubar(i,j) + u(i,j,k)*( z_w(i,j,k)+z_w(i-1,j,k)
     &                                  -z_w(i,j,k-1)-z_w(i-1,j,k-1))
          enddo
        enddo
        do i=1,Lm+1
          ubar(i,j)=ubar(i,j)/( z_w(i,j,N)+z_w(i-1,j,N)
     &                         -z_w(i,j,0)-z_w(i-1,j,0))
        enddo
      enddo
      end

      subroutine compute_vbar (Lm,Mm,N, z_w, v,vbar)
      implicit none
      integer Lm,Mm,N, i,j,k
      real(kind=8) z_w(0:Lm+1,0:Mm+1,0:N)
      real(kind=4) v(0:Lm+1,1:Mm+1,N), vbar(0:Lm+1,1:Mm+1)
      do j=1,Mm+1
        do i=0,Lm+1
          vbar(i,j)=0.0
        enddo
        do k=N,1,-1
          do i=0,Lm+1
            vbar(i,j)=vbar(i,j) + v(i,j,k)*( z_w(i,j,k)+z_w(i,j-1,k)
     &                                  -z_w(i,j,k-1)-z_w(i,j-1,k-1))
          enddo
        enddo
        do i=0,Lm+1
          vbar(i,j)=vbar(i,j)/( z_w(i,j,N)+z_w(i,j-1,N)
     &                         -z_w(i,j,0)-z_w(i,j-1,0))
        enddo
      enddo
      end
