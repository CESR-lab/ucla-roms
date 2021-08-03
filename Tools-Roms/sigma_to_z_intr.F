#define CUBIC

      module s_to_z_scratch
        integer(kind=4), allocatable, dimension(:) :: km
C$OMP THREADPRIVATE(km)
        real(kind=8), allocatable, dimension(:) :: zz
C$OMP THREADPRIVATE(zz)
#ifdef CUBIC
        real(kind=8), allocatable, dimension(:) :: dz,FC
C$OMP THREADPRIVATE(dz,FC)
#endif
      end module s_to_z_scratch


      subroutine sigma_to_z_thread(Lm,Mm,N,  z_r,z_w, rmask,  nz,z_lev,
     &                          imin,jmin,kmin, var,var_zlv, FillValue)

! Interpolate field "var" defined in sigma-space to a set of several
! fixed z-levels z_lev(1:nz).

      use s_to_z_scratch

      implicit none
      integer Lm,Mm,N, nz, imin,jmin,kmin
      real(kind=8) z_r(0:Lm+1,0:Mm+1,N), z_w(0:Lm+1,0:Mm+1,0:N),
     &                                     rmask(0:Lm+1,0:Mm+1)
      real(kind=4) var(imin:Lm+1,jmin:Mm+1,kmin:N), z_lev(nz),
     &             var_zlv(imin:Lm+1,jmin:Mm+1,nz), FillValue

!--> internal

      integer numthreads,trd, nsub_x,nsub_y, my_first,my_last,
     &           range, tile, istr,iend,jstr,jend,  size,siz3d
C$    integer omp_get_num_threads, omp_get_thread_num
      numthreads=1 ; trd=0
C$    numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()

      call set_tiles(Lm,Mm, nsub_x,nsub_y)
      range=(nsub_x*nsub_y +numthreads-1)/numthreads
      my_first=trd*range
      my_last=min(my_first + range-1, nsub_x*nsub_y-1)


! Allocate thread-private scratch space: (1) only once during execution
! time; (2) the scratch size may be too large to use automatic arrays.

      if (.not.allocated(km)) then
C$OMP CRITICAL(alloc_scr_cr_rgn)
        size=(Lm+2 +nsub_x-1)/nsub_x ; siz3d=size*(N+2)
        allocate(km(siz3d))
        size=size + 2*siz3d  !<-- actually allocated size expressed
        allocate(zz(siz3d))  !          in would be 4-byte numbers
      
#ifdef CUBIC
        allocate(dz(siz3d),FC(siz3d))
        size=size + 4*siz3d
#endif
        write(*,'(1x,A,F10.5,A,I3)')  'sigma_to_z_thread :: allocated',
     &  dble(size)/dble(512*512),' MBytes of scratch space, trd =',trd
C$OMP END CRITICAL(alloc_scr_cr_rgn)
      endif


      do tile=my_first,my_last
        call comp_tile_bounds(tile, Lm,Mm, nsub_x,nsub_y,  istr,iend,
     &                                                     jstr,jend)

        if (istr==1) istr=max(istr-1,imin) ; if (iend==Lm) iend=iend+1
        if (jstr==1) jstr=max(jstr-1,jmin) ; if (jend==Mm) jend=jend+1

        call sigma_to_z_tile( istr,iend,jstr,jend, Lm,Mm,N,
     &                  z_r,z_w,rmask, nz,z_lev, imin,jmin,kmin,
     &                           var,var_zlv, FillValue,  km,zz
#ifdef CUBIC
     &                                                  , dz,FC
#endif
     &                                                        )
      enddo
      end


      subroutine sigma_to_z_tile( istr,iend,jstr,jend, Lm,Mm,N,
     &                  z_r,z_w,rmask, nz,z_lev, imin,jmin,kmin,
     &                           var,var_zlv, FillValue,  km,zz
#ifdef CUBIC   
     &                                                  , dz,FC 
#endif
     &                                                        )
      implicit none
      integer istr,iend,jstr,jend, Lm,Mm,N, nz,  imin,jmin,kmin
      real(kind=8) z_r(0:Lm+1,0:Mm+1,N),   rmask(0:Lm+1,0:Mm+1),
     &             z_w(0:Lm+1,0:Mm+1,0:N), zlev
      real(kind=4) var(imin:Lm+1,jmin:Mm+1,kmin:N), z_lev(nz),
     &             var_zlv(imin:Lm+1,jmin:Mm+1,nz), FillValue
      integer km(istr:iend), i,j,k,m
      real(kind=8) zz(istr:iend,  0:N+1), dpth
#ifdef CUBIC
     &           , dz(istr:iend,kmin-1:N)
     &           , FC(istr:iend,kmin-1:N), p,q,cff
#endif

! Interpolate field "var" defined in sigma-space to a set of several
! fixed z-levels z_lev(1:nz).  Below imin,jmin,kmin assuming values of
! 0 or 1 each are used to determine variable placement on the staggered
! C-grid: 0,0,* are for horizontal RHO-points; 1,0,* for U-; 0,1,* for
! V-; and 1,1,* for vorticity; the third one *,*,1 for vertical RHO-
! and *,*,0 for W-points. 

      do j=jstr,jend
        if (imin == 0 .and. jmin == 0 .and. kmin == 0) then
          do k=0,N
            do i=istr,iend
              zz(i,k)=z_w(i,j,k)
            enddo
          enddo
        elseif (imin == 1 .and. jmin == 0 .and. kmin == 0) then
          do k=0,N
            do i=istr,iend
              zz(i,k)=0.5D0*(z_w(i,j,k)+z_w(i-1,j,k))
            enddo
          enddo
        elseif (imin == 0 .and. jmin == 1 .and. kmin == 0) then
          do k=0,N
            do i=istr,iend
              zz(i,k)=0.5*(z_w(i,j,k)+z_w(i,j-1,k))
            enddo
          enddo
        elseif (imin == 1 .and. jmin == 1 .and. kmin == 0) then
          do k=0,N
            do i=istr,iend
              zz(i,k)=0.25D0*( z_w(i,j,k)+z_w(i-1,j,k)
     &                    +z_w(i,j-1,k)+z_w(i-1,j-1,k))
            enddo
          enddo

        elseif (imin == 0 .and. jmin == 0 .and. kmin == 1) then
          do k=1,N
            do i=istr,iend
              zz(i,k)=z_r(i,j,k)
            enddo
          enddo
          do i=istr,iend
            zz(i,N+1)=z_w(i,j,N)
            zz(i,  0)=z_w(i,j,0)
          enddo
        elseif (imin == 1 .and. jmin == 0 .and. kmin == 1) then
          do k=1,N
            do i=istr,iend
              zz(i,k)=0.5D0*(z_r(i,j,k)+z_r(i-1,j,k))
            enddo
          enddo
          do i=istr,iend
            zz(i,N+1)=0.5D0*(z_w(i-1,j,N)+z_w(i,j,N))
            zz(i,  0)=0.5D0*(z_w(i-1,j,0)+z_w(i,j,0))
          enddo
        elseif (imin == 0 .and. jmin == 1 .and. kmin == 1) then
          do k=1,N
            do i=istr,iend
              zz(i,k)=0.5*(z_r(i,j,k)+z_r(i,j-1,k))
            enddo
          enddo
          do i=istr,iend
            zz(i,N+1)=0.5D0*(z_w(i,j,N)+z_w(i,j-1,N))
            zz(i,  0)=0.5D0*(z_w(i,j,0)+z_w(i,j-1,0))
          enddo
        elseif (imin == 1 .and. jmin == 1 .and. kmin == 1) then
          do k=1,N
            do i=istr,iend
              zz(i,k)=0.25D0*( z_r(i,j,k)+z_r(i-1,j,k)
     &                    +z_r(i,j-1,k)+z_r(i-1,j-1,k))
            enddo
          enddo
          do i=istr,iend
            zz(i,N+1)=0.25D0*( z_w(i,j,N)+z_w(i-1,j,N)
     &                      +z_w(i,j-1,N)+z_w(i-1,j-1,N))

            zz(i,  0)=0.25D0*( z_w(i,j,0)+z_w(i-1,j,0)
     &                      +z_w(i,j-1,0)+z_w(i-1,j-1,0))
          enddo
        endif

#ifdef CUBIC
        do k=kmin,N-1
          do i=istr,iend
            dz(i,k)=zz(i,k+1)-zz(i,k)
            FC(i,k)=var(i,j,k+1)-var(i,j,k)
          enddo
        enddo
        do i=istr,iend
          dz(i,kmin-1)=dz(i,kmin)
          FC(i,kmin-1)=FC(i,kmin)

          dz(i,N)=dz(i,N-1)
          FC(i,N)=FC(i,N-1)
        enddo
        do k=N,kmin,-1   !--> irreversible
          do i=istr,iend
            cff=FC(i,k)*FC(i,k-1)
            if (cff > 0.D0) then
              FC(i,k)=cff*(dz(i,k)+dz(i,k-1))/( (FC(i,k)+FC(i,k-1))
     &                                          *dz(i,k)*dz(i,k-1) )
            else
              FC(i,k)=0.D0
            endif
          enddo
        enddo
#endif
        do m=1,nz
          zlev=z_lev(m)
          if (kmin == 0) then                      !
            do i=istr,iend
              dpth=zz(i,N)-zz(i,0)
              if (rmask(i,j) < 0.5) then
                km(i)=-3          !--> masked out
              elseif (dpth*(zlev-zz(i,N)) > 0.) then
                km(i)=N+2         !<-- above surface
              elseif (dpth*(zz(i,0)-zlev) > 0.) then
                km(i)=-2          !<-- below bottom
              else
                km(i)=-1          !--> to search
              endif
            enddo
          else
            do i=istr,iend
              dpth=zz(i,N+1)-zz(i,0)
              if (rmask(i,j) < 0.5) then
                km(i)=-3          !--> masked out
              elseif (dpth*(zlev-zz(i,N+1)) > 0.) then
                km(i)=N+2         !<-- above surface

              elseif (dpth*(zlev-zz(i,N)) > 0.) then
                km(i)=N           !<-- below surface, but above z_r(N)
              elseif (dpth*(zz(i,0)-zlev) > 0.) then
                km(i)=-2          !<-- below bottom
              elseif (dpth*(zz(i,1)-zlev) > 0.) then
                km(i)=0           !<-- above bottom, but below z_r(1)
              else
                km(i)=-1          !--> to search
              endif
c       write(*,*) z_w(i,j,0), zz(i,0), zlev !, zz(i,N)
            enddo
          endif

c        write(*,*) 'm =',m, ' zlev =', zlev
c         do i=istr,iend
c           if (km(i) == -1) write(*,*) 'to search'
c         enddo


          do k=N-1,kmin,-1
            do i=istr,iend
              if (km(i) == -1) then
                if ((zz(i,k+1)-zlev)*(zlev-zz(i,k)) >= 0.) then
                  km(i)=k
c                 write(*,*) 'found, km(i)=',  km(i)
                endif
              endif
            enddo
          enddo

          do i=istr,iend
            if (km(i) == -3) then
              var_zlv(i,j,m)=0.             !<-- masked out
            elseif (km(i) == -2) then
              var_zlv(i,j,m)=FillValue      !<-- below bottom
            elseif (km(i) == N+2) then
              var_zlv(i,j,m)=-FillValue     !<-- above surface
            elseif (km(i) == N) then
              var_zlv(i,j,m)=var(i,j,N)     !-> R-point, above z_r(N)
#ifdef CUBIC
     &              +FC(i,N)*(zlev-zz(i,N))
#else
     &              +(zlev-zz(i,N))*(var(i,j,N)-var(i,j,N-1))
     &                                   /(zz(i,N)-zz(i,N-1))
#endif
            elseif (km(i) == kmin-1) then   !-> R-point below z_r(1),
              var_zlv(i,j,m)=var(i,j,kmin)  !     but above bottom
#ifdef CUBIC
     &           -FC(i,kmin)*(zz(i,kmin)-zlev)
#else
     &          -(zz(i,kmin)-zlev)*(var(i,j,kmin+1)-var(i,j,kmin))
     &                                  /(zz(i,kmin+1)-zz(i,kmin))
#endif
            else
              k=km(i)
#ifdef CUBIC
              cff=1.D0/(zz(i,k+1)-zz(i,k))
              p=zlev-zz(i,k) ; q=zz(i,k+1)-zlev

              var_zlv(i,j,m)=cff*( q*var(i,j,k) + p*var(i,j,k+1)
     &             -cff*p*q*( cff*(q-p)*(var(i,j,k+1)-var(i,j,k))
     &                                  +p*FC(i,k+1) -q*FC(i,k) )
     &                                                          )
#else
              var_zlv(i,j,m)=( var(i,j,k)*(zz(i,k+1)-zlev)
     &                        +var(i,j,k+1)*(zlev-zz(i,k))
     &                               )/(zz(i,k+1)-zz(i,k))
c               write(*,*) var_zlv(i,j,m)
#endif
#ifdef VERBOSE
              if ((zz(i,k+1)-zlev)*(zlev-zz(i,k)) < 0.D0) then
                write(*,'(/1x,2A/)') '### ERROR: sigma_to_z_tile ',
     &                               ':: Search algoritm failure.'
                stop
              endif
#endif
            endif
          enddo
        enddo ! <-- m
      enddo !<-- j
      end
