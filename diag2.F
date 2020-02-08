#include "cppdefs.h"
#ifdef BIOLOGY

      subroutine bio_diag (tile)

! An alternative diagnostic routine designed to measure quadratic
! variation of physical/biological tracers, mostly for the purpose
! of comparison of different runs.

      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call bio_diag_tile  (istr,iend,jstr,jend)
      end

      subroutine bio_diag_tile (istr,iend,jstr,jend)
      implicit none
      integer istr,iend,jstr,jend, i,j,k,itrc, iocheck,
     &               k1,k2
# include "param.h"
      real FC(PRIVATE_2D_SCRATCH_ARRAY,2),
     &     FX(PRIVATE_2D_SCRATCH_ARRAY),  cff1,
     &     FE(PRIVATE_2D_SCRATCH_ARRAY),  NSUB
      real*QUAD cff, my_sum(0:2*NT+1)
      character*8 echar
# include "grid.h"
# include "ocean3d.h"
# include "scalars.h"

c----# define  ETALON_CHECK
# ifdef ETALON_CHECK
      integer ncheck, nparam
      parameter (ncheck=10, nparam=4)
      integer icheck, check_point(ncheck)
      character*57 check_line, etalon_line(ncheck)
      real    A0(nparam), A1(nparam)
      integer P0(nparam), P1(nparam)

# include "etalon_data.h"
# endif

      if (mod(iic-1,ninfo)==0) then
        do itrc=1,NT
          my_sum(2*itrc  )=0.          ! <-- Reset local sums
          my_sum(2*itrc+1)=0.
          k2=1
          do k=1,N,+1
            k1=k2
            k2=3-k1
            if (k==1) then
              do j=jstr,jend
                do i=istr,iend
                  FC(i,j,k1)=0.
                enddo
              enddo
            endif
            if (k==N) then
              do j=jstr,jend
                do i=istr,iend
                  FC(i,j,k2)=0.
                enddo
              enddo
            else
              do j=jstr,jend
                do i=istr,iend
                 FC(i,j,k2)=(t(i,j,k+1,nstp,itrc)-t(i,j,k,nstp,itrc))
     &                                     /(z_r(i,j,k+1)-z_r(i,j,k))
                enddo
              enddo
            endif
            do j=jstr,jend
              do i=istr,iend+1
                FX(i,j)=(t(i,j,k,nstp,itrc)-t(i-1,j,k,nstp,itrc))
     &                                      *0.5*(pm(i,j)+pm(i-1,j))
# ifdef MASKING
     &                                                   *umask(i,j)
# endif
              enddo
            enddo
            do j=jstr,jend+1
              do i=istr,iend
                FE(i,j)=(t(i,j,k,nstp,itrc)-t(i,j-1,k,nstp,itrc))
     &                                      *0.5*(pn(i,j)+pn(i,j-1))
#ifdef MASKING
     &                                                   *vmask(i,j)
#endif
              enddo
            enddo
            do j=jstr,jend
              do i=istr,iend
              my_sum(itrc)=my_sum(itrc)+cff*(FX(i,j)**2+FX(i+1,j)**2
     &                                      +FE(i,j)**2+FE(i,j+1)**2
     &                                  +FC(i,j,k2)**2+FC(i,j,k1)**2)
              enddo
            enddo
          enddo         ! <-- k
        enddo         !<-- itrc








        do itrc=0,NT
          my_sum(2*itrc  )=0.          ! <-- Reset local sums
          my_sum(2*itrc+1)=0.
        enddo
        do j=jstr,jend
         do i=istr,iend
# ifdef MASKING
           cff1=rmask(i,j)/(pm(i,j)*pn(i,j))  ! <-- grid box area
# else
           cff1=1./(pm(i,j)*pn(i,j))
# endif
           do k=1,N
            cff=cff1*Hz(i,j,k)       ! <-- volume of grid box(i,j,k)
            my_sum(1)=my_sum(1)+cff  ! <-- accumulate volume
             do itrc=1,NT
              my_sum(2*itrc  )=my_sum(2*itrc  )+cff*t(i,j,k,nstp,itrc)
              my_sum(2*itrc+1)=my_sum(2*itrc+1)+cff*t(i,j,k,nstp,itrc)
     &                                             *t(i,j,k,nstp,itrc)
             enddo
           enddo
         enddo
        enddo


        if (SINGLE_TILE_MODE) then
          NSUB=1
        else
          NSUB=NSUB_X*NSUB_E
        endif

! Perform global summation: whoever gets first to the critical region
! resets global sums before global summation starts; after the global
! summation is completed, thread, which is the last one to enter the
! critical region, finalizes the computation of diagnostics and
! prints them out.

C$OMP CRITICAL (diag2_cr_rgn)

        if (tile_count==0) then
          do itrc=1,2*NT+1
            global_sum(itrc)=my_sum(itrc)
          enddo
        else
          do itrc=1,2*NT+1
            global_sum(itrc)=global_sum(itrc)+my_sum(itrc)
            enddo
        endif

        tile_count=tile_count+1

        if (tile_count==NSUB) then
          tile_count=0

          cff=1./global_sum(1)
          do itrc=1,NT
            global_sum(2*itrc  )=cff*global_sum(2*itrc )
            global_sum(2*itrc+1)=cff*global_sum(2*itrc+1)
     &                              -global_sum(2*itrc)**2
          enddo
          if (first_time==0) then
            first_time=1
            write(*,1) 'STEP', 'time[DAYS]', 'rms_PHYT+ZOO',
     &                 'rms____NO3','rms___NH4','rms__DETRIT','trd'
   1        format(1x,A4,3x,A10,1x,A12,4x,A10,4x,A10,4x,A10,3x,A3)
          endif

          do itrc=1,NT
            global_sum(2*itrc+1)=sqrt(global_sum(2*itrc+1))
          enddo

          cff=global_sum(13) + global_sum(15)

          write(*,2) iic-1,tdays, cff, global_sum(7),
     &                  global_sum(9), global_sum(11), proc(2)
   2      format(I6, F12.5, ES16.9, 3(ES14.7), I3)

# ifdef  ETALON_CHECK
          do icheck=1,ncheck
            if (iic-1==check_point(icheck)) then
              write(check_line,10,iostat=iocheck) cff,global_sum(7),
     &                                global_sum( 9),global_sum(11)
  10          format(ES15.9,3(ES14.7))
              if (check_line == etalon_line(icheck)) then
                write(*,*) 'PASSED_ETALON_CHECK'
              else
                read(check_line         ,11,iostat=iocheck) (A1(i),
     &                                           P1(i), i=1,nparam)
                read(etalon_line(icheck),11,iostat=iocheck) (A0(i),
     &                                           P0(i), i=1,nparam)
  11            format(f11.9,1x,i3,3(f10.7,1x,i3))
                do i=1,nparam
                  A1(i)=A1(i)-A0(i)*10.**float(P0(i)-P1(i))
                enddo
                write(check_line        ,12,iostat=iocheck) (A1(i)
     &                                                ,i=1,nparam)
  12            format(f14.9,3(1x,f13.7))
                j=0
                do i=2,56
                  if (check_line(i-1:i)=='0.') then
                    check_line(i:i)=':'
                    j=1
                  else if (j==1 .and. check_line(i:i)=='0') then
                    check_line(i:i)='.'
                  else
                    j=0
                  endif
                enddo
                write(*,15) ' difference:    ',check_line
  15            format(A16,A56)
                write(*,*) 'ETALON_CHECK: SOMETHING_IS_WRONG'
              endif
            endif
          enddo
# endif
        endif
C$OMP END CRITICAL (diag2_cr_rgn)
      endif
      end
#else
      subroutine bio_diag_empty
      end
#endif
