#include "cppdefs.h"
#ifdef BIOLOGY

      subroutine bio_diag (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call bio_diag_tile  (istr,iend,jstr,jend)
      end

      subroutine bio_diag_tile (istr,iend,jstr,jend)
      implicit none
      integer istr,iend,jstr,jend, i,j,k,itrc, iocheck
# include "param.h"
      real cff1, NSUB
      real*QUAD cff, my_sum(0:2*NT+1)


# include "grid.h"
# include "ocean3d.h"
# include "scalars.h"
# include "comm_vars.h"

# define  ETALON_CHECK
# ifdef ETALON_CHECK
      integer ncheck, nparam
      parameter (ncheck=10, nparam=4)
      integer icheck, check_point(ncheck)
      character*57 check_line, etalon_line(ncheck)
      character*1 symbol(57)
      pointer (isym,symbol)
      real    A0(nparam), A1(nparam)
      integer P0(nparam), P1(nparam)

# include "etalon_data.PACIFIC"
# endif

      if (mod(iic-1,ninfo)==0) then
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

! Perform global summation: whoever gets first to the critical region
! resets global sums before global summation starts; after the global
! summation is completed, thread, which is the last one to enter the
! critical region, finalizes the computation of diagnostics and
! prints them out.

        if (SINGLE_TILE_MODE) then
          NSUB=1
        else
          NSUB=NSUB_X*NSUB_E
        endif

C$OMP CRITICAL (bio_cr_rgn)
          if (bc_count==0) then
            do itrc=1,2*NT+1
              global_sum(itrc)=QuadZero   ! <-- Reset global sums
            enddo
          endif
                                       !-->  Perform global summation
          do itrc=1,2*NT+1
            global_sum(itrc)=global_sum(itrc)+my_sum(itrc)
          enddo

          bc_count=bc_count+1             ! This counter identifies
          if (bc_count==NSUB) then      ! the last thread, whoever
            bc_count=0                    ! it is, not always master.
# ifdef MPI
  Reduce sequence here
# endif

            cff=1./global_sum(1)
            do itrc=1,NT
              global_sum(2*itrc  )=cff*global_sum(2*itrc )
              global_sum(2*itrc+1)=cff*global_sum(2*itrc+1)
     &                                -global_sum(2*itrc)**2
            enddo
            if (first_time==0) then
              first_time=1
              write(*,1) 'STEP','time[DAYS]','rms_PHYT+ZOO',
     &                  'rms____NO3','rms___NH4','rms__DETRIT','trd'
   1          format(1x,A4,3x,A10,1x,A12,4x,A10,4x,A10,4x,A10,3x,A3)
            endif

            do itrc=1,NT
              global_sum(2*itrc+1)=global_sum(2*itrc+1)**0.5
            enddo

            cff=global_sum(13) + global_sum(15)

            write(*,2) iic-1,tdays, cff, global_sum(7),
     &                    global_sum(9), global_sum(11), proc(2)
   2        format(I6, F12.5, ES16.9, 3(ES14.7), I3)

# ifdef ETALON_CHECK
            do icheck=1,ncheck
              if (iic-1==check_point(icheck)) then
                write(check_line,3,iostat=iocheck) cff,global_sum(7),
     &                                  global_sum( 9),global_sum(11)
   3            format(ES15.9,3(ES14.7))
                if (check_line == etalon_line(icheck)) then
                  write(*,*) 'PASSED_ETALON_CHECK'
                else
                  read(check_line         ,4,iostat=iocheck) (A1(i),
     &                                             P1(i), i=1,nparam)
                  read(etalon_line(icheck),4,iostat=iocheck) (A0(i),
     &                                             P0(i), i=1,nparam)
   4              format(f11.9,1x,i3,3(f10.7,1x,i3))
                  do i=1,nparam
                    A1(i)=A1(i)-A0(i)*10.**float(P0(i)-P1(i))
                  enddo
                  write(check_line,5,iostat=iocheck) (A1(i)
     &                                            ,i=1,nparam)
   5              format(f14.9,3(1x,f13.7))
                  isym=loc(check_line)
                  j=0
                  do i=2,56
                    if (symbol(i)=='.'.and.symbol(i-1)=='0') then
                      symbol(i)=':'
                      j=1
                    else if (j==1 .and. symbol(i)=='0') then
                      symbol(i)='.'
                    else
                      j=0
                    endif
                  enddo
                  write(*,6) ' difference:    ',check_line
   6              format(A16,A56)
                  write(*,*) 'ETALON_CHECK: DIFFERENCE DISCOVERED'
                endif
              endif
            enddo
# endif
          endif
C$OMP END CRITICAL (bio_cr_rgn)
      endif
      end
#else
      subroutine bio_diag_empty
      end
#endif
