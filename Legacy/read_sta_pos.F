#include "cppdefs.h"
#ifdef STATIONS

      subroutine read_sta_pos (fname, ierr)

! Read i,j indices of station points from an input text file "fname".
! The file is expected to be organized as follows:  the first line is
! an integer number "nstation" the total number of station points; the
! second and the following lines contain i,j indices (two integers) of
! station points (one pair per line, hence nstation+1 total lines in
! the file).  Returned status "ierr" is assumed be initialized from
! the outside (normally to zero), so the code below just adds error
! counts to it if errors occur.

      implicit none
      character(len=*) fname
      integer ierr, iost, k, lstr, lenstr
      integer, parameter :: inpst=40
# include "param.h"
# include "ncvars.h"
      lstr=lenstr(fname)
      open(unit=inpst, file=fname(1:lstr), form='formatted',
     &                                 status='old', err=99)

      read(inpst,*,iostat=iost) nstation
      if (iost==0) then
        if (nstation<=NS) then
          mpi_master_only write(*,'(/4x,2A,I5)')    'Number of ',
     &        'stations to write into stations file:', nstation

          do k=1,nstation
            read(inpst,*,iostat=iost) ispos(k), jspos(k)
            if (iost==0) then
              if (ispos(k)<1 .or. ispos(k)>LLm  .or.
     &            jspos(k)<1 .or. jspos(k)>MMm) then
                mpi_master_only write(*,'(8x,A,I4,2x,A,2I6,1x,A)')
     &                '### ERROR: station #', k, 'i,j =', ispos(k),
     &                       jspos(k), 'is outside grid boundaries'
                ierr=ierr+1
              else
                mpi_master_only write(*,'(19x,A,I4,2x,A,2I6)')
     &             'station #', k, 'i,j =',  ispos(k), jspos(k)
              endif
            else
              write(*,'(1x,2A,I4,1x,3A)') '### ERROR: read_sta_pos ',
     &                      ':: Cannot read i,j-indices for station',
     &                       k, 'from ''', fname(1:lstr), ''' file.'
              ierr=ierr+1
            endif
          enddo
        else
          write(*,'(/1x,2A,I4,1x,2A,I4,A/)')         '### ERROR: ',
     &             'read_sta_pos :: parameter NS =', NS, 'is too ',
     &                 'small. Should be at least ', nstation, '.'
          ierr=ierr+1
        endif
      else
        write(*,'(/1x,4A)')    '### ERROR: read_sta_pos :: Cannot ',
     &        'read number of station from ''', fname(1:lstr), '''.'
        ierr=ierr+1
      endif
      close(inpst)
      return
  99  write(*,'(/1x,4A)')    '### ERROR: read_sta_pos :: Cannot ',
     &       'open station position file ''', fname(1:lstr), '''.'

      ierr=ierr+1
      end
#else
      subroutine read_sta_pos_empty
      end
#endif /* STATIONS */
