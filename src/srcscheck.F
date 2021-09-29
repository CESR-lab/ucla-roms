      program srcscheck

! PURPOSE: Read list of source code files SRCS (also KRNSRC, IOSRC)
! from Makefile and create file "check_srcs.F" which contains
! subroutine check_srcs to put the list of source codes into string
! "scrs". The string is subsequently written as a global attribute
! SCRS into model history/restart netCDF files to automatically
! document  the model configuration (that is composition of names of
! source files at the moment when model was compiled).

      use strings
      implicit none

      integer, parameter :: iin=11, iout=12, max_string=128
      character tab*1, backslash*1, string*128
      integer lstr, i, is,ie, indx, iocheck
      logical end_of_file, sources

      tab=char(9)
      backslash=char(92)
      sources=.false.
      end_of_file=.false.

      write(*,'(/1x,A,1x,A/)') 'This is SRCSCHECK:',
     &          'Creating new version of check_srcs.F.'

      open (unit=iout,file='check_srcs.F',form='formatted')
      write(iout,'(A/, /6x,A/, /3(A,1x,A/),A,29x,A/)')
     &     '#include "cppdefs.opt"',    'subroutine check_srcs',
     &     '!!!!!! WARNING: THIS IS A MACHINE GENERATED',
     &                              'CODE, DO NOT EDIT! !!!!!!',
     &     '!!!!!! This file needs to be updated only if',
     &                               'the new files     !!!!!!',
     &     '!!!!!! were introduced into or deleted from',
     &                              'the list of source !!!!!!',
     &     '!!!!!! codes SRCS in the Makefile.',       '!!!!!!'
      write(iout,'(4(/6x,A)/8x,A/6x,A)')
     &     'use strings', 'implicit none', 'integer i',
     &     'do i=1,max_opt_size',  'srcs(i:i)='' ''',   'enddo'

      indx=1
      open(unit=iin,file='Makefile')
   1  string=' '
       read(iin,'(A)',iostat=iocheck,end=2) string
        goto 3
   2    end_of_file=.true.
   3    lstr=max_string
   4    if (string(lstr:lstr)==' ') then
          lstr=lstr-1
          if (lstr>0) goto 4
        endif
        if (lstr==0 .and. .not.end_of_file) goto 1
        do i=1,lstr
          if (string(i:i)==tab .or.
     &        string(i:i)=='=') string(i:i)=' '
        enddo
        is=1
   5    if (string(is:is)==' ') then
          is=is+1
          if (is<lstr) goto 5
        endif

        if (string(is:is+5)=='KRNSRC') then
          sources=.true.
        elseif (string(is:is+4)=='IOSRC') then
          sources=.true.
        elseif (string(is:is+3)=='SRCS') then
          sources=.true.
        endif

        if (sources) then
          if (string(lstr:lstr)==backslash) then
            lstr=lstr-1
            if (lstr==0 .and. .not.end_of_file) goto 1
          else
            sources=.false.
          endif

   6      if (string(is:is)==' ') then
            is=is+1
            if (is<lstr) goto 6
          endif
          ie=is
   7      if (string(ie+1:ie+1)/=' ') then
            ie=ie+1
            if (ie<lstr) goto 7
          endif

          if (indx==1) then
            write(iout,'(6x,A,I4,A1,I4,3A)') 'srcs(', indx, ':',
     &                   indx+ie-is, ')=''', string(is:ie), ''''
            indx=indx+ie-is
          else
            write(iout,'(6x,A,I4,A1,I4,3A)') 'srcs(', indx, ':',
     &                indx+ie-is+1, ')='' ', string(is:ie), ''''
            indx=indx+ie-is+1
          endif
          if (indx>max_opt_size) then
            write(iout,'(A/2(A,1x,A/),A)') '*****', '*****   ERROR:',
     &           'parameter max_opt_size in file "strings" is not',
     &           '*****          sufficient to accomodate SCRS list',
     &           'from Makefile.', '*****'
            goto 99
          endif
          indx=indx+1
          is=ie+1
          if (is<=lstr) goto 6
        endif
       if (.not.end_of_file) goto 1
      close(unit=iin)
  99  write(iout,'(6x,A/6x,A)') 'return', 'end'
      close(iout)
      stop
      end

