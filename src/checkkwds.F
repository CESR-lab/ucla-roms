      program checkkwds

! Reads file "read_inp.F" as an input text file and creates file
! "setup_kwds.F" which contains subroutine setup_kwds which is adjoint
! to subroutine read_inp in the following sense:
!
!  1. CPP-control structure of "setup_kwds.F" is identical to that of
!     "read_inp.F"; and
!
!  2. within this CPP-control structure, for each code statement of
!     "read_inp.F",
!
!         if/elseif (keyword(1:kwlen)=='keyword_string') then
!
!     it creates a code segment in "setup_kwds.F" which records this
!     particular 'keyword_string' into a character array "kwds".
!
! When "setup_kwds.F" is compiled together with the rest of ROMS code
! and is executed, it pre-records all the keywords which should be
! present in the startup file to guarantee proper initialization of the
! model: "setup_kwds" along with two other subroutines, "cancel_kwd"
! and "check_kwds" provides a safety mechanism, to verify that all the
! keyword-blocks required by the particular model configuration are
! present in the startup file, and present only once.
! This works as follows: "setup_kwds" is called in the beginning of
! "read_inp" (before reading the startup file) and crates a character
! str which contains the set of keywords required by the model
! (depending on CPP-switches); then, every time when a keyword-block
! in identified and read from the startup file, call to "cancel_kwd"
! is made to remove (blank out) that keyword from the str.
! After the last reading of the startup file is complete, call to
! "check_kwds" is made to verify that the the whole keyword is blank;
! if not, print an error message about what keywords were not canceled
! and terminate the model run.
!
! Created by Alexander Shchepetkin <alex@atmos.ucla.edu>, June 2000.

      implicit none
      integer, parameter :: iin=11, iout=12, max_str_size=80
      integer  ierr, i, is,ie, lstr
      logical  cpp_cont_line
      character(len=max_str_size) str
      character(len=1) quote, backslash
      quote=char(39)
      backslash=char(92)

! Open output file and create subprogram title and declarations.

      write(*,'(/1x,A,1x,A/)') 'This is CHECKKWDS: Creating',
     &                         'new version of "setup_kwds.F".'
      open(iout,file='setup_kwds.F', form='formatted')
      write(iout,'(A/,/6x,A/3(/A,1x,A)/A,8x,A/A,1x,A/A,39x,A)')
     &  '#include "cppdefs.opt"',    'subroutine setup_kwds (ierr)',
     &  '!!!!!! WARNING: THIS IS A MACHINE GENERATED',
     &                                   'CODE, DO NOT EDIT! !!!!!!',
     &  '!!!!!! This file needs to be updated only if',
     &                                    'new keywords were !!!!!!',
     &  '!!!!!! introduced into "read_inp.F".',
     &                            'To create or refresh this !!!!!!',
     &  '!!!!!! file use compile and execute "checkkwds.F" as an',
     &  '!!!!!!', '!!!!!! independent program, or use',
     &                          'commands "make checkkwds"   !!!!!!',
     &  '!!!!!! or "make depend".',                         '!!!!!!'
      write(iout,'(5(/6x,A), /6x,A /8x,A, 2(/6x,A))')
     &  'use param', 'use scalars', 'use strings',
     &  'implicit none', 'integer ierr, is,ie', ! DevinD edited for param module '#include "param.h"',
     &  'do is=1,max_opt_size',  'kwds(is:is)='' ''',
     &  'enddo',    'is=1'

! Open input file and process it line-by-line.  At first, set "str"
! to all-blank, then read in a new line from input file. Find position
! of the last non-blank character.  Disregard lines which are empty,
! Fortran comment lines, and #include CPP-directives.  Set lstr=0 to
! signal that the line should be ignored.

      open(unit=iin, file='read_inp.F', status='old', form='formatted')
      cpp_cont_line=.false.
      ierr=0
      do while(ierr==0)
        do i=1,max_str_size
          str(i:i)=' '
        enddo
        read(iin,'(A)',iostat=ierr) str
        lstr=max_str_size
        do while (str(lstr:lstr)==' ' .and. lstr>1)
          lstr=lstr-1
        enddo
        if (str(lstr:lstr)==' ' .or. str(1:1)=='!' .or.
     &            str(1:1)=='C' .or. str(1:1)=='c') then
          lstr=0
        elseif (str(1:1)=='#') then
          is=2
          do while (str(is:is)/='i' .and. is<lstr-6)
            is=is+1
          enddo
          if (str(is:is+6)=='include') lstr=0
        endif

! Copy CPP-directives lines (other than  #include rejected above )into
! the target file.  Also set CPP-directive continuation switch to copy
! the next line.

        if (lstr>0) then
          if (str(1:1)=='#' .or. cpp_cont_line) then
            write(iout,'(A)') str(1:lstr)
            if (str(lstr:lstr)==backslash) then
              cpp_cont_line=.true.
            else
              cpp_cont_line=.false.
            endif
          else

! Recognize lines containing << keyword(1:kwlen) == 'keyword_string' >>
! and create the corresponding code to record the actual string inside
! the quotes...

            is=7
            do while (str(is:is)/='k' .and. is<lstr-20)
             is=is+1
            enddo
            if (str(is:is+6)=='keyword') then
              is=is+7
              do while (str(is:is)/='(' .and. is<lstr-13)
                is=is+1
              enddo
              do while (str(is:is)/='1' .and. is<lstr-12)
                is=is+1
              enddo
              do while (str(is:is)/=':' .and. is<lstr-11)
                is=is+1
              enddo
              do while (str(is:is)/='k' .and. is<lstr-10)
                is=is+1
              enddo
              if (str(is:is+4)=='kwlen') then
                is=is+5
                do while (str(is:is)/=')' .and. is<lstr-9)
                  is=is+1
                enddo
                is=is+1
                do while (str(is:is)==' ' .and. is<lstr-4)
                  is=is+1
                enddo
                if (str(is:is+1)=='==' .or. (str(is:is+1)=='.e'
     &                       .and. str(is+2:is+3)=='q.') ) then
                  is=is+2
                  do while (str(is:is)/=quote .and. is<lstr)
                    is=is+1
                  enddo
                  if (str(is:is)==quote) then
                    ie=is+1
                    do while (str(ie:ie)/=quote .and. ie<lstr)
                      ie=ie+1
                    enddo
                    if (str(ie:ie)==quote) then
                      write(iout,'(6x,A,I2/6x,A/6x,A,A/6x,A/6x,A)')
     &                           'ie=is +', ie-is-1,
     &                           'if (ie>=max_opt_size) goto 99',
     &                           'kwds(is:ie)=',     str(is:ie),
     &                           'kwds(ie+1:ie+1)='' ''', 'is=ie+2'
                    endif
                  endif
                endif
              endif
            endif
          endif
        endif
      enddo
      close(unit=iin)                 !--> Finalize the subprogram...

      write(iout,'(6x,A)') 'return'
      write(iout,'(2x,A/5x,A,2x,2A/5x,A,2x,A,2x,A)')
     &     '99  mpi_master_only write(*,''(/1x,2A/25x,A/)'')', '&',
     &          '''### ERROR: setup_kwds :: Insufficient size of ',
     &        'str "kwds" '',', '&','''in file "strings".'',',
     &                   '''Increase the size it and recompile.'''
      write(iout,'(6x,A,2(/6x,A))') 'ierr=ierr+1', 'return', 'end'
      close(iout)
      stop
      end
