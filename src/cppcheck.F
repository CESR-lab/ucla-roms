      program cppcheck

! PURPOSE: Scan all existing CPP-switches in file "cppdefs.opt" and
! automatically generate file "check_switches1.F" which contains
! subroutine "check_switches1". When later this file is compiled and
! executed as a part of ROMS code, it creates log of activated
! CPP-switches.
!
! Algorithm: this program reads line-by-line file "cppdefs.opt" and
! creates catalog of CPP-switches found there.  It does not matter,
! whether switches are in defined or undefined status, either way
! they are put into catalog. For the purpose of this algorithm
! CPP-switch (CPP-macro name) is a word which follows a command
! (reserved word) of C-preprocessor, such as "ifdef", "define" etc.
! Conversely, a word which follows another word which is not a
! CPP-command is not considered as a CPP-macro name.
!
! For the purpose of preceding paragraph "word" means a consecutive
! string of non-blank and non-special characters (i.e. letters, digits
! and underscore '_').
!
! The algorithm works as follows:
!---- --------- ----- -- --------
! 0. reset catalog: arrays of names of CPP-switches and their sizes;
! 1. read line from the file; ignore all lines which do not have #;
! 2. find non-trivial length of that line;
! 3. set all symbols within C-style comments to blank (' ');
! 4. set all special characters to blank;
! 5. identify words within the modified line;
! 6. identify words which are CPP commands and CPP-macro names.
! 7. for each CPP-macro name encountered in (6) check whether it
!     already listed in the catalog, and if not, place it there.
! 8. Once catalog is complete, generate code for checkdefs.
!
! Created by Alexander Shchepetkin <alex@atmos.ucla.edu> on May 2000.

      implicit none
      integer, parameter :: max_str_size=80, max_switches=1024,
     &                      iin=11, iout=12, nwords=16

      character(len=max_str_size) string
      character(len=64) fname(4), switch(max_switches)
      integer nfiles, nswtchs,   nexample,nstdpth, fswtch(8),
     &        ifile,  ierr, i,k, is, istr(nwords), size(max_switches),
     &        count,  lstr,ln,n, ie, iend(nwords), line(max_switches)
      logical word, new, status, state(max_switches),
     &        example,  comment, command, macro(nwords)

! Most modern Fortran compilers recognize "iargc" as an intrinsic
! function, so it does not need to be declared. IBM xlf95 is a notable
! exclusion. So do 7.x and earlier versions of Intel IFC compiler, but
! starting with version 8.x IFORT recognizes it as intrinsic. Portland
! Group PGF compiler also needs it to be declared.

#ifdef INT_IARGC
      integer iargc
#endif
      example=.true. ; nswtchs=0     ! Initialize catalog: reset/blank
      nexample=0     ; fswtch(1)=1   ! out all counters, arraysof macro
      do i=1,max_switches            ! names, and name sizes.
        size(i)=0
        switch(i)='                                '
      enddo     !  12345678901234567890123456789012

      nfiles=iargc()
      if (nfiles==0) then
        write(*,'(/1x,2A/24x,A)') '### ERROR: cppcheck :: No file ',
     &                  'name specified. At least one command-line',
     &                     'argument of must be given to cppcheck.'
        stop
      endif

       write(*,'(/1x,2A/)') 'This is cppcheck: Creating new ',
     &                      'version of "check_switches1.F".'
      do ifile=1,nfiles
        call getarg(ifile,fname(ifile))
        write(*,'(/1x,2A/)') 'Processing ', fname(ifile)
        open(unit=iin,file=fname(ifile),status='old',form='formatted')
        ierr=0
        count=0                 ! Reset error flag, line counter and,
        nstdpth=0               ! nesting depth counter of CPP #ifdefs.
        do while (ierr==0)      ! Blank-out "string", then read in a
          count=count+1         ! line from input file, ignore all
          do i=1,max_str_size   ! lines, which do not start with #,
           string(i:i)=' '      ! otherwise find length of the string
          enddo
          read(iin,'(A)',iostat=ierr) string
          if (string(1:1) == '#') then
            lstr=max_str_size
            do while(string(lstr:lstr) == ' ' .and. lstr > 1)
              lstr=lstr-1
            enddo               ! (position right-most non-blank)

! Suppress C-style comments and special characters.

            n=0     ! <-- reset counter of comments within the string.
            comment=.false.
            i=1
            do while(i < lstr)
              if (.not.comment .and. string(i:i+1)=='/*') then
                comment=.true. ; n=n+1 ; istr(n)=i
              elseif (comment  .and. string(i:i+1)=='*/') then
                comment=.false.        ; iend(n)=i+1
              endif
              i=i+1
            enddo

            if (comment) then             ! If string ends as an open
              lstr=istr(n)-1              ! C-style comment, restrict
              n=n-1                       ! lstr and disregard all
            endif                         ! symbols on the right.
            do k=1,n
              do i=istr(k),iend(k)
                string(i:i)=' '           ! blank out C-style comments
              enddo
            enddo                         ! Suppress special characters
            do i=1,lstr
c*            if (string(i:i) == '(' .or. string(i:i) == ')' .or.
c*   &            string(i:i) == '&' .or. string(i:i) == '|'  .or.
c*   &            string(i:i) == '!' .or. ichar(string(i:i)) == 9)
c*   &          string(i:i)=' '
                                          ! Character 9 is TaB symbol.
              k=ichar(string(i:i))
              if ( k < 48 .or. (k > 57 .and. k < 65) .or.
     &            (k > 90 .and. k < 95) .or. k == 96 .or.
     &                         k > 122 ) string(i:i)=' '
            enddo

! Identify words within the string, find positions of starting and
! ending characters of each word. Since all special characters have
! been removed, at this point word is a sequence of non-blank
! characters.

            n=0       ! <-- reset counter of words within the string.
            word=.false.
            i=1
            do while(i < lstr)
              i=i+1
              if (string(i:i) /= ' ' .and. .not.word) then
                word=.true. ;  n=n+1 ; istr(n)=i
              elseif (string(i:i) == ' ' .and.  word) then
                word=.false.         ; iend(n)=i-1
              endif
            enddo
            if (word) iend(n)=i

c**         write(*,'(/,I4,I4,/)') count, n   ! Print out words.
c**         do k=1,n
c**           write(*,'(10x,80A1)') (string(i:i), i=istr(k),iend(k))
c**         enddo
                                        ! Identify CPP-commands (i.e.
            command=.false.             ! reserved words) and CPP-macro
            do k=1,n                    ! names among the words of the
              macro(k)=.false.          ! line. Cancel example switch
              is=istr(k)                ! when encounter the first
              ie=iend(k)                ! conditional CPP-command.
              ln=ie-is+1
              status=.false.

              if (ln == 6 .and. string(is:ie) == 'define') then
                command=.true. ; status=.true.
              elseif (ln == 5 .and. string(is:ie) == 'undef') then
                command=.true.
              elseif (ln == 2 .and. string(is:ie) == 'if') then
                command=.true. ; example=.false.
                nstdpth=nstdpth+1
              elseif (ln == 5 .and. string(is:ie) == 'ifdef') then
                command=.true. ; example=.false.
                nstdpth=nstdpth+1
              elseif (ln == 7 .and. string(is:ie) == 'defined') then
                command=.true. ; example=.false.
              elseif (ln == 4 .and. string(is:ie) == 'elif') then
                command=.true. ; example=.false.
              elseif (ln == 4 .and. string(is:ie) == 'else') then

              elseif (ln == 5 .and. string(is:ie) == 'endif') then
                nstdpth=nstdpth-1
              elseif (ln == 7 .and. string(is:ie) == 'include') then

              elseif (command) then
                command=.false. ;  macro(k)=.true.
c**           elseif (string(istr(1):iend(1)) /= 'include') then
c**             write (*,'(6x,A,1x,A,1x,I4,A1/8x,A)')    '### ERROR: ',
c**  &        'Unknown CPP-command on line', count, ':', string(is:ie)
              endif
            enddo

c**         write(*,'(/2I4/)') count,n ! Print out CPP-macro names.
c**         do k=1,n
c**           if (macro(k)) then
c**             write(*,'(10x,80A1)') (string(i:i),i=istr(k),iend(k))
c**           endif
c**         enddo
                                          ! Scan catalog of previously
            do k=1,n                      ! discovered switches to find
              if (macro(k)) then          ! match with CPP-macro names
                is=istr(k)                ! found in the present line.
                ie=iend(k)                ! If no match is found, add
                ln=ie-is+1
                if (ifile > 1 .and. nstdpth > 0) then
                  new=.false.
                else                      ! the new switch to catalog.
                  new=.true.
                  do i=fswtch(ifile),nswtchs
                    if (ln == size(i)) then
                     if (string(is:ie) == switch(i)(1:ln)) new=.false.
                    endif
                  enddo                   ! Note that CPP-switches
                endif                     ! found prior to the first
                if (new) then             ! conditional CPP-command
                  nswtchs=nswtchs+1       ! correspond to predefined
                  size(nswtchs)=ln
                  switch(nswtchs)(1:ln)=string(is:ie)
                  state(nswtchs)=status
                  line(nswtchs)=count
                  if (example) nexample=nexample+1
                endif
              endif                       ! configurations (these
            enddo                         ! switches are considered
          endif  !<-- string(1:1) == '#'  ! to be main configuration
        enddo  !<-- while                 ! switches, so only one of
        close(unit=iin)                   ! the is allowed to be
        fswtch(ifile+1)=nswtchs+1         ! defined.)
      enddo  !<-- ifile

      write(*,'(/A,I4)')  'nswtchs =', nswtchs  ! Print out catalog.
      do ifile=1,nfiles
        write(*,'(/A)') fname(ifile)
        do i=fswtch(ifile),nswtchs
          ln=size(i)
          write(*,'(10x,I4,I4,2x,A)') line(i), ln, switch(i)(1:ln)
        enddo
      enddo

! Generate CPP-checking subroutine.

      open (unit=iout,file='check_switches1.F',form='formatted')
      write(iout,'(A/)')  '#include "cppdefs.opt"'
      write(iout,'(/6x,A/)') 'subroutine check_switches1(ierr)'
      write(iout,'(4(A,1x,A/),A,14x,A/A,1x,A,3x,A/A,14x,A/A,22x,A)')
     &  '!!!!!! WARNING: THIS IS A MACHINE GENERATED',
     &                                   'CODE, DO NOT EDIT! !!!!!!',
     &  '!!!!!! This file needs to be updated only if',
     &                                    'new CPP-switches  !!!!!!',
     &  '!!!!!! were introduced into "cppdefs.opt".',
     &                                ' NO ACTION IS NEEDED  !!!!!!',
     &  '!!!!!! if changes in "cppdefs.opt" are limited',
     &                                    'to activation or  !!!!!!',
     &  '!!!!!! deactivation of previously known switches.','!!!!!!',
     &  '!!!!!! To refresh this file compile and execute',
     &                                      '"cppcheck.F"', '!!!!!!',
     &  '!!!!!! as an independent program, or use commands','!!!!!!',
     &  '!!!!!! "make checkdefs" or "make depend".',        '!!!!!!'
      write(iout,'(A,20x,I3,1x,A/A,23x,I3,1x,A)')
     &  '!!!!!! Number of Configuration Choices:',nexample, '!!!!!!',
     &  '!!!!!! Total number of CPP-switches:', nswtchs,  '!!!!!!'
      write(iout,'(9(/6x,A))')
     &  'use param', 'use scalars', 'use strings', 'implicit none',
     &  'integer ierr, is,ie, iexample',
     &  'do is=1,max_opt_size', '  cpps(is:is)='' ''',   'enddo',
     &                          'iexample=0'

      do ifile=1,nfiles
        write(iout,'(A/A,1x,A/A)') '!', '!', fname(ifile), '!'

        ln=32
        do while(fname(ifile)(ln:ln) == ' ')
          ln=ln-1
        enddo
c**     write(iout,'(6x,4A)') 'mpi_master_only write(*,''(1x,A)'') ',
c**  &                              '''<', fname(ifile)(1:ln), '>'''

        if (ifile == 1) then
          if (ln > 8) then
             write(iout,'(6x,A,I2)')   'is=1 ; ie=is+', ln+1
          else
             write(iout,'(6x,A,I1)')   'is=1 ; ie=is+', ln+1
          endif
        else
          if (ln > 8) then
             write(iout,'(6x,A,I2)') 'is=ie+2 ; ie=is+', ln+1
          else
             write(iout,'(6x,A,I1)') 'is=ie+2 ; ie=is+', ln+1
          endif
        endif

        write(iout,'(6x,A/6x,3A/6x,A/6x,A)')
     &             'if (ie >= max_opt_size) goto 99',
     &             'cpps(is:ie)=''<', fname(ifile)(1:ln), '>''',
     &             'cpps(ie+1:ie+1)='' '''

        do i=fswtch(ifile),fswtch(ifile+1)-1
          ln=size(i)
          if (ifile == 1) write(iout,'(A,1x,A)') '#ifdef',
     &                                        switch(i)(1:ln)
          if (ifile == 1 .or. state(i)) then
            if (i <= nexample) write(iout,'(6x,A)')
     &                                 'iexample=iexample+1'
            if (ln > 10) then
              write(iout,'(6x,A,I2)') 'is=ie+2 ; ie=is+', ln-1
            else
              write(iout,'(6x,A,I1)') 'is=ie+2 ; ie=is+', ln-1
            endif
            write(iout,'(6x,A/6x,3A/6x,A)')
     &                 'if (ie >= max_opt_size) goto 99',
     &                 'cpps(is:ie)=''', switch(i)(1:ln), '''',
     &                 'cpps(ie+1:ie+1)='' '''
          endif
          if (ifile == 1) write(iout,'(A)')  '#endif'
        enddo
      enddo

      write(iout,'(6x,A/8x,A/5x,A1,1x,A/8x,A/6x,A)')
     &      'if (iexample == 0) then',
     &      'mpi_master_only write(*,''(/1x,A)'')', '&',
     & '''### ERROR in "cppdefs.opt": no configuration is specified.''',
     &      'ierr=ierr+1',  'elseif (iexample > 1) then'
      write(iout,'(8x,A/5x,A1,2x,A/5x,A1,13x,A/8x,A/6x,A/6x,A)')
     &      'mpi_master_only write(*,''(/1x,2A/)'')', '&',
     &      '''### ERROR: More than one main configuration '',', '&',
     &      '''switch is defined in "cppdefs.opt".''',   'ierr=ierr+1',
     &                                             'endif', 'return'
      write(iout,'(2x,A/5x,A1,2x,A,1x,A/5x,A1,2x,A,1x,A)')
     &      '99  mpi_master_only write(*,''(/1x,2A/14x,A)'')',  '&',
     &      '''### ERROR: checkdefs :: Insufficient size of string',
     &      '"cpps" '',', '&', '''in file "strings".'',',
     &                    '''Increase the size it and recompile.'''
      write(iout,'(6x,A,2(/6x,A))') 'ierr=ierr+1', 'return', 'end'
      close(unit=iout)
      end
