      program redefs
!
! A utility program to create global redefinition file "redefs.X"
! for the purpose of creating quasi-independent entities for multiply
! nested grids within the same executable file. The redefinition file
! consists of lines
!                     #define name name_suffix
!
! where 'name' is a name of a function, subroutine, common block, and
! /or fortran 90 module; and  _suffix = _1, _2, _3, etc, dipending on
!
! command line usage:
!                        redefs *.h *.F
! within the directory.

      implicit none
      integer, parameter :: iin=11, iout=12,
     &          max_string_size=72, maxwords=256
      logical not_end_of_file
      character suffix*2, fname*32,  bffr*(max_string_size),
     &  string*(16*max_string_size), small*(16*max_string_size)
      integer kA,kZ,ksa,ksz,k0,k9,ksft, line, istr,iend, length,
     &    nwords,is(maxwords),ie(maxwords), i,j,k,n, last_arg, iargc

! Derive suffix from the current directory name. The last symbol of
! the directory name is expected to be a digit 1, 2, 3, etc... where
! 1 corresponds to the parent grid and >1 to child grids. If the last
! symbol is not a digit, the directory is assumed to contain code for
! the parent grid and suffix is set to '_1'.

      do i=1,max_string_size
        bffr(i:i)=' '
      enddo
      call system ('pwd > pwd.d')
      open (unit=iin, file='pwd.d', form='formatted', status='old')
      read (iin,'(A)') bffr
      close (unit=iin, status='delete')
      i=1
      do while (bffr(i:i)==' ' .and. i<max_string_size)
        i=i+1
      enddo
      j=i
      do while (bffr(j:j)/=' ' .and. j<max_string_size)
        j=j+1
      enddo
      if (bffr(j:j)==' ') j=j-1
      if (bffr(j:j)>'1' .and. bffr(j:j)<='9') then
        suffix(1:1)='_'
        suffix(2:2)=bffr(j:j)
      else
        suffix='_1'
      endif

      kA=ichar('A')           ! Save ASCII code indices of some
      ksa=ichar('a')          ! symbols to use them later in case
      ksft=ksa-kA             ! filding algoritm and for name
      kZ=ichar('Z')           ! verification (it is assumed here
      ksz=ichar('z')          ! that valid fortran names may contain
      k0=ichar('0')           ! only letters, digits and underscore
      k9=ichar('9')           ! symbol).
!
! Main code starts here: read files one-by-lone and examine them for
! the presense of objects to be renamed. Note that auxiliary program
! "read_bffr" is called in delayed mode [i.e. at the beginning it
! copies content of array "bffr" left from reading of the previous
! line (hence preliminary call is needed to initialize) into "string"
! and call "read_bffr" again], which is needed in order to check
! whether the following line is a continuation line, and if is, it
! is atached to "string" (the process continues until a starting line
! is encountered, after which "string" is considered complete, and
! "bffr" is saved for the next time). This in needed to catch
! multiple common block declaration by a single common statement
! spread on several lines.
!
      last_arg=iargc()
      if (last_arg>0) then
        open (unit=iout, file='redefs.X', form='formatted')
        write(iout,'(A,1x,A)') '#define GRID_LEVEL', suffix(2:2)
        do n=1,last_arg
          call getarg(n,fname)
          i=1
          do while (fname(i:i)/=' ')
            i=i+1
          enddo
          i=i-1
          write(iout,'(A,1x,A,1x,A)') '/*', fname(1:i), '*/'
          open(unit=iin, file=fname(1:i), form='formatted',
     &                                 status='old', err=7)
          line=0
          not_end_of_file=.true.
          call read_bffr (fname, iin, line,  max_string_size,
     &                        bffr, istr,iend, not_end_of_file)
          do while (not_end_of_file)
            string(1:iend)=bffr(1:iend)
            length=iend
 1          call read_bffr (fname, iin, line,  max_string_size,
     &                          bffr, istr,iend, not_end_of_file)
            if (istr>0.and.iend>0 .and. bffr(6:6)/=' ') then
              string(length+1:length+1)=' '
              string(length+2:length+iend-istr+2)=bffr(istr:iend)
              length=length+iend-istr+2
              goto 1
            endif
            string(length+1:length+1)=' '
            if (length>0) then
              istr=1
              do while (string(istr:istr)==' '.and. istr<length)
                istr=istr+1
              enddo
!
! The following code segment forms an auxiliary string "small" which
! is a case-folded version of "string" with most of non-letter, non-
! digit, and non-underscore characters replaced with blank. This step
! is needed to allow case-blind identification of names (as words
! within the string), which occur in the filowing code segment.
!
              do i=istr,length
                k=ichar(string(i:i))
                if (k>=kA .and. k<=kZ) then
                  small(i:i)=char(k+ksft)
                elseif (k>=ksa .and. k<=ksz) then
                  small(i:i)=string(i:i)
                elseif (k>=k0 .and. k<=k9) then
                  small(i:i)=string(i:i)
                elseif (string(i:i)=='_') then
                  small(i:i)='_'
                elseif (string(i:i+1)=='/*' .or.
     &                  string(i:i+1)=='*/') then
                  small(i:i+1)='  '
                  string(i:i+1)='  '
                elseif (string(i:i)=='/') then
                  small(i:i)='/'
                else
                  small(i:i)=' '
                  string(i:i)=' '
                endif
              enddo

              nwords=0
              i=istr
              do while (i<length)
                do while (string(i:i)==' ')
                  i=i+1
                enddo
                if (nwords < maxwords) then
                  nwords=nwords+1
                  is(nwords)=i
                  ie(nwords)=0
                  do while (string(i:i)/=' '.and.i<length+1)
                    i=i+1
                  enddo
                  ie(nwords)=i-1
                else
                  write(*,*) 'WARNING: too many words on line',
     &                        line, '''', string(1:length), ''''
                  goto 2
                endif
              enddo
    2         if (ie(nwords)==0) ie(nwords)=length

!
! Identify lines which are declarations of objects of interest,
! find names of that objects (usually the following word, except
! common blocks, where it is a bit more complicated), and write
! out CPP-redefinition command in the output file.
!
              k=0
              do i=1,2
                if (small(is(i):ie(i)) == 'common' .or.
     &              small(is(i):ie(i)) == 'module'   .or.
     &              small(is(i):ie(i)) == 'subroutine' .or.
     &              small(is(i):ie(i)) == 'function') k=i
              enddo
              if (k==1 .and. small(is(1):ie(1))=='common') then
                k=0
                do i=is(2),ie(nwords)
                  if (string(i:i)=='/') then
                    if (k==0) then
                      k=i
                    else
                      write(iout,'(A,1x,A,1x,2A)') '#define',
     &                string(k+1:i-1), string(k+1:i-1), suffix
                      k=0
                    endif
                  endif
                enddo
              elseif ((k>0 .and. small(is(k):ie(k))=='function')
     &   .or. (k==1 .and. small(is(1):ie(1))=='subroutine')) then
                write(iout,'(A,1x,A,1x,2A)') '#define',
     &                         string(is(k+1):ie(k+1)),
     &                         string(is(k+1):ie(k+1)), suffix
              endif
            endif
          enddo
          close (iin)
        enddo
        close (iout)
      else
        write(*,'(/1x,A/)') 'ERROR in redefs: no files to process.'
      endif
      stop
  7   write(*,'(/3A/)') 'ERROR: Cannot open ''', fname(1:i), '''.'
      stop
      end

!
! Auxiliary code:  Read in one line of input file and perform
! preliminary analysis of it: strip out comment lines (of portion
! of the line which is comment), find starting and ending nonblanc
! characters of non-trivial portion of the line; increment line
! counter and set end-of-file flag.
!
      subroutine read_bffr (fname, iin, line, max_string_size,
     &                                    bffr, istr,iend, leof)
      implicit none
      integer iin, line, max_string_size, istr,iend, i, ierr
      character fname*(*),  bffr*(*),  quote*1
      logical leof, lswtch
                               ! Reset everything to all-blanc and
      quote=char(39)           ! zero status, then read in one line.
      iend=0                   ! Check whether it is an F77 style
      do i=1,max_string_size   ! comment, and if so, check further
        bffr(i:i)=' '          ! whether dollar sign is present, so
      enddo                    ! it is a directive or conditionally
      do while (leof .and. iend==0)
        line=line+1
        read(iin,'(A)',iostat=ierr,end=1) bffr
        goto 2
  1     leof=.false.
  2     if (bffr(1:1)=='C' .or. bffr(1:1)=='c' .or.
     &                            bffr(1:1)=='!') then
          istr=0
          do i=2,5
            if (bffr(i:i)=='$') istr=2
          enddo
        elseif (bffr(1:1)=='#') then
          istr=0
        else                   ! compiled line of code. Note that
          istr=1               ! variable "lstr" is used as signal
        endif                  ! flag here. If the line is not a
        if (istr>0) then    ! comment find positions of its first
          lswtch=.false.       ! and last non-blanc symbols.
          do i=istr,max_string_size
            if (bffr(i:i)==quote) then
              lswtch=.not.lswtch
            elseif (bffr(i:i)=='!' .and. .not.lswtch) then
              goto 3
            endif
            if (bffr(i:i)/=' ') iend=i
          enddo
c--#ifdef TRAP_SINGLE_QUOTES
          if (lswtch) then
            i=1
            do while (fname(i:i)/=' ')
              i=i+1
            enddo
            i=i-1
            write(*,'(/6x,2A,I4,1x,3A/)') '#### ERROR: unmatched ',
     &       'quote on line', line, 'in file ''', fname(1:i), '''.'
          endif
c--#endif
   3      do while (iend>1 .and. bffr(iend:iend)==' ')
            iend=iend-1
          enddo
          if (iend==1 .and. bffr(1:1)==' ') iend=0
        endif
      enddo
      if (iend>0) then
        istr=1
        if (bffr(6:6)/=' ') then
          istr=7
        else
          istr=1
        endif
        do while (bffr(istr:istr)==' ' .and. istr<iend)
          istr=istr+1
        enddo
      else
        istr=0
      endif
      do i=iend+1,max_string_size
        bffr(i:i)=' '
      enddo
c**   if (iend>0) write(*,*) bffr(istr:iend)
      end

