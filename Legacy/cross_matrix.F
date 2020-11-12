      program cross_matrix

! Purpose: Analyze file dependencies and build "Make.depend".
!---------
! Usage:
!--------       cross_matrix list_of_files
! for example
!               cross_matrix *.F
! or
!               cross_matrix *.F *.h
! or
!               cross_matrix $(SRCS)            inside a makefile

! Algorithm:
!-----------
! cross_matrix reads all files from the list arguments it is given,
! and searches for the CPP-style #include, Fortran-style INCLUDE, USE,
! and MODULE statements. If any of these found, it looks further and
! interprets the following word (enclosed or not into "...", '...') as
! the name of include file or F90 module. These new files will be read
! and searched for #include inside their own, hence "cross_matrix"
! investigates chained include statements, for example, if file1
! contains statement #include "file2", and file2 contains statement
! #include "file3", therefore file3 will be included into list of
! dependencies of "file1".

!--> Optionally (activated by CPP switch EXCLUDE), files which are
! not in the list of arguments AND which are not found in the present
! working directory can be excluded from the list of dependencies.
! A warning message will be printed about each file excluded from the
! analysis. However, WARNING: This practice is unsafe because, unlike
! the old IRIX SGI CPP, most Linux C-preprocessors are known to simply
! skip the #include command with little warning (or without and warning
! at all) if the file cannot be found by regardless of the reason (e.g.
! its name is misspelled in the source code being analyzed).  If the
! include file contains variable declarations (which is in most, but
! not all the cases) its absence triggers subsequent compilation error,
! so the mishap is discovered.  However, the file which contains purely
! executable code, the miss goes undetected resulting in wrong code
! without any visible signal during the compilation.

! The alternative CPP-switch EXCLUDE_SELECTED allows explicit
! specification of which files to be ignored in the dependency list.
! This option is safe against misspelled files.

! CPP-switch TWO_STAGE_UPCSF makes two-line dependency for .F-->.o
! corresponding to separate .o:.f and .f:.F dependencies, where the
! first one, .o:.f, contains all .mod and F77-style includes, but not
! CPP includes, while, conversely, .f:.F contains only dependencies
! from CPP includes.  This way each stage has correct minimal set of
! dependencies sufficient to execute it separately. Otherwise a single
! line containing all the dependencies from both stages is generated.

! CPP-switch SKIP_TRIVIAL suppresses printing dependency line where
! .o depends only on .f with same name and does not require any .mod or
! include files.   This also applies to .f:.F in the case of two-stage
! compilation.

! CPP-switch WARN_NONSTANDARD makes inspects the source codes for
! symbols beyond the standard Fortran character set and issues warning
! for each one found.

! Option TRAP_DOUBLE_SLASH inspects the source codes and prints warning
! about every double slash found.  This is useful because Fortran uses
! double slash as concatenation operator for character strings, however
! double slash is also known as comment under C++ syntax rules. Because
! of this, if a Fortran code is expected to go through CPP, CPP will
! remove all the symbols right from double slash (including double
! slash itself). This may result in wrong/damaged Fortran code without
! any compilation error or even warning.  The usual remedy against this
! is to put one blank space between slashes when double slash is needed
! as Fortran concatenation operator, however finding all instances of
! double slashes and making sure that all are properly separated may
! be elaborate, so this program automates the process.

! Neither WARN_NONSTANDARD, nor TRAP_DOUBLE_SLASH make any changes in
! any of the source codes being processed.

! Created and maintained by Alexander Shchepetkin old_galaxy@yahoo.com.
! All issues should be reported there.

c--#define EXCLUDE
#define EXCLUDE_SELECTED
#define TWO_STAGE_UPCSF
#define SKIP_TRIVIAL

#define WARN_NONSTANDARD
#define TRAP_DOUBLE_SLASH
c--#define VERBOSE

      implicit none
      integer, parameter :: max_names=561,  max_modules=128,  iin=11,
     &                max_string_size=128,  max_name_size=64, iout=12,
     &                                      max_width=80,     test=13
      character(len=max_string_size) string
      character(len=max_name_size) testword, fname(max_names),
     &                                       mod_name(max_modules)
      character(len=1) backslash, single_quote, double_quote,
     &                                       mtrx(max_names,max_names)
      integer lfnm(max_names), mod_src(max_names), mod_indx(max_names),
     &        nmax, last_arg, n_disc, nmods, itype, size, is, ie, ierr,
     &        line, lfort_all, cpps_all, dirs_all, comm_all, empty_all,
     &                                                        i,j,k,n
      logical new_name, not_end_of_file

! Most modern Fortran compilers recognize "iargc" as an intrinsic
! function, so it does not need to be declared. IBM xlf95 is a notable
! exclusion. So do 7.x and earlier versions of Intel IFC compiler, but
! starting with version 8.x IFORT recognizes it as intrinsic. Portland
! Group PGF compiler also needs it to be declared.

#if defined XLF || defined PGI
      integer iargc
#endif
#ifdef EXCLUDE
      integer, parameter :: exc_names=16
      integer nexc, exc_size(exc_names)
      character(len=32) exc_name(exc_names)
#endif
#ifdef EXCLUDE_SELECTED
      integer, parameter :: exc_names=2     ! Specify list of files
      integer exc_size(exc_names)           ! to be excluded from the
      character(len=16) exc_name(exc_names) ! dependency list.
      exc_name(1)='mpif.h'                  ! Typically this should be
      exc_size(1)=6                         ! be restricted to the 3rd-
      exc_name(2)='netcdf.inc'              ! party files which do not
      exc_size(2)=10                        ! not evolve with the code.
#endif
      single_quote=char(39)
      double_quote=char(34)
      backslash=char(92)
                             ! Reset everything...
      is=0 ; ie=0            !------ -------------
      lfort_all=0            ! Fortran code lines in all files together
      cpps_all=0             ! CPP-preprocessor directive lines
      dirs_all=0             ! compiler directive lines (Open MP, etc)
      comm_all=0             ! commentaries
      empty_all=0            ! empty lines

      n_disc=0               ! <-- counter of newly discovered files
#ifdef EXCLUDE
      nexc=0                 ! <-- counter of excluded files
#endif
      do j=1,max_names
        lfnm(j)=0            ! <-- array of sizes of filenames
        do i=1,max_name_size
          fname(j)(i:i)=' '  ! <-- character array of filenames
        enddo
        do i=1,max_names
          mtrx(i,j)=' '      ! <-- matrix of dependencies such that
        enddo                !         mtrx(i,j)==non-blank
      enddo                  !     means that file j depends on i.

      last_arg=iargc()
      if (last_arg==0) then
        write(*,'(/1x,A/)') '### ERROR: No files to process.'
        stop
      elseif (last_arg > max_names) then
        goto 98                               !--> ERROR
      endif

! The following block creates the list of given names and their sizes
! without actually opening and reading the files.

#ifdef VERBOSE
      write(*,'(1x,A)') 'Initial list:'
#endif
      do n=1,last_arg
        call getarg(n,fname(n))
        i=1
        do while (fname(n)(i:i)==' ' .and. i<max_name_size)
          i=i+1
        enddo
        j=i
        k=j+1
        do while (fname(n)(k:k)/=' ' .and. k<max_name_size)
          j=k
          k=j+1
        enddo
        if (i==1) then
          lfnm(n)=j
        else
          lfnm(n)=j-i+1
          fname(n)(1:j-i+1)=fname(n)(i:j)
        endif
#ifdef VERBOSE
        write(*,'(4I3,1x,3A)') n,i,j, lfnm(n), 'fname = ''',
     &                                   fname(n)(i:j), '''.'
#endif
      enddo
#ifdef VERBOSE
      write(*,*)
#endif
      nmax=last_arg
      nmods=0

! Create Make.depend file and write the upper portion of its
! header (all lines written here are comments starting with #).

      write(*,'(/1x,2A/)') 'This is cross_matrix: Generating ',
     &                       'a new version of "Make.depend".'
      open(unit=iout, file='Make.depend', form='formatted')
      write(iout,'(3(A,1x,A/),A/A/A/A/A/A)')
     &        '# Make.depend: list of dependencies generated by',
     &       'cross_matrix.', '# !!! WARNING: THIS IS A MACHINE',
     &         'GENERATED FILE: DO NOT EDIT !!!', '# To create',
     &         'or update this file use commands:', '#',
     &         '#        cross_matrix *.F',         '# or',
     &         '#        cross_matrix *.F *.h',     '# or',
     &         '#        make depend'


!      *****    *********    ******   *******    *********
!    ***   ***  *  ***  *   *   ***   ***   ***  *  ***  *
!    ***           ***     **   ***   ***   ***     ***
!      *****       ***    ***   ***   ***   **      ***
!          ***     ***    *********   ******        ***
!    ***   ***     ***    ***   ***   ***  **       ***
!      *****       ***    ***   ***   ***   ***     ***

! Reading and processing the input files starts here.

      n=0
  1    n=n+1     !<-- file number from list of args + newly discovered
        line=0     !<-- line counter within the file "n"
        not_end_of_file=.true.
#ifdef VERBOSE
        write(*,'(4x,3A)') 'Analyzing ''', fname(n)(1:lfnm(n)), '''...'
#endif
        open(unit=iin, file=fname(n)(1:lfnm(n)), form='formatted',
     &                                       status='old', err=7)
  2     string(1:1)=' '
         read(iin,'(A)',iostat=ierr,end=3) string
         line=line+1
         goto 4
  3      not_end_of_file=.false.

! Reset "itype" and determine position of the staring and ending
! position non-blank characters "is:ie" for current string. These
! will remain unchanged during the rest of processing of the line.

  4      if (ierr==0 ) then
           itype=0
           ie=max_string_size
           do while (ie>1 .and. string(ie:ie)==' ')
             ie=ie-1
           enddo
           if (string(ie:ie)==' ') ie=ie-1
           is=1
           do while (string(is:is)==' ' .and. is<ie)
             is=is+1
           enddo

! The following segment discriminates lines into several categories
! (empty lines, F77- and F90-style comments, CPP-lines starting with
! #, compiler directives, and regular Fortran statements). Once it is
! determined that the line is not empty, it checks whether the first
! word is "include" (with or without # in front), "use", or "module".
! If so, it sets flag "itype" to non-zero value to signal that the
! line needs further processing. Note that, unlike CPP-style #include,
! Fortran include/use/module can appear both in upper- and lower-case
! letters, so copy the first word of string into "testword", while
! performing folding to uppercase.  Note that each if-branch below
! increments exactly one line counter for the specific line type
! (the counters are used for reporting statistics only and do not
! affect computation of dependencies).

           if (ie==0) then
             empty_all=empty_all+1

           elseif (string(1:1)=='!' .or. string(1:1)=='C' .or.
     &             string(1:1)=='*' .or. string(1:1)=='c') then
             if (ie==1) then
               empty_all=empty_all+1
             elseif (string(2:2)=='$' .or. string(3:3)=='$' .or.
     &               string(4:4)=='$' .or. string(5:5)=='$') then
               dirs_all=dirs_all+1
             else
               comm_all=comm_all+1
             endif

           elseif (string(is:is)=='!') then
             if (is==ie) then
               empty_all=empty_all+1
             else
               comm_all=comm_all+1
             endif

           elseif (string(1:1)=='#') then
             cpps_all=cpps_all+1
             i=2
             do while (string(i:i)==' ' .and. i<ie)
               i=i+1
             enddo
             j=i
             do while (string(j:j)/=' ' .and. j<ie)
               j=j+1
             enddo
             if (string(j:j)==' ') j=j-1
             if (j-i==6 .and. string(i:j)=='include') then
               itype=1
             endif

           else
             lfort_all=lfort_all+1
             i=is
             j=is
             do while (string(j:j)/=' ' .and. j<ie)
               j=j+1
             enddo
             if (string(j:j)==' ') j=j-1

             size=j-i+1
             if (size<8) then
               testword=string(i:j)
               do k=1,size
                 if (testword(k:k)>='A' .and. testword(k:k)<='Z')
     &                     testword(k:k)=char(ichar(testword(k:k))+32)
               enddo
               if (size==7 .and. testword(1:size)=='include') then
                 itype=2
               elseif (size==3  .and. testword(1:size)=='use') then
                 itype=3
               elseif (size==6.and.testword(1:size)=='module') then
                 itype=4
               endif
             endif
           endif

! Once it is determined that the line begins with "[#]include", "use",
! or "module" statement, find bounds and read the second word, which is
! either name of include file or module. Copy it into "testword". Strip
! out double/single quotes, and append with suffix ".mod" if necessary.

           if (itype>0) then  !--> find the second word
             i=j+1
             do while (string(i:i)==' ' .and. i<ie)
               i=i+1
             enddo
             j=i
             do while (string(j:j)/=' ' .and. j<ie)
               j=j+1
             enddo
             if (string(j:j)==' ') j=j-1
             if (string(j:j)==',') j=j-1  !<- use mod_name, only smth

             if ( (string(i:i)=='<' .and. string(j:j)=='>') .or.
     &           ( string(i:i)==single_quote .and.
     &                             string(j:j)==single_quote ) .or.
     &           ( string(i:i)==double_quote .and.
     &                           string(j:j)==double_quote ) ) then
               i=i+1
               j=j-1                     ! Save the name of include
             endif                       ! file or module as "testword"
             size=j-i+1                  ! and, if module, fold all the
                                         ! letters into lowercase and
             if (size>0) then         ! append suffix ".mod".
               testword=string(i:j)
               if (itype==3  .or. itype==4) then
                 do k=1,size
                   if (testword(k:k)>='A' .and. testword(k:k)<='Z')
     &                      testword(k:k)=char(ichar(testword(k:k))+32)
                 enddo
                 testword(size+1:size+4)='.mod'
                 size=size+4
               endif

! At this moment "testword" contains name of an include file or module.
! Check whether this name already appeared in the list of files, and if
! not add it into the list. In both cases, "new" or "old", set element
! of dependency matrix.  Also note that F90-style module statement,
! itype=4, results in REVERSE dependency, since it is the resultant
! .mod file depends on the original source, and not the other way
! around, like in the case of include or use.

               new_name=.true.
#ifdef EXCLUDE_SELECTED
               do k=1,exc_names
                 if (size == exc_size(k)) then
                   if (exc_name(k)(1:exc_size(k)) ==
     &                         testword(1:size) ) then
                     new_name=.false.
# ifdef VERBOSE
                     write(*,'(8x,3A)')   'Excluded selected file ''',
     &               testword(1:size), ''' from the dependency list.'
# endif
                   endif
                 endif
               enddo
#endif
               if (new_name) then
                 do k=1,nmax
                   if (size == lfnm(k)) then
                     if (testword(1:size) == fname(k)(1:size)) then
                       new_name=.false.
                       if (itype==1) then
                         mtrx(k,n)='i'
#ifdef VERBOSE
                         write(*,'(8x,4A,I3,1x,A,I3,A)')     'Added ',
     &                  'previously found CPP-style #include file ''',
     &                                      fname(nmax)(1:lfnm(nmax)),
     &                                  ''' to mtrx(', k, ',', n, ')'
#endif
                       elseif (itype==2) then
                         mtrx(k,n)='j'
#ifdef VERBOSE
                         write(*,'(8x,4A,I3,1x,A,I3,A)')     'Added ',
     &                   'previously found F77-style include file ''',
     &                                      fname(nmax)(1:lfnm(nmax)),
     &                                  ''' to mtrx(', k, ',', n, ')'
#endif
                       elseif (itype==3) then
                         mtrx(k,n)='u'    !<--file "n" depends on "k"
#ifdef VERBOSE
                         write(*,'(8x,4A,I3,1x,A,I3,A)')     'Added ',
     &                     'previously found F90-style use module ''',
     &                                      fname(nmax)(1:lfnm(nmax)),
     &                                  ''' to mtrx(', k, ',', n, ')'
#endif
                       elseif (itype==4) then
                         mtrx(n,k)='M'    !<-- file "k" depends on "n"
                         nmods=nmods+1
                         mod_name(nmods)=testword(1:size)
                         mod_src(nmods)=n  !<-- save who creates module
                         mod_indx(nmods)=k !<-- index within  catalog
#ifdef VERBOSE
                         write(*,'(8x,4A,I3,3(A,I3))')    'Added ',
     &                         'previously found module definition ''',
     &                       fname(nmax)(1:lfnm(nmax)),
     &                                ''' to mtrx(', n, ',', k,
     &                         '), recorded mod_src(', nmods, ') =', n
#endif
                       endif
                     endif
                   endif
                 enddo
                 if (new_name) then
                   n_disc=n_disc+1
#ifdef EXCLUDE
                   if (itype==1 .or. itype==2) then
                     open(unit=test, file=testword(1:size),
     &                                status='old', err=5)
                     close(unit=test)
                   endif
#endif
                   nmax=nmax+1
                   if (nmax>max_names) goto 98  !--> ERROR
                   lfnm(nmax)=size
                   fname(nmax)(1:size)=testword(1:size)
                   if (itype==1) then
                     mtrx(nmax,n)='I'
#ifdef VERBOSE
                     write(*,'(8x,4A,I3,1x,A,I3)')    'Found new CPP-',
     &              'style #include file ''',fname(nmax)(1:lfnm(nmax)),
     &                                   ''', nmax =', nmax, 'n =', n
#endif
                   elseif (itype==2) then
                     mtrx(nmax,n)='J'
#ifdef VERBOSE
                     write(*,'(8x,4A,I3,1x,A,I3)')    'Found new F77-',
     &              'style include file ''', fname(nmax)(1:lfnm(nmax)),
     &                                   ''', nmax =', nmax, 'n =', n
#endif

                   elseif (itype==3) then
                     mtrx(nmax,n)='U'  !<--file "n" depends on "k"
#ifdef VERBOSE
                     write(*,'(8x,4A,I3,1x,A,I3)') 'Found new F90-',
     &               'style use module ''',  fname(nmax)(1:lfnm(nmax)),
     &                                   ''', nmax =', nmax, 'n =', n
#endif
                   elseif (itype==4) then
                     mtrx(n,nmax)='M'   !<--  "nmax" depends on "n"
                     nmods=nmods+1
                     mod_name(nmods)=testword(1:size)
                     mod_src(nmods)=n     !<-- save who creates module
                     mod_indx(nmods)=nmax !<-- index within the catalog
#ifdef VERBOSE
                     write(*,'(8x,4A,I3,2(1x,A,I3))')     'Found new ',
     &               'module definition ''', fname(nmax)(1:lfnm(nmax)),
     &                  ''', nmax =', nmax, 'n =', n, 'nmods =', nmods
#endif
                   endif
#ifdef EXCLUDE
                   goto 6
   5               do k=1,nexc
                     if (size == exc_size(k)) then
                      if (testword(1:size) == exc_name(k)(1:size))
     &                                              new_name=.false.
                     endif
                   enddo
                   if (new_name) then
                     nexc=nexc+1
                     if (nexc>exc_names) goto 99 !--> ERROR
                     exc_size(nexc)=size
                     exc_name(nexc)=testword(1:size)
# ifdef VERBOSE
                     write(*,'(I4,1x,3A)')   nexc, 'Exclude ''',
     &                   exc_name(nexc)(1:exc_size(nexc)), '''.'
# endif
                     write(iout,'(A/A,2x,A/2A)') '#',
     &                              '# WARNING: File is not found:',
     &                              exc_name(nexc)(1:exc_size(nexc)),
     &                              '# This file is excluded from ',
     &                               'the dependency list.'
                   endif
#endif
                 endif  !<-- new name
               endif  !<-- new name
             endif  ! <-- size >0
           endif  ! <-- itype>0
         endif ! <-- ierr==0
#ifdef TRAP_DOUBLE_SLASH
! Detect places where double slash occur // inside non-commented
! Fortran code and report line number.  Double slash is used as string
! concatenation operation in Fortran, however CPP interprets it as the
! beginning of C++-style comment, resulting in removal of the rest of
! the line including the slashes.  This is dangerous because it may
! result in correct Fortran syntax (hence no warning at compilation
! time), but altered and incorrect code.

         if ( string(is:is)/='!' .and. string(is:is)/='C' .and.
     &        string(is:is)/='*' .and. string(is:is)/='c' ) then
           do i=is,ie-1
             if (string(i:i)=='/' .and. string(i+1:i+1)=='/') then
               write(*,'(A/2A,I5,1x,3A/A/5(2A/),A)') '*****',
     &        '***** WARNING!!! Double slash /',        '/ on line',
     &         line,    'in file ''',  fname(n)(1:lfnm(n)),  '''.',
     &        '*****',  '***** Although it is legal in Fortran as ',
     &        'string concatenation operator, double',
     &        '***** slash is interpreted as C++ comment line ',
     &        'by GCC C-preprocessor leading to',
     &        '***** possible misinterpretation of Fortran code. ',
     &        'It shall be therefore policy',
     &        '***** of this code to eliminate all double slashes. ',
     &        'Use at least one blank ',
     &        '***** space between slashes in Fortran ',
     &        'concatenation operator.', '*****'
             endif
           enddo
         endif
#endif
#ifdef WARN_NONSTANDARD
! This part checks for presence of symbols which are not part of
! Fortran standard, i.e., tab-characters.  This applies only to .F
! and .h files, but not to Fortran modules, hence .mod are avoided.

         if (fname(n)(lfnm(n)-3:lfnm(n)) /= '.mod') then
           do i=is,ie
             k=ichar(string(i:i))
             if (string(i:i)==' ') then

             elseif (k==9) then
               write(*,'(/1x,4A,I5,1x,A,I3,A/)') '***** WARNING!!! ',
     &           'file ', fname(n)(1:lfnm(n)),  ' line',  line, ':',
     &          i, '  contains TAB character.'
             elseif (k<33 .or. k>126) then
               write(*,'(/1x,4A,I5,1x,A,I3,1x,A,I4,1X,2A/)') '***** ',
     &           'WARNING!!! file ',  fname(n)(1:lfnm(n)),  ' line',
     &            line,':',i, 'nonstandard symbol:',k, string(i:i),'.'
             endif
           enddo
       endif
#endif
#ifdef EXCLUDE
   6   continue
#endif
       if (not_end_of_file) goto 2
       close (unit=iin)
   7   if (n<nmax) goto 1


!          *******       ****       *****    *********
!          ***   ***   ***  ***   ***   ***  *  ***  *
!          ***   ***  ***    ***  ***           ***
!          ***   ***  ***    ***   *******      ***
!          *******    ***    ***        ***     ***
!          ***         ***  ***   ***   ***     ***
!          ***           ****       ****        ***


#ifdef VERBOSE
      write(*,'(2(/1x,A)/)') 'Post-processing...', '==== ============='
#endif

! Suppress circular dependencies associated with the F90-modules
! declared and used within the same file. In this case  the .mod file
! is created automatically at the compilation without the need for
! pre-existing .mod file.

      do k=1,nmods
        j=mod_src(k)     !<-- source code with module declaration
        i=mod_indx(k)    !<-- index of module name within the catalog
        if (mtrx(i,j)/=' ') then
          mtrx(i,j)=' '
#ifdef VERBOSE
        write(*,*)    'Removed self-dependency of ''',
     &              fname(i)(1:lfnm(i)), ''' from ''',
     &                     fname(j)(1:lfnm(j)), '''.'
#endif
        endif
      enddo

      do i=1,nmax                  ! Suppress main diagonal: needed to
        mtrx(i,i)=' '              ! avoid duplication of the file name
      enddo                        ! in the case of self-inclusion.

      i=0                          ! Investigate possible secondary
      j=1                          ! dependencies. This is equivalent
      do while(j>i)             ! to logical chain:
        i=j                                   ! if
        do n=1,nmax                           !     file N depends on K
          do k=1,nmax                         ! and
            if (mtrx(k,n) /= ' ') then        !     file K depends on J
              do j=1,nmax                     ! then
                if (mtrx(j,k) /= ' ' .and.    !     file N depends on J
     &              mtrx(j,k) /= 'M'  .and.   !
     &              mtrx(j,n) == ' ') then    ! however, do not follow
                 mtrx(j,n)=mtrx(k,n)          ! on "K depends on J" if
                endif                         ! K is an F90 module
              enddo                           ! declared inside J.
            endif
          enddo                    ! This is an iterative procedure,
        enddo                      ! because staged include statements
        j=0                        ! are possible. The number non-blank
        do n=1,nmax                ! elements in the matrix grows when
          do k=1,nmax              ! new dependencies are discovered.
            if (mtrx(k,n) /= ' ')  ! The procedure repeats itself until
     &                  j=j+1      ! until the next iteration does not
          enddo                    ! discover any new dependencies
        enddo                      ! (the termination is controlled by
      enddo                        ! the outermost "while" statement).


#ifdef VERBOSE
      if (nmax<=max_string_size) then
        write(*,'(/A)') 'Printing dependence matrix:'
        j=1
        do i=1,nmax
          if (j>9) j=0
          string(i:i)=char(j+48)
          j=j+1
        enddo
        write(*,'(32x,A)') string(1:nmax)
        do j=1,nmax
          do i=1,nmax
            if (mtrx(i,j) /= ' ') then
              if (i==j) then
                string(i:i)='W'
               else
                string(i:i)=mtrx(i,j)
              endif
            else
              string(i:i)='.'
            endif
          enddo
          write(*,'(A,I3,1x,A)') fname(j)(1:28), j, string(1:nmax)
        enddo
      else
        write(*,'(/1x,2A/)') 'Dependence matrix cannot be printed ',
     &                       'out because nmax > max_string_size.'
      endif
#endif

!        *******     ********    ***   *     ***  *********
!        ***   ***   ***   ***   ***   **    ***  *  ***  *
!        ***   ***   ***   ***   ***   ***   ***     ***
!        ***   ***   ***   **    ***   ***** ***     ***
!        *******     ******      ***   *** *****     ***
!        ***         ***  **     ***   ***   ***     ***
!        ***         ***   ***   ***   ***    **     ***

! Report statistics for all files:

      write(iout,'(A/A,12x,A/A,11x,A)')    '#',
     &   '#',   'P R O J E C T    S T A T I S T I C S',
     &   '#',  '---------------  ---------------------'

      write(iout,'(A/2A,13x,I4/A/A,26x,I4)') '#', '# Number of ',
     &     'files given for dependency analysis:', last_arg, '#',
     &     '# Number of newly discovered files:',  n_disc
#ifdef EXCLUDE
      write(iout,'(A/2A,9x,I4)') '#', '# Number of files ',
     &           'excluded from dependency analysis:', nexc
#endif

      write(iout,'(A,13x,A/A,13x,A,9x,I4/A/A)')                  '#',
     &    '---------------------------------------------------', '#',
     &    'total files analyzed for dependencies:',  nmax,  '#', '#'

      write(iout,'(A,6x,I6/2A/A/2A,5x,I6/A)')
     &  '# Total number of lines of Fortran code in all files:',
     &   lfort_all, '# (excluding CPP-commands, empty, and comment ',
     &  'lines)',  '#',  '# Total number of CPP-preprocessor lines ',
     &  'in all files:',   cpps_all, '#'

      write(iout,'(2A,3x,I3/A)')       '# Total number of compiler ',
     &              'directive lines, Open MP, etc...', dirs_all, '#'

      write(iout,'(A,14x,I6/A/A/A,16x,I6/A)')
     &  '# Total number of comment lines in all files:',    comm_all,
     &  '# (excluding empty and containing only "!", "C", "c", "*")',
     &  '#', '# Total number of empty lines in all files:',empty_all,
     &  '# (including lines containing only "!", "C", "c", and "*")'

      write(iout,'(A,38x,A/A,38x,A,7x,I7/)')
     &        '#', '--------------------------', '#', 'total lines:',
     &               lfort_all+cpps_all+dirs_all+comm_all+empty_all

! Write dependency lines associated with F90 modules. This needs to be
! done separately because of the special treatment: since .mod file is
! created after the compilation of source file with module statement,
! the suffix of the source file name (usually .F or .f, but sometimes
! having more than one letter like .f90, etc) needs to be changed into
! .o to initiate compilation of the source  code to create .mod.
! In the segment below it is assumed that the source file suffixes
! always starts with "F" or "f" (a warning is issued otherwise) but
! is not restricted to be a single letter.

      do k=1,nmods
        i=1
        do while (mod_name(k)(i:i)/=' ' .and. i<max_name_size)
          i=i+1
        enddo
        if (mod_name(k)(i:i)==' ') i=i-1
        string=fname(mod_src(k))
        size=lfnm(mod_src(k))
        j=size
        do while (string(j:j) /= '.' .and.  j>1)
          j=j-1
        enddo
        if (string(j:j)=='.') then
          j=j+1
          if (string(j:j)=='F' .or. string(j:j)=='f') then
            string(j:j)='o'
            write(iout,'(3A/)') mod_name(k)(1:i), ': ', string(1:j)
          else
            write(iout,'(6A)') '# WARNING: Cannot generate dependency',
     &                 ' for ', mod_name(k)(1:i), ': its source file ',
     &               fname(mod_src(k))(1:size), ' has unknown suffix.'
          endif
        else
          write(iout,'(6A)')   '# WARNING: Cannot generate dependency',
     &                 ' for ', mod_name(k)(1:i), ': its source file ',
     &                    fname(mod_src(k))(1:size), ' has no suffix.'
        endif
      enddo

! Print out list of dependencies:  Because compilation of .F files is
!------ --- ---- -- ------------   expected to be done two stages, at
! first .F --> .f via CPP  followed by .f --> .o, there are two lines
! of dependency corresponding to each .F, one is .o:.f containing only
! dependencies due to F90-style use and F77-style include statements,
! but not CPP-style #include; conversely f.:.F line has only #include
! dependencies. This way each line contains all necessary dependencies,
! but nothing extra, so, if desired, the first stage can be executed
! separately with all correct dependencies in place. All other types
! of files generate only a single dependency line. [This mechanism is
! activated by CPP-switch TWO_STAGE_UPCSF.]

! Suppress printing dependency line if the .o file depends only on
! a single source and does not need any .h or .mod. This is achieved
! by initializing j=0 for each dependency line and then changing it
! to j=1 (or 2 during the second stage if two-stage .o:.F rule) when
! a non-trivial dependency is discovered.  Do not print if j remains
! less that the treshhold 1 (or 2).

      do n=1,last_arg
        if (lfnm(n)>0) then
          if (n==1 .or. j>0) write(iout,*) !<-- just skip a line
#ifdef SKIP_TRIVIAL
          j=0
#else
          j=3
#endif
          i=lfnm(n)
          string(1:i)=fname(n)(1:i)
#ifdef TWO_STAGE_UPCSF
          if (fname(n)(i-1:i)=='.F') then
            string(i:i+2)='o: '
            i=i+2
            string(i+1:i+lfnm(n)-1)=fname(n)(1:lfnm(n)-1)
            i=i+lfnm(n)
            string(i:i)='f'                ! This segment prints .o:.f
                                           ! dependency, which does not
            do k=1,nmax                    ! need .h files, however F90
              if (mtrx(k,n) /= ' ') then   ! modules are required here.
                if (lfnm(k)>4) then
                  if (mtrx(k,n)=='U' .or. mtrx(k,n)=='u' .or.
     &                mtrx(k,n)=='J' .or. mtrx(k,n)=='j') then
  8                 if (i+lfnm(k) < max_width-2) then
                      j=1
                      i=i+1
                      string(i:i)=' '
                      string(i+1:i+lfnm(k))=fname(k)(1:lfnm(k))
                      i=i+lfnm(k)
                    else                          !--> this branch is
                      i=i+1                       !    to fold printed
                      string(i:i)=backslash       !    dependency line
                      write(iout,'(A)') string(1:i)
                      i=0                         !    if it is too
                      goto 8                      !    long to fit into
                    endif                         !    max_string_size.
                  endif
                endif
              endif
            enddo
            if (j>0) write(iout,'(A)') string(1:i)

            i=lfnm(n)                  ! Prepare the beginning of the
            string(1:i)=fname(n)(1:i)  ! string for printing .f:.F
            string(i:i+2)='f: '        ! dependency by the next segment
            i=i+2
            string(i+1:i+lfnm(n))=fname(n)(1:lfnm(n))
            i=i+lfnm(n)
          else
#endif
            do while (string(i:i)/='.' .and. i>1)
              i=i-1
            enddo                 ! This while is here because file
            string(i+1:i+3)='o: ' ! suffix may be more than one letter.
            i=i+3
            string(i+1:i+lfnm(n))=fname(n)(1:lfnm(n))
            i=i+lfnm(n)
#ifdef TWO_STAGE_UPCSF
          endif  !<-- =='.F'
#endif
          do k=1,nmax
            if (mtrx(k,n)/=' ') then
#ifdef TWO_STAGE_UPCSF
! .f:.F dependency requires only .h files included by CPP, but not
! for .mod and F77-style includes.

              if ( mtrx(k,n)=='I' .or. mtrx(k,n)=='i'  .or.
     &             fname(n)(lfnm(n)-1:lfnm(n) ) /= '.F' ) then
#endif
  9             if (i+lfnm(k) < max_width-2) then
                  j=2
                  i=i+1
                  string(i:i)=' '
                  string(i+1:i+lfnm(k))=fname(k)(1:lfnm(k))
                  i=i+lfnm(k)
                else                             !--> this branch is
                  i=i+1                          !    to fold printed
                  string(i:i)=backslash          !    dependency line
                  write(iout,'(A)') string(1:i)  !    if it is too
                  i=0                            !    long to fit into
                  goto 9                         !    max_string_size
                endif
#ifdef TWO_STAGE_UPCSF
              endif
#endif
            endif
          enddo
          if (j>1) write(iout,'(A)') string(1:i)
        endif
      enddo
      goto 100

  98  write(*,'(/1x,2A/12x,A/)')  '### ERROR: Too many files to ',
     &                   'process. Edit file ''cross_matrix.F'',',
     &          'increase parameter ''max_names'' and recompile.'
      goto 100
#ifdef EXCLUDE
  99  write(*,'(/1x,2A/12x,A/)')    '### ERROR: Encountered too ',
     &      'many files to be excluded. Edit ''cross_matrix.F'',',
     &          'increase parameter ''exc_names'' and recompile.'
#endif
 100  close (iout)
      stop
      end
