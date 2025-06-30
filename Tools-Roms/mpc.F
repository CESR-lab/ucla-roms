      program mpc

! Multi-functional pre-compiler to be inserted between CPP and Fortran
! compiler to:
!------------
! (1) Eliminate blank lines and garbage left by CPP as well as Fortran
!     comments, but keep compiler directives like C$OMP DO PARALLEL in
!     place. The essential purpose of this feature is to make the CPP-
!     pre-processed code to be clean and readable, without the intend
!     to change any of its functionality of results.

! (2) Normal Fortran syntax does not prohibit the practice when a
!     quotation is opened on one line and closes on the next line.
!     Doing so, however, may be considered as unsafe and sloppy from
!     a stylistic point of view; mpc enforces the rule when every
!     string inside quotation mark must be closed within the line.
! --> Optional:: activated by TRAP_UNMATCHED_QUOTES CPP-switch.

! (3) Convert default-size Fortran real type declarations into real*8.
!     Any mixture of lower- and upper-case letters is allowed in the
!     word real. Explicitly specified size declaration, such as real*4
!     or real*16 will be left unchanged.
! --> Optional:: activated by REAL_TO_REAL8 CPP-switch.

! (4) Enforce double precision accuracy for the real type constants
!     appearing anywhere in the program, for example, 1. --> 1.D0;
!          .1 --> .1D0;    9.81 --> 9.81D0;   .5e-8 --> .5D-8;
!         1.2e+14 --> 1.2D+14 etc.
! --> Safe. Will distinguish real type constants from look alike
!     combinations like, say format descriptors F8.4 or E12.5; or
!     logical expressions like 1.EQ.i etc.  Will fold line of the
!     code automatically according to Fortran syntax rules, when
!     this change causes line to increase beyond 72 characters.
! --> Optional: activated by DOUBLE_CONST CPP-switch.

! (5) Convert _default_ size Fortran integer type declarations into
!     integer*4, while at the same keeping explicitly declared sizes
!     like integer*2 or integer*8 unchanged.
! --> Optional: activated when both INT_TO_INT4 and REAL_TO_REAL8
!     CPP-switches are defined.

! (6) Convert parallelized loops over subdomains (tiles)
!
!          do tile=0,NSUB_X*NSUB_E-1
!
!     into two nested loops. The outer is loop over parallel threads
!     (to be parallelized), and the inner is over tiles to be processed
!     by the particular thread:
!
!          do trd=0,numthreads-1
!            do tile=forward_sweep
!     or
!          do trd=0,numthreads-1
!            do tile=backward_sweep
!
!     in such a way the inner loops in the subsequent parallel regions
!     are reversed, so that each thread "zig-zags" across the tiles it
!     is processing. Automatically append SHARED/PRIVATE lists with the
!     newly introduced variables.
!   --> Optional: activated by ZIG_ZAG CPP-switch.

! (7) Automatically insert !dir$ ivdep directives in front of each
!     i-loop which has bounds istr:iend or imin:imax, including their
!     combinations and modifications like istrR and others alike, as
!     well as additions/subtractions of bounds, like imin-1, imax+2.
!     Normally these are the innermost loops in ROMS, which are
!     expected to be vectorized. This is known to affect performance
!     to some degree, and the effect can be seen in Intel compiler
!     vectorization report when -qopt-report=5 -qopt-report-phase=vec
!     flags are invoked.  Activated by CPP-option VECTORIZE below.

! (8) Automatically insert calls to "sync_trap(barrier_number)" just
!     in front of each Open MP directive BARRIER.   This is merely to
!     detect loss of synchronization events and it does not affect
!     computational results provided (of course) that the code does
!     not do out of sync due to parallel bugs. Activated by CPP-option
!     TRAP_BARRIERS.

! More detailed description for each option is available below.

! In the case when after they were modified, lines of the Fortran
! program became longer than the prescribed width of 72 characters,
! they are automatically folded by mpc according to the Fortran
! syntax rules.

! Usage:  Overall similar to CPP:
!-------
!                         mpc source.file target.file
!                    or
!                         mpc source.file > target.file
!                    or
!                         mpc < source.file > target.file
!
! mpc is smart enough to recognize how many arguments are given, TWO,
! ONE or NONE and act accordingly.  If two arguments are present, the
! first one is input file name, while the second is output, if only
! one is present, it is input; output goes to standard output; if NONE
! mpc expects input from standard input, while output goes to standard
! output. The last option allows mpc to work as receiver in pipe with
! CPP:
!                /lib/cpp -P file.F | mpc > file.f
!
! Created and maintained by Alexander Shchepetkin, old_galaxy@yahoo.com


#define INT_TO_INT4
#define REAL_TO_REAL8
#define DOUBLE_CONST

#define F77_TO_F95
#define KIND_STYLE

#define TRAP_UNMATCHED_QUOTES

c-#define VECTORIZE

#define ZIG_ZAG
c--#define ZIG_ZAG_OLD
c--#define RUTGERS
c--#define TRAP_BARRIERS

      implicit none
      integer, parameter :: max_length=128, max_filename=32
      character(len=max_filename) file_in, file_out
      character(len=2*max_length) str
      character(len=max_length) scratch
#ifdef F77_TO_F95
     &                        , bffr
#endif
      character(len=1) tab, quote, double_quote, type
      character(len=1), parameter :: type_none=' ', char_type='1',
     &                               int_type='2',  real_type='3'
      logical not_end_of_file, omp_dir, dir_switch, lswtch
      integer iin, iout, last_arg, case_fold, ierr, line, istr,
     &  lstr, ibarr, i,j,k,m, is,ie, isft, ks0,ks1,ks2,ks3,ks4,ks5
#ifdef DOUBLE_CONST
      integer ndots, indx(max_length)
#endif

! Most modern Fortran compilers recognize "iargc" as an intrinsic
! function, so it does not need to be declared. IBM xlf95 is a notable
! exclusion. So do 7.x and earlier versions of Intel IFC compiler, but
! starting with version 8.x IFORT recognizes it as intrinsic. Portland
! Group PGF compiler also needs it to be declared.

#ifdef INT_IARGC
      integer iargc
#endif

#ifdef ZIG_ZAG_OLD
      logical par_region, lnest
      integer ip1,ip2, is1,is2
      lnest=.false.  ; par_region=.false.
#endif

! Note 2016: Fortran 2003 ISO Standard provides an intrinsic module
! called "iso_fortran_env" (recognized by a sufficiently new compiler
! without explicit declaration) to access standard I/O units, while
! the legacy stdin=5 and stdout=6 units are considered obsolescent and
! may be no longer recognized by the newer compilers in future.
! Because of this the use of units 5 and 6 is avoided within this
! program by replacing them with read(*, ) and write(*, ) respectively,
! while at the same time the use of module "iso_fortran_env" is avoided
! as well to be compatible with the legacy compilers.  Therefore the
! approach is to initialize both units specifiers, "iin" and "iout",
! to zero, then, check whether one or two arguments are present, open
! the corresponding files and assign the units to positive numbers.
! After that, every time when reading and writing is performed, check
! whether iin and/or iout are positive and use them as unit specifiers,
! otherwise use * for standard I/O.  This is somewhat cumbersome,
! especially for iout (logical checks are inserted in many occasions),
! but it is compatible with both old and new ways without any
! adjustment to the code and/or CPP-switches.

      iin=-1 ; iout=-1 ; ibarr=0
      tab=char(9) ; quote=char(39) ; double_quote=char(34)
      case_fold=ichar('a')-ichar('A')
      lswtch=.false. ;  dir_switch=.false.
      do i=1,max_filename
        file_in(i:i)=' ' ; file_out(i:i)=' '
      enddo
      last_arg=iargc()
      if (last_arg > 0) then
        iin=11 ; call getarg(1,file_in)
        open(iin, file=file_in, form='formatted', status='old')
      endif
      if (last_arg == 2) then
        iout=12 ; call getarg(2,file_out)
        open(iout, file=file_out, form='formatted', status='unknown')
      endif

      not_end_of_file=.true.
      line=0
  1    line=line+1                      ! Reset string length and
        lstr=0 ; omp_dir=.false.        ! blank out the string itself,
        do i=1,max_length               ! then read in a new string.
          str(i:i)=' '
        enddo
        if (iin > 0) then
          read(iin,'(A)',iostat=ierr,end=2) str
        else
          read(  *,'(A)',iostat=ierr,end=2) str
        endif
        goto 3
  2     not_end_of_file=.false.
  3     if ( str(1:1) =='!' .or. str(1:1) == '*' .or.
     &       str(1:1) =='C' .or. str(1:1) == 'c' ) then
          if ( str(2:2) == '$' .or. str(3:3) == '$' .or.
     &         str(4:4) == '$' .or. str(5:5) == '$' ) then
            omp_dir=.true.
          else               ! if the first str indicates that this
            goto 1           ! line is a Fortran comment, but a dollar
          endif              ! sign is present within positions 2,...5
        endif                ! then the line is a directive and should
                             ! be processed further.
        type=' '
        do i=1,max_length
          if (str(i:i) == quote) then        ! Meaning of symbols:
            if (type == ' ') then            ! 's' inside single quotes
              type='s'                       ! "d" inside double quotes
            elseif (type == 's') then        ! "'S'" inside single
              type=' '                       !           inside double
            elseif (type == 'd') then        ! '"D"' inside double
              type='S'                       !           inside single
            elseif (type == 'S') then        ! '"X'" illegal
              type='d'                       ! "'X"' illegal
            endif
          elseif (str(i:i) == double_quote) then
            if (type == ' ') then
              type='d'                       ! Exclamation sign !
            elseif (type == 'd') then        ! outside quotes starts
              type=' '                       ! Fortran 90 comment which
            elseif (type == 's') then        ! needs to be suppressed;
              type='D'                       ! "!" within a nest of
            elseif (type == 'D') then        ! quotes is just a symbol
              type='s'                       ! of a character string.
            endif
          elseif ( str(i:i) == '!' .and. type == ' ' .and.
     &                   .not.(omp_dir .and. i == 1) ) then
            goto 4
          endif
          if (str(i:i) /= ' ') lstr=i
        enddo
   4    continue

#ifdef TRAP_UNMATCHED_QUOTES
        if (type /= ' ') then
          if (iout > 0) then
            write(iout,'(/6x,2A,I4/)') '### ERROR: Unmatched ',
     &                                 'quote on line', line
          else
            write(   *,'(/6x,2A,I4/)') '### ERROR: Unmatched ',
     &                                 'quote on line', line
          endif
        endif
#endif

! Restrict "lstr" to the last non-blank and non-tab character of
! the actual Fortran code (that is ignoring the F90 commentaries
! starting with an exclamation sign outside quotes).

        if (lstr > 0) then
          do while( lstr > 1 .and. ( str(lstr:lstr) == ' ' .or.
     &                                 str(lstr:lstr) == tab ) )
            lstr=lstr-1
          enddo
          if (str(lstr:lstr) == ' ' .or. str(lstr:lstr) == tab) then
            lstr=lstr-1
          endif
        endif

! Find the starting non-blank and non-tab character of the line "istr".

        if (lstr > 0) then
          istr=1
          do while( istr < lstr .and. ( str(istr:istr) == ' ' .or.
     &                                  str(istr:istr) == tab ) )
            istr=istr+1
          enddo
        else
          if (not_end_of_file) then
            goto 1     ! empty line ==> proceed with the next one
          else
            goto 99    ! empty line and file ends ==> close everything
          endif        ! and quit (not having logical check here causes
        endif          ! "mpc" to terminate abnormally in the case when
                       ! the input file has an empty line at the end
                       ! because "istr=0" and the code is allowed to
                       ! proceed beyond this point.

c**     write(*,*) 'line=', line, ' istr=', istr,
c**  &      ' not_end_of_file =', not_end_of_file


! Recognize Fortran type declaration, REAL, INTEGER, CHARACTER (any
! combination of upper- and lowercase letters is permissible) and turn
! them default-precision REAL into REAL*8 or REAL(kind=8) depending on
! how relevant CPP-switches are set, however do not convert explicitly
! typed real*4 into real*8. Default-size integers may-be converted to
! integer*4 or integer(kind=4),

        i=istr ; type=type_none
#if defined REAL_TO_REAL8 || defined F77_TO_F95
        if (str(i  :i  )=='R' .or. str(i  :i  )=='r') then
         if (str(i+1:i+1)=='E' .or. str(i+1:i+1)=='e') then
          if (str(i+2:i+2)=='A' .or. str(i+2:i+2)=='a') then
           if (str(i+3:i+3)=='L' .or. str(i+3:i+3)=='l') then
             i=i+3
             type=real_type
           endif
          endif
         endif
        endif
#endif
#if defined INT_TO_INT4 || defined F77_TO_F95
        if (str(i  :i  )=='I' .or. str(i  :i  )=='i') then
         if (str(i+1:i+1)=='N' .or. str(i+1:i+1)=='n') then
          if (str(i+2:i+2)=='T' .or. str(i+2:i+2)=='t') then
           if (str(i+3:i+3)=='E' .or. str(i+3:i+3)=='e') then
            if (str(i+4:i+4)=='G' .or. str(i+4:i+4)=='g') then
             if (str(i+5:i+5)=='E' .or. str(i+5:i+5)=='e') then
              if (str(i+6:i+6)=='R' .or. str(i+6:i+6)=='r') then
                i=i+6
                type=int_type
              endif
             endif
            endif
           endif
          endif
         endif
        endif
#endif
#ifdef F77_TO_F95
        if (str(i  :i  )=='C' .or. str(i  :i  )=='c') then
         if (str(i+1:i+1)=='H' .or. str(i+1:i+1)=='h') then
          if (str(i+2:i+2)=='A' .or. str(i+2:i+2)=='a') then
           if (str(i+3:i+3)=='R' .or. str(i+3:i+3)=='r') then
            if (str(i+4:i+4)=='A' .or. str(i+4:i+4)=='a') then
             if (str(i+5:i+5)=='C' .or. str(i+5:i+5)=='c') then
              if (str(i+6:i+6)=='T' .or. str(i+6:i+6)=='t') then
               if (str(i+7:i+7)=='E' .or. str(i+7:i+7)=='e') then
                if (str(i+8:i+8)=='R' .or. str(i+8:i+8)=='r') then
                  i=i+8
                  type=char_type
                endif
               endif
              endif
             endif                      ! After this moment "i" is
            endif                       ! index of the last character
           endif                        ! of Fortran type declaration
          endif                         ! which may be either "real",
         endif                          ! "integer", or "character".
        endif                           ! Find is index "j" of the
#endif
        if (type /= type_none) then     ! first non-blanc symbol
          j=i+1
          do while (j < lstr .and. str(j:j)==' ')
            j=j+1
          enddo                         ! after type declaration.

! Once it is established that the line contains an "obsolescent" from
! Intel IFC compiler point of view style of explicit size declaration
! like real*16 or character*128, change it to Fortran 95 style,
!
!                   real*16  --> real(16)
!              character*128 --> character(len=128)
!
! which is then accepted by the compiler without warning message.
! To do so, first extract the size itself and save it into "bffr"
! (strip outer braces if any), and set initial value "isft" to be
! the number of deleted characters (negative), which is asterisk
! '*' itself and all blank spaces occurring between Fortran type
! specifier and size specifier.

          if (str(j:j) == '*') then
#ifdef F77_TO_F95
            is=j+1
            do while( str(is:is) == ' ' .and. is < lstr)
              is=is+1
            enddo
            ie=is
            do while(ie < lstr .and. str(ie+1:ie+1) /= ' '
     &                         .and. str(ie+1:ie+1) /= ',')
              ie=ie+1
            enddo
            isft=i+1-is
            if (str(is:is) == '(' .and. str(ie:ie) == ')') then
              is=is+1 ; ie=ie-1 ; isft=isft-2
            endif
            bffr(1:1+ie-is)=str(is:ie)

! Finally, modify the line of the code: add the number of newly
! inserted characters to "isft", move the rest of the line to the
! right in order to make room for the insertion, adjust the length
! of the line accordingly, then insert the new-style size
! specification.

            if (type == char_type) then
              isft=isft+6
              str(ie+1+isft:lstr+isft)=str(ie+1:lstr)
              lstr=lstr+isft
              str(i+1:i+5)='(len='
              str(i+6:i+ie-is+6)=bffr(1:ie-is+1)
              str(i+7+ie-is:i+7+ie-is)=')'

            elseif (type == int_type .or. type == real_type) then
# ifdef KIND_STYLE
              isft=isft+7
# else
              isft=isft+2
# endif
              str(ie+1+isft:lstr+isft)=str(ie+1:lstr)
              lstr=lstr+isft
# ifdef KIND_STYLE
              str(i+1:i+6)='(kind='
              str(i+7:i+7+ie-is)=bffr(1:1+ie-is)
              str(i+8+ie-is:i+8+ie-is)=')'

# else
              str(i+1:i+1)='('
              str(i+2:i+2+ie-is)=bffr(1:1+ie-is)
              str(i+3+ie-is:i+3+ie-is)=')'
# endif
            endif
#endif

! The following code segment deals with declaration of default-size
! real and integer type variables:  if so directed, insert explicit
! specification of the size.  Once again, two versions, F77 and F90
! are provided. The logical condition here is explained as follows:
! at first, reject declaration where the first symbol after type
! declaration is opening bracket (this is possible only if it is
! already a F90/95 style declaration with explicit size specification
! ==> no further action is required). If the symbol is not opening
! bracket, then three possibilities may occur: (i) it is a comma,
! separating type and F90 attribute (like pointer, dimension, etc.,);
! (ii) it is still a blanc character, which means that the line
! contains just one word, say REAL, and continues on the next line
! (typical situation in RUTGERS code); or (iii) there symbols other
! that '*' (already considered above) and '(', separated by at least
! one blank space: these are either variable names, or F90 double
! colon :: separator. In the all three cases (i -- iii) the line is
! considered as default type declaration, and should be processed.

#if defined REAL_TO_REAL8 || defined INT_TO_INT4 || defined F77_TO_F95
          elseif ( str(j:j) /= '(' .and. ( str(j:j) == ' ' .or.
     &                      j > i+1 .or. str(j:j) == ',' ) ) then

c--#if defined REAL_TO_REAL8 && defined INT_TO_INT4
c--            if (type==int_type .or. type==real_type) then
c--#elif defined REAL_TO_REAL8 || defined F77_TO_F95
c--            if (type==real_type) then
c--#elif defined INT_TO_INT4
c--            if (type==int_type) then
c--#endif

            if (type /= char_type) then
              if (type == int_type) then
                type='4'
              elseif (type == real_type) then
# ifdef REAL_TO_REAL8
                type='8'
# else
                type='4'
# endif
              endif
# ifdef F77_TO_F95
#  ifdef KIND_STYLE
              isft=8
#  else
              isft=3
#  endif
# else
              isft=2
# endif
              str(j+isft:lstr+isft)=str(j:lstr)

                                         ! Move the rest of the line
                                         ! isft symbols to the right
                                         ! to make room for explicit
              do k=j,j+isft-1            ! declaration of type size
                str(k:k)=' '             ! and increase lstr of
              enddo                      ! the line accordingly, then
              lstr=lstr+isft             ! insert size declaration.
# ifdef F77_TO_F95
#  ifdef KIND_STYLE
              str(i+1:i+6)='(kind='
              str(i+7:i+7)=type
              str(i+8:i+8)=')'
#  else
              str(i+1:i+1)='('
              str(i+2:i+2)=type
              str(i+3:i+3)=')'
#  endif
# else
              str(i+1:i+1)='*'
              str(i+2:i+2)=type
# endif
            endif
#endif
          endif
        endif

#ifdef DOUBLE_CONST

! Recognize numerical constants of real type in the text of the program
! and convert them into double precision constants, e.g.,
!
!           1. --> 1.D0;        .1 --> .1D0;      9.81 --> 9.81D0;
!        .5e-8 --> .5D-8;  1.2e+14 --> 1.2D+14  etc.
!
! Algorithm:
!-----------
! (1) Form list of indices of all dots within the line, except dots
!     which occur within quotations ' ... '. To do so, local logical
!     switch "lswtch" is used as a masking switch, it turns OFF when
!     meeting a quotation mark, when entering a region between quotes
!     and turns back ON again, when exiting.
!
! For each dot character in the list, starting from the last one,
! and moving from the right to the left:
!
! (2) Scan the characters adjacent from the _left_ to the dot in order
!     to find the first non-blank character which is not a digit (for
!     this purpose digits are ASCII symbols with numbers within the
!     range (48:57) inclusive). This search is terminated if
!     either
!               a non-blank non-digital symbol is found
!     or
!               the 7th position (the starting position for the
!               fixed format Fortran statements) has been reached.
!
!     During this search also save "is", which is position of the
!     leftmost digit among digits adjacent to the dot on the left
!     side (position of the dot itself, if there are no adjacent
!     digits on the left of it).
!
! (3) Check weather this symbol is .NOT. a letter, that is excluding
!     ASCII characters with numbers within the ranges of numbers
!     (65:90) or (97:122) inclusive, _OR_ weather the 7th position
!     has been reached (in this case it does not matter what the
!     symbol is). If either condition is true, continue processing,
!     otherwise terminate it.
!
!     It should be noted that in a legal Fortran code a constant
!     expression may be preceded by a mathematical operation symbol,
!     bracket, comma, dot, etc; but _NOT_ with a letter. If it
!     happens, the potential candidate for the numerical constant is
!     actually a format statement descriptor, like E16.8, and not a
!     constant of real type.  These are rejected at this moment.
!
! (4) Once condition (3) holds, scan the characters adjacent on the
!     right to the dot in order to find the first non-blanc character
!     which is non a digit. This search is limited by the length of
!     the line in the case when no such symbol is found (if it is the
!     case, it is then interpreted as a blanc symbol).   During this
!     search also save "ie", which is position of the last digit among
!     the digits adjacent to the dot on the right (it is set to the
!     position of the dot itself, if no adjacent digits are present).
!     Along with the previously saved "is" [see (2) above], "ie"
!     forms a logical expression "is<ie", indicating that there is at
!     least one digit adjacent to the dot, so further processing is
!     required.
!
! (5) Once this symbol is found, if any, or the search was terminated
!     (in this case index m is equal to lstr+1, so that that symbol
!     is blank, or ! (Fortran comment), this symbol may be
!     ether
!          'e', 'E', 'd' or 'D', so that it likely belongs to the
!          constant itself. In this case scan to the right, to verify
!          that this symbol is '+', '-' or a digit. If so, the
!          expression is a floating point real type constant to be
!          converted into double precision.
!     or
!          an underscore symbol '_', which may be associated with a
!          Fortran 90/95-style real-valued constant with explicit type
!          (kind), e.g., 1.2_4 (in which case no it is accepted "as
!          is" and no further action is taken); or a RUTGERS-style
!          constant -- e.g., 1.0_r8 or, 1.23_e8+18 -- in this case it
!          is converted into standard double precision syntax -->
!          1.0D0 and 1.23D+18.
!     or
!          any other character. In this case verify that it is not a
!          letter. (In a legal Fortran code a constant expression may
!          be followed by a mathematical/logical operation, bracket,
!          coma, etc, but NEVER a letter. If so, the expression is a
!          fixed point real valued constant to be converted into
!          double precision. Move the rest of the line two positions
!          to the right to make room and paste "D0" immediately after
!          its last digit.
!
! NOTE: if there are more than one real constant within the line in
! the code, the order of processing is from the right to the left.
! This is needed because in the case when 'D0' is pasted to the
! constant as in the step (5), the second case, the tail of the line
! is shifted to the right. Processing them in the forward order will
! also cause shift of the dots to be processed later. So that the
! indices indx(j+1:ndots) are no longer consistent with the actual
! position of the dots, if the dot indx(j) was found to be a fixed
! point real type constant as defined in (5), second case.
!
! Limitations:  1. Real-valued constants should not be continued
!-------------     to the next line of Fortran code.

# define k illegal
        ndots=0 ; lswtch=.true.                              ! Step (1)
        do i=7,lstr
          if (str(i:i) == quote) then
            lswtch=.not.lswtch
          elseif (str(i:i) == '.' .and. lswtch) then
            ndots=ndots+1 ; indx(ndots)=i
          endif
        enddo
        do j=ndots,1,-1  ! <-- REVERSED (see NOTE above)     ! Step (2)
          m=indx(j) ; is=m ; lswtch=.true.
          do while(m > 7 .and. lswtch)
            m=m-1
            if (str(m:m) >= '0' .and. str(m:m) <= '9') then
              is=m
            elseif (str(m:m) /= ' ') then
              lswtch=.false.
            endif
          enddo
                                                             ! Cond.(3)
          if ( lswtch .or. str(m:m) < 'A' .or.
     &        ('Z' < str(m:m) .and. str(m:m) < 'a')
     &                   .or. 'z' < str(m:m) ) then
            m=indx(j) ; ie=m ; lswtch=.true.                 ! Step (4)
            do while(m < lstr .and. lswtch)
              m=m+1
              if (str(m:m) >= '0' .and. str(m:m) <= '9') then
                ie=m
              elseif (str(m:m) /= ' ') then
                lswtch=.false.
              endif
            enddo
            if (lswtch) m=m+1

            if (is < ie) then                                ! Step (5)
              if (str(m:m) == 'e' .or. str(m:m) == 'E' .or.
     &            str(m:m) == 'd' .or. str(m:m) == 'D') then
                i=m+1
                do while(i<lstr .and. str(i:i)==' ')
                  i=i+1
                enddo
                if ( str(i:i) == '+' .or.  str(i:i) == '-' .or.
     &              (str(i:i) >= '0' .and. str(i:i) <= '9')
     &                                           ) str(m:m)='D'

              elseif (str(m:m) == '_') then

! Rutgers compatibility mode: the following code segment searches for
! suffices _r8 and _e8 attached to real valued constants and converts
! it into standard Fortran double precision constants, e.g.,
!               1.0_r8 --> 1.0D0;    2.4_e8+3 --> 2.4D+3.
! NOTE: there is still possibility for Fortran 90/95 real valued
! constant format with explicit type, such that 1.0_4 or 2.3_8.
! These are left unchanged.

# ifdef RUTGERS
                if (str(m+2:m+2) == '8') then
                  if (str(m+1:m+1) == 'r') then
                    str(m:m+1)='D0'
                    str(m+2:lstr-1)=str(m+3:lstr)
                    lstr=lstr-1
                  elseif (str(m+1:m+1) == 'e') then
                    str(m:m)='D'
                    str(m+1:lstr-2)=str(m+3:lstr)
                    lstr=lstr-2
                  endif
                endif
# endif
              elseif ( str(m:m) < 'A' .or.  ( 'Z' < str(m:m) .and.
     &                 str(m:m) < 'a' ) .or.  'z' < str(m:m) ) then
                do while(str(m-1:m-1) == ' ')
                  m=m-1
                enddo
                str(m+2:lstr+2)=str(m:lstr) ; lstr=lstr+2

# ifdef KIND_STYLE
                str(m:m+1)='_8'
# else
                str(m:m+1)='D0'
# endif
              endif
            endif
          endif
        enddo !<-- ndots
# undef k
#endif  /* DOUBLE_CONST */
#ifdef VECTORIZE

! Automatically place loop vectorization directives for i-loops

        do i=istr,lstr
          if ( 'A' <= str(i:i) .and. str(i:i) <= 'Z') then
            scratch(i:i)=char(ichar(str(i:i))+case_fold)
          else
            scratch(i:i)=str(i:i)
          endif
        enddo

        i=istr
        if (scratch(i:i+1) == 'do') then
          i=i+2
          do while(i < lstr .and. scratch(i:i) == ' ')
            i=i+1
          enddo
          if (scratch(i:i)=='i') then
            i=i+1
            do while(i < lstr .and. scratch(i:i) == ' ')
              i=i+1
            enddo
            if (scratch(i:i) == '=') then
              i=i+1
              do while(i < lstr .and. scratch(i:i) == ' ')
                i=i+1
              enddo

              if (scratch(i:i+3) == 'istr' .or.
     &            scratch(i:i+3) == 'imin') then
                i=i+4
                do while(i < lstr .and. scratch(i:i) /= ',')
                  i=i+1
                enddo
                if (i < lstr) then
                  i=i+1
                  if (scratch(i:i+3) ==' iend' .or.
     &                scratch(i:i+3) ==' imax') then
                    i=i+3
                    do while(i < lstr .and. scratch(i:i) /= ',')
                      i=i+1
                    enddo
                    if (i == lstr) then
                      if (iout > 0) then
                        write(iout,'(A)') '!dir$ ivdep'

! If the target processor is Itanium 2, suppress vectorization
! directive (available only for Intel IA32 architecture).
! Not doing so causes compiler warnings, but is harmless otherwise.

# if __IFC != IA64
c                       write(iout,'(A)') '!dir$ vector always'
# endif
c                       write(iout,'(A)') '!dir$ prefetch'
c                       write(iout,'(A)') '!dir$ unroll'
                      else
                        write(*,'(A)') '!dir$ ivdep'
# if __IFC != IA64
c                       write(*,'(A)') '!dir$ vector always'
# endif
c                       write(*,'(A)') '!dir$ prefetch'
c                       write(*,'(A)') '!dir$ unroll'
                      endif
                    endif
                  endif
                endif
              endif
            endif
          endif
        endif
#endif
#ifdef ZIG_ZAG

! Recognize and transform parallel loops:
!---------- --- --------- -------- ------
! The following code segment rearranges direction of consecutive
! 'do tile=my_first,my_last' loops into zig-zag order. This version
! is designed to work in conjunction with the OpenMP version of main.F
! (or any other code), so that no parallel directive is required in
! front of the loop in order for the rearrangement to occur: it is
! simply triggered by the occurrence of 'do tile=my_***' with possible
! blank spaces between 'tile', '=' and 'my_***', where 'my_***' may
! be anything like hand coded 'my_first,my_last' or 'my_tile_range'.

! Compared to the straightforward sequence, this measure eliminates
! secondary cache misses, because after passing a synchronization
! point, each thread proceeds with the same tile it was processing
! just before the synchronization point. It also reduces the
! probability of mutual cache_line invalidation (by factor of two)
! in a multiprocessor machine, if multiprocessing is allowed for the
! subdomains adjacent in the XI direction. This is because after a
! subdomain has been processed by a processor and a synchronization
! point has been reached, all cache_lines which are going across the
! subdivision partitioning are coherent with the cache of _that_
! processor (and invalid with respect to the one working on the
! adjacent subdomain). And it is _that_ processor (and not the
! adjacent one), who proceeds with this subdomain, and therefore
! enjoys all its cache_lines coherent with its cache at this moment.

        i=istr
        if (str(i:i+1) == 'do') then
          i=i+2
          do while(i < lstr .and. str(i:i) == ' ')
            i=i+1
          enddo
          if (str(i:i+3) == 'tile') then
            i=i+4
            do while(i < lstr .and. str(i:i) == ' ')
              i=i+1
            enddo
            if (str(i:i) == '=') then
              i=i+1
              do while(str(i:i) == ' ' .and. i < lstr)
                i=i+1
              enddo
              if (str(i:i+2) == 'my_') then
                dir_switch=.not.dir_switch
                if (dir_switch) then
                  scratch='my_first,my_last,+1'
                else
                  scratch='my_last,my_first,-1'
                endif
                str(i:i+18)=scratch(1:19)
                lstr=i+18
              endif
            endif
          endif
        endif
#endif
#ifdef ZIG_ZAG_OLD

! Similarly as above, but designed to work for version of main.F,
! where parallel directives are inserted in front of each parallel
! loop over subdomains "tiles", i.e., loops like
!
!      C$OMP PARALLEL DO PRIVATE(tile)
!          do tile=0,NSUB_X*NSUB_E-1
!
! are converted into two nested loops. The outer is loop over
! parallel threads (to be parallelized), and the inner is over
! tiles to be processed by the particular thread:
!
!          do trd=0,numthreads-1
!            do tile=forward_sweep
! or
!          do trd=0,numthreads-1
!            do tile=backward_sweep
!
! in such a way the inner loops in the subsequent parallel regions
! are reversed, so that each thread zig-zags across the tiles it is
! processing.
                              ! If the line is a parallel directive
        if (omp_dir) then     ! containing attribute LOCAL(..) or
          ip1=0               ! PRIVATE(..), and perhaps, SHARE(..),
          ip2=0               ! identify indices of the first and
          is1=0               ! the last symbol inside each bracket,
          is2=0               ! which are then used below to append
          i=1                 ! either list.
          do while(i < lstr)
            if (str(i:i+6)=='PRIVATE' .or. str(i:i+4)=='LOCAL' .or.
     &          str(i:i+6)=='private' .or. str(i:i+4)=='local') then
              i=i+5
              do while(i < lstr .and. str(i:i) /= '(')
                i=i+1
              enddo
              if (str(i:i)=='(') then
                ip1=i+1
                do while(i < lstr .and. str(i:i) /= ')')
                  i=i+1
                enddo
                if (str(i:i) == ')') ip2=i-1
              endif
            elseif (str(i:i+6) == 'SHARE' .or.
     &              str(i:i+6) == 'share') then
              i=i+5
              do while(i < lstr .and. str(i:) /= '(')
                i=i+1
              enddo
              if (str(i:i) == '(') then
                is1=i+1
                do while(i < lstr .and. str(i:i)/=')')
                  i=i+1
                enddo
                if (str(i:i) == ')') is2=i-1
              endif
           endif
                   ! Skip to the beginning of the next word
                     ! for further checking. NOTE that 'PRIVATE'
                     ! in context of 'DO PRIVATE' is recognized
                      ! (as it should be), but 'PRIVATE' within
                      ! 'THREADPRIVATE' is disregarded because
                            ! it not a complete word.
            i=i+1
            do while(str(i:i) /= ' ' .and. str(i:i) /= ',')
              i=i+1
            enddo
            do while(i < lstr .and. str(i:i) == ' ')
              i=i+1
            enddo
          enddo

          if (ip1 > 0 .and. ip2 == 0) then
            if (iout > 0) then
              write(iout,'(/6x,2A/)') '### ERROR: No closing bracket ',
     &                                     'in local or private list.'
            else
              write(   *,'(/6x,2A/)') '### ERROR: No closing bracket ',
     &                                     'in local or private list.'
            endif
          endif
          if (is1 > 0 .and. is2 == 0) then
            if (iout > 0) then
              write(iout,'(/6x,2A/)') '### ERROR: No closing bracket ',
     &                                               'in shared list.'
            else
              write(   *,'(/6x,2A/)') '### ERROR: No closing bracket ',
     &                                               'in shared list.'
            endif
          endif

! It should be noted that a valid DO PARALEL directive must have at
! least one variable in its LOCAL/PRIVATE list to privatize the index
! of do-loop which follows it. Hence, the presence of non-empty private
! list triggers further processing in the following code segment: paste
! additional variables "trd" and "subs" into private list and
! "numthreads" into shared, in order to allow loop rearrangement.
! There are three possibilities here: SHARED attribute goes before
! PRIVATE, after, and there is no SHARE at all.

                                            ! Save the whole string
          if (ip1 > 0 .and. ip2 > 0) then   ! into buffer array before
            bffr(1:lstr)=str(1:lstr)        ! start messing with it.

            if (ip1 > is2 .and. is2 >= is1 .and. is1 > 0) then
              str(is1:is1+11)='numthreads, '
              str(is1+12:ip1-11)=bffr(is1:ip1-1)
              str(ip1+12:ip1+21)='trd,subs, '
              str(ip1+22:lstr+22)=bffr(ip1:lstr)
              lstr=lstr+22
            else
              str(ip1:ip1+9)='trd,subs, '
              if (is1 > 0 .and. is2 > 0) then
                str(ip1+10:is1+9)=bffr(ip1:is1-1)
                str(is1+10:is1+21)='numthreads, '
                str(is1+22:lstr+22)=bffr(is1:lstr)
                lstr=lstr+22
              else
                str(ip1+10:ip2+11)=bffr(ip1:ip2+1)
                str(ip2+12:ip2+31)=', SHARED(numthreads)'
                str(ip2+32:lstr+30)=bffr(ip2+2:lstr)
                lstr=lstr+30
              endif
            endif

! Write out the modified compiler directive, which is obviously longer
! than the original one, check its lstr and if it exceeds 72 symbols,
! split it into two. For the aesthetic purposes the split goes along a
! natural divider str, such as a intentionally placed blank
! character or a comma.

            if (lstr <= 72) then
              if (iout > 0) then
                write(iout,'(A)') str(1:lstr)
              else
                write(   *,'(A)') str(1:lstr)
              endif
            else    !--> split line longer than 72 characters
              ks1=72
              do while(str(ks1:ks1) /= ' ' .and. ks1 > 0)
                ks1=ks1-1
              enddo
              ks2=72
              do while (str(ks2:ks2) /= ',' .and. ks2 > 0)
                ks2=ks2-1
              enddo
              if (ks1 > 54) then
                k=ks1
              elseif (ks2 > 6) then
                k=ks2
              else
                if (iout > 0) then
                  write(iout,'(/6x,2A/)') '### ERROR: mpc :: Cannot ',
     &                                      'split OpenMP directive.'
                else
                  write(   *,'(/6x,2A/)') '### ERROR: mpc :: Cannot ',
     &                                      'split OpenMP directive.'
                endif
              endif
              if (iout > 0) then
                write(iout,'(A)') str(1:k)
              else
                write(   *,'(A)') str(1:k)
              endif
              m=6
              str(m:m)='&'
              do i=m+1,k
                str(i:i)=' '
              enddo
              str(k+1-lstr+k:k)=str(k+1:lstr)
              if (iout > 0) then
                write(iout,'(A)') str(1:k)
              else
                write(   *,'(A)') str(1:k)
              endif
            endif
            par_region=.true.
            goto 1
          endif
        endif  !<--  omp_dir

! Once a parallel region is detected and the compiler directive above
! was modified, transform the parallel loop over subdomains (tiles)
! into a set of nested loops over threads (outer loop) and subdomains
! within the work zone of each thread (inner loop), in such a way
! that the direction of the inner loop is always reversed with respect
! to the direction of similar loop in the previous parallel region
! (zig-zag sequence).

        if (par_region) then
          is=7
          do while(str(is:is) == ' ' .and. is < lstr)
            is=is+1
          enddo
          if (str(is:is+1) == 'do') then
            is=is+2
            do while(is < lstr  .and. str(is:is) == ' ')
              is=is+1
            enddo
            if (str(is:is+4) == 'tile=') then
              if (str(is+5:is+21) /= '0,NSUB_X*NSUB_E-1') then
                lswtch=.true.
              else
                lswtch=.false.
              endif
              if (lswtch) then
                lnest=.true.
                if (iout > 0) then
                  write(iout,'(7x,A/8x,A)') 'do trd=0,numthreads-1',
     &                              'subs=NSUB_X*NSUB_E/numthreads'
                else
                  write(   *,'(7x,A/8x,A)') 'do trd=0,numthreads-1',
     &                              'subs=NSUB_X*NSUB_E/numthreads'
                endif
                dir_switch=.not.dir_switch
                if (dir_switch) then
                  if (iout > 0) then
                    write(iout,'(9x,A)')
     &               'do tile=subs*trd,subs*(trd+1)-1,+1'
                  else
                    write(   *,'(9x,A)')
     &               'do tile=subs*trd,subs*(trd+1)-1,+1'
                  endif
                else
                  if (iout > 0) then
                    write(iout,'(9x,A)')
     &               'do tile=subs*(trd+1)-1,subs*trd,-1'
                  else
                    write(   *,'(9x,A)')
     &               'do tile=subs*(trd+1)-1,subs*trd,-1'
                  endif
                endif
                goto 1
              endif
            endif
          elseif (lnest) then
            if (str(is:is+4) == 'enddo') then
              lnest=.false.
              par_region=.false.
              if (iout > 0) then
                write(iout,'(9x,A/7x,A)') 'enddo', 'enddo'
              else
                write(   *,'(9x,A/7x,A)') 'enddo', 'enddo'
              endif
              goto 1
            endif
          endif
        endif  !<--  par_region
#endif
#ifdef TRAP_BARRIERS
! instrument barrier directives.

        if (str(1:1) =='C' .or. str(1:1)=='c' .or. str(1:1)=='!') then
          if (str(2:5) == '$OMP' .or. str(2:5) == '$omp') then
            is=6
            do while(is < lstr .and. str(is:is) == ' ')
              is=is+1
            enddo
            if (is < lstr-5  .and. ( str(is:is+6)=='BARRIER' .or.
     &                               str(is:is+6)=='barrier' )) then
              ibarr=ibarr+1
              write(scratch,'(I4)') ibarr
              is=1; ie=4
              do while(is < 4 .and. scratch(is:is) == ' ')
                is=is+1
              enddo
              if (iout > 0) then
                write(iout,'(A,4x,3A)') 'C$', 'call sync_trap(',
     &                                      scratch(is:ie), ')'
              else
                write(   *,'(A,4x,3A)') 'C$', 'call sync_trap(',
     &                                      scratch(is:ie), ')'
              endif
            endif
          endif
        endif
#endif


! Write the modified line of code into the output file.
!------ --- -------- ---- -- ---- ---- --- ------ -----
! Because its lstr after modification may exceed the standard of 72
! symbols, it may be necessary to split the line into two. In this case
! several attempts are made to find a good (from an aesthetic point of
! view) splitting point. First attempt is made by searching for the
! first blank symbol starting from position 72 and moving to the left;
! if such is not found, then the search starts from the beginning and
! attempt is made to split along a comma; if this is not successful, an
! assignment operator is being searched; then a mathematical operation.

        if (lstr > 0 .and. lstr <= 72) then
          if (iout > 0) then
            write(iout,'(A)') str(1:lstr)
          else
            write(   *,'(A)') str(1:lstr)
          endif
        elseif (lstr > 72) then

          ks0=0        ! Find an appropriate place within the
          ks1=0        ! line where it can be split into two.
          ks2=0        ! This is merely aesthetic matter: split
          ks3=0        ! can be made where there is a natural
          ks4=0        ! divider: semicolon, blank character,
          ks5=0        ! comma, or arithmetic operation

          do k=7,72
            if (str(k:k) == ';') then
              ks0=k
            elseif (str(k:k) == ' ') then
              ks1=k
            elseif (str(k:k) == ',') then
              ks2=k
            elseif (str(k:k) == '=') then
              ks3=k
            elseif (str(k:k) == '/') then
              ks4=k
            elseif (str(k:k) == '+' .or.
     &              str(k:k) == '-'  .or.
     &              str(k:k) /= '*') then
              ks5=k
            endif
          enddo

          if (lstr-ks1 > 66) ks1=0
          if (lstr-ks2 > 66) ks2=0
          if (lstr-ks3 > 66) ks3=0
          if (lstr-ks4 > 66) ks4=0
          if (lstr-ks5 > 66) ks5=0

! Make the final decision about the line split: there is nothing
! special here, it is just a matter of aesthetics to decide which
! of the possible breaking points (if more than one are available)
! is the most appropriate; the logical sequence below is designed
! as hierarchy of preferences.

          if (ks0 > 6) then
            k=ks0
          elseif (ks1 > 34) then
            k=ks1
          elseif (ks4 > 6) then
            k=ks4-1
          elseif (ks2 > 54) then
            k=ks2
          elseif (ks3 > 60) then
            k=ks3
          elseif (ks5 > 6) then
            k=ks5-1
          else
            if (iout > 0) then
              write(iout,'(/6x,2A/)') '### ERROR: mpc :: Cannot ',
     &                                              'split line.'
            else
              write(   *,'(/6x,2A/)') '### ERROR: mpc :: Cannot ',
     &                                              'split line.'
            endif
          endif

! Write out the line.  First write the first part of the line; then
! create a continuation line resetting the starting symbols to blank,
! and them moving the tail to the left to center the mathematical
! expression in the continuation line around the center of the area
! limited by positions 7...72.

          if (iout > 0) then
            write(iout,'(A)') str(1:k)
          else
            write(   *,'(A)') str(1:k)
          endif
          if (omp_dir) then     ! Blank out already printed portion
            do i=7,k            ! of the string, however if the line
              str(i:i)=' '      ! is an Open MP directive, the first
            enddo               ! five symbols are its sentinel, so
          else                  ! it must appear on the continuation
            do i=1,k            ! line as well.  On the other hand,
              str(i:i)=' '      ! if it is just a normal line, it may
            enddo               ! contain label which should not be
          endif                 ! repeated.
          str(6:6)='&'
          m=(lstr+k-79)/2       ! <-- center the continuation line
          if (lstr-m > 72) m=lstr-72
          str(k+1-m:lstr-m)=str(k+1:lstr)
          if (iout > 0) then
            write(iout,'(A)') str(1:lstr-m)
          else
            write(   *,'(A)') str(1:lstr-m)
          endif
        endif
       if (not_end_of_file) goto 1
  99  if (iout > 0) close(iout)
      if (iin  > 0) close(iin)
      end
