      program cpp_clean

! A basic tool to clean up empty lines and Fortran comments after CPP.
! Note,this tool leaves CPP hash-lines intact, so it can be applied to
! .F files as well without changing functionality.

#define TRAP_UNMATCHED_QUOTES

      implicit none
      integer, parameter :: max_length=128, max_filename=32
      character(len=max_filename) file_in, file_out
      character(len=7) fmt_in, fmt_out
      character(len=1) tab,quote,double_quote,type,symbol(max_length)
      logical not_end_of_file, omp_dir
      integer last_arg, iin, iout, i, line, ierr, length

! Function "iargc" is viewed as intrinsic by most modern compilers
! and does not need to be declared. IBM xlf95 is a notable exclusion.
! So do 7.x and earlier versions of Intel IFC compiler, but starting
! with version 8.x IFORT recognize it as intrinsic.

#ifdef INT_IARGC
      integer iargc
#endif
      tab=char(9)
      quote=char(39)
      double_quote=char(34)

      file_in='                                '
      file_out='                                '
      last_arg=iargc()
      if (last_arg>0) then
        iin=11
        call getarg(1,file_in)
        open(iin, file=file_in, form='formatted', status='old')
      else
        iin=5                   !<-- read from standard input
      endif
      if (last_arg==2) then
        iout=12
        call getarg(2,file_out)
        open(iout, file=file_out, form='formatted', status='unknown')
      else
        iout=6                  !--> write into standard output
      endif
      write(fmt_in,'(A,I3,A)') '(', max_length, 'A1)'

      not_end_of_file=.true.
      line=0
  1    line=line+1
        length=0                ! Reset the string length and blank
        omp_dir=.false.         ! out the string itself. Then read
        do i=1,max_length       ! in a new string.
          symbol(i)=' '
        enddo
        read(iin,fmt_in, iostat=ierr,end=2) (symbol(i), i=1,max_length)
        goto 3
  2     not_end_of_file=.false.
  3     if ( symbol(1)=='!' .or. symbol(1)=='*' .or.
     &       symbol(1)=='C' .or. symbol(1)=='c' ) then
          if ( symbol(2)=='$' .or. symbol(3)=='$' .or.
     &         symbol(4)=='$' .or. symbol(5)=='$' ) then
            omp_dir=.true.
          else                   ! if the first symbol indicates that
            goto 1               ! this line is a fortran comment, but
          endif                  ! if a dollar sign is present within
        endif                    ! positions 2,...5, the line is an
                                 ! directive and should be processed
        type=' '                 ! further.
        do i=1,max_length
          if (symbol(i)==quote) then     ! Meaning of symbols:
            if (type==' ') then          ! 's' inside single quotes
              type='s'                     ! "d" inside double quotes
            elseif (type=='s') then      ! "'S'" inside single
              type=' '                     !            inside double
            elseif (type=='d') then      ! '"D"' inside double
              type='S'                     !            inside single
            elseif (type=='S') then      ! '"X'" illegal
              type='d'                     ! "'X"' illegal
            endif
          elseif (symbol(i)==double_quote) then
            if (type==' ') then
              type='d'                     ! Exclamation sign !
            elseif (type=='d') then      ! outside quotes starts
              type=' '                     ! Fortran 90 comment which
            elseif (type=='s') then      ! needs to be suppressed;
              type='D'                     ! "!" within a nest of
            elseif (type=='D') then      ! quotes is just a symbol
              type='s'                     ! of a character string.
            endif
          elseif ( symbol(i)=='!' .and. type==' '
     &          .and. .not.(omp_dir .and. i==1) ) then
            goto 4
          endif
          if (symbol(i)/=' ') length=i
        enddo
   4    continue
#ifdef TRAP_UNMATCHED_QUOTES
        if ( type/=' ') write(iout,'(/6x,A,I4,1x,A/)')
     &      '### ERROR: Unmatched quote on line', line
#endif

! Restrict "length" to the last non-blank and non-tab (ASCII code 9)
        if (length>0) then
          do while( length>0 .and. ( symbol(length)==' '
     &                          .or. symbol(length)==tab ))
            length=length-1
          enddo
        endif
        if (length>0) then
          write(fmt_out,'(A,I3,A)') '(', length, 'A1)'
          write(iout,fmt_out) (symbol(i),i=1,length)
        endif
       if (not_end_of_file) goto 1
  99  if (iout>6) close(iout)
      if (iin>5) close(iin)
      stop
      end
