      program chkindent

! Check indenting of nested #ifdef-s and report line numbers with
! miss-aligned CPP directives.     Traditionally one-position step to
! the right is expected for each inner #ifdef in the nesting hierarchy.
! This corresponds to undefined setting of switch FLEX_OFFSET below.
! Defining it make this code to accept formatting with any offset as
! long as it is used consistently within each nest. In this case the
! offset is determined automatically as the offset of the first CPP
! directive within each #ifdef, and all subsequent directives are
! expected to be aligned with it.
! This code produces no output on the screen if no miss-alignment is
! found and does not attempt to modify the file being checked.

c--#define FLEX_OFFSET

      implicit none
      integer, parameter :: iin=11,   max_string_size=512,
     &                      max_dpth=16, max_fname_size=80
      integer nargs, line, ndpth, indx, ierr, i,is,ie,
     &                              lstr, ipos(0:max_dpth)
      character(len=max_string_size) :: fname
      character(len=max_string_size) :: string


      nargs=iargc()
      if (nargs==1) then
        call getarg(nargs,fname)
        ie=max_fname_size
        do while (fname(ie:ie)==' ' .and. ie>1)
          ie=ie-1
        enddo
        is=1
        do while (fname(is:is)==' ' .and. is<ie)
          ie=ie-1
        enddo
        open(unit=iin,  file=fname(is:ie),   form='formatted',
     &                              status='old', iostat=ierr)
        write(*,*) 'checking ''', fname(is:ie), '''....'
        line=0
        ndpth=0
        do i=0,max_dpth
#ifdef FLEX_OFFSET
          ipos(i)=0
#else
          ipos(i)=i+2  !<-- 1-position offset step
#endif
        enddo

        do while (ierr==0)
          do i=1,max_string_size
             string(i:i)=' '
          enddo
          read(iin,'(A)',iostat=ierr) string
          if (ierr==0) then
            line=line+1
            if (string(1:1)=='#') then
              lstr=max_string_size
              do while(string(lstr:lstr)==' ' .and. lstr>1)
                lstr=lstr-1
              enddo
              if (lstr>0) then
                is=0
                ie=0
                i=1
                do while (ie==0 .and. i<lstr)
                  i=i+1
                  if (is==0) then
                    if (string(i:i)/=' ') is=i
                  else
                    if (string(i:i)==' ') ie=i-1
                  endif
                enddo
                if (ie==0) ie=lstr

                if (string(is:is+1)=='el'  .or.
     &              string(is:is+2)=='end') then
                  indx=ndpth-1
                else
                  indx=ndpth
#ifdef FLEX_OFFSET
                  if (ipos(ndpth)==0) ipos(ndpth)=is
#endif
                endif

                if (is/=ipos(indx)) then
                  write(*,'(A,I5,2x,A,2(2x,A,I3))') 'line =', line,
     &                   string(is:ie), 'should start at position',
     &                                ipos(indx), 'instead of', is
                endif

                if (string(is:is+1) == 'if') then
                  ndpth=ndpth+1
#ifdef FLEX_OFFSET
                  ipos(ndpth)=0
#endif
                elseif (string(is:is+2) == 'end') then
                  ndpth=ndpth-1
                endif

              endif ! <-- lstr>0
            endif  ! <-- string(1:1)=='#'
          endif  ! <--  ierr==0
        enddo  ! <-- while
        close(unit=iin)
      else
        write(*,'(/1x,3A/1x,A/24x,A/ /1x,A,2(/1x,2A)/)') 'This is ',
     &   '"chkindent" - a basic tool to check indenting of CPP ',
     &   'directives.',   'Usage should be:',  'chkindent filename',
     &   'with a single argument - name of the file to be checked.',
     &   '"Chkindent" only reports the line numbers where ',
     &   'misalignments occur ',   'without making an attempt to ',
     &                                         'change the file.'
      endif
      stop
      end


