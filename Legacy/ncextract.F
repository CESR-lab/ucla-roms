      program ncextract

! Purpose: to exstract a subset of declarations from "netcdf.inc",
! which are relevant to ROMS code. To do so, first search through the
! whole ROMS code (that is all *.F and *.h files) and build catalog
! of all "nf_" words which occur in the code. Once the catalog is
! built, it read the standard "netcdf.inc" file line-by-line and save
! only lines, which contains nf_words which are present in the
! catalog.

! NOTE: This measure is provided purely for convenience just to speed
! up things a little bit by shortening the "netcdf.inc" include file
! (typically from 1300+ lines of the standard "netcdf.inc" file down
! to 140 lines, which leads to more compact and easier to search .f
! files).  It does not affect the result of compilation of the entire
! model. Some manual editing of the resultant "netcdf.inc.ext" file
! may ne needed in order to eliminate conflicts during subsequent
! compilation --- typically eliminate usage of few obsolecent ncxxxxx
! parameters, which are provided by UNIDATA for compatibility with
! older versions of netCDF library.

      implicit none
      integer, parameter :: max_string=72,  max_size=32,
     &                      max_found=512, input=11, iout=12
      character(len=max_string) string, buff, nfword(max_found)
      character(len=1) quote, double_quote


      integer nlen(max_found), nfound, last_arg, n,
     &    case_fold, i,j, ie,is, istrt,length, line,iocheck
      logical not_end_of_file, lswtch, lprint

#ifdef INT_IARGC
      integer iargc
#endif

      nfound=0
      quote=char(39)
      double_quote=char(34)
      case_fold=ichar('a')-ichar('A')

      last_arg=iargc()
      if (last_arg > 0) then
        do n=1,last_arg+1
          if (n < last_arg+1) then
            call getarg(n, string)
          else
            string='/usr/local/include/netcdf.inc'
            open(unit=iout, file='netcdf.inc.ext',form='formatted')
          endif
          open(unit=input, file=string, form='formatted')
          line=0
          not_end_of_file=.true.
          do while (not_end_of_file)
            line=line+1
            string(1:1)=' '
            read(input,'(A)',iostat=iocheck,end=1) string
            goto 2
   1        not_end_of_file=.false.
   2        length=0
            lswtch=.false.
            do i=1,max_string,+1
              if (string(i:i)==quote) then
                lswtch=.not.lswtch
              elseif (string(i:i)=='!' .and. .not.lswtch) then
                goto 3
              endif
              if (string(i:i)/=' ') length=i
            enddo
   3        continue
            do while(length > 1 .and. string(length:length) == ' ')
              length=length-1
            enddo
            if (string(length:length) == ' ') length=length-1
            istrt=1
            do while(istrt < length .and. string(istrt:istrt) == ' ')
              istrt=istrt+1
            enddo

! Save the non-modified copy of the string into buffer (only for
! "netcdf.inc" file, then to simplify comparison and extraction fold
! cases to lower case and suppress all symbols except letters, digits
! and underscore.

            if (n == last_arg+1) buff=string(1:length)
            do i=istrt,length
              if (string(i:i) >= 'A' .and. string(i:i) <= 'Z') then
                string(i:i)=char(ichar(string(i:i))+case_fold)
              elseif (string(i:i) >= 'a'.and.string(i:i) <= 'z') then

              elseif (string(i:i) >= '0'.and.string(i:i) <= '9') then

              elseif (string(i:i) == '_') then

              else
                string(i:i)=' '
              endif
            enddo

! Search for nf_words inside the string; once found, check whether it
! is already in the catalog, and if it is not, put it in (all passes,
! except the last, where input file is "netcdf.inc") during the last
! pass, if, as opposite, nf_word is in catalog, set switch to print
! the whole non-modified version of the line into the new file.

            lprint=.false.
            i=istrt
            do while (i < length)
              if (string(i:i+2)=='nf_') then
                is=i
              else
                is=0
              endif
              i=i+3
              do while (string(i:i) /= ' ')
                i=i+1
              enddo
              if (is>0) then
                ie=i-1
                lswtch=.true.
                do j=1,nfound
                  if (ie-is+1==nlen(j)) then
                    if (string(is:ie)==nfword(j)(1:nlen(j))) then
                      lswtch=.false.
                    endif
                  endif
                enddo
                if (lswtch .and. n/=last_arg+1) then
                  nfound=nfound+1
                  nlen(nfound)=ie-is+1
                  nfword(nfound)=string(is:ie)
                  write(*,'(3A)') '''',
     &              nfword(nfound)(1:nlen(nfound)), ''''
                elseif (.not.lswtch .and. n==last_arg+1) then
                  lprint=.true.
                endif
              endif
              do while(string(i:i)==' ' .and. i<length)
                i=i+1
              enddo
            enddo
            if (lprint) write(iout,'(A)') buff(1:length)
          enddo
          close (unit=input)
        enddo
        close (unit=iout)
      else
        write(*,'(/1x,2A/, /8x,A/)') '#### ERROR: No filenames ',
     &          'are specified. Usage of ncextract should be',
     &          'ncextract *.F *.h'
      endif
      end
