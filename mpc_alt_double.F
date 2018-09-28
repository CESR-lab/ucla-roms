This part is extracted from "mpc.F". It is kept for reference
purposes only and is considered obsolete, since a better
performing version is present in mpc under

#ifdef DOUBLE_CONST_ALT
!
! Recognize numerical constants of real type in the text of the
! program and convert them into double precision constants, for
! example, 1. --> 1.D0;  .1 --> .1D0;  9.81 --> 9.81D0;
! .5e-8 --> .5D-8;  1.2e+14 --> 1.2D+14 etc.
!
!  Alternative Algorithm:
!------------------------
!    (1) Form list of indices of all dot symbols within the line
!        [dots which occur within quotations ' ... ' are ignored,
!        to do so , variable m is used as masking switch, it turns
!        off when meeting a quotation mark, when entering a region
!        between quotes and turns ond back again, when exiting];
!
!    (2) scan this list in order to recognize fortran logical
!        relation operators, such as ==, /=, .or., >, <,
!        >=, <=, /=, .and., .not.; exclude the dots bounding
!        these expressions from the list of indices. Here setting of
!        indx(j,j+1)=0 signals that the associated dots should be
!        ignored during the subsequent processing; Note that it is
!        not necessary to recognize and exclude dots surrounding
!        logical constants .true. and .false., since there is no
!        possibility in a legal fortran code that there is a digit
!        in a position immediately adjacent to the left or right
!        from  these dots.
!
!    (3) Among remaining dots, exclude those for which none of the
!        characters adjacent on the left or on the right
!        is a digit [i.e. symbol with  ASCII table number between 48
!        to 57 (inclusive)]; It should be noted that (2) and (3) are
!        not interchangeable: (2) must be done before (3), because
!        algorithm (2) relies on parity [opening and closing dot
!        for the logical operator]. Algorithm (3) may destroy the
!        parity, and, as the result, some of the closing dots may
!        not be excluded in (2), if (3) is applied before (2).
!
!    (4) Process the remaining dots: Once at least one of the symbols
!        adjacent to the dot on the left or on the right is a digit
!        [i.e. a symbol with  ASCII table number between 48 to 57
!        (inclusive)], then this dot is associated with
!        either
!               a real type constant expression.
!        or
!               a specifier in operator format like F8.4 for example.
!
!        To exclude dots associated with format specifiers, search
!        for the first nonblank non-digital character to the left
!        from the dot. If this symbol is 'f', 'F', 'e', 'E', 'd' and
!        'D', the dot should be excluded. [For the sake of simplicity
!        the algorithm is designed to exclude all letter characters,
!        i.e. ASCII symbols within the ranges (65:90) and (97:122)
!        inclusive]. So that things like B12.5 are still interpreted
!        as format specifiers, even thought they are not legal
!        fortran expressions: who cares, they should not appear
!        anyway in a correct fortran program.
!
!
!          If the dot is not excluded as a format specifier, scan
!          the symbols to the right from the dot in order to find
!          the first non-blank non-digital symbol. This symbol may
!          be ether
!             'e' or 'E', so that it needs to be converted into 'D';
!          or
!             'd' or 'D', so that the constant is already double
!                         precision and no action is needed (it is
!                         converted to the upper case 'D', however,
!                         for stylistic reasons);.
!          or
!               any other symbol, which is not a digit. In this case
!                         'D0' is placed just after it, while the
!          rest of the line is moved two positions to the right, and
!          the length of the line in increased by two every time when
!          the expression 'D0' is being placed.
!
!        Finally, the alternative part, which prints the error
!        message is placed here just to test bullet proofness of the
!        algorithm: in principle, it should never happen, since all
!        possible possibilities are _positively_ identified in the
!        steps presented above.
!
! Limitations:
!-------------
! i.  This algorithm implies that all fortran logical relation
!     operators are written in dense mode, that is no blank spaces
!     are allowed within fortran logical operations, e.g.,
!     .a n d. is not an acceptable format for .and.;
!
! ii. Real type numerical constant may be written in nondense format,
!     for example  3. 141 5926 will be recognized and interpreted
!     correctly (that is converted into 3. 141 5926D0), however, at
!     least one digit should be adjacent (without leaving blank space
!     between) to the decimal point, e.g., 3 . 1415926 will not be
!     recognized as a real type constant.
!
        ndots=0                                          ! Step (1)
        m=1
        do i=7,length
          if (symbol(i)==quote) m=-m
          if (m>0 .and. symbol(i)=='.') then
            ndots=ndots+1
            indx(ndots)=i
          endif
        enddo
        do j=1,ndots                                     ! Step (2)
         i=indx(j)
          if (i>0) then
            if (symbol(i+3)=='.') then
              if (symbol(i+2)=='r' .or. symbol(i+2)=='R') then
                if (symbol(i+1)=='o' .or. symbol(i+1)=='O') then
                  indx(j)=0
                  indx(j+1)=0
                endif
              elseif(symbol(i+2)=='q' .or. symbol(i+2)=='Q') then
                if ( symbol(i+1)=='e' .or. symbol(i+1)=='E') then
                  indx(j)=0
                  indx(j+1)=0
                endif
              elseif(symbol(i+2)=='t' .or. symbol(i+2)=='T') then
                if ( symbol(i+1)=='g' .or. symbol(i+1)=='G' .or.
     &               symbol(i+1)=='l' .or. symbol(i+1)=='L') then
                  indx(j)=0
                  indx(j+1)=0
                endif
              elseif(symbol(i+2)=='e' .or. symbol(i+2)=='E') then
                if ( symbol(i+1)=='g' .or. symbol(i+1)=='G' .or.
     &               symbol(i+1)=='l' .or. symbol(i+1)=='L' .or.
     &               symbol(i+1)=='n' .or. symbol(i+1)=='N') then
                  indx(j)=0
                  indx(j+1)=0
                endif
              endif
            elseif (symbol(i+4)=='.') then
              if (symbol(i+1)=='a' .or. symbol(i+1)=='A') then
                if (symbol(i+2)=='n' .or.symbol(i+2)=='N') then
                  if (symbol(i+3)=='d'.or.symbol(i+3)=='D') then
                    indx(j)=0
                    indx(j+1)=0
                  endif
                endif
              endif
              if (symbol(i+1)=='n' .or. symbol(i+1)=='N') then
                if (symbol(i+2)=='o'.or. symbol(i+2)=='O') then
                  if (symbol(i+3)=='t'.or.symbol(i+3)=='T') then
                    indx(j)=0
                    indx(j+1)=0
                  endif
                endif
              endif
            endif
          endif
        enddo

        do j=1,ndots                                     ! Step (3)
         i=indx(j)
          if (i>0) then
            m=ichar(symbol(i-1))
            k=ichar(symbol(i+1))
            if ((m<48 .or. m>57) .and.
     &          (k<48 .or. k>57)) then
              indx(j)=0
            endif
          endif
        enddo

        do j=ndots,1,-1   !<-- REVERSED !!!              ! Step (4)
          if (indx(j)>0) then
            m=ichar(symbol(i-1))
            k=ichar(symbol(i+1))
            if ((m>=48 .and. m<=57) .or.
     &          (k>=48 .and. k<=57)) then
              m=indx(j)
  8            m=m-1
                if (symbol(m)==' ') go to 8
               k=ichar(symbol(m))
               if (k>=48 .and. k<=57) then
                 if (m>7) goto 8
               endif
              if (k<65 .or. (k>90. and. k<97) .or.
     &                                        k>122) then
                m=indx(j)
  9              m=m+1
                  if (symbol(m)==' ') go to 9
                 k=ichar(symbol(m))
                if (k>=48 .and. k<=57) then
                   if (m<=length)  goto 9
                endif
                if (symbol(m)=='e' .or. symbol(m)=='E') then
                  symbol(m)='D'
                elseif (symbol(m)=='d' .or. symbol(m)=='D') then
                  symbol(m)='D'
                else
  91              m=m-1
                   if (symbol(m)==' ') goto 91
                  m=m+1
                  do k=length,m,-1
                    symbol(k+2)=symbol(k)
                  enddo
                  symbol(m)='D'
                  symbol(m+1)='0'
                  length=length+2
                endif
              endif
            else
              write(iout,*)'MPC ERROR: cannot recognize expression'
              write(iout,*)'  associated with (.) in the next line.'
            endif
          endif
        enddo
#endif


