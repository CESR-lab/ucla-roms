! CASE = undefined, 2, 3 are for integer kind=4,2,8 
! CASE = 4 or 8 for real type of kind=4 or 8 respectively; 

#ifdef CASE
# if CASE == 4
#  define read_str_TYPE read_str_real
#  define var_TYPE real(kind=4)
# elif CASE == 8
#  define read_str_TYPE read_str_dble
#  define var_TYPE real(kind=8)
# elif CASE == 2
#  define read_str_TYPE read_str_int2
#  define var_TYPE integer(kind=2)
# elif CASE == 3
#  define read_str_TYPE read_str_int8
#  define var_TYPE integer(kind=8)
# endif
#else
# define read_str_TYPE read_str_int
# define var_TYPE integer
# define VERBOSE
#endif

      function read_str_TYPE(str, value)

! Reads numerical value from a character string while avoiding the use
! of FORTRAN internal read.  Literally, this is a functional equivalent
! of
!           read(str,*,iostat=ierr) value

! where variable "value" may be either integer or floating-point,
! single or double precision type.  This is needed as a workaround
! for Intel compiler issue #548345, which is a failure to return
! non-zero iostat error status for certain strings if the operation
! fails [it is known that strings starting with symbols , / F T f t
! cause read operation to yield value=-1.000 and iostat=ierr=0,
! instead of returning non-zero errors status].

! input:  str -- a character string of unknown length;
! output: read_str_TYPE -- integer error status, normally 0; if an
!                  error occurs, then the position of the leftmost
!                  character in the str, which causes the error;
!         value -- an integer- or real-kind, single/double precision
!                  variable (the type should match name of the
!                  function itself); note that if an error occurs
!                  then variable "value" REMAINS UNTOUCHED.

      implicit none
      integer read_str_TYPE
      character(len=*) str
      var_TYPE value
      integer(kind=8) mnts, mnts_bak  !<-- to prevent/detect overflows
      integer is,ie, ierr, isgn
      
#if CASE > 3
      logical dcm_point, exp_present
      integer xpsgn, k,n
#else
      var_TYPE itst  !<-- test value to detect integer overflow
#endif
      ierr=0 ; ie=len(str)                        ! This logic is to
      if (ie > 0) then                            ! avoid hitting zero
        do while(ie > 1 .and. str(ie:ie) == ' ')  ! position when input
          ie=ie-1                                 ! is a zero-size
        enddo                                     ! string.
        if (str(ie:ie) == ' ') ie=ie-1

        isgn=1     ! sign of the whole number
        mnts=0     ! abs.val. mantissa (integer,ignoring decimal point)
#if CASE > 3
        exp_present=.false. ! presence of exponent
        xpsgn=1             ! sign of exponent
        dcm_point=.false.   ! presence of decimal period, logical
        n=0                 ! number of decimal places after period
        k=0                 ! decimal exponent
#endif
        if (str(1:1) == '+') then       ! Determine sign of the whole
          is=2                          ! number, check if here is
        elseif(str(1:1) == '-') then    ! either + or -  sign present
          is=2 ; isgn=-1                ! in the fist character; no
        else                            ! sign means positive; shift
          is=1                          ! "is" to the second position,
        endif                           ! if the sign is present.

        is=is-1
        do while( is < ie .and.
#if CASE > 3
     &     .not.exp_present .and.
#endif
     &               ierr == 0 )
          is=is+1

! Normally mantissa is expected to grow as digits are being read from
! left to right. If there is a decimal point then the number of digits
! on the right from it is counted as well, so it can be divided later;
! however there is a possibility that the number has too many digits
! past decimal point, so integer(kind=8) mantissa may overflow.
! To prevent this, save the previous value of the mantissa, and revert
! back if it no longer grows (algorithm to detect this situation
! utilizes the fact that over-limiting positive integer results in
! negative integer).  This means that the digits beyond this point
! are ignored and the number is rounded off.

          if (str(is:is) >= '0' .and. str(is:is) <= '9') then
            mnts_bak=mnts ; mnts=mnts*10 + ichar(str(is:is))-48
#if CASE > 3
            if (mnts > mnts_bak) then !--> normal growth of mantissa
              if (dcm_point) n=n+1
            elseif (mnts == 0) then   !--> still going through zeros
              if (dcm_point) n=n+1    !          in front of the first
            elseif (dcm_point) then   !          significant digit
              mnts=mnts_bak           !<-- revert back and ignore the
            else                      !    subsequent digits (overflow)
#else
            if (mnts < mnts_bak) then
#endif
              ierr=1000
#ifdef VERBOSE
              write(*,*) '### ERROR: read_string :: mantissa overflow'
#endif
            endif
#if CASE > 3
          elseif (str(is:is) == '.' .and. .not.dcm_point) then
            dcm_point=.true.
          elseif (str(is:is) == 'E' .or. str(is:is) == 'D' .or.
     &            str(is:is) == 'e' .or. str(is:is) == 'd') then
            exp_present=.true. ; is=is+1
          elseif (str(is:is) == '+' .or. str(is:is) == '-') then
            exp_present=.true.
#endif
          else                          ! The above search may or may
            ierr=is                     ! not yield the presence of
#ifdef VERBOSE
            write(*,*) '### ERROR: read_string :: unexpected non-digit'
#endif
          endif                         ! exponent, 'E' or 'E+' or just
        enddo                           ! + or - in the middle followed
                                        ! followed by digits. If found,
                                        ! record its magnitude as "k"
        if (ierr == 0) then             ! and sign as "xpsgn".
#if CASE > 3
          if (exp_present .and. is <= ie) then
            if (str(is:is) == '+') then
              is=is+1
            elseif (str(is:is) == '-') then
              is=is+1 ; xpsgn=-1
            endif
            do while(is <= ie .and. ierr == 0)
              if (str(is:is) >= '0' .and. str(is:is) <= '9') then
                k=k*10 + ichar(str(is:is))-48
              else
                ierr=is   !--> ERROR, not a digit
              endif
              is=is+1
            enddo
          endif

c**       write(*,*) 'exp_present =',exp_present, 'xpsgn =',xpsgn
c**       write(*,*) 'mnts=', mnts, 'k,n =',k,n

          if (ierr == 0) then           ! Finalize computation of
            value=isgn*dble(mnts)       ! floating-point number: the
            k=xpsgn*k-n                 ! 5 incoming variables are all
            if (k > 0) then             ! integers:
              xpsgn=1                   !   "mnts" (>0) mantissa;
            else                        !   "isgn" sign of mantissa;
              xpsgn=-1 ; k=-k           !   "n" the number of digits
            endif                       !       after decimal point;
            do while(k > 31)            !   "k" (>0) decimal exponent 
              k=k-32                    !       as it is read from 
              if (xpsgn > 0) then       !       the input string; 
                value=value*1.0D+32     !   "xpsgn" sign of exponent.
              else                      !
                value=value*1.0D-32     ! The code on the left is
              endif                     ! basically to compute
            enddo                       !
            do while(k > 15)            !  isgn*mnts * 10^(xpsgn*k-n)
              k=k-16                    !
              if (xpsgn > 0) then       ! however a special care is 
                value=value*1.0D+16     ! needed to avoid overflows 
              else                      ! within the intermediate 
                value=value*1.0D-16     ! operations and/or round-off 
              endif                     ! errors due to underflows 
            enddo                       ! in the middle.
            if (k > 7) then
              k=k-8
              if (xpsgn > 0) then
                value=value*1.0D+8
              else
                value=value*1.0D-8
              endif
            endif
            if (k > 3) then
              k=k-4
              if (xpsgn > 0) then
                value=value*1.0D+4
              else
                value=value*1.0D-4
              endif
            endif
            if (k > 1) then
              k=k-2
              if (xpsgn > 0) then
                value=value*1.0D+2
              else
                value=value*1.0D-2
              endif
            endif
            if (k > 0) then
              k=k-1
              if (xpsgn > 0) then
                value=value*10.D0
              else
                value=value*0.1D0
              endif
            endif
          endif
#else
          mnts=isgn*mnts             ! This silly-looking code is to
          itst=mnts ; mnts_bak=itst  ! detect integer(kind<8) overflow:
          if (mnts == mnts_bak) then ! here "mnts" and "mnts_bak" are
            value=mnts               ! integer(kind=8) type, but "itst"
          else                       ! and "value" are kind=2 or 4, so
            ierr=+32768  !--> ERROR  ! if "itst" overflows, copying it
                                     ! back to will cause difference.
# ifdef VERBOSE
            write(*,*) '### ERROR: read_string :: integer overflow'
# endif
          endif
#endif
        endif
      else
        ierr=1   !--> ERROR, zero-size str
#ifdef VERBOSE
        write(*,*) '### ERROR: read_string :: zero-size string'
#endif
      endif
      read_str_TYPE=ierr
      end

#undef var_TYPE
#undef read_str_TYPE
#ifndef CASE
# define CASE 2
#elif CASE == 2
# undef CASE
# define CASE 3
#elif CASE == 3
# undef CASE
# define CASE 4
#elif CASE == 4
# undef CASE
# define CASE 8
#elif CASE == 8
# undef CASE
#endif
#ifdef CASE
# include "read_string.F"
#endif
