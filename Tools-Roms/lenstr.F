      function lenstr (string)

! Find the position of the last non-blank character in input string
! after removing all leading blank characters, if any.     At first,
! find the length of input string using intrinsic function "len" and
! starting from this, search for the last and the first non-blank
! character, "ie" and "is".    If there are leading blanks (is > 1),
! move the whole string to the beginning.     Return "lenstr" is the
! position of the last non-blanc character of the modified string.
!
! WARNING: if there are leading blank characters, user must ensure
! that the string is "writable", i.e., there is a character variable
! in the calling program which holds the string: otherwise call to
! lenstr results in segmentation fault, i.e. passing directly typed
! argument like
!
!                     lstr=lenstr(' x...')
!
! is not allowed, however
!
!                      lstr=lenstr('x...')
!
! is OK because lenstr makes no attempt to shift the string.

      implicit none                       ! In the code below there
      integer lenstr, is,ie               ! are two possible outcomes
      character string*(*)                ! of the search for the first
      ie=len(string)
      do while (string(ie:ie)==' ' .and. ie>1)
        ie=ie-1
      enddo                               ! non-blank character "is":
      is=1
      do while (string(is:is)==' ' .and. is<ie)
        is=is+1
      enddo                               ! it either finds one, or
      if (string(is:is)/=' ') then      ! the search is terminated
        if (is>1) string=string(is:ie) ! by reaching the condition
        lenstr=ie-is+1                    ! (is==ie), while the
      else                                ! character is still blank,
        lenstr=0                          ! which  means that the
      endif                               ! whole string consists of
      return                              ! blank characters only.
      end
