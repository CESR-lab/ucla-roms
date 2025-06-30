c--#define TEST_INSERT  ! <-- self testing

#ifdef TIME_INDEX
# ifdef EXTRACT_INDEX
      subroutine extract_time_index(name, lstr, indx, ierr)
# else
      subroutine   insert_time_index(name, lstr, indx, ierr)
# endif
#else
      subroutine insert_node(name, lstr, node, nnodes, ierr)
#endif

! Insert MPI-node number "node" into character string "name" and
! adjust the length of the string "lstr" to accommodate the inserted
! number. The node number is coded using as many digits as necessary
! to accommodate numbers from 0 to nnodes-1, which is determined
! internally in this subprogram. The number is placed just before
! the suffix, if suffix is present, e.g.: '.nc', or in the end of
! the name, if it does not have suffix.
!
! input:  name   character string
!         lstr   length of the string
!      indx/node time Index/MPI-node number (MPI-rank, within
!                   the from 0 to nnodes-1) to be inserted into
!                    the string.
!        nnodes  total number of MPI nodes, the same as parameter
!                    NNODES in "param.h", but since this routine is
!                    also used in "partit.F", it is passed as an
!                    argument, rather than parameter in include file;
!                    in this code nnodes is used only to determine
!                    how many digits are needed to accommodate all
!                    possible MPI node numbers.
!         ierr   input/output error counter incremented by 1 for each
!                    error occurring here; if no errors it remains
!                    untouched.  Note that "ierr" is ASSUMED TO BE
!                    INITIALIZED FROM THE OUTSIDE, normally to 0.
!
! output: name    modified string with node number inserted.
!         lstr    length of the modified string
!         ierr    same as at entry, if no errors; incremented, if any

      implicit none
      character(len=*) name
      integer lstr, ierr, i,j,k, lsffx, ndots, idot(3)
#ifdef TIME_INDEX
     &                                , indx
#else
     &                                , node,  nnodes
#endif
#ifndef EXTRACT_INDEX
     &                                , power
      character(len=16) sffx
#endif
      integer, parameter :: digits=5
      logical leading_dots

      ndots=0                              ! Determine how many dots
      leading_dots=.true.                  ! are present in the string.
      do i=1,lstr                          ! Do not allow more than 3
        if (name(i:i) == '.') then         ! dots: signal about error,
         if (.not.leading_dots) then       ! if it occurs, however, do
            if (ndots < 3) then            ! not count leading dots,
              ndots=ndots+1                ! so that that names like
              idot(ndots)=i                ! "../grid.00.nc" are OK.
            else
              write(*,'(/1x,4A/)')
#ifdef TIME_INDEX
# ifdef EXTRACT_INDEX
     &                         '### ERROR: extract_time_index :: ',
# else
     &                          '### ERROR: insert_time_index :: ',
# endif
#else
     &                                '### ERROR: insert_node :: ',
#endif
     &         'too many dots in file name ''', name(1:lstr), '''.'
              ierr=ierr+1
              return
            endif
          endif
        else                               ! Determine whether the
          leading_dots=.false.             ! name has a suffix: the
        endif                              ! last segment of the name
      enddo                                ! is considered to be a
                                           ! suffix if it has least one
      lsffx=0                              ! character which is not a
      if (ndots > 0) then                  ! digit or UNIX wildcard.
        i=idot(ndots)+1
   1    if ( (name(i:i) < '0' .or. name(i:i) > '9') .and.
     &      name(i:i) /= '*' .and. name(i:i) /= '?' ) then
          lsffx=lstr-idot(ndots)+1
        elseif (i < lstr) then             ! If suffix is present, it
          i=i+1                            ! will be saved and added
          goto 1                           ! back to the name later.
        endif
      endif                                ! Inspect segments between
                                           ! the dots. These segments
      do j=1,ndots-1                       ! must contain only digits
        i=idot(j)+1
   2    if (name(i:i) < '0' .or. name(i:i) > '9') then
          if (name(i:i) /= '*' .and. name(i:i) /= '?') then
            write(*,'(/1x,2A/20x,3A/)')
#ifdef TIME_INDEX
# ifdef EXTRACT_INDEX
     &                          '### ERROR: extract_time_index :: ',
# else
     &                           '### ERROR: insert_time_index :: ',
# endif
#else
     &                                 '### ERROR: insert_node :: ',
#endif
     &                    'a non-digital character found in index ',
     &                    'segment of name ''',  name(1:lstr), '''.'
            ierr=ierr+1
          endif                            ! or UNIX wild cards '*'
        elseif (i < idot(j+1)-1) then      ! and '?'. These segments
          i=i+1                            ! are used to store Time
          goto 2                           ! Index or MPI-node number.
        endif                              ! Set "ierr" signal the
      enddo                                ! presence of an illegal
      if (ierr /= 0) return                ! symbol in this area.


#ifdef EXTRACT_INDEX
      if (ndots == 1 .and. lsffx == 0) then   ! Read the digital
        i=idot(1)+1 ; j=lstr                  ! segment in the file
      elseif (ndots > 1) then                 ! name which contains
        i=idot(1)+1 ; j=idot(2)-1             ! time index (identified
      else                                    ! as the leftmost segment
        i=0 ; j=0                             ! which has length larger
      endif                                   ! or equal to the value
      indx=0                                  ! of parameter "digits"
      if (j-i+1 >= digits) then               ! specified above).
        do k=i,j
          indx=10*indx + ichar(name(k:k))-48
        enddo
      endif
#else

! Determine where to put Time Index and/or MPI-node number (rank).
! Since the string name may or may not contain digital segments,
! three possibilities exist:
!
! (1) there are no digital segments: either there are no dots
!     (hence no suffix), or there is only one dot (which separates
!     the suffix from the root name). In this case an ne digital
!     segment is created for either Time Index or MPI-node (rank).
!
! (2) only ONE digital segment exist. In this case it has to be
!     determined whether it is to be used as Time Index or MPI-node.
!     The determination is made based upon the length of the segment:
!
!     --> if the segment length is greater or equal than parameter
!     "digits" specified above, then it will be interpreted as the
!     place to store time index. (if MPI-node needs to be inserted,
!     a new digital segment adjacent to the right from the existing
!     one will be created in this case.)
!
!     --> if, in the other hand, the segment length is smaller than
!     "digits", then it will be interpreted as the MPI-node number.
!     A a new digital segment adjacent to the LEFT from the existing
!     one will be created to place Time Index.
!
! (3) There are already TWO digital segments in string "name". The
!     left one will be used for time index, the right for MPI-node.
!
! In the code segment below, "i" is the starting dot of digital
! segment to be inserted,  while "j" is the starting dot of the tail
! of the string "name", i.e. name(j:lstr) contains either suffix of
! the string "name" (including starting dot), or the right digital
! segment (if there is one, and there is no suffix); or both segment
! and suffix.

      if (ndots == 0) then
        i=lstr+1 ; j=lstr+1
        name(i:i)='.'
      else
# ifdef TIME_INDEX
        i=idot(1)
        if (ndots == 1) then
          if (lsffx > 0 .or. lstr-idot(1) < digits) then
            j=idot(1)
          else
            j=lstr+1
          endif
        elseif (ndots == 2 .and. idot(2)-idot(1) <= digits) then
          j=idot(1)
        else
          j=idot(2)
        endif
# else
        if (ndots == 1) then
          i=idot(1)
        elseif (ndots == 2) then
          if (idot(2)-idot(1) <= digits) then
            i=idot(1)
          else
            i=idot(2)
          endif
        else
          i=idot(ndots-1)
        endif
        if (lsffx > 0) then
          j=idot(ndots)
        else
          j=lstr+1
        endif
# endif
      endif
      lsffx=lstr+1-j
      if (lsffx > 0) sffx(1:lsffx)=name(j:lstr)

! Load Time Index or MPI-node (rank) into temporal variable "k".
! This variable will be written into digital segment. Also specify
! the maximum allowed number, which sets the number of digits in the
! segment.

# ifdef TIME_INDEX
      k=indx
      power=10**digits
# else
      k=node
      power=10                           ! Determine how many digits
      do while(nnodes > power)           ! are needed to accommodate
        power=10*power                   ! the largest possible MPI-
      enddo                              ! node number (rank).
      if (power >= 10**digits) then
        write(*,'(/1x,2A/6x,2A/6x,A/)')
#ifdef TIME_INDEX
# ifdef EXTRACT_INDEX
     &                          '### ERROR: extract_time_index :: ',
# else
     &                            '### ERROR: insert_time_index :: ',
# endif
#else
     &                                  '### ERROR: insert_node :: ',
#endif
     &   'Possible ambiguity between MPI-node segment',    'length ',
     &   'and time index segment length. To fix: increase parameter',
     &   '''digits'' in file "insert_node.F" and recompile.'
        ierr=ierr+1
        return
      endif
# endif
      do while(power > 9)              ! Insert time index or MPI node
        power=power/10 ; j=k/power     ! number (rank) into the string,
        i=i+1          ; k=k-j*power   ! then attach suffix, if needed.
        name(i:i)=char(48+j)
      enddo
      if (lsffx > 0) name(i+1:i+lsffx)=sffx(1:lsffx)
      lstr=i+lsffx
#endif  /* EXTRACT_INDEX */
      end

#ifndef TIME_INDEX
# define TIME_INDEX
# include "insert_node.F"
#else
# ifndef EXTRACT_INDEX
#  define EXTRACT_INDEX
#  include "insert_node.F"
# else
#  ifdef TEST_INSERT
      implicit none                       ! Self-testing program
      character(len=64) fname             ! for the subroutines above
      integer lstr, i,iout, nnodes, ierr
      ierr=0
      fname='his_00_*/.*'
      fname='../dir/root_name.000.nc'
      fname='../dir/root_name.*.*.nc'
      fname='../dir/root_name.0'

      fname='../dir/root_name.123.3459.nc'

      nnodes=23

      lstr=1
      do while (lstr < 64 .and. fname(lstr:lstr) /= ' ')
        lstr=lstr+1
      enddo
      if (fname(lstr:lstr) == ' ') lstr=lstr-1

      write(*,'(/1x,A/1x,A)') 'Testing insert_node...',
     &                        '------------------------'
      do i=0,nnodes-1
        call insert_node(fname, lstr, i, nnodes, ierr)
        write(*,'(I4,1x,A,1x,A,I3)') i, fname(1:lstr), 'len =', lstr
      enddo

      write(*,'(/1x,A/1x,A)') 'Testing insert_time_index...',
     &                        '------------------------------'
      do i=0,1000,80
        call insert_time_index (fname, lstr, i, ierr)
        write(*,'(1x,A,I4,2x,3A,I3)') 'tindx =', i, 'fname = ''',
     &                          fname(1:lstr), '''  len =', lstr
      enddo

      write(*,'(/1x,A/1x,A)') 'Testing extract_time_index...',
     &                        '-------------------------------'
      do i=0,1000,75
        call insert_time_index(fname, lstr, i, ierr)
        call extract_time_index(fname, lstr, iout, ierr)
        write(*,'(1x,3A,I8)') 'fname = ''',fname(1:lstr),
     &                              '''  tindx =', iout
      enddo
      end
#  endif
# endif
#endif
