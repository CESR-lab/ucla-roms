module mg_tictoc
  ! intrinsec fortran function

  implicit none

  integer(kind=4), parameter :: st=4, lg=8

  integer(kind=st), parameter :: levmax=32, submax=32

  integer(kind = lg) , dimension(levmax,submax) :: ntic
  integer(kind = lg) , dimension(levmax,submax) :: ntoc
  real(kind = lg)    , dimension(levmax,submax) :: time_tictoc
  real(kind = lg)    , dimension(levmax,submax) :: time_tictoc_sort
  integer(kind=st)   , dimension(levmax,submax) :: calls
  integer(kind=st)   , dimension(levmax,submax) :: calls_sort
  character(len=32)  , dimension(submax)        :: subname
  character(len=32)  , dimension(submax)        :: subname_sort
  integer(kind=st)                              :: nblev = 0
  integer(kind=st)                              :: nbsub = 0

contains

  !------------------------------------------------
  subroutine tic(lev, string)
    integer(kind=st), intent(in) :: lev
    character(len=*), intent(in) :: string

    integer(kind=st) :: ns
    logical :: flag 

    if (nbsub > 0) then

       flag = .true.

       !- Search if subroutine is already timed
       !- if yes -> cpu_time(tic)
       do ns=1, nbsub
          if (TRIM(string) == subname(ns)) then
             !call cpu_time(ntic(lev,ns))
             call system_clock(ntic(lev,ns))
             flag = .false.
             exit
          endif
       end do

       !- New subroutine to time
       !- Add its name to "subname"
       !- cpu_time(tic)
       if (flag) then
          nbsub = nbsub + 1
          subname(nbsub)=TRIM(string)
!          call cpu_time(ntic(lev,nbsub))
          call system_clock(ntic(lev,nbsub))
          time_tictoc(lev,nbsub)  = 0._lg
          calls(lev,nbsub) = 0
       endif

    else
       !- First subroutine to time
       !- add its name to "subname"
       !- cpu_time(tic)
       nbsub = 1
       subname(nbsub)=TRIM(string)
!       call cpu_time(ntic(lev,nbsub))
       call system_clock(ntic(lev,nbsub))

       time_tictoc(lev,nbsub)  = 0._lg
       calls(lev,nbsub) = 0
    endif

    if (lev > nblev) nblev = lev

  end subroutine tic

  !------------------------------------------------
  subroutine toc(lev, string)
    integer(kind=st), intent(in) :: lev
    character(len=*), intent(in) :: string

    integer(kind=st) :: ns
    logical :: flag 

    real(kind=lg)   :: rate

    call system_clock(count_rate=rate)
!    call system_clock(count_max=cm)

    if (nbsub > 0) then

       flag = .true.

       do ns=1, nbsub
          if (TRIM(string) == subname(ns)) then
!             call cpu_time(ntoc(lev,ns))
             call system_clock(ntoc(lev,ns))
             time_tictoc(lev,ns) = time_tictoc(lev,ns) + real(ntoc(lev,ns) - ntic(lev,ns),kind=lg)/rate
             calls(lev,ns) = calls(lev,ns) + 1
             if (lev > nblev) nblev = lev
             flag = .false.
             exit
          endif
       end do

       if (flag) then
          write(*,*)'Error: tictoc: a toc is calling before a tic !'
          write(*,*)'Error: check if a tic exist for:', TRIM(string)
       endif

    else
       write(*,*)'Error: tictoc: a toc is calling before a tic !'
       write(*,*)'Error: check if a tic exist for:', TRIM(string)
    endif

  end subroutine toc

  !------------------------------------------------
  subroutine get_tictoc(val, lev, string)
    real(kind = lg) , intent(out) :: val
    integer(kind=st), intent(in)  :: lev
    character(len=*), intent(in)  :: string

    integer(kind=st) :: ns
    logical :: flag 

    if (nbsub > 0) then

       flag = .true.

       do ns=1, nbsub
          if (TRIM(string) == subname(ns)) then
             val = time_tictoc(lev,ns)
             flag = .false.
             exit
          endif
       end do

       if (flag) then
          write(*,*)'Error: tictoc: no information for this name !'
          write(*,*)'Error: check if a tictoc calls exist for:', TRIM(string)
       endif

    else
       write(*,*)'Error: tictoc: no information for this name !'
       write(*,*)'Error: check if a tictoc calls exist for:', TRIM(string)
    endif

  end subroutine get_tictoc

  !------------------------------------------------
  subroutine print_tictoc(myrank)
    integer(kind=st), optional, intent(in)::myrank

    integer(kind=st)  :: lev  ! level
    integer(kind=st)  :: ii
    integer(kind=st)  :: lun  ! logical unit number
    integer(kind=st)  :: mr   ! my rank
    CHARACTER(len=32) :: filename
    CHARACTER(len=14) :: cmftf, cmfti
    character(len=32) :: ctmp
    real(kind=lg)   , dimension(1:nblev) :: tmp
    integer(kind=st), dimension(1:nblev) :: itmp
    real(kind=lg)     :: sum1, sum2
    logical           :: end_sort

    if (present(myrank)) then
       mr  = myrank
       lun = myrank + 10
    else
       mr  = 9999
       lun =  10
    endif

! Sort results
    time_tictoc_sort(:,:) = time_tictoc(:,:)
    subname_sort(:) = subname(:)
    calls_sort(:,:) = calls(:,:)

    do                
       end_sort = .true.

       do ii=2,nbsub

          sum1=sum(time_tictoc_sort(1:nblev,ii  ))
          sum2=sum(time_tictoc_sort(1:nblev,ii-1))

          if (sum1 > sum2) then

             end_sort = .false.

             tmp(:) = time_tictoc_sort(1:nblev,ii-1)
             time_tictoc_sort(1:nblev,ii-1) = time_tictoc_sort(1:nblev,ii)
             time_tictoc_sort(1:nblev,ii) = tmp(:)

             ctmp   = trim(subname_sort(ii-1))
             subname_sort(ii-1) = trim(subname_sort(ii))
             subname_sort(ii) = trim(ctmp)

             itmp(:) = calls_sort(1:nblev,ii-1)
             calls_sort(1:nblev,ii-1) = calls_sort(1:nblev,ii)
             calls_sort(1:nblev,ii) = itmp(:)

          end if

       end do

       if (end_sort) exit
    end do

 write(filename, '(A14,I4.4,A4)')'nhmg_prof_rank',mr,'.log'

 open(unit=lun,file=trim(filename),form='FORMATTED')

 WRITE(cmftf , 1000) nblev
 WRITE(cmfti , 1001) nblev
1000 FORMAT('(', I3, '(x,E9.3))')
1001 FORMAT('(', I3, '(x,I9))')

 write(lun,'(t22)', ADVANCE="no")
 write(lun,'(A10)', ADVANCE="no") 'Total'
 do lev=1, nblev
    write(lun,'(x,I9)', ADVANCE="no") lev
 enddo

 write(lun,'(x)', ADVANCE="yes")

 do ii=1, nbsub
    write(lun,'(x,A20)' , ADVANCE="no" ) TRIM(subname_sort(ii))
    write(lun,'(x,E9.3)', ADVANCE="no" ) sum(time_tictoc_sort(1:nblev,ii))
    write(lun,FMT=cmftf , ADVANCE="no" ) time_tictoc_sort(1:nblev,ii)
    write(lun,'(x)'     , ADVANCE="yes")
    write(lun,'(t22)'   , ADVANCE="no" )
    write(lun,'(x,I9)'  , ADVANCE="no" ) sum(calls_sort(1:nblev,ii))
    write(lun,FMT=cmfti , ADVANCE="no" ) calls_sort(1:nblev,ii)
    write(lun,'(x)'     , ADVANCE="yes")
 end do

!!$write(lun,*)""
!!$ do ii=1, nbsub
!!$    write(lun,'(x,A20)' , ADVANCE="no" ) TRIM(subname(ii))
!!$    write(lun,'(x,E9.3)', ADVANCE="no" ) sum(time_tictoc(1:nblev,ii))
!!$    write(lun,FMT=cmftf , ADVANCE="no" ) time_tictoc(1:nblev,ii)
!!$    write(lun,'(x)'     , ADVANCE="yes")
!!$    write(lun,'(t22)'   , ADVANCE="no" )
!!$    write(lun,'(x,I9)'  , ADVANCE="no" ) sum(calls(1:nblev,ii))
!!$    write(lun,FMT=cmfti , ADVANCE="no" ) calls(1:nblev,ii)
!!$    write(lun,'(x)'     , ADVANCE="yes")
!!$ end do

 close(lun)

end subroutine print_tictoc

  !------------------------------------------------
  subroutine print_tictoc_old(myrank)
    integer(kind=st), optional, intent(in)::myrank

    integer(kind=st)  :: lev
    integer(kind=st)  :: ii
    integer(kind=st)  :: lun
    CHARACTER(len=14) :: cmftf, cmfti

    if (present(myrank)) then
       lun = myrank + 10
    else
       lun = 10
    endif

    WRITE(cmftf , 1000) nblev
    WRITE(cmfti , 1001) nblev

1000 FORMAT('(', I3, '(x,E9.3))')
1001 FORMAT('(', I3, '(x,I9))')

    write(lun,'(t22)', ADVANCE="no")
    write(lun,'(A10)', ADVANCE="no") 'Total'
    do lev=1, nblev
       write(lun,'(x,I9)', ADVANCE="no") lev
    enddo

    write(lun,'(x)', ADVANCE="yes")

    do ii=1, nbsub
       write(lun,'(x,A20)' , ADVANCE="no" ) TRIM(subname(ii))
       write(lun,'(x,E9.3)', ADVANCE="no" ) sum(time_tictoc(1:nblev,ii))
       write(lun,FMT=cmftf , ADVANCE="no" ) time_tictoc(1:nblev,ii)
       write(lun,'(x)'     , ADVANCE="yes")
       write(lun,'(t22)'   , ADVANCE="no" )
       write(lun,'(x,I9)'  , ADVANCE="no" ) sum(calls(1:nblev,ii))
       write(lun,FMT=cmfti , ADVANCE="no" ) calls(1:nblev,ii)
       write(lun,'(x)'     , ADVANCE="yes")
    end do

  end subroutine print_tictoc_old

end module mg_tictoc
