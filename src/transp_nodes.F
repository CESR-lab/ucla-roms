      program transp_nodes

! A program to automatically generate MPI machine file machines.LINUX
! for the optimal placement of MPI ranks to the actual hardware nodes
! (hosts) on a cluster made of multi-processor nodes. In takes an input
! file (usually  $PBS_NODEFILE or $LSB_HOSTS) which contains the list
! of hostnames of compute nodes allocated to the batch job submitted
! into the queue (obviously there is no way to know the hostnames until
! the job actually starts) and rearranges it such a way that MPI-nodes
! adjacent in ETA-direction are placed on the same multi-processor
! compute node.  The goal is to minimize sending MPI messages across
! the network interconnect (Infiniband or whatever) by directing as
! much communication as possible to be internal within each multi-
! processor compute node.

! For example, optimal placement for 6 quad-processor nodes (hence
! 24 CPUs total) configuration is with 3x8 sub-domain decomposition is

!   ==============     ==============     ==============     0 HOST_1
!  |              |   |              |   |              |    1 HOST_2
!  |   rank 21    |   |   rank 22    |   |   rank 23    |    2 HOST_3
!  |              |   |              |   |              |    3 HOST_1
!   --------------     --------------     --------------     4 HOST_2
!  |              |   |              |   |              |    5 HOST_3
!  |   rank 18    |   |   rank 19    |   |   rank 20    |    6 HOST_1
!  |              |   |              |   |              |    7 HOST_2
!   --- HOST_4 ---     --- HOST_5 ---     --- HOST_6 ---     8 HOST_3
!  |              |   |              |   |              |    9 HOST_1
!  |   rank 15    |   |   rank 16    |   |   rank 17    |   10 HOST_2
!  |              |   |              |   |              |   11 HOST_3
!   --------------     --------------     --------------    12 HOST_4
!  |              |   |              |   |              |   13 HOST_5
!  |   rank 12    |   |   rank 13    |   |   rank 14    |   14 HOST_6
!  |              |   |              |   |              |   15 HOST_4
!   ==============     ==============     ==============    16 HOST_5
!                                                           17 HOST_6
!   ==============     ==============     ==============    18 HOST_4
!  |              |   |              |   |              |   19 HOST_5
!  |    rank 9    |   |    rank 10   |   |    rank 11   |   20 HOST_6
!  |              |   |              |   |              |   21 HOST_4
!   --------------     --------------     ---------- ---    22 HOST_5
!  |              |   |              |   |              |   23 HOST_6
!  |    rank 6    |   |    rank 7    |   |    rank 8    |
!  |              |   |              |   |              |
!   --- HOST_1 ---     --- HOST_2 ---     --- HOST_3 ----
!  |              |   |              |   |              |
!  |    rank 3    |   |    rank 4    |   |    rank 5    |
!  |              |   |              |   |              |
!   ---------- ---     ---------- ---     --------------
!  |              |   |              |   |              |
!  |    rank 0    |   |    rank 1    |   |    rank 2    |
!  |              |   |              |   |              |
!   ==============     ==============     ==============

! where HOST_1, ..., HOST_6 are names of quad-processor compute hosts,
! while rank 0, ..., rank 23 are MPI-ranks working on individual MPI
! subdomains.  The table on the right shows the placement of each MPI
! rank (left column) to the physical node (right column, which is also
! the actual output file "machines.LINUX" of this program).
!
! This program is expected to be used in batch script submitted to
! queue as follows:
!
!        transp_nodes $PBS_NODEFILE
!        mpirun_rsh -np 256 -hostfile machines.LINUX roms_exe
! or
!        echo $LSB_HOSTS > lsb_hosts
!        transp_nodes lsb_hosts
!        mpiexec -np N -machinefile machines.LINUX roms_exe
!
! depending on the specific queueing system: some like PBS may provide
! a file which already looks like machines.LINUX (listing one hostname
! per line); others like LSB provide just an environmental variable --
! a long character string containing all the node names allocated to
! the current batch job.  Either way "transp_nodes" takes the list
! provided by the queuing system and creates file "machines.LINUX"
! specifically suitable for grid partitioning geometry NP_XI,NP_ETA
! (thus "transp_nodes" is expected to be compiled together with roms
! executable -- it needs the same "param.h"). The number of CPUs
! (cores) per host is determined from the number of repetitions of
! each individual hostname in the list


#define VERBOSE
#define VISUALIZE

#include "cppdefs.opt"
#ifdef MPI

      use param
      implicit none
      integer, parameter :: iin=11, iout=12
      character(len=14), parameter :: machinefile = 'machines.LINUX'
      character(len=10) fmt
      character(len=64) fname, bffr
      character(len=:), allocatable :: str, host(:)
      integer,          allocatable :: lhost(:), jhost(:), jrpts(:)
      integer nhosts, nunique, npnode, max_str_size, max_name_size,
     &        is,ie, ierr, inode, jnode,  n_xi_rpt, n_eta_rpt, jinc,
     &                        i,j, ii,jj, lstr, lbfr, lfnm, lenstr
      logical lswitch

! Function "iargc" is viewed as intrinsic by most modern compilers
! and does not need to be declared. IBM xlf95 is a notable exclusion.
! So do 7.x and earlier versions of Intel IFC compiler, but starting
! with version 8.x IFORT recognize it as intrinsic.

#ifdef INT_IARGC
      integer iargc
#endif

! Open hostname file, inquire about its size in bytes, allocate
! a character string of sufficient length to accommodate the whole
! file and try to read it. Because different type queueing system
! may have supply it in two different formats -- (i) either a long
! character string listing all the nodes, or (ii) a file listing
! one node per line -- the code below is designed to handle both
! contingencies automatically: at first it attempts to read the file
! as would be a single string, checks whether the string contains
! blanc characters separating the names: if not keep reading the
! file line-by-line and add new names into the string -- the result
! of this operation will be essentially the same for both cases --
! a long string.

      nhosts=0
      if (iargc() == 1) then
        call getarg(1,fname)
        lfnm=lenstr(fname)
        open(unit=iin, file=fname(1:lfnm), status='old', iostat=ierr)
        if (ierr == 0) then
          inquire(unit=iin, size=max_str_size)
          allocate(character(len=max_str_size) :: str)
          read(iin,'(A)',iostat=ierr) str
          if (ierr == 0) then
            ie=lenstr(str)
            nhosts=1
            do i=2,ie-1
              if (str(i:i) == ' ') nhosts=nhosts+1
            enddo
            if (nhosts == 1) then
              do while(ierr == 0)
                read(iin,'(A)',iostat=ierr) bffr
                if (ierr == 0) then
                  lbfr=lenstr(bffr)
                  is=ie+2 ; ie=is+lbfr-1
                  if (ie <= max_str_size) then
                    str(is-1:is-1)=' '
                    str(is:ie)=bffr(1:lbfr)
                    nhosts=nhosts+1
                  else
                    write(*,'(/2A/)') '### ERROR: transp_nodes :: ',
     &                      'insufficient length of buffer string.'
                    stop
                  endif
                endif
              enddo
            endif
          else
            write(*,'(/4A/)')    '### ERROR: transp_nodes :: Cannot ',
     &                  'read host name file ''', fname(1:lfnm), '''.'
            stop
          endif
          close(unit=iin)
        else
          write(*,'(/4A/)') '### ERROR: transp_nodes :: Cannot open ',
     &                       'host name file ''', fname(1:lfnm), '''.'
          stop
        endif
        lstr=lenstr(str)
# ifdef VERBOSE
        write(*,'(3A,I7,1x,A,I5)') 'processing = ''', fname(1:lfnm),
     &                        ''' size =', lstr, 'nhosts =', nhosts
        write(*,'(/3A/)')   'hostnames = ''',  str(1:lstr), ''''
# endif
      else
        write(*,'(/2A/)') '### ERROR: transp_nodes requires ',
     &            'a single argument -- filename is missing.'
        stop
      endif                             ! Preprocessing: find how many
                                        ! hostnames (words separated by
      nhosts=0                          ! blank characters) are in the
      max_name_size=0                   ! string and find the maximum
      lswitch=.true.                    ! length among all words; use
      do i=1,lstr                       ! these two numbers to allocate
        if (lswitch) then               ! storage with appropriate size
          if (str(i:i) /= ' ') then     ! for further processing...
            is=i
            lswitch=.false.
          endif
        else
          if (str(i:i) == ' ' .or. i == lstr) then
            if (i == lstr) then
              ie=lstr
            else
              ie=i-1
            endif
            lswitch=.true.
            nhosts=nhosts+1
            max_name_size=max(max_name_size, ie-is+1)
          endif
        endif
      enddo

      allocate(character(len=max_name_size) :: host(nhosts))
      allocate(lhost(nhosts),  jhost(nhosts), jrpts(nhosts))

      lswitch=.true.                        ! Copy individual hostnames
      j=0                                   ! into indexed character
      do i=1,lstr                           ! array "host" and record
        if (lswitch) then                   ! length of each hostname
          if (str(i:i) /= ' ') then
            is=i
            lswitch=.false.
          endif
        else
          if (str(i:i) == ' ' .or. i == lstr) then
            if (i == lstr) then
              ie=lstr
            else
              ie=i-1
            endif
            lswitch=.true.
            j=j+1
            host(j)=str(is:ie)
            lhost(j)=ie-is+1
          endif
        endif
      enddo
      do j=1,nhosts
        do i=lhost(j)+1,max_name_size
         host(j)(i:i)=' '
        enddo
      enddo         !--> discard str

#define str illegal
                                            ! Post-processing: find
      nunique=1                             ! unique names among all
      jhost(1)=1                            ! the hostnames and
      do i=1,nhosts                         ! record index of the
        lswitch=.true.                      ! first appearance...
        do j=1,nunique
          if ( lhost(i) == lhost(jhost(j)) .and.
     &          host(i) ==  host(jhost(j)) ) then
            lswitch=.false.
          endif
        enddo
        if (lswitch) then
          nunique=nunique+1
          jhost(nunique)=i
        endif
      enddo
                                             ! determine how many
      do j=1,nunique                         ! times each unique
        jrpts(j)=0                           ! hostname is repeated
        do i=1,nhosts
          if ( lhost(i) == lhost(jhost(j)) .and.
     &          host(i) ==  host(jhost(j)) ) then
            jrpts(j)=jrpts(j)+1
          endif
        enddo                                ! verify that that all
      enddo                                  ! unique names appear
                                             ! exactly the same number
      lswitch=.true.                         ! of times
      do j=2,nunique
        if (jrpts(j)/=jrpts(1)) lswitch=.false.
      enddo

# ifdef VERBOSE
      write(*,'(1x,2A,I5/23x,A,I5)')      'SUMMARY: Total number ',
     &                      'of hostnames in the string:', nhosts,
     &                      'number of unique hostnames:', nunique
      if (lswitch) then
        write(*,'(16x,2A,I5,1x,A)')  'All the hostnames repeated ',
     &                               'equally', jrpts(1), 'times.'
      endif
      write(*,'(1x,2A)') '---------------------------------',
     &                   '---------------------------------'
      do j=1,nunique
        write(*,'(4x,3A,I5,1x,A,I4,1x,A)')          'hostname = ''',
     &    host(jhost(j))(1:lhost(jhost(j))), ''' first occurrence =',
     &                     jhost(j), ' repeated', jrpts(j), 'times'
      enddo
      write(*,'(1x,2A)') '---------------------------------',
     &                   '---------------------------------'
# endif

                                      ! At this stage the processing
      ierr=0                          ! of hostnames in LSB string is
      if (lswitch) then               ! complete. Perform various
        npnode=jrpts(1)               ! checks and create machinefile..
      else
        write(*,'(1x,4A)') '### ERROR :: Cannot determine number of ',
     &     'CPUs per node (cores per node) from ', fname(1:lfnm), '.'
        ierr=ierr+1
      endif
      if (nunique*npnode < NP_XI*NP_ETA) then
        write(*,'(1x,4A,I5,1x,A,I5,A)')    '### ERROR :: Not enough ',
     &             'hostnames in file ''', fname(1:lfnm), ''': found',
     &                   nunique*npnode, 'needed ', NP_XI*NP_ETA, '.'
        ierr=ierr+1
      endif

      if (NP_ETA < npnode) then
        n_xi_rpt=npnode/NP_ETA
        if (n_xi_rpt*NP_ETA /= npnode) then
          write(*,'(1x,2A,I5,2A,I5)')   '### ERROR :: Number of CPUs ',
     &            '(cores) per node,',  npnode,  ', cannot be evenly ',
     &                                 'divided by  NP_ETA =', NP_ETA
          ierr=ierr+1
        endif
        n_eta_rpt=1
      else
        n_xi_rpt=1
        n_eta_rpt=NP_ETA/npnode
        if (n_eta_rpt*npnode /= NP_ETA) then
          write(*,'(1x,A,I5,1x,2A,I5)')      '### ERROR :: NP_ETA = ',
     &           NP_ETA, 'cannot be evenly divided by the number of ',
     &                               'CPUs (cores) per node', npnode
          ierr=ierr+1
        endif
      endif

      jinc=nunique/n_eta_rpt
      if (jinc*n_eta_rpt /= nunique) then
        write(*,'(1x,2A,I4,A,I4,A)')  '### ERROR :: Number of unique ',
     &  'hostnames, ', nunique, ', is not divisible by', n_eta_rpt, '.'
        ierr=ierr+1
      endif

      if (ierr == 0) then
        open(unit=iout, file=machinefile, form='formatted')
        do jnode=0,NP_ETA-1
          jj=jnode/npnode
          do inode=0,NP_XI-1
            ii=inode/n_xi_rpt
            j=1+ii  +jj*jinc
c**         write(*,*) j, jhost(j)
            write(iout,'(A)') host(jhost(j))(1:lhost(jhost(j)))
          enddo
        enddo                         ! At this stage machinefile is
        close(unit=iout)              ! complete; the following part
                                      ! is merely for visual checking:
# ifdef VISUALIZE
        open(unit=iin, file=machinefile, form='formatted',status='old')
        do i=1,NP_XI*NP_ETA
          read(iin,'(A)') host(i)     ! read back the hostnames from
          lhost(i)=lenstr(host(i))    ! machinefile and print them out
        enddo                         ! in format of a 2D NP_XI*NP_ETA
        close(unit=iin)               ! matrix to visualize mapping.

        write(*,'(1x,A,I4,2(3x,A,I4))') 'NP_XI =', NP_XI, 'NP_ETA =',
     &           NP_ETA, 'Number of CPUs (cores) per node =', npnode
        write(*,'(1x,2A)') '---------------------------------',
     &                     '---------------------------------'
        write(fmt,'(A,I2,A)'), '(', NP_XI, '(1x,A))'
        do jnode=0,NP_ETA-1
          jj=jnode*NP_XI
          write(*,fmt) (host(i+jj)(1:lhost(i+jj)), i=1,NP_XI)
        enddo
        write(*,'(1x,2A)') '---------------------------------',
     &                     '---------------------------------'
# endif
      endif
      end
# undef str
# include "lenstr.F"
#else
      end
#endif   /* MPI */
