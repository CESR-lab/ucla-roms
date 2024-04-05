#ifndef EXCLUDE_MAIN
# define EXCLUDE_MAIN

      program ncjoin

      ! Creator: Sasha (Alex S.)
      ! Edits:
      ! - 2021-Jan: DevinD changed nf_create and nf_def_var to nf90
      !             routines to create netcdf-4 files to allow for
      !             file compression.
      !             Change compression level with 'deflate_level='.

! Generic netCDF assembly tool: reads ROMS partial netCDF files and
! assembles them into a file for the whole physical grid.  It is the
! inverse operation to partit.  Usage:
!
!                  ncjoin files.???.nc 
! or
!                  ncjoin -d files.???.nc
! or
!                  ncjoin --delete files.???.nc
!
! where files.???.nc matches a complete set of partial files (or
! several complete sets) and "-d" or "--delete" forces deletion of
! partial files upon successful joining. 


! This code uses self-inclusion in order to generate several similar
! subroutines from the same source code. CPP-macro EXCLUDE_MAIN makes
! main program appear only once after CPP is executed: it is assumed 
! to be undefined initially.  All other CPP-switches below are user
! configurable. 

! CPP-switch "FAST_MODE" activates mode in which partial arrays are
! assembled into a big array in operating memory, and the whole data
! is written into the disk as a one long record.  If not defined, the
! code relies on the netCDF capability to write rectangular sub-arrays
! via nf_put_vara_XXX routines.  Because in hardware terms the latter
! results in writing multiple small records (with record length equal
! to the first fortran dimension of sub-array), it is considerably
! slower.

! CPP-switch "DEL_PART_FILES" enables the user-activated option to
! delete partial files just after they are joined.  If compiled with
! this switch is defined, ncjoin becomes sensitive to whether its first
! argument is "-d" or "--delete".  The deletion occurs one-set-at-at-
! a-time and ONLY if the joining is successful, i.e., if an error of
! any kind occurs, so the set cannot be joined for whatever reason (a
! missing file, or one file has fewer records than others, or reading
! error; etc), the files stay intact.  If ncjoin is interrupted, the
! partial files stay intact.
 
! CPP-switch "KEEP_CLOSED" forces regime in which to keep no more
! than one partial file for MPI-nodes 1:nnnodes-1 be opened at any
! given time (file belonging to MPI-node=0 is kept permanently open,
! so does the target file). This regime is needed to avoid netCDF
! internal memory limitation associated with opening too many files
! simultaneously.

! CPP-switch "AUTORENICE" forces automatic lowering priority of ncjoin
! to 19 (the lowest possible) to run it in "background" regime.



#define FAST_MODE
#define DEL_PART_FILES
c--#define KEEP_CLOSED
#define TIMING
c--#define AUTORENICE
c--#define VERBOSE

      use netcdf ! 2021 addition for nf90_create & nf90_def_var

      implicit none
#include "netcdf.inc"

      ! START USER INPUT
      integer :: deflate_level=1 ! Compression level for joined file
      logical :: shuffle=.true.  ! Shuffle on for extra conmpression
      ! END USER INPUT
      
      integer, dimension(4) :: chunk_sizes ! DevinD: chunksize debugging

      integer, parameter :: maxdims=32, maxvars=220 ! BGC has 201 vars in restart file
      integer, parameter :: part_XI=1,  part_ETA=2, part_2D=3

      logical complete, clean_set,  digit, var_mask,    lnewvar
      integer nargs, nnodes,  size_XI,  XI_rho, id_xi_rho,  id_xi_u,
     &        arg,   node,    size_ETA, ETA_rho,id_eta_rho, id_eta_v,
     &        ierr,  maxnodes,size_S,   tsize,  unlimdimid, rec,
     &        ntest, nctarg,  ndims,    size,   code_size,  lvar,
     &        nvars, ngatts,  varatts,  size1,  code_size_bak,
     &        i,j,k, is,ie,   lncn,     ltrg,   lstr, lbak, lenstr

      character(len=8) sffx, sffx_bak
      character(len=32) vname(maxvars), dimname(maxdims) 
      character(len=128) :: nctestname, nctargname, root, root_bak, string
      character(len=128), dimension(:), allocatable :: ncname
      character(len=1) :: buff_char ! DevinD added as Expanse complained about putting
                                    ! spherical character in numerical 'buff' array...

      integer, dimension(:), allocatable:: ncid, xi_start, eta_start
      integer, dimension(:,:), allocatable :: vid, dimsize
      logical, dimension(:), allocatable :: western_edge,
     &         eastern_edge, southern_edge, northern_edge

      logical series(maxvars)
      integer, dimension(maxvars) :: varid,   vnode,  vdims,
     &                               vartype, part_type
      integer, dimension(maxdims) :: dimid,   ldim,   ibuff,
     &                               start,   count,  start1
      integer, dimension(maxdims,maxvars) :: dimids
      integer*8 max_buff_size, alloc_buff_size
      real*8, allocatable, dimension(:) :: buff
#ifdef FAST_MODE
      integer count1(maxdims)
      integer*8 max_bfr_out, alloc_bfr_out
      real*8, allocatable, dimension(:) :: bfr_out
#endif
#ifdef DEL_PART_FILES
      logical del_part_files
      character(len=128) rmcmd
#endif
#ifdef TIMING
      real*4 tstart, RUN_time, CPU_time(2)
      integer iclk(2), nclk, clk_rate, clk_max, iclk_init
      integer*8 net_read_size, net_wrt_size, net_fcrt_clk,
     &          net_read_clk,  net_wrt_clk,  net_assm_clk,
     &          net_sync_clk,  net_gray_clk, inc_clk
      real*8 ReadSize, ReadTime, WrtSize,  WrtTime, 
     &       FcrtTime, AssmTime, SyncTime, GrayTime
# ifdef __IFC
      real*4 etime
# endif
# ifdef DEL_PART_FILES
      integer*8 net_rmcmd_clk
# endif
#endif             /* ! Function "iargc" is viewed as intrinsic by  */
#ifdef AIX         /* ! most modern compilers and does not need to   */
      integer iargc   ! be declared. IBM xlf95 is a notable exclusion.
#endif             /* ! So do 7.x and earlier versions of Intel IFC; */
                      ! 8.x and later IFORT recognize it as intrinsic.
#ifdef AUTORENICE
      integer getpid, pid              ! Sometimes it makes sense to
      character(len=32) cmd            ! run ncjoin in the background
      pid=getpid()                     ! mode with lowered priority to
      write(cmd,'(I8)') pid            ! not interfere with running MPI
      lstr=lenstr(cmd)                 ! job. This code segment catches
      cmd(11:lstr+10)=cmd(1:lstr)      ! pid of its own process and
      lstr=lstr+10                     ! executes re-nice command.
      cmd(1:10)='renice 19 '
      write(*,'(/3A/)') 'Autorenice: executing ''',cmd(1:lstr),'''.'
      call system(cmd(1:lstr))
      write(*,*)
#endif
#ifdef TIMING
# ifdef __IFC
      tstart = etime(CPU_time)
# else
      call etime(CPU_time, tstart)
# endif
      nclk=1
      call system_clock (iclk(nclk), clk_rate, clk_max)
      iclk_init=iclk(nclk)
      net_read_clk=0
      net_read_size=0
      net_wrt_size=0
      net_wrt_clk=0
      net_sync_clk=0
      net_assm_clk=0
      net_gray_clk=0
      net_fcrt_clk=0 ! DevinD added this was missing and caused seg fault.
#endif
#ifdef DEL_PART_FILES
      del_part_files=.false.
# ifdef TIMING
      net_rmcmd_clk=0
# endif
#endif
      ntest=-1                         ! initialize sizes of buffer 
      maxnodes=-1                      ! arrays to be allocated. Here
      max_buff_size=0                  ! "max_*" means the needed size, 
      alloc_buff_size=0                ! and "alloc_*" is the size of 
#ifdef FAST_MODE
      max_bfr_out=0                    ! the actually allocated array.
      alloc_bfr_out=0
#endif
      nargs=iargc()
      arg=0


! Extract a set of files which cover the whole physical grid.

  1   nnodes=-1                        ! Reset variables which
      root_bak(1:1)=' '                ! define the file set.
      sffx_bak(1:1)=' '
      code_size_bak=-1
 
      ! DD: will cycle through argument list of names of partitioned
      !     files by coming back here with 'goto 2'
  2   arg=arg+1
        ! DD: get the arguements from calling ncjoin arg1 ...
        !     e.g. ncjoin his.0000.*.nc will actually send in
        !     ncjoin his.0000.00.nc his.0000.01.nc ... and so on.
        call getarg(arg,nctestname)
        lncn=lenstr(nctestname)
#ifdef DEL_PART_FILES
        if (arg.eq.1 .and. (lncn.eq.2 .and. nctestname(1:2).eq.'-d'
     &     .or. lncn.eq.8 .and. nctestname(1:8).eq.'--delete') ) then
           write(*,'(/1x,2A/)') '>>>> Flag to delete partial files ',
     &                                                  'is raised.'
           del_part_files=.true.
          goto 2
        endif
#endif
        if (ntest.ne.-1) then
          ierr=nf_close(ntest)
          ntest=-1
        endif
        ierr=nf_open (nctestname, nf_nowrite, ntest)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_att (ntest, nf_global, 'partition', i, lvar)
          ! DD: e.g. :partition = 0, 6, 1, 1 ;
          ! i=attribute type, lvar= number of values stored in attribute.
          if (ierr .eq. nf_noerr) then
            if (i.eq.nf_int .and. lvar.eq.4) then
              ierr=nf_get_att_int (ntest,nf_global, 'partition', ibuff)
              ! DD: ibuff array stores 4 integers of 'partition'.
              if (ierr .eq. nf_noerr) then

                if (nnodes.eq.-1) then
                  nnodes=ibuff(2)
                  if (nnodes.gt.maxnodes) then
                    maxnodes=nnodes
                    if (allocated(ncid)) then
                      deallocate(dimsize)
                      deallocate(vid)
                      deallocate(northern_edge)
                      deallocate(southern_edge)
                      deallocate(eastern_edge)
                      deallocate(western_edge)
                      deallocate(eta_start)
                      deallocate(xi_start)
                      deallocate(ncname)
                      deallocate(ncid)
                    endif
                    allocate (ncid(0:nnodes-1))
                    allocate (ncname(0:nnodes-1))
                    allocate (xi_start(0:nnodes-1))
                    allocate (eta_start(0:nnodes-1))
                    allocate (western_edge(0:nnodes-1))
                    allocate (eastern_edge(0:nnodes-1))
                    allocate (southern_edge(0:nnodes-1))
                    allocate (northern_edge(0:nnodes-1))
                    allocate (vid(maxvars,0:nnodes-1))
                    allocate (dimsize(maxdims,0:nnodes))
                  endif
                                          ! Reset variables defining
                  complete=.false.        ! a complete set of partitial
                  do node=0,nnodes-1      ! files. These variables will
                    ncid(node)=-1         ! receive meaningful values
                    xi_start(node)=-1     ! from data read from netCDF
                    eta_start(node)=-1    ! file headers, subsequently
                  enddo                   ! be used to verify the set
                                          ! completeness.

                elseif (nnodes.ne.ibuff(2)) then
                  write(*,'(/1x,2A,I4/14x,3A/14x,A,I4,4x,A/)')
     &                 '### WARNING: Number of MPI nodes in global ',
     &                 'attribute ''partition'', nnodes =', ibuff(2),
     &                 'in netCDF file ''',       nctestname(1:lncn),
     &                 ''' contradicts that from the initial',
     &                 'file in the sequence, nnodes =',      nnodes,
     &                                   ' ==> The file is ignored.'
                  arg=arg-1
                  goto 5
                endif

                node=ibuff(1)
                if (ncid(node).ne.-1) then 
                  write(*,'(/1x,2A,I4,1x,A)') '### ERROR: netCDF ID ',
     &                         'for file corresponding to MPI-node =',
     &                                   node,  'is already in use.'
                  stop
                endif


                if (ncid(node).eq.-1 .and. xi_start(node).eq.-1
     &                         .and.  eta_start(node).eq.-1) then
                  ncid(node)=ntest 
                  ncname(node)=nctestname
                  xi_start(node)=ibuff(3)
                  eta_start(node)=ibuff(4)

#define ntest illegal


! Lexical analysis of the file name: It is assumed that name of the
! file consists of root name (eg.: "history"); integer number which
! contains MPI node number (eg.: 03) and; suffix (eg.: ".nc").
! Files which belong to the same set normally have the same (1) root
! and (2) suffix names; the same (3) number of digits in the MPI node
! number segment in the filename and; (4) MPI node number from the
! file name should match the number determined from global attribute
! "partition".
                                               ! Determine positions
                  digit=.false.                ! of starting and ending
                  is=0                         ! characters of MPI node
                  ie=0                         ! segment (is:ie)
                  i=lncn+1
                  do while (is.eq.0 .and. i.gt.1)
                    i=i-1
                    if (nctestname(i:i).ge.'0' .and.
     &                  nctestname(i:i).le.'9') then 
                      if (.not.digit) then
                        if (i.lt.lncn) then
                          if (nctestname(i+1:i+1).eq.'.') then
                            ie=i
                            digit=.true.         ! check that node
                          endif                  ! segment and suffix
                        else                     ! are separated by '.'
                          ie=i                   ! no-suffix case
                          digit=.true.
                        endif
                      endif
                    elseif (digit .and. nctestname(i:i).eq.'.') then
                      digit=.false.
                      is=i+1 
                    endif
                  enddo

                  if (is.gt.0 .and. ie.ge.is) then
                    root=nctestname(1:is-1)
                    if (ie.lt.lncn) then         ! Extract common
                      sffx=nctestname(ie+1:lncn) ! part of file names,
                    else                         ! MPI node number
                      sffx(1:1)=' '              ! and suffix (if any)
                    endif
                    k=0
                    do i=is,ie
                      k=10*k + ichar(nctestname(i:i))-48
                    enddo
                    code_size=ie-is+1
                  else
                    write(*,'(/1x,3A/)')        '### ERROR: Cannot ',
     &                  'determine MPI node number from filename ''',
     &                                     nctestname(1:lncn), '''.'
                  endif
# ifdef VERBOSE
                  write(*,'(1x,3A,I3,1x,A,I4,1x,A,I4,3x,A,2I4)')
     &           'fname = ''', nctestname(1:lncn),  ''' code_size =',
     &            code_size,  'code =', k, 'node =', node, 'i,jSW =',
     &                              xi_start(node), eta_start(node)
# endif

! Checking consistency of root name with previously found.

                  ierr=nf_noerr
                  if (root_bak(1:1).eq.' ') then
                    root_bak=root
                  else
                    lvar=lenstr(root)
                    lbak=lenstr(root_bak)
                    if (lvar.ne.lbak .or. root.ne.root_bak) then
                      ierr=ierr+1
                      write(*,'(/8x,6A/17x,3A/)') 'WARNING: file ''',
     &                     nctestname(1:lncn),   ''' has different ',
     &                    'root name ''',  root(1:lvar),   ''' than',
     &                    'previously found root name ''',
     &                     root_bak(1:lbak), ''' from the same set.'
                    endif
                  endif

! Checking consistency of suffix with previously found..

                  if (sffx_bak(1:1).eq.' ') then
                    sffx_bak=sffx
                  else
                    lvar=lenstr(sffx)
                    lbak=lenstr(sffx_bak)
                    if (lvar.ne.lbak .or. sffx.ne.sffx_bak) then
                      ierr=ierr+1
                      write(*,'(/8x,7A/17x,3A/)')       'WARNING: ',
     &                  'file ''',  nctestname(1:lncn),   ''' has ',
     &                  'different suffix name ''',    sffx(1:lvar),
     &                  ''' than','previously found suffix name ''',
     &                  sffx_bak(1:lbak),   ''' from the same set.'
                    endif
                  endif

! Checking consistency of length of node number segment

                  if (code_size_bak.eq.-1) then
                    code_size_bak=code_size
                  elseif (code_size .ne. code_size_bak) then
                    ierr=ierr+1
                    write(*,'(/8x,A,I2,1x,A/17x,3A,I2,A/)')
     &              'WARNING: number of digits in MPI node segment',
     &               code_size, 'in filename', '''',
     &               nctestname(1:lncn),
     &                ''' is different than previously determined',
     &                                     code_size_bak, '.'
                  endif
 
! Checking consistency of node number with the file name.

                  if (k.ne.node) then
                    ierr=ierr+1
                    write(*,'(/8x,3A,I3/17x,2A/17x,A,I3,A/)')
     &                   'WARNING: file ''', nctestname(1:lncn),
     &                   ''' belongs to different MPI node',   node,
     &                   '(as determined from its global attribute',
     &                   '''partition'')', 'than node', k,
     &                   ' determined from to the file name.'
                  endif

! Stop, if something is wrong.

                  if (ierr.ne.nf_noerr) goto 97
                else
                  arg=arg-1
                  goto 5
                endif
              else
                write(*,'(/1x,2A/14x,3A/)')   '### WARNING: Cannot ',
     &           'aquire global attribute ''partition'' from netCDF',
     &                                 'file ''', nctestname(1:lncn),
     &                               '''. ==> This file is ignored.'
              endif
            else
              write(*,'(/1x,2A/14x,3A/)')  '### WARNING: Wrong type ',
     &                'or size of global attribute ''partition'' in ',
     &                           'netCDF file ''', nctestname(1:lncn),
     &                                '''. ==> This file is ignored.'
            endif
          else
            write(*,'(/1x,3A/)') '### WARNING: ''', nctestname(1:lncn),
     &      ''' is not a partial netCDF file: ==> The file is ignored.'
          endif
        else
          write(*,'(/1x,4A/14x,A/)')   '### WARNING: Cannot open ''',
     &                   nctestname(1:lncn), ''' as a netCDF file: ',
     &                nf_strerror(ierr), ' ==> The file is ignored.'
        endif

#define nctestname illegal

! Verify, whether ncname(0:nnodes-1) and ncid(0:nnodes-1) > 0 (i.e.,
! successfully opened for reading) comprise a complete set of partial
! files.  Keep searching, if not.

   5    continue
        if (nnodes.gt.0) then
          complete=.true.
          do node=0,nnodes-1
            if (ncid(node).lt.0) complete=.false.
          enddo
        endif
      if (.not.complete .and. arg.lt.nargs) goto 2  !--> next file

#ifdef VERBOSE
      write(*,*) ' line 433, complete =', complete, ' nnodes =', nnodes 
#endif

! Once a complete set is identified, print the finenames. 

      if (complete) then
        lncn=lenstr(ncname(0))
        write(*,'(2(1x,A,I4),1x,A,2x,A,2I5)')   'Processing set of ',
     &                         nnodes, 'files', 0, ncname(0)(1:lncn),
     &                          'i,jSW =', xi_start(0), eta_start(0) 
        do node=1,nnodes-1
          if (node.lt.16 .or. (nnodes.gt.16 .and.
     &                         node.eq.nnodes-1 )) then
            write(*,'(29x,I4,1x,A,2x,A,2I5)') node,
     &                  ncname(node)(1:lncn), 'i,jSW =',
     &                  xi_start(node), eta_start(node)
          elseif (nnodes.gt.16 .and. node.lt.18) then
            write(*,'(24x,A)') '.................................'
          endif
        enddo

#undef ntest
        if (ntest.ne.-1) then       ! Thus far netCDF file id array
          ierr=nf_close(ntest)      ! "ncid(0:nnodes-1)" was used just
          ntest=-1                  ! to signal that a complete set of
        endif                       ! partitioned files has been
        do node=0,nnodes-1          ! identified, but all the files are
          ncid(node)=-1             ! actually closed at this moment.
        enddo                       ! Reset the ids accordingly. 

      elseif (arg.lt.nargs) then
        goto 1
      else
        write(*,*) 'stop at 466'
        stop
      endif



!      *****    *********    ******   *******    *********
!    ***   ***  *  ***  *   **  ***   ***   ***  *  ***  *
!    ***           ***     **   ***   ***   ***     *** 
!      *****       ***    ***   ***   ***   **      ***
!          ***     ***    *********   ******        ***
!    ***   ***     ***    ***   ***   ***  **       ***
!      *****       ***    ***   ***   ***   ***     ***


! At this moment a set of files recorded as ncname(0:nnodes-1),
! xi_start(0:nnodes-1), eta_start(0:nnodes-1) comprise a complete
! set, but the files are actually in "closed" state and all netCDF
! IDs are reset to -1.

! Verify that ndims, ngatts, unlimdimid are the same for all nodes,
! however, note that different files may store different composition
! of variables, and netCDF variable IDs for the same variable (with
! the same name) may be different across the set of files.

#ifdef TIMING
        nclk=3-nclk
        call system_clock (iclk(nclk), clk_rate,clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        net_gray_clk=net_gray_clk+inc_clk
#endif
        do node=0,nnodes-1
          lncn=lenstr(ncname(node))
          if (ncid(node).eq.-1) ierr=nf_open (ncname(node),
     &                               nf_nowrite, ncid(node))
          if (ierr .eq. nf_noerr) then

            ierr=nf_inq (ncid(node), ibuff(1), ibuff(2),
     &                               ibuff(3), ibuff(4))
            if (ierr .ne. nf_noerr) then
              write(*,'(/1x,4A/12x,A/)')  '### ERROR: Cannot make ',
     &                        'general inquiry into netCDF file ''',
     &               ncname(node)(1:lncn), '''.', nf_strerror(ierr)
              goto 97
            elseif (ibuff(1) .gt. maxdims) then
              write(*,'(/1x,2A,I4,1x,3A/12x,2A/)')   '### ERROR: ',
     &         'number of dimensions', ibuff(1), 'in netCDF file ''',
     &          ncname(node)(1:lncn),    '''',   'exceeds limit.  ',
     &             'Increase parameter maxdims in file "ncjoin.F".'
              goto 97
            elseif (ibuff(2) .gt. maxvars) then
              write(*,'(/1x,2A,I4,1x,3A/12x,2A/)')   '### ERROR: ',
     &         'number of variables',  ibuff(2), 'in netCDF file ''',
     &          ncname(node)(1:lncn),  '''',      'exceeds limit. ',
     &             'Increase parameter maxvars in file "ncjoin.F".'
              goto 97
            elseif (node.eq.0) then
              ndims=ibuff(1)
c**           nvars=ibuff(2) 
              ngatts=ibuff(3)
              unlimdimid=ibuff(4)
            else
              if (ibuff(1) .ne. ndims) then
                write(*,'(/4x,4A/15x,3A/)')     '### ERROR: netCDF ',
     &                    'file ''', ncname(node)(1:lncn), ''' has ',
     &                    'different number of dimensions than ''',
     &                                       ncname(0)(1:lstr), '''.'
                ierr=ierr+1
              endif
c**           if (ibuff(2) .ne. nvars) then
c**             write(*,'(/4x,4A/15x,3A/)')     '### ERROR: netCDF ',
c**  &                    'file ''', ncname(node)(1:lncn), ''' has ',
c**  &                       'different number of variables than ''',
c**  &                                      ncname(0)(1:lstr), '''.'
c**             ierr=ierr+1
c**           endif
              if (ibuff(3) .ne. ngatts) then
                write(*,'(/4x,4A/15x,3A/)')     '### ERROR: netCDF ',
     &                   'file ''',  ncname(node)(1:lncn), ''' has ',
     &               'different number of global attributes than ''',
     &                                       ncname(0)(1:lstr),'''.'
                ierr=ierr+1
              endif 
              if (ibuff(4) .ne. unlimdimid) then
                write(*,'(/4x,4A/15x,3A/)')     '### ERROR: netCDF ',
     &                    'file ''', ncname(node)(1:lncn), ''' has ',
     &                'different ID for unlimited dimension than ''',
     &                                      ncname(0)(1:lstr), '''.'
                ierr=ierr+1
              endif
              if (ierr .ne. nf_noerr) goto 97 
            endif


! Verify that the sequence of dimension names is consistent
! throughout the entire set of variables.

#define i ilegal
            do j=1,ibuff(1)
              ierr=nf_inq_dim (ncid(node),j,string,dimsize(j,node))
              if (ierr .eq. nf_noerr) then
                lstr=lenstr(string)
                if (node.eq.0) then
                  ldim(j)=lstr
                  dimname(j)=string(1:lstr)
                elseif (lstr.ne.ldim(j) .or. string(1:lstr).ne.
     &                                 dimname(j)(1:ldim(j)) )then
                  write(*,'(/1x,2A,I3,3A/12x,6A/12x,3A/)')    '### ',
     &                 'ERROR: Name of dimension #', j, ', named ''',
     &                  string(1:lstr),   ''' in netCDF file',  '''',
     &                  ncname(node)(1:lncn),   ''' does not match ',
     &                 'name ''',  dimname(j)(1:ldim(j)), ''' with ',
     &                 'the corresponding name from netCDF file ''',
     &                                  ncname(0)(1:lncn),    '''.'
                  goto 97
                endif
              else
                write(*,'(/1x,2A,I3/12x,3A/12x,A)')    '### ERROR: ',
     &            'Cannot determine name and size of dimension #', j,
     &            'in netCDF file ''',   ncname(node)(1:lncn), '''.',
     &                                             nf_strerror(ierr)
                goto 97
              endif
            enddo
#undef i

! Create catalog of variable names, IDs, and ranks (number of
! dimensions) throughout the entire set. The meaning of the arrays
! defined here is as follows:
!
!     nvars    -- total number of variables discovered;
!     vname(i), where i=1:nvars -- variable name; 
!     vid(i,node) -- netCDF ID for that variable in netCDF file for
!                      MPI node "node", where node=0:nnodes-1;
!     vdims(i) -- "rank", i.e. namber of dimensions of that variable;
!     vnode(i) -- the file index for the lowest MPI node where the
!                 variable has been found, or has achieved first time
!                 its full rank (this is needed to differentiate
!                 between the true variable and its proxy "dummy
!                 scalar", if partit creates one, e.g., in the case
!                 of boundary forcing variable in inner MPI node,
!                 where it is not needed).   

            if (node.eq.0) nvars=0 
            do i=1,ibuff(2)
              ierr=nf_inq_varname (ncid(node), i, string)
              if (ierr .eq. nf_noerr) then
                lstr=lenstr(string)
                ierr=nf_inq_varndims (ncid(node), i, k)
                if (ierr .eq. nf_noerr) then
                  lnewvar=.true.
                  do j=1,nvars
                    lvar=lenstr(vname(j))
                    if (lstr.eq.lvar .and. string(1:lstr)
     &                               .eq.vname(j)(1:lvar)) then
                      lnewvar=.false.
                      vid(j,node)=i
                      if (k.gt.vdims(j)) then
                        vdims(j)=k
                        vnode(j)=node
                      endif
                    endif
                  enddo
                  if (lnewvar) then
                    nvars=nvars+1
                    vname(nvars)=string(1:lstr)
                    vid(nvars,node)=i
                    vnode(nvars)=node
                    vdims(nvars)=k
                  endif
                else
                  write(*,'(/1x,3A,I3/12x,5A/12x,A)')  '### ERROR: ',
     &                      'Cannot determine number of dimensions ',
     &                      'for variable with id =', i,  'named ''',
     &                       string(1:lstr),  ''' in netCDF file ''',
     &                       ncname(node)(1:lncn), '''.',
     &                                             nf_strerror(ierr)
                  goto 97
                endif
              else
                write(*,'(/1x,2A,I3/12x,3A/12x,A)')    '### ERROR: ',
     &                 'Cannot determine name of variable with id =',
     &                  i, 'in netCDF file ''', ncname(node)(1:lncn),
     &                                      '''.', nf_strerror(ierr)
                goto 97
              endif
            enddo
#ifdef KEEP_CLOSED
            if (node.gt.0) then           ! Close all the files,
              ierr=nf_close(ncid(node))   ! except for MPI rank=0.
              ncid(node)=-1
            endif
#endif
          else
            write(*,'(/1x,A,1x,3A/14x,A)')    '### ERROR: Cannot ',
     &                 'open netCDF file ''', ncname(node)(1:lncn),
     &                                    '''.', nf_strerror(ierr)
            goto 97
          endif
        enddo  !<-- node=0,nnodes-1 

#ifdef VERBOSE
        write(*,'(/1x,A,I3)') 'Inventory of variables: nvars =',nvars
        do i=1,nvars 
          lvar=lenstr(vname(i))
          write(*,'(4I4,2x,3A)') i, vid(i,vnode(i)), vnode(i),
     &                 vdims(i), '''', vname(i)(1:lvar), ''''
        enddo
        write(*,*) '...............................'
#endif

! Determine sizes of dimensions for combined file: For partitionable
! dimensions 'xi_rho', 'xi_u', 'eta_rho' and 'eta_v' determine the
! extent of the physical grid in each direction as the maximum over all
! subdomains of the dimension of each partial file combined with its
! starting index, "xi_start" or "eta_start".   This is straghtforward
! for RHO-points, but for U- and V-dimensions it requires to take into
! account the fact that the subdomains adjacent to eastern or southern
! edge have one point less than the corresponding RHO-dimension.
! Consequently, all subsequent subdomains receive one-point shift.
! For all other dimensions, verify that the sizes are the same for
! all nodes.  Also find size of unlimited dimension, if it is present.
! Note that variable "tsize" is set to its default value 1 (meaning
! one record), and may or may not be overwritten by the actual size of
! unlimited dimension (if it exists). If the unlimited dimension does
! not exist, it retains its value of 1 so that the loop over records
! is still executed, but only once.

        id_xi_rho=0
        id_eta_rho=0
        id_xi_rho=0
        id_eta_v=0

        XI_rho=0
        ETA_rho=0

        size_XI=1
        size_ETA=1
        size_S=1
        tsize=1

        do i=1,ndims
          dimsize(i,nnodes)=0
          lvar=lenstr(dimname(i))
          if (lvar.eq.6 .and. dimname(i)(1:lvar).eq.'xi_rho') then
            id_xi_rho=i
            do node=0,nnodes-1
              dimsize(i,nnodes)=max( dimsize(i,nnodes),
     &               dimsize(i,node) +xi_start(node)-1 )
              size_XI=max(size_XI,dimsize(i,node))
              XI_rho=max(XI_rho, dimsize(i,nnodes))
            enddo
           
          elseif (lvar.eq.4 .and.dimname(i)(1:lvar).eq.'xi_u') then
            id_xi_u=i
            do node=0,nnodes-1
              if (xi_start(node).gt.1) then
                dimsize(i,nnodes)=max( dimsize(i,nnodes),
     &                 dimsize(i,node) +xi_start(node)-2 )
              else
                dimsize(i,nnodes)=max( dimsize(i,nnodes),
     &                                  dimsize(i,node) )
              endif
              size_XI=max(size_XI,dimsize(i,node))
              XI_rho=max(XI_rho, dimsize(i,nnodes)+1)
            enddo

          elseif (lvar.eq.7.and. dimname(i)(1:lvar).eq.'eta_rho') then
            id_eta_rho=i
            do node=0,nnodes-1
              dimsize(i,nnodes)=max( dimsize(i,nnodes),
     &              dimsize(i,node) +eta_start(node)-1 )
              size_ETA=max(size_ETA,dimsize(i,node))
              ETA_rho=max(ETA_rho, dimsize(i,nnodes))
            enddo

          elseif (lvar.eq.5 .and. dimname(i)(1:lvar).eq.'eta_v') then
            id_eta_v=i
            do node=0,nnodes-1
              if (eta_start(node).gt.1) then
                dimsize(i,nnodes)=max( dimsize(i,nnodes),
     &                 dimsize(i,node) +eta_start(node)-2 )
              else
                dimsize(i,nnodes)=max( dimsize(i,nnodes),
     &                                   dimsize(i,node))
              endif
              size_ETA=max(size_ETA,dimsize(i,node))
              ETA_rho=max(ETA_rho, dimsize(i,nnodes)+1)
            enddo

          else
            dimsize(i,nnodes)=dimsize(i,0)
            do node=1,nnodes-1
              if (dimsize(i,0).ne.dimsize(i,node)) then
                lncn=lenstr(ncname(node))
                write(*,'(/1x,A,I3,3A,I4,1x,A/12x,4A/12x,3A,I4,A/)')
     &                 '### ERROR: Nonpartitionable dimension #',  i,
     &                 ' named ''', dimname(i)(1:lvar), ''', size =',
     &                  dimsize(i,node),   'in netCDF',    'file ''',
     &                  ncname(node)(1:lncn),    ''' has different ',
     &                 'size than the corresponding', 
     &                 'dimension from file ''',   ncname(0)(1:lncn),
     &                     ''', which has size =', dimsize(i,0), '.'
                goto 97
              endif
            enddo
            if (lvar.eq.5 .and. dimname(i)(1:lvar).eq.'s_rho') then
              size_S=max(size_S, dimsize(i,0))
            elseif (lvar.eq.3.and.dimname(i)(1:lvar).eq.'s_w') then
              size_S=max(size_S, dimsize(i,0))
            endif
          endif
          if (i.eq. unlimdimid) then
            tsize=dimsize(i,nnodes)
            dimsize(i,nnodes)=nf_unlimited
          endif
        enddo ! <-- i loop over dimensions


#ifdef VERBOSE
        write(*,'(1x,A)') 'Identifying presense of boundary edges:'
#endif
        do node=0,nnodes-1
          western_edge(node)=.true.
          eastern_edge(node)=.true.
          southern_edge(node)=.true.
          northern_edge(node)=.true.

          if (xi_start(node).gt.1) then
            western_edge(node)=.false.
          endif
          if (id_xi_rho.gt.0) then
            if ( xi_start(node)+dimsize(id_xi_rho,node)
     &          .lt.XI_rho ) eastern_edge(node)=.false. 
          endif
          if (id_xi_u.gt.0) then
            if ( xi_start(node)+dimsize(id_xi_u,node)
     &          .lt.XI_rho ) eastern_edge(node)=.false.
          endif
          if (eta_start(node).gt.1) then
            southern_edge(node)=.false.
          endif
          if (id_eta_rho.gt.0) then
            if ( eta_start(node)+dimsize(id_eta_rho,node)
     &           .lt.ETA_rho ) northern_edge(node)=.false.
          endif
          if (id_eta_v.gt.0) then
            if ( eta_start(node)+dimsize(id_eta_v,node)
     &          .lt.ETA_rho ) northern_edge(node)=.false.
          endif

#ifdef VERBOSE
          if (node.eq.0) then
            write(*,'(8x,A,I4,4(2x,A,L1))') 'node =', node,
     &     'WST=', western_edge(node),   'EST=', eastern_edge(node),
     &     'SOU=', southern_edge(node),  'NOR=', northern_edge(node)
          else
              write(*,'(14x,I4,4(6x,L1))') node,
     &                     western_edge(node),   eastern_edge(node),
     &                     southern_edge(node),  northern_edge(node)
          endif
#endif
        enddo


! set buffer size for needed to accommodate data for the largest
! possible MPI-subdomain (i.e., partial file).
! 
        max_buff_size=size_XI*size_ETA*size_S
!
! Create combined netCDF file:   Once the completeness of the set of
!------- -------- ------ -----   partial files have been established
! and dimensions survive consistency check, create the combined file,
! define its dimensions and copy global attributes.

        i=lenstr(root_bak)
        if (sffx_bak(1:1).ne.' ') then
          j=lenstr(sffx_bak)
          if (root_bak(i:i).eq.'.' .and. sffx_bak(1:1).eq.'.') then
            nctargname=root_bak(1:i)/ /sffx_bak(2:j)
          else
            nctargname=root_bak(1:i)/ /sffx_bak(1:j)
          endif
        else
          nctargname=root_bak(1:i)
        endif
        ltrg=lenstr(nctargname)

        j=0
        do i=1,ltrg
          if (nctargname(i:i).eq.'/') j=i+1
        enddo
        if (j.gt.0) then
          nctargname=nctargname(j:ltrg)
          ltrg=ltrg-j+1
        endif

!--> create...

        ierr=nf90_create (nctargname(1:ltrg), NF90_NETCDF4, ! old verson: nf_clobber+nf_64bit_offset
     &                                                          nctarg)
        if (ierr .eq. nf_noerr) then
          write(*,'(/1x,3A)')  'Created netCDF file ''',
     &                        nctargname(1:ltrg), '''.'
          write(*,*)  '-> compression level (deflate_level) = ',
     &                        deflate_level
          write(*,*)  '-> compression shuffle = ',
     &                        shuffle
        else
          write(*,'(/1x,4A/12x,A/)')     '### ERROR: Cannot create ',
     &                          'netCDF file ''', nctargname(1:ltrg),
     &                                      '''.', nf_strerror(ierr)
          goto 97
        endif

! Define dimensions: also compute the size of buffer needed to
! accommodate the largest array.

# ifdef VERBOSE
        write(*,'(/1x,A,3x,A,1x,A,1x,A)') 'Dimensions:', 'id',
     &                                           'size', 'name'
# endif
        size_XI=1
        size_ETA=1
        size_S=1

        do i=1,ndims
          lvar=lenstr(dimname(i))
          ierr=nf_def_dim (nctarg, dimname(i)(1:lvar),
     &                     dimsize(i,nnodes), dimid(i))
          if (ierr .eq. nf_noerr) then
            if (dimid(i) .eq. i) then
# ifdef VERBOSE
              write(*,'(14x,I3,I5,1x,3A)') dimid(i),
     &                               dimsize(i,nnodes), '''',
     &                               dimname(i)(1:lvar), ''''
# endif
              if (dimname(i)(1:3) .eq. 'xi_') then
                size_XI=max(size_XI, dimsize(i,nnodes))
              elseif (dimname(i)(1:4) .eq. 'eta_') then
                size_ETA=max(size_ETA, dimsize(i,nnodes))
              elseif (dimname(i)(1:5) .eq. 's_rho' .or.
     &                dimname(i)(1:3) .eq. 's_w') then
                size_S=max(size_S, dimsize(i,nnodes))
              endif
            else
              write(*,'(/1x,2A,I3,1x,5A/12x,2A,I3,A/)')  '### ERROR: ',
     &        'id =', dimid(i), 'for dimension ''', dimname(i)(1:lvar),
     &        ''' from netCDF file ''',    nctargname(1:ltrg),    '''',
     &                    'differs from ', 'the original id =', i, '.'
              goto 97
            endif
          else
           write(*,'(/1x,4A/12x,A/)')      '### ERROR: Cannot define ',
     &    'dimension ''', dimname(i)(1:lvar), '''.', nf_strerror(ierr)
              goto 97
          endif 
        enddo ! <-- i loop over dimensions
# ifdef FAST_MODE
        max_bfr_out=size_XI*size_ETA*size_S 
# endif

! Copy all global attributes, except 'partition'.

# ifdef VERBOSE
        write(*,'(1x,A)') 'Copying global attributes:'
# endif
        lncn=lenstr(ncname(0))
        if (ncid(0).eq.-1) ierr=nf_open (ncname(0), nf_nowrite,
     &                                                 ncid(0))
        if (ierr .eq. nf_noerr) then

          do i=1,ngatts
            ierr=nf_inq_attname (ncid(0), nf_global, i, string)
            if (ierr. eq. nf_noerr) then
              lvar=lenstr(string)
              if (string(1:lvar) .ne. 'partition') then
                ierr=nf_copy_att (ncid(0), nf_global, string(1:lvar),
     &                                             nctarg, nf_global)
                if (ierr .ne. nf_noerr) then
                  write(*,'(/1x,4A/12x,3A/12x,A)')     '### ERROR: ',
     &             'Cannot copy global attribute ''', string(1:lvar),
     &             ''' into netCDF',  'file ''',  nctargname(1:ltrg),
     &                                     '''.',  nf_strerror(ierr)
                  goto 97
                endif
# ifdef VERBOSE
                write(*,'(20x,3A)') '''', string(1:lvar), ''''
# endif      
              endif
            else
            write(*,'(/1x,2A,I3/12x,3A/12x,A/)') '### ERROR: Cannot',
     &                    ' determine name of global attribute #', i,
     &                     'from netCDF file ''',  ncname(0)(1:lncn),
     &                                     '''.',  nf_strerror(ierr)
              goto 97
            endif
          enddo
        else
          write(*,'(/1x,A,1x,3A/14x,A)')     '### ERROR: Cannot open ',
     &   'netCDF file ''', ncname(0)(1:lncn), '''.', nf_strerror(ierr)
          goto 97
        endif

! Define variables and copy their attributes.

#ifdef VERBOSE
        write(*,'(1x,2A)') 'Variables, their dimensions and ',
     &                                           'attributes:'
#endif
        do i=1,nvars
          node=vnode(i)
          lncn=lenstr(ncname(node))
#ifdef KEEP_CLOSED
          if (ncid(node).eq.-1) ierr=nf_open  (ncname(node),
     &                               nf_nowrite, ncid(node))
          if (ierr .eq. nf_noerr) then
#endif
            ierr=nf_inq_var (ncid(node), vid(i,node), vname(i),
     &                vartype(i), vdims(i), dimids(1,i),  varatts)
            if (ierr .eq. nf_noerr) then
              lvar=lenstr(vname(i))
              ! OLD F77 method:
!              ierr=nf_def_var (nctarg, vname(i)(1:lvar),vartype(i),
!     &                            vdims(i), dimids(1,i), varid(i))
              ! NEW F90 method:
              ! - nf90 routine no longer requires number of dimensions for
              !   the variable (hence vdims(i) not required).
              ! - However, dimids needs to be passed as an array.
              ! - Default compression is set to 1.
              !   Change to 'deflate_level = 0' for no compression.
              ! - Default is to turn shuffle on, this is what nc3to4z had.
              !   In tests it gave an extra ~10% compression over just d1 compression.
              if(vdims(i) > 2) then  ! 2D & 3D vars
                ierr=nf90_def_var (nctarg, vname(i)(1:lvar),vartype(i),
     &                             dimids(1:vdims(i),i), varid(i),
     &                             deflate_level=deflate_level, shuffle=shuffle)
!         ierr=nf90_inquire_variable(nctarg,varid(i),chunksizes=chunk_sizes) ! DevinD debug
!         print *, 'var: ',vname(i)(1:lvar),' chunksizes: ',chunk_sizes      ! DevinD debug
              else ! 1D and grd vars
                ierr=nf90_def_var (nctarg, vname(i)(1:lvar),vartype(i),
     &                             dimids(1:vdims(i),i), varid(i))
              endif

              if (ierr .eq. nf_noerr) then
#ifdef VERBOSE
                write(*,'(8x,3A,8I3)')   '''', vname(i)(1:lvar),
     &              ''',  dimids =', (dimids(j,i), j=1,vdims(i))
#endif
                do j=1,varatts
                  ierr=nf_inq_attname (ncid(node), vid(i,node),
     &                                              j, string)
                  if (ierr .eq. nf_noerr) then
                    lstr=lenstr(string)
                    ierr=nf_copy_att (ncid(node), vid(i,node),
     &                       string(1:lstr), nctarg, varid(i))
                    if (ierr. ne. nf_noerr) then
                      write(*,'(/1x,2A,I3,3A/12x,4A)')   '### ERROR: ',
     &                 'Cannot copy attribute #', j,' for variable ''',
     &                  vname(i)(1:lvar),  ''' into netCDF', 'file ''',
     &                  nctargname(1:ltrg), '''.  ', nf_strerror(ierr)
                      goto 97
                    endif
#ifdef VERBOSE
                     write(*,'(16x,3A)') '''', string(1:lstr), ''''
#endif
                  else
                    write(*,'(/1x,2A,I3/12x,3A/12x,A/)') '### ERROR: ',
     &                             'Cannot get name of attribute #', j,
     &                      'for variable ''', vname(i)(1:lvar), '''.',
     &                                               nf90_strerror(ierr)
                    goto 97
                  endif
                enddo
              else
                write(*,'(/8x,4A/)') 'ERROR: Cannot define ',
     &                  'variable ''', vname(i)(1:lvar), '''.'
                goto 97
              endif
            else

            write(*,'(/8x,2A/15x,A,I3,1x,3A/)')  '### ERROR: Cannot ',
     &        'determine name, type and attributes for variable #', i,
     &            'from netCDF file ''', ncname(node)(1:lncn),  '''.'
              goto 97
            endif

! Determine whether partitionable dimensions or unlimited dimension
! are present for this variable: the convention adopted here is:
!           part_type = 0 -- non-partitionable array;
!                     = 1 -- has partitionable XI-dimension only;
!                     = 2 -- has partitionable ETA-dimension only;
!                     = 3 -- partitionable in both XI and ETA.

            series(i)=.false.
            part_type(i)=0
            do j=1,vdims(i)
              if (dimids(j,i).eq.id_xi_rho .or.
     &            dimids(j,i).eq.id_xi_u) then
                part_type(i)=part_type(i)+1
              elseif (dimids(j,i).eq.id_eta_rho .or.
     &                dimids(j,i).eq.id_eta_v) then
                part_type(i)=part_type(i)+2
              elseif (dimids(j,i).eq.unlimdimid) then
                series(i)=.true.
              endif
            enddo

#ifdef KEEP_CLOSED
            if (node.gt.0) then
              ierr=nf_close(ncid(node))
              ncid(node)=-1
            endif
          else
            write(*,'(/1x,A,1x,3A/12x,A)')  '### ERROR: Cannot open ',
     &                  'netCDF file ''', ncname(node)(1:lncn), '''.',
     &                                              nf_strerror(ierr)
            goto 97
          endif
#endif
        enddo       ! <-- i=1,nvars, variable IDs.

! Leave definition mode

        ierr=nf_enddef (nctarg)
# ifdef VERBOSE
        write(*,'(/1x,A)') 'Leaving definition mode.'
# endif





!     **     *     ***  *******    ***  *********  ********
!      *    ***   ***   ***   ***  ***  *  ***  *  ***    *
!       *   ***   ***   ***   ***  ***     ***     ***
!       *  *** * ***    ***   **   ***     ***     ******
!        * **  * **     ******     ***     ***     ***
!        ***    ***     ***  **    ***     ***     ***    *
!         *     **      ***   ***  ***     ***     ********
!

! Allocate necessary buffer arrays: here "max_*_size" is needed size
! of buffer array, as determined by processing of current file, while
! "alloc_*_size" is size of array already allocated.  Both are
! initialized to zero at the beginning.   Basically this anticipates
! the possibility of gradual increase of needed size of buffer array,
! so that if it is allocated at this moment, but insufficient, it is
! deallocated first. 

#ifdef FAST_MODE
        if (max_bfr_out .gt. alloc_bfr_out) then
          if (allocated(bfr_out)) deallocate(bfr_out)
        endif
#endif
        if (max_buff_size .gt. alloc_buff_size) then
          if (allocated(buff)) deallocate(buff)
          allocate(buff(max_buff_size))
          alloc_buff_size=max_buff_size
          write(*,*) 'allocated "buff" with  max_buff_size =',
     &                                         max_buff_size
        endif
#ifdef FAST_MODE
        if (max_bfr_out .gt. alloc_bfr_out) then
          allocate(bfr_out(max_bfr_out))
          alloc_bfr_out=max_bfr_out
          write(*,*) 'allocated "bfr_out" with ',
     &               'max_bfr_out =', max_bfr_out
        endif
#endif
#ifdef TIMING
        nclk=3-nclk
        call system_clock (iclk(nclk), clk_rate,clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        net_fcrt_clk=net_fcrt_clk+inc_clk
#endif

! Transfer variables into combined file.

        do rec=1,tsize
          if (tsize.gt.1) then
#ifdef TIMING
            nclk=3-nclk
            call system_clock (iclk(nclk), clk_rate,clk_max)
            inc_clk=iclk(nclk)-iclk(3-nclk)
            net_gray_clk=net_gray_clk+inc_clk

            write(*,'(F8.1,1x,A,I8,1x,A,I8,1x,A)')
     &         dble(iclk(nclk)-iclk_init)/dble(clk_rate),
#else
            write(*,'(4x,A,I8,1x,A,I8,1x,A)')
#endif
     &     'Processing record', rec, 'out of', tsize,  '...'
          endif

          do i=1,nvars
            if (rec.eq.1 .or. series(i)) then
              if (part_type(i).eq.0 .and. .not.series(i)) then

! Scalar (zero-dimensional) variables:

                lvar=lenstr(vname(i))
                write(*,'(16x,3A)') 'Copy scalar variable: ''',
     &                               vname(i)(1:lvar), '''...'

                if (vartype(i) .eq. nf_char) then
                  ierr=nf_get_var_text  (ncid(0), vid(i,0), buff_char)
                elseif (vartype(i) .eq. nf_byte) then
                  ierr=nf_get_var_int1   (ncid(0), vid(i,0), buff)
                elseif (vartype(i) .eq. nf_short) then
                  ierr=nf_get_var_int2   (ncid(0), vid(i,0), buff)
                elseif (vartype(i) .eq. nf_int) then
                  ierr=nf_get_var_int    (ncid(0), vid(i,0), buff)
                elseif (vartype(i) .eq. nf_float) then
                  ierr=nf_get_var_real   (ncid(0), vid(i,0), buff)
                elseif (vartype(i) .eq. nf_double) then
                  ierr=nf_get_var_double (ncid(0), vid(i,0), buff)
                else
                  lvar=lenstr(vname(i))
                  write(*,'(/8x,4A/)') '### ERROR: scalar variable ',
     &              '''', vname(i)(1:lvar), ''' has unknown type.'
                  goto 97
                endif
                if (ierr .eq. nf_noerr) then
                  if (vartype(i) .eq. nf_char) then
                    ierr=nf_put_var_text  (nctarg, varid(i), buff_char)
                  elseif (vartype(i) .eq. nf_byte) then
                    ierr=nf_put_var_int1   (nctarg, varid(i), buff)
                  elseif (vartype(i) .eq. nf_short) then
                    ierr=nf_put_var_int2   (nctarg, varid(i), buff)
                  elseif (vartype(i) .eq. nf_int) then
                    ierr=nf_put_var_int    (nctarg, varid(i), buff)
                  elseif (vartype(i) .eq. nf_float) then
                    ierr=nf_put_var_real   (nctarg ,varid(i), buff)
                  elseif (vartype(i) .eq. nf_double) then
                    ierr=nf_put_var_double (nctarg, varid(i), buff)
                  endif
                  if (ierr .ne. nf_noerr) then
                    lvar=lenstr(vname(i))
                    write(*,'(/1x,4A/12x,3A/12x,A)')  '### ERROR: ',
     &                            'Cannot write scalar variable ''',
     &                           vname(i)(1:lvar), ''' into netCDF',
     &                                'file ''', nctargname(1:ltrg),
     &                                    '''.',  nf_strerror(ierr)
                    goto 97
                  endif
                else
                  lvar=lenstr(vname(i))
                  write(*,'(/1x,4A/12x,A/)')  '### ERROR: Cannot ',
     &                 'read scalar variable ''', vname(i)(1:lvar),
     &                                    '''.', nf_strerror(ierr)
                  goto 97
                endif
              elseif (part_type(i).eq.0) then

! Non-partitionable array.

                lvar=lenstr(vname(i))
                write(*,'(16x,3A)') 'Copy non-partitioned array: ''',
     &                                    vname(i)(1:lvar), '''...'
                size=1
                do j=1,vdims(i)
                  if (dimids(j,i).eq.unlimdimid) then
                    start(j)=rec
                    count(j)=1
                  else
                    start(j)=1
                    count(j)=dimsize(dimids(j,i),0)
                  endif
                  size=size*count(j)
                enddo
                if (vartype(i) .eq. nf_char .or.
     &              vartype(i) .eq. nf_byte) then
                  size=size*1
                elseif (vartype(i) .eq. nf_short) then
                  size=size*2
                elseif (vartype(i) .eq. nf_int .or.
     &                  vartype(i) .eq. nf_float) then
                  size=size*4
                elseif (vartype(i) .eq. nf_double) then
                  size=size*8
                else
                  lvar=lenstr(vname(i))
                  write(*,'(/8x,3A/)')  '### ERROR: variable ''',
     &                  vname(i)(1:lvar), ''' has unknown type.'
                  goto 97
                endif
                if (size .gt. 8*max_buff_size) then
                  if (allocated(buff)) deallocate(buff)
                  max_buff_size=(size+7)/8
                  allocate(buff(max_buff_size))
                  write(*,*) 'allocated "buff" with max_buff_size =',
     &                                              max_buff_size
                endif
 
                if (vartype(i) .eq. nf_char) then
                  ierr=nf_get_vara_text  (ncid(0), vid(i,0),
     &                                    start,count, buff)
                elseif (vartype(i) .eq. nf_byte) then
                  ierr=nf_get_vara_int1   (ncid(0), vid(i,0),
     &                                    start,count, buff)
                elseif (vartype(i) .eq. nf_short) then
                  ierr=nf_get_vara_int2   (ncid(0), vid(i,0),
     &                                    start,count, buff)
                elseif (vartype(i) .eq. nf_int) then
                  ierr=nf_get_vara_int    (ncid(0), vid(i,0),
     &                                    start,count, buff)
                elseif (vartype(i) .eq. nf_float) then
                  ierr=nf_get_vara_real   (ncid(0), vid(i,0),
     &                                    start,count, buff)
                elseif (vartype(i) .eq. nf_double) then
                  ierr=nf_get_vara_double (ncid(0), vid(i,0),
     &                                    start,count, buff)
                endif
                if (ierr .eq. nf_noerr) then
                  if (vartype(i) .eq. nf_char) then
                    ierr=nf_put_vara_text  (nctarg, varid(i),
     &                                    start,count, buff)
                  elseif (vartype(i) .eq. nf_byte) then
                    ierr=nf_put_vara_int1  (nctarg, varid(i),
     &                                    start,count, buff)
                  elseif (vartype(i) .eq. nf_short) then
                    ierr=nf_put_vara_int2  (nctarg, varid(i),
     &                                    start,count, buff)
                  elseif (vartype(i) .eq. nf_int) then
                    ierr=nf_put_vara_int   (nctarg, varid(i),
     &                                    start,count, buff)
                  elseif (vartype(i) .eq. nf_float) then
                    ierr=nf_put_vara_real  (nctarg, varid(i),
     &                                    start,count, buff)
                  elseif (vartype(i) .eq. nf_double) then
                    ierr=nf_put_vara_double(nctarg, varid(i),
     &                                    start,count, buff)
                  endif
                  if (ierr .ne. nf_noerr) then
                    lvar=lenstr(vname(i))
                    write(*,'(/8x,4A,I3/15x,3A/)')    '### ERROR: ',
     &                 'Cannot write variable ''', vname(i)(1:lvar),
     &              ''' for time record',rec, 'into netCDF file ''',
     &                  nctargname(1:ltrg),'''.', nf_strerror(ierr)
                    goto 97
                  endif
                else
                  lvar=lenstr(vname(i))
                  write(*,'(/8x,4A,I3,A/15x,A/)')     '### ERROR: ',
     &              'Cannot read variable ''',     vname(i)(1:lvar),
     &              ''' for time record',rec,'.', nf_strerror(ierr)
                  goto 97
                endif
              elseif (part_type(i).gt.0) then

! Partitioned array:

                lvar=lenstr(vname(i))
                write(*,'(16x,2A,I3,1x,3A)')  'Assembly partitioned ',
     &                       'array type', part_type(i),  'name = ''',
     &                                         vname(i)(1:lvar), ''''

                do node=0,nnodes-1
                  var_mask=.false.
                  if (part_type(i).eq.1 .and. lvar.gt.6) then 
                    if (vname(i)(lvar-5:lvar).eq.'_south' 
     &                        .and. southern_edge(node)) then
                      var_mask=.true.
#ifdef VERBOSE
                      write(*,'(3x,A,I4,1x,4A)')  'node =', node,
     &               'identified XI-partitioned southern ',
     &               'boundary array ''', vname(i)(1:lvar), ''''
#endif
                    elseif (vname(i)(lvar-5:lvar).eq.'_north'
     &                           .and. northern_edge(node)) then
                      var_mask=.true.
#ifdef VERBOSE
                      write(*,'(3x,A,I4,1x,4A)')   'node =', node,
     &               'identified XI-partitioned northern ',
     &               'boundary array ''', vname(i)(1:lvar), ''''
#endif
                    endif
                  elseif (part_type(i).eq.2 .and. lvar.gt.5) then
                    if (vname(i)(lvar-4:lvar).eq.'_west'
     &                                .and. western_edge(node)) then
                      var_mask=.true.
#ifdef VERBOSE
                      write(*,'(3x,A,I4,1x,4A)')     'node =', node,
     &               'identified ETA-partitioned western boundary ',
     &               'array ''', vname(i)(1:lvar), ''''
#endif
                    elseif (vname(i)(lvar-4:lvar).eq.'_east'
     &                          .and. eastern_edge(node)) then
                      var_mask=.true.
#ifdef VERBOSE
                      write(*,'(3x,A,I4,1x,4A)')     'node =', node,
     &               'identified ETA-partitioned eastern boundary ',
     &               'array ''', vname(i)(1:lvar), ''''
#endif
                    endif
                  elseif (part_type(i).eq.3) then
                    var_mask=.true.
#ifdef VERBOSE
                    write(*,'(3x,A,I4,1x,4A)')   'node =',  node,
     &                      'identified 2D-partitioned array ''',
     &                       vname(i)(1:lvar), ''''
#endif

                  endif
                  if (var_mask) then
                    size=1
                    size1=1
                    do j=1,vdims(i)
                      k=dimids(j,i)
                      if (k.eq.id_xi_rho  .or. k.eq.id_xi_u  .or.
     &                    k.eq.id_eta_rho .or. k.eq.id_eta_v) then
                        start(j)=1
                        count(j)=dimsize(k,node)

                        if (k.eq.id_xi_rho) then
                          start1(j)=xi_start(node)
#ifdef FAST_MODE
                          count1(j)=XI_rho
#endif
                        elseif (k.eq.id_xi_u) then
                          start1(j)=max(xi_start(node)-1,1)
#ifdef FAST_MODE
                          count1(j)=XI_rho-1
#endif
                        elseif (k.eq.id_eta_rho) then
                          start1(j)=eta_start(node)
#ifdef FAST_MODE
                          count1(j)=ETA_rho
#endif
                        elseif (k.eq.id_eta_v) then
                          start1(j)=max(eta_start(node)-1,1)
#ifdef FAST_MODE
                          count1(j)=ETA_rho-1
#endif
                        endif

                      elseif (k.eq.unlimdimid) then
                        start(j)=rec
                        count(j)=1
                        start1(j)=rec
#ifdef FAST_MODE
                        count1(j)=1
#endif
                      else
                        start(j)=1
                        count(j)=dimsize(k,nnodes)
                        start1(j)=1
#ifdef FAST_MODE
                        count1(j)=count(j)
#endif
                      endif
                      size=size*count(j)
#ifdef FAST_MODE
                      size1=size1*count1(j)
#else
                      size1=size*count(j)
#endif
                    enddo

! Convert sizese to Bytes

                    if (vartype(i) .eq. nf_char .or.
     &                  vartype(i) .eq. nf_byte) then
                      size=size*1
                      size1=size1*1
                    elseif (vartype(i) .eq. nf_short) then
                      size=size*2
                      size1=size1*2
                    elseif (vartype(i) .eq. nf_int .or.
     &                      vartype(i) .eq. nf_float) then
                      size=size*4
                      size1=size1*4
                    elseif (vartype(i) .eq. nf_double) then
                      size=size*8
                      size1=size1*8
                    else
                      lvar=lenstr(vname(i))
                      write(*,'(/8x,4A/)') '### ERROR: variable ''',
     &                     vname(i)(1:lvar), ''' has unknown type.'
!dd                      goto 97
                    endif

!                   print *,size,max_buff_size*8

                    if (size .gt. 8*max_buff_size) then
                      if (allocated(buff)) deallocate(buff)
                      max_buff_size=(size+7)/8
                      allocate(buff(max_buff_size))
                      write(*,*) 'allocated "buff" with ',
     &                   'max_buff_size =', max_buff_size
                    endif

#ifdef KEEP_CLOSED
                    if (ncid(node).eq.-1) ierr=nf_open (ncname(node),
     &                                       nf_nowrite, ncid(node))
                    if (ierr.eq.nf_noerr) then
#endif
#ifdef VERBOSE
                      write(*,'(3x,A,I4,2x,A,I7,2x,A,I4,1x,A,I4,1x,A)')
     &                     'node =', node, 'ncid=',ncid(node),
     &                     'xi_start =', xi_start(node),
     &                     'eta_start =', eta_start(node), 'reading...'
#endif

# ifdef TIMING
                      nclk=3-nclk
                      call system_clock (iclk(nclk), clk_rate,clk_max)
                      inc_clk=iclk(nclk)-iclk(3-nclk)
                      net_gray_clk=net_gray_clk+inc_clk
# endif
                      ! DevinD: nf90_get_var is one call, it doesn't need
                      ! to know the data type.
                      if (vartype(i) .eq. nf_char) then
                        ! DevinD: This call is technically illegal, and will cause
                        ! ncjoin to abort in debug mode because you cannot read a
                        ! 'character' into a 'real'.
                        ! This final result still seems correct though, hence not
                        ! previously noticed. (Only for 'spherical' which is not
                        ! commonly looked at.
                        ierr=nf_get_vara_text (ncid(node), vid(i,node),
     &                                             start, count, buff)
                      elseif (vartype(i) .eq. nf_byte) then
                        ierr=nf_get_vara_int1 (ncid(node), vid(i,node),
     &                                             start, count, buff)
                      elseif (vartype(i) .eq. nf_short) then
                        ierr=nf_get_vara_int2 (ncid(node), vid(i,node),
     &                                             start, count, buff)
                      elseif (vartype(i) .eq. nf_int) then
                        ierr=nf_get_vara_int  (ncid(node), vid(i,node),
     &                                             start, count, buff)
                      elseif (vartype(i) .eq. nf_float) then
!                        print *, 'READING FLOAT: '
                        ierr=nf_get_vara_real (ncid(node), vid(i,node),
     &                                             start, count, buff)
                      elseif (vartype(i) .eq. nf_double) then
                        ierr=nf_get_vara_double(ncid(node),vid(i,node),
     &                                             start, count, buff)
                      endif
# ifdef TIMING
                      if (ierr.eq.nf_noerr) then
                        net_read_size=net_read_size+size
                        nclk=3-nclk
                        call system_clock (iclk(nclk),clk_rate,clk_max)
                        inc_clk=iclk(nclk)-iclk(3-nclk)
                        net_read_clk=net_read_clk+inc_clk
                      else
# else
                      if (ierr.ne.nf_noerr) then
# endif
                        lvar=lenstr(vname(i))
                        lncn=lenstr(ncname(node))
                        write(*,'(/1x,4A,I3/15x,3A/15x,A/)')  '### ',
     &                              'ERROR: Cannot read variable ''',
     &                   vname(i)(1:lvar), ''' for time record', rec,
     &                  'from netCDF file ''',  ncname(node)(1:lncn),
     &                                      '''.', nf_strerror(ierr)
!dd                        goto 97
                      endif

#ifdef KEEP_CLOSED
                      if (node.gt.0) then
                        ierr=nf_close(ncid(node))
                        ncid(node)=-1
                      endif
                    else
                      lncn=lenstr(ncname(node))
                      write(*,'(/1x,A,1x,3A/14x,A)')   '### ERROR: ',
     &                                  'Cannot open netCDF file ''',
     &                ncname(node)(1:lncn), '''.', nf_strerror(ierr)
!dd                      goto 97
                    endif
#endif
# ifdef VERBOSE
                    write(*,'(1x,A)')  'copying ...'
# endif

! In the code segment below there are two strategies, resulting in
! equivalent outcome, but different performance: CPP-switch FAST_MODE
! activated code in which the partitioned array is first assembled
! into the intermediate buffer array "bfr_out" horizontal dimensions
! of which correspond to the whole physical grid. After that, the
! whole array is defined, it is written once as a single record [this
! occurs immediately after the reading of the last node]. Note that
! all values of "count1(k)" [except for "k" corresponding to
! unlimited dimension] are equal to the actual dimensions of the
! variable in netCDF file, while all all "start(k)" [except unlimited
! dimension] are equal to 1.

! Alternatively, partial rectangular blocks are written immediately
! after they read from the partial netCDF files, relying on sub-array
! writing capability of "nf_put_vara_XXX". This results physical
! write into a large number of small records of size count(1).


#ifdef FAST_MODE
                    if (size1 .gt. 8*max_bfr_out) then
                      if (allocated(bfr_out)) deallocate(bfr_out)
                      max_bfr_out=(size1+7)/8
                      allocate(bfr_out(max_bfr_out))
                      write(*,*) 'allocated "bfr_out" with ',
     &                           'max_bfr_out =', max_bfr_out 
                    endif
# ifdef TIMING
                    nclk=3-nclk
                    call system_clock (iclk(nclk), clk_rate, clk_max)
                    inc_clk=iclk(nclk)-iclk(3-nclk)
                    net_gray_clk=net_gray_clk+inc_clk
# endif
                    if (vartype(i) .eq. nf_char) then
                      call assembly_text (buff, count, bfr_out,
     &                                 start1,count1, vdims(i))
                    elseif (vartype(i) .eq. nf_byte) then
                      call assembly_byte (buff, count, bfr_out, 
     &                                 start1,count1, vdims(i))
                    elseif (vartype(i) .eq. nf_short) then
                      call assembly_int2 (buff, count, bfr_out,
     &                                 start1,count1, vdims(i))
                    elseif (vartype(i) .eq. nf_int) then
                      call assembly_int  (buff, count, bfr_out,
     &                                 start1,count1, vdims(i))
                    elseif (vartype(i) .eq. nf_float) then
                      call assembly_real (buff, count, bfr_out,
     &                                 start1,count1, vdims(i))
                    elseif (vartype(i) .eq. nf_double) then
                      call assembly_double (buff,count,bfr_out,
     &                                 start1,count1, vdims(i))
                    endif
# ifdef TIMING
                    nclk=3-nclk
                    call system_clock (iclk(nclk), clk_rate, clk_max)
                    inc_clk=iclk(nclk)-iclk(3-nclk)
                    net_assm_clk=net_assm_clk+inc_clk
# endif
#else
# ifdef TIMING
                    nclk=3-nclk
                    call system_clock (iclk(nclk), clk_rate, clk_max)
                    inc_clk=iclk(nclk)-iclk(3-nclk)
                    net_gray_clk=net_gray_clk+inc_clk
# endif

                    if (vartype(i) .eq. nf_char) then
                      ierr=nf_put_vara_text   (nctarg, varid(i),
     &                                     start1, count, buff)
                    elseif (vartype(i) .eq. nf_byte) then
                      ierr=nf_put_vara_int1   (nctarg, varid(i),
     &                                     start1, count, buff)
                    elseif (vartype(i) .eq. nf_short) then
                      ierr=nf_put_vara_int2   (nctarg, varid(i),
     &                                     start1, count, buff)
                    elseif (vartype(i) .eq. nf_int) then
                      ierr=nf_put_vara_int    (nctarg, varid(i),
     &                                     start1, count, buff)
                    elseif (vartype(i) .eq. nf_float) then
                       ierr=nf_put_vara_real  (nctarg, varid(i),
     &                                     start1, count, buff)
                    elseif (vartype(i) .eq. nf_double) then
                       ierr=nf_put_vara_double(nctarg, varid(i),
     &                                     start1, count, buff)
                    endif
# ifdef TIMING
                    if (ierr.eq.nf_noerr) then
                      net_wrt_size=net_wrt_size+size

                      nclk=3-nclk
                      call system_clock(iclk(nclk), clk_rate,clk_max)
                      inc_clk=iclk(nclk)-iclk(3-nclk)
                      net_wrt_clk=net_wrt_clk+inc_clk
                    else
# else
                    if (ierr.ne.nf_noerr) then
# endif
                      lvar=lenstr(vname(i))
                      lncn=lenstr(vname(i))
                      write(*,'(/1x,3A,I3/12x,3A/12x,A/)')
     &                         '### ERROR: Cannot write variable ''',
     &                     vname(i)(1:lvar),''' for time record',rec,
     &                     'into netCDF file ''', nctargname(1:ltrg),
     &                                     '''.',  nf_strerror(ierr)
!dd                      goto 97
                    endif
#endif /* FAST_MODE */
                  endif  ! <-- var_mask 

#ifdef FAST_MODE
                  if (node .eq. nnodes-1) then
# ifdef TIMING
                    nclk=3-nclk
                    call system_clock (iclk(nclk), clk_rate, clk_max)
                    inc_clk=iclk(nclk)-iclk(3-nclk)
                    net_gray_clk=net_gray_clk+inc_clk
# endif

                    if (vartype(i) .eq. nf_char) then
                      ierr=nf_put_vara_text   (nctarg, varid(i),
     &                                  start, count1, bfr_out)
                    elseif (vartype(i) .eq. nf_byte) then
                      ierr=nf_put_vara_int1   (nctarg, varid(i),
     &                                  start, count1, bfr_out)
                    elseif (vartype(i) .eq. nf_short) then
                      ierr=nf_put_vara_int2   (nctarg, varid(i),
     &                                  start, count1, bfr_out)
                    elseif (vartype(i) .eq. nf_int) then
                      ierr=nf_put_vara_int    (nctarg, varid(i),
     &                                  start, count1, bfr_out)
                    elseif (vartype(i) .eq. nf_float) then
                      ierr=nf_put_vara_real   (nctarg, varid(i),
     &                                  start, count1, bfr_out)
                    elseif (vartype(i) .eq. nf_double) then
                      ierr=nf_put_vara_double (nctarg, varid(i),
     &                                  start, count1, bfr_out)
                    endif
# ifdef TIMING
                    if (ierr.eq.nf_noerr) then
                      net_wrt_size=net_wrt_size+size1

                      nclk=3-nclk
                      call system_clock(iclk(nclk), clk_rate,clk_max)
                      inc_clk=iclk(nclk)-iclk(3-nclk)
                      net_wrt_clk=net_wrt_clk+inc_clk
                    else
# else
                    if (ierr.ne.nf_noerr) then
# endif
                      lvar=lenstr(vname(i))
                      lncn=lenstr(vname(i))
                      write(*,'(/1x,3A,I3/12x,3A/12x,A/)')
     &                         '### ERROR: Cannot write variable ''',
     &                     vname(i)(1:lvar),''' for time record',rec,
     &                    'into netCDF file ''',  nctargname(1:ltrg),
     &                                     '''.',  nf_strerror(ierr)
!dd                      goto 97
                    endif
                  endif    !<--  node.eq.nnodes-1
#endif /* FAST_MODE */
                enddo   !<-- node=0,nnodes-1

              endif    !<-- part_type .eq./.gt. 0 switch
            endif    !<-- rec.eq.1 .or. series(i) switch

c        ierr=nf_sync (nctarg)

          enddo    !<-- i=1,nvars, loop over variables

! Use intermediate one-per-record nf_sync of the target file only in
! verbose mode because it slows down the execution speed.   Note that
! time spent by nf_sync is counted as writing time, and so does time
! spent to close the target file. 

c--#ifdef VERBOSE
# if defined TIMING
          nclk=3-nclk
          call system_clock(iclk(nclk), clk_rate, clk_max)
          inc_clk=iclk(nclk)-iclk(3-nclk)
          net_gray_clk=net_gray_clk+inc_clk
# endif
          ierr=nf_sync (nctarg)

c          write(*,*) 'close' 
c          ierr=nf_close (nctarg) 
c          write(*,*) 'reopen'
c          ierr=nf_open (nctargname(1:ltrg), nf_write, nctarg)
# if defined TIMING
          nclk=3-nclk
          call system_clock(iclk(nclk), clk_rate, clk_max)
          inc_clk=iclk(nclk)-iclk(3-nclk)
          net_sync_clk=net_sync_clk+inc_clk
# endif
c--#endif 

        enddo    !<--- rec, loop over records

        if (ierr.eq.nf_noerr) then
          clean_set=.true.
          goto 98
        endif
  97    clean_set=.false. 

! Close all files                        ! At this moment open/closed
                                         ! status of partial files
  98    write(*,*) 'closing files...'    ! depends on the state of CPP
        do node=0,nnodes-1               ! switch KEEP_CLOSED.   If the
          if (ncid(node).ne.-1) then     ! switch is defined, then only
            ierr=nf_close(ncid(node))    ! node=0 file is expected to
            ncid(node)=-1                ! be opened here.   Otherwise
          endif                          ! the entire set is opened and
        enddo                            ! needs to be closed.   Either
        write(*,*) '...........input'    ! way, as ncid(node).eq/ne.-1
#if defined TIMING
        nclk=3-nclk
        call system_clock(iclk(nclk), clk_rate, clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        net_gray_clk=net_gray_clk+inc_clk
#endif
        ierr=nf_close (nctarg)           ! is used as flag indicating
#if defined TIMING
        nclk=3-nclk
        call system_clock(iclk(nclk), clk_rate, clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        net_sync_clk=net_sync_clk+inc_clk
#endif
        write(*,*) '...........output'   ! status of each file. 

#ifdef DEL_PART_FILES
        if (del_part_files) then
          if (clean_set) then
# ifdef TIMING
            nclk=3-nclk
            call system_clock (iclk(nclk), clk_rate, clk_max)
            inc_clk=iclk(nclk)-iclk(3-nclk)
            net_gray_clk=net_gray_clk+inc_clk
# endif
            write(*,'(/1x,A)') 'Deleting partial files...'
            do node=0,nnodes-1
              rmcmd='/bin/rm -f '/ /ncname(node)
              lstr=lenstr(rmcmd)
              if (node.lt.16 .or. (nnodes.gt.16 .and.
     &                         node.eq.nnodes-1 )) then
                write(*,'(27x,3A)')  '''', rmcmd(1:lstr), ''''

              elseif (nnodes.gt.16 .and. node.lt.18) then
                write(*,'(24x,A)') '.................................'
              endif
              call system (rmcmd(1:lstr))
            enddo
            write(*,*)
# ifdef TIMING
            nclk=3-nclk
            call system_clock (iclk(nclk), clk_rate, clk_max)
            inc_clk=iclk(nclk)-iclk(3-nclk)
            net_rmcmd_clk=net_rmcmd_clk+inc_clk
# endif
          else
            write(*,'(/1x,2A/)')  '### ERROR: Not removing ',
     &                     'partial files because of errors.'
          endif
        endif
#endif
      if (arg .lt. nargs)  goto 1     !--> next set of partial files.

#ifdef TIMING
# ifdef __IFC
      RUN_time=etime(CPU_time)
# else
      call etime(CPU_time, RUN_time)
# endif
      RUN_time=RUN_time-tstart
      

      write(*,'(/3(1x,A,F11.2,1x))') 'CPU_time:  run =', RUN_time,
     &                   'usr =', CPU_time(1),  'sys =', CPU_time(2)

      if (clk_rate.gt.0) then
        ReadSize=1.0D-6*net_read_size
        WrtSize=1.0D-6*net_wrt_size  
        ReadTime=net_read_clk/dble(clk_rate)
        AssmTime=net_assm_clk/dble(clk_rate)
        WrtTime = net_wrt_clk/dble(clk_rate)
        SyncTime=net_sync_clk/dble(clk_rate)
        FcrtTime=net_fcrt_clk/dble(clk_rate)

        write(*,'(/1x,A,22x,F12.2,1x,A)') 'Analysis/file creation :',
     &                                               FcrtTime, 'sec'
        write(*,'(8x,A,F12.2,1x,A,F12.2,1x,A,F8.2,1x,A)')
     &         'Total data read :', ReadSize, 'MBytes in',  ReadTime,
     &                          'sec (', ReadSize/ReadTime, 'MB/sec)'
        write(*,'(5x,A,F12.2,1x,A,F12.2,1x,A,F8.2,1x,A)')
     &      'Total data written :', WrtSize,  'MBytes in',   WrtTime,
     &                          'sec (',  WrtSize/WrtTime,  'MB/sec)'
#ifdef FAST_MODE
        write(*,'(5x,A,22x,F12.2,1x,A)')      'Data assembly time :',
     &                                               AssmTime, 'sec'
#endif
        write(*,'(2x,A,22x,F12.2,1x,A)')   'Output file sync time :',
     &                                               SyncTime, 'sec'

#ifdef DEL_PART_FILES
        if (del_part_files) then
          write(*,'(1x,A,22x,F12.2,1x,A)') 'Removing partial files :',
     &                          net_rmcmd_clk/dble(clk_rate), 'sec'
        endif
#endif
        nclk=3-nclk
        call system_clock (iclk(nclk), clk_rate, clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        net_gray_clk=net_gray_clk+inc_clk
        GrayTime=net_gray_clk/dble(clk_rate)
        write(*,'(14x,A,22x,F12.2,1x,A)') 'Gray time :',GrayTime,'sec'
        inc_clk=iclk(nclk)-iclk_init
        write(*,'(47x,A/12x,A,11x,F12.2,1x,A/)') '------------------',
     &   'Elapsed wall-clock time:', inc_clk/dble(clk_rate), 'sec'
      endif
#endif
      stop
      end

#endif   /* EXCLUDE_MAIN */





! The following four subroutines differ only by name and type of the
! first and the third argument, but their inner codes are identically
! the same. Use self-inclusion to generate multiple versions. 

#ifdef FAST_MODE
# if F_TYPE == 1
      subroutine assembly_byte (buff,  count, bfr_out,
     &                          start1, count1, vdims)
      implicit none
      integer(kind=1) ::  buff(*), bfr_out(*)

# elif F_TYPE == 2
      subroutine assembly_int2 (buff, count,  bfr_out,
     &                          start1, count1, vdims)
      implicit none
      integer(kind=2) :: buff(*), bfr_out(*)

# elif F_TYPE == 3   
      subroutine assembly_int  (buff, count,  bfr_out,
     &                          start1, count1, vdims)
      implicit none
      integer(kind=4) :: buff(*), bfr_out(*)

# elif F_TYPE == 4   
      subroutine assembly_real (buff, count,  bfr_out,
     &                          start1, count1, vdims)
      implicit none
      real(kind=4) :: buff(*), bfr_out(*)

# elif F_TYPE == 5   
      subroutine assembly_double (buff,count, bfr_out,
     &                          start1, count1, vdims)
      implicit none
      real(kind=8) :: buff(*), bfr_out(*)


# else
      subroutine assembly_text (buff, count,  BFr_out,
     &                          start1, count1, vdims)
      implicit none
      character(len=1) buff(*), bfr_out(*)
# endif
      integer vdims, count(vdims),  start1(vdims), count1(vdims),
     &        ndims, i, istr,imax,imax1, j,js,js1, jstr,jmax,jmax1,
     &        k,ks,ks1, kstr,kmax,kmax1, l,ls,ls1, lstr,lmax,lmax1

! Check whether the last dimension is unlimited dimension and if so,
! ignore it, since records corresponding to the different indices along
! the unlimited dimension are always written one-by-one.

      if (count1(vdims).eq.1) then
        ndims=vdims-1
      else
        ndims=vdims
      endif
      imax=count(1)                ! WARNING: This code is restricted
      istr=start1(1)               ! for partitioned arrays to have no
      imax1=count1(1)              ! more than 4 dimensions, not
      if (ndims.gt.1) then         ! counting unlimited dimension.
        jmax=count(2)
        jstr=start1(2)             ! Furthermore, it is assumed that
        jmax1=count1(2)            ! unlimited dimension is always the
      else                         ! last one [this is the standard
        jstr=1                     ! practice in ROMS, however netCDF
        jmax=1                     ! is not restricted to do so].
        jmax1=1
      endif
      if (ndims.gt.2) then
        kmax=count(3)
        kstr=start1(3)
        kmax1=count1(3)
      else
        kstr=1
        kmax=1
        kmax1=1
      endif
      if (ndims.gt.3) then
        lmax=count(4)
        lstr=start1(4)
        lmax1=count1(4)
      else
        lstr=1
        lmax=1
        lmax1=1
      endif
      if (ndims.gt.4) then
        write(*,'(/1x,2A/12x,A/)')   '### ERROR: Exceeding limit of ',
     &                           '4 dimensions for partitioned array',
     &                        '[unlimited dimension does not count].'
        stop
      endif

c*    write(*,'(1x,A,2I3,3(3x,A,4I4))')  'ndims,vdims =', ndims,
c*   &   vdims, 'imax1,jmax1,kmax1,lmax1 =', imax1,jmax1,kmax1,lmax1,
c*   &                  'imax,jmax,kmax,lmax =', imax,jmax,kmax,lmax,
c*   &                  'istr,jstr,kstr,lstr =', istr,jstr,kstr,lstr

      do l=1,lmax
        ls=l-1
        ls1=l+lstr-2
        do k=1,kmax
          ks=k-1 +ls*kmax
          ks1=k+kstr-2 + ls1*kmax1
          do j=1,jmax
            js=j-1 +ks*jmax
            js1=j+jstr-2 + ks1*jmax1
            do i=1,imax
              bfr_out(i+istr-1 + js1*imax1)=buff(i + js*imax)
            enddo
          enddo
        enddo
      enddo
      return
      end

# if F_TYPE == 1
#  undef F_TYPE
#  define F_TYPE 2
# elif F_TYPE == 2
#  undef F_TYPE
#  define F_TYPE 3
# elif F_TYPE == 3
#  undef F_TYPE
#  define F_TYPE 4
# elif F_TYPE == 4
#  undef F_TYPE
#  define F_TYPE 5
# elif F_TYPE == 5
#  undef F_TYPE
# else
#  define F_TYPE 1
# endif
# ifdef F_TYPE
#  include "ncjoin.F"
# endif
#endif /* FAST_MODE */
