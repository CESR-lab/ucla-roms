      program partit

! Generic netCDF partitioning tool: reads netCDF files corresonding to
! the whole physical grid and prepares multiple files which hold data
! corresponding to different subdomains.  These files can be then read
! in parallel by different MPI processes.
!
! Usage:  partit NP_XI NP_ETA ncname1 ... ncnameN
! ------  where NP_XI  number of subdomains along XI-direction
!               NP_ETA number of subdomains along ETA-direction
!               ncname1 ... ncnameN  names of netCDF files
!
! Non-partitionable objects of netCDF files, such as scalar variables
! and attributes (both global and attributes to variables) are copied
! redundantly into the partitioned files, while partitionable array
! data is subdivided into subdomains and distributed among the
! partitioned files in such a manner that all files contain individual
! data without any overlap or redundantly stored data.
!
! The partitioning algorithm works as follows:  The partitionable
! dimensions ('xi_rho', 'xi_u', 'eta_rho' and 'eta_v') are identified
! by their names, then then their values are read and compared in pairs
! to detect if any of the directions have periodicity.  It is assumed
! that ghost points corresponding to physical boundaries are stored in
! the file, but computational margins (including periodic margins) are
! not.  Consequently, if 'xi_rho' and 'xi_u' are equal to each other
! then XI-direction is periodic, and if they differ by one, it is not.
! ETA-direction is treated similarly. Once periodicity type is
! determined, the internal number of internal points in each
! direction (i.e. excluding ghost points corresponding to physical
! boundaries) id divided by the number of subdomains in that
! direction and then physical boundary points are attached to
! subdomains which are adjacent to the boundaries. This results in
! slightly different dimension sizes of netCDF files corresponding
! to diffeent subdomains.
!
! Once all dimensions are sorted out, data corresponding to
! subdomains is extracted from the source file and copied into
! partial files.

#define AUTORENICE
#define TIMING
c--#define VERBOSE


      implicit none
#include "netcdf.inc"
      integer, parameter ::  maxdims=32,  maxvars=64
      character(len=80) src_name,         string
      character(len=32) dimname(maxdims), varname(maxvars)

      integer narg,  NP_XI,  nnodes, xi_rho, id_xi_rho, id_xi_psi,
     &        arg,   NP_ETA, node,   xi_u,   id_xi_u,   id_xi_v,
     &               eta_rho,        id_eta_rho,        id_eta_psi,
     &        ierr,  eta_v,          id_eta_v,          id_eta_u,
     &        ndims, nvars,  ngatts, tsize,  unlimdimid,varatts,
     &        i,j,k, ncsrc,  lvar,   size,  rec, lstr, lsrc,  lenstr

      character(len=80), allocatable, dimension(:) :: ncname
      integer, allocatable, dimension(:,:) :: vid
      integer, allocatable, dimension(:) :: ncid, xi_start, xi_size,
     &                                           eta_start, eta_size
      logical, allocatable, dimension(:) ::  western_edge,
     &          eastern_edge, southern_edge, northern_edge


      integer max_buff_size,  vartype(maxvars),   vardims(maxvars),
     &        dimid(maxdims), ibuff(maxdims),dimids(maxdims,maxvars),
     &        dimsize(maxdims), start(maxdims),   count(maxdims),
     &                                            start1(maxdims)
      logical part_switch(maxvars),                 series(maxvars)
      common /partit_int/ max_buff_size,   vartype, vardims,
     &                    dimid,  ibuff,   dimids,  dimsize, start,
     &                    count,  start1,  part_switch,      series

      real*8, allocatable, dimension(:) ::  buff

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
#endif
                    ! Function "iargc" is viewed as intrinsic by
#ifdef AIX       /* ! most modern compilers and does not need to  */
      integer iargc ! be declared. IBM xlf95 is a notable exclusion.
#endif           /* ! So do 7.x and earlier versions of Intel IFC, */
                    ! but 8.x, 9.x IFORT recognize it as intrinsic.
#ifdef AUTORENICE
      integer getpid, pid            ! Sometimes it makes sense to
      character(len=32) cmd          ! run partit in the background
      pid=getpid()                   ! mode with lowered priority to
      write(cmd,'(I8)') pid          ! not interfere with running MPI
      lstr=lenstr(cmd)               ! job. This code segment catches
      cmd(11:lstr+10)=cmd(1:lstr)    ! pid of its own process and
      lstr=lstr+10                   ! executes re-nice command.
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
#endif

! Check how many arguments are given, complain about the error,
! if too few, otherwise extract NP_X and NP_E from the first two
! arguments.

      narg=iargc()
      if (narg .lt. 3) then
        write(*,'(/1x,A/32x,A/)')      'Usage of partit should be:',
     &                       'partit NP_X NP_E ncname1 ... ncnameN'
        stop
      endif
 
      call getarg(1,string)
      lstr=lenstr(string)
      NP_XI=0
      do i=1,lstr
        j=ichar(string(i:i))-48
        if (j.ge.0 .and. j.le.9) then
          NP_XI=10*NP_XI+j
        else
          write(*,'(/8x,3A/)')     '### ERROR: First argument ',
     &           string(1:lstr), ', must be an integer number.'
          stop
        endif
      enddo
 
      call getarg(2,string)
      lstr=lenstr(string)
      NP_ETA=0
      do i=1,lstr
        j=ichar(string(i:i))-48
        if (j.ge.0 .and. j.le.9) then
          NP_ETA=10*NP_ETA+j
        else
          write(*,'(/8x,3A/)')    '### ERROR: Second argument ',
     &           string(1:lstr), ', must be an integer number.'
          stop
        endif
      enddo
      nnodes=NP_XI*NP_ETA
      write(*,'(/1x,2(4x,A,I3)/)') 'NP_XI =',  NP_XI,
     &                             'NP_ETA =', NP_ETA

      allocate(ncid(0:nnodes-1))
      allocate(ncname(0:nnodes-1))
      allocate(xi_size(0:nnodes-1))
      allocate(eta_size(0:nnodes-1))
      allocate(xi_start(0:nnodes-1))
      allocate(eta_start(0:nnodes-1))
      allocate(western_edge(0:nnodes-1))
      allocate(eastern_edge(0:nnodes-1))
      allocate(southern_edge(0:nnodes-1))
      allocate(northern_edge(0:nnodes-1))
      allocate(vid(maxvars,0:nnodes-1))

      max_buff_size=16384
      allocate(buff(max_buff_size))



! Process netCDF files: open, determine if it is already
! a partitioned file, then make general inquiry. Complain
! about error if it already partitioned, or if number of
! variables and/or dimensions exceeds specified limits.

      do arg=3,narg
        call getarg(arg,src_name)
        lsrc=lenstr(src_name)
        ierr=nf_open (src_name(1:lsrc), nf_nowrite, ncsrc)
        if (ierr .eq. nf_noerr) then
          ierr=nf_inq_att (ncsrc, nf_global, 'partition', i,j)
          if (ierr .eq. nf_noerr) then
            write(*,'(/1x,4A/14x,2A/)')'### WARNING: netCDF file ''',
     &        src_name(1:lsrc),  ''' is already ', 'a partial file ',
     &       'and cannot be partitioned any further ==> ignoring it.'
            goto 97     !--> next file
          endif
        else
          write(*,'(/1x,4A/14x,A)')      '### WARNING: Cannot open ',
     &    'netCDF file ''',src_name(1:lsrc),'''.', nf_strerror(ierr)
          goto 97     !--> next file
        endif

 
        write(*,'(1x,3A)') 'Processing netCDF file ''',
     &                        src_name(1:lsrc), '''...'
 
        ierr=nf_inq (ncsrc, ndims, nvars, ngatts, unlimdimid)
        if (ierr.eq.nf_noerr) then
          if (ndims .gt. maxdims) then
            write(*,'(/1x,2A,I4,1x,4A/12x,A,I4,2A/)')  '### ERROR: ',
     &        'Number of dimensions',      ndims,       'in netCDF ',
     &        'file ''', src_name(1:lsrc), ''' exceeds the limit of',
     &        'maxdims =',     maxdims,      '. Increase parameter ',
     &                      '"maxdims" in "partit.F" and recompile.'
            ierr=ierr+1
          endif
          if (nvars .gt. maxvars) then
            write(*,'(/1x,2A,I4,1x,3A/12x,A,I4,2A/)')  '### ERROR: ',
     &           'Number of variables',  nvars,  'in netCDF file ''',
     &            src_name(1:lsrc),       ''' exceeds the limit of ',
     &           'maxvars =',     maxvars,   '. Increase parameter ',
     &                      '"maxvars" in "partit.F" and recompile.'
            ierr=ierr+1
          endif
        else
          write(*,'(/1x,4A/12x,A/)')       '### ERROR: Cannot make ',
     &                         'general inquiry into netCDF file ''',
     &                    src_name(1:lsrc), '''.', nf_strerror(ierr)
        endif
        if (ierr .ne. nf_noerr) stop

! Sort out dimensions: For each dimension find and save its name and
! size. Then check whether all partitionable dimensions (identified
! by names 'xi_rho', 'xi_u', 'eta_rho' and 'eta_v')  are present and
! save their IDs and sizes.

        tsize=1      ! <-- default value.
        do i=1,ndims
          ierr=nf_inq_dim (ncsrc, i, dimname(i), dimsize(i))
          if (ierr.eq.nf_noerr) then
            if (i.eq. unlimdimid) then
              tsize=dimsize(i)
              dimsize(i)=nf_unlimited
            endif
          else
            write(*,'(/1x,2A,I4/12x,3A/12x,A/)')       '### ERROR: ',
     &              'Cannot determine name and size for dimension #',
     &              i,  'in netCDF file ''', src_name(1:lsrc), '''.',
     &                                             nf_strerror(ierr)
             goto 97     !--> next file
          endif
        enddo

! Determine IDs and sizes of partitionable dimensions, 'xi_rho',
! 'xi_u', 'eta_rho' and 'eta_v'. Also save IDs of obsolete dimensions
! 'xi_psi', 'xi_v', 'eta_psi' and and 'eta_u'. These are used to
! readress obsolete (redundant) dimensions according to the rules:

#ifdef VERBOSE
        write(*,'(4x,A)') 'Identifying dimensions:'
#endif
                                  ! redundant dimension mapping:
        xi_rho=0                  !
        id_xi_rho=0               !          xi_psi  --> xi_u
        xi_u=0                    !
        id_xi_u=0                 !          xi_v    --> xi_rho
        eta_rho=0                 !
        id_eta_rho=0              !          eta_psi --> eta_v
        eta_v=0                   !
        id_eta_v=0                !          eta_u   --> eta_rho
        id_xi_psi=0
        id_xi_v=0
        id_eta_psi=0
        id_eta_u=0

        do i=1,ndims
          lvar=lenstr(dimname(i))
          if (lvar.eq.6 .and. dimname(i)(1:lvar).eq.'xi_rho') then
            id_xi_rho=i
            xi_rho=dimsize(i)
          elseif (lvar.eq.4 .and. dimname(i)(1:lvar).eq.'xi_u') then
            id_xi_u=i
            xi_u=dimsize(i)
          elseif (lvar.eq.7.and.dimname(i)(1:lvar).eq.'eta_rho') then
            id_eta_rho=i
            eta_rho=dimsize(i)
          elseif (lvar.eq.5 .and. dimname(i)(1:lvar).eq.'eta_v') then
            id_eta_v=i
            eta_v=dimsize(i)
          elseif (lvar.eq.6 .and.dimname(i)(1:lvar).eq.'xi_psi') then
            id_xi_psi=i
          elseif (lvar.eq.4 .and. dimname(i)(1:lvar).eq.'xi_v') then
            id_xi_v=i
          elseif (lvar.eq.7.and.dimname(i)(1:lvar).eq.'eta_psi') then
            id_eta_psi=i
          elseif (lvar.eq.5 .and. dimname(i)(1:lvar).eq.'eta_u') then
            id_eta_u=i
          endif
#ifdef VERBOSE
          write(*,'(8x,I3,1x,A,T24,I3)') i, dimname(i)(1:lvar),
     &                                              dimsize(i)
#endif
        enddo

        if (id_xi_rho.ne.0 .and.  id_xi_u.eq.0) then
          xi_u=xi_rho-1
        elseif (id_xi_rho.eq.0 .and.  id_xi_u.ne.0) then
          xi_rho=xi_u+1
        endif
        if (id_eta_rho.ne.0 .and. id_eta_v.eq.0) then
          eta_v=eta_rho-1 
        elseif (id_eta_rho.eq.0 .and. id_eta_v.eq.0) then
          eta_rho=eta_v+1
        endif
        if (xi_rho.eq.0  .or. xi_u.eq.0 .or.
     &      eta_rho.eq.0 .or. eta_v.eq.0) then
          write(*,'(/8x,2A/15x,3A/)')       '### ERROR: not all ',
     &            'partitionable dimensions are found in netCDF ',
     &                         'file ''', src_name(1:lsrc), '''.'
          goto 97     !--> next file
        endif






! Set horizontal dimensions for each subdomain
!---- ---------- ---------- --- ---- ---------

        call mpi_setup (NP_XI,NP_ETA, xi_rho,eta_rho,
     &           xi_start,xi_size, eta_start,eta_size,
     &                    western_edge, eastern_edge,
     &                   southern_edge, northern_edge)




! Create partitioned files:
!======= =========== ======

        do node=0,nnodes-1
          lsrc=lenstr(src_name)
          ncname(node)=src_name(1:lsrc)
          lstr=lsrc
          ierr=0
          call insert_node (ncname(node), lstr, node, nnodes, ierr)

          if (ierr. eq. 0) then 
            ierr=nf_create (ncname(node)(1:lstr), nf_clobber+nf_64bit_offset,
     &                                                           ncid(node))
            if (ierr .eq. nf_noerr) then
              write(*,'(4x,3A)')      'Created partitioned file ''',
     &                                  ncname(node)(1:lstr), '''.'
            else
              write(*,'(/8x,A,1x,3A/8x,A)')    '### ERROR: cannot ',
     &                'create netCDF file ''', ncname(node)(1:lstr),
     &                                    '''.',  nf_strerror(ierr)
            endif
          endif
          if (ierr. ne. 0) goto 97     !--> next file to process

! Define dimensions for the partitioned files:
!------- ---------- --- --- ----------- -----


          do i=1,ndims
            if (i .eq. id_xi_rho) then
              size=xi_size(node)
              if (western_edge(node)) size=size+1
              if (eastern_edge(node)) size=size+1

            elseif (i .eq. id_xi_u) then
              size=xi_size(node)
              if (eastern_edge(node)) size=size+1

            elseif (i .eq. id_eta_rho) then
              size=eta_size(node)
              if (southern_edge(node)) size=size+1
              if (northern_edge(node)) size=size+1

            elseif (i .eq. id_eta_v) then
              size=eta_size(node)
              if (northern_edge(node)) size=size+1
            else
              size=dimsize(i)
            endif

            dimid(i)=0
            lvar=lenstr(dimname(i))
            if (i.ne.unlimdimid .and. size.eq.0) then
              if (node.eq.0) write(*,'(4x,4A)')     'Suppressing ',
     &          'zero-size dimension ''', dimname(i)(1:lvar), '''.'
            elseif  (i.eq.unlimdimid .and. tsize.eq.0) then 
              if (node.eq.0) write(*,'(4x,4A)')     'Suppressing ',
     &                          'zero-size unlimited dimension ''',
     &                                 dimname(i)(1:lvar), '''.'
            elseif (i.eq.id_xi_psi) then
              if (node.eq.0) write(*,'(4x,2A)')     'Suppressing ',
     &                            'obsolete dimension ''xi_psi''.'
            elseif (i.eq.id_xi_v) then
              if (node.eq.0) write(*,'(4x,2A)')     'Suppressing ',
     &                              'obsolete dimension ''xi_v''.'
            elseif (i.eq.id_eta_psi) then
              if (node.eq.0) write(*,'(4x,2A)')     'Suppressing ',
     &                           'obsolete dimension ''eta_psi''.'
            elseif (i.eq.id_eta_u)  then
              if (node.eq.0) write(*,'(4x,2A)')     'Suppressing ',
     &                              'obsolete dimension ''xi_u''.'
            else
              ierr=nf_def_dim (ncid(node), dimname(i)(1:lvar),
     &                                           size, dimid(i))
#ifdef VERBOSE
              if (ierr .eq. nf_noerr) then
                write(*,'(8x,2A,I2,2x,2A,I4,T50,A,I4)') 'Defined ',
     &                  'dimension #',dimid(i), dimname(i)(1:lvar),
     &                                 ' =', size, 'node =', node
              else
#else
              if (ierr .ne. nf_noerr) then
#endif
                write(*,'(/1x,4A/12x,A)')    '### ERROR: Cannot ',
     &                  'define dimension ''', dimname(i)(1:lvar),
     &                                   '''.', nf_strerror(ierr)
              endif
            endif
          enddo

! After this moment array dimid(1:ndims) contains the set of NEW
! dimension IDs.  If the original file contains any of the four
! obsolete dimensions, 'xi_psi', 'eta_psi', 'xi_v', and 'eta_u'
! which have been eliminated, then dimid(i) does not correspond to
! the set of dimension IDs of the original file [which would be
! just dimid(i)=i]. Consequently, array dimid(1:ndims) will later
! be used later to remap old dimension IDs into new ones.
!
! Put global attributes:    The new attribute 'partition' identifies
!---- ------ -----------    identifies positon of each individual
!                           subdomain within the processor grid
          ibuff(1)=node
          ibuff(2)=nnodes
          ibuff(3)=xi_start(node)
          ibuff(4)=eta_start(node)
          ierr=nf_put_att_int (ncid(node), nf_global, 'partition',
     &                                          nf_int, 4, ibuff)
        enddo

! Copy global attributes

        do i=1,ngatts
          ierr=nf_inq_attname (ncsrc, nf_global, i, string)
          if (ierr. eq. nf_noerr) then
            lvar=lenstr(string)
            do node=0,nnodes-1
              ierr=nf_copy_att (ncsrc, nf_global, string(1:lvar),
     &                                     ncid(node), nf_global)
              if (ierr. ne. nf_noerr) then
                lstr=lenstr(ncname(node))
                write(*,'(/1x,7A/12x,A)')  '### ERROR: Cannot copy ',
     &             'global attribute ''', string(1:lvar), ''' into ',
     &             'netCDF file ''',   ncname(node)(1:lstr),   '''.',
     &                                             nf_strerror(ierr)
                goto 97
              endif
            enddo
          else
            lstr=lenstr(ncname(0))
            write(*,'(/1x,2A,I3,1x,4A/12x,A/)') '### ERROR: Cannot ',
     &       'determine name of global attribute #', i, 'in netCDF ',
     &       'file ''', src_name(1:lsrc), '''.',  nf_strerror(ierr)
            goto 97
          endif
        enddo

! Define variables and their attributes in the partitioned files.
!------- --------- --- ----- ---------- -- --- ----------- ------

        do i=1,nvars
          ierr=nf_inq_var (ncsrc,   i, varname(i),  vartype(i),
     &                     vardims(i), dimids(1,i),   varatts)
          if (ierr.ne.nf_noerr) then
             write(*,'(/1x,2A,I3/12x,3A/12x,A/)')     '### ERROR: ',
     &            'Cannot make general inquiry about variable ID =',
     &             i, 'in netCDF file ''',  src_name(1:lsrc), ''',',
     &                                            nf_strerror(ierr)

            goto 97
          endif

! Readress obsolete dimensions, if any:

          do j=1,vardims(i)
            if (dimids(j,i).eq.id_xi_psi) then
              dimids(j,i)=id_xi_u
            elseif (dimids(j,i).eq.id_xi_v) then
              dimids(j,i)=id_xi_rho
            elseif (dimids(j,i).eq.id_eta_psi) then
              dimids(j,i)=id_eta_v
            elseif (dimids(j,i).eq.id_eta_u) then
              dimids(j,i)=id_eta_rho
            endif
          enddo

! Determine whether partitionable dimensions or unlimited dimension
! are present for this variable.

          series(i)=.false.
          part_switch(i)=.false.
          do j=1,vardims(i)
            if (dimids(j,i).eq.id_xi_rho .or.
     &          dimids(j,i).eq.id_xi_u    .or.
     &          dimids(j,i).eq.id_eta_rho .or.
     &          dimids(j,i).eq.id_eta_v) then
              part_switch(i)=.true.
            elseif (dimids(j,i).eq.unlimdimid) then
              series(i)=.true.
            endif
          enddo

          if (tsize.gt.0 .or. .not.series(i)) then

! WARNING: Since dimids(1:vardims(i),i) contains dimension IDs
! corresponding to the set of IDs of the ORIGINAL file, and since
! some of the original dimensions were eliminated (merged), the
! set of dimension IDs in the NEW definitions is obtained by
! inverse mapping of dimids(j,i) onto ibuff(j) using dimid(k) as
! a mapping array.

            do j=1,vardims(i)
              do k=1,ndims
                if (dimids(j,i).eq.k) ibuff(j)=dimid(k)
              enddo
            enddo
c**         write(*,*) 'old_dimids:', (dimids(j,i),j=1,vardims(i))
c**         write(*,*) 'new_dimids:',    (ibuff(j),j=1,vardims(i))
 
            lvar=lenstr(varname(i))
            do node=0,nnodes-1
              vid(i,node)=0
              if (lvar.gt.5) then
                if (varname(i)(lvar-4:lvar).eq.'_west' .and.
     &              .not.western_edge(node)) vid(i,node)=-1000

                if (varname(i)(lvar-4:lvar).eq.'_east'  .and.
     &              .not.eastern_edge(node)) vid(i,node)=-1000
              endif
              if (lvar.gt.6) then 
                if (varname(i)(lvar-5:lvar).eq.'_south' .and.
     &              .not.southern_edge(node)) vid(i,node)=-1000

                if (varname(i)(lvar-5:lvar).eq.'_north'  .and.
     &              .not.northern_edge(node)) vid(i,node)=-1000
              endif
              if (vid(i,node).eq.0) then 
                ierr=nf_def_var (ncid(node), varname(i)(1:lvar),
     &                vartype(i), vardims(i), ibuff, vid(i,node))
                if (ierr.ne.nf_noerr) then
                  write(*,'(/1x,4A,I4/12x,A/)') '### ERROR: Cannot ',
     &                    'create variable ''',   varname(i)(1:lvar),
     &                  ''', node =',   node,    nf_strerror(ierr)
#ifdef VERBOSE
                else
                 write(*,'(I3,3A,T36,2(A,I3,1x),2(1x,A,L1),1x,A,I4)')
     &                 i, ' Created variable ''', varname(i)(1:lvar),
     &                 '''',  'varid =',   vid(i,node),  'vardims =',
     &                    vardims(i), 'part_switch=', part_switch(i),
     &                          'series=',series(i), 'node =', node
#endif
                endif

#ifdef DUMMY_BRY_VAR
              else
                ierr=nf_def_var (ncid(node), varname(i)(1:lvar),
     &                           vartype(i), 0, ibuff, vid(i,node))
# ifdef VERBOSE
                write(*,'(I3,4A,T84,A,I4)')     i, ' Created dummy ',
     &                'scalar variable ''', varname(i)(1:lvar), '''',
     &                                              'node =', node
# endif
#else
# ifdef VERBOSE
              else
                write(*,'(I3,4A,T84,A,I4)')      i, ' Suppressed ',
     &                     'variable ''', varname(i)(1:lvar), '''',
     &                                              'node =', node
# endif
#endif
              endif
            enddo
            do j=1,varatts
              ierr=nf_inq_attname (ncsrc, i, j, string)
              lvar=lenstr(string)
              do node=0,nnodes-1
                if (vid(i,node).gt.0) then
                  ierr=nf_copy_att (ncsrc, i, string(1:lvar),
     &                              ncid(node), vid(i,node))
                  if (ierr.ne.nf_noerr) then
                    write(*,'(/1x,4A,I4/12x,A/)') '### ERROR: ',
     &               'Cannot copy attribute ''', string(1:lvar),
     &               ''' for variable ''',   varname(i)(1:lvar),
     &               ''', node =',   node,    nf_strerror(ierr)

                  endif
                endif
              enddo
            enddo
          endif
        enddo    ! <--- i, loop over variables

! Leave definition mode:
!------ ---------- -----

        do node=0,nnodes-1
c**          ierr=nf_set_fill (ncid(node), nf_nofill, i)
          ierr=nf_enddef(ncid(node))
        enddo
#ifdef VERBOSE
        write(*,'(/4x,A/)') 'Left definition mode.'
#endif
#ifdef TIMING
        nclk=3-nclk
        call system_clock (iclk(nclk), clk_rate,clk_max)
        inc_clk=iclk(nclk)-iclk(3-nclk)
        net_gray_clk=net_gray_clk+inc_clk
#endif



!     **     *     ***  *******    ***  *********  ********
!      *    ***   ***   ***   ***  ***  *  ***  *  ***    *
!       *   ***   ***   ***   ***  ***     ***     ***
!       *  *** * ***    ***   **   ***     ***     ******
!        * **  * **     ******     ***     ***     ***
!        ***    ***     ***  **    ***     ***     ***    *
!         *     **      ***   ***  ***     ***     ********



! Transfer variables into newly created files.
!========= ========= ==== ===== ======= ======

        do rec=1,max(tsize,1)
          if (tsize.gt.1) write(*,'(8x,2(1x,A,I5))')
     &      'Processing record', rec, 'out of', tsize
          do i=1,nvars
#ifdef VERBOSE
            if (tsize.eq.1) then
              lvar=lenstr(varname(i))
              write(*,'(16x,I4,1x,A,I4,1x,3A)') i, 'of', nvars,
     &           'processing variable ''', varname(i)(1:lvar), '''.'
            endif
#endif
            if ((rec.eq.1 .and. .not.series(i)) .or.
     &          (series(i) .and. tsize.gt.0)) then
              if (.not.part_switch(i) .and. .not.series(i)) then

! Scalar (zero-dimensional) variables:

                if (vartype(i) .eq. nf_char) then
                  ierr=nf_get_var_text (ncsrc, i, buff)
                elseif (vartype(i) .eq. nf_byte) then
                  ierr=nf_get_var_int1   (ncsrc, i, buff)
                elseif (vartype(i) .eq. nf_short) then
                  ierr=nf_get_var_int2   (ncsrc, i, buff)
                elseif (vartype(i) .eq. nf_int) then
                  ierr=nf_get_var_int    (ncsrc, i, buff)
                elseif (vartype(i) .eq. nf_float) then
                  ierr=nf_get_var_real   (ncsrc, i, buff)
                elseif (vartype(i) .eq. nf_double) then
                  ierr=nf_get_var_double (ncsrc, i, buff)
                else
                  lvar=lenstr(varname(i))
                  write(*,'(/8x,4A/)') '### ERROR: Scalar variable ',
     &              '''', varname(i)(1:lvar), ''' has unknown type.'
                  stop 
                endif
                if (ierr .eq. nf_noerr) then
                  do node=0,nnodes-1
                    if (vid(i,node).gt.0) then
                      if (vartype(i) .eq. nf_char) then
                       ierr=nf_put_var_text (ncid(node), vid(i,node),
     &                                                         buff)
                      elseif (vartype(i) .eq. nf_byte) then
                       ierr=nf_put_var_int1 (ncid(node), vid(i,node),
     &                                                         buff)
                      elseif (vartype(i) .eq. nf_short) then
                       ierr=nf_put_var_int2 (ncid(node), vid(i,node),
     &                                                         buff)
                      elseif (vartype(i) .eq. nf_int) then
                       ierr=nf_put_var_int  (ncid(node), vid(i,node),
     &                                                         buff)
                      elseif (vartype(i) .eq. nf_float) then
                       ierr=nf_put_var_real (ncid(node), vid(i,node),
     &                                                         buff)
                      elseif (vartype(i) .eq. nf_double) then
                       ierr=nf_put_var_double(ncid(node),vid(i,node),
     &                                                         buff)
                      endif
                      if (ierr .ne. nf_noerr) then
                        lvar=lenstr(varname(i))
                        lstr=lenstr(ncname(node))
                        write(*,'(/1x,3A/12x,3A/12x,A/)')
     &                  '### ERROR: Cannot write scalar variable ''',
     &                     varname(i)(1:lvar),     ''' into netCDF ',
     &                  'file ''',    ncname(node)(1:lstr),    '''.',
     &                                             nf_strerror(ierr)
                        goto 97
                      endif
                    endif
                  enddo
                else
                  lvar=lenstr(varname(i))
                  write(*,'(/1x,3A/12x,3A/12x,A/)')      
     &                  '### ERROR: Cannot read scalar variable ''',
     &                   varname(i)(1:lvar),      ''' from netCDF ',
     &                  'file ''',     src_name(1:lsrc),      '''.',
     &                                            nf_strerror(ierr)
                  goto 97
                endif
              elseif (.not.part_switch(i)) then

! Non-partitionable array.

                size=1
                do j=1,vardims(i)
                  if (dimids(j,i).eq.unlimdimid) then
                    start(j)=rec
                    count(j)=1
                  else
                    start(j)=1
                    count(j)=dimsize(dimids(j,i))
                  endif
                  size=size*count(j)
                enddo
                if (vartype(i) .eq. nf_char .or.
     &              vartype(i) .eq. nf_byte) then
                  size=size*1
                elseif (vartype(i) .eq. nf_short) then
                  size=size*2
                elseif (vartype(i) .eq. nf_int) then
                  size=size*4
                elseif (vartype(i) .eq. nf_float) then
                  size=size*4
                elseif (vartype(i) .eq. nf_double) then
                  size=size*8
                else
                  lvar=lenstr(varname(i))
                  write(*,'(/8x,3A/)')   '### ERROR: variable ''',
     &                 varname(i)(1:lvar), ''' has unknown type.'
                  goto 97
                endif
                if (size .gt. 8*max_buff_size) then

                  if (allocated(buff)) deallocate(buff)
                  max_buff_size=(size+7)/8
                  allocate (buff(max_buff_size))

#ifdef VERBOSE
                  write(*,*) 'Allocated buffer size =', max_buff_size
#endif
                endif
 
                if (vartype(i) .eq. nf_char) then
                  ierr=nf_get_vara_text (ncsrc, i, start,
     &                                        count, buff)
                elseif (vartype(i) .eq. nf_byte) then
                  ierr=nf_get_vara_int1  (ncsrc, i, start,
     &                                         count, buff)
                elseif (vartype(i) .eq. nf_short) then
                  ierr=nf_get_vara_int2  (ncsrc, i, start,
     &                                         count, buff)
                elseif (vartype(i) .eq. nf_int) then
                  ierr=nf_get_vara_int    (ncsrc, i, start,
     &                                         count, buff)
                elseif (vartype(i) .eq. nf_float) then
                  ierr=nf_get_vara_real   (ncsrc, i, start,
     &                                         count, buff)
                elseif (vartype(i) .eq. nf_double) then
                  ierr=nf_get_vara_double (ncsrc, i, start,
     &                                         count, buff)
                endif
                if (ierr .eq. nf_noerr) then
                  do node=0,nnodes-1
                    if (vid(i,node).gt.0) then
                      if (vartype(i) .eq. nf_char) then
                        ierr=nf_put_vara_text (ncid(node),
     &                       vid(i,node), start, count, buff)
                      elseif (vartype(i) .eq. nf_byte) then
                        ierr=nf_put_vara_int1 (ncid(node),
     &                       vid(i,node), start, count, buff)
                      elseif (vartype(i) .eq. nf_short) then
                        ierr=nf_put_vara_int2 (ncid(node),
     &                       vid(i,node), start, count, buff)
                      elseif (vartype(i) .eq. nf_int) then
                        ierr=nf_put_vara_int  (ncid(node),
     &                       vid(i,node), start, count, buff)
                      elseif (vartype(i) .eq. nf_float) then
                        ierr=nf_put_vara_real (ncid(node),
     &                       vid(i,node), start, count, buff)
                      elseif (vartype(i) .eq. nf_double) then
                        ierr=nf_put_vara_double (ncid(node),
     &                       vid(i,node), start, count, buff)
                      endif
                      if (ierr .ne. nf_noerr) then
                        lvar=lenstr(varname(i))
                        lstr=lenstr(ncname(node))
                        write(*,'(/1x,A,I4,1x,3A/12x,3A/12x,A/)')
     &                    '### ERROR: Cannot write time record =',
     &                     rec,   'for nonpartitionable array ''',
     &                     varname(i)(1:lvar),   '''into netCDF ',
     &                    'file ''', ncname(node)(1:lstr),  '''.',
     &                                          nf_strerror(ierr)
                        goto 97
                      endif
                    endif
                  enddo
                else
                  lvar=lenstr(varname(i))
                  write(*,'(/1x,A,I4,1x,3A/12x,3A/12x,A/)')
     &                    '### ERROR: Cannot read time record =',
     &                     rec,  'for nonpartitionable array ''',
     &                     varname(i)(1:lvar),  ''' from netCDF',
     &                    'file ''',   src_name(1:lsrc),   '''.',
     &                                         nf_strerror(ierr)
                  goto 97
                endif
              elseif (part_switch(i)) then

! Partitioned array:

                do node=0,nnodes-1
                  if (vid(i,node).gt.0) then
#ifdef VERBOSE
                    lvar=lenstr(varname(i))
                    write(*,'(24x,A,I4,4x,3A,I3)')  'part var id =',
     &                           i, 'name = ''', varname(i)(1:lvar),
     &                                           ''' node =', node
#endif
                    size=1
                    do j=1,vardims(i)
                      start1(j)=1
                      if (dimids(j,i).eq.id_xi_rho) then
                        start(j)=xi_start(node)
                        count(j)=xi_size(node)
                        if (western_edge(node)) then
                          count(j)=count(j)+1
                        endif
                        if (eastern_edge(node)) then
                          count(j)=count(j)+1
                        endif
 
                      elseif (dimids(j,i).eq.id_xi_u) then
                        start(j)=xi_start(node)
                        count(j)=xi_size(node)

                        if (.not.western_edge(node)) then
                          start(j)=start(j)-1
                        endif
                        if (eastern_edge(node)) then
                          count(j)=count(j)+1
                        endif
 
                      elseif (dimids(j,i).eq.id_eta_rho) then
                        start(j)=eta_start(node)
                        count(j)=eta_size(node)
                        if (southern_edge(node)) then
                          count(j)=count(j)+1
                        endif
                        if (northern_edge(node)) then
                          count(j)=count(j)+1
                        endif
 
                      elseif (dimids(j,i).eq.id_eta_v) then
                        start(j)=eta_start(node)
                        count(j)=eta_size(node)
                        if (.not.southern_edge(node)) then
                          start(j)=start(j)-1
                        endif
                        if (northern_edge(node)) then
                          count(j)=count(j)+1
                        endif
 
                      elseif (dimids(j,i).eq.unlimdimid) then
                        start(j)=rec
                        count(j)=1
                        start1(j)=rec
                       else
                        start(j)=1
                        count(j)=dimsize(dimids(j,i))

                      endif
                      size=size*count(j)
                    enddo
c**               write(*,*) 'dimids:', (dimids(j,i),j=1,vardims(i))
c**               write(*,*) ' start:',    (start(j),j=1,vardims(i))
c**               write(*,*) ' count:',    (count(j),j=1,vardims(i))
 
 
                    if (vartype(i) .eq. nf_char .or.
     &                  vartype(i) .eq. nf_byte) then
                      size=size*1
                    elseif (vartype(i) .eq. nf_short) then
                      size=size*2
                    elseif (vartype(i) .eq. nf_int) then
                      size=size*4
                    elseif (vartype(i) .eq. nf_float) then
                      size=size*4
                    elseif (vartype(i) .eq. nf_double) then
                      size=size*8
                    else
                      lvar=lenstr(varname(i))
                      write(*,'(/1x,4A/)') '### ERROR: variable ''',
     &                   varname(i)(1:lvar), ''' has unknown type.'
                      goto 97
                    endif

                    if (size .gt. 8*max_buff_size) then
                      if (allocated(buff)) deallocate(buff)

                      max_buff_size=(size+7)/8
                      allocate (buff(max_buff_size))
#ifdef VERBOSE
                 write(*,*) 'Allocated buffer size =',max_buff_size
#endif
                    endif
 
                    if (vartype(i) .eq. nf_char) then
                      ierr=nf_get_vara_text (ncsrc, i, start,
     &                                            count, buff)
                    elseif (vartype(i) .eq. nf_byte) then
                      ierr=nf_get_vara_int1   (ncsrc, i, start,
     &                                           count, buff)
                    elseif (vartype(i) .eq. nf_short) then
                      ierr=nf_get_vara_int2   (ncsrc, i, start,
     &                                           count, buff)
                    elseif (vartype(i) .eq. nf_int) then
                      ierr=nf_get_vara_int    (ncsrc, i, start,
     &                                           count, buff)
                    elseif (vartype(i) .eq. nf_float) then
                      ierr=nf_get_vara_real   (ncsrc, i, start,
     &                                            count, buff)
                    elseif (vartype(i) .eq. nf_double) then
                      ierr=nf_get_vara_double (ncsrc, i, start,
     &                                           count, buff)
                    endif
 
                    if (ierr .eq. nf_noerr) then
                      if (vartype(i) .eq. nf_char) then
                        ierr=nf_put_vara_text (ncid(node),
     &                       vid(i,node), start1, count, buff)
                      elseif (vartype(i) .eq. nf_byte) then
                        ierr=nf_put_vara_int1 (ncid(node),
     &                       vid(i,node), start1, count, buff)
                      elseif (vartype(i) .eq. nf_short) then
                        ierr=nf_put_vara_int2 (ncid(node),
     &                       vid(i,node), start1, count, buff)
                      elseif (vartype(i) .eq. nf_int) then
                        ierr=nf_put_vara_int  (ncid(node),
     &                       vid(i,node), start1, count, buff)
                      elseif (vartype(i) .eq. nf_float) then
                        ierr=nf_put_vara_real (ncid(node),
     &                       vid(i,node), start1, count, buff)
                      elseif (vartype(i) .eq. nf_double) then
                        ierr=nf_put_vara_double (ncid(node),
     &                       vid(i,node), start1, count, buff)
                      endif
                      if (ierr .ne. nf_noerr) then
                        lvar=lenstr(varname(i))
                        lstr=lenstr(ncname(node))
                        write(*,'(/1x,A,I4,1x,3A/12x,3A/12x,A)')
     &                 '### ERROR: Cannot write time record =',  rec,
     &                 'of partitioned array ''', varname(i)(1:lvar),
     &                 '''',  'into netCDF file ''',
     &                 ncname(node)(1:lstr), '''.', nf_strerror(ierr)
                        goto 97
                      endif
                    else
                      lvar=lenstr(varname(i))
                      write(*,'(/1x,A,I4,1x,3A/12x,3A/12x,A)')
     &                 '### ERROR: Cannot read time record =',   rec,
     &                 'of partitioned array ''', varname(i)(1:lvar),
     &                 '''', 'from netCDF file ''', src_name(1:lsrc),
     &                                     '''.',  nf_strerror(ierr)
                      write(*,*) 'start =', start
                      write(*,*) 'count =', count
                      goto 97
                    endif
                  endif
                enddo       ! <-- node=0,nnodes-1
              endif       ! <-- part_switch
            endif       ! <--series(i) .or. rec.eq.1
          enddo       ! <-- i=1,nvars
        enddo       ! <-- rec=1,tsize

! Close all netCDF files

  97    write(*,*) 'closing files...'
        ierr=nf_close (ncsrc)
        write(*,*) '...........input'
        do node=0,nnodes-1
         ierr=nf_close (ncid(node))
        enddo
        write(*,*) '..........output'
      enddo       ! <-- arg=3,narg

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



! Setup horizontal dimensions and associated variables for each
! subdomain. The following code is extracted into a separate entity,
! which is mathematically consistent with the actual mpi_setup, and
! compute_starts_and_counts.h.  In essense, this subroutine receives 
! four integer numbers: dimension of the grid (including boundary
! rows on the side) -- xi_rho,eta_rho; and number of partitions in
! each direction, NP_XI,NP_ETA.  These four are translated into
! dimensions of subdomains, Lm,Mm, bounds of used portions of arrays
! iwest,ieast,jsouth,jnorth and global-to-relative index shift
! translations, SW_corn,jSW_corn (exactly the same way as in
! mpi_setup.F), which are then further translated into
! xi_start(node),xi_size(node), and eta_start(node),eta_size(node),
! which have meaning of starting netCDF indices for RHO-point sub-
! array in netCDF file belonging to each individual subdomain, and
! the sizes of subarrays.  The other four variables, western_,
! eastern_, southern_, and northern_edge are logical flags to
! identify the proximity of side boundary on each side for each
! subdomain. 

! Note that the code above --- the main part of "partit" --- is
! written in such a way that it makes no assumption about the
! structure of "processor grid" (the arrangement of subdomains
! corresponding to MPI nodes relatively to the physical grid), but
! rather relies exclussively on the eight variables defined in the
! code below. 


      subroutine mpi_setup (NP_XI,NP_ETA, xi_rho,eta_rho,
     &               xi_start,xi_size, eta_start,eta_size,
     &                        western_edge, eastern_edge,
     &                       southern_edge, northern_edge)
      implicit none
      integer, intent(in) :: NP_XI,NP_ETA, xi_rho,eta_rho
! out
      integer, intent(out), dimension(0:NP_XI*NP_ETA-1) ::
     &              xi_start, xi_size, eta_start, eta_size
      logical, intent(out), dimension(0:NP_XI*NP_ETA-1) ::
     & western_edge, eastern_edge, southern_edge, northern_edge
! internal
      integer LLm,Lm, MMm,Mm, nnodes,node, inode,jnode,
     &        iwest,ieast,jsouth,jnorth, iSW_corn,jSW_corn,
     &                                   off_XI,off_ETA


      nnodes=NP_XI*NP_ETA 
      LLm=xi_rho-2                   !
      MMm=eta_rho-2                  !
      Lm=(LLm+NP_XI-1)/NP_XI         !
      Mm=(MMm+NP_ETA-1)/NP_ETA       !
#ifdef VERBOSE
        write(*,'(4x,2A,2I5,4x,A,2I5)') 'mpi_setup: found grid ',
     &          'sizes: LLm,MMm =',  LLm,MMm,  'Lm,Mm =', Lm,Mm
#endif


      do node=0,nnodes-1
        jnode=node/NP_XI             ! the following segment maps
        inode=node-jnode*NP_XI       ! exactly onto "mpi_setup" in
                                     ! the actual ROMS code
        off_XI=NP_XI*Lm-LLm
        iSW_corn=inode*Lm-off_XI/2
        if (inode.eq.0) then
          iwest=1+off_XI/2
        else
          iwest=1
        endif
        if (inode.lt.NP_XI-1) then
          ieast=Lm
        else
          ieast=Lm -(off_XI+1)/2
        endif

        off_ETA=NP_ETA*Mm-MMm
        jSW_corn=jnode*Mm-off_ETA/2
        if (jnode.eq.0) then
          jsouth=1+off_ETA/2
        else
          jsouth=1
        endif
        if (jnode.lt.NP_ETA-1) then
          jnorth=Mm
        else
          jnorth=Mm -(off_ETA+1)/2
        endif

        xi_size(node)=ieast - iwest+1
        eta_size(node)=jnorth - jsouth+1

        if (inode.eq.0) then
          xi_start(node)=iSW_corn+iwest
        else
          xi_start(node)=iSW_corn+iwest+1
        endif

        if (jnode.eq.0) then
          eta_start(node)=jSW_corn+jsouth
        else
          eta_start(node)=jSW_corn+jsouth+1
        endif


        if (inode.eq.0) then
          western_edge(node)=.true.
        else
          western_edge(node)=.false.
        endif
        if (inode.lt.NP_XI-1) then
          eastern_edge(node)=.false. 
        else
          eastern_edge(node)=.true.
        endif
        if (jnode.eq.0) then
          southern_edge(node)=.true.
        else
          southern_edge(node)=.false.
        endif
        if (jnode.lt.NP_ETA-1) then
          northern_edge(node)=.false. 
        else
          northern_edge(node)=.true.
        endif

#ifdef VERBOSE
        write(*,'(A,3I5,2x,A,4(1x,L1))') 'node =', node, inode,jnode,
     &                    'western,eastern,southern,northern_edge =',
     &                      western_edge(node),   eastern_edge(node),
     &                      southern_edge(node), northern_edge(node)
#endif
      enddo   !--> discard iwest,ieast,jsouth,jnorth
      return
      end 
