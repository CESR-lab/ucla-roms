      module bry_extract_child   ! Interpolate and extract online the boundary files
                                 ! needed to run child (nested) simulation

      ![ INFO:
      !
      ! coded by: Devin Dollery & Jeroen Molemaker - 2021 March
      !
      ! STEPS:
      !
      ! 1) READ CHILD BOUNDARY LOCATIONS:
      !
      !    subroutine assign_child_bry_to_parent
      !
      !    The location of the child boundary is stored in one full-domain file
      !    and is read in by every mpi node.
      !    The child boundary file consists of fractional i'th and j'th positions of the
      !    child boundary points relative to the parent's index domain.
      !    The i/j value is -1 over a land mask and -2 if not in parent domain.
      !    Each i'th/j'th data consists of 2 rows per boundary side. They are orientated
      !    as per child grid. I.e. west_bry row 1 is western most. north_bry row 1 is south side.
      !
      !    Seen in DIAGRAM (A) below, each child point (a-f) will have a input variable
      !    giving it's i'th position, i.e. variable 'i_west', which for point b = 1.4.
      !    Another variable will hold the j'th position, i.e. 'j_west' which for point b = 1.2.
      !    If the boundary is closed then i and j variables will not exist for that boundary.
      !
      !         DIAGRAM (A)
      !
      !   C  ^                ^  D
      ! (1,2)|                |(2,2)
      !      X---o---o--------X-->         A-D = parent grid points (X)
      !      |   e   f        |
      ! N    |                |            a-f = child boundary grid points (West boundary)
      ! ^    |   o   o        |
      ! |    |   c   d        |            Fractional i/j coordinates of child
      ! |    |                |            relative to parent grid:
      !      |   o   o        |              e = (1.2, 2.0)  f = (1.4, 2.0)
      !      |   a   b        |              c = (1.2, 1.6)  d = (1.4, 1.6)
      !      X----------------X-->           a = (1.2, 1.2)  b = (1.4, 1.2)
      ! (1,1)                  (2,1)
      !   A                      B
      !
      !
      ! 2) ASSIGN CHILD BRY TO PARENT NODE:
      !
      !    subroutine assign_child_bry_to_parent
      !
      !    As mentioned each parent node will read in the full child boundary data.
      !    It then needs to calculate which (if any) of the child boundary fall within
      !    the parent node's range.
      !    Since the child bry is given in i/j relative to the parent's global grid points,
      !    one can simply compare the parent's GP's against the child's i/j number.
      !    However, the child i/j values are relative to the global parent domain.
      !    Typically roms uses 1:Lm (istr:iend) for local xi ranges, and 1:Mm for local
      !    eta ranges of a parent node.  These are local to the parent node and
      !    thus you can't directly compare the child bry i/j to (istr:iend).
      !
      !    We need the global start and end i/j of the parent node, stored in variables:
      !    Compare child values to parent global value using glob_istr, glob_iend.
      !    Note, glob_istr=0 for western edge and 1 for rest. This is to catch child
      !    points in the ghost cell.
      !
      !    The northern grid points of a sub-domain might have j=50, then the sub-domain
      !    to the north will start with j=51. If the child grid point has a j value
      !    of 50.41 in terms of the parent domain, we decide to use the midpoint to
      !    chose which parent sub-domain to assign the child point to. For this example
      !    the midpoint is 50.5. So the point j=50.41 would fall in the sub-domain to
      !    the south.
      !
      !    For further explanation of local and global grid point numbering see:
      !    Documentation/readme-specifics/README.grid.2021.03
      !
      !
      ! 3) INTERPOLATE CHILD BRY AT REQUIRED OUTPUT RATE:
      !
      !    When required the child bry is interpolated using the surrounding 4 parent
      !    grid points.
      !    The 4 parent grid points are calculated on the fly based on the floor and
      !    ceiling of the fractional i/j child numbers. I.e. child point b is interpolated
      !    from parent points A-D in DIAGRAM (A) using bilinear interpolation.
      !
      !    However, since the parent's variable arrays use local indexing,
      !    e.g. i=1:Lm and j=1:Mm, we need to convert the child's global i/j position
      !    to a local i/j position in the parent node. This way we can access the
      !    correct entries in the local parent variable array to interpolate from.
      !
      !    Need to catch the fact that glob_istr=0 on western edge and 1 everywhere else.
      !
      !
      ! 4) WRITE RESULTS TO FILE:
      !
      !    subroutine wrt_child_bry
      !
      !
      ! 5) JOIN PARTITIONED RESULTS FILES:
      !
      !    This probably requires its own tool to join.
      !
      !
      ! X) ADDITIONAL INFO:
      !
      !    The ascii diagrams so far ignore that fact that nodes so not actually share
      !    a boundary. If you zoom in on an intersection of 4 node corners, it looks like:
      !
      !      C    |    |    D
      !           | e  |
      !           |    |
      ! 41  ------      ------
      !                    d
      !        a    b
      !
      ! 40  ------      ------
      !           |    |
      !           |  c |
      !      A    |    |    B
      !           50   51
      !
      !    We find that the child boundary might have an i value between 50-51
      !    and a j value between 40-41.
      !    In such case we use the mid point (i=40.5 and j=50.5) to divide the parent
      !    node ranges. Must use <= for one and > for other to catch exactly on mid point.
      !    For example, store:
      !    child point a, b -> parent node A
      !    child point    c -> parent node B
      !    child point    e -> parent node C
      !    child point    d -> parent node D
      !
      !
      !
      !]

#include "cppdefs.h"
#ifdef BRY_EXTRACT_CHILD

      use dimensions  ! has dimensions, a list of forcing files, and rmask, glob_istr, etc
      use read_write
      use netcdf
      use mpi         ! debug temporary for mpi abort
      use tracers, only: t, t_vname, t_lname, t_units  ! need to get names of tracers

      implicit none
      private

#include "ocean2d.h"
! scalars.h for 'knew' var
#include "scalars.h"

      ! ****************************************************************
      ! USER INPUTS
      integer :: rec_rate      = 1     ! output record rate
      integer :: recs_per_file = 10    ! total recs per file
      ! END USER INPUTS
      ! ****************************************************************

      ! DECLARE VARIABLES: ![

      character(len=9)               :: child_bry_name=' '                     ! root name, append e/w/s/n
      character(len=50)              :: c_bry_file='../sample_child_ij_bry.nc' ! child boundary input file
      character(len=5), dimension(4) :: bry_name = (/ 'south', 'west',
     &                                                'north', 'east'  /)

      ! INPUT VARS:
        ! Would be cleaner to have variables 1-4 for W/S/E/N rather than loads of variables.

                                                                             ! CHANGE: since i/j same length can be one var e.g. ij_west
      real, dimension(:,:), allocatable :: i_west, i_south, i_east, i_north  ! i position of child bry in parent grid
      real, dimension(:,:), allocatable :: j_west, j_south, j_east, j_north  ! j position of child bry in parent grid

      integer, dimension(2,4)           :: rho_tot=0                         ! Number of child points within parent on each boundary. 4 = s/w/n/e, 2 = 2 rows of i/j input
     &
      integer, dimension(2,4)           :: rho_c_str=0                       ! start indices of child boundary within parent
     &

      integer :: glob_c_s, glob_c_w, glob_c_n, glob_c_e                      ! Total number of global child points along child boundary

      ! INTERPOLATION VARS:
      real,    dimension(:,:), allocatable :: vwest, vsouth, veast, vnorth   ! array storing interpolated variable (2D or 3D)
                                                                             ! allocated for each var so only need one.
      real,    dimension(:,:), allocatable :: west_r_cffs, south_r_cffs,     ! 4 interpolation coefficients to interp from 4 parent points surronding child point.
     &                                        east_r_cffs, north_r_cffs      ! i=each bry point. j=1-4 (4 interpolants). Will need for u_cff, v_cff. 4 interpolation coefficients
      integer, dimension(:), allocatable   :: iw_sw, jw_sw, is_sw, js_sw,    ! sw parent interpolant i and j position
     &                                        ie_sw, je_sw, in_sw, jn_sw     ! could merge with coefficients.
      real                                 :: rt2                            ! square root of 2 stored for efficiency

      ! NETCDF:
      integer            :: ncid=-1                     ! in file for child bry positions
      integer            :: ncid_out                    ! output file
      integer            :: ierr=0                      ! debug temporary
      integer            :: dimid_t                     ! time dimension id in output file
      integer            :: output_rec = 0              ! record number of output. 0 indicates we need new file
      integer            :: total_output_rec = 0        ! total flux output recordings so far
      character (len=60) :: file_name = '_child_bry.nc' ! suffix to append to simulation output rootname

      public set_child_bry

      ! auxiliary tool to join partial files created here: Tools-Roms/join_child_bry.F

      !] <- end declare variables

      contains

! ----------------------------------------------------------------------
      subroutine set_child_bry  ![
      implicit none             !  Interpolate child boundary at desired intervals.

      ! local
      integer :: itrc  ! tracer index

      ! Use this to only set at required output rate. Rather have it here than in main
      ! to keep main less cluttered? Although wrt_his has it in main? Ask JM
      ! only plot if timestep is divisble by record rate
!      if (mod(iic-1,rec_rate) == 0) then ! 1 = ntstart

      if (ncid==-1) then                               ! first time only
        call read_child_bry                            ! not per variable, just r,u,v types
        call calc_child_interpolation_coeffs
      endif


      if(output_rec == 0) then                         ! if file doesn't yet exist or need new file
        call create_child_output_file
      endif

      ierr=nf90_open(file_name, nf90_write, ncid_out)
      output_rec = output_rec + 1                      ! current file record
      total_output_rec = total_output_rec + 1          ! total number of recs written

      ! zeta:
      call interpolate_child_bry( zeta(:,:,knew), 1 )  ! , rp_var var, nz, grd_type
      call wrt_child_bry( indxZ, 1, .false. )

      ! ubar:
      ! vbar:
      ! u:
      ! v:

      ! tracers:
      do itrc=1,NT
        call interpolate_child_bry( t(START_2D_ARRAY,1,nstp,itrc), N )
        call wrt_child_bry( itrc, N, .true. )
      enddo


      ierr=nf90_close(ncid_out)
      if (output_rec == recs_per_file) then  ! Reset output_rec to create new file
        output_rec = 0
      endif

!      endif ! <- mod(iic-1,rec_rate)

      end subroutine set_child_bry  !]

! ----------------------------------------------------------------------
      subroutine read_child_bry  ![
      ! Assign 4 parent GP from which to interpolate 1 child BP. Done only once.

      implicit none

      ! Do for 4 sides, if it has i then it must have j:
      call assign_child_bry_to_parent(bry_name(1), i_south, j_south, rho_tot(:,1),
     &                                rho_c_str(:,1), glob_c_s )
      call MPI_Barrier(MPI_COMM_WORLD, ierr)  ! debug

      call assign_child_bry_to_parent(bry_name(2), i_west,  j_west,  rho_tot(:,2),
     &                                rho_c_str(:,2) , glob_c_w )
      call MPI_Barrier(MPI_COMM_WORLD, ierr)  ! debug

      call assign_child_bry_to_parent(bry_name(3), i_north, j_north, rho_tot(:,3),
     &                                rho_c_str(:,3), glob_c_n )
      call MPI_Barrier(MPI_COMM_WORLD, ierr)  ! debug

      call assign_child_bry_to_parent(bry_name(4), i_east,  j_east,  rho_tot(:,4),
     &                                rho_c_str(:,4) , glob_c_e )
      call MPI_Barrier(MPI_COMM_WORLD, ierr)  ! debug

      end subroutine read_child_bry  !]

! ----------------------------------------------------------------------
      subroutine assign_child_bry_to_parent( side_name,   i_bry, j_bry,     ![
     &                                       side_totals, c_str, dim_size )
      ! Read in N/W/E/S arrays of child BP's. (in i,j coordinates)

      ! All procs read the child bc file and determine if and what work they need to do.
      ! If work, determine range along the boundary  Loop through all child BP's and
      ! compare to corners of parent node.

      implicit none

      ! input
      character(len=*),                  intent(in)    :: side_name
      real, dimension(:,:), allocatable, intent(out)   :: i_bry, j_bry  ! Final lists local to parent.
      integer, dimension(2),             intent(inout) :: side_totals   ! Number of child points within parent
      integer, dimension(2),             intent(out)   :: c_str         ! full child boundary start index falling in parent
      integer,                           intent(out)   :: dim_size      ! Size of global child boundary

      ! local
      integer               :: varid, ij, l, counter
      character(len=30)     :: vname
      integer, dimension(2) :: dimids

      real, dimension(:,:), allocatable :: i_in_bry, j_in_bry ! Store read in child boundary

      if (ncid==-1) then  ! File not yet opened.
        ierr=nf90_open(c_bry_file, nf90_nowrite, ncid) ! var will later be in grdname file
        if(ierr/=0) call handle_ierr(ierr,
     &             'read_child_bry:: Error opening file with child bry!')
      endif


      vname='i_' / / side_name
!      if(mynode==0) print *, 'child bry name: ', vname
      ierr=nf90_inq_varid(ncid, vname, varid)          ! Check if contains bry
      if (ierr==nf90_noerr) then
!        if(mynode==0) write (*,*) 'Found child boundary in file for: ', vname

        ierr=nf90_inquire_variable(ncid, varid, dimids = dimids)
        ierr=nf90_inquire_dimension(ncid, dimids(1), len=dim_size)
!        if(mynode==0) write (*,*) 'dim_size=', dim_size

        allocate(i_in_bry(dim_size,2))                 ! Allocate temporary array for all bry data.
        ierr=nf90_get_var(ncid, varid, i_in_bry)       ! read in data
        if(ierr/=0) print *, 'bad read of bry', mynode

        allocate(j_in_bry(dim_size,2))
        vname='j_' / / side_name
        ierr=nf90_inq_varid(ncid, vname, varid)
        ierr=nf90_get_var(ncid, varid, j_in_bry)

                            ! CHECK IF CHILD CONTAINED WITHIN PARENT NODE:
        do l=1,2            ! Need to check both rows as might slightly differ e.g. i=9.95 vs i=10.05
          counter=0
          do ij=1,dim_size
                                                               ! Note for boundary edges, western edge glob_istr=0
                                                               ! eatern edge and glob_iend=LLm+1. This catches child points in the ghost cell.
            if (glob_istr-0.5  <= i_in_bry(ij,l) .and.         ! Child point within parent node
     &          glob_jstr-0.5  <= j_in_bry(ij,l) .and.
     &          i_in_bry(ij,l) <  glob_iend+0.5  .and.         ! Catch +/- 0.5 to catch between parent nodes
     &          j_in_bry(ij,l) <  glob_jend+0.5        ) then  ! <= and > to catch child exactly on mid-points between parents

              counter=counter+1

            endif

            if (counter==1) c_str(l)=ij                        ! Store starting child index

          enddo
          side_totals(l)=counter
        enddo

        if (side_totals(1) > 0 .or. side_totals(2) > 0) then
          allocate(i_bry(maxval(side_totals),2))               ! Store both rows/columns in one
          allocate(j_bry(maxval(side_totals),2))               ! array using larger of two

          i_bry(1:side_totals(1),1)=i_in_bry(c_str(1):c_str(1)+side_totals(1)-1,1) ! -1 for closed interval counting
          i_bry(1:side_totals(2),2)=i_in_bry(c_str(2):c_str(2)+side_totals(2)-1,2)
          j_bry(1:side_totals(1),1)=j_in_bry(c_str(1):c_str(1)+side_totals(1)-1,1)
          j_bry(1:side_totals(2),2)=j_in_bry(c_str(2):c_str(2)+side_totals(2)-1,2)

!          write(*,'(/1x,I2,A)') mynode, ' Child points within parent node' ! debug
!          print *, mynode, 'side_totals=', side_totals               ! debug
!          print *, mynode, 'c_str(1)=',c_str(1),'c_str(2)=',c_str(2) ! debug
!          print *, mynode, 'i_bry(str,1)= ', i_bry(1,1)              ! debug
!          print *, mynode, 'i_bry(end,1)= ', i_bry(side_totals(1),1) ! debug
!          print *, mynode, 'j_bry(str,1)= ', j_bry(1,1)              ! debug
!          print *, mynode, 'j_bry(end,1)= ', j_bry(side_totals(1),1) ! debug
        endif

      else     ! Side not contained
        if(mynode==0)
     &    write (*,'(/1x,3A)') 'No child boundary in parent domain ',
     &                         'for: ', side_name
        dim_size=0
      endif

      end subroutine assign_child_bry_to_parent  !]

! ----------------------------------------------------------------------
      subroutine calc_child_interpolation_coeffs  ![

      ! Get interpolation coefficients from fraction of i/j positons.
      ! Calculate interpolation coefficients once only.

      implicit none

      integer :: ij

      if (allocated(i_south)) then                                       ! 2 rows are orientated as per child grid. i.e. row 1 = west/south edge
        call calc_bry_coeffs(rho_tot(1,1), i_south(:,1), j_south(:,1),
     &                       south_r_cffs, is_sw, js_sw)                 ! index 1 since uses southern ghost points.
      endif
      if (allocated(i_west))  then
        call calc_bry_coeffs(rho_tot(1,2),  i_west(:,1),  j_west(:,1),
     &                       west_r_cffs,  iw_sw, jw_sw)                 ! index 1 since uses wester   ghost points.
      endif
      if (allocated(i_north)) then
        call calc_bry_coeffs(rho_tot(2,3), i_north(:,2), j_north(:,2),
     &                       north_r_cffs, in_sw, jn_sw)                 ! index 2 since uses northern ghost points.
      endif
      if (allocated(i_east))  then
        call calc_bry_coeffs(rho_tot(2,4),  i_east(:,2),  j_east(:,2),
     &                       east_r_cffs,  ie_sw, je_sw)                 ! index 2 since uses eastern  ghost points.
      endif

      rt2=sqrt(2.0)  ! set interpolation coefficient once

      end subroutine calc_child_interpolation_coeffs  !]

! ----------------------------------------------------------------------
      subroutine calc_bry_coeffs(bry_tot, i_bry, j_bry, rho_cffs, i_sw, j_sw)  ![

      ![Since we have the position of the child boundary in terms of i/j of the
      ! parent sub-domain, we can simply interpolate using the i/j fraction.
      !
      ! For this we use bilinear interpolation: (en.wikipedia.org/wiki/Bilinear_interpolation)
      !
      !   f(x,y) = 1/[(x_2-x_1)(y_2-y_1)] | f(Q_11)(x_2-x  )(y_2-y  ) + |      (1)
      !                                   | f(Q_21)(x  -x_1)(y_2-y  ) + |
      !                                   | f(Q_21)(x_2-x  )(y  -y_1) + |
      !                                   | f(Q_21)(x  -x_1)(y  -y_1)   |
      !
      ! In this routine we pre-compute all the interpolation coeffs
      !
      !   cff(1) = (x_2-x  )(y_2-y  ) / [(x_2-x_1)(y_2-y_1)]                   (2)
      !
      ! Since we are interpolating from 4 parent nodes that are exactly 1 point apart in i/j,
      ! the demoninator in (2) is just 1 and can be ignored, giving:
      !
      !   cff(1) = (x_2-x  )(y_2-y  )                                          (3)
      !
      ! and so on, such that:
      !
      !   f(x,y) = f(Q_11)*cff(1) + f(Q_21)*cff(2) + f(Q_12)*cff(3) + f(Q_22)*cff(4)
      !
      ! These interpolation coefficients can be used for all vertical layers.
      !
      !]

      implicit none

      ! hidden_mpi_vars.h needed for inode & jnode to establish if boundary edge
# include "hidden_mpi_vars.h"

      ! input/output
      integer,                              intent(in)  :: bry_tot     ! total child bry points in parent
      real,    dimension(bry_tot),          intent(in)  :: i_bry
      real,    dimension(bry_tot),          intent(in)  :: j_bry
      real,    dimension(:,:), allocatable, intent(out) :: rho_cffs    ! 4 parent interpolants for each child point
      integer, dimension(:),   allocatable, intent(out) :: i_sw, j_sw  ! i and j of sw parent interpolant

      ! local
      integer :: ij
      real    :: ri, rj                       ! position fractions

      allocate( rho_cffs(4,bry_tot) )         ! bry_tot for each bry points in parent. 4 for all four interpolants
      allocate( i_sw(bry_tot) )               ! i_sw is i grid point of sw parent interpolant
      allocate( j_sw(bry_tot) )

!      if(mynode==4) then
!        print *,'calcB FI-glob_istr=',floor(i_bry(1)),glob_istr,mynode
!        print *,'calcB FJ-glob_jstr=',floor(j_bry(1)),glob_jstr,mynode
!      endif

      do ij=1,bry_tot                                  ! Since glob_istr=0 on western_edge to catch child in ghost row in assign_child_interpolation_coeffs
                                                       ! we need to treat W & S edges diferently.
        if (inode==0) then                             ! western_edge
          i_sw(ij) = floor(i_bry(ij)) - glob_istr      ! local coords, set the SW parent i position. western edge glob_istr=0 not 1 hence no +1
        else
          i_sw(ij) = floor(i_bry(ij)) - glob_istr + 1  ! +1 since e.g. glob_jstr=51, j_bry=51.23 want sw corner to be 1 not 0
        endif

        if (jnode==0) then                             ! eastern_edge
          j_sw(ij) = floor(j_bry(ij)) - glob_jstr      ! local coords, set the SW parent j position
        else
          j_sw(ij) = floor(j_bry(ij)) - glob_jstr + 1
        endif

        ri = i_bry(ij) - floor(i_bry(ij))             ! i fraction (e.g. 23.47 - 23 = 0.47)
        rj = j_bry(ij) - floor(j_bry(ij))             ! j fraction

        rho_cffs(1,ij) = (1-ri)*(1-rj)        ! sw parent interp
        rho_cffs(2,ij) =    ri *(1-rj)        ! se parent interp
        rho_cffs(3,ij) = (1-ri)*   rj         ! nw parent interp
        rho_cffs(4,ij) =    ri *   rj         ! ne parent interp

!        if(mynode==5) then ! debug
!          if(ij<3) then
!            print *, 'ij=',ij,'i_bry(ij)=',i_bry(ij),'j_bry(ij)=',j_bry(ij)
!            print *, 'ri=',ri,'rj=',rj
!            print *, 'cff_sum=', rho_cffs(1,ij)+rho_cffs(2,ij)+rho_cffs(3,ij)+rho_cffs(4,ij)
!          endif
!        endif

      enddo

      end subroutine calc_bry_coeffs  !]

! ----------------------------------------------------------------------
      subroutine interpolate_child_bry( var, nmax )  ![

      ! Interpolate all boundaries if applicable

      implicit none

      integer, intent(in):: nmax
      real               :: var(GLOBAL_2D_ARRAY,nmax)  ! variable to interpolate.
                                                            ! wasn't able to use dimension(:,:,:), but vars start from -1 index so maybe this is safer.

      if (allocated(i_south))
     &  call interpolate_one_child_bry(rho_tot(1,1), south_r_cffs, is_sw, js_sw,
     &                                 nmax, var, vsouth)

      if (allocated(i_west))
     &  call interpolate_one_child_bry(rho_tot(1,2),  west_r_cffs, iw_sw, jw_sw,
     &                                 nmax, var, vwest)

!      if(mynode==1) then  ! debug
!        print *, 'vwest: rho_tot=', rho_tot(1,2),'mynode',mynode
!      print *,'i/j(end)=',i_west(rho_tot(1,2),1),j_west(rho_tot(1,2),1)
!        print *,'i/jsw(end)=',iw_sw(rho_tot(1,2)),jw_sw(rho_tot(1,2))
!        print *,'vwest(end)=',vwest(rho_tot(1,2),1)
!      endif
!      if(mynode==4) then  ! debug
!        print *, 'vwest: rho_tot=', rho_tot(1,2),'mynode',mynode
!        print *,'i/j(1)=',i_west(1,1),j_west(1,1)
!        print *,'i/jsw(1)=',iw_sw(1),jw_sw(1)
!        print *,'vwest(1)=',vwest(1,1)
!      endif

      if (allocated(i_north))
     &  call interpolate_one_child_bry(rho_tot(2,3), north_r_cffs, in_sw, jn_sw,  ! Note (2) as using northern most row of data
     &                                 nmax, var, vnorth)

      if (allocated(i_east))
     &  call interpolate_one_child_bry(rho_tot(2,4),  east_r_cffs, ie_sw, je_sw,  ! Note (2) as using eastern  most row of data
     &                                 nmax, var, veast)


      end subroutine interpolate_child_bry  !]

! ----------------------------------------------------------------------
      subroutine interpolate_one_child_bry(bry_tot, rho_cffs, i_sw, j_sw, nmax, var, vbry)  ![

      ! Interpolation coefficients already calculated in calc_child_interpolation_coeffs

      implicit none

      ! inputs
      integer,                              intent(in)  :: bry_tot     ! total child bry points in parent
      real,    dimension(:,:),              intent(in)  :: rho_cffs    ! 4 parent interpolants for each child point
      integer, dimension(:),                intent(in)  :: i_sw, j_sw  ! i and j of sw parent interpolant
      integer,                              intent(in)  :: nmax        ! vertical layers of var (1 if 2D)
      real,                                 intent(in)  :: var(GLOBAL_2D_ARRAY,nmax)  ! parent var to interpolate from
      real,    dimension(:,:), allocatable, intent(out) :: vbry        ! bry to interpolated 2D/3D var


      ! local
      integer :: ij, k
      real, dimension(nmax) :: var_sw, var_se, var_nw, var_ne  ! temporary copies of columns of data
      real(kind=8), parameter :: spv_set=1.D+33

      if (.not. allocated(vbry)) allocate( vbry(bry_tot, N))   ! allocate once for 3D var and reuse for all 2D/3D vars
                                                               ! maybe move this to initialization to avoid if statement repetition


      ! Normally would have k as outer loop for cache coherence but can avoid many
      ! if statements by swapping order here:
      do ij=1,bry_tot

        if(i_sw(ij) >= 0) then  ! child point not masked. >=0 as mask value is -1. i_sw can be 0 if point is between nodes.

          if ( rmask( i_sw(ij)  , j_sw(ij)  ) == 0) then
            call interp_masked_var( i_sw(ij)  , j_sw(ij)  , nmax, var, var_sw )
          else
            var_sw = var(i_sw(ij)  , j_sw(ij)  , :)
          endif

          if ( rmask( i_sw(ij)+1, j_sw(ij)  ) == 0) then
            call interp_masked_var( i_sw(ij)+1, j_sw(ij)  , nmax, var, var_se )
          else
            var_se = var(i_sw(ij)+1, j_sw(ij)  , :)
          endif

          if ( rmask( i_sw(ij)  , j_sw(ij)+1) == 0) then
            call interp_masked_var( i_sw(ij)  , j_sw(ij)+1, nmax, var, var_nw )
          else
            var_nw = var(i_sw(ij)  , j_sw(ij)+1, :)
          endif

          if ( rmask( i_sw(ij)+1, j_sw(ij)+1) == 0) then
            call interp_masked_var( i_sw(ij)+1, j_sw(ij)+1, nmax, var, var_ne )
          else
            var_ne = var(i_sw(ij)+1, j_sw(ij)+1, :)
          endif


          do k=1,nmax

            vbry(ij,k) = rho_cffs(1,ij) * var_sw(k)   ! interp from sw parent point.      Note even though 'k' index doesn't exist for 2D var,
     &                 + rho_cffs(2,ij) * var_se(k)   ! interp from se parent point.      fortran won't complain. This is done in ncdf_write_mod() too.
     &                 + rho_cffs(3,ij) * var_nw(k)   ! interp from nw parent point.
     &                 + rho_cffs(4,ij) * var_ne(k)   ! interp from ne parent point.

          enddo

        else  ! child point masked

          vbry(ij,:) = spv_set

        endif ! <- i_sw(ij)>0

      enddo   ! <- ij

      end subroutine interpolate_one_child_bry  !]

! ----------------------------------------------------------------------
      subroutine create_child_output_file  ![

      implicit none

      ! local
      integer               :: prev_fill_mode
      integer, dimension(5) :: rho_dimid_bry   ! dimension ID for boundaries
      integer               :: itrc            ! tracer index

      if(total_output_rec == 0) then                            ! Add root output name to
        file_name = trim(output_root_name) / / trim(file_name)  ! file name (only done once)
      endif

      ! For now decided each file still creates file even if not child contained in node.
      ! so that ncjoin can still check files are correct.
      ! Might need to then include all dimensions even if no child bry, for ncjoin to work
      call ncdf_create_blank_file(file_name, ncid_out,  prev_fill_mode,
     &                              recs_per_file, total_output_rec       )

      ! DEFINE DIMENSIONS & VARIABLES:

      ierr=nf90_put_att(ncid_out, NF90_GLOBAL, 'global_south_size', glob_c_s)
      ierr=nf90_put_att(ncid_out, NF90_GLOBAL, 'global_west_size' , glob_c_w)
      ierr=nf90_put_att(ncid_out, NF90_GLOBAL, 'global_north_size', glob_c_n)
      ierr=nf90_put_att(ncid_out, NF90_GLOBAL, 'global_east_size' , glob_c_e)
      ierr=nf90_def_dim( ncid_out, 'rec_time', 0, dimid_t)

      ! rho dims:
      ierr=nf90_def_dim( ncid_out, 'xi_rho_south', max( rho_tot(1,1), 1 ), rho_dimid_bry(1) )  ! Set size as 1 if size==0 to avoid unlimited dimension (could be confused with time in ncjoin)
      ierr=nf90_def_dim( ncid_out, 'eta_rho_west', max( rho_tot(1,2), 1 ), rho_dimid_bry(2) )
      ierr=nf90_def_dim( ncid_out, 'xi_rho_north', max( rho_tot(2,3), 1 ), rho_dimid_bry(3) )
      ierr=nf90_def_dim( ncid_out, 'eta_rho_east', max( rho_tot(2,4), 1 ), rho_dimid_bry(4) )

      ierr=nf90_def_dim( ncid_out, 's_rho', N, rho_dimid_bry(5) )

      ! u dims:

      ! v dims:

      ! define boundary even if not contained in parent for joining tool.
      call define_child_bry( indxZ, rho_tot, rho_c_str, rho_dimid_bry, 1, .false. )

      ! tracers:
      do itrc=1,NT
        call define_child_bry( itrc , rho_tot, rho_c_str, rho_dimid_bry, N, .true.  )
      enddo


      ierr=nf90_enddef(ncid_out)  ! end definition to write changes to disk
      if (mynode == 0) then
        write(*,'(/1x,2A/)')
     &     'bry_extract_child :: created new netCDF file ', file_name
      endif

      end subroutine create_child_output_file  !]

! ----------------------------------------------------------------------
      subroutine wrt_child_bry ( indx_v, nmax, trc_yn )  ![

      ! Each sub-domain writes out its portion of child boundary if applicable.
      ! all procs, make child bc netcdf file, only relevant procs define bc variables (east, north, etc)

      implicit none

      ! inputs
      integer, intent(in)   :: indx_v  ! index of primary var or itrc
      integer, intent(in)   :: nmax    ! vertical levels
      logical, intent(in)   :: trc_yn  ! is variable a tracer y/n

      ! local
      integer               :: varid, ndims
      integer               :: dimid, dim_len ! debug only
      integer, dimension(3) :: start_out, count_out, count_dummy=1
      character(len=40)     :: bname, rname


      if (trc_yn) then                              ! handle primary variable or tracers name
        rname=t_vname(indx_v)
      else
        rname=vname(1,indx_v)
      endif

      if (nmax == 1) then                           ! handle 2D or 3D variable dimids
        ndims=2
        start_out(1:ndims) = (/1, output_rec/)
        count_out(1:ndims) = (/1, 1         /)
      else
        ndims=3
        start_out(1:ndims) = (/1, 1   , output_rec/)
        count_out(1:ndims) = (/1, nmax, 1         /)
      endif


      bname = trim(rname) / /'_south'
      ierr  = nf90_inq_varid(ncid_out, bname, varid )
      count_out(1) = rho_tot(1,1)
      if (allocated(vsouth)) then
        ierr=nf90_put_var( ncid_out, varid , vsouth(:,:),         ! put dummy value so ncjoin doesn't complain
     &                     start=start_out(1:ndims), count=count_out(1:ndims) )
        if(ierr/=0) write(*,*) 'wrts_bry ierr= ',nf90_strerror(ierr), mynode
      else

        ierr=nf90_put_var( ncid_out, varid, (/1.0/),
     &                     start=start_out(1:ndims), count=count_dummy(1:ndims) )
      endif

      bname = trim(rname) / /'_west'
      ierr=nf90_inq_varid(ncid_out, bname, varid )
      count_out(1) = rho_tot(1,2)
      if (allocated(i_west )) then
        ierr=nf90_put_var( ncid_out, varid    , vwest(:,:) ,         ! put dummy value
     &                     start=start_out(1:ndims), count=count_out(1:ndims) )
        if(ierr/=0) write(*,*) 'wrtw_bry ierr= ',nf90_strerror(ierr), mynode
      else
        ierr=nf90_put_var( ncid_out, varid, (/1.0/),
     &                     start=start_out(1:ndims), count=count_dummy(1:ndims) )
      endif

      bname = trim(rname) / /'_north'
      ierr=nf90_inq_varid(ncid_out, bname, varid )
      count_out(1) = rho_tot(2,3)
      if (allocated(i_north)) then
        ierr=nf90_put_var( ncid_out, varid    , vnorth(:,:),         ! put dummy value
     &                     start=start_out(1:ndims), count=count_out(1:ndims) )
        if(ierr/=0) write(*,*) 'wrtn_bry ierr= ',nf90_strerror(ierr), mynode
      else
        ierr=nf90_put_var( ncid_out, varid, (/1.0/),
     &                     start=start_out(1:ndims), count=count_dummy(1:ndims) )
      endif

      bname = trim(rname) / /'_east'
      ierr=nf90_inq_varid(ncid_out, bname, varid )
      count_out(1) = rho_tot(2,4)
      if (allocated(i_east )) then
        ierr=nf90_put_var( ncid_out, varid    , veast(:,:) ,         ! put dummy value
     &                     start=start_out(1:ndims), count=count_out(1:ndims) )
        if(ierr/=0) write(*,*) 'wrte_bry ierr= ',nf90_strerror(ierr), mynode
      else
        ierr=nf90_put_var( ncid_out, varid, (/1.0/),
     &                     start=start_out(1:ndims), count=count_dummy(1:ndims) )
      endif

      end subroutine wrt_child_bry  !]

! ----------------------------------------------------------------------
      subroutine define_child_bry( indx_v, sides_tot, sides_c_str,
     &                             dimid_bry, nmax, trc_yn )        ![ Define child bry variable in output netCDF file.
      implicit none

      ! input/outputs
      integer,                 intent(in) :: indx_v
      integer, dimension(2,4), intent(in) :: sides_tot
      integer, dimension(2,4), intent(in) :: sides_c_str
      integer, dimension(5),   intent(in) :: dimid_bry   ! dimension ID for boundaries
      integer,                 intent(in) :: nmax
      logical,                 intent(in) :: trc_yn

      ! local
      integer               :: varid, bry, ndims
      character(len=40)     :: bname, rname
      character(len=60)     :: lname, uname
      integer, dimension(4) :: ghost_p=(/ 1, 1, 2, 2 /)  ! Ghost point position of 2 rows of input i/j
                                                         ! for S & W row 1, for N & E row 2
      integer, dimension(3) :: vdimids


      if (trc_yn) then           ! handle primary variable or tracers name
        rname=t_vname(indx_v)
        lname=t_lname(indx_v)
        uname=t_units(indx_v)
      else
        rname=vname(1,indx_v)
        lname=vname(2,indx_v)
        uname=vname(3,indx_v)
      endif

      do bry=1,4                                                 ! loop through 4 boundaries:

        if (nmax == 1) then                                      ! handle 2D or 3D variable dimids
          ndims=2
          vdimids(1:ndims) = (/dimid_bry(bry), dimid_t/)
        else
          ndims=3
          vdimids(1:ndims) = (/dimid_bry(bry), dimid_bry(5), dimid_t/)
        endif


        bname=trim(rname) / /'_'/ /trim( bry_name(bry) )

        call nc_define_var(ncid_out, bname, lname, uname,        ! longname, units
     &                     vdimids(1:ndims), ierr )


        ierr=nf90_inq_varid(ncid_out, bname, varid)

        ierr=nf90_put_att(ncid_out, varid, 'start', sides_c_str( ghost_p(bry), bry ))  ! add attribute with global child start value (i.e. south_c_str)
        ierr=nf90_put_att(ncid_out, varid, 'count', sides_tot(   ghost_p(bry), bry ))  ! and its count (i.e. south_tot). This will make it easy for ncjoin.
        ierr=nf90_put_att(ncid_out, varid, 'child_bry', 'True')                        ! 'child_bry' for ncjoin to identify variable. 'True' just a filler.

      enddo

      end subroutine define_child_bry  !]

! ----------------------------------------------------------------------
      subroutine interp_masked_var( i, j, nmax, var, vinterp )  ![

      ![ Since we interpolate child point from 4 parent points, there are
      !  situations where 1, 2 or 3 of the parent points are in the land mask.
      !
      !  This routine interpolates those specific points using the 8 surrounding
      !  surrounding grip points, using a inverse distance weighting method:
      !
      !    f(x) = sum( w_i(x) f_i ) / sum( w_i(x) )
      !
      !  where:
      !
      !    w_i(x) = 1 / d(x,x_i)
      !
      !  If any of the 8 surrounding points are also in the mask, we set their
      !  weighting contribution to zero. This way the total contribution of the
      !  weightings is still 1. This is a crude approximation, but we are favoring
      !  simplicity here.
      !
      !  Since we are interpolating from +/- 1 in i/j, the weights are 1 or 1/sqrt(2)
      !  as seen in the following diagram.
      !
      !  Diagram: X point interpolated from 8 surrouring points A-H
      !
      !       F   G   H               weights:
      !                               X to middles B,D,E,G = 1
      !       D   X   E               X to corners A,C,F,H = 1/sqrt(2)
      !
      !       A   B   C
      !]

      implicit none

      ! input/output
      integer,               intent(in)  :: i, j                       ! indices of masked point to interpolate
      integer,               intent(in)  :: nmax                       ! vertical layers of variable
      real,                  intent(in)  :: var(GLOBAL_2D_ARRAY,nmax)  ! parent variable array
      real, dimension(nmax), intent(out) :: vinterp                    ! interpolated scalar (2D var) or column (3D var) at masked point

      ! local
      real    :: ca, cc, cf, ch, mb, md, me, mg  ! middle & corner interpolants as per diagram above
      real    :: sum_cffs                        ! sum of weightings (denominator in interpolation)
      integer :: k

      ca = rt2*vmask(i-1,j-1)  ! corner point weightings
      cc = rt2*vmask(i+1,j-1)  ! (calculate once outside of vertical loop)
      cf = rt2*vmask(i-1,j+1)
      ch = rt2*vmask(i+1,j+1)

      mb = vmask(i  ,j-1)      ! middle point weightings
      md = vmask(i-1,j  )
      me = vmask(i+1,j  )
      mg = vmask(i  ,j+1)

      sum_cffs = ca+cc+cf+ch+mb+md+me+mg

      do k=1,nmax

        vinterp(k) = (  var(i-1,j-1,k)*ca
     &                + var(i+1,j-1,k)*cc
     &                + var(i-1,j+1,k)*cf
     &                + var(i+1,j+1,k)*ch
     &                + var(i  ,j-1,k)*mb
     &                + var(i-1,j  ,k)*md
     &                + var(i+1,j  ,k)*me
     &                + var(i  ,j+1,k)*mg )
     &                / sum_cffs

      enddo ! <- k

      end subroutine interp_masked_var  !]

! ----------------------------------------------------------------------

#endif /* BRY_EXTRACT_CHILD */
      end module bry_extract_child
