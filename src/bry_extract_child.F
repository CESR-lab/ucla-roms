      module bry_extract_child   ! Interpolate and extract online the boundary files
                                 ! needed to run child (nested) simulation

      ![ INFO:

      ! coded by: Devin Dollery & Jeroen Molemaker - 2021 March

      ! STEPS:

      ! 1) READ CHILD BOUNDARY LOCATIONS:
      !
      !    The location of the child boundary is stored in one full-domain file
      !    and is read in by every mpi node.
      !    The child boundary file consists of fractional i'th and j'th positions of the
      !    child boundary points relative to the parent's index domain.
      !    Each i'th/j'th data consists of 2 rows per boundary side.
      !    Seen in DIAGRAM (A) below, each child point (a-f) will have a input variable
      !    giving it's i'th position, i.e. variable 'i_west', which for point b = 1.4.
      !    Another variable will hold the j'th position, i.e. 'j_west' which for point b = 1.2.
      !    If the boundary is closed then i and j variables will not exist for that boundary.

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


      ! 2) ASSIGN CHILD BRY TO PARENT NODE:
      !
      !    As mentioned each parent node will read in the full child boundary data.
      !    It then needs to calculate which (if any) of the child boundary fall within
      !    the parent node's range.
      !    Since the child bry is given in i/j relative to the parent grid points,
      !    one can simply compare the parent's node number against the child's i/j number.
      !    However, the child i/j values are relative to the global parent domain.
      !    Typically roms uses 1:Lm (istr:iend) for local xi ranges, and 1:Mm for local
      !    eta ranges of a parent node.  These are local to the parent node and
      !    thus you can't directly compare the child bry i/j to (istr:iend).
      !
      !    We need the global start and end i/j of the parent node, stored in variables:
      !    insert names here!!!
      !
      !    For further explanation of local and global grid point numbering see:
      !    Documentation/readme-specifics/README.grid.2021.03


      ! 3) INTERPOLATE CHILD BRY AT REQUIRED OUTPUT RATE:
      !
      !    When required the child bry is interpolated using the surrounding 4 parent
      !    grid points.
      !    The 4 parent grid points are calculated on the fly based on the floor and
      !    ceiling of the fractional i/j child numbers. I.e. child point b is interpolated
      !    from parent points A-D in DIAGRAM (A) using bilinear interpolation.


      ! 4) WRITE RESULTS TO FILE:
      !
      !


      ! 5) JOIN PARTITIONED RESULTS FILES:
      !
      !    This probably requires its own tool to join.
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

      !]

#include "cppdefs.h"
#if defined BRY_EXTRACT_CHILD

      use dimensions !! has dimensions, a list of forcing files, and rmask, glob_istr, etc
      use read_write
      use netcdf
      use mpi ! debug temporary for mpi abort

      implicit none
      private

      ! scalars.h needed for FIRST_TIME_STEP
#include "scalars.h"

      character(len=9)  :: child_bry_name=' '                     ! root name, append e/w/s/n
      character(len=50) :: c_bry_file='../sample_child_ij_bry.nc' ! child boundary input file

      real, dimension(:,:), allocatable :: i_west, i_south, i_east, i_north  ! i position of child bry in parent grid
      real, dimension(:,:), allocatable :: j_west, j_south, j_east, j_north  ! j position of child bry in parent grid

      integer, dimension(2) :: west_tot=-1, south_tot=-1, east_tot=-1, north_tot=-1  ! Number of child points within parent

      integer :: ncid=-1
      integer :: ierr=0 ! debug temporary

      public set_child_bry

      ! need: method to join results into 1 file after simulation.

      contains

! ----------------------------------------------------------------------
      subroutine set_child_bry
      ! Interpolate child boundary at desired intervals.

      implicit none

      if (ncid==-1) then     ! First time only
        call read_child_bry
      endif

      ! DO INTERPOLATION:


      ! WRT CHILD BOUNDARY VALUES:
      call wrt_child_bry

      call MPI_Barrier(MPI_COMM_WORLD, ierr)  ! debug
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr) ! debug

      end subroutine set_child_bry

! ----------------------------------------------------------------------
      subroutine read_child_bry
      ! Assign 4 parent GP from which to interpolate 1 child BP. Done only once.

      ! a) Loop through contained child BP's from part 1) and find 4 enclosing parent GP's.
      ! b) Calulate (once) interpolation coefficients for contribution of 4 parent GP's to 1 child BP.

      implicit none

      ! Do for 4 sides, if it has i then it must have j:
      call assign_child_bry_to_parent('south', i_south, j_south, south_tot )
      call MPI_Barrier(MPI_COMM_WORLD, ierr)  ! debug
      call assign_child_bry_to_parent('west',  i_west,  j_west,  west_tot  )
      call MPI_Barrier(MPI_COMM_WORLD, ierr)  ! debug
      call assign_child_bry_to_parent('north', i_north, j_north, north_tot )
      call MPI_Barrier(MPI_COMM_WORLD, ierr)  ! debug
      call assign_child_bry_to_parent('east',  i_east,  j_east,  east_tot  )
      call MPI_Barrier(MPI_COMM_WORLD, ierr)  ! debug

      end subroutine read_child_bry

! ----------------------------------------------------------------------
      subroutine assign_child_bry_to_parent(side_name, i_bry, j_bry,
     &                                      side_totals              )
      ! Read in N/W/E/S arrays of child BP's. (in i,j coordinates)

      ! all procs read the child bc file and determine if and what work they need to do.
      !If work, determine range along the boundary  Loop through all child BP's and
      ! compare to corners of parent node.
!      use mpi ! debug temporary for mpi abort
      implicit none

      ! input
      character(len=*), intent(in)                   :: side_name
      real, dimension(:,:), allocatable, intent(out) :: i_bry, j_bry  ! Final lists local to parent.
      integer, dimension(2), intent(inout)           :: side_totals   ! Number of child points within parent


      ! local
      integer               :: varid, dim_size, ij, l, counter
      character(len=30)     :: vname
      integer, dimension(2) :: dimids, c_str  ! child_start indices for those within parent

      real, dimension(:,:), allocatable :: i_in_bry, j_in_bry ! Store read in child boundary

      if (ncid==-1) then  ! File not yet opened.
        ierr=nf90_open(c_bry_file, nf90_nowrite, ncid) ! var will later be in grdname file
        if(ierr/=0) call handle_ierr(ierr,
     &             'read_child_bry:: Error opening file with child bry!')
      endif


      vname='i_' / / side_name
      if(mynode==0) print *, 'child bry name: ', vname
      ierr=nf90_inq_varid(ncid, vname, varid)          ! Check if contains bry
      if (ierr==nf90_noerr) then
        if(mynode==0) write (*,*) 'Found child boundary in file for: ', vname
        ierr=nf90_inquire_variable(ncid, varid, dimids = dimids)
        ierr=nf90_inquire_dimension(ncid, dimids(1), len=dim_size)
        if(mynode==0) write (*,*) 'dim_size=', dim_size

        allocate(i_in_bry(dim_size,2))                 ! Allocate temporary array for all bry data.
        ierr=nf90_get_var(ncid, varid, i_in_bry)       ! read in data
        if(ierr/=0) print *, 'bad read of bry'

        allocate(j_in_bry(dim_size,2))
        vname='j_' / / side_name
        ierr=nf90_inq_varid(ncid, vname, varid)
        ierr=nf90_get_var(ncid, varid, j_in_bry)

                            ! CHECK IF CHILD CONTAINED WITHIN PARENT NODE:
        do l=1,2            ! Need to check both rows as might slightly differ e.g. i=9.95 vs i=10.05
          counter=0
          do ij=1,dim_size

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

          write(*,'(/1x,I2,A)') mynode, ' Child points within parent node' ! debug
          print *, mynode, 'side_totals=', side_totals               ! debug
          print *, mynode, 'c_str(1)=',c_str(1),'c_str(2)=',c_str(2) ! debug
          print *, mynode, 'i_bry(str,1)= ', i_bry(1,1)              ! debug
          print *, mynode, 'i_bry(end,1)= ', i_bry(side_totals(1),1) ! debug
          print *, mynode, 'j_bry(str,1)= ', j_bry(1,1)              ! debug
          print *, mynode, 'j_bry(end,1)= ', j_bry(side_totals(1),1) ! debug
        endif

      else     ! Side not contained
        if(mynode==0) write (*,*) 'No child boundary for: ', vname
      endif

      end subroutine assign_child_bry_to_parent

! ----------------------------------------------------------------------
      subroutine wrt_child_bry
      ! Each sub-domain writes out its portion of child boundary if applicable.
      ! all procs, make child bc netcdf file, only relevant procs define bc variables (east, north, etc)

      implicit none



      ierr=nf90_create(fname, NF90_NETCDF4, ncid)

      end subroutine wrt_child_bry

! ----------------------------------------------------------------------

#endif /* BRY_EXTRACT_CHILD */
      end module bry_extract_child
