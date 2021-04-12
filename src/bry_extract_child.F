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
      !    Each i'th/j'th data consists of 2 rows per boundary side.
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

      implicit none
      private

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

      real, dimension(:,:), allocatable :: i_west, i_south, i_east, i_north  ! i position of child bry in parent grid
      real, dimension(:,:), allocatable :: j_west, j_south, j_east, j_north  ! j position of child bry in parent grid

      integer, dimension(2)             :: west_tot=0,   south_tot=0,        ! Number of child points within parent
     &                                     east_tot=0,   north_tot=0
      integer, dimension(2)             :: west_c_str=0, south_c_str=0,      ! start indices of child boundary within parent
     &                                     east_c_str=0, north_c_str=0

      integer :: glob_c_s, glob_c_w, glob_c_n, glob_c_e                      ! Total number of global child points along child boundary

      ! INTERPOLATION VARS:
      real,    dimension(:,:), allocatable :: vwest, vsouth, veast, vnorth   ! array storing interpolated variable (2D or 3D)
      real,    dimension(:,:), allocatable :: west_r_cffs, south_r_cffs,     ! 4 interpolation coefficients to interp from 4 parent points surronding child point.
     &                                        east_r_cffs, north_r_cffs      ! i=each bry point. j=1-4 (4 interpolants). Will need for u_cff, v_cff. 4 interpolation coefficients
      integer, dimension(:), allocatable   :: iw_sw, jw_sw, is_sw, js_sw,    ! sw parent interpolant i and j position
     &                                        ie_sw, je_sw, in_sw, jn_sw

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

      ! Use this to only set at required output rate. Rather have it here than in main
      ! to keep main less cluttered? Although wrt_his has it in main? Ask JM
      ! only plot if timestep is divisble by record rate
!      if (mod(iic-1,rec_rate) == 0) then ! 1 = ntstart

      if (ncid==-1) then                       ! First time only
        call read_child_bry
        call calc_child_interpolation_coeffs
      endif

      ! DO INTERPOLATION:
      call interpolate_child_bry

      ! WRT CHILD BOUNDARY VALUES:
      call wrt_child_bry

!      call MPI_Barrier(MPI_COMM_WORLD, ierr)  ! debug
!      call MPI_Abort(MPI_COMM_WORLD, 1, ierr) ! debug

!      endif ! <- mod(iic-1,rec_rate)

      end subroutine set_child_bry  !]

! ----------------------------------------------------------------------
      subroutine read_child_bry  ![
      ! Assign 4 parent GP from which to interpolate 1 child BP. Done only once.

      ! a) Loop through contained child BP's from part 1) and find 4 enclosing parent GP's.
      ! b) Calulate (once) interpolation coefficients for contribution of 4 parent GP's to 1 child BP.

      implicit none

      ! Do for 4 sides, if it has i then it must have j:
      call assign_child_bry_to_parent(bry_name(1), i_south, j_south, south_tot,
     &                                south_c_str, glob_c_s )
      call MPI_Barrier(MPI_COMM_WORLD, ierr)  ! debug

      call assign_child_bry_to_parent(bry_name(2), i_west,  j_west,  west_tot ,
     &                                west_c_str , glob_c_w )
      call MPI_Barrier(MPI_COMM_WORLD, ierr)  ! debug

      call assign_child_bry_to_parent(bry_name(3), i_north, j_north, north_tot,
     &                                north_c_str, glob_c_n )
      call MPI_Barrier(MPI_COMM_WORLD, ierr)  ! debug

      call assign_child_bry_to_parent(bry_name(4), i_east,  j_east,  east_tot ,
     &                                east_c_str , glob_c_e )
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
      integer,                           intent(out)   :: dim_size      ! Size of global child boundary

      ! local
      integer               :: varid, ij, l, counter
      character(len=30)     :: vname
      integer, dimension(2) :: dimids, c_str  ! full child boundary start index falling in parent

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
        if(mynode==0) write (*,*) 'No child boundary for: ', vname
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
        call calc_bry_coeffs(south_tot(1), i_south(:,1), j_south(:,1),
     &                       south_r_cffs, is_sw, js_sw)                 ! index 1 since uses southern ghost points.
      endif
      if (allocated(i_west))  then
        call calc_bry_coeffs(west_tot(1),  i_west(:,1),  j_west(:,1),
     &                       west_r_cffs,  iw_sw, jw_sw)                 ! index 1 since uses wester   ghost points.
      endif
      if (allocated(i_north)) then
        call calc_bry_coeffs(north_tot(2), i_north(:,2), j_north(:,2),
     &                       north_r_cffs, in_sw, jn_sw)                 ! index 2 since uses northern ghost points.
      endif
      if (allocated(i_east))  then
        call calc_bry_coeffs(east_tot(2),  i_east(:,2),  j_east(:,2),
     &                       east_r_cffs,  ie_sw, je_sw)                 ! index 2 since uses eastern  ghost points.
      endif

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

      do ij=1,bry_tot

        i_sw(ij) = floor(i_bry(ij)) - glob_istr          ! set the SW parent i position
        j_sw(ij) = floor(j_bry(ij)) - glob_jstr          ! set the SW parent j position

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
      subroutine interpolate_child_bry  ![

      ! interpolate all boundaries if applicable

      implicit none

#include "ocean2d.h"
! scalars.h for 'knew' var
#include "scalars.h"

      if (allocated(i_south))
     &  call interpolate_one_child_bry(south_tot(1), south_r_cffs, is_sw, js_sw,
     &                                 1, zeta(:,:,knew), vsouth)

      if (allocated(i_west))
     &  call interpolate_one_child_bry(west_tot(1),  west_r_cffs,  iw_sw, jw_sw,
     &                                 1, zeta(:,:,knew), vwest)

      if (allocated(i_north))
     &  call interpolate_one_child_bry(north_tot(1), north_r_cffs, in_sw, jn_sw,
     &                                 1, zeta(:,:,knew), vnorth)

      if (allocated(i_east))
     &  call interpolate_one_child_bry(east_tot(1),  east_r_cffs,  ie_sw, je_sw,
     &                                 1, zeta(:,:,knew), veast)


      end subroutine interpolate_child_bry  !]
! ----------------------------------------------------------------------
      subroutine interpolate_one_child_bry(bry_tot, rho_cffs, i_sw, j_sw, nz, var, vbry)  ![

      ! Interpolation coefficients already calculated in calc_child_interpolation_coeffs

      implicit none

      ! inputs
      integer,                              intent(in)  :: bry_tot     ! total child bry points in parent
      real,    dimension(:,:),              intent(in)  :: rho_cffs    ! 4 parent interpolants for each child point
      integer, dimension(:),                intent(in)  :: i_sw, j_sw  ! i and j of sw parent interpolant
      integer,                              intent(in)  :: nz          ! vertical layers of var (1 if 2D)
      real,                                 intent(in)  :: var(GLOBAL_2D_ARRAY,nz)  ! parent var to interpolate from
      real,    dimension(:,:), allocatable, intent(out) :: vbry        ! bry to interpolated 2D/3D var


      ! local
      integer :: ij, k

      if (.not. allocated(vbry)) allocate( vbry(bry_tot, N))  ! allocate once for 3D var and reuse for all 2D/3D vars
                                                              ! maybe move this to initialization to avoid if statement repetition
      do ij=1,bry_tot
        do k=1,nz
          vbry(ij,k) = rho_cffs(1,ij) * var(i_sw(ij)  , j_sw(ij)  , k)   ! interp from sw parent point.      Note even though 'k' index doesn't exist for 2D var,
     &               + rho_cffs(2,ij) * var(i_sw(ij)+1, j_sw(ij)  , k)   ! interp from se parent point.      fortran won't complain. This is done in ncdf_write_mod() too.
     &               + rho_cffs(3,ij) * var(i_sw(ij)  , j_sw(ij)+1, k)   ! interp from se parent point.
     &               + rho_cffs(4,ij) * var(i_sw(ij)+1, j_sw(ij)+1, k)   ! interp from se parent point.
        enddo
      enddo


      end subroutine interpolate_one_child_bry  !]

! ----------------------------------------------------------------------
      subroutine wrt_child_bry  ![

      ! Each sub-domain writes out its portion of child boundary if applicable.
      ! all procs, make child bc netcdf file, only relevant procs define bc variables (east, north, etc)

      ! For simplicity have only one variable based on the longest of the 2 rows.
      ! Can add row length as attribute for ncjoin to know how far to go.
      ! Can still write based on length

      implicit none

      ! local
      integer               :: prev_fill_mode, varid
      integer               :: dimid, dim_len ! debug only
      integer, dimension(2) :: start

      if(total_output_rec == 0) then                            ! Add root output name to
        file_name = trim(output_root_name) / / trim(file_name)  ! file name (only done once)
      endif

      if(output_rec == 0) then  ! file doesn't yet exist or need new file

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

!        if (allocated(i_south)) then  ! contains south child boundary ! define boundary even if not contained in parent for joining tool.
        call define_child_bry( 'zeta_south', 'xi_rho_south', south_tot, south_c_str, indxZ )
!        endif
        call define_child_bry( 'zeta_west' , 'eta_rho_west', west_tot , west_c_str,  indxZ )
        call define_child_bry( 'zeta_north', 'xi_rho_north', north_tot, north_c_str, indxZ )
        call define_child_bry( 'zeta_east' , 'eta_rho_east', east_tot , east_c_str,  indxZ )

        ierr=nf90_enddef(ncid_out)  ! end definition to write changes to disk
        if (mynode == 0) then
          write(*,'(/1x,2A/)')
     &     'bry_extract_child :: created new netCDF file ', file_name
        endif

      endif ! <- (output_rec == 0; needed to create a file)


      ! WRITE VARIABLES TO OUTPUT:

      ierr=nf90_open(file_name, nf90_write, ncid_out)

      output_rec = output_rec + 1              ! Current file record
      total_output_rec = total_output_rec + 1  ! total number of recs written
      start=(/1,1/) ! (/1,output_rec/)

      ierr=nf90_inq_varid(ncid_out, 'zeta_south', varid )
      if (allocated(vsouth)) then
        ierr=nf90_put_var(  ncid_out, varid    , vsouth(:,:),
     &                      start=start, count=(/south_tot(1),1/) )
        if(ierr/=0) write(*,*) 'wrts_bry ierr= ',nf90_strerror(ierr), mynode
      else
        ! Put in dummy value so that ncjoin doesn't complain
        ierr=nf90_put_var(  ncid_out, varid, (/1.0/), start=start, count=(/1,1/) )
      endif

      ierr=nf90_inq_varid(ncid_out, 'zeta_west' , varid )
      if (allocated(i_west )) then
        ierr=nf90_put_var(  ncid_out, varid    , vwest(:,:) ,
     &                      start=start, count=(/west_tot(1),1/) )
        if(ierr/=0) write(*,*) 'wrtw_bry ierr= ',nf90_strerror(ierr), mynode
      else
        ! Put in dummy value so that ncjoin doesn't complain
        ierr=nf90_put_var(  ncid_out, varid, (/1.0/), start=start, count=(/1,1/) )
      endif

      ierr=nf90_inq_varid(ncid_out, 'zeta_north', varid )
      if (allocated(i_north)) then
        ierr=nf90_put_var(  ncid_out, varid    , vnorth(:,:),
     &                      start=start, count=(/north_tot(1),1/) )
        if(ierr/=0) write(*,*) 'wrtn_bry ierr= ',nf90_strerror(ierr), mynode
      else
        ! Put in dummy value so that ncjoin doesn't complain
        ierr=nf90_put_var(  ncid_out, varid, (/1.0/), start=start, count=(/1,1/) )
      endif

      ierr=nf90_inq_varid(ncid_out, 'zeta_east' , varid )
      if (allocated(i_east )) then
        ierr=nf90_put_var(  ncid_out, varid    , veast(:,:) ,
     &                      start=start, count=(/east_tot(1),1/) )
        if(ierr/=0) write(*,*) 'wrte_bry ierr= ',nf90_strerror(ierr), mynode
      else
        ! Put in dummy value so that ncjoin doesn't complain
        ierr=nf90_put_var(  ncid_out, varid, (/1.0/), start=start, count=(/1,1/) )
      endif


      ierr=nf90_close(ncid_out)

      if (output_rec == recs_per_file) then  ! Reset output_rec to create new file
        output_rec = 0
      endif

      end subroutine wrt_child_bry  !]

! ----------------------------------------------------------------------
      subroutine define_child_bry( bry_name, dim_name, side_tot, side_c_str, indx_v ) ![ Define child bry variable in output netCDF file.
      implicit none

      ! input/outputs
      character(len=*),      intent(in) :: bry_name
      character(len=*),      intent(in) :: dim_name
      integer, dimension(2), intent(in) :: side_tot
      integer, dimension(2), intent(in) :: side_c_str
      integer,               intent(in) :: indx_v

      ! local
      integer, dimension(4)             :: dimid_bry   ! dimension ID for boundaries
      integer                           :: varid
      character(len=40)                 :: lname

      ! Possibly need to check if already exists rather?
      ierr=nf90_def_dim( ncid_out, dim_name, max( side_tot(1), 1 ), dimid_bry(1) )  ! Set size as 1 if size==0 to avoid unlimited dimension (could be confused with time in ncjoin)

      lname=trim(vname(2,indx_v))                                 ! / / ' of child bry from parent simulation'
      call nc_define_var(ncid_out, bry_name, lname,
     &      vname(3,indx_v),                                      ! units
     &       (/dimid_bry(1), dimid_t/), ierr )


      ierr=nf90_inq_varid(ncid_out, bry_name, varid)              ! ADD ATTRIBUTE WITH GLOBAL CHILD START VALUE (i.e. south_c_str)
      ierr=nf90_put_att(ncid_out, varid, 'start', side_c_str(1))  ! and its count (i.e. south_tot). This will make it easy for ncjoin.
      ierr=nf90_put_att(ncid_out, varid, 'count', side_tot(1))
      ierr=nf90_put_att(ncid_out, varid, 'child_bry', 'True')     ! True is not required

      end subroutine define_child_bry  !]

! ----------------------------------------------------------------------

#endif /* BRY_EXTRACT_CHILD */
      end module bry_extract_child
