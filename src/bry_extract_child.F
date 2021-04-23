      module bry_extract_child   ! Interpolate and extract online the boundary files
                                 ! needed to run child (nested) simulation

      ! Auxiliary tools required:
      ! - Tools-Roms/bry_child_input_gen.m   -> creates input file of child bry i/j values used here
      ! - Tools-Roms/join_child_bry.F        -> joins partial files created here

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
      !    giving it's i'th position, i.e. variable 'ir_west', which for point b = 1.4.
      !    Another variable will hold the j'th position, i.e. 'jr_west' which for point b = 1.2.
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
      !    Tools-Roms/join_child_bry.F (program)
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
      integer           :: rec_rate      = 1                           ! output record rate
      integer           :: recs_per_file = 10                          ! total recs per file
      character(len=50) :: c_bry_file    = '../sample_child_ij_bry.nc' ! child boundary input file
      ! END USER INPUTS
      ! ****************************************************************

      ! DECLARE VARIABLES: ![

      ! Would be cleaner to have variables 1-4 for W/S/E/N rather than loads of variables.

      integer :: glob_c_rs, glob_c_rw, glob_c_rn, glob_c_re                         ! total rho number of global child points along child boundary
      integer :: glob_c_us, glob_c_uw, glob_c_un, glob_c_ue                         ! total u   ...
      integer :: glob_c_vs, glob_c_vw, glob_c_vn, glob_c_ve                         ! total v   ...
                                                                                    ! not necessary to store all of these but easier for subroutine arguements
      character(len=9),       dimension(4) :: bry_name = (/ 'south', 'west',
     &                                                      'north', 'east'  /)


      ! RHO-VARS:
      ! - input vars:
      character(len=9),       dimension(4) :: rbry_name = (/ 'rho_south', 'rho_west',
     &                                                       'rho_north', 'rho_east'  /)

      real,    dimension(:),   allocatable :: ir_south, ir_west, ir_north, ir_east  ! i position of child bry in parent grid
      real,    dimension(:),   allocatable :: jr_south, jr_west, jr_north, jr_east  ! j position of child bry in parent grid

      integer, dimension(4)                :: r_tot=0                               ! number of child points within parent on each boundary. 4 = s/w/n/e, 2 = 2 rows of i/j input
      integer, dimension(4)                :: r_c_str=0                             ! start indices of child boundary within parent
      ! - interp vars:
      real,    dimension(:,:), allocatable :: south_r_cffs, west_r_cffs,            ! 4 interpolation coefficients to interp from 4 parent points surronding child point.
     &                                        north_r_cffs, east_r_cffs             ! i=each bry point. j=1-4 (4 interpolants).

      integer, dimension(:),   allocatable :: isr_sw, jsr_sw, iwr_sw, jwr_sw,       ! sw parent interpolant i and j position
     &                                        inr_sw, jnr_sw, ier_sw, jer_sw        ! could merge with coefficients.

      ! U-VARS:
      ! - input vars:
      character(len=9),       dimension(4) :: ubry_name = (/ 'u_south', 'u_west',
     &                                                       'u_north', 'u_east'  /)

      real,    dimension(:),   allocatable :: iu_south, iu_west, iu_north, iu_east  ! i position of child bry in parent grid
      real,    dimension(:),   allocatable :: ju_south, ju_west, ju_north, ju_east  ! j position of child bry in parent grid

      integer, dimension(4)                :: u_tot=0                               ! number of child points within parent on each boundary. 4 = s/w/n/e, 2 = 2 rows of i/j input
      integer, dimension(4)                :: u_c_str=0                             ! start indices of child boundary within parent
      ! - interp vars:
      real,    dimension(:,:), allocatable :: south_u2u_cffs, west_u2u_cffs,        ! 4 interpolation coefficients to interp from 4 parent u-points surronding child u-point.
     &                                        north_u2u_cffs, east_u2u_cffs         ! i=each bry point. j=1-4 (4 interpolants).

      real,    dimension(:,:), allocatable :: south_v2u_cffs, west_v2u_cffs,        ! 4 interpolation coefficients to interp from 4 parent v-points surronding child u-point.
     &                                        north_v2u_cffs, east_v2u_cffs         ! i=each bry point. j=1-4 (4 interpolants).

      real,    dimension(:,:), allocatable :: south_r2u_cffs, west_r2u_cffs,        ! 4 interpolation coefficients to interp from 4 parent r-points surronding child u-point. (angle)
     &                                        north_r2u_cffs, east_r2u_cffs

      integer, dimension(:),   allocatable :: isu2u_sw, jsu2u_sw, iwu2u_sw, jwu2u_sw, ! uparent 2 uchild sw parent interpolant i and j position
     &                                        inu2u_sw, jnu2u_sw, ieu2u_sw, jeu2u_sw

      integer, dimension(:),   allocatable :: isv2u_sw, jsv2u_sw, iwv2u_sw, jwv2u_sw, ! vparent 2 uchild sw parent interpolant i and j position
     &                                        inv2u_sw, jnv2u_sw, iev2u_sw, jev2u_sw

      integer, dimension(:),   allocatable :: isr2u_sw, jsr2u_sw, iwr2u_sw, jwr2u_sw, ! vparent 2 uchild sw parent interpolant i and j position
     &                                        inr2u_sw, jnr2u_sw, ier2u_sw, jer2u_sw

      real,    dimension(:),   allocatable :: us_angcos, uw_angcos, un_angcos, ue_angcos, ! calculated rotation required at child bry point to get parent u/v to child u/v
     &                                        us_angsin, uw_angsin, un_angsin, ue_angsin  ! difference between parent grid angle at child point compared to child grid angle

      real,    dimension(:,:), allocatable :: usouth,   uwest,    unorth,   ueast,    ! final interpolated & rotated variable (2D or 3D)
     &                                        us_uprnt, uw_uprnt, un_uprnt, ue_uprnt, ! u value of parent at child u location
     &                                        us_vprnt, uw_vprnt, un_vprnt, ue_vprnt  ! v value of parent at child u location

      ! INTERPOLATION:
      real,    dimension(:,:), allocatable :: vsouth, vwest, vnorth, veast   ! array storing interpolated variable (2D or 3D)
      real                                 :: rt2                            ! square root of 2 stored for efficiency
                                                                             ! allocated for each var so only need one.
      ! NETCDF:
      integer            :: ncid=-1                     ! in file for child bry positions
      integer            :: ncid_out                    ! output file
      integer            :: ierr=0                      ! debug temporary
      integer            :: dimid_t                     ! time dimension id in output file
      integer            :: output_rec = 0              ! record number of output. 0 indicates we need new file
      integer            :: total_output_rec = 0        ! total flux output recordings so far
      character (len=60) :: file_name = '_child_bry.nc' ! suffix to append to simulation output rootname


      public set_child_bry

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
      call interpolate_child_bry( zeta(:,:,knew), 1 )  ! var, nz
      call wrt_child_bry( indxZ,  1, .false., vsouth, vwest, vnorth, veast, r_tot )

      ! ubar:
      call interpolate_u_child_bry( ubar(:,:,knew), vbar(:,:,knew), 1 )
      call wrt_child_bry( indxUb, 1, .false., usouth, uwest, unorth, ueast, u_tot )

      ! vbar:
      ! u:
      ! v:

      ! tracers:
      do itrc=1,NT
        call interpolate_child_bry( t(START_2D_ARRAY,1,nstp,itrc), N )
        call wrt_child_bry( itrc, N, .true., vsouth, vwest, vnorth, veast, r_tot )
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

      call assign_child_bry_to_parent(rbry_name(1), ir_south,  jr_south, r_tot(1),
     &                                r_c_str(1),   glob_c_rs, rp_var )
      call assign_child_bry_to_parent(rbry_name(2), ir_west,   jr_west,  r_tot(2),
     &                                r_c_str(2) ,  glob_c_rw, rp_var )
      call assign_child_bry_to_parent(rbry_name(3), ir_north,  jr_north, r_tot(3),
     &                                r_c_str(3),   glob_c_rn, rp_var )
      call assign_child_bry_to_parent(rbry_name(4), ir_east,   jr_east,  r_tot(4),
     &                                r_c_str(4) ,  glob_c_re, rp_var )

      call assign_child_bry_to_parent(ubry_name(1), iu_south,  ju_south, u_tot(1),
     &                                u_c_str(1),   glob_c_us, up_var )
      call assign_child_bry_to_parent(ubry_name(2), iu_west,   ju_west,  u_tot(2),
     &                                u_c_str(2) ,  glob_c_uw, up_var )
      call assign_child_bry_to_parent(ubry_name(3), iu_north,  ju_north, u_tot(3),
     &                                u_c_str(3),   glob_c_un, up_var )
      call assign_child_bry_to_parent(ubry_name(4), iu_east,   ju_east,  u_tot(4),
     &                                u_c_str(4) ,  glob_c_ue, up_var )

!      call MPI_Barrier(MPI_COMM_WORLD, ierr) ! debug
!      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)

      end subroutine read_child_bry  !]

! ----------------------------------------------------------------------
      subroutine assign_child_bry_to_parent( side_name,  i_bry, j_bry,            ![
     &                                       side_total, c_str, dim_size, vtype )
      ! Read in N/W/E/S arrays of child BP's. (in i,j coordinates)

      ! All procs read the child bc file and determine if and what work they need to do.
      ! If work, determine range along the boundary  Loop through all child BP's and
      ! compare to corners of parent node.

      implicit none

      ! input
      character(len=*),                intent(in)    :: side_name
      real, dimension(:), allocatable, intent(out)   :: i_bry, j_bry  ! final lists local to parent.
      integer,                         intent(inout) :: side_total    ! number of child bry points within parent
      integer,                         intent(out)   :: c_str         ! full child bry start index falling in parent
      integer,                         intent(out)   :: dim_size      ! size of global child boundary
      integer,                         intent(in)    :: vtype         ! r, u or v type

      ! local
      integer               :: varid, ij, l, counter
      real                  :: sx, sy                                 ! shift in x & y for different var type
      character(len=30)     :: vname
      integer, dimension(2) :: dimids

      real, dimension(:), allocatable :: i_in_bry, j_in_bry           ! store read in child boundary

      if (ncid==-1) then                                              ! file not yet opened.
        ierr=nf90_open(c_bry_file, nf90_nowrite, ncid)                ! var will later be in grdname file
        if(ierr/=0) call handle_ierr(ierr,
     &             'read_child_bry:: Error opening file with child bry!')
      endif


      vname='i_' / / side_name
!      if(mynode==0) print *, 'child bry name: ', vname
      ierr=nf90_inq_varid(ncid, vname, varid)                         ! Check if contains bry
      if (ierr==nf90_noerr) then
!        if(mynode==0) write (*,*) 'Found child boundary in file for: ', vname

        ierr=nf90_inquire_variable(ncid, varid, dimids = dimids)
        ierr=nf90_inquire_dimension(ncid, dimids(1), len=dim_size)    ! bry vars only have 1 dim.
!        if(mynode==0) write (*,*) 'dim_size=', dim_size

        allocate(i_in_bry(dim_size))                                  ! allocate temporary array for all bry data.
        ierr=nf90_get_var(ncid, varid, i_in_bry)                      ! read in data
        if(ierr/=0) print *, 'bad read of bry', mynode

        allocate(j_in_bry(dim_size))
        vname='j_' / / side_name
        ierr=nf90_inq_varid(ncid, vname, varid)
        ierr=nf90_get_var(ncid, varid, j_in_bry)
                                                             ! check if in sub-domain relative to parent rho-coords
        sx=0.; sy=0.                                         ! reset var type shift
        if (vtype==up_var) sx=-0.5                           ! shift i_u to i_rho
        if (vtype==vp_var) sy=-0.5                           ! shift j_v to j_rho
                                                             ! shift remains 0 for rho vars
        counter=0
        do ij=1,dim_size                                     ! check if child bry point in parent sub-domain
                                                             ! note for boundary edges, western edge glob_istr=0
                                                             ! eatern edge and glob_iend=LLm+1. This catches child points in the ghost cell.
          if (glob_istr-0.5   <= i_in_bry(ij)+sx .and.       ! child point within parent node
     &        glob_jstr-0.5   <= j_in_bry(ij)+sy .and.
     &        i_in_bry(ij)+sx <  glob_iend+0.5   .and.       ! Catch +/- 0.5 to catch between parent nodes
     &        j_in_bry(ij)+sy <  glob_jend+0.5        ) then ! <= and > to catch child exactly on mid-points between parents

            counter=counter+1
          endif

          if (counter==1) c_str=ij                           ! Store starting child index

        enddo
        side_total=counter

        if (side_total > 0) then
          allocate(i_bry( side_total ))                      ! Store both rows/columns in one
          allocate(j_bry( side_total ))                      ! array using larger of two

          i_bry(1:side_total) = i_in_bry( c_str:c_str+side_total-1 ) ! -1 for closed interval counting
          j_bry(1:side_total) = j_in_bry( c_str:c_str+side_total-1 )

!          write(*,'(/1x,I2,A)') mynode, ' Child points within parent node' ! debug
!          print *, mynode, 'side_total=', side_total         ! debug
!          print *, mynode, 'c_str=',c_str                    ! debug
!          print *, mynode, 'i_bry(str)= ', i_bry(1)          ! debug
!          print *, mynode, 'i_bry(end)= ', i_bry(side_total) ! debug
!          print *, mynode, 'j_bry(str)= ', j_bry(1)          ! debug
!          print *, mynode, 'j_bry(end)= ', j_bry(side_total) ! debug
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

      ! rho cffs:
      if (allocated(ir_south)) then
        call calc_bry_coeffs(r_tot(1), ir_south, jr_south, south_r_cffs,
     &                       isr_sw,   jsr_sw,   rp_var,   rp_var        )
      endif
      if (allocated(ir_west))  then
        call calc_bry_coeffs(r_tot(2), ir_west,  jr_west,  west_r_cffs,
     &                       iwr_sw,   jwr_sw,   rp_var,   rp_var        )
      endif
      if (allocated(ir_north)) then
        call calc_bry_coeffs(r_tot(3), ir_north, jr_north, north_r_cffs,
     &                       inr_sw,   jnr_sw,   rp_var,   rp_var        )
      endif
      if (allocated(ir_east))  then
        call calc_bry_coeffs(r_tot(4), ir_east,  jr_east,  east_r_cffs,
     &                       ier_sw,   jer_sw,   rp_var,   rp_var        )
      endif

      ! u cffs:
      if (allocated(iu_south)) then
        call calc_bry_coeffs(u_tot(1), iu_south, ju_south, south_u2u_cffs,        ! cffs to interp parent u to child u
     &                       isu2u_sw, jsu2u_sw, up_var,   up_var        )
        call calc_bry_coeffs(u_tot(1), iu_south, ju_south, south_v2u_cffs,        ! cffs to interp parent v to child u
     &                       isv2u_sw, jsv2u_sw, vp_var,   up_var        )
        call calc_bry_coeffs(u_tot(1), iu_south, ju_south, south_r2u_cffs,        ! cffs to interp parent r to child u (angle)
     &                       isr2u_sw, jsr2u_sw, rp_var,   up_var        )

        call calc_parent_vs_child_grid_angles(   u_tot(1), u_c_str(1),            ! don't need south_r2u after this, should deallocate
     &       south_r2u_cffs, isr2u_sw, jsr2u_sw, us_angcos,us_angsin, up_var, 1 ) ! 1=S
      endif
      if (allocated(iu_west)) then
        call calc_bry_coeffs(u_tot(2), iu_west,  ju_west,  west_u2u_cffs,
     &                       iwu2u_sw, jwu2u_sw, up_var,   up_var        )
        call calc_bry_coeffs(u_tot(2), iu_west,  ju_west,  west_v2u_cffs,
     &                       iwv2u_sw, jwv2u_sw, vp_var,   up_var        )
        call calc_bry_coeffs(u_tot(2), iu_west,  ju_west,  west_r2u_cffs,
     &                       iwr2u_sw, jwr2u_sw, rp_var,   up_var        )
        call calc_parent_vs_child_grid_angles(   u_tot(2), u_c_str(2),
     &       west_r2u_cffs,  iwr2u_sw, jwr2u_sw, uw_angcos,uw_angsin, up_var, 2 ) ! 2=W
      endif
      if (allocated(iu_north)) then
        call calc_bry_coeffs(u_tot(3), iu_north, ju_north, north_u2u_cffs,
     &                       inu2u_sw, jnu2u_sw, up_var,   up_var        )
        call calc_bry_coeffs(u_tot(3), iu_north, ju_north, north_v2u_cffs,
     &                       inv2u_sw, jnv2u_sw, vp_var,   up_var        )
        call calc_bry_coeffs(u_tot(3), iu_north, ju_north, north_r2u_cffs,
     &                       inr2u_sw, jnr2u_sw, rp_var,   up_var        )
        call calc_parent_vs_child_grid_angles( u_tot(3),   u_c_str(3),
     &       north_r2u_cffs, inr2u_sw, jnr2u_sw, un_angcos,un_angsin, up_var, 3 )
      endif
      if (allocated(iu_east)) then
        call calc_bry_coeffs(u_tot(4), iu_east,  ju_east,  east_u2u_cffs,
     &                       ieu2u_sw, jeu2u_sw, up_var,   up_var        )
        call calc_bry_coeffs(u_tot(4), iu_east,  ju_east,  east_v2u_cffs,
     &                       iev2u_sw, jev2u_sw, vp_var,   up_var        )
        call calc_bry_coeffs(u_tot(4), iu_east,  ju_east,  east_r2u_cffs,
     &                       ier2u_sw, jer2u_sw, rp_var,   up_var        )
        call calc_parent_vs_child_grid_angles(   u_tot(4), u_c_str(4),
     &       east_r2u_cffs,  ier2u_sw, jer2u_sw, ue_angcos,ue_angsin, up_var, 4 )
      endif

      rt2=sqrt(2.0)  ! set interpolation coefficient once

!      call MPI_Barrier(MPI_COMM_WORLD, ierr) ! debug
!      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)

      end subroutine calc_child_interpolation_coeffs  !]

! ----------------------------------------------------------------------
      subroutine calc_bry_coeffs( bry_tot, i_bry, j_bry,   cffs,      ![
     &                            i_sw,    j_sw,  p_vtype, c_vtype )

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
      ! For u & v boundaries, we need to get both parent's u, v and angle to the child
      ! u or v point in order to rotate u and v to extract the child u & v (since axis
      ! orientation is different).
      ! Therefore we also calculate coefficients for:
      ! - u vars:
      ! -- u parent to u child
      ! -- v parent to u child
      ! -- r parent to u child (for parent angle at child u)
      ! - v vars:
      ! -- v parent to v child
      ! -- u parent to v child
      ! -- r parent to v child (for parent angle at child v)
      !
      ! In order to calculate all coefficients using this same routine,
      ! We take the the child ij values (e.g. south_u2u_cffs), which are in terms of
      ! parent rho, u or v coordinates, and shift them in terms of their 4 surrounding
      ! rho, u or v parent points using sx and sy by +/- 0.5.
      !
      ! If the parent type is the same as child type then sx=0 and sy=0 (no shift).
      !
      !]

      implicit none

      ! hidden_mpi_vars.h needed for inode & jnode to establish if boundary edge
# include "hidden_mpi_vars.h"

      ! input/output
      integer,                              intent(in)  :: bry_tot     ! total child bry points in parent
      real,    dimension(bry_tot),          intent(in)  :: i_bry
      real,    dimension(bry_tot),          intent(in)  :: j_bry
      real,    dimension(:,:), allocatable, intent(out) :: cffs        ! 4 parent interpolants for each child point
      integer, dimension(:),   allocatable, intent(out) :: i_sw, j_sw  ! i and j of sw parent interpolant
      integer,                              intent(in)  :: p_vtype     ! parent interp var type
      integer,                              intent(in)  :: c_vtype     ! child bry point type

      ! local
      integer :: ij
      real    :: ri, rj                       ! position fractions
      real    :: sx, sy                       ! shift for getting u-ij in terms of v-ij

      allocate( cffs(4,bry_tot) )             ! bry_tot for each bry points in parent. 4 for all four interpolants
      allocate( i_sw(bry_tot) )               ! i_sw is i grid point of sw parent interpolant
      allocate( j_sw(bry_tot) )

!      if(mynode==4) then
!        print *,'calcB FI-glob_istr=',floor(i_bry(1)),glob_istr,mynode
!        print *,'calcB FJ-glob_jstr=',floor(j_bry(1)),glob_jstr,mynode
!      endif
                                                           ! Shift child type to parent type
      sx=0.; sy=0.                                         ! if parent & child same type no shift
      if     (p_vtype==up_var .and. c_vtype==vp_var) then
        sx=+0.5                                            ! shift i_rho to i_u
        sy=-0.5                                            ! shift j_v   to j_rho
      elseif (p_vtype==vp_var .and. c_vtype==up_var) then
        sx=-0.5                                            ! shift i_u   to i_rho
        sy=+0.5                                            ! shift j_rho to j_v
      elseif (p_vtype==rp_var .and. c_vtype==up_var) then
        sx=-0.5                                            ! shift i_u   to i_rho
        sy= 0.0                                            ! both  j_rho no shift
      elseif (p_vtype==rp_var .and. c_vtype==vp_var) then
        sx= 0.0                                            ! both  i_rho no shift
        sy=-0.5                                            ! shift j_v   to j_rho
      endif

      do ij=1,bry_tot                                      ! Since glob_istr=0 on western_edge to catch child in ghost row in assign_child_interpolation_coeffs
                                                           ! we need to treat W & S edges diferently.
        if (inode==0) then                                 ! western_edge
          i_sw(ij) = floor(i_bry(ij)+sx) - glob_istr       ! local coords, set the SW parent i position. western edge glob_istr=0 not 1 hence no +1
        else
          i_sw(ij) = floor(i_bry(ij)+sx) - glob_istr + 1   ! +1 since e.g. glob_jstr=51, j_bry=51.23 want sw corner to be 1 not 0
        endif

        if (jnode==0) then                                 ! eastern_edge
          j_sw(ij) = floor(j_bry(ij)+sy) - glob_jstr       ! local coords, set the SW parent j position
        else
          j_sw(ij) = floor(j_bry(ij)+sy) - glob_jstr + 1
        endif

        ri = i_bry(ij)+sx - floor(i_bry(ij)+sx)            ! i fraction (e.g. 23.47 - 23 = 0.47)
        rj = j_bry(ij)+sy - floor(j_bry(ij)+sy)            ! j fraction

        cffs(1,ij) = (1-ri)*(1-rj)                         ! sw parent interp
        cffs(2,ij) =    ri *(1-rj)                         ! se parent interp
        cffs(3,ij) = (1-ri)*   rj                          ! nw parent interp
        cffs(4,ij) =    ri *   rj                          ! ne parent interp

!        if(mynode==1) then ! debug
!          if(ij<3) then
!            print *, 'ij=',ij,'i_bry(ij)=',i_bry(ij),'j_bry(ij)=',j_bry(ij)
!            print *, 'i_sw(ij)=',i_sw(ij),'j_sw(ij)=',j_sw(ij)
!            print *, 'ri=',ri,'rj=',rj
!            print *, 'cff_sum=', cffs(1,ij)+cffs(2,ij)+cffs(3,ij)+cffs(4,ij)
!          endif
!        endif

      enddo

      end subroutine calc_bry_coeffs  !]

! ----------------------------------------------------------------------
      subroutine calc_parent_vs_child_grid_angles( bry_tot, c_str,  cffs,     i_sw, j_sw,    ![
     &                                             anglcos, anglsin, c_vtype, bry   )

      ! Calcluate the angle between the parent and child and grid in order to
      ! rotate u/v in parent axes to u/v in child axes.
      ! Can use 1D array since doesn't change through vertical layers.
      ! Store once to avoid expensive trig functions.

      implicit none

      ! input/output
      integer,                              intent(in)  :: bry_tot     ! total child bry points in parent
      integer,                              intent(in)  :: c_str       ! full child bry start index falling in parent
      integer, dimension(:),                intent(in)  :: i_sw, j_sw  ! i and j of sw parent interpolant
      real,    dimension(:,:), allocatable, intent(in)  :: cffs        ! 4 parent interpolants for each child point
      real   , dimension(:),   allocatable, intent(out) :: anglcos     ! interpolated child angle at child point
      real   , dimension(:),   allocatable, intent(out) :: anglsin
      integer,                              intent(in)  :: c_vtype     ! child bry point type
      integer,                              intent(in)  :: bry         ! 1/2/3/4 = S/W/N/E

      ! local
      integer                         :: ij, varid
      real, dimension(:), allocatable :: c_angle                       ! store read in child angles
      real                            :: d_angle
      character(len=30)               :: vname


      allocate( c_angle( bry_tot ))
      allocate( anglcos( bry_tot ))
      allocate( anglsin( bry_tot ))

      if     (c_vtype==up_var) then
        vname='angle_u_' / / bry_name(bry)                             ! append with S/W/E/N
      elseif (c_vtype==vp_var) then
        vname='angle_v_' / / bry_name(bry)
      endif

      ierr=nf90_inq_varid(ncid, vname, varid)                          ! Check if contains bry
      if (ierr==nf90_noerr) then
        ierr=nf90_get_var( ncid, varid, c_angle, start=(/c_str/), count=(/bry_tot/) ) ! only read child points in subdomain
        if (ierr/=0) call handle_ierr(ierr, 'Bad read of ', vname)
      else
        write (*,'(/1x,3A)') 'ERROR: No child boundary in parent ',    ! angle not contained but should be because
     &                       'domain for: ', vname                     ! side has been found already
        error stop
      endif


      do ij=1,bry_tot                                                  ! angle = p_angl - c_angl
        d_angle = cffs(1,ij) * angler(i_sw(ij)  , j_sw(ij)  )         ! interp from sw parent point.
     &          + cffs(2,ij) * angler(i_sw(ij)+1, j_sw(ij)  )         ! interp from se parent point.
     &          + cffs(3,ij) * angler(i_sw(ij)  , j_sw(ij)+1)         ! interp from nw parent point.
     &          + cffs(4,ij) * angler(i_sw(ij)+1, j_sw(ij)+1)         ! interp from ne parent point.
     &          - c_angle(ij)                                         ! child angle is at point (no interp needed)

        anglcos(ij) = cos(d_angle)
        anglsin(ij) = sin(d_angle)

!        if (ij==1 .and. mynode==5) then
!          print *, 'angler(i_sw,j_sw)',angler(i_sw(ij),j_sw(ij))
!          print *, 'ij=',ij,'c_angle(ij)',c_angle(ij)
!          print *, 'angl(ij)=',angl(ij)
!          print *, 'anglcos(ij)',anglcos(ij)
!          print *, 'anglsin(ij)',anglsin(ij)
!        endif

      enddo

      end subroutine calc_parent_vs_child_grid_angles  !]

! ----------------------------------------------------------------------
      subroutine interpolate_child_bry( var, nmax )  ![

      ! Interpolate all boundaries if applicable

      implicit none

      integer, intent(in) :: nmax
      real   , intent(in) :: var(GLOBAL_2D_ARRAY,nmax)  ! variable to interpolate.
                                                        ! wasn't able to use dimension(:,:,:), but vars start from -1 index so maybe this is safer.

      if (allocated(ir_south))
     &  call interpolate_one_child_bry(r_tot(1), south_r_cffs, isr_sw, jsr_sw,
     &                                 nmax,     var,          vsouth, rmask  )

      if (allocated(ir_west))
     &  call interpolate_one_child_bry(r_tot(2),  west_r_cffs, iwr_sw, jwr_sw,
     &                                 nmax,      var,         vwest,  rmask  )

!      if(mynode==1) then  ! debug
!        print *, 'vwest: r_tot=', r_tot(2),'mynode',mynode
!      print *,'i/j(end)=',ir_west(r_tot(2)),jr_west(r_tot(2))
!        print *,'i/jsw(end)=',iwr_sw(r_tot(2)),jwr_sw(r_tot(2))
!        print *,'vwest(end)=',vwest(r_tot(2),1)
!      endif
!      if(mynode==4) then  ! debug
!        print *, 'vwest: r_tot=', r_tot(2),'mynode',mynode
!        print *,'i/j(1)=',ir_west(1),jr_west(1)
!        print *,'i/jsw(1)=',iwr_sw(1),jwr_sw(1)
!        print *,'vwest(1)=',vwest(1,1)
!      endif

      if (allocated(ir_north))
     &  call interpolate_one_child_bry(r_tot(3), north_r_cffs, inr_sw, jnr_sw,
     &                                 nmax,     var,          vnorth, rmask  )

      if (allocated(ir_east))
     &  call interpolate_one_child_bry(r_tot(4), east_r_cffs,  ier_sw, jer_sw,
     &                                 nmax,     var,          veast,  rmask  )


      end subroutine interpolate_child_bry  !]

! ----------------------------------------------------------------------
      subroutine interpolate_u_child_bry( uvar, vvar, nmax )  ![
      ! Interpolate all u boundaries if applicable

      implicit none


      integer, intent(in) :: nmax
      real   , intent(in) :: uvar(GLOBAL_2D_ARRAY,nmax)  ! variable to interpolate.
      real   , intent(in) :: vvar(GLOBAL_2D_ARRAY,nmax)  ! variable to interpolate.


      if (allocated(iu_south)) then
        call interpolate_one_child_bry(u_tot(1), south_u2u_cffs, isu2u_sw,  jsu2u_sw,
     &                                 nmax,     uvar,           us_uprnt,  umask   )
        call interpolate_one_child_bry(u_tot(1), south_v2u_cffs, isv2u_sw,  jsv2u_sw,
     &                                 nmax,     vvar,           us_vprnt,  umask   )

        call interpolate_one_u_bry(    u_tot(1), nmax, usouth,   us_uprnt,  us_vprnt,
     &                                 us_angcos,us_angsin                          )
      endif

      if (allocated(iu_west)) then
        call interpolate_one_child_bry(u_tot(2), west_u2u_cffs,  iwu2u_sw,  jwu2u_sw,
     &                                 nmax,     uvar,           uw_uprnt,  umask   )
        call interpolate_one_child_bry(u_tot(2), west_v2u_cffs,  iwv2u_sw,  jwv2u_sw,
     &                                 nmax,     vvar,           uw_vprnt,  umask   )

        call interpolate_one_u_bry(    u_tot(2), nmax, uwest,    uw_uprnt,  uw_vprnt,
     &                                 uw_angcos,uw_angsin                          )
      endif

      if (allocated(iu_north)) then
        call interpolate_one_child_bry(u_tot(3), north_u2u_cffs, inu2u_sw,  jnu2u_sw,
     &                                 nmax,     uvar,           un_uprnt,  umask   )
        call interpolate_one_child_bry(u_tot(3), north_v2u_cffs, inv2u_sw,  jnv2u_sw,
     &                                 nmax,     vvar,           un_vprnt,  umask   )

        call interpolate_one_u_bry(    u_tot(3), nmax, unorth,   un_uprnt,  un_vprnt,
     &                                 un_angcos,un_angsin                          )
      endif

      if (allocated(iu_east)) then
        call interpolate_one_child_bry(u_tot(4), east_u2u_cffs,  ieu2u_sw,  jeu2u_sw,
     &                                 nmax,     uvar,           ue_uprnt,  umask   )
        call interpolate_one_child_bry(u_tot(4), east_v2u_cffs,  iev2u_sw,  jev2u_sw,
     &                                 nmax,     vvar,           ue_vprnt,  umask   )

        call interpolate_one_u_bry(    u_tot(4), nmax, ueast,    ue_uprnt,  ue_vprnt,
     &                                 ue_angcos,ue_angsin                          )
      endif

      end subroutine interpolate_u_child_bry  !]

! ----------------------------------------------------------------------
      subroutine interpolate_one_child_bry(bry_tot, rho_cffs, i_sw, j_sw,  ![
     &                                     nmax,    var,      vbry, mask )

      ! Interpolation coefficients already calculated in calc_child_interpolation_coeffs

      implicit none

      ! inputs
      integer,                              intent(in)  :: bry_tot     ! total child bry points in parent
      real,    dimension(:,:),              intent(in)  :: rho_cffs    ! 4 parent interpolants for each child point
      integer, dimension(:),                intent(in)  :: i_sw, j_sw  ! i and j of sw parent interpolant
      integer,                              intent(in)  :: nmax        ! vertical layers of var (1 if 2D)
      real,                                 intent(in)  :: var(GLOBAL_2D_ARRAY,nmax)  ! parent var to interpolate from
      real,    dimension(:,:), allocatable, intent(out) :: vbry        ! bry to interpolated 2D/3D var
      real,                                 intent(in)  :: mask(GLOBAL_2D_ARRAY)      ! grid type r/u/v mask

      ! local
      integer                 :: ij, k
      real, dimension(nmax)   :: var_sw, var_se, var_nw, var_ne  ! temporary copies of columns of data
      real(kind=8), parameter :: spv_set=1.D+33

      if (.not. allocated(vbry)) allocate( vbry(bry_tot, N))     ! allocate once for 3D var and reuse for all 2D/3D vars
                                                                 ! maybe move this to initialization to avoid if statement repetition

      do ij=1,bry_tot           ! Normally would have k as outer loop for cache coherence but can avoid many if statements by swapping order here:

        if(i_sw(ij) >= 0) then  ! child point not masked. >=0 as mask value is -1. i_sw can be 0 if point is between nodes.

          if ( mask( i_sw(ij)  , j_sw(ij)  ) == 0) then
            call interp_masked_var( i_sw(ij)  , j_sw(ij)  , nmax, var, var_sw, mask )
          else
            var_sw = var(i_sw(ij)  , j_sw(ij)  , :)
          endif

          if ( mask( i_sw(ij)+1, j_sw(ij)  ) == 0) then
            call interp_masked_var( i_sw(ij)+1, j_sw(ij)  , nmax, var, var_se, mask )
          else
            var_se = var(i_sw(ij)+1, j_sw(ij)  , :)
          endif

          if ( mask( i_sw(ij)  , j_sw(ij)+1) == 0) then
            call interp_masked_var( i_sw(ij)  , j_sw(ij)+1, nmax, var, var_nw, mask )
          else
            var_nw = var(i_sw(ij)  , j_sw(ij)+1, :)
          endif

          if ( mask( i_sw(ij)+1, j_sw(ij)+1) == 0) then
            call interp_masked_var( i_sw(ij)+1, j_sw(ij)+1, nmax, var, var_ne, mask )
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
      subroutine interpolate_one_u_bry( bry_tot, nmax, ubry, uparent, vparent,  ![
     &                                  angcos,  angsin )

      implicit none

      ! inputs
      integer,                              intent(in)  :: bry_tot     ! total child bry points in parent
      integer,                              intent(in)  :: nmax        ! vertical layers of var (1 if 2D)
      real,    dimension(:,:), allocatable, intent(out) :: ubry        ! bry to interpolated 2D/3D var
      real,    dimension(:,:),              intent(in)  :: uparent     ! u-value of parent grid at child u-point
      real,    dimension(:,:),              intent(in)  :: vparent     ! v-value of parent grid at child u-point
      real,    dimension(:),                intent(in)  :: angcos      ! cos of angle difference between parent-child grids
      real,    dimension(:),                intent(in)  :: angsin      ! sin of ...

      ! local
      integer :: ij, k

      if (.not. allocated(ubry)) allocate( ubry(bry_tot, N))   ! allocate once for 3D var and reuse for all 2D/3D vars

      do ij=1,bry_tot
        do k=1,nmax
          ubry(ij,k) = uparent(ij,k) * angcos(ij) - vparent(ij,k) * angsin(ij)
        enddo
      enddo

      end subroutine interpolate_one_u_bry  !]

! ----------------------------------------------------------------------
      subroutine create_child_output_file  ![

      ! For now decided each file still creates file even if not child contained in node.
      ! so that ncjoin can still check files are correct.
      ! define boundary even if not contained in parent for joining tool.
      !
      ! On a boundary, the r, u and v points of the child are in different positions (obviously),
      ! This means that the number of r, u and v points in a sub-domain might not be the
      ! same, since the child bry is likely at an angle to the parent boundary.
      ! This means we need dimensions for u and v even though they might sit on the
      ! same points in the joined file.
      !
      !   <-----u------
      !       r        |      Assume rho & u are from the child west boundary.
      !            u   |      They would both has j coordinates in eta_rho.
      !          r     |      However, there are only 2 u's in the parent sub-domain
      !               u|      and 3 rho's. Therefore in the output netcdf file, they
      !             r  |      cannot share an 'eta_rho_west' sub-domain, as their lengths
      !                       are not the same.
      !
      ! Therefore, we need to give u and v bry variables their own xi_rho, and eta_rho
      ! dimensions, e.g. eta_urho_west.

      implicit none

      ! local
      integer               :: prev_fill_mode
      integer, dimension(5) :: r_dimid_bry,  u_dimid_bry,  v_dimid_bry  ! dimension ID for boundaries
      integer               :: itrc                                     ! tracer index

      if(total_output_rec == 0) then                                    ! Add root output name to
        file_name = trim(output_root_name) / / trim(file_name)          ! file name (only done once)
      endif


      call ncdf_create_blank_file(file_name, ncid_out,  prev_fill_mode,
     &                              recs_per_file, total_output_rec       )

      ! DEFINE DIMENSIONS & VARIABLES:

      ierr=nf90_put_att(ncid_out, NF90_GLOBAL, 'global_south_size', glob_c_rs)
      ierr=nf90_put_att(ncid_out, NF90_GLOBAL, 'global_west_size' , glob_c_rw)
      ierr=nf90_put_att(ncid_out, NF90_GLOBAL, 'global_north_size', glob_c_rn)
      ierr=nf90_put_att(ncid_out, NF90_GLOBAL, 'global_east_size' , glob_c_re)
      ierr=nf90_def_dim( ncid_out, 'rec_time', 0, dimid_t)

      ! rho dims:
      ierr=nf90_def_dim( ncid_out, 'xi_rho_south', max( r_tot(1), 1 ), r_dimid_bry(1) )  ! Set size as 1 if size==0 to avoid unlimited dimension (could be confused with time in ncjoin)
      ierr=nf90_def_dim( ncid_out, 'eta_rho_west', max( r_tot(2), 1 ), r_dimid_bry(2) )
      ierr=nf90_def_dim( ncid_out, 'xi_rho_north', max( r_tot(3), 1 ), r_dimid_bry(3) )
      ierr=nf90_def_dim( ncid_out, 'eta_rho_east', max( r_tot(4), 1 ), r_dimid_bry(4) )
      ierr=nf90_def_dim( ncid_out, 's_rho'       , N                 , r_dimid_bry(5) )

      ! u dims: (see preamble for eta_urho_west & eta_urho_east info)
      ierr=nf90_def_dim( ncid_out, 'xi_u_south',   max( u_tot(1), 1 ), u_dimid_bry(1) )
      ierr=nf90_def_dim( ncid_out, 'eta_urho_west',max( u_tot(2), 1 ), u_dimid_bry(2) )
      ierr=nf90_def_dim( ncid_out, 'xi_u_north',   max( u_tot(3), 1 ), u_dimid_bry(3) )
      ierr=nf90_def_dim( ncid_out, 'eta_urho_east',max( u_tot(4), 1 ), u_dimid_bry(4) )
      u_dimid_bry(5) = r_dimid_bry(5)  ! k dimension of u-point is same as rho

      ! v dims:


      call define_child_bry( indxZ,  r_tot, r_c_str, r_dimid_bry, 1, .false. )
      call define_child_bry( indxUb, u_tot, u_c_str, u_dimid_bry, 1, .false. )

      ! tracers:
      do itrc=1,NT
        call define_child_bry( itrc, r_tot, r_c_str, r_dimid_bry, N, .true.  )
      enddo


      ierr=nf90_enddef(ncid_out)  ! end definition to write changes to disk
      if (mynode == 0) then
        write(*,'(/1x,2A/)')
     &     'bry_extract_child :: created new netCDF file ', file_name
      endif

      end subroutine create_child_output_file  !]

! ----------------------------------------------------------------------
      subroutine wrt_child_bry ( indx_v, nmax, trc_yn, sbry, wbry, nbry, ebry, bry_tot )  ![

      ! Each sub-domain writes out its portion of child boundary if applicable.
      ! all procs, make child bc netcdf file, only relevant procs define bc variables (east, north, etc)

      implicit none

      ! inputs
      integer,                           intent(in) :: indx_v  ! index of primary var or itrc
      integer,                           intent(in) :: nmax    ! vertical levels
      logical,                           intent(in) :: trc_yn  ! is variable a tracer y/n
      real, dimension(:,:), allocatable, intent(in) :: sbry    ! interpolated south bry values
      real, dimension(:,:), allocatable, intent(in) :: wbry    ! interpolated west  bry values
      real, dimension(:,:), allocatable, intent(in) :: nbry    ! interpolated north bry values
      real, dimension(:,:), allocatable, intent(in) :: ebry    ! interpolated east  bry values
      integer, dimension(4),             intent(in) :: bry_tot ! total child bry points in parent

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


      bname = trim(rname) / /'_south'               ! Move this inside if allocated...
      ierr  = nf90_inq_varid(ncid_out, bname, varid )
      count_out(1) = bry_tot(1)
      if (allocated(sbry)) then
        ierr=nf90_put_var( ncid_out, varid , sbry(:,:),
     &                     start=start_out(1:ndims), count=count_out(1:ndims) )
        if(ierr/=0) write(*,*) 'wrts_bry ierr= ',nf90_strerror(ierr), mynode
      else

        ierr=nf90_put_var( ncid_out, varid, (/1.0/),                               ! put dummy value so ncjoin doesn't complain
     &                     start=start_out(1:ndims), count=count_dummy(1:ndims) )
      endif

      bname = trim(rname) / /'_west'
      ierr=nf90_inq_varid(ncid_out, bname, varid )
      count_out(1) = bry_tot(2)
      if (allocated(wbry )) then
        ierr=nf90_put_var( ncid_out, varid    , wbry(:,:) ,
     &                     start=start_out(1:ndims), count=count_out(1:ndims) )
        if(ierr/=0) write(*,*) 'wrtw_bry ierr= ',nf90_strerror(ierr), mynode
      else
        ierr=nf90_put_var( ncid_out, varid, (/1.0/),
     &                     start=start_out(1:ndims), count=count_dummy(1:ndims) )
      endif

      bname = trim(rname) / /'_north'
      ierr=nf90_inq_varid(ncid_out, bname, varid )
      count_out(1) = bry_tot(3)
      if (allocated(nbry)) then
        ierr=nf90_put_var( ncid_out, varid    , nbry(:,:),
     &                     start=start_out(1:ndims), count=count_out(1:ndims) )
        if(ierr/=0) write(*,*) 'wrtn_bry ierr= ',nf90_strerror(ierr), mynode
      else
        ierr=nf90_put_var( ncid_out, varid, (/1.0/),
     &                     start=start_out(1:ndims), count=count_dummy(1:ndims) )
      endif

      bname = trim(rname) / /'_east'
      ierr=nf90_inq_varid(ncid_out, bname, varid )
      count_out(1) = bry_tot(4)
      if (allocated(ebry )) then
        ierr=nf90_put_var( ncid_out, varid    , ebry(:,:) ,
     &                     start=start_out(1:ndims), count=count_out(1:ndims) )
        if(ierr/=0) write(*,*) 'wrte_bry ierr= ',nf90_strerror(ierr), mynode
      else
        ierr=nf90_put_var( ncid_out, varid, (/1.0/),
     &                     start=start_out(1:ndims), count=count_dummy(1:ndims) )
      endif

      end subroutine wrt_child_bry  !]

! ----------------------------------------------------------------------
      subroutine define_child_bry( indx_v, sides_tot, sides_c_str,  ![
     &                             dimid_bry, nmax, trc_yn )

      ! Define child bry variable in output netCDF file.

      implicit none

      ! input/outputs
      integer,               intent(in) :: indx_v
      integer, dimension(4), intent(in) :: sides_tot
      integer, dimension(4), intent(in) :: sides_c_str
      integer, dimension(5), intent(in) :: dimid_bry     ! dimension ID for boundaries
      integer,               intent(in) :: nmax
      logical,               intent(in) :: trc_yn

      ! local
      integer               :: varid, bry, ndims
      character(len=40)     :: bname, rname
      character(len=60)     :: lname, uname
      integer, dimension(3) :: vdimids


      if (trc_yn) then                                   ! handle primary variable or tracers name
        rname=t_vname(indx_v)
        lname=t_lname(indx_v)
        uname=t_units(indx_v)
      else
        rname=vname(1,indx_v)
        lname=vname(2,indx_v)
        uname=vname(3,indx_v)
      endif

      do bry=1,4                                         ! loop through 4 boundaries:

        if (nmax == 1) then                              ! handle 2D or 3D variable dimids
          ndims=2
          vdimids(1:ndims) = (/dimid_bry(bry), dimid_t/)
        else
          ndims=3
          vdimids(1:ndims) = (/dimid_bry(bry), dimid_bry(5), dimid_t/)
        endif


        bname=trim(rname) / /'_'/ /trim( bry_name(bry) )

        call nc_define_var(ncid_out, bname, lname, uname, ! longname, units
     &                     vdimids(1:ndims), ierr )


        ierr=nf90_inq_varid(ncid_out, bname, varid)

        ierr=nf90_put_att(ncid_out, varid, 'start', sides_c_str( bry ))  ! add attribute with global child start value (i.e. south_c_str)
        ierr=nf90_put_att(ncid_out, varid, 'count', sides_tot(   bry ))  ! and its count (i.e. south_tot). This will make it easy for ncjoin.
        ierr=nf90_put_att(ncid_out, varid, 'child_bry', 'True')                        ! 'child_bry' for ncjoin to identify variable. 'True' just a filler.

      enddo

      end subroutine define_child_bry  !]

! ----------------------------------------------------------------------
      subroutine interp_masked_var( i, j, nmax, var, vinterp, mask )  ![

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
      real,                  intent(in)  :: mask(GLOBAL_2D_ARRAY)      ! grid type r/u/v mask

      ! local
      real    :: ca, cc, cf, ch, mb, md, me, mg  ! middle & corner interpolants as per diagram above
      real    :: sum_cffs                        ! sum of weightings (denominator in interpolation)
      integer :: k

      ca = rt2*mask(i-1,j-1)  ! corner point weightings BUG THIS SHOULD BE RMASK
      cc = rt2*mask(i+1,j-1)  ! (calculate once outside of vertical loop)
      cf = rt2*mask(i-1,j+1)
      ch = rt2*mask(i+1,j+1)

      mb = mask(i  ,j-1)      ! middle point weightings
      md = mask(i-1,j  )
      me = mask(i+1,j  )
      mg = mask(i  ,j+1)

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
