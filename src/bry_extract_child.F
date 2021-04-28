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
      !    The i/j value is -1 over a land mask and -2 if not in parent domain (for rho points).
      !
      !    Seen in DIAGRAM (A) below, each child point (a-f) will have a input variable
      !    giving it's i'th position, i.e. variable 'rp%i_west', which for point b = 1.4.
      !    Another variable will hold the j'th position, i.e. 'rp%j_west' which for point b = 1.2.
      !    If the boundary is closed then i and j variables will not exist for that boundary.
      !
      !         DIAGRAM (A)
      !
      !   C  ^        u       ^  D
      ! (1,2)|              v |(2,2)
      !      X-------rc-------X-->         A-D = parent grid points (X)
      !      |                |
      !      |      ue        |            a-f = child boundary grid points (West boundary)
      ! N    |            vg  |
      ! ^    |     rb         |
      ! |    |                |            Fractional i/j coordinates of child
      ! |    |    ud          |            relative to parent grid:
      !      |          vf    |            - r-vars:
      !      |   ra           |              c = (1.2 , 2.0)   b = (1.3 , 1.6)  a = (1.4, 1.2)
      !      |                |
      !      X--u-------------X-->         - u vars:
      ! (1,1)                  (2,1)         d = (1.25, 2.0)   e = (1.35, 1.6)
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
      !]

#include "cppdefs.h"
#ifdef BRY_EXTRACT_CHILD

      use dimensions                                   ! has dimensions, a list of forcing files, and rmask, glob_istr, etc
      use read_write
      use netcdf
      use mpi                                          ! debug temporary for mpi abort
      use tracers, only: t, t_vname, t_lname, t_units  ! need to get names of tracers

      implicit none
      private

#include "ocean2d.h"
#include "ocean3d.h"
! scalars.h for 'knew' var
#include "scalars.h"


      ! ****************************************************************
      ! USER INPUTS
      integer            :: rec_rate      = 2                                            ! output record rate for interpolated child boundaries
      integer            :: recs_per_file = 2                                            ! total recs per output file
      integer, parameter :: nc_grds       = 2                                            ! number of child grids required
      character(len=50)  :: c_bry_files(nc_grds) = (/ 'sample_child_ij_bry.nc',
     &                                                'sample_child_ij_bry_2nd.nc' /)    ! child boundary input file
      character(len=60)  :: file_names(nc_grds)  = (/ '_child_bry.nc',
     &                                                '_child_bry_2.nc' /)               ! suffix to append to simulation output rootname
      ! END USER INPUTS
      ! ****************************************************************

      ! DECLARE VARIABLES: ![

      type grid_type_vars  ! contains all interpolation record keeping variables for rho, u and v points.
                           ! would be cleaner to have variables 1-4 for W/S/E/N rather than loads of variables.

        ! input vars:
        integer                              :: glob_cs, glob_cw, glob_cn, glob_ce    ! total rho,u or v number of global child points along child boundary
                                                                                      ! not necessary to store all of these but easier for subroutine arguements
        integer, dimension(4)                :: tot=0                                 ! number of child points within parent on each boundary. 4 = s/w/n/e, 2 = 2 rows of i/j input
        integer, dimension(4)                :: c_str=0                               ! start index of each child boundary within parent
        real,    dimension(:),   allocatable :: i_south, i_west, i_north, i_east      ! i position of child bry in parent grid coords
        real,    dimension(:),   allocatable :: j_south, j_west, j_north, j_east      ! j position of child bry in parent grid coords

        ! interp vars:
        ! -- r2r, u2u, or v2v:
        real,    dimension(:,:), allocatable :: s_cffs,  w_cffs, n_cffs,  e_cffs      ! 4 interpolation coefficients to interp from 4 parent points to corresponding (r2r, u2u, v2v) surronding child point.
                                                                                      ! i=each bry point. j=1-4 (4 interpolants).
        integer, dimension(:),   allocatable :: is_sw,   js_sw,  iw_sw,   jw_sw,      ! sw parent interpolant i and j position
     &                                          in_sw,   jn_sw,  ie_sw,   je_sw       ! could merge with coefficients.

        ! -- d2d (dim 2 dim) for v2u (u-child) and u2v (v-child):
        real,    dimension(:,:), allocatable :: s_d2d_cffs, w_d2d_cffs,                   ! 4 interpolation coefficients to interp from 4 parent v/u points surronding child v/u point.
     &                                          n_d2d_cffs, e_d2d_cffs                    ! i=each bry point. j=1-4 (4 interpolants).
        integer, dimension(:),   allocatable :: isd2d_sw,   jsd2d_sw, iwd2d_sw, jwd2d_sw, ! v or u parent 2 u or v child sw parent interpolant i and j position
     &                                          ind2d_sw,   jnd2d_sw, ied2d_sw, jed2d_sw

        ! -- r2d (rho 2 dim) for r2u (u-child) and r2v (v-child)::
        real,    dimension(:,:), allocatable :: s_r2d_cffs, w_r2d_cffs,                   ! 4 interpolation coefficients to interp from 4 parent r-points surronding child u or v point. (angle)
     &                                          n_r2d_cffs, e_r2d_cffs
        integer, dimension(:),   allocatable :: isr2d_sw,   jsr2d_sw, iwr2d_sw, jwr2d_sw,
     &                                          inr2d_sw,   jnr2d_sw, ier2d_sw, jer2d_sw

        ! -- trig of angle difference between parent and child grids (up & vp only):
        real,    dimension(:),   allocatable :: s_angcos,   w_angcos, n_angcos, e_angcos, ! calculated rotation required at child bry point to get parent u/v to child u/v
     &                                          s_angsin,   w_angsin, n_angsin, e_angsin  ! difference between parent grid angle at child point compared to child grid angle

        ! -- store parent u or v at child u or v to rotate:
        real,    dimension(:,:), allocatable :: s_uprnt,    w_uprnt,  n_uprnt,  e_uprnt,  ! u value of parent at child u location
     &                                          s_vprnt,    w_vprnt,  n_vprnt,  e_vprnt   ! v value of parent at child u location

        ! final bry values:
        real,    dimension(:,:), allocatable :: south,      west,     north,    east      ! final interpolated variable (2D or 3D)

      end type grid_type_vars

      type (grid_type_vars) :: rp_vars(nc_grds) ! rho variables
      type (grid_type_vars) :: up_vars(nc_grds) ! u   variables
      type (grid_type_vars) :: vp_vars(nc_grds) ! v   variables


      character(len=9), dimension(4) ::  bry_name =     (/ 'south',     'west',
     &                                                     'north',     'east'  /)
      character(len=9), dimension(4) :: rbry_name = (/ 'rho_south', 'rho_west',
     &                                                 'rho_north', 'rho_east'  /)
      character(len=9), dimension(4) :: ubry_name = (/   'u_south',   'u_west',
     &                                                   'u_north',   'u_east'  /)
      character(len=9), dimension(4) :: vbry_name = (/   'v_south',   'v_west',
     &                                                   'v_north',   'v_east'  /)


      ! INTERPOLATION:
      real                :: rt2                         ! square root of 2 stored for efficiency
                                                         ! allocated for each var so only need one.
      ! NETCDF:
      integer             :: ierr=0                      ! debug temporary
      type nc_vars
        integer           :: ncid=-1                     ! in file for child bry positions
        integer           :: ncid_out                    ! output file
        integer           :: dimid_t                     ! time dimension id in output file
        integer           :: output_rec = 0              ! record number of output. 0 indicates we need new file
        integer           :: total_output_rec = 0        ! total flux output recordings so far
        character(len=50) :: c_bry_file                  ! input  child boundary file name
        character(len=60) :: file_name                   ! output child boundary file name
      end type nc_vars

      type (nc_vars)      :: ncvars(nc_grds)             ! netcdf variables for each child grd

      public set_child_bry

      !] <- end declare variables

      contains

! ----------------------------------------------------------------------
      subroutine set_child_bry  ![
      implicit none             !  Interpolate child boundary at desired intervals.

      ! local
      integer :: itrc, c  ! tracer and child grid indices


      if (mod(iic-1,rec_rate) == 0) then  ! only plot if timestep is divisble by record rate. 1 replaces ntstart?

        do c=1,nc_grds

          if (ncvars(c)%ncid==-1) then                                                            ! first time only
            ncvars(c)%c_bry_file = c_bry_files(c)                                       ! set to user input value
            ncvars(c)%file_name  = file_names(c)
            call read_child_bry(       ncvars(c), rp_vars(c), up_vars(c), vp_vars(c) )  ! change to rp1 once working                            ! not per variable, just r,u,v types
            call calc_child_interpolation_coeffs( rp_vars(c), up_vars(c), vp_vars(c), ncvars(c) )
          endif


          if(ncvars(c)%output_rec == 0) then                                                  ! if file doesn't yet exist or need new file
            call create_child_output_file( rp_vars(c), up_vars(c), vp_vars(c), ncvars(c) )
          endif

          ierr=nf90_open(ncvars(c)%file_name, nf90_write, ncvars(c)%ncid_out)
          ncvars(c)%output_rec       = ncvars(c)%output_rec       + 1                                                   ! current file record
          ncvars(c)%total_output_rec = ncvars(c)%total_output_rec + 1                                       ! total number of recs written

          ! zeta:
          call interpolate_child_bry(   zeta(:,:,knew), 1,  rp_vars(c)            )
          call wrt_child_bry(           indxZ,  1, .false., rp_vars(c), ncvars(c) )

          ! ubar:
          call interpolate_u_child_bry( ubar(:,:,knew), vbar(:,:,knew), 1, up_vars(c) )
          call wrt_child_bry(           indxUb, 1, .false., up_vars(c),     ncvars(c) )

          ! vbar:
          call interpolate_v_child_bry( ubar(:,:,knew), vbar(:,:,knew), 1, vp_vars(c) )
          call wrt_child_bry(           indxVb, 1, .false., vp_vars(c),     ncvars(c) )

          ! u:
          call interpolate_u_child_bry( u(START_2D_ARRAY,1,nstp),
     &                                  v(START_2D_ARRAY,1,nstp), N,    up_vars(c) )
          call wrt_child_bry(           indxU,  N, .false., up_vars(c),  ncvars(c) )

          ! v:
          call interpolate_v_child_bry( u(START_2D_ARRAY,1,nstp),
     &                                  v(START_2D_ARRAY,1,nstp), N, vp_vars(c) )
          call wrt_child_bry( indxV,  N, .false., vp_vars(c),         ncvars(c) )


          ! tracers:
          do itrc=1,NT
            call interpolate_child_bry( t(START_2D_ARRAY,1,nstp,itrc), N, rp_vars(c) )
            call wrt_child_bry(         itrc, N, .true., rp_vars(c),       ncvars(c) )
          enddo

          if (mynode == 0) then
            write(*,'(9x,A,1x,F11.4,2x,A,I7,1x,A,I4,A,I4,1x,A,I3)')
     &       'wrote child grid, tdays =', tdays,
     &       'step =', iic-1, 'rec =',    ncvars(c)%output_rec, '/',
     &       ncvars(c)%total_output_rec
          endif

          ierr=nf90_close(ncvars(c)%ncid_out)
          if (ncvars(c)%output_rec == recs_per_file) then                              ! reset output_rec to create new file
            ncvars(c)%output_rec = 0
          endif

        end do ! <- c

!        write(*,'(9x,A,A12,5x,A,G14.6,1x,A,I4)')   ! 4x is spaces, A is string. / at end adds blank line
!     &  'set_frc :: ',         var_name,         ! 5x is 5 spaces, A is string
!     &  'input time (days) =', var_time,         ! G is number for time
!     &  'rec =',               var_irec MYID     ! I is integer

      endif ! <- mod(iic-1,rec_rate)

      end subroutine set_child_bry  !]

! ----------------------------------------------------------------------
      subroutine read_child_bry(ncv, rp, up, vp )  ![
      ! Assign 4 parent GP from which to interpolate 1 child BP. Done only once.

      implicit none

      ! inputs
      type(nc_vars)        :: ncv         ! netcdf variables for child
      type(grid_type_vars) :: rp, up, vp  ! rho, u, v-point variable types

      call assign_child_bry_to_parent(rbry_name(1), rp%i_south, rp%j_south, rp%tot(1),
     &                                rp%c_str(1),  rp%glob_cs, rp_var,     ncv      )
      call assign_child_bry_to_parent(rbry_name(2), rp%i_west,  rp%j_west,  rp%tot(2),
     &                                rp%c_str(2) , rp%glob_cw, rp_var,     ncv      )
      call assign_child_bry_to_parent(rbry_name(3), rp%i_north, rp%j_north, rp%tot(3),
     &                                rp%c_str(3),  rp%glob_cn, rp_var,     ncv      )
      call assign_child_bry_to_parent(rbry_name(4), rp%i_east,  rp%j_east,  rp%tot(4),
     &                                rp%c_str(4) , rp%glob_ce, rp_var,     ncv      )

      call assign_child_bry_to_parent(ubry_name(1), up%i_south, up%j_south, up%tot(1),
     &                                up%c_str(1),  up%glob_cs, up_var,     ncv      )
      call assign_child_bry_to_parent(ubry_name(2), up%i_west,  up%j_west,  up%tot(2),
     &                                up%c_str(2) , up%glob_cw, up_var,     ncv      )
      call assign_child_bry_to_parent(ubry_name(3), up%i_north, up%j_north, up%tot(3),
     &                                up%c_str(3),  up%glob_cn, up_var,     ncv      )
      call assign_child_bry_to_parent(ubry_name(4), up%i_east,  up%j_east,  up%tot(4),
     &                                up%c_str(4) , up%glob_ce, up_var,     ncv      )

      call assign_child_bry_to_parent(vbry_name(1), vp%i_south, vp%j_south, vp%tot(1),
     &                                vp%c_str(1),  vp%glob_cs, vp_var,     ncv      )
      call assign_child_bry_to_parent(vbry_name(2), vp%i_west,  vp%j_west,  vp%tot(2),
     &                                vp%c_str(2) , vp%glob_cw, vp_var,     ncv      )
      call assign_child_bry_to_parent(vbry_name(3), vp%i_north, vp%j_north, vp%tot(3),
     &                                vp%c_str(3),  vp%glob_cn, vp_var,     ncv      )
      call assign_child_bry_to_parent(vbry_name(4), vp%i_east,  vp%j_east,  vp%tot(4),
     &                                vp%c_str(4) , vp%glob_ce, vp_var,     ncv      )

      end subroutine read_child_bry  !]

! ----------------------------------------------------------------------
      subroutine assign_child_bry_to_parent( side_name,  i_bry, j_bry,                 ![
     &                                       side_total, c_str, dim_size, vtype, ncv )
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
      type(nc_vars)                                  :: ncv           ! netcdf variables for child

      ! local
      integer               :: varid, ij, l, counter
      real                  :: sx, sy                                 ! shift in x & y for different var type
      character(len=30)     :: vname
      integer, dimension(2) :: dimids

      real, dimension(:), allocatable :: i_in_bry, j_in_bry           ! store read in child boundary

      if (ncv%ncid==-1) then                                          ! file not yet opened.
        ierr=nf90_open(ncv%c_bry_file, nf90_nowrite, ncv%ncid)        ! var will later be in grdname file
        if(ierr/=0) call handle_ierr(ierr,
     &              'Error opening child bry file:', ncv%c_bry_file)
      endif


      vname='i_' / / side_name
      ierr=nf90_inq_varid(ncv%ncid, vname, varid)                         ! Check if contains bry
      if (ierr==nf90_noerr) then
        ierr=nf90_inquire_variable( ncv%ncid, varid,     dimids = dimids)
        ierr=nf90_inquire_dimension(ncv%ncid, dimids(1), len=dim_size)    ! bry vars only have 1 dim.

        allocate(i_in_bry(dim_size))                                      ! allocate temporary array for all bry data.
        ierr=nf90_get_var(ncv%ncid, varid, i_in_bry)                      ! read in data
        if(ierr/=0) print *, 'bad read of bry', mynode

        allocate(j_in_bry(dim_size))
        vname='j_' / / side_name
        ierr=nf90_inq_varid(ncv%ncid, vname, varid)
        ierr=nf90_get_var(  ncv%ncid, varid, j_in_bry)
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

        endif

      else
        if(mynode==0)
     &    write (*,'(/1x,3A)') 'No child boundary in parent domain ',
     &                         'for: ', side_name
        dim_size=0
      endif

      end subroutine assign_child_bry_to_parent  !]

! ----------------------------------------------------------------------
      subroutine calc_child_interpolation_coeffs( rp, up, vp, ncv )  ![

      ! Get interpolation coefficients from fraction of i/j positons.
      ! Calculate interpolation coefficients once only.

      implicit none

      ! inputs
      type(grid_type_vars) :: rp, up, vp  ! rho, u, v-point variable type
      type(nc_vars)        :: ncv

      ! rho cffs:
      if (allocated(rp%i_south)) then
        call calc_bry_coeffs(rp%tot(1),   rp%i_south,  rp%j_south, rp%s_cffs,
     &                       rp%is_sw,    rp%js_sw,    rp_var,     rp_var   )
      endif
      if (allocated(rp%i_west))  then
        call calc_bry_coeffs(rp%tot(2),   rp%i_west,   rp%j_west,  rp%w_cffs,
     &                       rp%iw_sw,    rp%jw_sw,    rp_var,     rp_var   )
      endif
      if (allocated(rp%i_north)) then
        call calc_bry_coeffs(rp%tot(3),   rp%i_north,  rp%j_north, rp%n_cffs,
     &                       rp%in_sw,    rp%jn_sw,    rp_var,     rp_var   )
      endif
      if (allocated(rp%i_east))  then
        call calc_bry_coeffs(rp%tot(4),   rp%i_east,   rp%j_east,  rp%e_cffs,
     &                       rp%ie_sw,    rp%je_sw,    rp_var,     rp_var   )
      endif

      ! u cffs:
      if (allocated(up%i_south)) then
        call calc_bry_coeffs(up%tot(1),   up%i_south,  up%j_south, up%s_cffs,            ! cffs to interp parent u to child u
     &                       up%is_sw,    up%js_sw,    up_var,     up_var   )
        call calc_bry_coeffs(up%tot(1),   up%i_south,  up%j_south, up%s_d2d_cffs,        ! cffs to interp parent v to child u
     &                       up%isd2d_sw, up%jsd2d_sw, vp_var,     up_var       )
        call calc_bry_coeffs(up%tot(1),   up%i_south,  up%j_south, up%s_r2d_cffs,        ! cffs to interp parent r to child u (angle)
     &                       up%isr2d_sw, up%jsr2d_sw, rp_var,     up_var       )

        call calc_parent_vs_child_grid_angles(         up%tot(1),  up%c_str(1),          ! don't need south_r2u after this, should deallocate
     &        up%s_r2d_cffs, up%isr2d_sw, up%jsr2d_sw, up%s_angcos,up%s_angsin, up_var, 1, ncv ) ! 1=S
      endif
      if (allocated(up%i_west)) then
        call calc_bry_coeffs(up%tot(2),   up%i_west,   up%j_west,  up%w_cffs,
     &                       up%iw_sw,    up%jw_sw,    up_var,     up_var   )
        call calc_bry_coeffs(up%tot(2),   up%i_west,   up%j_west,  up%w_d2d_cffs,
     &                       up%iwd2d_sw, up%jwd2d_sw, vp_var,     up_var       )
        call calc_bry_coeffs(up%tot(2),   up%i_west,   up%j_west,  up%w_r2d_cffs,
     &                       up%iwr2d_sw, up%jwr2d_sw, rp_var,     up_var       )
        call calc_parent_vs_child_grid_angles(         up%tot(2),  up%c_str(2),
     &        up%w_r2d_cffs, up%iwr2d_sw, up%jwr2d_sw, up%w_angcos,up%w_angsin, up_var, 2, ncv ) ! 2=W
      endif
      if (allocated(up%i_north)) then
        call calc_bry_coeffs(up%tot(3),   up%i_north,  up%j_north, up%n_cffs,
     &                       up%in_sw,    up%jn_sw,    up_var,     up_var   )
        call calc_bry_coeffs(up%tot(3),   up%i_north,  up%j_north, up%n_d2d_cffs,
     &                       up%ind2d_sw, up%jnd2d_sw, vp_var,     up_var       )
        call calc_bry_coeffs(up%tot(3),   up%i_north,  up%j_north, up%n_r2d_cffs,
     &                       up%inr2d_sw, up%jnr2d_sw, rp_var,     up_var       )
        call calc_parent_vs_child_grid_angles(         up%tot(3),  up%c_str(3),
     &        up%n_r2d_cffs, up%inr2d_sw, up%jnr2d_sw, up%n_angcos,up%n_angsin, up_var, 3, ncv )
      endif
      if (allocated(up%i_east)) then
        call calc_bry_coeffs(up%tot(4),   up%i_east,   up%j_east,  up%e_cffs,
     &                       up%ie_sw,    up%je_sw,    up_var,     up_var   )
        call calc_bry_coeffs(up%tot(4),   up%i_east,   up%j_east,  up%e_d2d_cffs,
     &                       up%ied2d_sw, up%jed2d_sw, vp_var,     up_var       )
        call calc_bry_coeffs(up%tot(4),   up%i_east,   up%j_east,  up%e_r2d_cffs,
     &                       up%ier2d_sw, up%jer2d_sw, rp_var,     up_var       )
        call calc_parent_vs_child_grid_angles(         up%tot(4),  up%c_str(4),
     &        up%e_r2d_cffs, up%ier2d_sw, up%jer2d_sw, up%e_angcos,up%e_angsin, up_var, 4, ncv )
      endif

      ! v cffs:
      if (allocated(vp%i_south)) then
        call calc_bry_coeffs(vp%tot(1),   vp%i_south,  vp%j_south, vp%s_cffs,            ! cffs to interp parent u to child u
     &                       vp%is_sw,    vp%js_sw,    vp_var,     vp_var   )
        call calc_bry_coeffs(vp%tot(1),   vp%i_south,  vp%j_south, vp%s_d2d_cffs,        ! cffs to interp parent v to child u
     &                       vp%isd2d_sw, vp%jsd2d_sw, up_var,     vp_var       )
        call calc_bry_coeffs(vp%tot(1),   vp%i_south,  vp%j_south, vp%s_r2d_cffs,        ! cffs to interp parent r to child u (angle)
     &                       vp%isr2d_sw, vp%jsr2d_sw, rp_var,     vp_var       )

        call calc_parent_vs_child_grid_angles(         vp%tot(1),  vp%c_str(1),          ! don't need south_r2u after this, should deallocate
     &        vp%s_r2d_cffs, vp%isr2d_sw, vp%jsr2d_sw, vp%s_angcos,vp%s_angsin, vp_var, 1, ncv ) ! 1=S
      endif
      if (allocated(vp%i_west)) then
        call calc_bry_coeffs(vp%tot(2),   vp%i_west,   vp%j_west,  vp%w_cffs,
     &                       vp%iw_sw,    vp%jw_sw,    vp_var,     vp_var   )
        call calc_bry_coeffs(vp%tot(2),   vp%i_west,   vp%j_west,  vp%w_d2d_cffs,
     &                       vp%iwd2d_sw, vp%jwd2d_sw, up_var,     vp_var       )
        call calc_bry_coeffs(vp%tot(2),   vp%i_west,   vp%j_west,  vp%w_r2d_cffs,
     &                       vp%iwr2d_sw, vp%jwr2d_sw, rp_var,     vp_var       )
        call calc_parent_vs_child_grid_angles(         vp%tot(2),  vp%c_str(2),
     &        vp%w_r2d_cffs, vp%iwr2d_sw, vp%jwr2d_sw, vp%w_angcos,vp%w_angsin, vp_var, 2, ncv ) ! 2=W
      endif
      if (allocated(vp%i_north)) then
        call calc_bry_coeffs(vp%tot(3),   vp%i_north,  vp%j_north, vp%n_cffs,
     &                       vp%in_sw,    vp%jn_sw,    vp_var,     vp_var   )
        call calc_bry_coeffs(vp%tot(3),   vp%i_north,  vp%j_north, vp%n_d2d_cffs,
     &                       vp%ind2d_sw, vp%jnd2d_sw, up_var,     vp_var       )
        call calc_bry_coeffs(vp%tot(3),   vp%i_north,  vp%j_north, vp%n_r2d_cffs,
     &                       vp%inr2d_sw, vp%jnr2d_sw, rp_var,     vp_var       )
        call calc_parent_vs_child_grid_angles(         vp%tot(3),  vp%c_str(3),
     &        vp%n_r2d_cffs, vp%inr2d_sw, vp%jnr2d_sw, vp%n_angcos,vp%n_angsin, vp_var, 3, ncv )
      endif
      if (allocated(vp%i_east)) then
        call calc_bry_coeffs(vp%tot(4),   vp%i_east,   vp%j_east,  vp%e_cffs,
     &                       vp%ie_sw,    vp%je_sw,    vp_var,     vp_var   )
        call calc_bry_coeffs(vp%tot(4),   vp%i_east,   vp%j_east,  vp%e_d2d_cffs,
     &                       vp%ied2d_sw, vp%jed2d_sw, up_var,     vp_var       )
        call calc_bry_coeffs(vp%tot(4),   vp%i_east,   vp%j_east,  vp%e_r2d_cffs,
     &                       vp%ier2d_sw, vp%jer2d_sw, rp_var,     vp_var       )
        call calc_parent_vs_child_grid_angles(         vp%tot(4),  vp%c_str(4),
     &        vp%e_r2d_cffs, vp%ier2d_sw, vp%jer2d_sw, vp%e_angcos,vp%e_angsin, vp_var, 4, ncv )
      endif

      rt2=sqrt(2.0)  ! set interpolation coefficient once

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
      ! We take the the child ij values (e.g. up%s_cffs), which are in terms of
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
      subroutine calc_parent_vs_child_grid_angles( bry_tot, c_str,   cffs,    i_sw, j_sw,    ![
     &                                             anglcos, anglsin, c_vtype, bry,  ncv )

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
      type(nc_vars)                                     :: ncv

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

      ierr=nf90_inq_varid(ncv%ncid, vname, varid)                          ! Check if contains bry
      if (ierr==nf90_noerr) then
        ierr=nf90_get_var(ncv%ncid, varid, c_angle, start=(/c_str/), count=(/bry_tot/) ) ! only read child points in subdomain
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
      subroutine interpolate_child_bry( var, nmax, rp )  ![ ! CHANGE NAME TO INTERP RHO

      ! Interpolate all boundaries if applicable

      implicit none

      ! inputs
      type(grid_type_vars) :: rp  ! rho-point variable type

      integer, intent(in) :: nmax
      real   , intent(in) :: var(GLOBAL_2D_ARRAY,nmax)  ! variable to interpolate.
                                                        ! wasn't able to use dimension(:,:,:), but vars start from -1 index so maybe this is safer.

      if (allocated(rp%i_south))
     &  call interpolate_one_child_bry(rp%tot(1), rp%s_cffs, rp%is_sw, rp%js_sw,
     &                                 nmax,      var,       rp%south, rmask   )

      if (allocated(rp%i_west))
     &  call interpolate_one_child_bry(rp%tot(2), rp%w_cffs, rp%iw_sw, rp%jw_sw,
     &                                 nmax,      var,       rp%west,  rmask   )

      if (allocated(rp%i_north))
     &  call interpolate_one_child_bry(rp%tot(3), rp%n_cffs, rp%in_sw, rp%jn_sw,
     &                                 nmax,      var,       rp%north, rmask   )

      if (allocated(rp%i_east))
     &  call interpolate_one_child_bry(rp%tot(4), rp%e_cffs, rp%ie_sw, rp%je_sw,
     &                                 nmax,      var,       rp%east,  rmask   )


      end subroutine interpolate_child_bry  !]

! ----------------------------------------------------------------------
      subroutine interpolate_u_child_bry( uvar, vvar, nmax, up )  ![
      ! Interpolate all u boundaries if applicable

      implicit none

      ! inputs
      integer,  intent(in) :: nmax
      real   ,  intent(in) :: uvar(GLOBAL_2D_ARRAY,nmax)  ! variable to interpolate.
      real   ,  intent(in) :: vvar(GLOBAL_2D_ARRAY,nmax)  ! variable to interpolate.
      type(grid_type_vars) :: up                          ! u-point variable type


      if (allocated(up%i_south)) then
        call interpolate_one_child_bry(up%tot(1), up%s_cffs,     up%is_sw,    up%js_sw,
     &                                 nmax,      uvar,          up%s_uprnt,  umask   )
        call interpolate_one_child_bry(up%tot(1), up%s_d2d_cffs, up%isd2d_sw, up%jsd2d_sw,
     &                                 nmax,      vvar,          up%s_vprnt,  vmask      )

        call rotate_one_u_bry(         up%tot(1), nmax, up%south, up%s_uprnt, up%s_vprnt,
     &                                 up%s_angcos, up%s_angsin                         )
      endif

      if (allocated(up%i_west)) then
        call interpolate_one_child_bry(up%tot(2), up%w_cffs,     up%iw_sw,    up%jw_sw,
     &                                 nmax,      uvar,          up%w_uprnt,  umask   )
        call interpolate_one_child_bry(up%tot(2), up%w_d2d_cffs, up%iwd2d_sw, up%jwd2d_sw,
     &                                 nmax,      vvar,          up%w_vprnt,  vmask      )

        call rotate_one_u_bry(         up%tot(2), nmax, up%west, up%w_uprnt,  up%w_vprnt,
     &                                 up%w_angcos, up%w_angsin                         )
      endif

      if (allocated(up%i_north)) then
        call interpolate_one_child_bry(up%tot(3), up%n_cffs,     up%in_sw,    up%jn_sw,
     &                                 nmax,      uvar,          up%n_uprnt,  umask   )
        call interpolate_one_child_bry(up%tot(3), up%n_d2d_cffs, up%ind2d_sw, up%jnd2d_sw,
     &                                 nmax,      vvar,          up%n_vprnt,  vmask      )

        call rotate_one_u_bry(         up%tot(3), nmax, up%north, up%n_uprnt, up%n_vprnt,
     &                                 up%n_angcos, up%n_angsin                          )
      endif

      if (allocated(up%i_east)) then
        call interpolate_one_child_bry(up%tot(4), up%e_cffs,     up%ie_sw,    up%je_sw,
     &                                 nmax,      uvar,          up%e_uprnt,  umask   )
        call interpolate_one_child_bry(up%tot(4), up%e_d2d_cffs, up%ied2d_sw, up%jed2d_sw,
     &                                 nmax,      vvar,          up%e_vprnt,  vmask      )

        call rotate_one_u_bry(         up%tot(4), nmax, up%east, up%e_uprnt,  up%e_vprnt,
     &                                 up%e_angcos, up%e_angsin                         )
      endif

      end subroutine interpolate_u_child_bry  !]

! ----------------------------------------------------------------------
      subroutine interpolate_v_child_bry( uvar, vvar, nmax, vp )  ![
      ! Interpolate all v boundaries if applicable
      ! Shouldn't really need this, rather have uv function with a type as an arguement!!!

      implicit none

      integer,  intent(in) :: nmax
      real   ,  intent(in) :: uvar(GLOBAL_2D_ARRAY,nmax)  ! variable to interpolate.
      real   ,  intent(in) :: vvar(GLOBAL_2D_ARRAY,nmax)  ! variable to interpolate.
      type(grid_type_vars) :: vp                          ! v-point variable type


      if (allocated(vp%i_south)) then
        call interpolate_one_child_bry(vp%tot(1), vp%s_cffs,     vp%is_sw,    vp%js_sw,
     &                                 nmax,      vvar,          vp%s_vprnt,  vmask   )
        call interpolate_one_child_bry(vp%tot(1), vp%s_d2d_cffs, vp%isd2d_sw, vp%jsd2d_sw,
     &                                 nmax,      uvar,          vp%s_uprnt,  umask      )

        call rotate_one_v_bry(         vp%tot(1), nmax, vp%south, vp%s_uprnt, vp%s_vprnt,
     &                                 vp%s_angcos, vp%s_angsin                         )
      endif

      if (allocated(vp%i_west)) then
        call interpolate_one_child_bry(vp%tot(2), vp%w_cffs,     vp%iw_sw,    vp%jw_sw,
     &                                 nmax,      vvar,          vp%w_vprnt,  vmask   )
        call interpolate_one_child_bry(vp%tot(2), vp%w_d2d_cffs, vp%iwd2d_sw, vp%jwd2d_sw,
     &                                 nmax,      uvar,          vp%w_uprnt,  umask      )

        call rotate_one_v_bry(         vp%tot(2), nmax, vp%west, vp%w_uprnt,  vp%w_vprnt,
     &                                 vp%w_angcos, vp%w_angsin                         )
      endif

      if (allocated(vp%i_north)) then
        call interpolate_one_child_bry(vp%tot(3), vp%n_cffs,     vp%in_sw,    vp%jn_sw,
     &                                 nmax,      vvar,          vp%n_vprnt,  vmask   )
        call interpolate_one_child_bry(vp%tot(3), vp%n_d2d_cffs, vp%ind2d_sw, vp%jnd2d_sw,
     &                                 nmax,      uvar,          vp%n_uprnt,  umask      )

        call rotate_one_v_bry(         vp%tot(3), nmax, vp%north, vp%n_uprnt, vp%n_vprnt,
     &                                 vp%n_angcos, vp%n_angsin                          )
      endif

      if (allocated(vp%i_east)) then
        call interpolate_one_child_bry(vp%tot(4), vp%e_cffs,     vp%ie_sw,    vp%je_sw,
     &                                 nmax,      vvar,          vp%e_vprnt,  vmask   )
        call interpolate_one_child_bry(vp%tot(4), vp%e_d2d_cffs, vp%ied2d_sw, vp%jed2d_sw,
     &                                 nmax,      uvar,          vp%e_uprnt,  umask      )

        call rotate_one_v_bry(         vp%tot(4), nmax, vp%east, vp%e_uprnt,  vp%e_vprnt,
     &                                 vp%e_angcos, vp%e_angsin                         )
      endif

      end subroutine interpolate_v_child_bry  !]

! ----------------------------------------------------------------------
      subroutine interpolate_one_child_bry(bry_tot, cffs, i_sw,   j_sw,  ![
     &                                     nmax,    var,  varbry, mask )

      ! Interpolate 4 parent points to 1 child bry point (r2r, u2u, or v2v)
      ! Interpolation coefficients already calculated in calc_child_interpolation_coeffs

      implicit none

      ! inputs
      integer,                              intent(in)  :: bry_tot                    ! total child bry points in parent
      real,    dimension(:,:),              intent(in)  :: cffs                       ! 4 parent interpolants for each child point
      integer, dimension(:),                intent(in)  :: i_sw, j_sw                 ! i and j of sw parent interpolant
      integer,                              intent(in)  :: nmax                       ! vertical layers of var (1 if 2D)
      real,                                 intent(in)  :: var(GLOBAL_2D_ARRAY,nmax)  ! parent var to interpolate from
      real,    dimension(:,:), allocatable, intent(out) :: varbry                     ! bry to interpolated 2D/3D var
      real,                                 intent(in)  :: mask(GLOBAL_2D_ARRAY)      ! grid type r/u/v mask

      ! local
      integer                 :: ij, k
      real, dimension(nmax)   :: var_sw, var_se, var_nw, var_ne  ! temporary copies of columns of data
      real(kind=8), parameter :: spv_set=1.D+33

      if (.not. allocated(varbry)) allocate( varbry(bry_tot, N)) ! allocate once for 3D var and reuse for all 2D/3D vars
                                                                 ! maybe move this to initialization to avoid if statement repetition

      do ij=1,bry_tot                                            ! Normally would have k as outer loop for cache coherence but can avoid many if statements by swapping order here:

        if(i_sw(ij) >= 0) then                                   ! child point not masked. >=0 as mask value is -1. i_sw can be 0 if point is between nodes.

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

            varbry(ij,k) = cffs(1,ij) * var_sw(k)   ! interp from sw parent point.      Note even though 'k' index doesn't exist for 2D var,
     &                   + cffs(2,ij) * var_se(k)   ! interp from se parent point.      fortran won't complain. This is done in ncdf_write_mod() too.
     &                   + cffs(3,ij) * var_nw(k)   ! interp from nw parent point.
     &                   + cffs(4,ij) * var_ne(k)   ! interp from ne parent point.

          enddo

        else  ! child point masked

          varbry(ij,:) = spv_set

        endif ! <- i_sw(ij)>0

      enddo   ! <- ij

      end subroutine interpolate_one_child_bry  !]

! ----------------------------------------------------------------------
      subroutine rotate_one_u_bry( bry_tot, nmax, ubry, uparent, vparent,  ![
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

!        if (mynode==5) then
!          if( up%js_sw(ij) > 45 ) then  ! debug
!            print *, 'up%js_sw(ij)',up%js_sw(ij), mynode
!            print *, 'uparent(ij,1)',uparent(ij,1), mynode
!            print *, 'vparent(ij,1)',vparent(ij,1), mynode
!            print *, 'angcos(ij)',angcos(ij),'angsin(ij)',angsin(ij)
!          endif
!        endif

      enddo

      end subroutine rotate_one_u_bry  !]

! ----------------------------------------------------------------------
      subroutine rotate_one_v_bry( bry_tot, nmax, vbry, uparent, vparent,  ![
     &                                  angcos,  angsin )

      implicit none

      ! inputs
      integer,                              intent(in)  :: bry_tot     ! total child bry points in parent
      integer,                              intent(in)  :: nmax        ! vertical layers of var (1 if 2D)
      real,    dimension(:,:), allocatable, intent(out) :: vbry        ! bry to interpolated 2D/3D var
      real,    dimension(:,:),              intent(in)  :: uparent     ! u-value of parent grid at child u-point
      real,    dimension(:,:),              intent(in)  :: vparent     ! v-value of parent grid at child u-point
      real,    dimension(:),                intent(in)  :: angcos      ! cos of angle difference between parent-child grids
      real,    dimension(:),                intent(in)  :: angsin      ! sin of ...

      ! local
      integer :: ij, k

      if (.not. allocated(vbry)) allocate( vbry(bry_tot, N))   ! allocate once for 3D var and reuse for all 2D/3D vars

      do ij=1,bry_tot
        do k=1,nmax
          vbry(ij,k) = uparent(ij,k) * angsin(ij) + vparent(ij,k) * angcos(ij)
        enddo
      enddo

      end subroutine rotate_one_v_bry  !]

! ----------------------------------------------------------------------
      subroutine create_child_output_file( rp, up, vp, ncv )  ![

      ![For now decided each file still creates file even if not child contained in node.
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
      !]

      implicit none

      ! inputs
      type(grid_type_vars)  :: rp, up, vp                                   ! rho,u,v-point variable types
      type(nc_vars)         :: ncv                                          ! netcdf variables for child

      ! local
      integer               :: prev_fill_mode
      integer, dimension(5) :: r_dimid_bry,  u_dimid_bry,  v_dimid_bry      ! dimension ID for boundaries
      integer               :: itrc                                         ! tracer index

      if(ncv%total_output_rec == 0) then                                    ! Add root output name to
        ncv%file_name = trim(output_root_name) / / trim(ncv%file_name)      ! file name (only done once)
      endif


      call ncdf_create_blank_file(ncv%file_name, ncv%ncid_out,  prev_fill_mode,
     &                            recs_per_file, ncv%total_output_rec )

      ! DEFINE DIMENSIONS & VARIABLES:

      ierr=nf90_put_att( ncv%ncid_out, NF90_GLOBAL, 'global_south_size', rp%glob_cs )
      ierr=nf90_put_att( ncv%ncid_out, NF90_GLOBAL, 'global_west_size' , rp%glob_cw )
      ierr=nf90_put_att( ncv%ncid_out, NF90_GLOBAL, 'global_north_size', rp%glob_cn )
      ierr=nf90_put_att( ncv%ncid_out, NF90_GLOBAL, 'global_east_size' , rp%glob_ce )
      ierr=nf90_def_dim( ncv%ncid_out, 'rec_time', 0, ncv%dimid_t)

      ! rho dims:
      ierr=nf90_def_dim( ncv%ncid_out, 'xi_rho_south', max( rp%tot(1), 1 ), r_dimid_bry(1) )  ! Set size as 1 if size==0 to avoid unlimited dimension (could be confused with time in ncjoin)
      ierr=nf90_def_dim( ncv%ncid_out, 'eta_rho_west', max( rp%tot(2), 1 ), r_dimid_bry(2) )
      ierr=nf90_def_dim( ncv%ncid_out, 'xi_rho_north', max( rp%tot(3), 1 ), r_dimid_bry(3) )
      ierr=nf90_def_dim( ncv%ncid_out, 'eta_rho_east', max( rp%tot(4), 1 ), r_dimid_bry(4) )
      ierr=nf90_def_dim( ncv%ncid_out, 's_rho'       , N                  , r_dimid_bry(5) )

      ! u dims: (see preamble for eta_urho_west & eta_urho_east info)
      ierr=nf90_def_dim( ncv%ncid_out, 'xi_u_south',   max( up%tot(1), 1 ), u_dimid_bry(1) )
      ierr=nf90_def_dim( ncv%ncid_out, 'eta_urho_west',max( up%tot(2), 1 ), u_dimid_bry(2) )
      ierr=nf90_def_dim( ncv%ncid_out, 'xi_u_north',   max( up%tot(3), 1 ), u_dimid_bry(3) )
      ierr=nf90_def_dim( ncv%ncid_out, 'eta_urho_east',max( up%tot(4), 1 ), u_dimid_bry(4) )
      u_dimid_bry(5) = r_dimid_bry(5)  ! k dimension of u-point is same as rho

      ! v dims:
      ierr=nf90_def_dim( ncv%ncid_out, 'xi_vrho_south',max( vp%tot(1), 1 ), v_dimid_bry(1) )
      ierr=nf90_def_dim( ncv%ncid_out, 'eta_v_west',   max( vp%tot(2), 1 ), v_dimid_bry(2) )
      ierr=nf90_def_dim( ncv%ncid_out, 'xi_vrho_north',max( vp%tot(3), 1 ), v_dimid_bry(3) )
      ierr=nf90_def_dim( ncv%ncid_out, 'eta_v_east',   max( vp%tot(4), 1 ), v_dimid_bry(4) )
      v_dimid_bry(5) = r_dimid_bry(5)  ! k dimension of v-point is same as rho

      ! 2D vars:
      call define_child_bry( indxZ,  rp%tot, rp%c_str, r_dimid_bry, 1, .false., ncv )
      call define_child_bry( indxUb, up%tot, up%c_str, u_dimid_bry, 1, .false., ncv )
      call define_child_bry( indxVb, vp%tot, vp%c_str, v_dimid_bry, 1, .false., ncv )

      ! 3D vars:
      call define_child_bry( indxU,  up%tot, up%c_str, u_dimid_bry, N, .false., ncv )
      call define_child_bry( indxV,  vp%tot, vp%c_str, v_dimid_bry, N, .false., ncv )

      ! tracers:
      do itrc=1,NT
        call define_child_bry( itrc, rp%tot, rp%c_str, r_dimid_bry, N, .true. , ncv )
      enddo


      ierr=nf90_enddef(ncv%ncid_out)  ! end definition to write changes to disk
      if (mynode == 0) then
        write(*,'(/1x,2A/)')
     &     'bry_extract_child :: created new netCDF file ', ncv%file_name
      endif

      end subroutine create_child_output_file  !]

! ----------------------------------------------------------------------
      subroutine wrt_child_bry ( indx_v, nmax, trc_yn, ap, ncv )  ![

      ! Each sub-domain writes out its portion of child boundary if applicable.
      ! all procs, make child bc netcdf file, only relevant procs define bc variables (east, north, etc)

      implicit none

      ! inputs
      integer,              intent(in) :: indx_v    ! index of primary var or itrc
      integer,              intent(in) :: nmax      ! vertical levels
      logical,              intent(in) :: trc_yn    ! is variable a tracer y/n
      type(grid_type_vars), intent(in) :: ap        ! any variable type
      type(nc_vars)                    :: ncv

      ! local
      integer               :: varid, ndims
      integer, dimension(3) :: start_out, count_out, count_dummy=1
      character(len=40)     :: bname, rname


      if (trc_yn) then                              ! handle primary variable or tracers name
        rname=t_vname(indx_v)
      else
        rname=vname(1,indx_v)
      endif

      if (nmax == 1) then                           ! handle 2D or 3D variable dimids
        ndims=2
        start_out(1:ndims) = (/1, ncv%output_rec/)
        count_out(1:ndims) = (/1, 1             /)
      else
        ndims=3
        start_out(1:ndims) = (/1, 1   , ncv%output_rec/)
        count_out(1:ndims) = (/1, nmax, 1             /)
      endif


      bname = trim(rname) / /'_south'
      ierr  = nf90_inq_varid( ncv%ncid_out, bname, varid )
      count_out(1) = ap%tot(1)
      if (allocated(ap%south)) then
        ierr=nf90_put_var( ncv%ncid_out, varid, ap%south(:,:),
     &                     start=start_out(1:ndims), count=count_out(1:ndims) )
        if(ierr/=0) write(*,*) 'wrts_bry ierr= ',nf90_strerror(ierr), mynode
      else
        ierr=nf90_put_var( ncv%ncid_out, varid, (/1.0/),                               ! put dummy value so ncjoin doesn't complain
     &                     start=start_out(1:ndims), count=count_dummy(1:ndims) )
      endif

      bname = trim(rname) / /'_west'
      ierr  = nf90_inq_varid( ncv%ncid_out, bname, varid )
      count_out(1) = ap%tot(2)
      if (allocated(ap%west)) then
        ierr=nf90_put_var( ncv%ncid_out, varid, ap%west(:,:) ,
     &                     start=start_out(1:ndims), count=count_out(1:ndims) )
        if(ierr/=0) write(*,*) 'wrtw_bry ierr= ',nf90_strerror(ierr), mynode
      else
        ierr=nf90_put_var( ncv%ncid_out, varid, (/1.0/),
     &                     start=start_out(1:ndims), count=count_dummy(1:ndims) )
      endif

      bname = trim(rname) / /'_north'
      ierr  = nf90_inq_varid( ncv%ncid_out, bname, varid )
      count_out(1) = ap%tot(3)
      if (allocated(ap%north)) then
        ierr=nf90_put_var( ncv%ncid_out, varid, ap%north(:,:),
     &                     start=start_out(1:ndims), count=count_out(1:ndims) )
        if(ierr/=0) write(*,*) 'wrtn_bry ierr= ',nf90_strerror(ierr), mynode
      else
        ierr=nf90_put_var( ncv%ncid_out, varid, (/1.0/),
     &                     start=start_out(1:ndims), count=count_dummy(1:ndims) )
      endif

      bname = trim(rname) / /'_east'
      ierr  = nf90_inq_varid( ncv%ncid_out, bname, varid )
      count_out(1) = ap%tot(4)
      if (allocated(ap%east )) then
        ierr=nf90_put_var( ncv%ncid_out, varid, ap%east(:,:) ,
     &                     start=start_out(1:ndims), count=count_out(1:ndims) )
        if(ierr/=0) write(*,*) 'wrte_bry ierr= ',nf90_strerror(ierr), mynode
      else
        ierr=nf90_put_var( ncv%ncid_out, varid, (/1.0/),
     &                     start=start_out(1:ndims), count=count_dummy(1:ndims) )
      endif

      end subroutine wrt_child_bry  !]

! ----------------------------------------------------------------------
      subroutine define_child_bry( indx_v, sides_tot, sides_c_str,  ![
     &                             dimid_bry, nmax, trc_yn, ncv  )

      ! Define child bry variable in output netCDF file.

      implicit none

      ! input/outputs
      integer,               intent(in) :: indx_v
      integer, dimension(4), intent(in) :: sides_tot
      integer, dimension(4), intent(in) :: sides_c_str
      integer, dimension(5), intent(in) :: dimid_bry     ! dimension ID for boundaries
      integer,               intent(in) :: nmax
      logical,               intent(in) :: trc_yn
      type(nc_vars)                     :: ncv           ! netcdf variables for child

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
          vdimids(1:ndims) = (/dimid_bry(bry), ncv%dimid_t/)
        else
          ndims=3
          vdimids(1:ndims) = (/dimid_bry(bry), dimid_bry(5), ncv%dimid_t/)
        endif


        bname=trim(rname) / /'_'/ /trim( bry_name(bry) )

        call nc_define_var( ncv%ncid_out,     bname, lname, uname, ! longname, units
     &                      vdimids(1:ndims), ierr )


        ierr=nf90_inq_varid(ncv%ncid_out, bname, varid)

        ierr=nf90_put_att(  ncv%ncid_out, varid, 'start',     sides_c_str( bry ))  ! add attribute with global child start value (i.e. south_c_str)
        ierr=nf90_put_att(  ncv%ncid_out, varid, 'count',     sides_tot(   bry ))  ! and its count (i.e. south_tot). This will make it easy for ncjoin.
        ierr=nf90_put_att(  ncv%ncid_out, varid, 'child_bry', 'True')                        ! 'child_bry' for ncjoin to identify variable. 'True' just a filler.

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

      ca = rt2*mask(i-1,j-1)                     ! corner point weightings BUG THIS SHOULD BE RMASK
      cc = rt2*mask(i+1,j-1)                     ! (calculate once outside of vertical loop)
      cf = rt2*mask(i-1,j+1)
      ch = rt2*mask(i+1,j+1)

      mb = mask(i  ,j-1)                         ! middle point weightings
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
