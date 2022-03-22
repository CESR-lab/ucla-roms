      module extract_data  
      ! Extract data at various positions and frequencies
      ! only for things that need interpolating

      ! We're retooling this to be more generally a way to output a
      ! variety of things at a range of locations at various frequencies

      ! Auxiliary tools required:
      ! - Tools-Roms/bry_child_input_gen.m   -> creates input file of child bry i/j values used here
      ! - Tools-Roms/join_child_bry.F        -> joins partial files created here

      !
      ! STEPS:
      !
      ! 0) CREATE INPUT FILE OF EXTRACTION OBJECTS
      !
      !    - Positions (in fractional i,j locations of the grid)
      !    - Name of the object ('child1_south','mooring2', etc, etc,...
      !    - Output frequency of each object
      !    - Variables to be output
      !    - Object name, positions, var type (rho/vel), variables, frequency
      !
      ! 1) Assign positions to subdomains
      !    - if ispresent(iobj) then see above
      !    - keep track of placement in global object arrays
      !    The ROMS parent simulation needs to know where the child boundary sits in its
      !    parent domain. To do this the child bry points are given i and j coords in
      !    terms of the parent grid. We created a matlab tool to do this found here:
      !    Tools-Roms/bry_child_input_gen.m

      ! 2) create a single file per subdomain containing all the object
      !    for which ispresent(iobj) is true
      !
      ! 3) Loop through all objects and all vars, and write when needed
      !    Add averaging capability at some point
      !
      !    vel points always need to interpolate both u, and v. in order
      !    to rotate the vector to a desired angle
      !
      !]

#include "cppdefs.opt"
#ifdef EXTRACT_DATA

      use dimensions
      use grid
      use nc_read_write
      use netcdf
      use mpi                                          ! debug temporary for mpi abort
      use tracers, only: t, t_vname, t_lname, t_units  ! need to get names of tracers
      use ocean2d
      use ocean3d
      use scalars

      implicit none
      private

#include "bry_extract_child.opt"

      ! scalar and vectors will always be in their own object
      ! for now, lets do all rho vars for scalar objects (zeta,T,S and other
      ! tracers), and all velocities for vector objects (u,v,ubar,vbar)
      ! next step is to make that selectable
      ! next step to add averaging

      character(len=60) :: extract_file = 'sample_edata.nc'
      integer :: nobj  ! number of extraction objects
 
      type extract_object  ! contains all information for a data_extraction object

        ! needed as input
        character(len=60) :: obj_name                   ! name of object
        logical                            :: scalar    ! scalar or vector
        integer                            :: period    ! output period
        integer                            :: np        ! local number of locations
        real,   dimension(:)  ,allocatable :: ipos,jpos ! fractional index locations
        real,   dimension(:)  ,allocatable :: ang       ! desired angle for velocities
        integer                            :: start_idx ! starting position in the obj array

        real                               :: otime=0   ! time since last output
        integer                            :: record=0  ! record number

        real,   dimension(:,:),pointer     :: vari      ! data

        ! these are only for scalars
        integer,dimension(:)  ,pointer     :: ip,jp     ! lower left hand point index
        real   ,dimension(:,:),pointer     :: coef      ! interpolation coefficients

        ! these are only for velocities
        real,   dimension(:)  ,allocatable :: cosa,sina ! for rotation of vectors
        integer,dimension(:)  ,pointer     :: ipu,jpu   ! only for vectors
        integer,dimension(:)  ,pointer     :: ipv,jpv   ! only for vectors
        real   ,dimension(:,:),pointer     :: cfu,cfv   ! only for vectors
        real   ,dimension(:,:),pointer     :: ui,vi     ! only for vectors 

      end type extract_object

      type(extract_object),dimension(:),allocatable :: obj

      !! examples of objects: 'grid1_south_r', mooring1_r'

      interface interpolate
        module procedure  interpolate_2D, interpolate_3D,
     &                  interpolate_vec_2D,interpolate_vec_3D
      end interface

      public do_extract_data

      contains
! ----------------------------------------------------------------------
      subroutine init_extract_data  ![
      ! Allocate space and compute interpolation coefficients
      ! for rho,u,and v variables.
      implicit none

      !local
      integer :: i,np

      call read_extraction_objects

      do i = 1,nobj
        np = obj(i)%np
        allocate(obj(i)%vari(np,nz))
        allocate(obj(i)%coef(np,4))
        allocate(obj(i)%ip(np))
        allocate(obj(i)%jp(np))

        ! from absolute index to rho-index
        obj(i)%ipos = obj(i)%ipos+0.5
        obj(i)%jpos = obj(i)%jpos+0.5
        call compute_coef(obj(i)%ipos,obj(i)%jpos,
     &              obj(i)%coef,obj(i)%ip,obj(i)%jp,rmask)

        if (.not.obj(i)%scalar) then
!          call compute_rotate(obj)  ! DevinD commented out to get to compile. JM must fix.

          allocate(obj(i)%ui(np,nz))
          allocate(obj(i)%cfu(np,4))
          allocate(obj(i)%ipu(np))
          allocate(obj(i)%jpu(np))

          ! from rho-index to u-index 
          obj(i)%ipos = obj(i)%ipos+0.5
          call compute_coef(obj(i)%ipos,obj(i)%jpos,
     &              obj(i)%cfu,obj(i)%ipu,obj(i)%jpu,umask)

          allocate(obj(i)%vi(obj(i)%np,nz))
          allocate(obj(i)%cfv(obj(i)%np,4))
          allocate(obj(i)%ipv(obj(i)%np))
          allocate(obj(i)%jpv(obj(i)%np))

          ! from u-index to v-index
          obj(i)%ipos = obj(i)%ipos-0.5
          obj(i)%jpos = obj(i)%jpos+0.5
          call compute_coef(obj(i)%ipos,obj(i)%jpos,
     &              obj(i)%cfv,obj(i)%ipv,obj(i)%jpv,vmask)

        endif
      enddo

      end subroutine init_extract_data !]
! ----------------------------------------------------------------------
      subroutine read_extraction_objects  ![
      ! Read all objects from file and determine local ranges 
      implicit none

      ! local
      integer               :: iobj,ncid
      integer,dimension(2)  :: dimids
      character(len=30) :: objname
      real,dimension(:,:),allocatable :: object
      integer               :: n1,n2,i0,i1
      integer req(2), status(MPI_STATUS_SIZE,2), ierr


      if (mynode==0) then
        print *, 'read objects: ',extract_file
      endif
      ierr = nf90_open(extract_file,nf90_nowrite,ncid)
      ierr = nf90_inquire(ncid,nVariables=nobj)

      allocate(obj(nobj))

      ! Read all objects from file. 
      do iobj = 1,nobj
        ierr = nf90_inquire_variable(ncid,iobj,name=objname,dimids=dimids)
        ierr = nf90_inquire_dimension(ncid,dimids(1),len=n1)
        ierr = nf90_inquire_dimension(ncid,dimids(2),len=n2)

        obj(iobj)%obj_name = objname
        if (n2==2) obj(iobj)%scalar = .true.

        allocate(object(n1,n2))
        call ncread(ncid,objname,object)

        call find_local_points(obj(iobj),object(:,1),object(:,2))

        if (n2==3) then
          allocate(obj(iobj)%ang(obj(iobj)%np))
          i0 = obj(iobj)%start_idx
          i1 = i0+obj(iobj)%np-1
          obj(iobj)%ang = object(i0:i1,3)
        endif

        deallocate(object)
      enddo
      ierr = nf90_close(ncid)

      end subroutine read_extraction_objects !]
! ----------------------------------------------------------------------
      subroutine find_local_points(obj,gobj_i,gobj_j) ![

      ! Find object index locations that are within the subdomain
      ! Assign start and lenght of the local points in the global array
      ! of the object
      ! Translate global index locations to local ones
      implicit none
      ! import/export
      type(extract_object),intent(inout) ::obj
      real,dimension(:)   ,intent(inout) ::gobj_i,gobj_j ! global indices

      ! local
      integer :: i,start_idx,end_idx

      gobj_i = gobj_i-iSW_corn + 1
      gobj_j = gobj_j-jSW_corn + 1

      ! Assume that local ranges of objects are contiguous
      start_idx = 0
      do i = 1,nobj
        if ( gobj_i(i)>=0.and.gobj_i(i)<nx .and.
     &       gobj_j(i)>=0.and.gobj_j(i)<ny ) then
          start_idx = i
        endif
      enddo
      if (start_idx>0) then
        do i = start_idx,nobj
          if (gobj_i(i)<0.or.gobj_i(i)>=nx.or.
     &        gobj_j(i)<0.or.gobj_j(i)>=ny ) then
            end_idx = i-1
          endif
        enddo
      else ! object not in local range
        obj%np = 0
      endif
      
      obj%start_idx = start_idx
      obj%np = end_idx-start_idx+1
      allocate(obj%ipos(obj%np))
      allocate(obj%jpos(obj%np))
      obj%ipos = gobj_i(start_idx:end_idx)
      obj%jpos = gobj_j(start_idx:end_idx)

      end subroutine find_local_points  !]
! ----------------------------------------------------------------------
      subroutine compute_coef(ipos,jpos,coef,ip,jp,mask)  ![
      ! compute interpolation coefficients
      implicit none

      ! inport/export
      real   ,dimension(:)  ,intent(in) :: ipos,jpos
      real   ,dimension(:,:),intent(out):: coef
      integer,dimension(:)  ,intent(out):: ip,jp
      real   ,dimension(-1:nx+1,-1:ny+1),intent(in) :: mask

      ! local
      integer :: i,np
      real :: cfx,cfy

      np = size(ip,1)
      do i = 1,np
        ip(i)  = floor(ipos(i))
        cfx    = ipos(i)-ip(i)
        jp(i)  = floor(jpos(i))
        cfy    = jpos(i)-jp(i)
        coef(i,1) = (1-cfx)*(1-cfy)*mask(ip(i)  ,jp(i)  ) 
        coef(i,2) = cfx    *(1-cfy)*mask(ip(i)+1,jp(i)  ) 
        coef(i,3) = (1-cfx)*   cfy *mask(ip(i)  ,jp(i)+1)
        coef(i,4) =    cfx *   cfy *mask(ip(i)+1,jp(i)+1)
        coef(i,:) = coef(i,:)/sum(coef(i,:))
      enddo

      end subroutine compute_coef !]
! ----------------------------------------------------------------------
      subroutine interpolate_2D(var,vari,coef,ip,jp)  ![
      ! Interpolate a scalar variable
      use dimensions
      implicit none

      ! inputs
      real   ,dimension(-1:nx+1,-1:ny+1),intent(in) :: var  ! assumed size arrays always start at 1.
      real   ,dimension(:)              ,intent(out):: vari
      real   ,dimension(:,:)            ,intent(in) :: coef
      integer,dimension(:)              ,intent(in) :: ip
      integer,dimension(:)              ,intent(in) :: jp

      ! local
      integer :: i,k,np

      np = size(ip,1)
      do i = 1,np
          vari(i) = var(ip(i)  ,jp(i)  )*coef(i,1) + 
     &              var(ip(i)+1,jp(i)  )*coef(i,2) + 
     &              var(ip(i)  ,jp(i)+1)*coef(i,3) + 
     &              var(ip(i)+1,jp(i)+1)*coef(i,4)
      enddo

      end subroutine interpolate_2D  !]
! ----------------------------------------------------------------------
      subroutine interpolate_3D(var,vari,coef,ip,jp)  ![
      ! Interpolate a variable
      use dimensions
      implicit none

      ! inputs
      real   ,dimension(-1:nx+1,-1:ny+1,nz),intent(in) :: var  ! assumed size arrays would always start at 1.
      real   ,dimension(:,:)               ,intent(out):: vari
      real   ,dimension(:,:)               ,intent(in) :: coef
      integer,dimension(:)                 ,intent(in) :: ip
      integer,dimension(:)                 ,intent(in) :: jp

      ! local
      integer :: i,k,np

      np = size(ip,1)
      do i = 1,np
        do k = 1,nz
          vari(i,k) = var(ip(i)  ,jp(i)  ,k)*coef(i,1) + 
     &                var(ip(i)+1,jp(i)  ,k)*coef(i,2) + 
     &                var(ip(i)  ,jp(i)+1,k)*coef(i,3) + 
     &                var(ip(i)+1,jp(i)+1,k)*coef(i,4)
        enddo
      enddo

      end subroutine interpolate_3D  !]
! ----------------------------------------------------------------------
      subroutine interpolate_vec_2D(u,v,ui,vi,cfu,ipu,jpu,cfv,ipv,jpv)  ![
      ! Interpolate a vector on a c-grid and rotate to angle
      use dimensions
      implicit none

      ! inputs
      real   ,dimension(-1:nx+1,-1:ny+1),intent(in) :: u    ! assumed size arrays always start at 1.
      real   ,dimension(-1:nx+1,-1:ny+1),intent(in) :: v    ! assumed size arrays always start at 1.
      real   ,dimension(:)              ,intent(out):: ui
      real   ,dimension(:)              ,intent(out):: vi
!     real   ,dimension(:)              ,intent(in) :: cosa,sina
      real   ,dimension(:,:)            ,intent(in) :: cfu,cfv
      integer,dimension(:)              ,intent(in) :: ipu,ipv
      integer,dimension(:)              ,intent(in) :: jpu,jpv

      ! local
      integer :: i,np

      call interpolate(u,ui,cfu,ipu,jpu)
      call interpolate(v,vi,cfv,ipv,jpv)

!     np = size(ipu,1)
!     do i = 1,np
!        ui(i) = ui(i)*cosa(i) - vi(i)*sina(i)
!        vi(i) = vi(i)*cosa(i) + ui(i)*sina(i)
!     enddo

      end subroutine interpolate_vec_2D  !]
! ----------------------------------------------------------------------
      subroutine interpolate_vec_3D(u,v,ui,vi,cosa,sina,cfu,ipu,jpu,cfv,ipv,jpv)  ![
      ! Interpolate a vector on a c-grid and rotate to angle
      use dimensions
      implicit none

      ! inputs
      real   ,dimension(-1:nx+1,-1:ny+1,nz),intent(in) :: u    ! assumed size arrays always start at 1.
      real   ,dimension(-1:nx+1,-1:ny+1,nz),intent(in) :: v    ! assumed size arrays always start at 1.
      real   ,dimension(:,:)               ,intent(out):: ui
      real   ,dimension(:,:)               ,intent(out):: vi
      real   ,dimension(:)                 ,intent(in) :: cosa,sina
      real   ,dimension(:,:)               ,intent(in) :: cfu,cfv
      integer,dimension(:)                 ,intent(in) :: ipu,ipv
      integer,dimension(:)                 ,intent(in) :: jpu,jpv

      ! local
      integer :: i,k,np

      call interpolate(u,ui,cfu,ipu,jpu)
      call interpolate(v,vi,cfv,ipv,jpv)

      np = size(ipu,1)
      do k = 1,nz
        do i = 1,np
          ui(i,k) = ui(i,k)*cosa(i) - vi(i,k)*sina(i)
          vi(i,k) = vi(i,k)*cosa(i) + ui(i,k)*sina(i)
        enddo
      enddo

      end subroutine interpolate_vec_3D  !]
! ----------------------------------------------------------------------
      subroutine do_extract_data ![
      ! extract data for all objects, for all vars
      ! and write to file
      implicit none

      ! local
      integer :: i,itrc,record
      character(len=30) :: obj_name
      real, dimension(:,:),pointer :: vi
      real, dimension(:,:),pointer :: ui
      real, dimension(:,:),pointer :: coef, cfu, cfv
      integer,dimension(:),pointer :: ip,jp,ipu,ipv,jpu,jpv

      ! For now, extract zeta and all tracers for rho objects
      ! and ubar,vbar,u, and v for velocity objects

      ! different record rates will have to go in different output files....

      if (.not.allocated(obj)) call init_extract_data

!     if (mod(total_rec_diag,nrpf_diag)==0) then
!       call create_diagnostics_file(total_rec_diag,fname)
!       record = 0
!     endif 
!     total_rec_diag  = total_rec_diag +1
!     record = record + 1

!     ierr=nf90_open(fname,nf90_write,ncid)
!     ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)

      do i = 1,nobj
        obj(i)%otime = obj(i)%otime + dt
        if (obj(i)%otime > obj(i)%period ) then
          obj_name = obj(i)%obj_name
          record = obj(i)%record + 1


          if (obj(i)%scalar) then

!           call interpolate(obj(i),zeta(:,:,knew))
!           call write_data(obj(i)%vi(:,1),'zeta',record)


            coef => obj(i)%coef
            ip   => obj(i)%ip
            jp   => obj(i)%jp

            call interpolate(zeta(:,:,knew),obj(i)%vari(:,1),coef,ip,jp)
!           call write_data(vari(:,1),'zeta',record)
            do itrc=1,nt
              call interpolate(t(:,:,:,nstp,itrc),obj(i)%vari,coef,ip,jp)
!             call write_data(vari,obj_name,t_vname(itrc),record)
            enddo
            
          else
            ui => obj(i)%ui
            vi => obj(i)%vi
            cfu=> obj(i)%cfu
            cfv=> obj(i)%cfv
            ipu=> obj(i)%ipu
            ipv=> obj(i)%ipv
            jpu=> obj(i)%jpu
            jpv=> obj(i)%jpv

            call interpolate(ubar(:,:,knew),ui(:,1),cfu,ipu,jpu)
            call interpolate(vbar(:,:,knew),vi(:,1),cfv,ipv,jpv)
!            if (obj(i)%u) then ! DevinD %u (not legal), JM must fix.
              obj(i)%vari(:,1) = obj(i)%cosa*ui(:,1)
!            else
!              vari(:,1) = obj(i)%cosa*vi(:,1)
!            endif

!           call interpolate(ui(:,1),vi(:,1),ubar(:,:,knew),
!    &            vbar(:,:,knew),obj(i))
!           call write_data(vari(:,1),obj_name,'ubar',record)

!           call interpolate(ui,vi,u(:,:,:,knew),v(:,:,:,knew),obj(i))
!           call rotate_vel(vari,ui,vi,cosa,sina)
!           call write_data(ui,'u',record)

!           call interpolate(ui(:,1),vi(:,1),ubar(:,:,knew),
!    &            vbar(:,:,knew),cosa,sina,cfu,ipu,jpu,cfv,ipv,jpv)
!           call write_data(vi(:,1),'vbar',record)
       
          endif
          obj(i)%record = record
          obj(i)%otime = 0
        endif
      enddo

      end subroutine do_extract_data  !]
! ----------------------------------------------------------------------

#endif /* EXTRACT_DATA */
      end module extract_data
