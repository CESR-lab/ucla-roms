      module extract_data  
      ! Extract data at various positions and frequencies
      ! Only for things that need interpolating


      ! Auxiliary tools required:
      ! - Tools-Roms/scripts/pre/add_object.m -> add a data extraction object to a netcdf file
      ! - Tools-Roms/scripts/pre/bry_extract_obj.m -> Example showing
      !             how to add objects for a child grid boundary forcing file

      !
      ! STEPS:
      !
      ! 0) CREATE INPUT FILE OF EXTRACTION OBJECTS
      !
      !    - Positions (in fractional i,j locations of the grid)
      !    - Name of the object ('child1_south','mooring2', etc, etc,...
      !    - Output frequency of each object
      !    - Variables to be output
      !  
      !
      ! 1) Assign positions to subdomains
      !    - if ispresent(iobj) then see above
      !    - keep track of placement in global object arrays
      !    The ROMS parent simulation needs to know where the child boundary sits in its
      !    parent domain. To do this the child bry points are given i and j coords in
      !    terms of the parent grid. We created a matlab tool to do this found here:
      !    Tools-Roms/scripts/pre/

      ! 2) Create a single file per subdomain containing all the objects
      !    for which ispresent(iobj) is true
      !
      ! 3) Loop through all objects and all vars, and write when needed
      !    Add averaging capability at some point
      !
      !    Vel points always need to interpolate both u, and v in order
      !    to rotate the vector to a desired angle
      !
      !]

#include "cppdefs.opt"

      use dimensions
      use grid
      use nc_read_write
      use roms_read_write
      use netcdf
      use tracers, only: t, t_vname, t_lname, t_units  ! need to get names of tracers
      use ocean_vars
      use basic_output, vn=>vname
      use scalars
      use diagnostics


      ! TODO: add averaging


      implicit none
      private

#include "extract_data.opt"


      integer :: total_rec=0 ! total records
      integer :: nobj  ! number of extraction objects

      logical                         :: extend_up !! flag to extend up,vp
      real,dimension(:,:),allocatable :: upe,vpe  !! buffer filled versions of up,vp
 
      type extract_object  ! contains all information for a data_extraction object

        ! needed as input
        character(len=60) :: obj_name                   ! name of object
        character(len=20) :: set                        ! name of set that the object belongs to
        character(len=20) :: bnd                        ! name of boundary (for bry type data)
        character(:),allocatable :: pre                 ! preamble for vars and dims
        character(:),allocatable :: dname               ! dimension name for object
        integer                  :: dsize               ! dimension of object
        logical                            :: scalar    ! scalar or vector
        real                               :: period    ! output period (seconds)
        integer                            :: np        ! local number of locations
        real,   dimension(:)  ,allocatable :: ipos,jpos ! fractional index locations
        real,   dimension(:)  ,allocatable :: ang       ! desired angle for velocities
        integer                            :: start_idx ! starting position in the obj array


        !! Initializing record at nrpf will trigger the making of a file
        !! at the first time step
        integer                        :: record=nrpf   ! record number in file
        real                               :: otime=0   ! time since last output

        real,   dimension(:,:),pointer     :: vari      ! data for output

        ! these are only for scalars
        integer,dimension(:)  ,pointer     :: ip,jp     ! lower left hand point index
        real   ,dimension(:,:),pointer     :: coef      ! interpolation coefficients

        ! these are only for velocities
        real,   dimension(:)  ,allocatable :: cosa,sina ! for rotation of vectors
        integer,dimension(:)  ,pointer     :: ipu,jpu   ! only for vectors
        integer,dimension(:)  ,pointer     :: ipv,jpv   ! only for vectors
        real   ,dimension(:,:),pointer     :: cfu,cfv   ! only for vectors
        real   ,dimension(:,:),pointer     :: ui,vi     ! only for vectors 

        ! These logicals determine which variables are desired for an
        ! object. Only the ones listed below are currently functional
        ! False by default
        logical :: zeta = .false.
        logical :: ubar = .false.
        logical :: vbar = .false.
        logical :: u    = .false.
        logical :: v    = .false.
        logical :: temp = .false.
        logical :: salt = .false.
        logical :: up   = .false.
        logical :: vp   = .false.

      end type extract_object

      type(extract_object),dimension(:),allocatable :: obj

      interface interpolate
        module procedure  interpolate_2D, interpolate_3D
      end interface

      public do_extract_data
      contains
! ----------------------------------------------------------------------
      subroutine init_extract_data  ![
      ! Allocate space and compute interpolation coefficients
      ! for rho,u,and v variables.
      implicit none

      !local
      character(len=30) :: preamb
      integer :: i,np,ierr,lpre
      real,dimension(:),allocatable :: angp  ! grid angle

      call read_extraction_objects

      if (diag_pflx) then
        allocate( upe(GLOBAL_2D_ARRAY) ); upe = 0
        allocate( vpe(GLOBAL_2D_ARRAY) ); vpe = 0
      endif

      do i = 1,nobj
        np = obj(i)%np
        if (np>0) then

          preamb = trim(obj(i)%obj_name)
          lpre = len(trim(preamb))-1
          allocate(character(len=lpre) :: obj(i)%pre)
          obj(i)%pre = preamb(1:lpre)

          allocate(obj(i)%vari(np,nz))
          allocate(obj(i)%coef(np,4))
          allocate(obj(i)%ip(np))                          ! DON'T NEED THESE FOR U or V objects
          allocate(obj(i)%jp(np))

          ! from absolute index to rho-index
          ! ipos,jpos are in 'absolute' index space: [0,nx]x[0,ny]
          obj(i)%ipos = obj(i)%ipos+0.5
          obj(i)%jpos = obj(i)%jpos+0.5
          call compute_coef(obj(i)%ipos,obj(i)%jpos,
     &                obj(i)%coef,obj(i)%ip,obj(i)%jp,rmask)

          if (.not.obj(i)%scalar) then
            allocate(obj(i)%cosa(np))
            allocate(obj(i)%sina(np))

            allocate(angp(np))
            call interpolate(angler,angp,obj(i)%coef,obj(i)%ip,obj(i)%jp)
            obj(i)%cosa = cos(angp-obj(i)%ang)
            obj(i)%sina = sin(angp-obj(i)%ang)

            deallocate(angp)

            allocate(obj(i)%ui(np,nz))
            allocate(obj(i)%cfu(np,4))
            allocate(obj(i)%ipu(np))
            allocate(obj(i)%jpu(np))

            ! from rho-index to u-index 
            obj(i)%ipos = obj(i)%ipos+0.5
            call compute_coef(obj(i)%ipos,obj(i)%jpos,
     &                obj(i)%cfu,obj(i)%ipu,obj(i)%jpu,umask)

            allocate(obj(i)%vi(obj(i)%np,nz))
            allocate(obj(i)%cfv(obj(i)%np,4))
            allocate(obj(i)%ipv(obj(i)%np))
            allocate(obj(i)%jpv(obj(i)%np))

            ! from u-index to v-index
            obj(i)%ipos = obj(i)%ipos-0.5
            obj(i)%jpos = obj(i)%jpos+0.5
            call compute_coef(obj(i)%ipos,obj(i)%jpos,
     &                obj(i)%cfv,obj(i)%ipv,obj(i)%jpv,vmask)

          endif
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
      character(len=20) :: dname
      character(len=30) :: objname
      character(len=80) :: variables
      real,dimension(:,:),allocatable :: object
      integer               :: n1,n2,i0,i1,lstr
      integer ierr,sidx,sidx2
      real :: period

      if (mynode==0) then
        write(*,'(7x,2A)')
     &    'extract_data :: read objects: ',extract_file
      endif
      ierr = nf90_open(extract_file,nf90_nowrite,ncid)
      if (ierr/=0)
     &  call handle_ierr(ierr,'extract_data :: read objects file', extract_file)
      ierr = nf90_inquire(ncid,nVariables=nobj)
      if (ierr/=0) call handle_ierr(ierr,'read extr ')

      allocate(obj(nobj))

      ! Read all objects from file. 
      do iobj = 1,nobj
        ierr = nf90_inquire_variable(ncid,iobj,name=objname,dimids=dimids)
        if (ierr/=0) call handle_ierr(ierr,'read extr ')
        ierr = nf90_inquire_dimension(ncid,dimids(1),name=dname,len=n1)
        if (ierr/=0) call handle_ierr(ierr,'read extr ')
        ierr = nf90_inquire_dimension(ncid,dimids(2),len=n2)
        if (ierr/=0) call handle_ierr(ierr,'read extr ')

        lstr = len(trim(dname))
        allocate(character(len=lstr) :: obj(iobj)%dname)
        obj(iobj)%dname = trim(dname)
        obj(iobj)%dsize = n1
        
        ! scalar objects have i,j. Vector obj also have ang
        if (n2==2) obj(iobj)%scalar = .true.

        allocate(object(n1,n2))
        call ncread(ncid,objname,object)

        call find_local_points(obj(iobj),object(:,1),object(:,2))

        obj(iobj)%obj_name = objname
        lstr = len(trim(objname))
        if (findstr(objname,'_',sidx) ) then
          obj(iobj)%set = objname(1:sidx-1)
          if (findstr(objname(sidx+1:lstr),'_',sidx2) ) then
            obj(iobj)%bnd = objname(sidx:sidx+sidx2-1)
          else
            obj(iobj)%bnd = ' '
          endif
        else
          obj(iobj)%set = objname
        endif

!       if (mynode==84) then
!         print *,i,obj(iobj)%obj_name,obj(iobj)%set,obj(iobj)%bnd
!         print *,trim(obj(iobj)%set)//'_'//vname//'_'//trim(obj(iobj)%bnd)
!         stop
!       endif

        ierr = nf90_get_att(ncid,iobj,'output_period',period)
        if (ierr/=0) call handle_ierr(ierr,'read extr ')
        obj(iobj)%period   = period

        if (obj(iobj)%np>0) then
          !! only for objects with a presences in this subdomain

          if (n2==3) then
            allocate(obj(iobj)%ang(obj(iobj)%np))
            i0 = obj(iobj)%start_idx
            i1 = i0+obj(iobj)%np-1
            obj(iobj)%ang = object(i0:i1,3)
          endif

          ! Figure out which variables to output for this object
          ierr = nf90_get_att(ncid,iobj,'output_vars',variables)
          if (ierr/=0) call handle_ierr(ierr,'read extr ')

          if (findstr(variables,'zeta') ) obj(iobj)%zeta = .True.
          if (findstr(variables,'temp') ) obj(iobj)%temp = .True.
          if (findstr(variables,'salt') ) obj(iobj)%salt = .True.
          if (findstr(variables,'ubar') ) obj(iobj)%ubar = .True.
          if (findstr(variables,'vbar') ) obj(iobj)%vbar = .True.
          if (findstr(variables,'u'   ) ) obj(iobj)%u    = .True.
          if (findstr(variables,'v'   ) ) obj(iobj)%v    = .True.
          if (findstr(variables,'up'  ) ) obj(iobj)%up   = .True.
          if (findstr(variables,'vp'  ) ) obj(iobj)%vp   = .True.

          ! DIAGNOSTICS needs to be activated and diag_pflx=.true.
          if (obj(iobj)%up.and. .not.diag_pflx) then
            print *,'Fatal ERROR: up extraction, diag_pflx is not set'
            stop
          endif
          if (obj(iobj)%vp.and. .not.diag_pflx) then
            print *,'Fatal ERROR: vp extraction, diag_pflx is not set'
            stop
          endif

          if (obj(iobj)%up.or.obj(iobj)%vp) extend_up = .true.

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
      integer :: np
      
      np = size(gobj_i)
      gobj_i = gobj_i-iSW_corn
      gobj_j = gobj_j-jSW_corn

      ! Assume that local ranges of objects are contiguous
      start_idx = 0
      do i = 1,np
        if ( gobj_i(i)>=0.and.gobj_i(i)<nx .and.
     &       gobj_j(i)>=0.and.gobj_j(i)<ny ) then
          start_idx = i
          exit
        endif
      enddo

      end_idx = np
      if (start_idx>0) then
        do i = start_idx,np
          if (gobj_i(i)<0.or.gobj_i(i)>=nx.or.
     &        gobj_j(i)<0.or.gobj_j(i)>=ny ) then
            end_idx = i-1
            exit
          endif
        enddo
        obj%np = end_idx - start_idx + 1
      else ! object not in local range
        obj%np = 0
      endif

      if (obj%np>0) then
        obj%start_idx = start_idx
        obj%np = end_idx-start_idx+1
        allocate(obj%ipos(obj%np))
        allocate(obj%jpos(obj%np))
        obj%ipos = gobj_i(start_idx:end_idx)
        obj%jpos = gobj_j(start_idx:end_idx)
      endif

      end subroutine find_local_points  !]
! ----------------------------------------------------------------------
      subroutine compute_coef(ipos,jpos,coef,ip,jp,mask)  ![
      ! compute interpolation coefficients
      implicit none

      ! inport/export
      real   ,dimension(:)  ,intent(in) :: ipos,jpos
      real   ,dimension(:,:),intent(out):: coef
      integer,dimension(:)  ,intent(out):: ip,jp
      real   ,dimension(-1:nx+2,-1:ny+2),intent(in) :: mask

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
        !! possibly check for all masked ....
        coef(i,:) = coef(i,:)/sum(coef(i,:))
      enddo

      end subroutine compute_coef !]
! ----------------------------------------------------------------------
      subroutine interpolate_2D(var,vari,coef,ip,jp)  ![
      ! Interpolate a scalar variable
      use dimensions
      implicit none

      ! inputs
      real   ,dimension(-1:nx+2,-1:ny+2),intent(in) :: var  ! assumed size arrays always start at 1.
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
      real   ,dimension(-1:nx+2,-1:ny+2,nz),intent(in) :: var  ! assumed size arrays would always start at 1.
      real   ,dimension(:,:)               ,intent(out):: vari
      real   ,dimension(:,:)               ,intent(in) :: coef
      integer,dimension(:)                 ,intent(in) :: ip
      integer,dimension(:)                 ,intent(in) :: jp

      ! local
      integer :: i,k,np


!     if (mynode==1) then
!       print *,'interp'
!       print *,shape(var)
!       print *, ip(10),jp(10)
!       print *, nx,ny
!       print *, var(ip(10),jp(10),nz)
!     endif
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
      subroutine do_extract_data ![
      ! extract data for all objects, for all vars
      ! and write to file
      implicit none

      ! local
      integer :: i,itrc,ierr,ncid,k,record
      character(len=30) :: obj_name
      character(len=99),save :: fname
      character(len=20)              :: tname
      character(len=40) :: oname
      integer :: lpre
      real, dimension(:,:),pointer :: vi
      real, dimension(:,:),pointer :: ui
      real, dimension(:,:),pointer :: coef, cfu, cfv
      integer,dimension(:),pointer :: ip,jp,ipu,ipv,jpu,jpv
      real,dimension(:),allocatable :: dummy
      integer,dimension(3) :: start2D
      integer,dimension(2) :: start1D

      if (.not.allocated(obj)) then
        call init_extract_data
        return
      endif

      if (obj(1)%record==nrpf) then
          call create_edata_file(fname)
          obj(:)%record = 0
      endif

      if (diag_pflx) then
        upe(1:nx,1:ny) = up
        vpe(1:nx,1:ny) = vp
# ifdef EXCHANGE   
        call exchange_xxx(upe,vpe)
# endif
      endif

      ierr=nf90_open(fname,nf90_write,ncid)
!     ierr=nf90_set_fill(ncid, nf90_nofill, prev_fill_mode)

      !! We have to update the object records regardless of whether
      !! there are points in the sub-domain to ensure correct file
      !! names for all
      do i = 1,nobj
        obj(i)%otime = obj(i)%otime + dt
        if (obj(i)%otime>=obj(i)%period) then
          obj(i)%record = obj(i)%record+1
          obj(i)%otime = 0
          if (i==1) total_rec = total_rec+1
          if (obj(i)%np>0) then

            record = obj(i)%record
            start1D = (/1, record/)
            start2D = (/1,1,record/)

            coef => obj(i)%coef
            ip   => obj(i)%ip
            jp   => obj(i)%jp

            if (.not.obj(i)%scalar) then
              ui => obj(i)%ui
              vi => obj(i)%vi
              cfu=> obj(i)%cfu
              cfv=> obj(i)%cfv
              ipu=> obj(i)%ipu
              ipv=> obj(i)%ipv
              jpu=> obj(i)%jpu
              jpv=> obj(i)%jpv
            endif

            tname = trim(obj(i)%set)//'_time'
            call ncwrite(ncid,tname,(/time/),(/record/))
            !if (mynode==0) print *,'writing extract: ',time,mynode,tname

            if (obj(i)%zeta) then
              oname = trim(obj(i)%set)//'_zeta'//trim(obj(i)%bnd)
              call interpolate(zeta(:,:,knew),obj(i)%vari(:,1),coef,ip,jp)
              call ncwrite(ncid,oname,obj(i)%vari(:,1),start1D)
            endif
            if (obj(i)%temp) then
              oname = trim(obj(i)%set)//'_temp'//trim(obj(i)%bnd)
              call interpolate(t(:,:,:,nstp,itemp),obj(i)%vari,coef,ip,jp)
              call ncwrite(ncid,oname,obj(i)%vari,start2D)
            endif
#ifdef SALINITY
            if (obj(i)%salt) then
              oname = trim(obj(i)%set)//'_salt'//trim(obj(i)%bnd)
              call interpolate(t(:,:,:,nstp,isalt),obj(i)%vari,coef,ip,jp)
              call ncwrite(ncid,oname,obj(i)%vari,start2D)
            endif
#endif

            if (obj(i)%ubar) then
              oname = trim(obj(i)%set)//'_ubar'//trim(obj(i)%bnd)
              call interpolate(ubar(:,:,knew),ui(:,1),cfu,ipu,jpu)
              call interpolate(vbar(:,:,knew),vi(:,1),cfv,ipv,jpv)
              obj(i)%vari(:,1) = obj(i)%cosa*ui(:,1) - obj(i)%sina*vi(:,1)
              call ncwrite(ncid,oname,obj(i)%vari(:,1),start1D)
            endif
            if (obj(i)%vbar) then
              oname = trim(obj(i)%set)//'_vbar'//trim(obj(i)%bnd)
              call interpolate(ubar(:,:,knew),ui(:,1),cfu,ipu,jpu)
              call interpolate(vbar(:,:,knew),vi(:,1),cfv,ipv,jpv)
              obj(i)%vari(:,1) = obj(i)%sina*ui(:,1) + obj(i)%cosa*vi(:,1)
              call ncwrite(ncid,oname,obj(i)%vari(:,1),start1D)
            endif
            if (obj(i)%u) then
              call interpolate(u(:,:,:,nstp),ui,cfu,ipu,jpu)
              call interpolate(v(:,:,:,nstp),vi,cfv,ipv,jpv)
              obj(i)%vari = ui
              do k=1,nz
                obj(i)%vari(:,k) = obj(i)%cosa*ui(:,k) - obj(i)%sina*vi(:,k)
              enddo
              oname = trim(obj(i)%set)//'_u'//trim(obj(i)%bnd)
              call ncwrite(ncid,oname,obj(i)%vari,start2D)
            endif
            if (obj(i)%v) then
              call interpolate(u(:,:,:,nstp),ui,cfu,ipu,jpu)
              call interpolate(v(:,:,:,nstp),vi,cfv,ipv,jpv)
              do k=1,nz
                obj(i)%vari(:,k) = obj(i)%sina*ui(:,k) + obj(i)%cosa*vi(:,k)
              enddo
              oname = trim(obj(i)%set)//'_v'//trim(obj(i)%bnd)
              call ncwrite(ncid,oname,obj(i)%vari,start2D)
            endif
            if (obj(i)%up) then
              call interpolate(upe,ui(:,1),cfu,ipu,jpu)
              call interpolate(vpe,vi(:,1),cfv,ipv,jpv)
              obj(i)%vari(:,1) = obj(i)%cosa*ui(:,1) - obj(i)%sina*vi(:,1)
              oname = trim(obj(i)%set)//'_up'//trim(obj(i)%bnd)
              call ncwrite(ncid,oname,obj(i)%vari(:,1),start1D)
            endif
            if (obj(i)%vp) then
              call interpolate(upe,ui(:,1),cfu,ipu,jpu)
              call interpolate(vpe,vi(:,1),cfv,ipv,jpv)
              obj(i)%vari(:,1) = obj(i)%sina*ui(:,1) + obj(i)%cosa*vi(:,1)
              oname = trim(obj(i)%set)//'_vp'//trim(obj(i)%bnd)
              call ncwrite(ncid,oname,obj(i)%vari(:,1),start1D)
            endif
          endif
        endif
      enddo

      ierr=nf90_close(ncid)

      end subroutine do_extract_data  !]
! ----------------------------------------------------------------------
      subroutine create_edata_file(fname) ![
      implicit none

      !input/output
      character(len=99),intent(out) :: fname

      !local
      integer :: ncid,ierr

      call create_file('_ext',fname)

      ierr=nf90_open(fname,nf90_write,ncid)

      call create_edata_vars(ncid)

      ierr = nf90_close(ncid)

      end subroutine create_edata_file !]
! ----------------------------------------------------------------------
      subroutine create_edata_vars(ncid)  ![
      ! Add edata variables to an opened netcdf file

      use nc_read_write
      implicit none

      !import/export
      integer, intent(in) :: ncid

      !local
      character(len=20)              :: vname
      character(len=20)              :: tname
      integer                        :: varid,ierr
      integer                        :: it,i,lpre
      character(len=3) :: label
      character(len=20),dimension(2) :: dname ! dimension names
      integer,          dimension(2) :: dsize ! dim lengths
      character(len=20),dimension(3) :: dname3 ! dimension names
      integer,          dimension(3) :: dsize3 ! dim lengths

      do i = 1,nobj
        if (obj(i)%np>0) then
          write(label,'(I0.3)') i
          dname = (/ 'np'//label, 'time_'//trim(obj(i)%set) /)
          dsize = (/ obj(i)%np ,  0/)
          dname3(1) = 'np'//label
          dname3(2) = 's_rho'
          dname3(3) = 'time_'//trim(obj(i)%set)
          dsize3 = (/ obj(i)%np , nz, 0/)

          tname = trim(obj(i)%set)//'_time'
          ierr=nf90_inq_varid(ncid,tname,varid)    
          if (ierr/=0) then   ! Only create if not already present
            varid = nccreate(ncid,tname,(/dname(2)/),(/0/),nf90_double)
            ierr = nf90_put_att(ncid,varid,'long_name', 
     &                                     'Time since 2000')
            ierr = nf90_put_att(ncid,varid,'units','second' )
          endif

          if (obj(i)%zeta) call create_var(ncid,obj(i),'zeta',dname,dsize,indxZ)
          if (obj(i)%ubar) call create_var(ncid,obj(i),'ubar',dname,dsize,indxUb)
          if (obj(i)%vbar) call create_var(ncid,obj(i),'vbar',dname,dsize,indxVb)
          if (obj(i)%temp) call create_var(ncid,obj(i),'temp',dname3,dsize3,indxT)
#ifdef SALINITY
          if (obj(i)%salt) call create_var(ncid,obj(i),'salt',dname3,dsize3,indxS)
#endif
          if (obj(i)%u)    call create_var(ncid,obj(i),'u',dname3,dsize3,indxU)
          if (obj(i)%v)    call create_var(ncid,obj(i),'v',dname3,dsize3,indxV)
          if (obj(i)%up)   call create_var(ncid,obj(i),'up',dname,dsize)
          if (obj(i)%vp)   call create_var(ncid,obj(i),'vp',dname,dsize)
        endif
      enddo

      end subroutine create_edata_vars  !]
! ----------------------------------------------------------------------
      subroutine create_var(ncid,obj,vname,dname,dsize,idx) ![
      implicit none

      ! import/export
      integer :: ncid
      type(extract_object),intent(in) ::obj
      character(len=*)              ,intent(in) :: vname
      character(len=20),dimension(:),intent(in) :: dname ! dimension names
      integer,          dimension(:),intent(in) :: dsize ! dim lengths
      integer,optional,              intent(in) :: idx

      ! local
      integer :: varid,ierr
      character(len=40) :: oname

      oname = trim(obj%set)//'_'//vname//trim(obj%bnd)

      varid = nccreate(ncid,oname,dname,dsize,edat_prec)
      ierr = nf90_put_att(ncid,varid,'start',obj%start_idx)
      ierr = nf90_put_att(ncid,varid,'count',obj%np)
      ierr = nf90_put_att(ncid,varid,'dname',obj%dname)
      ierr = nf90_put_att(ncid,varid,'dsize',obj%dsize)
      if (ierr/=0) call handle_ierr(ierr,'extract create var: '//vname)

      if (present(idx)) then
        ierr = nf90_put_att(ncid,varid,'long_name',vn(2,idx))
        ierr = nf90_put_att(ncid,varid,'units',vn(3,idx))
      else
        ierr = nf90_put_att(ncid,varid,'long_name','pressure flux')
        ierr = nf90_put_att(ncid,varid,'units','W/m^2')
      endif

      end subroutine create_var  !]

! ----------------------------------------------------------------------

      end module extract_data
