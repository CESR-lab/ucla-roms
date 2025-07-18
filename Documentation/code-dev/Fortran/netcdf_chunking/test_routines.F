      module test_routines
      
      ! ROMS functionality to facilitate partioning and joining of files

      use mpi
      use nc_read_write
      use netcdf

      private

      ! mpi variables
      integer,public :: part_comm,mynode,nnodes  ! mpi vars

      ! misc.
      real,dimension(:,:,:),allocatable,public :: temp
      integer,dimension(:),allocatable,public :: cs
      integer,public :: surplus_y,loc_y,ny,jstart,chunks_per_node
      integer,dimension(:),allocatable,public :: j0,j1
      integer,dimension(:),allocatable,public :: csizes
      integer,public                          :: gnx,gny,gnz
      character(len=128),public               :: file

      public :: load_data
      public :: create_file
      public :: bychunk
      public :: byrow
      public :: write_file
      
      contains
! ----------------------------------------------------------------------
      subroutine load_data(temp)   ![

      implicit none

      integer :: ncid,ierr,dimid
      real,dimension(:,:,:),allocatable,intent(out) :: temp
      character(len=30) :: vname


      ! read in rst file (pachug preferably)
      ierr = nf90_open(adjustl(trim(file)),nf90_nowrite,ncid)
       if (ierr/=nf90_noerr) call handle_ierr(ierr)

      ! get global dimensions for reading
      ierr = nf90_inq_dimid(ncid,'xi_rho',dimid)
       if (ierr/=nf90_noerr) call handle_ierr(ierr,'getting dimid: ','xi_rho')
      ierr = nf90_inquire_dimension(ncid,dimid,len=gnx)
       if (ierr/=nf90_noerr) call handle_ierr(ierr,'dimension: ','xi_rho')

      ierr = nf90_inq_dimid(ncid,'eta_rho',dimid)
       if (ierr/=nf90_noerr) call handle_ierr(ierr,'getting dimid: ','eta_rho')
      ierr = nf90_inquire_dimension(ncid,dimid,len=gny)
       if (ierr/=nf90_noerr) call handle_ierr(ierr,'dimension: ','eta_rho')

      ierr = nf90_inq_dimid(ncid,'s_rho',dimid)
       if (ierr/=nf90_noerr) call handle_ierr(ierr,'getting dimid: ','s_rho')
      ierr = nf90_inquire_dimension(ncid,dimid,len=gnz)
       if (ierr/=nf90_noerr) call handle_ierr(ierr,'dimension: ','s_rho')

      allocate(temp(gnx,gny,gnz))

      ! grab data
      call ncread(ncid,'temp',temp)

      ierr=nf90_close(ncid)

      end subroutine load_data   !]
! ----------------------------------------------------------------------
      subroutine create_file   ![

      implicit none
      ! local
      integer :: npid,ierr,d1,d2,d3
      integer :: cs,varid,part
      integer,dimension(3) :: dids

      ierr = nf90_create_par('output.nc',nf90_netcdf4,part_comm,MPI_Info_Null,npid)
       if (ierr/=nf90_noerr) call handle_ierr(ierr,'creating file')

      ierr=nf90_def_dim(npid,'xi_rho',gnx,d1)
      ierr=nf90_def_dim(npid,'eta_rho',gny,d2)
      ierr=nf90_def_dim(npid,'s_rho',gnz,d3)
 
      dids=(/d1,d2,d3/)
      csizes=(/gnx,100,100/)

      ierr=nf90_def_var(npid,'temp',NF90_REAL,dids,varid)

      ierr=nf90_def_var_chunking(npid,varid,nf90_chunked,csizes)

      ierr=nf90_var_par_access(npid, nf90_global, nf90_collective)

      ierr=nf90_close(npid)

      end subroutine create_file   !]
! ---------------------------------------------------------------------
      subroutine byrow   ![
 
      implicit none
      ! local
      integer :: ncid,ierr,part

      allocate(j0(nnodes))
      allocate(j1(nnodes))

      ny = ceiling(1.0*(gny-2)/nnodes)
      surplus_y = ny*nnodes - (gny-2)
      loc_y = 1
      do part=0,nnodes-1
        if (part==0) then ! bottom-most subdomain
          jstart = ny - surplus_y/2
        elseif (part==nnodes) then ! upper-most subdomain
          jstart = ny - (surplus_y+1)/2
        else
          jstart = ny
        endif
        j0(part) = loc_y
        j1(part) = loc_y + jstart-1
        loc_y = loc_y + jstart
      enddo


      end subroutine   !]
! ----------------------------------------------------------------------
      subroutine bychunk   ![
 
      implicit none
      ! local
      integer :: ncid,ierr,chunks,part

      chunks = ceiling(1.0*(gny-2)/csizes(2))
      chunks_per_node = ceiling(1.0*chunks/nnodes)

      allocate(j0(chunks))
      allocate(j1(chunks))

      if (mynode==0) print*,'Chunks: ',chunks
      if (mynode==0) print*,'Nodes per Chunk',chunks_per_node

      loc_y = 1
      do part=0,nnodes-1
        if (part==nnodes-1) then ! upper-most subdomain
          jstart = gny-loc_y+1
        elseif (part==0) then
          jstart = chunks_per_node*csizes(2)
        else
          jstart = chunks_per_node*csizes(2)
        endif
        j0(part) = loc_y
        j1(part) = loc_y + jstart-1
        loc_y = loc_y + jstart
      enddo           
    

      end subroutine bychunk   !]
! ----------------------------------------------------------------------
      subroutine write_file(temp)   ![

      implicit none
      integer :: ncid,ierr
      real,dimension(:,:,:),intent(in) :: temp

      ierr = nf90_open_par('output.nc',nf90_write,part_comm,MPI_Info_Null,ncid)
       if (ierr/=nf90_noerr) call handle_ierr(ierr,'opening file')
      call ncwrite(ncid,'temp',temp(:,j0(mynode):j1(mynode),:),(/1,j0(mynode),1,1/))

      call MPI_Barrier(part_comm,ierr)
      ierr = nf90_close(ncid)
     

      end subroutine write_file!]
! ----------------------------------------------------------------------
      end module test_routines
