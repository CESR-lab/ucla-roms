      program ncjoin

      ! To distinguish between the partitions and the nodes of this 
      ! code, we'll use mynode,nnodes for the mpi-nodes that this
      ! application uses and npart,part for the number of partitions

      ! ----- Module -----
      use netcdf
      use mpi
      use roms_part

      ! ----- Allocation ----- 
      implicit none

      ! Local
!     type(ncstructure) :: ncinfo            ! structure to access netcdf info  
!     integer :: ierr,ncgrd,count,dimid ! netcdf 
!     integer :: gnx,gny                     ! Global domain dimensions
!     integer :: loc_x,loc_y,lnx,lny         ! partition lengths
!     integer :: ifile

      integer :: ierr


      call MPI_Init(ierr)

      ! Set up nodes 
      part_comm = MPI_COMM_WORLD
      call MPI_Comm_rank(part_comm, mynode, ierr)
      call MPI_Comm_size(part_comm, nnodes, ierr)

      ! Handle command line arguments
      call cmd_args_ncjoin

      if (mynode==0) print *, 'Joining partial files'

      ! Make joined file
      call make_joined_file

      call MPI_Barrier(part_comm, ierr)
     
      if (mynode==0) print *, 'Writing joined file'

      call write_joined_file

      if (mynode==0) print *, 'Done writing joined file'

      call MPI_Barrier(part_comm, ierr)
      call MPI_Finalize(ierr)

      end program ncjoin

