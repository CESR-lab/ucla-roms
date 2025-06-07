      program ncjoin

      ! To distinguish between the partitions and the nodes of this 
      ! code, we'll use mynode,nnodes for the mpi-nodes that this
      ! application uses and npart,part for the number of partitions

      ! ----- Module -----
      use netcdf
      use mpi
      use roms_part_join

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

!     call MPI_Bcast(nparts,1,MPI_INTEGER,0,part_comm,ierr)

      print *, 'Joining partial files: ',mynode,nparts

      ! Make joined file
      call make_joined_file

      call MPI_Barrier(part_comm, ierr)
     
      print *, 'Writing joined file', mynode

      call write_joined_file

      call MPI_Finalize(ierr)

      end program ncjoin

