      program partit

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

!     integer :: i,nargs 

      integer :: npartx,nparty
      integer :: ierr
      integer :: ifile

      ! Handle system arguments
      call cmd_args_partit(npartx,nparty)

      call MPI_Init(ierr)

      ! Set up nodes 
      part_comm = MPI_COMM_WORLD
      call MPI_Comm_rank(part_comm, mynode, ierr)
      call MPI_Comm_size(part_comm, nnodes, ierr)


      if (mynode==0) then
      ! possibly: all nodes do this, making their on set of partial files
        
        ! figures how many relevant nodes, their location, sizes and neighbors
        call partition_mask(npartx,nparty)

        if (nparts<nnodes) then
          stop 'more nodes than partial files'
        endif

      endif

      call MPI_Bcast(nparts,1,MPI_INTEGER,0,part_comm,ierr)

      print *, 'Making partial files: ',mynode,nparts

!     ! do ifile = 1,nfiles
      ifile = 1
      srcfile = srcfiles(ifile)
      fullsrcfile = trim(spath) // trim(srcfile)

      ! Make partial files 
      if (mynode==0) then
        ! possibly: all nodes do this, making their on set of partial files
        call make_partial_files
      endif

      call MPI_Barrier(part_comm, ierr)

      print *, 'Writing partial files', mynode


      call write_partial_files
      ! read required block of data for partitionable
      ! vars from srcfile and write to partial files

      call MPI_Finalize(ierr)

      end program partit

