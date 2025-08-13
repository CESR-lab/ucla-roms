      program test

      use netcdf
      use mpi
      use test_routines

      implicit none

      integer :: ierr,i
      real :: time,start,kill
      character(len=128) :: arg,option

      call MPI_Init(ierr)

      ! Set up nodes 
      part_comm = MPI_COMM_WORLD
      call MPI_Comm_rank(part_comm, mynode, ierr)
      call MPI_Comm_size(part_comm, nnodes, ierr)

      ! Handle command line arguments
      call get_command_argument(1,arg)
      read(arg,*) file
      call get_command_argument(2,arg)
      read(arg,*) option


      ! call routine
      if (mynode==0) print*,'Loading data'

      call load_data(temp)

      if (mynode==0) print*,'Creating file'
      call create_file

      call MPI_Barrier(part_comm,ierr)

      if (option=='row') then
        if (mynode==0) print*,'Testing option a: by row'
        call byrow
      else
        if (mynode==0) print*,'Testing option b: by chunk'
        call bychunk
      endif

      if (mynode==0) print*,'Writing...'
      print*,mynode,j0(mynode),j1(mynode)

      call cpu_time(start)
      ! loop over 5 time to produce some measureable run time
      do i=1,5
        call write_file(temp)
      enddo
      call cpu_time(kill)

      print*,'Time to write: ',kill-start

      call MPI_Finalize(ierr)

      end program test
