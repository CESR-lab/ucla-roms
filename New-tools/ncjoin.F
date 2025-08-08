      program ncjoin

      ! To distinguish between the partitions and the nodes of this 
      ! code, we'll use mynode,nnodes for the mpi-nodes that this
      ! application uses and npart,part for the number of partitions

      ! ----- Module -----
      use netcdf
      use mpi
      use roms_part
      use IFPORT

      ! ----- Allocation ----- 
      implicit none

      ! Local
!     type(ncstructure) :: ncinfo            ! structure to access netcdf info  
!     integer :: ierr,ncgrd,count,dimid ! netcdf 
!     integer :: gnx,gny                     ! Global domain dimensions
!     integer :: loc_x,loc_y,lnx,lny         ! partition lengths
!     integer :: ifile

      integer                        :: ierr,files_removed,i,idx
      integer                        :: stridx,endidx
      character(len=:), allocatable  :: command,file,newfile
      type(filestruct)               :: joinlist


      call MPI_Init(ierr)

      ! Set up nodes 
      part_comm = MPI_COMM_WORLD
      call MPI_Comm_rank(part_comm, mynode, ierr)
      call MPI_Comm_size(part_comm, nnodes, ierr)

      ! Handle command line arguments
      call cmd_args_ncjoin

      call MPI_Barrier(part_comm, ierr)
      ! loop through all files that need to be joined
      do i=1,size(file_indices)
      if (MOD(i,2)/=0) then
      stridx=file_indices(i); endidx=file_indices(i+1)
      allocate(partfiles(0:endidx-stridx))
      
      !print*,'number of files: ', size(partfiles)
      ! with new set of files, set definitions
      partfiles=allfiles(stridx:endidx)
      !print*,'File Ranges: ',partfiles
      nparts   =size(partfiles)
      part0    =0  ! always zero, but placed here for clarity
      part1    =nparts-1

      if (mynode==0) print *, 'Joining partial files'

      ! Make joined file
      call make_joined_file

      call MPI_Barrier(part_comm, ierr)
     
      if (mynode==0) print *, 'Writing joined file'

      call write_joined_file

      if (mynode==0) print *, 'Done writing joined file'

      deallocate(partfiles)
      endif
      call MPI_Barrier(part_comm, ierr)
      enddo  ! <-- file loop

      call MPI_Barrier(part_comm, ierr)
      call MPI_Finalize(ierr)

      if (mynode==0) then  ! avoid multiple cores deleting files
      if (del_part_files) then
       print*,'Deleting partial files'
       do i=0,size(partfiles)-1
         command = 'rm -r -- "' // trim(partfiles(i)) // '"'
         call execute_command_line(command)
       enddo
       print*,'Removed ',size(partfiles),' files'
      endif
      endif  ! <-- mynode=0

      end program ncjoin

