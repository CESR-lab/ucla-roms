module mg_mpi

  use mg_cst
  use mpi

  implicit none

  integer(kind=ip) :: myrank
  integer(kind=ip) :: nprocs

contains
  !---------------------------------------------------------------------
  subroutine mg_mpi_init()

    call mpi_myrank()

    call mpi_nprocs()

  end subroutine mg_mpi_init

  !---------------------------------------------------------------------
  subroutine mpi_myrank()

    integer(kind=ip) :: ierr

    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

  end subroutine mpi_myrank

  !---------------------------------------------------------------------
  subroutine mpi_nprocs()

    integer(kind=ip) :: ierr

    call mpi_comm_size(mpi_comm_world, nprocs, ierr)

  end subroutine mpi_nprocs

end module mg_mpi
