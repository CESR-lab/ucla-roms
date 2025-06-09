      program optimal_part

      ! Provide a range of 'optimal' partitions that skips
      ! subdomains that are fully masked and ends up with a
      ! number of active subdomains that is close to the amount
      ! of available cores.

      use netcdf
      use mpi
      use roms_part_join

      implicit none


!     integer :: i,nargs 

      integer :: npartx,nparty
      integer :: ierr
      integer :: ifile

      ! Get available cores and grid file through command line argument
      call cmd_args_optimal

  
      ar = gnx/gny  ! aspect ratio grid:
      ocean = sum(msk)/(gnx*gny) ! percentage of grid that is not masked

      ! first guess of x,y partition that gives roughly square subdomains
      ! and should exceed the number of available cores by a bit.
      npx1 = ceiling(sqrt(ncores)*ar/ocean)
      npy1 = ceiling(sqrt(ncores)/ar/ocean)

      ! 16 guesses ok?
      npx0 = max(npx1-4,1)
      npy0 = max(npy1-4,1)
        
      i = 0
      do npy = npy0,npy1
        do npx = npx0,npx1
          call partition_mask(npx,npy)
          i = i+1
          active(i,1) = 
          active(i,2) = npx
          active(i,3) = npy
        enddo
      enddo

      sort(active)


      end program optimal_part

