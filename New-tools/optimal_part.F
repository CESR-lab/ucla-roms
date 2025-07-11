      program optimal_part

      ! Provide a range of 'optimal' partitions that skips
      ! subdomains that are fully masked and ends up with a
      ! number of active subdomains that is close to the amount
      ! of available cores.

      use netcdf
      use nc_read_write
      use mpi
      use roms_part

      implicit none

      integer :: npartx,nparty
      integer :: ifile

      integer :: ncid,npid,ierr,       ! netcdf error handling
     &           dimid,                ! netcdf dimensions
!    &           gnx,gny,
     &           npx0,npx1,npy0,npy1,  ! loop ranges
     &           npx,npy,i,k
      real                              :: cores,ocean,ar
      real,dimension(:,:),allocatable   :: msk,active

      character(len=7) :: vname
      character(len=7),dimension(2) :: dname
      integer,dimension(2)          :: dsize
      integer                       :: varid

      ! Get available cores and grid file through command line argument
      call cmd_args_optimal(cores,grdfile)

      ! read in data for global domains
      ierr = nf90_open(trim(adjustl(grdfile)), nf90_nowrite, ncid)
       if (ierr/=0) call handle_ierr(ierr,'opening: ',trim(adjustl(grdfile)))

      ! global dimensions
      ierr = nf90_inq_dimid(ncid,'xi_rho',dimid)
       if (ierr/=0) call handle_ierr(ierr,'getting dimid: ','xi_rho')
      ierr = nf90_inquire_dimension(ncid, dimid, len=gnx)
       if (ierr/=0) call handle_ierr(ierr,'getting dimension: ','xi_rho')
        
      ierr = nf90_inq_dimid(ncid,'eta_rho',dimid)
       if (ierr/=0) call handle_ierr(ierr,'getting dimid: ','eta_rho')
      ierr = nf90_inquire_dimension(ncid, dimid, len=gny)
       if (ierr/=0) call handle_ierr(ierr,'getting dimension: ','eta_rho')

      gnx = gnx-2
      gny = gny-2

      ! mask rho
      allocate(msk(0:gnx+1,0:gny+1)); msk=1 
      call ncread(ncid,'mask_rho',msk)

      ierr = nf90_close(ncid)
       if (ierr/=0) call handle_ierr(ierr,'closing: ',trim(adjustl(grdfile)))

      ! make new netcdf to visalize parition
      ierr = nf90_create('part.nc',nf90_netcdf4,npid)
       if (ierr/=0) call handle_ierr(ierr,'creating: ','part.nc')

      ! add variables
      vname='p'
      dname=(/'xi_rho','eta_rho'/)
      dsize=(/gnx+2,gny+2/)
      varid=nccreate(npid,vname,dname,dsize)

      ar    = real(gnx)/real(gny)   ! aspect ratio grid:
      ocean = (sum(msk)/(gnx*gny))  ! fraction of the grid that is not masked

      ! first guess of x,y partition that gives roughly square subdomains
      ! and should exceed the number of available cores by a bit.
      npx1 = ceiling(sqrt(cores*ar/ocean))
      npy1 = ceiling(sqrt(cores/ar/ocean))


      npx0 = max(npx1-7,1)
      npy0 = max(npy1-7,1)

      allocate(active(1:56,1:3))
      i = 0
      do npy = npy0,npy1
        do npx = npx0,npx1
          if ((npx*npy>=cores).and.(npx*npy<=cores/ocean)) then
            call check_mask(npx,npy,gnx,gny,msk)
            i = i+1
            if (nparts<=cores) then
              active(i,1) = nparts
              active(i,2) = npx
              active(i,3) = npy
            endif
          endif
        enddo
      enddo

      ierr = nf90_close(npid)
                     
      call sort(active)

      do k = 1,10
        print *,'---',int(active(k,2)),'x',int(active(k,3)),'---'
        print *,int(active(k,1)),'cores out of',int(cores),'used'
      enddo

      deallocate(active)

      contains

! ----------------------------------------------------------------------
      subroutine sort(active)  ![

      real,dimension(:,:),intent(inout) :: active

      ! local
      integer :: i, j
      real,dimension(size(active, 2)) :: key

      do i  = 2, size(active,1)
        key = active(i,:)
        j   = i - 1

        ! Shift elements to the right
        do while ((j >= 1))
          if (active(j,1) < key(1)) then
            active(j+1,:) = active(j,:)
            j = j - 1
          else
            exit
          endif
        enddo

        ! Insert the element
        active(j+1,:) = key

      enddo
      end subroutine sort  !]
! ----------------------------------------------------------------------

      end program optimal_part
