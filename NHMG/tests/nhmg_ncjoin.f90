!-----------------------------------------------------------------------
! A simple program to join NHMG netcdf files 
!
! Arguments :
!   1- the number of procs in x
!   2- the number of procs in y
!   3- Number of dimension and order [ij, ijk, ji, kji, dkji]
!   4- ghost points : before 1 in x, after nx in x, before 1 in y, 
!      after ny in y [0000, 1111, 0111,1101, 2222, ...]
!   5..n - the netcdf file names in good order in term of rank
!
! Examples : 
! nhmg_join 2 2  ijk 1111 nc0.nc nc1.nc nc2.nc nc3.nc
! or
! nhmg_join 2 2   ij 1111 h*.nc
! or
! nhmg_join 2 2  ijk 0111 u*.nc
! or
! nhmg_join 2 2  ijk 1101 v*.nc
! or
! nhmg_join 2 2  ijk 1111 w*.nc
! or
! nhmg_join 2 2  kji 2222 zr*.nc
! or
! nhmg_join 2 2 dkji 1111 cA*.nc  <----- matrix coefficients
!
!--------------------------------------------------------------------------
!- Nicolas.Grima@univ-brest.fr january 2017
!--------------------------------------------------------------------------

program nhmg_ncjoin

  use netcdf

  implicit none

  character(:)                   , allocatable :: a_single_argument
  character(len=64), dimension(:), allocatable :: ncfiles

  character(len = 64) :: final_netcdf_file_name

  character(len = 32) :: vname

  character(len =  4) :: dimsorder  !-> ij, ijk, ji, kji, dkji

  character(len =  1) :: cdum

  real(kind = 8), dimension(:,:)    , allocatable :: a2D
  real(kind = 8), dimension(:,:,:)  , allocatable :: a3D
  real(kind = 8), dimension(:,:,:,:), allocatable :: a4D

  integer(kind = 4), dimension(:)   , allocatable :: sstart
  integer(kind = 4), dimension(:)   , allocatable :: scount
  integer(kind = 4), dimension(:,:) , allocatable :: gstart

  integer(kind = 4) :: sz
  integer(kind = 4) :: stat
  integer(kind = 4) :: ii, jj
  integer(kind = 4) :: countargs
  integer(kind = 4) :: ierr
  integer(kind = 4) :: ncid, gncid

  integer(kind = 4) :: nbdim
  integer(kind = 4) :: dim1ID, dim2ID, dim3ID, dim4ID
  integer(kind = 4) :: dim1, dim2, dim3, dim4
  integer(kind = 4) :: gdim1, gdim2, gdim3, gdim4
  integer(kind = 4) :: stii, stjj
  integer(kind = 4) :: xs, xe, ys, ye

  integer(kind = 4) :: counter

  integer(kind = 4) :: varID
  integer(kind = 4) :: gvarID

  integer(kind = 4) :: npx, npy
  integer(kind = 4) :: npxy

  integer(kind = 4) :: parg = 4

  !=================!
  != Get arguments =!
  !=================!
  countargs = command_argument_count()
  print "('Number of arguments:',I2)", countargs

  if (countargs > parg) then

     do ii=1,parg

        call get_command_argument(ii, LENGTH=sz, STATUS=stat)
        allocate(character(sz) :: a_single_argument)
        call get_command_argument(ii, value=a_single_argument, STATUS=stat)

        if (ii == 1) then !- Arg1 = npx

           call str2int(a_single_argument,npx,ierr)
           print "('Nb proc in x:',I2)",npx

        elseif (ii == 2) then !- Arg2 = npy

           call str2int(a_single_argument,npy,ierr)
           print "('Nb proc in y:',I2)",npx
           npxy = npx * npy

        elseif (ii == 3) then !- Arg3 = dimension order

           dimsorder = trim(a_single_argument)

        elseif (ii == 4) then !- Arg4 = ghost points in x and y

           cdum = a_single_argument(1:1)
           call str2int(cdum,xs,ierr)
           cdum = a_single_argument(2:2)
           call str2int(cdum,xe,ierr)
           cdum = a_single_argument(3:3)
           call str2int(cdum,ys,ierr)
           cdum = a_single_argument(4:4)
           call str2int(cdum,ye,ierr)

           print "('ghost points:',I2,' ',I2,' ',I2,' ',I2)",xs,xe,ys,ye

           !---------------------------------!
           !- Reading the netcdf file names -!
           !---------------------------------!
           if ((countargs-parg) == npxy) then
              allocate(ncfiles(npxy))
              do jj = parg+1, countargs
                 call get_command_argument(jj, LENGTH=sz, STATUS=stat)
                 call get_command_argument(jj, value=ncfiles(jj-parg), STATUS=stat)
                 print "('Netcdf file ',I2,': ',A)",jj-parg,trim(ncfiles(jj-parg))
              enddo
           else
              stop "Sorry, but the number of netcdf files is not complet"
           endif

        endif

        deallocate(a_single_argument)

     enddo
  else
     stop "Please, read the help in the header of the nhmg_ncjoin.f90 file."
  endif

  !=============================!
  !- Read subdomain dimensions -!
  !=============================!
  !- depend of ijk or kji order !

  if (trim(dimsorder)=='ij') then

     write(*,*)'-> IJ dims order <-'

     ierr = nf90_open(trim(ncfiles(1)), nf90_nowrite, ncid)
     if (ierr /= nf90_noerr) stop "Error to open Netcdf file"

     !- dim 1 (x)
     ierr = nf90_inq_dimid(ncid, "dim1", Dim1ID)
     if (ierr /= nf90_noerr) stop "Error to read dim1 ID!"
     ierr = nf90_inquire_dimension(ncid, Dim1ID, len = dim1)
     dim1 = dim1 - xs - xe
     gdim1 = dim1 * npx

     !- dim 2 (y)
     ierr = nf90_inq_dimid(ncid, "dim2", Dim2ID)
     if (ierr /= nf90_noerr) stop "Error to read dim2 ID!"
     ierr = nf90_inquire_dimension(ncid, Dim2ID, len = dim2)
     dim2 = dim2 - ys - ye
     gdim2 = dim2 * npy

     nbdim = 2

     allocate(sstart(nbdim))
     allocate(scount(nbdim))

     sstart=[xs+1,ys+1]
     scount=[dim1,dim2]

     allocate(gstart(npxy,nbdim))
     counter = 1
     do jj=1,npy
        stjj = (jj-1) * dim2  + 1
        do ii=1,npx
           stii = (ii-1) * dim1  + 1
           gstart(counter,:) = [stii,stjj]
           counter = counter + 1
        enddo
     enddo

  elseif (trim(dimsorder)=='ijk') then

     write(*,*)'-> IJK dims order <-'

     ierr = nf90_open(trim(ncfiles(1)), nf90_nowrite, ncid)
     if (ierr /= nf90_noerr) stop "Error to open Netcdf file"

     !- dim 1 (x)
     ierr = nf90_inq_dimid(ncid, "dim1", Dim1ID)
     if (ierr /= nf90_noerr) stop "Error to read dim1 ID!"
     ierr = nf90_inquire_dimension(ncid, Dim1ID, len = dim1)
     dim1 = dim1 - xs - xe
     gdim1 = dim1 * npx

     !- dim 2 (y)
     ierr = nf90_inq_dimid(ncid, "dim2", Dim2ID)
     if (ierr /= nf90_noerr) stop "Error to read dim2 ID!"
     ierr = nf90_inquire_dimension(ncid, Dim2ID, len = dim2)
     dim2 = dim2 - ys - ye
     gdim2 = dim2 * npy

     !- dim 3 (z if available)
     ierr = nf90_inq_dimid(ncid, "dim3", Dim3ID)
     if (ierr /= nf90_noerr) stop "Error to read dim3 ID!"
     ierr = nf90_inquire_dimension(ncid, Dim3ID, len = dim3)
     gdim3 = dim3

     nbdim = 3

     allocate(sstart(nbdim))
     allocate(scount(nbdim))

     sstart=[xs+1,ys+1,1]
     scount=[dim1,dim2,dim3]

     allocate(gstart(npxy,nbdim))
     counter = 1
     do jj=1,npx
        stjj = (jj-1) * dim1  + 1
        do ii=1,npy
           stii = (ii-1) * dim2  + 1
           gstart(counter,:) = [1,stjj,stii]
           counter = counter + 1
        enddo
     enddo

  elseif ((trim(dimsorder)=='ji')) then

     write(*,*)'-> JI dims order <-'

     ierr = nf90_open(trim(ncfiles(1)), nf90_nowrite, ncid)
     if (ierr /= nf90_noerr) stop "Error to open Netcdf file"

     !- dim 1 (y)
     ierr = nf90_inq_dimid(ncid, "dim1", Dim1ID)
     if (ierr /= nf90_noerr) stop "Error to read dim1 ID!"
     ierr = nf90_inquire_dimension(ncid, Dim1ID, len = dim1)
     dim1 = dim1
     gdim1 = dim1

     !- dim 2 (x)
     ierr = nf90_inq_dimid(ncid, "dim2", Dim2ID)
     if (ierr /= nf90_noerr) stop "Error to read dim2 ID!"
     ierr = nf90_inquire_dimension(ncid, Dim2ID, len = dim2)
     dim2 = dim2 - ys - ye
     gdim2 = dim2 * npy

     nbdim = 2

     allocate(sstart(nbdim))
     allocate(scount(nbdim))

     sstart=[ys+1,xs+1]
     scount=[dim1,dim2]

     allocate(gstart(npxy,nbdim))
     counter = 1
     do jj=1,npy
        stjj = (jj-1) * dim1  + 1
        do ii=1,npx
           stii = (ii-1) * dim2  + 1
           gstart(counter,:) = [stii,stjj]
           counter = counter + 1
        enddo
     enddo

  elseif ((trim(dimsorder)=='kji')) then

     write(*,*)'-> KJI dims order <-'

     ierr = nf90_open(trim(ncfiles(1)), nf90_nowrite, ncid)
     if (ierr /= nf90_noerr) stop "Error to open Netcdf file"

     !- dim 1 (z)
     ierr = nf90_inq_dimid(ncid, "dim1", Dim1ID)
     if (ierr /= nf90_noerr) stop "Error to read dim1 ID!"
     ierr = nf90_inquire_dimension(ncid, Dim1ID, len = dim1)
     dim1 = dim1
     gdim1 = dim1

     !- dim 2 (y)
     ierr = nf90_inq_dimid(ncid, "dim2", Dim2ID)
     if (ierr /= nf90_noerr) stop "Error to read dim2 ID!"
     ierr = nf90_inquire_dimension(ncid, Dim2ID, len = dim2)
     dim2 = dim2 - ys - ye
     gdim2 = dim2 * npy

     !- dim 3 (x)
     ierr = nf90_inq_dimid(ncid, "dim3", Dim3ID)
     if (ierr /= nf90_noerr) stop "Error to read dim3 ID!"
     ierr = nf90_inquire_dimension(ncid, Dim3ID, len = dim3)
     dim3 = dim3 - xs -xe
     gdim3 = dim3 * npx

     nbdim = 3

     allocate(sstart(nbdim))
     allocate(scount(nbdim))

     sstart=[1,ys+1,xs+1]
     scount=[dim1,dim2,dim3]

     allocate(gstart(npxy,nbdim))
     counter = 1

     do jj=1,npy
        stjj = (jj-1) * dim2  + 1
        do ii=1,npx
           stii = (ii-1) * dim3  + 1
           gstart(counter,:) = [1,stjj,stii]
           counter = counter + 1
        enddo
     enddo

  elseif (trim(dimsorder)=='dkji') then

     write(*,*)'-> DKJI dims order <-'

     ierr = nf90_open(trim(ncfiles(1)), nf90_nowrite, ncid)
     if (ierr /= nf90_noerr) stop "Error to open Netcdf file"

     !- dim 1 (z)
     ierr = nf90_inq_dimid(ncid, "dim1", Dim1ID)
     if (ierr /= nf90_noerr) stop "Error to read dim1 ID!"
     ierr = nf90_inquire_dimension(ncid, Dim1ID, len = dim1)
     dim1 = dim1
     gdim1 = dim1

     !- dim 2 (z)
     ierr = nf90_inq_dimid(ncid, "dim2", Dim2ID)
     if (ierr /= nf90_noerr) stop "Error to read dim2 ID!"
     ierr = nf90_inquire_dimension(ncid, Dim2ID, len = dim2)
     dim2 = dim2
     gdim2 = dim2

     !- dim 3 (y)
     ierr = nf90_inq_dimid(ncid, "dim3", Dim3ID)
     if (ierr /= nf90_noerr) stop "Error to read dim3 ID!"
     ierr = nf90_inquire_dimension(ncid, Dim3ID, len = dim3)
     dim3 = dim3 - ys - ye
     gdim3 = dim3 * npy

     !- dim 4 (x)
     ierr = nf90_inq_dimid(ncid, "dim4", Dim4ID)
     if (ierr /= nf90_noerr) stop "Error to read dim4 ID!"
     ierr = nf90_inquire_dimension(ncid, Dim4ID, len = dim4)
     dim4 = dim4 - xs -xe
     gdim4 = dim4 * npx

     nbdim = 4

     allocate(sstart(nbdim))
     allocate(scount(nbdim))

     sstart=[   1,   1,ys+1,xs+1]
     scount=[dim1,dim2,dim3,dim4]

     allocate(gstart(npxy,nbdim))
     counter = 1
     do jj=1,npy
        stjj = (jj-1) * dim3  + 1
        do ii=1,npx
           stii = (ii-1) * dim4  + 1
           gstart(counter,:) = [1,1,stjj,stii]
           counter = counter + 1
        enddo
     enddo

  else
     stop 'Error, dimension order is not valid.'
  endif

  ierr = nf90_inquire_variable(ncid, 1, name = vname)
  if (ierr /= nf90_noerr) stop

  ierr = nf90_close(ncid)

  !=============================!
  !- Create global NetCDF file -!
  !=============================!
  final_netcdf_file_name = 'global_netcdf_file.nc'
  ierr = nf90_create                       ( &
       path  = trim(final_netcdf_file_name), & ! NetCDF filename created.
       cmode = NF90_64BIT_OFFSET           , & ! 64_bits record, Large File
       ncid  = gncid                       )
  if (ierr /= nf90_noerr) stop

  !------------------!
  !- Set Dimensions -!
  !------------------!
  ierr=NF90_def_dim   ( &
       ncid  = gncid  , &
       name  = 'dim1' , &
       len   = gdim1  , &
       dimid = dim1id )
  if (ierr /= nf90_noerr) stop

  ierr=NF90_def_dim   ( &
       ncid  = gncid  , &
       name  = 'dim2' , &
       len   = gdim2  , &
       dimid = dim2id )
  if (ierr /= nf90_noerr) stop

  !- 2D -!
  if (nbdim == 2) then

     ierr=NF90_def_var              ( &
          ncid   = gncid            , &
          name   = trim(vname)      , &
          xtype  = NF90_DOUBLE      , &
          dimids = (/dim1id,dim2id/), &
          varid  = gvarid           )

     if (ierr /= nf90_noerr) stop

  !- 3D -!
  elseif (nbdim == 3) then

     ierr=NF90_def_dim  ( &
          ncid  = gncid , &
          name  = 'dim3', &
          len   = gdim3 , &
          dimid = dim3id)

     if (ierr /= nf90_noerr) stop

     ierr=NF90_def_var                     ( &
          ncid   = gncid                   , &
          name   = trim(vname)             , &
          xtype  = NF90_DOUBLE             , &
          dimids = (/dim1id,dim2id,dim3id/), &
          varid  = gvarid                  )

     if (ierr /= nf90_noerr) stop

  !- 4D -!
  elseif (nbdim == 4) then

     ierr=NF90_def_dim  ( &
          ncid  = gncid , &
          name  = 'dim3', &
          len   = gdim3 , &
          dimid = dim3id)

     if (ierr /= nf90_noerr) stop

     ierr=NF90_def_dim  ( &
          ncid  = gncid , &
          name  = 'dim4', &
          len   = gdim4 , &
          dimid = dim4id)

     if (ierr /= nf90_noerr) stop

     ierr=NF90_def_var                            ( &
          ncid   = gncid                          , &
          name   = trim(vname)                    , &
          xtype  = NF90_DOUBLE                    , &
          dimids = (/dim1id,dim2id,dim3id,dim4id/), &
          varid  = gvarid                         )

     if (ierr /= nf90_noerr) stop

  endif

  !------------------------------!
  !- Close the header defintion -!
  !------------------------------!
  ierr =  NF90_enddef(ncid = gncid)
  if (ierr /= nf90_noerr) stop

  !==========================================================================!
  !- Reading subdomain data and Writing this data in the global NetCDF file -!
  !==========================================================================!
  ! depend of the number of dimensions !

  if (nbdim == 2) then

     allocate(a2D(dim2,dim1))

     counter = 1

     do ii = 1,npxy

        print "('Reading rank:',I2)", counter-1

        ierr = nf90_open(trim(ncfiles(counter)), nf90_nowrite, ncid)
        if (ierr /= nf90_noerr) stop "Error to open Netcdf file"

        ierr = nf90_inq_varid(ncid, trim(vname), varID)
        if (ierr /= nf90_noerr) stop

        ierr = nf90_get_var(ncid, varID, a2D, start = sstart, count=scount)
        if (ierr /= nf90_noerr) stop "Error reading data"

        ierr = nf90_close(ncid)

        !print "('Global indices:',I2,' ',I2,' ',I2)",counter,stii,stjj

        ierr =  nf90_put_var(gncid, gvarID, a2D, start = gstart(ii,:), count = scount )
        if (ierr /= nf90_noerr) stop "Error writing data 2D in the global netcdf file"

        counter = counter + 1

     enddo

     deallocate(a2D)

  elseif (nbdim == 3) then

     allocate(a3D(dim3,dim2,dim1))

     counter = 1

     do ii = 1,npxy

        print "('Reading rank:',I2)", counter-1

        ierr = nf90_open(trim(ncfiles(counter)), nf90_nowrite, ncid)
        if (ierr /= nf90_noerr) stop "Error to open Netcdf file"

        ierr = nf90_inq_varid(ncid, trim(vname), varID)
        if (ierr /= nf90_noerr) stop

        ierr = nf90_get_var(ncid, varID, a3D, start = sstart, count=scount)
        if (ierr /= nf90_noerr) stop "Error reading data"

        ierr = nf90_close(ncid)

        ierr =  nf90_put_var(gncid, gvarID, a3D, start = gstart(ii,:), count=scount )
        if (ierr /= nf90_noerr) stop "Error writing 3D data in the global netcdf file"

        counter = counter + 1

     enddo

     deallocate(a3D)

  elseif (nbdim == 4) then

     allocate(a4D(dim4,dim3,dim2,dim1))

     counter = 1

     do ii=1, npxy

        print "('Reading rank:',I2)", counter-1

        ierr = nf90_open(trim(ncfiles(counter)), nf90_nowrite, ncid)
        if (ierr /= nf90_noerr) stop "Error to open Netcdf file"

        ierr = nf90_inq_varid(ncid, trim(vname), varID)
        if (ierr /= nf90_noerr) stop

        ierr = nf90_get_var(ncid, varID, a4D, start = sstart, count=scount)
        if (ierr /= nf90_noerr) stop "Error reading data"

        ierr = nf90_close(ncid)

        ierr =  nf90_put_var(gncid, gvarID, a4D, start = gstart(ii,:), count=scount )
        if (ierr /= nf90_noerr) stop "Error writing data 4D in the global netcdf file"

        counter = counter + 1

     enddo

     deallocate(a4D)

  endif

  !--------------------------------!
  !- Close the global netcdf file -!
  !--------------------------------!
  ierr = nf90_close(ncid = gncid)

  if (ierr /= nf90_noerr) stop

contains

  !---------------------------------!
  !- Transfort a string to integer -!
  !---------------------------------!
  elemental subroutine str2int(str,int,stat)

    character(len=*) , intent(in)  :: str
    integer(kind = 4), intent(out) :: int
    integer(kind = 4), intent(out) :: stat

    read(str,*,iostat=stat)  int

  end subroutine str2int

end program nhmg_ncjoin
