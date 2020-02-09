      implicit none      ! This code just prints
      integer i          ! ASCII characters agains 
      character A*1      ! their numbers
      do i=0,255
        A=char(i)
        write(*,'(I4,1x,3A1)') i, '''',  A, ''''
      enddo
      stop
      end
 
 
