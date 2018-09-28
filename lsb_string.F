      program lsb_string

! Create nodename string as in supercomputer center environmental
! variable $LSB_HOSTS: a single string containing hostnames of compute
! nodes repeated as many times as number of CPU cores per node.

      implicit none
      character(len=2), parameter :: family='tn'
      integer, parameter :: first=1, last=32, ncores=8

c      character(len=2), parameter :: family='nd'
c      integer, parameter :: first=1, last=16, ncores=2

      character(len=16) name
      character(len=2000) string
      integer i, is,ie, node, lnm, lenstr

      is=1
      do node=first,last
        if (node<10) then
          write(name,'(A,I1)') family, node
        else
          write(name,'(A,I2)') family, node
        endif
        lnm=lenstr(name)
        do i=1,ncores
          ie=is+lnm-1
          string(is:ie)=name
          string(ie+1:ie+1)=' '
          is=ie+2
        enddo
      enddo

      write(*,'(3A)') 'string = ''', string(1:ie), ''''
      stop
      end

#include "lenstr.F"



