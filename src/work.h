! This is "work.h": declaration of utility work array.
!
      real work(GLOBAL_2D_ARRAY,0:N)
CSDISTRIBUTE_RESHAPE work(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /work3a/ work
 
