! Character strings to hold names of activated cpp-switches and list
! of source-code filenames to keep track which ones are activated/used.
! Basically the names are pasted together consecutively into long
! strings and written as global attributes into output netCDF files.

      integer, parameter :: max_opt_size=2048
      character*(max_opt_size) cpps, srcs, kwds
      common /strings/ cpps, srcs, kwds
 
