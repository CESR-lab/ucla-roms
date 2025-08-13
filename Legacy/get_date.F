#include "cppdefs.opt"

      subroutine get_date (date_str)
!
! Get todays date, day of the week and time called (default month
! and weekday are December & Saturday respectively).
! By default this subroutine uses SUN intrinsic date routine.
!
! Output:  date_str  Concatenated string for the day of the
!                    week, date (month,day,year), and time
!                    (12hr clock) of day (hour:minute:sec)
!
! Calls: day_code (vax/cray only, code is provided in this file);
!        none, except cray/vax;
!
! NOTE: This file is shared with plottin package.
!
! Copyright (c) 1996 Rutgers University
!
      implicit none
      character*(*) date_str
      integer year,hour, minute,sec, half, iday,imon, dstat,tstat,
     &        nday, lmonth(12), lday(31), len1, len2, len3, lenstr
      character*3  ampm(0:1)
      character*9 day(0:6), month(12)
      data lmonth/7,8,5,5,3,4,4,6,9,7,8,8/ ampm/' AM',' PM'/
     &    lday /9*1,22*2/ day /'Sunday   ', 'Monday   ', 'Tuesday  ',
     &            'Wednesday', 'Thursday ', 'Friday   ', 'Saturday '/
     &      month/'January  ', 'February ', 'March    ', 'April    ',
     &            'May      ', 'June     ', 'July     ', 'August   ',
     &            'September', 'October  ', 'November ', 'December '/
      character*11 ctime*11, today*18, fmt*20, wkday*44

#if defined cray
      integer century
      parameter (century=1900)
      character*8  tstring
#elif defined sun || defined sgi || defined DECALPHA
      character*3  day3, mon
      character*28 fdate, tmpday
#elif defined XLF
      character*3  day3, mon
      character*28 tmpday
#endif
!
! Get weekday, date and time in short format, then extract this
! information.
!
#if defined vax
      dstat=0
      call idate (imon,nday,year)
      year=year+century
      call time (tstring)
      read(tstring,'(i2,1x,i2,1x,i2)',iostat=tstat) hour, minute, sec
      if(tstat/=0) ctime=tstring
#elif defined cray
      write(tstring,'(a8)') date()
      read(tstring,'(i2,1x,i2,1x,i2)',iostat=dstat) imon, nday, year
      year=year+century
      if(dstat/=0) then
        wkday=tstring
        today=' '
      endif
      write(tstring,'(a8)') clock()
      read(tstring,'(i2,1x,i2,1x,i2)',iostat=tstat) hour, minute, sec
      if(tstat/=0) ctime=tstring
#elif defined sun || defined sgi || defined DECALPHA
      tmpday=fdate()
      read(tmpday,'(a3,1x,a3,1x,i2)',iostat=dstat) day3, mon, nday
      read(tmpday,'(11x,i2,1x,i2,1x,i2)',iostat=tstat)hour,minute,sec
      tstat=max(abs(dstat),abs(tstat))
      read(tmpday,'(20x,i4)',iostat=dstat) year
      if(dstat/=0 .or. tstat/=0) then
        dstat=1
        tstat=1
        wkday=tmpday
        today=' '
        ctime=' '
      endif
#elif defined XLF
      call fdate_(tmpday)
      read(tmpday,'(a3,1x,a3,1x,i2)',iostat=dstat) day3, mon, nday
      read(tmpday,'(11x,i2,1x,i2,1x,i2)',iostat=tstat)hour,minute,sec
      tstat=max(abs(dstat),abs(tstat))
      read(tmpday,'(20x,i4)',iostat=dstat) year
      if((dstat/=0).or.(tstat/=0)) then
        dstat=1
        tstat=1
        wkday=tmpday
        today=' '
        ctime=' '
      endif
#else
      hour=0
      minute=0
      sec=0
      nday=1
c      mon='Jan'
      dstat=1
      tstat=1
      wkday=' '
      today=' '
      ctime=' '
#endif
      if (tstat==0) then            ! Convert from 24 hour clock
        half=hour/12                  ! clock to 12 hour AM/PM clock.
        hour=hour-half*12
        if (hour==0) hour=12
        if (half==2) half=0
      endif
      if (dstat==0) then
#if defined vax || defined cray
                                             ! Get index for
        call day_code (imon,nday,year,iday)  ! the day of the week.

#elif defined sun || defined sgi || defined XLF || defined DECALPHA
                          ! Loop to find full day name by comparing
                          ! DAY3 with the first 3 letters of day.
        iday=0
        do while ((day3/=day(iday)(1:3)).and.(iday<6))
          iday=iday+1
        enddo             ! Loop to find full month name by comparing
                          ! MON with first 3 letters of month.
        imon=1
        do while ((mon/=month(imon)(1:3)).and.(imon<12))
          imon=imon+1
        enddo
#endif
!
! Construct date, time and day of the week output string.
!
        write(fmt,10) lmonth(imon), lday(nday)
  10    format('(a',i1,',1x,i',i1,',1h,,1x,i4)')
        write(today,fmt) month(imon),nday,year
        wkday=day(iday)
      endif
      if(tstat==0) then
        write(ctime,20) hour, minute, sec, ampm(half)
  20    format(i2,':',i2.2,':',i2.2,a3)
      endif

      len1=lenstr(wkday)                  ! Concatenate date string.
      len2=lenstr(today)
      len3=lenstr(ctime)
      date_str=wkday(1:len1)
      if (len2>0) then
        len1=lenstr(date_str)
        date_str=date_str(1:len1)/ /' - '/ /today(1:len2)
      endif
      if (len3>0) then
        len1=lenstr(date_str)
        date_str=date_str(1:len1)/ /' - '/ /ctime(1:len3)
      endif
      return
      end


#if defined vax || defined cray

      subroutine day_code (month,day,year,code)
!
! Compute a code for the day of the week, given the date. This code
! is good for date after January 1, 1752 AD, the year the Gregorian
! calander was adopted in Britian and the American colonies.
!
! Input:
!     month     The month, 1=January, 2=February, ... (integer).
!     day       The day of the month (integer).
!     year      The year, including the century (integer).
!
! Output:
!     code      A code for the corresponding day of the week
!               (integer):
!               code = 0  =>  Sunday
!               code = 1  =>  Monday
!               code = 2  =>  Tuesday
!               code = 3  =>  Wednesday
!               code = 4  =>  Thursday
!               code = 5  =>  Friday
!               code = 6  =>  Saturday
!
      implicit none
      logical leap_year
      integer month, day, year, code, base_cen,   base_qcen, i,
     &        base_qyear, base_year,  bym1_dec31, feb_end,   leap,
     &        no_day, no_yr,  nqy,    nyc, nyqc,  month_day(12)
      parameter (base_cen=1700,   base_qcen=1600, base_qyear=1748,
     &           base_year=1752,  bym1_dec31=5,   feb_end=59)
      data month_day /31,28,31,30,31,30,31,31,30,31,30,31/
!
! Compute the number of years since the base year, the number of
! years since the beginning of the base century and the number of
! years since the beginning of the base 400 year.
!
      no_yr=year-base_year
      nqy=year-base_qyear
      nyc=year-base_cen
      nyqc=year-base_qcen
!
! Compute the number of leapdays in that time.
! Determine if this is a leap year.
!
      leap=nqy/4-nyc/100+nyqc/400
      leap_year=(mod(nqy,4)==0 .and. mod(nyc,100)/=0)
     &                            .or. mod(nyqc,400)==0
!
! Compute the number of days this year.  The leap year corrections
! are:   Jan. 1 - Feb. 28   Have not had the leap day counted above.
!        Feb.29             Counting leap day twice.
!
      no_day=day
      do i=1,month-1
        no_day=no_day+month_day(i)
      enddo
      if (leap_year .and. no_day<=feb_end)  no_day=no_day-1
      if (leap_year .and. month==2 .and. day==29) no_day=no_day-1
!
! Compute the total number of days since Jan. 1 of the base year,
! exclusive of the 364 day per year which represent an even 52
! weeks.  Actually, only need to do the addition mod 7.
!
      no_day=mod(no_day,7)+mod(leap,7)+mod(no_yr,7)+bym1_dec31
!
! Get the day of the week code.
!
      code=mod(no_day,7)
      return
      end
#endif

