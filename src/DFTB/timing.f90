module timing
!real :: timetab(20), timetab_start(20), timetab_stop(20)
real    (kind=8) :: timetab(20)
integer (kind=8) :: timetab_start(20), timetab_stop(20)
integer (kind=8) :: counts, clock_rate,  count_max
public :: timetab, timetab_start,  timetab_stop


interface time_stamp
  module procedure  time_stamp_only
  module procedure  time_stamp_mesg
end interface time_stamp


interface time_start 
  module procedure  start_time
end interface time_start 


interface time_stop
  module procedure  stop_time
end interface time_stop



!---------------------
contains

  subroutine start_time(item)
     integer item
   !  real, save ::  timetab_start(20)
   !  real, save ::  timetab_stop(20)
   !  real ::  timetab_start(20)
   !  real ::  timetab_stop(20)
     !real       :: t0
     !call cpu_time(t0)
     call system_clock (counts, clock_rate, count_max )
   if (item.le.20) timetab_start(item)= counts
  end subroutine


  subroutine stop_time(item)
    integer item
  !  real, save ::  timetab_start(20)
  !  real, save ::  timetab_stop(20)
  !  real ::  timetab_start(20)
  !  real ::  timetab_stop(20)
    !real       :: t0
    !call cpu_time(t0)
    call system_clock (counts, clock_rate, count_max )
    if (item.le.20)  then
        !timetab_stop(item)= t0
        timetab_stop(item)= counts
        timetab(item) =  timetab(item) + real(timetab_stop(item)-timetab_start(item))/real(clock_rate)
    endif
  end subroutine


  subroutine  time_stamp_only( )
      !character* ,intent(in)   :: mesg
      character(8)  :: date
      character(10) :: time
      character(5)  :: zone
      integer,dimension(8) :: values
    ! using keyword arguments
    ! VALUES is INTENT(OUT) and provides the following:
    !       VALUE(1):       The year
    !       VALUE(2):       The month
    !       VALUE(3):       The day of the month
    !       VALUE(4):       Time difference with UTC in minutes
    !       VALUE(5):       The hour of the day
    !       VALUE(6):       The minutes of the hour
    !       VALUE(7):       The seconds of the minute
    !       VALUE(8):       The milliseconds of the second 
    ! DATE    (Optional) The type shall be CHARACTER(LEN=8) or larger, and of default kind.
    ! TIME    (Optional) The type shall be CHARACTER(LEN=10) or larger, and of default kind.
    ! ZONE    (Optional) The type shall be CHARACTER(LEN=5) or larger, and of default kind.
    ! VALUES  (Optional) The type shall be INTEGER(8). 
    !
    !call date_and_time(date,time,zone,values)
    !call date_and_time(DATE=date,ZONE=zone)
    !call date_and_time(TIME=time)
    !call date_and_time(VALUES=values)
    !print '(a,2x,a,2x,a)', date, time, zone
    !print '(8i5)', values


    call date_and_time(date,time,zone,values)
    write(*,1000)values
    !1000 format ('DATE: 'I4,'-',I2,'-',I2,' (UTC:',I4,')   TIME: ',I2,':',I2.2,':',I2.2,'.',I3)      
    1000 format ('DATE: 'I4,'-',I2,'-',I2,' (UTC:',I4,')   TIME: ',I2,':',I0.2,':',I0.2,'.',I3)
  end subroutine


  subroutine  time_stamp_mesg(mesg )
      character(len=*) ,intent(in)   :: mesg
      character(8)  :: date
      character(10) :: time
      character(5)  :: zone
      integer,dimension(8) :: values
    ! using keyword arguments
    ! VALUES is INTENT(OUT) and provides the following:
    !       VALUE(1):       The year
    !       VALUE(2):       The month
    !       VALUE(3):       The day of the month
    !       VALUE(4):       Time difference with UTC in minutes
    !       VALUE(5):       The hour of the day
    !       VALUE(6):       The minutes of the hour
    !       VALUE(7):       The seconds of the minute
    !       VALUE(8):       The milliseconds of the second 
    ! DATE    (Optional) The type shall be CHARACTER(LEN=8) or larger, and of default kind.
    ! TIME    (Optional) The type shall be CHARACTER(LEN=10) or larger, and of default kind.
    ! ZONE    (Optional) The type shall be CHARACTER(LEN=5) or larger, and of default kind.
    ! VALUES  (Optional) The type shall be INTEGER(8). 
    !
    !call date_and_time(date,time,zone,values)
    !call date_and_time(DATE=date,ZONE=zone)
    !call date_and_time(TIME=time)
    !call date_and_time(VALUES=values)
    !print '(a,2x,a,2x,a)', date, time, zone
    !print '(8i5)', values

    call date_and_time(date,time,zone,values)
    write(*,1000)trim(mesg),values
    !1000 format ('DATE: 'I4,'-',I2,'-',I2,' (UTC:',I4,')   TIME: ',I2,':',I2.2,':',I2.2,'.',I3)      
    1000 format (A,' DATE: 'I4,'-',I2,'-',I2,' (UTC:',I4,')   TIME: ',I2,':',I0.2,':',I0.2,'.',I3)
  end subroutine




end module timing 
