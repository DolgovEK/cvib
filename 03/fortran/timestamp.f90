subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 May 2014 by E. Dolgov
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  use iso_fortran_env

  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 ) time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5 ) zone
  real (kind = 8)     ::wall_current,cpu_current
  real (kind = 8),save::cpu_last,wall_last,cpu_start,wall_start
  integer(kind=4),save::n_start,s_start,h_start,d_start
  

  call date_and_time ( date, time, zone, values )
  call cpu_time(cpu_current)

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)
  
wall_current=((24.0D0*dble(d)+dble(h))*60.0D0+dble(n))*60.0D0+dble(s)+dble(mm)/1000.0D0 ! will be zeroed out in 1 month


  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if
  
!  write (*,*)
!  write (*,*) ' Current date and time:'
!  write ( *, '(2x,a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
!    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )
  if (wall_start/=0.0) then
     write (*,*)
     write (*,'(A35,F15.3)') 'Step wall clock time in seconds:  ', wall_current-wall_last
     write (*,'(A35,F15.3)') 'Step CPU clock time in seconds:   ', cpu_current-cpu_last
     write (*,'(A35,F15.3)') 'Step CPU utilization in percent:  ', (cpu_current-cpu_last)/(wall_current-wall_last)*100.0D0 
     write (*,*)
     

     d=int((wall_current-wall_start)/(60*60*24))
     h=int((wall_current-wall_start)/(60*60))-24*d
     n=int((wall_current-wall_start)/(60))-24*60*d-60*h
     s=(wall_current-wall_start)-24*60*60*d-60*60*h-60*n

     write (*,'(A30,I2,A6,I2,A7,I2,A9,I2,A9)') &
      'Total wall clock time passed: ',d,' days ',h,' hours ',n,' minutes ',s,' seconds '
     write (*,*)
     write (*,'(A35,F15.3)') 'Total wall clock time in seconds: ', wall_current-wall_start
     write (*,'(A35,F15.3)') 'Total CPU clock time in seconds:  ', cpu_current-cpu_start
     write (*,'(A35,F15.3)') 'Total CPU utilization in percent: ', (cpu_current-cpu_start)/(wall_current-wall_start)*100.0D0 
  else
     cpu_start=cpu_current
     wall_start=wall_current
  end if
  write (*,*)
  
  cpu_last=cpu_current
  wall_last=wall_current
  
  flush(output_unit);
  
  return
end
