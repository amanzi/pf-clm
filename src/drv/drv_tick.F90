!#include <misc.h>

subroutine drv_tick (drv)

!=========================================================================
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!=========================================================================
! DESCRIPTION:
!  Advance the time 1 timestep.
!
! REVISION HISTORY:
!  15 Oct 1999: Paul Houser; Initial Version
!=========================================================================
! $Id: drv_tick.F90,v 1.1.1.1 2006/02/14 23:05:52 kollet Exp $
!=========================================================================

  use clm_precision
  use drv_type_module      ! Land Driver 1-D variables
  implicit none

!=== Arguments ===========================================================

  type (drv_type) drv              

!=== Local Variables =====================================================

  integer days(13)
  data days/31,28,31,30,31,30,31,31,30,31,30,31,31/

!=== End Variable List ===================================================

  drv%pda=drv%da  !used to determine end-of-day restart writes

  drv%ss=drv%ss+drv%ts
!  write(*,*) "Seconds", drv%ss

  do while(drv%ss.gt.59)
     drv%ss=drv%ss-60
     drv%mn=drv%mn+1
  enddo

  do while(drv%mn.gt.59)
     drv%mn=drv%mn-60
     drv%hr=drv%hr+1
  enddo

  do while(drv%hr.ge.24)
     drv%hr=drv%hr-24
     drv%da=drv%da+1
  enddo
!  write(*,*)"hours:",drv%hr

  if ((mod(drv%yr,4).eq.0.AND.mod(drv%yr,100).ne.0)  & !correct for leap year
       .OR.(mod(drv%yr,400).eq.0))then                  !correct for Y2K
     days(2)=29                  
  else
     days(2)=28
  endif

  do while(drv%da.gt.days(drv%mo))
     drv%da=drv%da-days(drv%mo)
     drv%mo=drv%mo+1
  enddo

  do while(drv%mo.gt.12)
     drv%mo=drv%mo-12
     drv%yr=drv%yr+1
  enddo

!=== Update DRV current model TIME Variable
  call drv_date2time(drv%time,drv%doy,drv%day,drv%gmt, &
                     drv%yr,drv%mo,drv%da,drv%hr,drv%mn,drv%ss)
  write(*,24)'CLM-DRV Time: ',drv%mo,'/',drv%da,'/', &
       drv%yr,drv%hr,':',drv%mn,':',drv%ss

 24 format(a15,i2,a1,i2,a1,i4,1x,i2,a1,i2,a1,i2)
  return
end subroutine drv_tick


!=========================================================================
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!=========================================================================
! DESCRIPTION:
!  Determines time in years, based on year, month, day hour etc..
!   OR reverse (date2time).  
!
! REVISION HISTORY:
!  15 Oct 1999: Paul Houser; Initial Code
!=========================================================================

subroutine drv_date2time(time,     doy,     day,     gmt,       &
                         yr,       mo,      da,      hr,        &
                         mn,       ss)

  use clm_precision
  implicit none
  integer yr,mo,da,hr,mn,ss,yrdays,doy,days(13),k
  real*8 time
  real(r8) gmt,day
  data days /31,28,31,30,31,30,31,31,30,31,30,31,30/

  if ((mod(yr,4).eq.0.AND.mod(yr,100).ne.0)  & !correct for leap year
       .OR.(mod(yr,400).eq.0))then              !correct for Y2K
     days(2) = 29
     yrdays=366                  
  else
     yrdays=365
  endif
  doy=0
  do k=1,(mo-1)
     doy=doy+days(k)
  enddo
  doy=doy+da-1

  time=dfloat(yr)+((((((dfloat(ss)/60.d0)+dfloat(mn))/60.d0)+ &
       dfloat(hr))/24.d0)+dfloat(doy-1))/dfloat(yrdays)

  gmt=(((float(ss)/60.0) +float(mn)) /60.0)+float(hr)
  day=float(doy-1)+ &
       ((((float(ss)/60.0)+float(mn))/60.0)+float(hr))/24.0 

  return
end subroutine drv_date2time


subroutine drv_time2date(time,     doy,     day,     gmt,       &
                         yr,       mo,      da,      hr,        &
                         mn, ss)
  use clm_precision
  implicit none
  integer,intent(out) :: yr,mo,da,hr,mn,ss,doy
  real(r8),intent(out) :: gmt,day
  real(r8),intent(in) :: time
  
  ! locals
  real(r8) :: tmp
  integer :: yrdays, tmpdoy, days(13)
  data days /31,28,31,30,31,30,31,31,30,31,30,31,30/

  print*, "initing at time = ", time
  
  yr  = dint(time)
  tmp =     (time) 
  print*, " year = ", yr
  
  if ((mod(yr,4).eq.0.AND.mod(yr,100).ne.0)  & !correct for leap year
       .OR.(mod(yr,400).eq.0))then              !correct for Y2K
     yrdays=366                  
  else
     yrdays=365
  endif

  if (yrdays.eq.366) days(2)=29

  doy  = dint((tmp-yr)*dfloat(yrdays)) 
  tmp =      ((tmp-yr)*dfloat(yrdays)) 
  print*, " doy = ", doy

  hr  = dint((tmp-doy)*24.d0) 
  tmp =     ((tmp-doy)*24.d0) 
  print*, " hr = ", hr

  mn  = dint((tmp-hr)*60.d0) 
  tmp =     ((tmp-hr)*60.d0) 
  print*, " mn = ", mn

  ss  = dint((tmp-hr)*60.d0) 
  print*, " ss = ", ss

  tmpdoy = doy
  mo=0
  do while (tmpdoy.ge.0)
     mo=mo+1
     tmpdoy=tmpdoy-days(mo)
  enddo
  da=tmpdoy+days(mo)+1
  print*, " mo = ", mo
  print*, " da = ", da

  gmt=(((float(ss)/60.0)+float(mn))/60.0)+float(hr)
  day=float(doy)+ &
       ((((float(ss)/60.0)+float(mn))/60.0)+float(hr))/24.0 

  print*, " gmt = ", gmt
  print*, " day = ", day
  write(*,24)'CLM-DRV Time: ',mo,'/',da,'/', &
       yr,hr,':',mn,':',ss
 24 format(a15,i2,a1,i2,a1,i4,1x,i2,a1,i2,a1,i2)
  
  return
end subroutine drv_time2date











