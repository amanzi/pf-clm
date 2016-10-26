module drv_type_module 
!=========================================================================
!
!  CLMCLMCLMCLMCLMCLMCLMCLMCL  A community developed and sponsored, freely   
!  L                        M  available land surface process model.  
!  M --COMMON LAND MODEL--  C  
!  C                        L  CLM WEB INFO: http://clm.gsfc.nasa.gov
!  LMCLMCLMCLMCLMCLMCLMCLMCLM  CLM ListServ/Mailing List: 
!
!=========================================================================
! drv_module.f90: 
!
! DESCRIPTION:
!  Module for 1-D land model driver variable specification.
!
! REVISION HISTORY:
!  15 Jan 2000: Paul Houser; Initial code
!   3 Mar 2000: Jon Radakovich; Revision for diagnostic output
!=========================================================================     
! $Id: drv_module.F90,v 1.1.1.1 2006/02/14 23:05:52 kollet Exp $
!=========================================================================

  use clm_precision
  use clm_infnan
  implicit none

  public :: drv_init, drv_destroy
  
  type, public :: drv_type

!=== Driver User-Defined Parameters ======================================
     real(r8) ::               &
          mina,             & !Min grid area for tile (%)
          udef                !Undefined value

     character*40 ::       &
          vegtf,            & !Vegetation Tile Specification File                   
          vegpf,            & !Vegetation Type Parameter Values
          poutf1d,          & !CLM 1D Parameter Output File
          metf1d,           & !Meterologic input file
          outf1d,           & !CLM output file
          rstf                !CLM active restart file

!=== CLM Parameters ======================================================
     integer :: nch          !actual number of tiles

!=== Driver Parameters ====================================================
     integer :: nc           !Number of Columns in Grid
     integer :: nr           !Number of Rows in Grid
     integer :: nt           !Number of Vegetation Types 
     integer :: startcode    !0=restart date, 1=card date
     integer :: sss          !Starting Second 
     integer :: sdoy         !Starting Day of Year 
     integer :: smn          !Starting Minute 
     integer :: shr          !Starting Hour 
     integer :: sda          !Starting Day 
     integer :: smo          !Starting Month 
     integer :: syr          !Starting Year  
     integer :: ess          !Ending Second
     integer :: emn          !Ending Minute
     integer :: edoy         !Ending Day of Year
     integer :: ehr          !Ending Hour
     integer :: eda          !Ending Day
     integer :: emo          !Ending Month
     integer :: eyr          !Ending Year
     integer :: ts           !Timestep (seconds) 
     integer :: ts_old
     real(r8) :: writeintc   !CLM Output Interval (hours)
     integer :: maxt         !Maximum tiles per grid  

!=== Timing Variables ==========
     real*8  :: time                  !CLM Current Model Time in Years
     real*8  :: etime                 !CLM End Time in Years
     integer :: pda                   !CLM Previous Timestep Day
     integer :: doy,yr,mo,da,hr,mn,ss !CLM Current Model Timing Variables   
     integer :: endtime               !CLM Stop (0=continue time looping)
     real(r8):: day,gmt,eday,egmt,sgmt

!=== Arguments ==========================================================
     real(r8) :: ctime                !CLM Restart Time 
     integer :: cyr,cmo,cda       !Restart Model Timing Variables
     integer :: chr,cmn,css       !Restart Model Timing Variables

!=== Initial CLM conditions =============================================
     real(r8) :: t_ini                !Initial temperature [K] 
     real(r8) :: h2osno_ini              !Initial snow cover, water equivalent [mm] 
     real(r8) :: sw_ini               !Initial average soil volumetric water&ice content [m3/m3] 

!=== CLM diagnostic parameters ==========================================
     integer :: surfind      !Number of surface diagnostic variables
     integer :: soilind      !Number of soil layer diagnostic variables
     integer :: snowind      !Number of snow layer diagnostic variables

     integer :: vclass            !Vegetation Classification Scheme (1=UMD,2=IGBP,etc.) NOT the index 
     integer :: clm_ic            !CLM Initial Condition Source
  
!=== CLM.PF varibales
     integer  :: sat_flag         ! 0: enough storage in the domain; 1: too little storage in the domain, full saturation
!     real(r8) :: dx,dy,dz         
     real(r8) :: begwatb, endwatb ! beg and end water balance over domain      
     
!=== End Variable List ===================================================
  end type drv_type

contains

  subroutine drv_init(drv)
    type(drv_type) drv

    drv%mina = 0d0
    drv%udef = NaN

    drv%vegtf = ''
    drv%vegpf = ''
    drv%poutf1d = ''
    drv%metf1d = ''
    drv%outf1d = ''
    drv%rstf = ''

    drv%nch = -1
    drv%nc = -1
    drv%nr = -1
    drv%nt = 18 ! default 18 IGBP land cover classes
    drv%startcode = 0
    drv%sss = -1
    drv%sdoy = -1
    drv%smn = -1
    drv%shr = -1
    drv%sda = -1
    drv%smo = -1
    drv%syr = -1
    drv%ess = -1
    drv%edoy = -1
    drv%emn = -1
    drv%ehr = -1
    drv%eda = -1
    drv%emo = -1
    drv%eyr = -1
    drv%ss = -1
    drv%doy = -1
    drv%mn = -1
    drv%hr = -1
    drv%da = -1
    drv%mo = -1
    drv%yr = -1

    drv%ts = -1
    drv%ts_old = -1

    drv%writeintc = -1
    drv%maxt = -1

    drv%time = NaN
    drv%etime = NaN
    drv%pda = -1
    drv%endtime = -1
    drv%day = NaN
    drv%gmt = NaN
    drv%eday = NaN
    drv%egmt = NaN
    drv%sgmt = NaN

    drv%ctime = NaN
    drv%cyr = -1
    drv%cmo = -1
    drv%cda = -1
    drv%chr = -1
    drv%cmn = -1
    drv%css = -1

    drv%t_ini = NaN
    drv%h2osno_ini = NaN
    drv%sw_ini = NaN

    drv%surfind = -1
    drv%soilind = -1    
    drv%snowind = -1
    drv%vclass = -1
    drv%clm_ic = 0
    drv%sat_flag = -1

    drv%begwatb = 0.d0
    drv%endwatb = 0.d0
    return
  end subroutine drv_init

  subroutine drv_destroy(drv)
    type(drv_type) drv
    return
  end subroutine drv_destroy
  
  
end module drv_type_module
