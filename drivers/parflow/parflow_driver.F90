subroutine clm_lsm(pressure,saturation,evap_trans,topo,porosity,pf_dz_mult,istep_pf,dt,time,           &
start_time,pdx,pdy,pdz,ix,iy,nx,ny,nz,nx_f,ny_f,nz_f,nz_rz,ip,npp,npq,npr,gnx,gny,rank,sw_pf,lw_pf,    &
prcp_pf,tas_pf,u_pf,v_pf,patm_pf,qatm_pf,lai_pf,sai_pf,z0m_pf,displa_pf,                               &
eflx_lh_pf,eflx_lwrad_pf,eflx_sh_pf,eflx_grnd_pf,                                                     &
qflx_tot_pf,qflx_grnd_pf,qflx_soi_pf,qflx_eveg_pf,qflx_tveg_pf,qflx_in_pf,swe_pf,t_g_pf,               &
t_soi_pf,clm_dump_interval,clm_1d_out,clm_forc_veg,clm_output_dir,clm_output_dir_length,clm_bin_output_dir,         &
write_CLM_binary,beta_typepf,veg_water_stress_typepf,wilting_pointpf,field_capacitypf,                 &
res_satpf,irr_typepf, irr_cyclepf, irr_ratepf, irr_startpf, irr_stoppf, irr_thresholdpf,               &
qirr_pf,qirr_inst_pf,irr_flag_pf,irr_thresholdtypepf,soi_z,clm_next,clm_write_logs,                    &
clm_last_rst,clm_daily_rst)

  use clm_precision
  use clm_io_config, only : MAX_FILENAME_LENGTH
  use clm_type_module
  use clm_host
  use clm1d_varpar, only : nlevsoi
  use grid_type_module
  use tile_type_module
  use clm1d_type_module
  implicit none

  !=== interface variables =====================================================
  integer,intent(in) :: nx,ny,nz,nx_f,ny_f,nz_f,nz_rz
  integer,intent(in) :: soi_z                    ! Specify layer shold be used for reference temperature
  real(r8),intent(in) :: pressure((nx+2)*(ny+2)*(nz+2))     ! pressure head, from parflow on grid w/ ghost nodes for current proc
  real(r8),intent(in) :: saturation((nx+2)*(ny+2)*(nz+2))   ! saturation from parflow, on grid w/ ghost nodes for current proc
  real(r8),intent(out) :: evap_trans((nx+2)*(ny+2)*(nz+2))   ! ET flux from CLM to ParFlow on grid w/ ghost nodes for current proc
  real(r8),intent(in) :: topo((nx+2)*(ny+2)*(nz+2))         ! mask from ParFlow 0 for inactive, 1 for active, on grid w/ ghost nodes for current proc
  real(r8),intent(in) :: porosity((nx+2)*(ny+2)*(nz+2))     ! porosity from ParFlow, on grid w/ ghost nodes for current proc
  real(r8),intent(in) :: pf_dz_mult((nx+2)*(ny+2)*(nz+2))   ! dz multiplier from ParFlow on PF grid w/ ghost nodes for current proc
  real(r8),intent(in) :: dt                                 ! parflow dt in parflow time units not CLM time units
  real(r8),intent(in) :: time                               ! parflow time in parflow units
  real(r8),intent(in) :: start_time                         ! starting time in parflow units
  real(r8),intent(in) :: pdx,pdy,pdz                        ! parflow DX, DY and DZ in parflow units
  integer,intent(in)  :: istep_pf                           ! istep, now passed from PF
  integer,intent(in)   :: ix                                 ! parflow ix, starting point for local grid on global grid
  integer,intent(in)   :: iy                                 ! parflow iy, starting point for local grid on global grid
  integer,intent(in)   :: ip                               
  integer,intent(in)   :: npp,npq,npr                        ! number of processors in x,y,z
  integer,intent(in)   :: gnx, gny                           ! global grid, nx and ny
  integer,intent(in)   :: rank                               ! processor rank, from ParFlow

  integer,intent(in)  :: clm_next                           ! NBE: Passing flag to sync outputs
  integer,intent(in)  :: clm_write_logs                     ! NBE: Enable/disable writing of the log files
  integer,intent(in)  :: clm_last_rst                       ! NBE: Write all the CLM restart files or just the last one
  integer,intent(in)  :: clm_daily_rst                      ! NBE: Write daily restart files or hourly

  ! surface fluxes & forcings
  real(r8),intent(out)  :: eflx_lh_pf((nx+2)*(ny+2)*3)        ! e_flux   (lh)    output var to send to ParFlow, on grid w/ ghost nodes for current proc but nz=1 (2D)
  real(r8),intent(out) :: eflx_lwrad_pf((nx+2)*(ny+2)*3)     ! e_flux   (lw)    output var to send to ParFlow, on grid w/ ghost nodes for current proc but nz=1 (2D)
  real(r8),intent(out) :: eflx_sh_pf((nx+2)*(ny+2)*3)        ! e_flux   (sens)  output var to send to ParFlow, on grid w/ ghost nodes for current proc but nz=1 (2D)
  real(r8),intent(out) :: eflx_grnd_pf((nx+2)*(ny+2)*3)      ! e_flux   (grnd)  output var to send to ParFlow, on grid w/ ghost nodes for current proc but nz=1 (2D)
  real(r8),intent(out) :: qflx_tot_pf((nx+2)*(ny+2)*3)       ! h2o_flux (total) output var to send to ParFlow, on grid w/ ghost nodes for current proc but nz=1 (2D)
  real(r8),intent(out) :: qflx_grnd_pf((nx+2)*(ny+2)*3)      ! h2o_flux (grnd)  output var to send to ParFlow, on grid w/ ghost nodes for current proc but nz=1 (2D)
  real(r8),intent(out) :: qflx_soi_pf((nx+2)*(ny+2)*3)       ! h2o_flux (soil)  output var to send to ParFlow, on grid w/ ghost nodes for current proc but nz=1 (2D)
  real(r8),intent(out) :: qflx_eveg_pf((nx+2)*(ny+2)*3)      ! h2o_flux (veg-e) output var to send to ParFlow, on grid w/ ghost nodes for current proc but nz=1 (2D)
  real(r8),intent(out) :: qflx_tveg_pf((nx+2)*(ny+2)*3)      ! h2o_flux (veg-t) output var to send to ParFlow, on grid w/ ghost nodes for current proc but nz=1 (2D)
  real(r8),intent(out) :: qflx_in_pf((nx+2)*(ny+2)*3)        ! h2o_flux (infil) output var to send to ParFlow, on grid w/ ghost nodes for current proc but nz=1 (2D)
  real(r8),intent(out) :: swe_pf((nx+2)*(ny+2)*3)            ! swe              output var to send to ParFlow, on grid w/ ghost nodes for current proc but nz=1 (2D)
  real(r8),intent(out) :: t_g_pf((nx+2)*(ny+2)*3)            ! t_grnd           output var to send to ParFlow, on grid w/ ghost nodes for current proc but nz=1 (2D)
  real(r8),intent(out) :: t_soi_pf((nx+2)*(ny+2)*(nlevsoi+2))!tsoil             output var to send to ParFlow, on grid w/ ghost nodes for current proc, but nz=10 (3D)
  real(r8),intent(in) :: sw_pf((nx+2)*(ny+2)*3)             ! SW rad, passed from PF
  real(r8),intent(in) :: lw_pf((nx+2)*(ny+2)*3)             ! LW rad, passed from PF
  real(r8),intent(in) :: prcp_pf((nx+2)*(ny+2)*3)           ! Precip, passed from PF
  real(r8),intent(in) :: tas_pf((nx+2)*(ny+2)*3)            ! Air temp, passed from PF
  real(r8),intent(in) :: u_pf((nx+2)*(ny+2)*3)              ! u-wind, passed from PF
  real(r8),intent(in) :: v_pf((nx+2)*(ny+2)*3)              ! v-wind, passed from PF
  real(r8),intent(in) :: patm_pf((nx+2)*(ny+2)*3)           ! air pressure, passed from PF
  real(r8),intent(in) :: qatm_pf((nx+2)*(ny+2)*3)           ! air specific humidity, passed from PF
  real(r8),intent(in) :: lai_pf((nx+2)*(ny+2)*3)            ! BH: lai, passed from PF
  real(r8),intent(in) :: sai_pf((nx+2)*(ny+2)*3)            ! BH: sai, passed from PF
  real(r8),intent(in) :: z0m_pf((nx+2)*(ny+2)*3)            ! BH: z0m, passed from PF
  real(r8),intent(in) :: displa_pf((nx+2)*(ny+2)*3)         ! BH: displacement height, passed from PF
  real(r8),intent(in) :: irr_flag_pf((nx+2)*(ny+2)*3)       ! irrigation flag for deficit-based scheduling -- 1 = irrigate, 0 = no-irrigate
  real(r8),intent(inout) :: qirr_pf((nx+2)*(ny+2)*3)           ! irrigation applied above ground -- spray or drip (2D)
  real(r8),intent(inout) :: qirr_inst_pf((nx+2)*(ny+2)*(nlevsoi+2))! irrigation applied below ground -- 'instant' (3D)

  ! output keys
  real(r8),intent(in) :: clm_dump_interval                  ! dump inteval for CLM output, passed from PF, always in interval of CLM timestep, not time   !FIXME bug -- should be integer --etc
  integer,intent(in)  :: clm_1d_out                         ! whether to dump 1d output 0=no, 1=yes
  integer,intent(in)  :: clm_forc_veg                       ! BH: whether vegetation (LAI, SAI, z0m, displa) is being forced 0=no, 1=yes
  integer,intent(in)  :: clm_output_dir_length              ! for output directory
  integer,intent(in)  :: clm_bin_output_dir                 ! output directory
  integer,intent(in)  :: write_CLM_binary                   ! whether to write CLM output as binary 
  character(LEN=clm_output_dir_length),intent(in) :: clm_output_dir ! output dir location

  ! ET keys
  integer,intent(in)  :: beta_typepf                        ! beta formulation for bare soil Evap 0=none, 1=linear, 2=cos
  integer,intent(in)  :: veg_water_stress_typepf            ! veg transpiration water stress formulation 0=none, 1=press, 2=sm
  real(r8),intent(in) :: wilting_pointpf                    ! wilting point in m if press-type, in saturation if soil moisture type
  real(r8),intent(in) :: field_capacitypf                   ! field capacity for water stress same as units above
  real(r8),intent(in) :: res_satpf                          ! residual saturation from ParFlow

  ! irrigation keys
  integer,intent(in)  :: irr_typepf                         ! irrigation type flag (0=none,1=spray,2=drip,3=instant)
  integer,intent(in)  :: irr_cyclepf                        ! irrigation cycle flag (0=constant,1=deficit)
  real(r8),intent(in) :: irr_ratepf                         ! irrigation application rate for spray and drip [mm/s]
  real(r8),intent(in) :: irr_startpf                        ! irrigation daily start time for constant cycle
  real(r8),intent(in) :: irr_stoppf                         ! irrigation daily stop tie for constant cycle
  real(r8),intent(in) :: irr_thresholdpf                    ! irrigation threshold criteria for deficit cycle (units of soil moisture content)
  integer,intent(in)  :: irr_thresholdtypepf                ! irrigation threshold criteria type -- top layer, bottom layer, column avg

  
  !=== local variables =====================================================
  type(clm_type),save :: clm            ! note use of SAVE to stay consistent with parflow "1 point of entry"
  type(host_type),save :: host

  character(len=MAX_FILENAME_LENGTH) :: output_dir
  real(r8) :: junk2d((nx+2)*(ny+2)*3)            ! placeholder junk data on a 2d surface
  real(r8) :: junk3d((nx+2)*(ny+2)*(nz+2))       ! placeholder junk data on a 3d surface
  integer :: d_stp                              ! NBE: Dummy for CLM restart

  real(r8) :: latlon((nx+2)*(ny+2)*3,2)    ! latitude,longitude [degrees]
  real(r8) :: sand((nx+2)*(ny+2)*(nz+2))          ! percent sand FIXME: 0-1 or 0-100? --etc
  real(r8) :: clay((nx+2)*(ny+2)*(nz+2))          ! percent clay FIXME: 0-1 or 0-100? --etc
  integer :: color_index((nx+2)*(ny+2)*3)  ! color index FIXME: document! --etc
  real(r8) :: fractional_ground((nx+2)*(ny+2)*3, 18) ! fraction of land surface of type t
  
  !=== begin code ============================================================

  if (time == start_time) then
  
     !--- initialize the clm master object
     call clm_init(clm, rank, nx, ny)

     !--- open log files, set up io options
     write(output_dir,*) trim(adjustl(clm_output_dir))
     call io_open(clm%io, output_dir, rank, clm_write_logs)
     clm%io%restart_last = clm_last_rst
     clm%io%restart_daily = clm_daily_rst
     clm%io%dump_interval = clm_dump_interval 
     clm%io%dump_current = 0
     clm%io%output_1d = clm_1d_out
     clm%io%write_binaries = write_CLM_binary
  
     !--- initialize the host object, pushing grid info into it
     call host_info_init(host, nx, ny, nz, topo)
     if (clm%io%log /= 0) call host_write_to_log(host, clm%io%log)
     if (clm%io%ranked_log /= 0) call host_write_to_log(host, clm%io%ranked_log)

     !--- initalize the clm master object's grid
     clm%grid => grid_create_2d(nx, ny, clm%drv%nt)

     !--- read clm input files, data
     ! FIXME -- can we remove this and use a setter? --etc
     ! FIXME -- fix io to not pass clm_write_logs, instead iounit --etc
     call drv_readclmin(clm%drv, clm%grid, rank, clm_write_logs)

     if (clm%drv%startcode == 0) stop
     if (clm%drv%clm_ic == 0) stop

     !--- initialize the clm master objects tiles, 1d columns
     clm%tile => tile_create_n(clm%ntiles)
     clm%clm => clm1d_create_n(clm%ntiles, clm%drv%surfind, clm%drv%soilind, clm%drv%snowind)

     !--- set up some clm1d options
     clm%istep = istep_pf
     clm%clm%soi_z = soi_z

     !--- setup -- set dz, ground stuff which is then pushed out to clm1d and tiles in clm_setup
     call parflow_read_ground(host, clm%drv, ix,iy, gnx,gny, latlon, sand,clay, &
          color_index, fractional_ground)
     call host_to_clm_ground_properties(host, clm, latlon, sand, clay, color_index, fractional_ground)
     
     !--- setup -- transfers grid to tile/clm1d properties
     call clm_setup(clm)

     !--- setup -- set data, other parameters FIXME: move these above clm_setup if possible for cleanliness --etc
     call host_to_clm_dz(host, pdz, pf_dz_mult, clm)
     call host_to_clm_et_controls(host, beta_typepf, veg_water_stress_typepf, wilting_pointpf, &
          field_capacitypf, res_satpf, clm)
     call host_to_clm_irrigation(host, irr_typepf, irr_cyclepf, irr_ratepf, irr_startpf, &
          irr_stoppf, irr_thresholdpf, irr_thresholdtypepf, clm)
     call host_to_clm_wc(host, porosity, saturation, clm)
     call host_to_clm_tksat_from_porosity(host, porosity, clm)
     if (clm_forc_veg == 1) then
        call host_to_clm_forced_vegetation(host, lai_pf, sai_pf, z0m_pf, displa_pf, clm)
     end if

     !--- potential clobber stuff with restart data
     call clm_restart(1, clm, -1)
     
  else
     !--- if not intial time, just open logfiles and set water data
     call host_to_clm_wc(host, porosity, saturation, clm)
  end if

  !--- set pressure, get forcing met data, get ready to take the step
  call host_to_clm_pressure(host, pressure, 1000., clm) ! note unit conversion from m to mm
  call host_to_clm_met_data(host, sw_pf, lw_pf, prcp_pf, tas_pf, qatm_pf, u_pf, v_pf, patm_pf, clm)

  !--- advance the time step
  call clm_advance_time(clm, host, istep_pf, time*3600, dt*3600) ! note unit conversion from hrs to seconds
  
  !--- potentially write 2d output
  if ((clm%io%dump_current == 1) .or. &
       (mod(clm%istep, clm%io%dump_interval) == 0)) then
     if (clm%io%write_binaries /= 0) then
        ! Call subroutine to open (2D-) output files
        call parflow_open_files(clm%clm, clm%drv,clm%rank, ix,iy,istep_pf,clm_output_dir, &
             clm_output_dir_length,clm_bin_output_dir) 

        ! Call subroutine to write 2D output
        call drv_2dout(clm%drv,clm%grid,clm%clm)

        ! Call to subroutine to close (2D-) output files
        call parflow_close_files(clm%clm,clm%drv)
        
     end if
  end if

  !--- Copy values from 2D CLM arrays to PF arrays for printing from PF (as Silo)
  ! NOTE the use of junk arrays which are tossed.  FIXME -- make these pointers with NULL? --etc
  call clm_to_host_total_energy_fluxes(host, clm, eflx_lh_pf, eflx_sh_pf, eflx_lwrad_pf, eflx_grnd_pf)
  call clm_to_host_mass_fluxes(host, clm, qflx_tot_pf, qflx_grnd_pf, qflx_soi_pf, &
       qflx_eveg_pf, qflx_tveg_pf, qflx_in_pf, qirr_pf, qirr_inst_pf, irr_flag_pf, junk3d)
  call clm_to_host_diagnostics(host, clm, swe_pf, junk2d, junk2d, t_g_pf, junk2d, t_soi_pf)

  !--- write restart files?
  if (clm_next == 1) then
     if (clm_last_rst == 1) then
        d_stp = 0
     else
        d_stp = istep_pf
     end if

     if (clm_daily_rst == 1) then
        if ( (clm%drv%gmt==0.0).or.(clm%drv%endtime==1) ) then
           if (rank==0) then
              open(393, file="clm_restart.tcl",action="write")
              write(393,*) "set istep ",istep_pf
              close(393)
           end if  !  write istep corresponding to restart step

           call clm_restart(2, clm, d_stp)
        end if
     else
        call clm_restart(2,clm,d_stp)
     endif
  end if

  !--- copy back final forcing data
  call clm_to_host_total_mass_fluxes_combined(host, clm, evap_trans)
  evap_trans(:) = evap_trans(:) * 3600. / 1000. ! units conversion from mm/s to m/hr

  !--- perform mass balance check
  call parflow_check_mass_balance(clm%drv, clm%clm, clm%tile, evap_trans, &
       saturation, pressure, porosity, nx,ny,nz,host%j_incr,host%k_incr,ip,istep_pf)

  !--- write spatially-averaged BCs and ICs to file for diagnostics
  if (clm_write_logs == 1) then
     if (istep_pf == 1) call drv_pout(clm%drv, clm%tile, clm%clm, rank)
  end if

  !--- if at end of simulation...
  if (clm%drv%endtime == 1) then
     !--- close log files
     call io_close(clm%io)

     !--- destroy memory
     call clm_host_destroy(host)
     call clm_destroy(clm)
  end if
  
end subroutine clm_lsm
     


  
  
  
