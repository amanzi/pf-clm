! ------------------------------------------------------------------------------
! ATS Driver module
!
! Effectively a singleton module for ATS points of entry.
!
! Author: Ethan Coon (ecoon _at_ lanl.gov)
! ------------------------------------------------------------------------------

module ats_clm
  use clm_precision
  use clm_type_module
  use clm_host
  use clm_host_transfer
  use clm1d_type_module
  use io_type_module
  implicit none

  private

! namespace-local data to hold the singleton clm and host instances  
  type(clm_type) :: clm
  type(host_type) :: host
  real(r8) :: year0
  
  public :: ats_to_clm_ground_properties, &
       ats_to_clm_dz, &
       ats_to_clm_met_data, &
       ats_to_clm_forced_vegetation, &
       ats_to_clm_pressure, &
       ats_to_clm_wc, &
       ats_to_clm_tksat_from_porosity, &
       ats_to_clm_irrigation, &
       ats_to_clm_et_controls, &
       clm_to_ats_ground_energy_fluxes, &
       clm_to_ats_total_energy_fluxes, &
       clm_to_ats_mass_fluxes, &
       clm_to_ats_total_mass_fluxes, &
       clm_to_ats_total_mass_fluxes_combined, &
       clm_to_ats_diagnostics, &
       ats_clm_init, &
       ats_clm_setup_begin, &
       ats_clm_setup_end, &
       ats_clm_advance_time

contains
  

  !
  ! sets the node depths based on dz
  ! ------------------------------------------------------------------
  subroutine ats_to_clm_dz(dz) bind(C)
    implicit none
    real(r8),intent(in) :: dz(host%ncells_g)
    call host_to_clm_dz(host, 1.d0, dz, clm)
  end subroutine ats_to_clm_dz
  
  !
  ! Set soil properties
  ! ------------------------------------------------------------------
  subroutine ats_to_clm_ground_properties(latlon, sand, clay, &
       color_index, fractional_ground) bind(C)
    implicit none
    real(r8),intent(in) :: latlon(2,host%ncolumns_g)    ! latitude,longitude [decimal degrees]
    real(r8),intent(in) :: sand(host%ncells_g)          ! percent sand FIXME: 0-1 or 0-100? --etc
    real(r8),intent(in) :: clay(host%ncells_g)          ! percent clay FIXME: 0-1 or 0-100? --etc
    integer(i4),intent(in) :: color_index(host%ncolumns_g)  ! color index FIXME: document! --etc
    real(r8),intent(in) :: fractional_ground(clm%drv%nt,host%ncolumns_g) ! fraction of land surface of type t
    call host_to_clm_ground_properties(host, latlon, sand, clay, color_index, fractional_ground, clm)
  end subroutine ats_to_clm_ground_properties

  
  !
  ! Sets the Meterological data forcing from the host code
  ! ------------------------------------------------------------------
  subroutine ats_to_clm_met_data(eflx_swin, eflx_lwin, precip, &
       air_temp, air_spec_hum, wind_x, wind_y, patm) bind(C)
    implicit none
    real(r8),intent(in) :: eflx_swin(host%ncolumns_g) ! shortwave incoming radiation [W/m^2]
    real(r8),intent(in) :: eflx_lwin(host%ncolumns_g) ! longwave incoming radiation [W/m^2]
    real(r8),intent(in) :: precip(host%ncolumns_g) ! precipitation rate [mm/s]
    real(r8),intent(in) :: air_temp(host%ncolumns_g) ! air temperature [K]
    real(r8),intent(in) :: air_spec_hum(host%ncolumns_g) ! air specific humidity [kg/kg]
    real(r8),intent(in) :: wind_x(host%ncolumns_g) ! wind speed, eastward direction [m/s]
    real(r8),intent(in) :: wind_y(host%ncolumns_g) ! wind speed, northward direction [m/s]
    real(r8),intent(in) :: patm(host%ncolumns_g) ! atmospheric pressure [Pa]
    call host_to_clm_met_data(host, eflx_swin, eflx_lwin, precip, &
         air_temp, air_spec_hum, wind_x, wind_y, patm, clm)
  end subroutine ats_to_clm_met_data

  !
  ! Sets the Meterological data forcing from the host code
  !
  ! Expected units: pressure [mm]
  ! ------------------------------------------------------------------
  subroutine ats_to_clm_forced_vegetation(lai, sai, z0m, displacement_ht)
    implicit none
    real(r8),intent(in) :: lai(host%ncolumns_g) ! exposed leaf area index [-]
    real(r8),intent(in) :: sai(host%ncolumns_g) ! exposed stem area index [-]
    real(r8),intent(in) :: z0m(host%ncolumns_g) ! aerodynamic roughness length [m]
    real(r8),intent(in) :: displacement_ht(host%ncolumns_g) ! displacement height [m]
    call host_to_clm_forced_vegetation(host, lai, sai, z0m, displacement_ht, clm)
  end subroutine ats_to_clm_forced_vegetation
  
  !
  ! Copies pressure from a host array to a clm tiled array.
  !
  ! Expected units: pressure [Pa]
  ! ------------------------------------------------------------------
  subroutine ats_to_clm_pressure(pressure, p_atm) bind(C)
    use clm1d_varcon, only : denh2o, grav
    implicit none
    real(r8),intent(in) :: pressure(host%ncells_g)  ! [Pa]
    real(r8),intent(in) :: p_atm
    
    ! local
    real(r8) :: pressure_adj(host%ncells_g)  ! [mm] 
    pressure_adj(:) = (pressure(:) - p_atm) / (denh2o*grav) * 1e3
    call host_to_clm_pressure(host, pressure_adj, 1.d0, clm)
  end subroutine ats_to_clm_pressure
  
  !
  ! Copies porosity from a host array to a clm tiled array.
  ! ------------------------------------------------------------------
  subroutine ats_to_clm_wc(porosity, saturation) bind(C)
    implicit none
    real(r8),intent(in) :: porosity(host%ncells_g)      ! [-]
    real(r8),intent(in) :: saturation(host%ncells_g)    ! [-]
    call host_to_clm_wc(host, porosity, saturation, clm)
  end subroutine ats_to_clm_wc


  !
  ! Calculates thermal conductivity of the saturated soil from a model
  ! based on porosity and thermal conductivity of the grain material.
  ! ------------------------------------------------------------------
  subroutine ats_to_clm_tksat_from_porosity(poro) bind(C)
    implicit none
    real(r8),intent(in) :: poro(host%ncells_g)  ! [-]
    call host_to_clm_tksat_from_porosity(host, poro, clm)
  end subroutine ats_to_clm_tksat_from_porosity


  !
  ! Sets irrigiation data
  !
  ! Expected units: 
  ! ------------------------------------------------------------------
  subroutine ats_to_clm_irrigation(irr_type, irr_cycle, &
       irr_rate, irr_start, irr_stop, irr_threshold, irr_thresholdtype)
    implicit none
    ! irrigation keys
    integer(i4), intent(in) :: irr_type            ! irrigation type flag (0=none,1=spray,2=drip,3=instant)
    integer(i4), intent(in) :: irr_cycle           ! irrigation cycle flag (0=constant,1=deficit)
    real(r8), intent(in) :: irr_rate           ! irrigation application rate for spray and drip [mm/s]
    real(r8), intent(in) :: irr_start          ! irrigation daily start time for constant cycle
    real(r8), intent(in) :: irr_stop           ! irrigation daily stop tie for constant cycle
    real(r8), intent(in) :: irr_threshold      ! irrigation threshold criteria for deficit cycle
                                                 ! (units of soil moisture content)
    integer(i4), intent(in)  :: irr_thresholdtype  ! irrigation threshold criteria type -- top layer,
                                                 ! bottom layer, column avg
    call host_to_clm_irrigation(host, irr_type, irr_cycle, &
         irr_rate, irr_start, irr_stop, irr_threshold, irr_thresholdtype, clm)
  end subroutine ats_to_clm_irrigation

  !
  ! Sets ET controls
  !
  ! Expected units: 
  ! ------------------------------------------------------------------
  subroutine ats_to_clm_et_controls(beta_type, veg_water_stress_type, wilting_point, &
       field_capacity, res_sat) bind(C)
    implicit none

    ! ET controls
    integer(i4), intent(in) :: beta_type              ! beta formulation for bare soil Evap 0=none, 1=linear, 2=cos
    integer(i4), intent(in) :: veg_water_stress_type  ! veg transpiration water stress formulation
                                                  ! 0=none, 1=press, 2=sm
    real(r8), intent(in) :: wilting_point         ! wilting point in m if press-type, in saturation
                                                  ! if soil moisture type
    real(r8), intent(in) :: field_capacity        ! field capacity for water stress same as units above
    real(r8), intent(in) :: res_sat               ! residual saturation
    call host_to_clm_et_controls(host, beta_type, veg_water_stress_type, wilting_point, &
         field_capacity, res_sat, clm)
  end subroutine ats_to_clm_et_controls



  !
  ! Gets surface energy balance on the ground.
  !
  ! Expected units: all returned in W/m^2
  ! ------------------------------------------------------------------
  subroutine clm_to_ats_ground_energy_fluxes(eflx_lh, eflx_sh, eflx_lwrad_out, eflx_soil)
    implicit none
    real(r8),intent(inout) :: eflx_lh(host%ncolumns_g)          ! latent heat flux [W/m^2]
    real(r8),intent(inout) :: eflx_sh(host%ncolumns_g)          ! sensible heat from ground [W/m^2]
    real(r8),intent(inout) :: eflx_lwrad_out(host%ncolumns_g)   ! outgoing long-wave radiation from ground [W/m^2]
    real(r8),intent(inout) :: eflx_soil(host%ncolumns_g)        ! flux conducted to ground [W/m^2]
    call clm_to_host_ground_energy_fluxes(host, clm, eflx_lh, eflx_sh, eflx_lwrad_out, eflx_soil)
  end subroutine clm_to_ats_ground_energy_fluxes

  

  !
  ! Gets surface energy balance total fluxes (canopy + ground)
  !
  ! Expected units: all returned in W/m^2
  ! ------------------------------------------------------------------
  subroutine clm_to_ats_total_energy_fluxes(eflx_lh, eflx_sh, eflx_lwrad_out, eflx_soil)
    implicit none
    real(r8),intent(inout) :: eflx_lh(host%ncolumns_g)          ! latent heat flux [W/m^2]
    real(r8),intent(inout) :: eflx_sh(host%ncolumns_g)          ! sensible heat from ground [W/m^2]
    real(r8),intent(inout) :: eflx_lwrad_out(host%ncolumns_g)   ! outgoing long-wave radiation from ground [W/m^2]
    real(r8),intent(inout) :: eflx_soil(host%ncolumns_g)        ! flux conducted to ground [W/m^2]
    call clm_to_host_total_energy_fluxes(host, clm, eflx_lh, eflx_sh, eflx_lwrad_out, eflx_soil)
  end subroutine clm_to_ats_total_energy_fluxes


  !
  ! Gets mass fluxes due to surface energy balance relative to
  ! the ground.
  !
  ! Expected units: all returned in [mm/s]
  ! ------------------------------------------------------------------
  subroutine clm_to_ats_mass_fluxes(qflx_evap_tot, qflx_evap_ground, qflx_evap_soil, &
       qflx_evap_veg, qflx_tran_veg, qflx_infl, qflx_irr, qflx_irr_inst, irr_flag, qflx_tran_soil)
    implicit none
    real(r8),intent(inout) :: qflx_evap_tot(host%ncolumns_g)    ! total evaporation [mm/s]
    real(r8),intent(inout) :: qflx_evap_ground(host%ncolumns_g) ! ground evaporation (does not
                                                                !  include snow sublimation) [mm/s]
    real(r8),intent(inout) :: qflx_evap_soil(host%ncolumns_g)   ! soil evaporation [mm/s]
    real(r8),intent(inout) :: qflx_evap_veg(host%ncolumns_g)    ! canopy evaporation [mm/s]
    real(r8),intent(inout) :: qflx_tran_veg(host%ncolumns_g)    ! canopy transpiration [mm/s]
    real(r8),intent(inout) :: qflx_infl(host%ncolumns_g)        ! net infiltration [mm/s]
    real(r8),intent(inout) :: qflx_irr(host%ncolumns_g)         ! irrigation sources [mm/s]
    real(r8),intent(inout) :: qflx_irr_inst(host%ncells_g)      ! instantaneous irrigation sources [mm/s]
    real(r8),intent(inout) :: irr_flag(host%ncolumns_g)         ! flag for irrigation type
    real(r8),intent(inout) :: qflx_tran_soil(host%ncells_g)     ! transpiration, distributed via roots [mm/s]
    call clm_to_host_mass_fluxes(host, clm, qflx_evap_tot, qflx_evap_ground, qflx_evap_soil, &
         qflx_evap_veg, qflx_tran_veg, qflx_infl, qflx_irr, qflx_irr_inst, irr_flag, qflx_tran_soil)
  end subroutine clm_to_ats_mass_fluxes


  !
  ! Gets the aggregated mass fluxes -- this is what a host code should
  ! feel as a source.
  !
  ! Expected units: all returned in [mm/s]
  ! ------------------------------------------------------------------
  subroutine clm_to_ats_total_mass_fluxes(qflx_surface, qflx_subsurface)
    implicit none
    real(r8),intent(inout) :: qflx_surface(host%ncolumns_g)  ! total mass flux to/from the surface [mm/s]
    real(r8),intent(inout) :: qflx_subsurface(host%ncells_g) ! total mass flux to/from the subsurface [1/s]
    call clm_to_host_total_mass_fluxes(host, clm, qflx_surface, qflx_subsurface)
  end subroutine clm_to_ats_total_mass_fluxes


  !
  ! Gets the aggregated mass fluxes -- this is what a host code should
  ! feel as a source.  Puts the surface fluxes into the top cell of
  ! the subsurface.
  !
  ! Expected units: all returned in [mm/s]
  ! ------------------------------------------------------------------
  subroutine clm_to_ats_total_mass_fluxes_combined(qflx_subsurface)
    implicit none
    real(r8),intent(inout) :: qflx_subsurface(host%ncells_g) ! total mass flux to/from the subsurface [1/s]
    call clm_to_host_total_mass_fluxes_combined(host, clm, qflx_subsurface)
  end subroutine clm_to_ats_total_mass_fluxes_combined
  

  !
  ! Gets surface energy balance total fluxes (canopy + ground)
  !
  ! Expected units: all returned in W/m^2
  ! ------------------------------------------------------------------
  subroutine clm_to_ats_diagnostics(swe, snow_depth, canopy_storage, T_skin, T_veg, T_soil)
    implicit none
    real(r8),intent(inout) :: swe(host%ncolumns_g) ! snow-water equivalent (mass) [kg/m^2]
    real(r8),intent(inout) :: canopy_storage(host%ncolumns_g) ! canopy storage (mass) [kg/m^2]
    real(r8),intent(inout) :: snow_depth(host%ncolumns_g)! snow depth [m]
    real(r8),intent(inout) :: T_skin(host%ncolumns_g) ! skin temperature [K]
    real(r8),intent(inout) :: T_veg(host%ncolumns_g) ! leaf temperature [K]
    real(r8),intent(inout) :: T_soil(host%ncells_g)   ! soil temperature [K]
    call clm_to_host_diagnostics(host, clm, swe, snow_depth, canopy_storage, T_skin, T_veg, T_soil)
  end subroutine clm_to_ats_diagnostics



  subroutine ats_clm_init(ncells, ncolumns, col_inds, startcode, rank, verbosity) bind(C)
    use clm_io_config,only : MAX_FILENAME_LENGTH
    implicit none
    integer(i4),intent(in) :: ncells, ncolumns, startcode, rank, verbosity
    integer(i4),intent(in) :: col_inds(ncolumns, 2)

    ! locals
    character(len=MAX_FILENAME_LENGTH) :: output_dir = "clm/"
    
    !--- initialize the clm master object
    call clm_init(clm, rank, ncolumns, 1, 18, verbosity)

    !--- open log files, set up io options
    print*, "Opening logfile: ", output_dir
    call io_open(clm%io, output_dir, rank, 0)
    clm%io%restart_last = 1
    clm%io%restart_daily = 0
    clm%io%dump_interval = -1
    clm%io%dump_current = 0
    clm%io%output_1d = 0
    clm%io%write_binaries = 0

    !--- initialize the host object, pushing grid info into it
    call host_init(host, ncells, ncells, ncolumns, ncolumns, col_inds)
    if (io_ok(clm%io, VERBOSITY_LOW)) call host_write_to_log(host, clm%io%log)

    !--- read clm input files, data
    ! FIXME -- can we remove this and use a setter? --etc
    ! FIXME -- fix io to not pass clm_write_logs, instead iounit --etc
    ! FIXME -- hard-coded 1 PFT per cell
    clm%drv%maxt = 1
    clm%drv%startcode = startcode
    clm%drv%clm_ic = startcode
    
    if (io_ok(clm%io, VERBOSITY_LOW)) then
       write(clm%io%log,*) "  CLM startcode for date (1=restart, 2=defined):", clm%drv%startcode
       write(clm%io%log,*) "  CLM IC (1=restart, 2=defined):", clm%drv%clm_ic
       if (clm%drv%startcode == 0) then
          write(clm%io%log,*) "ERROR: startcode = 0"
          stop
       end if
       if (clm%drv%clm_ic == 0) then
          write(clm%io%log,*) "ERROR: clm_ic = 0"
          stop
       end if
    end if
    
    clm%clm => clm1d_create_n(clm%ntiles, clm%drv%surfind, clm%drv%soilind, clm%drv%snowind)
  end subroutine ats_clm_init

  subroutine ats_clm_setup_begin() bind(C)
    implicit none
    call clm_setup_begin(clm)
  end subroutine ats_clm_setup_begin

  subroutine ats_clm_setup_end() bind(C)
    implicit none
    call clm_setup_end(clm)
  end subroutine ats_clm_setup_end

  subroutine ats_clm_advance_time(istep, time, dtime) bind(C)
    implicit none
    integer(i4),intent(in) :: istep
    real(r8),intent(in) :: time, dtime
    if ((istep.eq.0).or.(clm%drv%time < 0)) then
       ! tick to set initial time
       clm%drv%ts = nint(time)
       call drv_tick(clm%drv)
    end if
    call clm_advance_time(clm, host, istep, time, dtime)
  end subroutine ats_clm_advance_time

  subroutine ats_clm_zero_time(zero_year) bind(C)
    ! zero time, in years
    implicit none
    real(r8),intent(in) :: zero_year
    year0 = zero_year
    call drv_time2date(year0, clm%drv%doy, clm%drv%day, clm%drv%gmt,       &
         clm%drv%yr, clm%drv%mo, clm%drv%da, clm%drv%hr,        &
         clm%drv%mn, clm%drv%ss)
    call drv_time2date(year0, clm%drv%sdoy, clm%drv%day, clm%drv%sgmt,       &
         clm%drv%syr, clm%drv%smo, clm%drv%sda, clm%drv%shr,        &
         clm%drv%smn, clm%drv%sss)
  end subroutine ats_clm_zero_time
 end module ats_clm
