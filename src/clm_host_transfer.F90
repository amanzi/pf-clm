! ------------------------------------------------------------------------------
! CLM-to-host data transfer functions
!
! This provides an interface for setting and getting data to/from CLM,
! and is used in writing a driver.
!
! This set of functions should not be edited by the user.  Instead, a
! hostfile should be written (see clm_host_structured or
! clm_host_unstructured for examples), and these functions should then
! be used.
!
!
! Author: Ethan Coon (ecoon _at_ lanl.gov)
! ------------------------------------------------------------------------------


module clm_host_transfer
  use clm_precision, only : r8,i4
  use clm_infnan, only : invalid
  use clm_host
  use clm_type_module 
  implicit none

  private

  public :: &
! host_to_clm: setters
       host_to_clm_ground_properties, &
       host_to_clm_dz, &
       host_to_clm_met_data, &
       host_to_clm_forced_vegetation, &
       host_to_clm_pressure, &
       host_to_clm_wc, &
       host_to_clm_tksat, &
       host_to_clm_tksat_from_porosity, &
       host_to_clm_irrigation, &
       host_to_clm_et_controls, &
  ! clm_to_host: getters
       clm_to_host_ground_energy_fluxes, &
       clm_to_host_total_energy_fluxes, &
       clm_to_host_mass_fluxes, &
       clm_to_host_total_mass_fluxes, &
       clm_to_host_total_mass_fluxes_combined, &
       clm_to_host_diagnostics  
  
contains

  !
  ! sets the node depths based on dz
  ! ------------------------------------------------------------------
  subroutine host_to_clm_dz(host, dz_base, dz_mult, clm)
    use clm1d_varpar, only : nlevsoi
    type(host_type),intent(in) :: host
    real(r8),intent(in) :: dz_mult(host%ncells_g)  ! this oddity of having a scalar dz and a vector
    real(r8),intent(in) :: dz_base                 ! multiplier is thanks to ParFlow history, should
                                                   ! get removed eventually
                                                   ! dz_base is in [m], dz_mult in [-]
                                                   ! or alternatively dz_base is 1 and dz_mult in [m]
    type(clm_type) :: clm

    ! locals
    integer :: i,j,k,l,t
    
    ! Set up variable DZ over root column
    ! -- Copy dz multipliers for root zone cells from host grid to 1D
    !    array, computing zi (interface centroid) and z (cell
    !    centroid)
    ! -- Reset rooting fractions based on this dz

    ! FIXME Should probably push these dz/z/zi onto a grid attribute,
    ! and then in g2clm push to clm attributes? --etc
    do t = 1,clm%drv%nch

       i = clm%tile(t)%col
       j = clm%tile(t)%row

       clm%clm(t)%zi(0) = 0.   

       if (clm%grid(i,j)%rootfr /= clm%drv%udef) then 
          ! use rooting fraction assigned by user in clmin file.
          ! NOTE: this is stupid and likely unphysical --etc
          do k=1,nlevsoi
             clm%clm(t)%rootfr(k)=clm%grid(i,j)%rootfr    
          enddo
       else
       
          if (host%planar_mask(t) == 1) then
             ! loop over CLM-active cells setting dz, zi
             do k = 1, nlevsoi
                l = host_cell_index(host, i,j,k)
                clm%clm(t)%dz(k) = dz_base * dz_mult(l)
                clm%clm(t)%z(k) = clm%clm(t)%zi(k-1) + 0.5 * clm%clm(t)%dz(k)
                clm%clm(t)%zi(k) = clm%clm(t)%zi(k-1) + clm%clm(t)%dz(k)
             enddo

             do k = 1, nlevsoi-1
                ! evaluates the integral of the rooting fraction function
                ! from zi(k-1) to zi(k)
                clm%clm(t)%rootfr(k) = .5*( exp(-clm%tile(t)%roota*clm%clm(t)%zi(k-1))  &
                     + exp(-clm%tile(t)%rootb*clm%clm(t)%zi(k-1))  &
                     - exp(-clm%tile(t)%roota*clm%clm(t)%zi(k))  &
                     - exp(-clm%tile(t)%rootb*clm%clm(t)%zi(k)) )
             enddo
             ! evalutes the integral of the rooting fraction function
             ! from zi(nlevsoi-1) to zi=inf to ensure the sum of all
             ! rooting fractions is 1.
             clm%clm(t)%rootfr(nlevsoi)=.5*( exp(-clm%tile(t)%roota*clm%clm(t)%zi(nlevsoi-1))&
                  + exp(-clm%tile(t)%rootb*clm%clm(t)%zi(nlevsoi-1)))
          end if
       end if
    end do
  end subroutine host_to_clm_dz
  
  !
  ! Set soil properties
  ! ------------------------------------------------------------------
  subroutine host_to_clm_ground_properties(host, latlon, sand, clay, &
       color_index, fractional_ground, clm)
    use clm1d_varpar, only : nlevsoi
    implicit none
    type(host_type),intent(in) :: host
    type(clm_type), intent(inout) :: clm
    real(r8),intent(in) :: latlon(2,host%ncolumns_g)    ! latitude,longitude [decimal degrees]
    real(r8),intent(in) :: sand(host%ncells_g)          ! percent sand FIXME: 0-1 or 0-100? --etc
    real(r8),intent(in) :: clay(host%ncells_g)          ! percent clay FIXME: 0-1 or 0-100? --etc
    integer(i4),intent(in) :: color_index(host%ncolumns_g)  ! color index FIXME: document! --etc
    real(r8),intent(in) :: fractional_ground(clm%drv%nt,host%ncolumns_g) ! fraction of land surface of type t

    ! local
    integer :: i,j,k,l,m

    ! push data into the grid
    do i=1,clm%drv%nc
       do j=1,clm%drv%nr
          l = host_column_index(host, i,j)

          clm%grid(i,j)%latdeg = latlon(1,l)
          clm%grid(i,j)%londeg = latlon(2,l)
          clm%grid(i,j)%isoicol = color_index(l)

          do m = 1,clm%drv%nt
             clm%grid(i,j)%fgrd(m) = fractional_ground(m,l)
          end do

          do k = 1,nlevsoi
             l = host_cell_index(host, i,j,k)
             clm%grid(i,j)%sand(k) = sand(l)
             clm%grid(i,j)%clay(k) = clay(l)
          end do
       end do
    end do
  end subroutine host_to_clm_ground_properties

  
  !
  ! Sets the Meterological data forcing from the host code
  ! ------------------------------------------------------------------
  subroutine host_to_clm_met_data(host, eflx_swin, eflx_lwin, precip, &
       air_temp, air_spec_hum, wind_x, wind_y, patm, clm)
    use clm1d_varcon, only : tfrz, tcrit
    implicit none
    type(host_type),intent(in) :: host
    real(r8),intent(in) :: eflx_swin(host%ncolumns_g) ! shortwave incoming radiation [W/m^2]
    real(r8),intent(in) :: eflx_lwin(host%ncolumns_g) ! longwave incoming radiation [W/m^2]
    real(r8),intent(in) :: precip(host%ncolumns_g) ! precipitation rate [mm/s]
    real(r8),intent(in) :: air_temp(host%ncolumns_g) ! air temperature [K]
    real(r8),intent(in) :: air_spec_hum(host%ncolumns_g) ! air specific humidity [kg/kg]
    real(r8),intent(in) :: wind_x(host%ncolumns_g) ! wind speed, eastward direction [m/s]
    real(r8),intent(in) :: wind_y(host%ncolumns_g) ! wind speed, northward direction [m/s]
    real(r8),intent(in) :: patm(host%ncolumns_g) ! atmospheric pressure [Pa]
    type(clm_type) :: clm

    ! locals
    integer :: t,i,j,l

    do t=1,clm%drv%nch
       i=clm%tile(t)%col
       j=clm%tile(t)%row
       l = host_column_index(host, i,j)

       clm%clm(t)%forc_lwrad = eflx_lwin(l)
       clm%clm(t)%forc_t = air_temp(l)
       clm%clm(t)%forc_u = wind_x(l)
       clm%clm(t)%forc_v = wind_y(l)
       clm%clm(t)%forc_pbot = patm(l)
       clm%clm(t)%forc_q = air_spec_hum(l)

       ! air density
       clm%clm(t)%forc_rho = clm%clm(t)%forc_pbot / (clm%clm(t)%forc_t*2.8704e2)

       ! solar SW, partitioned
       clm%clm(t)%forc_solad(1)   = eflx_swin(l)*35./100.   ! direct/beam solar radiation
       clm%clm(t)%forc_solad(2)   = eflx_swin(l)*35./100.   ! 
       clm%clm(t)%forc_solai(1)   = eflx_swin(l)*15./100.   ! indirect/diffuse solar radiation
       clm%clm(t)%forc_solai(2)   = eflx_swin(l)*15./100.   ! 

       ! precip partitioning
       if (precip(l) > 0.) then
          if (clm%clm(t)%forc_t > (tfrz+tcrit)) then
             clm%clm(t)%itypprc = 1
             clm%clm(t)%forc_rain = precip(l)
             clm%clm(t)%forc_snow = 0.
          else
             clm%clm(t)%itypprc = 2
             clm%clm(t)%forc_rain = 0.
             clm%clm(t)%forc_snow = precip(l)
          end if
       else
          clm%clm(t)%itypprc = 0
          clm%clm(t)%forc_rain = 0.
          clm%clm(t)%forc_snow = 0.
       end if       
    end do
  end subroutine host_to_clm_met_data

  !
  ! Sets the vegetation data from the host code, for use in case of
  ! prescribed vegetation.
  !
  ! Expected units: pressure [mm]
  ! ------------------------------------------------------------------
  subroutine host_to_clm_forced_vegetation(host, lai, sai, z0m, displacement_ht, clm)
    implicit none
    type(host_type),intent(in) :: host
    real(r8),intent(in) :: lai(host%ncolumns_g) ! exposed leaf area index [-]
    real(r8),intent(in) :: sai(host%ncolumns_g) ! exposed stem area index [-]
    real(r8),intent(in) :: z0m(host%ncolumns_g) ! aerodynamic roughness length [m]
    real(r8),intent(in) :: displacement_ht(host%ncolumns_g) ! displacement height [m]
    type(clm_type) :: clm

    ! locals
    integer :: t,i,j,l

    do t=1,clm%drv%nch
       i=clm%tile(t)%col
       j=clm%tile(t)%row
       l = host_column_index(host, i,j)

       clm%clm(t)%elai = lai(l)
       clm%clm(t)%esai = sai(l)
       clm%clm(t)%z0m = z0m(l)
       clm%clm(t)%displa = displacement_ht(l)
    end do
  end subroutine host_to_clm_forced_vegetation
  
  !
  ! Copies pressure from a host array to a clm tiled array.
  !
  ! Expected units: pressure [mm]
  ! ------------------------------------------------------------------
  subroutine host_to_clm_pressure(host, pressure, unit_conversion, clm)
    use clm1d_varpar, only : nlevsoi
    implicit none
    type(host_type),intent(in) :: host
    real(r8),intent(in) :: pressure(host%ncells_g)  ! [mm] if unit_conversion is 1
    real(r8),intent(in) :: unit_conversion
    type(clm_type) :: clm

    ! locals
    integer :: t,i,j,k,l

    do t=1,clm%drv%nch
       if(host%planar_mask(t) == 1) then
          i=clm%tile(t)%col
          j=clm%tile(t)%row

          do k = 1,nlevsoi
             l = host_cell_index(host, i,j,k)
             clm%clm(t)%pf_press(k) = pressure(l) * unit_conversion
          end do
       end if
    end do
  end subroutine host_to_clm_pressure
  
  !
  ! Copies porosity from a host array to a clm tiled array.
  ! ------------------------------------------------------------------
  subroutine host_to_clm_wc(host, porosity, saturation, clm)
    use clm1d_varpar, only : nlevsoi
    use clm1d_varcon, only : denh2o
    implicit none
    type(host_type),intent(in) :: host
    real(r8),intent(in) :: porosity(host%ncells_g)      ! [-]
    real(r8),intent(in) :: saturation(host%ncells_g)    ! [-]
    type(clm_type) :: clm

    ! locals
    integer :: t,i,j,k,l

    do t=1,clm%drv%nch
       if(host%planar_mask(t) == 1) then
          i=clm%tile(t)%col
          j=clm%tile(t)%row

          do k = 1,nlevsoi
             l = host_cell_index(host, i,j,k)
             clm%clm(t)%watsat(k) = porosity(l)
             clm%clm(t)%pf_vol_liq(k) = porosity(l) * saturation(l)

             ! KNOWN BUG FIX:
             ! The old code read:
             ! clm%clm(t)%h2osoi_liq(k) = porosity(l) * saturation(l)*clm%clm(t)%dz(1)*denh2o
             ! THIS SHOULD READ:
             clm%clm(t)%h2osoi_liq(k) = porosity(l) * saturation(l)*clm%clm(t)%dz(k)*denh2o
             ! --etc
          end do
       end if
    end do
  end subroutine host_to_clm_wc


  !
  ! Copies thermal conductivity for the saturated soil from a host
  ! array to a clm tiled array.
  ! ------------------------------------------------------------------
  subroutine host_to_clm_tksat(host, tksat, clm)
    use clm1d_varpar, only : nlevsoi
    implicit none
    type(host_type),intent(in) :: host
    real(r8),intent(in) :: tksat(host%ncells_g) ! thermal conductivity of saturated soil [W / m-K]
    type(clm_type) :: clm

    ! locals
    integer :: t,i,j,k,l

    do t=1,clm%drv%nch
       if(host%planar_mask(t) == 1) then
          i=clm%tile(t)%col
          j=clm%tile(t)%row

          do k = 1,nlevsoi
             l = host_cell_index(host, i,j,k)
             clm%clm(t)%tksatu(k) = tksat(l)
          end do
       end if
    end do
  end subroutine host_to_clm_tksat


  !
  ! Calculates thermal conductivity of the saturated soil from a model
  ! based on porosity and thermal conductivity of the grain material.
  ! ------------------------------------------------------------------
  subroutine host_to_clm_tksat_from_porosity(host, poro, clm)
    use clm1d_varpar, only : nlevsoi
    implicit none
    type(host_type),intent(in) :: host
    real(r8),intent(in) :: poro(host%ncells_g)  ! [-]
    type(clm_type) :: clm

    ! locals
    integer :: t,i,j,k,l

    do t=1,clm%drv%nch
       if(host%planar_mask(t) == 1) then
          i=clm%tile(t)%col
          j=clm%tile(t)%row

          do k = 1,nlevsoi
             l = host_cell_index(host, i,j,k)
             clm%clm(t)%tksatu(k) = clm%clm(t)%tkmg(k)*0.57**poro(l)
          end do
       end if
    end do
  end subroutine host_to_clm_tksat_from_porosity


  !
  ! Sets irrigiation data
  !
  ! Expected units: 
  ! ------------------------------------------------------------------
  subroutine host_to_clm_irrigation(host, irr_type, irr_cycle, &
       irr_rate, irr_start, irr_stop, irr_threshold, irr_thresholdtype, clm)
    implicit none
    type(host_type),intent(in) :: host

    ! irrigation keys
    integer, intent(in) :: irr_type            ! irrigation type flag (0=none,1=spray,2=drip,3=instant)
    integer, intent(in) :: irr_cycle           ! irrigation cycle flag (0=constant,1=deficit)
    real(r8), intent(in) :: irr_rate           ! irrigation application rate for spray and drip [mm/s]
    real(r8), intent(in) :: irr_start          ! irrigation daily start time for constant cycle
    real(r8), intent(in) :: irr_stop           ! irrigation daily stop tie for constant cycle
    real(r8), intent(in) :: irr_threshold      ! irrigation threshold criteria for deficit cycle
                                                 ! (units of soil moisture content)
    integer, intent(in)  :: irr_thresholdtype  ! irrigation threshold criteria type -- top layer,
                                                 ! bottom layer, column avg
    type(clm_type) :: clm

    ! locals
    integer :: t

    do t=1,clm%drv%nch
       if(host%planar_mask(t) == 1) then
          ! for irrigation
          clm%clm(t)%irr_type           = irr_type
          clm%clm(t)%irr_cycle          = irr_cycle
          clm%clm(t)%irr_rate           = irr_rate
          clm%clm(t)%irr_start          = irr_start
          clm%clm(t)%irr_stop           = irr_stop
          clm%clm(t)%irr_threshold      = irr_threshold
          clm%clm(t)%threshold_type     = irr_thresholdtype
       end if
    end do
  end subroutine host_to_clm_irrigation

  !
  ! Sets ET controls
  !
  ! Expected units: 
  ! ------------------------------------------------------------------
  subroutine host_to_clm_et_controls(host, beta_type, veg_water_stress_type, wilting_point, &
       field_capacity, res_sat, clm)
    implicit none
    type(host_type),intent(in) :: host

    ! ET controls
    integer, intent(in) :: beta_type              ! beta formulation for bare soil Evap 0=none, 1=linear, 2=cos
    integer, intent(in) :: veg_water_stress_type  ! veg transpiration water stress formulation
                                                  ! 0=none, 1=press, 2=sm
    real(r8), intent(in) :: wilting_point         ! wilting point in m if press-type, in saturation
                                                  ! if soil moisture type
    real(r8), intent(in) :: field_capacity        ! field capacity for water stress same as units above
    real(r8), intent(in) :: res_sat               ! residual saturation

    type(clm_type) :: clm

    ! locals
    integer :: t

    do t=1,clm%drv%nch
       if(host%planar_mask(t) == 1) then
          ! for beta and veg stress formulations
          clm%clm(t)%beta_type          = beta_type
          clm%clm(t)%vegwaterstresstype = veg_water_stress_type
          clm%clm(t)%wilting_point      = wilting_point
          clm%clm(t)%field_capacity     = field_capacity
          clm%clm(t)%res_sat            = res_sat
       end if
    end do
  end subroutine host_to_clm_et_controls



  !
  ! Gets surface energy balance on the ground.
  !
  ! Expected units: all returned in W/m^2
  ! ------------------------------------------------------------------
  subroutine clm_to_host_ground_energy_fluxes(host, clm, eflx_lh, eflx_sh, eflx_lwrad_out, eflx_soil)
    implicit none
    type(host_type),intent(in) :: host
    type(clm_type), intent(in) :: clm
    real(r8),intent(inout) :: eflx_lh(host%ncolumns_g)          ! latent heat flux [W/m^2]
    real(r8),intent(inout) :: eflx_sh(host%ncolumns_g)          ! sensible heat from ground [W/m^2]
    real(r8),intent(inout) :: eflx_lwrad_out(host%ncolumns_g)   ! outgoing long-wave radiation from ground [W/m^2]
    real(r8),intent(inout) :: eflx_soil(host%ncolumns_g)        ! flux conducted to ground [W/m^2]

    ! local
    integer :: t,i,j,l
    
    do t=1,clm%drv%nch
       i = clm%tile(t)%col
       j = clm%tile(t)%row
       l = host_column_index(host, i,j)
       if (host%planar_mask(t) == 1) then
          eflx_lh(l) = clm%clm(t)%eflx_lh_grnd
          eflx_sh(l) = clm%clm(t)%eflx_sh_grnd
          eflx_lwrad_out(l) = clm%clm(t)%eflx_lwrad_grnd
          eflx_soil(l) = clm%clm(t)%eflx_soil_grnd
       else
          eflx_lh(l) = invalid
          eflx_sh(l) = invalid
          eflx_lwrad_out(l) = invalid
          eflx_soil(l) = invalid
       end if
    end do
  end subroutine clm_to_host_ground_energy_fluxes

  

  !
  ! Gets surface energy balance total fluxes (canopy + ground)
  !
  ! Expected units: all returned in W/m^2
  ! ------------------------------------------------------------------
  subroutine clm_to_host_total_energy_fluxes(host, clm, eflx_lh, eflx_sh, eflx_lwrad_out, eflx_soil)
    implicit none
    type(host_type),intent(in) :: host
    type(clm_type), intent(in) :: clm
    real(r8),intent(inout) :: eflx_lh(host%ncolumns_g)          ! latent heat flux [W/m^2]
    real(r8),intent(inout) :: eflx_sh(host%ncolumns_g)          ! sensible heat from ground [W/m^2]
    real(r8),intent(inout) :: eflx_lwrad_out(host%ncolumns_g)   ! outgoing long-wave radiation from ground [W/m^2]
    real(r8),intent(inout) :: eflx_soil(host%ncolumns_g)        ! flux conducted to ground [W/m^2]

    ! local
    integer :: t,i,j,l
    
    do t=1,clm%drv%nch
       i = clm%tile(t)%col
       j = clm%tile(t)%row
       l = host_column_index(host, i,j)
       if (host%planar_mask(t) == 1) then
          eflx_lh(l) = clm%clm(t)%eflx_lh_tot
          eflx_sh(l) = clm%clm(t)%eflx_sh_tot
          eflx_lwrad_out(l) = clm%clm(t)%eflx_lwrad_out
          eflx_soil(l) = clm%clm(t)%eflx_soil_grnd
       else
          eflx_lh(l) = invalid
          eflx_sh(l) = invalid
          eflx_lwrad_out(l) = invalid
          eflx_soil(l) = invalid
       end if          
    end do
  end subroutine clm_to_host_total_energy_fluxes


  !
  ! Gets mass fluxes due to surface energy balance relative to
  ! the ground.
  !
  ! Expected units: all returned in [mm/s]
  ! ------------------------------------------------------------------
  subroutine clm_to_host_mass_fluxes(host, clm, qflx_evap_tot, qflx_evap_ground, qflx_evap_soil, &
       qflx_evap_veg, qflx_tran_veg, qflx_infl, qflx_irr, qflx_irr_inst, irr_flag, qflx_tran_soil)
    use clm1d_varpar, only : nlevsoi
    implicit none
    type(host_type),intent(in) :: host
    type(clm_type), intent(in) :: clm
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

    ! local
    integer :: t,i,j,k,l
    
    do t=1,clm%drv%nch
       i = clm%tile(t)%col
       j = clm%tile(t)%row
       l = host_column_index(host, i,j)
       if (host%planar_mask(t) == 1) then
          qflx_evap_tot(l) = clm%clm(t)%qflx_evap_tot
          qflx_evap_ground(l) = clm%clm(t)%qflx_evap_grnd
          qflx_evap_soil(l) = clm%clm(t)%qflx_evap_soi
          qflx_evap_veg(l) = clm%clm(t)%qflx_evap_veg 
          qflx_tran_veg(l) = clm%clm(t)%qflx_tran_veg
          qflx_infl(l) = clm%clm(t)%qflx_infl 
          qflx_irr(l) = clm%clm(t)%qflx_qirr
          irr_flag(l) = clm%clm(t)%irr_flag
          do k = 1,nlevsoi
             l = host_cell_index(host, i,j,k)
             qflx_irr_inst(l) = clm%clm(t)%qflx_qirr_inst(k)
             qflx_tran_soil(l) = clm%clm(t)%qflx_tran_veg * clm%clm(t)%rootfr(k)
          end do

       else
          qflx_evap_tot(l) = invalid
          qflx_evap_ground(l) = invalid
          qflx_evap_soil(l) = invalid
          qflx_evap_veg(l) = invalid
          qflx_tran_veg(l) = invalid
          qflx_infl(l) = invalid
          qflx_irr(l) = invalid
          irr_flag(l) = invalid
          do k = 1,nlevsoi
             l = host_cell_index(host, i,j,k)
             qflx_irr_inst(l) = invalid
             qflx_tran_soil(l) = invalid
          end do
          
       end if
    end do
  end subroutine clm_to_host_mass_fluxes


  !
  ! Gets the aggregated mass fluxes -- this is what a host code should
  ! feel as a source.
  !
  ! Expected units: all returned in [mm/s]
  ! ------------------------------------------------------------------
  subroutine clm_to_host_total_mass_fluxes(host, clm, qflx_surface, qflx_subsurface)
    use clm1d_varpar, only : nlevsoi
    implicit none
    type(host_type),intent(in) :: host
    type(clm_type), intent(in) :: clm
    real(r8),intent(inout) :: qflx_surface(host%ncolumns_g)  ! total mass flux to/from the surface [mm/s]
    real(r8),intent(inout) :: qflx_subsurface(host%ncells_g) ! total mass flux to/from the subsurface [mm/m/s]

    ! local
    integer :: t,i,j,k,l
    
    do t=1,clm%drv%nch
       i = clm%tile(t)%col
       j = clm%tile(t)%row
       l = host_column_index(host, i,j)
       if (host%planar_mask(t) == 1) then
          qflx_surface(l) = clm%clm(t)%qflx_infl + clm%clm(t)%qflx_dew_grnd
          do k = 1,nlevsoi
             l = host_cell_index(host, i,j,k)
             qflx_subsurface(l) = (-clm%clm(t)%qflx_tran_veg*clm%clm(t)%rootfr(k) + clm%clm(t)%qflx_qirr_inst(k))&
                  / clm%clm(t)%dz(k)
          end do

       else
          qflx_surface(l) = 0.
          do k = 1,nlevsoi
             l = host_cell_index(host, i,j,k)
             qflx_subsurface(l) = 0.
          end do

       end if
    end do
  end subroutine clm_to_host_total_mass_fluxes


  !
  ! Gets the aggregated mass fluxes -- this is what a host code should
  ! feel as a source.  Puts the surface fluxes into the top cell of
  ! the subsurface.
  !
  ! Expected units: all returned in [mm/s]
  ! ------------------------------------------------------------------
  subroutine clm_to_host_total_mass_fluxes_combined(host, clm, qflx_subsurface)
    use clm1d_varpar, only : nlevsoi
    implicit none
    type(host_type),intent(in) :: host
    type(clm_type), intent(in) :: clm
    real(r8),intent(inout) :: qflx_subsurface(host%ncells_g) ! total mass flux to/from the subsurface [mm/m/s]

    ! local
    integer :: t,i,j,k,l
    
    do t=1,clm%drv%nch
       i = clm%tile(t)%col
       j = clm%tile(t)%row
       l = host_column_index(host, i,j)
       if (host%planar_mask(t) == 1) then
          do k = 1,nlevsoi
             l = host_cell_index(host, i,j,k)
             if (k == 1) then 
                qflx_subsurface(l) = (-clm%clm(t)%qflx_tran_veg*clm%clm(t)%rootfr(k) &
                     + clm%clm(t)%qflx_qirr_inst(k) + clm%clm(t)%qflx_infl) + clm%clm(t)%qflx_dew_grnd &
                     / clm%clm(t)%dz(k)
             else
                qflx_subsurface(l) = (-clm%clm(t)%qflx_tran_veg*clm%clm(t)%rootfr(k) &
                     + clm%clm(t)%qflx_qirr_inst(k)) &
                     / clm%clm(t)%dz(k)
             end if
          end do

       else
          do k = 1,nlevsoi
             l = host_cell_index(host, i,j,k)
             qflx_subsurface(l) = 0.
          end do

       end if
    end do
  end subroutine clm_to_host_total_mass_fluxes_combined
  

  !
  ! Gets surface energy balance total fluxes (canopy + ground)
  !
  ! Expected units: all returned in W/m^2
  ! ------------------------------------------------------------------
  subroutine clm_to_host_diagnostics(host, clm, swe, snow_depth, canopy_storage, T_skin, T_veg, T_soil)
    use clm1d_varpar, only : nlevsoi
    implicit none
    type(host_type),intent(in) :: host
    type(clm_type), intent(in) :: clm
    real(r8),intent(inout) :: swe(host%ncolumns_g) ! snow-water equivalent (mass) [mm]
    real(r8),intent(inout) :: canopy_storage(host%ncolumns_g) ! canopy storage (mass) [mm]
    real(r8),intent(inout) :: snow_depth(host%ncolumns_g)! snow depth [m]
    real(r8),intent(inout) :: T_skin(host%ncolumns_g) ! skin temperature [K]
    real(r8),intent(inout) :: T_veg(host%ncolumns_g) ! leaf temperature [K]
    real(r8),intent(inout) :: T_soil(host%ncells_g)   ! soil temperature [K]

    ! local
    integer :: t,i,j,k,l
    
    do t=1,clm%drv%nch
       i = clm%tile(t)%col
       j = clm%tile(t)%row
       l = host_column_index(host, i,j)
       if (host%planar_mask(t) == 1) then
          swe(l) = clm%clm(t)%h2osno
          snow_depth(l) = clm%clm(t)%snowdp
          canopy_storage(l) = clm%clm(t)%h2ocan
          T_skin(l) = clm%clm(t)%t_grnd
          T_veg(l) = clm%clm(t)%t_veg
          do k = 1,nlevsoi
             l = host_cell_index(host, i,j,k)
             T_soil(l) = clm%clm(t)%t_soisno(k)
          end do
       else
          swe(l) = invalid
          canopy_storage(l) = invalid
          snow_depth(l) = invalid
          T_skin(l) = invalid
          T_veg(l) = invalid
          do k = 1,nlevsoi
             l = host_cell_index(host, i,j,k)
             T_soil(l) = invalid
          end do
       end if          
    end do
  end subroutine clm_to_host_diagnostics


end module clm_host_transfer
