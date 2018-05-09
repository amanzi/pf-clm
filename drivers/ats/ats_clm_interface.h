/*----------------------------------------------------------------------------*
 * ATS
 *
 * License: see $ATS_DIR/COPYRIGHT
 *
 * C-to-F90 wrapper for CLM
 *----------------------------------------------------------------------------*/


#ifndef ATS_CLM_INTERFACE_HH_
#define ATS_CLM_INTERFACE_HH_

namespace ATS {
namespace CLM {


//
// sets the node depths based on dz
// ------------------------------------------------------------------
int32_t ats_to_clm_dz(const std::vector<double>& dz);

//
// Set soil properties
// ------------------------------------------------------------------
int32_t ats_to_clm_ground_properties(double *latlon, double *sand,
				     double *clay, double *color_index,
				     double *fractional_ground);

// //
// // Sets the Meterological data forcing from the host code
// // ------------------------------------------------------------------
// int32_t ats_to_clm_met_data(double *eflx_swin, double *eflx_lwin, double *precip,
// 			    double *air_temp, double *air_spec_hum,
// 			    double *wind_x, double *wind_y,
// 			    double *patm);

// //
// // Sets the Meterological data forcing from the host code
// //
// // Expected units: pressure [mm]
// // ------------------------------------------------------------------
// int32_t ats_to_clm_forced_vegetation(double *lai, double *sai,
// 				     double *z0m, double *displacement_ht);

// //
// // Copies pressure from a host array to a clm tiled array.
// //
// // Expected units: pressure [mm]
// // ------------------------------------------------------------------
// int32_t ats_to_clm_pressure(double *pressure, double *p_atm);

// //
// // Copies porosity from a host array to a clm tiled array.
// // ------------------------------------------------------------------
// int32_t ats_to_clm_wc(double *porosity, double *saturation);
  
// //
// // Calculates thermal conductivity of the saturated soil from a model
// // based on porosity and thermal conductivity of the grain material.
// // ------------------------------------------------------------------
// int32_t ats_to_clm_tksat_from_porosity(double *poro);


// //
// // Sets irrigiation data
// //
// // Expected units: 
// // ------------------------------------------------------------------
// int32_t ats_to_clm_irrigation(int32_t *irr_type, int32_t *irr_cycle,
// 			      double *irr_rate, double *irr_start, double *irr_stop,
// 			      double *irr_threshold, int32_t *irr_thresholdtype);


// //
// // Sets ET controls
// //
// // Expected units: 
// // ------------------------------------------------------------------
// int32_t ats_to_clm_et_controls(int32_t* beta_type, int32_t* veg_water_stress_type,
// 			       double *wilting_point, double f*ield_capacity,
// 			       double *res_sat);
    
  
// //
// // Gets surface energy balance on the ground.
// //
// // Expected units: all returned in W/m^2
// // ------------------------------------------------------------------
// int32_t clm_to_ats_ground_energy_fluxes(double *eflx_lh, double *eflx_sh,
// 					double *eflx_lwrad_out, double *eflx_soil);


// //
// // Gets surface energy balance total fluxes (canopy + ground)
// //
// // Expected units: all returned in W/m^2
// // ------------------------------------------------------------------
// int32_t clm_to_ats_total_energy_fluxes(double *eflx_lh, double *eflx_sh,
// 				       double *eflx_lwrad_out, double *eflx_soil);


// //
// // Gets mass fluxes due to surface energy balance relative to
// // the ground.
// //
// // Expected units: all returned in [mm/s]
// // ------------------------------------------------------------------
// int32_t clm_to_ats_mass_fluxes(double *qflx_evap_tot, double *qflx_evap_ground,
// 			       double *qflx_evap_soil, double *qflx_evap_veg,
// 			       double *qflx_tran_veg, double *qflx_infl,
// 			       double *qflx_irr, double *qflx_irr_inst,
// 			       double *irr_flag, double *qflx_tran_soil);


// //
// // Gets the aggregated mass fluxes -- this is what a host code should
// // feel as a source.
// //
// // Expected units: all returned in [mm/s]
// // ------------------------------------------------------------------
// int32_t clm_to_ats_total_mass_fluxes(double *qflx_surface, double *qflx_subsurface);


// //
// // Gets the aggregated mass fluxes -- this is what a host code should
// // feel as a source.  Puts the surface fluxes into the top cell of
// // the subsurface.
// //
// // Expected units: all returned in [mm/s]
// // ------------------------------------------------------------------
// int32_t clm_to_ats_total_mass_fluxes_combined(double *qflx_subsurface);


// //
// // Gets surface energy balance total fluxes (canopy + ground)
// //
// // Expected units: all returned in W/m^2
// // ------------------------------------------------------------------
// int32_t clm_to_ats_diagnostics(double *swe, double *snow_depth,
// 			       double *canopy_storage, double *T_skin,
// 			       double *T_veg, double *T_soil);



int32_t ats_clm_init(int32_t ncells, int32_t ncolumns, int32_t col_inds,
		     int32_t rank, int32_t verbosity);
    

} // namespace CLM
} // namespace ATS
  

#endif
