#ifndef ATS_CLM_INTERFACE_PRIVATE_HH_
#define ATS_CLM_INTERFACE_PRIVATE_HH_


extern "C" {
  void ats_to_clm_dz(const double* dz);
  void ats_to_clm_ground_properties(const double* latlon,
          const double* sand, const double* clay, const int* color_index,
          const double* fractional_ground);

  void ats_clm_init(int* ncells, int* ncolumns, int* col_inds, int* startcode,
                    int* rank, int* verbosity);
  void ats_clm_setup_begin();
  void ats_clm_setup_end();
  void ats_clm_advance_time(int* step, double* time, double* dt);


  void ats_to_clm_et_controls(int* beta_type, int* veg_water_stress_type,
          double* wilting_point, double* field_capacity,
          double* res_sat);
  void ats_to_clm_wc(double* porosity, double* sat);
  void ats_to_clm_pressure(double* pressure, double* patm);
  void ats_to_clm_tksat_from_porosity(double* porosity);
  void ats_to_clm_met_data(double* qSW, double* qLW, double* precip,
                           double* air_temp, double* rel_hum,
                           double* wind_x, double* wind_y, double* patm);

  void ats_clm_zero_time(double* zero_year);
}


#endif

