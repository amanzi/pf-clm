/*----------------------------------------------------------------------------*
 * ATS
 *
 * License: see $ATS_DIR/COPYRIGHT
 *
 * C-to-F90 wrapper for CLM
 *----------------------------------------------------------------------------*/


#include "ats_clm_interface.hh"
#include "ats_clm_interface_private.hh"

namespace ATS {
namespace CLM {


int init(int ncells, int ncolumns, int startcode, int rank, int verbosity) {
  int col_inds[ncolumns][2];
  int count = 0;
  int ncells_per = ncells / ncolumns;
  for (int i=0; i!=ncolumns; ++i) {
    col_inds[i][0] = count+1;
    count += ncells_per;
    col_inds[i][1] = count;
  }

  ats_clm_init(&ncells, &ncolumns, &col_inds[0][0], &startcode, &rank, &verbosity);
  return 0;
}

//
// sets the node depths based on dz
// ------------------------------------------------------------------
int set_dz(const std::vector<double>& dz) {
  ats_to_clm_dz(dz.data());
  return 0;
}

int set_ground_properties(double* latlon,
        const std::vector<double>& sand, const std::vector<double>& clay,
        const std::vector<int>& color_index,
        double* fractional_ground) {
  ats_to_clm_ground_properties(latlon, sand.data(), clay.data(),
          color_index.data(), fractional_ground);
  return 0;
}

int set_wc(Epetra_MultiVector& porosity, Epetra_MultiVector& saturation) {
  ats_to_clm_wc(porosity[0], saturation[0]);
  return 0;
}

int set_tksat_from_porosity(Epetra_MultiVector& porosity) {
  ats_to_clm_tksat_from_porosity(porosity[0]);
  return 0;
}

int set_pressure(Epetra_MultiVector& pressure, double patm) {
  ats_to_clm_pressure(pressure[0], &patm);
  return 0;
}

int set_et_controls(int beta_type, int veg_water_stress_type,
                    double wilting_point, double field_capacity,
                    double res_sat) {
  ats_to_clm_et_controls(&beta_type, &veg_water_stress_type,
                         &wilting_point, &field_capacity, &res_sat);
  return 0;
}

int set_met_data(Epetra_MultiVector& qSW, Epetra_MultiVector& qLW, Epetra_MultiVector& pRain,
                 Epetra_MultiVector& pSnow, Epetra_MultiVector& air_temp,
                 Epetra_MultiVector& rel_hum, Epetra_MultiVector& wind_u,
                 double patm) {
  Epetra_MultiVector precip(pRain);
  precip.Update(1000.,pSnow,1000.); // converts m/s --> mm/s
  

  Epetra_MultiVector wind_y(wind_u);
  wind_y.PutScalar(0.);
  Epetra_MultiVector patm_v(rel_hum);
  patm_v.PutScalar(patm);
  
  ats_to_clm_met_data(qSW[0], qLW[0], precip[0], air_temp[0], rel_hum[0], wind_u[0], wind_y[0], patm_v[0]);
  return 0;
}



int setup_begin() {
  ats_clm_setup_begin();
  return 0;
}

int setup_end() {
  ats_clm_setup_end();
  return 0;
}

int advance_time(int step, double time, double dt) {
  ats_clm_advance_time(&step, &time, &dt);
  return 0;
}

int set_zero_time(double time) {
  // expects zero time in years
  ats_clm_zero_time(&time);
  return 0;
}
  



} // namespace CLM
} // namespace ATS
