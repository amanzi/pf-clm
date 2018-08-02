/*----------------------------------------------------------------------------*
 * ATS
 *
 * License: see $ATS_DIR/COPYRIGHT
 *
 * C-to-F90 wrapper for CLM
 *----------------------------------------------------------------------------*/


#ifndef ATS_CLM_INTERFACE_HH_
#define ATS_CLM_INTERFACE_HH_

#include <cstdint>
#include <vector>

#include "Epetra_MultiVector.h"

#define NUM_LC_CLASSES 18

namespace ATS {
namespace CLM {


//
// sets the node depths based on dz
// ------------------------------------------------------------------
int set_dz(const std::vector<double>& dz);

//
// Set soil properties
// ------------------------------------------------------------------
int set_ground_properties(double* latlon,
        const std::vector<double>& sand, const std::vector<double>& clay,
        const std::vector<int>& color_index,
        double* fractional_ground);

//
// Copies porosity from a host array to a clm tiled array.
// ------------------------------------------------------------------
int set_wc(Epetra_MultiVector& porosity, Epetra_MultiVector& saturation);
int set_tksat_from_porosity(Epetra_MultiVector& porosity);
int set_pressure(Epetra_MultiVector& pressure, double patm); // this requires more thought in frozen cases potentially
int set_met_data(Epetra_MultiVector& qSW, Epetra_MultiVector& qLW, Epetra_MultiVector& pRain,
                 Epetra_MultiVector& pSnow, Epetra_MultiVector& air_temp,
                 Epetra_MultiVector& rel_hum, Epetra_MultiVector& wind_u,
                 double patm);

//
// Sets ET controls
//
// Expected units: 
// ------------------------------------------------------------------
int set_et_controls(int beta_type, int veg_water_stress_type,
                    double wilting_point, double field_capacity,
                    double res_sat);    

int init(int ncells, int ncolumns, int startcode,
         int rank, int verbosity);
int set_zero_time(double zero_time);
    
int setup_begin();
int setup_end();
int advance_time(int step, double time, double dt);

} // namespace CLM
} // namespace ATS
  

#endif
