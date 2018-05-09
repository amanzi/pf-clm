/*----------------------------------------------------------------------------*
 * ATS
 *
 * License: see $ATS_DIR/COPYRIGHT
 *
 * C-to-F90 wrapper for CLM
 *----------------------------------------------------------------------------*/


#include "ats_clm_interface.hh"

namespace ATS {
namespace CLM {

//
// sets the node depths based on dz
// ------------------------------------------------------------------
int32_t ats_to_clm_dz(const double* const dz) {
  ats_to_clm_dz(dz);
  return 0;
}

//
// Set soil properties
// ------------------------------------------------------------------
int32_t ats_to_clm_ground_properties(
  
  double *latlon, double *sand,
				     double *clay, double *color_index,
				     double *fractional_ground) {

}



int32_t ats_clm_init(int32_t ncells, int32_t ncolumns, int32_t col_inds,
		     int32_t rank, int32_t verbosity);


} // namespace CLM
} // namespace ATS
