#include <iostream>

#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"

#include "ats_clm_interface.hh"


int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  Epetra_MpiComm comm(MPI_COMM_SELF);
  
  std::cout << "Calling init()" << std::endl;
  int ncells_per = 15;
  int ncols = 1;
  int ncells = ncols * ncells_per;
  int startcode = 2; // new simulation
  int rank = 0; // logging only
  int verbosity = 3; // verbose!
  int npft = 1;
  Epetra_Map cell_map(ncells, 0, comm);
  Epetra_Map col_map(ncols, 0, comm);

  // init and allocate space
  ATS::CLM::init(ncells, ncols,
                 startcode, rank, verbosity);
  ATS::CLM::set_zero_time(2004);
  ATS::CLM::set_initial_state(300., 0.);

  // set the soil
  double latlon[ncols][2];
  for (int i=0; i!=ncols; ++i) {
    latlon[i][0] = 36.0;
    latlon[i][1] = -84.0;
  }
  Epetra_MultiVector sand(cell_map,1); sand.PutScalar(0.2);
  Epetra_MultiVector clay(cell_map,1); clay.PutScalar(0.4);
  std::vector<int> color_index(ncols, 1);
  double fractional_ground[ncols][NUM_LC_CLASSES];
  for (int c=0; c!=ncols; ++c) {
    for (int p=0; p!=NUM_LC_CLASSES; ++p) {
      fractional_ground[c][p] = 0.;
    }
    fractional_ground[c][3] = 1.0; // NOTE this is IGBP class 4, deciduous broadleaf
  }
  ATS::CLM::set_ground_properties(&latlon[0][0], sand, clay, color_index, &fractional_ground[0][0]);

  // begin setup (distributes to columns
  ATS::CLM::setup_begin();

  // set the grid spacing as 1m cells
  Epetra_MultiVector dz(cell_map, 1); dz.PutScalar(1.);
  ATS::CLM::set_dz(dz);

  // set up a bunch of control parameters
  ATS::CLM::set_et_controls(1, 2, 0.1, 1.0, 0.1);

  // not setting irrigation at this point
  //ATS::CLM::set_irrigation_controls();

  // finalize setup, move data to columns
  ATS::CLM::setup_end();

  // re-clobber, this is poor design and needs fixed
  ATS::CLM::set_dz(dz);

  // initial conditions
  Epetra_MultiVector poro(cell_map, 1);
  poro.PutScalar(0.45);
  Epetra_MultiVector sat(cell_map, 1);
  sat.PutScalar(1.0);
  ATS::CLM::set_wc(poro,sat);
  ATS::CLM::set_tksat_from_porosity(poro);

  // timestep
  Epetra_MultiVector pressure(cell_map, 1);
  pressure.PutScalar(101325.);
  ATS::CLM::set_pressure(pressure, 101325.);

  Epetra_MultiVector qSW(col_map,1);
  qSW.PutScalar(150.);
  Epetra_MultiVector qLW(col_map,1);
  qLW.PutScalar(150.);
  Epetra_MultiVector pRain(col_map,1);
  pRain.PutScalar(1.e-7);
  Epetra_MultiVector pSnow(col_map,1);
  pSnow.PutScalar(0.);
  Epetra_MultiVector air_temp(col_map,1);
  air_temp.PutScalar(274.15);
  Epetra_MultiVector rel_hum(col_map,1);
  rel_hum.PutScalar(0.9);
  Epetra_MultiVector wind(col_map,1);
  wind.PutScalar(3.);
  ATS::CLM::set_met_data(qSW, qLW, pRain, pSnow, air_temp, rel_hum, wind, 101325.);

  // set the start time, endtime
  ATS::CLM::advance_time(0, 43200.0, 3600.0); //units in seconds


  // get diagnostics for vis
  Epetra_MultiVector latent_heat(col_map, 1);
  Epetra_MultiVector sensible_heat(col_map, 1);
  Epetra_MultiVector lw_out(col_map, 1);
  Epetra_MultiVector conducted_e(col_map, 1);
  ATS::CLM::get_total_energy_fluxes(latent_heat, sensible_heat, lw_out, conducted_e);

  Epetra_MultiVector evap_total(col_map, 1);
  Epetra_MultiVector evap_ground(col_map, 1);
  Epetra_MultiVector evap_soil(col_map, 1);
  Epetra_MultiVector evap_canopy(col_map, 1);
  Epetra_MultiVector tran_veg(col_map, 1);
  Epetra_MultiVector influx(col_map, 1);
  Epetra_MultiVector irrigation(col_map, 1);
  Epetra_MultiVector inst_irrigation(cell_map, 1);
  Epetra_MultiVector irrigation_flag(col_map, 1);
  Epetra_MultiVector tran_soil(cell_map, 1);
  ATS::CLM::get_mass_fluxes(evap_total, evap_ground, evap_soil, evap_canopy, tran_veg,
                            influx, irrigation, inst_irrigation, irrigation_flag, tran_soil);

  Epetra_MultiVector swe(col_map, 1);
  Epetra_MultiVector snow_depth(col_map, 1);
  Epetra_MultiVector canopy_storage(col_map, 1);
  Epetra_MultiVector Tskin(col_map, 1);
  Epetra_MultiVector Tveg(col_map, 1);
  Epetra_MultiVector Tsoil(cell_map, 1);
  ATS::CLM::get_diagnostics(swe, snow_depth, canopy_storage, Tskin, Tveg, Tsoil);

  // get source terms for richards/overland
  Epetra_MultiVector qW_surf(col_map, 1);
  Epetra_MultiVector qW_subsurf(cell_map, 1);
  ATS::CLM::get_total_mass_fluxes(qW_surf, qW_subsurf);  

  qW_surf.Print(std::cout);
  qW_subsurf.Print(std::cout);
}
