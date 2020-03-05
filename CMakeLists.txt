cmake_minimum_required(VERSION 3.9)

set(CLM_LIB clm)
project(${CLM_LIB} LANGUAGES Fortran)

option(UNSTRUCTURED_HOST "ON/OFF builds CLM library for unstructured/structured host driver [ON]" ON)

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  message(STATUS "Setting CLM installation prefix as: ${CMAKE_SOURCE_DIR}" )
  set (CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}" CACHE PATH "default install path" FORCE )
else()
  message(STATUS "Installing CLM in provided installation prefix: ${CMAKE_INSTALL_PREFIX}" )
endif()

if (UNSTRUCTURED_HOST)
  message(STATUS "Configuring CLM library for use with unstructured host driver")
  set(HOST_SOURCE clm_host_unstructured)
else()
  message(STATUS "Configuring CLM library for use with unstructured host driver")
  set(HOST_SOURCE clm_host_structured)
endif()

if (BUILD_SHARED_LIBS)
  message(STATUS "Building CLM as a shared library" )
else()
  message(STATUS "Building CLM as a static library" )
endif()

set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/include )

set (SRC_FILES
src/utils/clm_precision.F90
src/utils/clm_infnan.F90
src/utils/clm_abort.F90
src/utils/clm_io_config.F90
src/clm1d/clm1d_varpar.F90
src/clm1d/clm1d_varcon.F90
src/clm1d/clm1d_type.F90
src/clm1d/clm1d_typini.F90
src/clm1d/clm1d_balchk.F90
src/clm1d/clm1d_lake.F90
src/clm1d/clm1d_subdiv.F90
src/clm1d/clm1d_combin.F90
src/clm1d/clm1d_leaftem.F90
src/clm1d/clm1d_surfalb.F90
src/clm1d/clm1d_combo.F90
src/clm1d/clm1d_main.F90
src/clm1d/clm1d_surfrad.F90
src/clm1d/clm1d_compact.F90
src/clm1d/clm1d_meltfreeze.F90
src/clm1d/clm1d_thermal.F90
src/clm1d/clm1d_condch.F90
src/clm1d/clm1d_obuini.F90
src/clm1d/clm1d_thermalk.F90
src/clm1d/clm1d_condcq.F90
src/clm1d/clm1d_obult.F90
src/clm1d/clm1d_coszen.F90
src/clm1d/clm1d_psi.F90
src/clm1d/clm1d_tridia.F90
src/clm1d/clm1d_dynvegpar.F90
src/clm1d/clm1d_qsadv.F90
src/clm1d/clm1d_twostream.F90
src/clm1d/clm1d_hydro_canopy.F90
src/clm1d/clm1d_snowage.F90
src/clm1d/clm1d_hydro_irrig.F90
src/clm1d/clm1d_snowalb.F90
src/clm1d/clm1d_hydro_snow.F90
src/clm1d/clm1d_soilalb.F90
src/clm1d/clm1d_hydro_soil.F90
src/clm1d/clm1d_hydro_wetice.F90
src/clm1d/clm1d_stomata.F90
src/drv/grid_type.F90
src/drv/tile_type.F90
src/drv/io_type.F90
src/drv/drv_type.F90
src/drv/drv_almaout.F90
src/drv/drv_g2clm.F90
src/drv/drv_readclmin.F90
src/drv/drv_t2g.F90
src/drv/drv_astp.F90
src/drv/drv_getforce.F90
src/drv/drv_readvegpf.F90
src/drv/drv_tick.F90
src/drv/drv_1dout.F90
src/drv/drv_clm2g.F90
src/drv/drv_gridave.F90
src/drv/drv_g2tile.F90
src/drv/drv_2dout.F90
src/drv/drv_clmini.F90
src/drv/drv_pout.F90
src/drv/drv_restart.F90
src/${HOST_SOURCE}.F90
src/clm_type.F90
src/clm_host_transfer.F90
)

include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

add_library(${CLM_LIB} ${SRC_FILES})

install (TARGETS ${CLM_LIB} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )
