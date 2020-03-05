cmake_minimum_required(VERSION 3.9)

set(CLM_LIB clm)
project(${CLM_LIB} LANGUAGES Fortran)

option(UNSTRUCTURED_HOST "ON/OFF builds CLM library for unstructured/structured host [ON]." ON)

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}" CACHE PATH "default install path" FORCE )
endif()

if (UNSTRUCTURED_HOST)
  set(HOST_SOURCE clm_host_unstructured)
else()
  set(HOST_SOURCE clm_host_structured)
endif()

set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/include )

set (SRC_FILES
utils/clm_precision.F90
utils/clm_infnan.F90
utils/clm_abort.F90
utils/clm_io_config.F90
clm1d/clm1d_varpar.F90
clm1d/clm1d_varcon.F90
clm1d/clm1d_type.F90
clm1d/clm1d_typini.F90
clm1d/clm1d_balchk.F90
clm1d/clm1d_lake.F90
clm1d/clm1d_subdiv.F90
clm1d/clm1d_combin.F90
clm1d/clm1d_leaftem.F90
clm1d/clm1d_surfalb.F90
clm1d/clm1d_combo.F90
clm1d/clm1d_main.F90
clm1d/clm1d_surfrad.F90
clm1d/clm1d_compact.F90
clm1d/clm1d_meltfreeze.F90
clm1d/clm1d_thermal.F90
clm1d/clm1d_condch.F90
clm1d/clm1d_obuini.F90
clm1d/clm1d_thermalk.F90
clm1d/clm1d_condcq.F90
clm1d/clm1d_obult.F90
clm1d/clm1d_coszen.F90
clm1d/clm1d_psi.F90
clm1d/clm1d_tridia.F90
clm1d/clm1d_dynvegpar.F90
clm1d/clm1d_qsadv.F90
clm1d/clm1d_twostream.F90
clm1d/clm1d_hydro_canopy.F90
clm1d/clm1d_snowage.F90
clm1d/clm1d_hydro_irrig.F90
clm1d/clm1d_snowalb.F90
clm1d/clm1d_hydro_snow.F90
clm1d/clm1d_soilalb.F90
clm1d/clm1d_hydro_soil.F90
clm1d/clm1d_hydro_wetice.F90
clm1d/clm1d_stomata.F90
drv/grid_type.F90
drv/tile_type.F90
drv/io_type.F90
drv/drv_type.F90
drv/drv_almaout.F90
drv/drv_g2clm.F90
drv/drv_readclmin.F90
drv/drv_t2g.F90
drv/drv_astp.F90
drv/drv_getforce.F90
drv/drv_readvegpf.F90
drv/drv_tick.F90
drv/drv_1dout.F90
drv/drv_clm2g.F90
drv/drv_gridave.F90
drv/drv_g2tile.F90
drv/drv_2dout.F90
drv/drv_clmini.F90
drv/drv_pout.F90
drv/drv_restart.F90
${HOST_SOURCE}.F90
clm_type.F90
clm_host_transfer.F90
)

include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

add_library(${CLM_LIB} ${SRC_FILES})

install (TARGETS ${CLM_LIB} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )
