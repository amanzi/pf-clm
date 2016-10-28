SRCDIR	      =	.
OBJECT        = .


HOST_SOURCE   = clm_host_structured
#HOST_SOURCE   = clm_host_unstructured

default: all

include $(OBJECT)/config/Makefile.config
FCFLAGS += -I $(CLM_UTILS_DIR) -I $(CLM_CLM1D_DIR) -I $(CLM_DRV_DIR)

# Order is important due to .mod file creation
OBJS = \
	utils/clm_precision.o \
	utils/clm_infnan.o \
	utils/clm_abort.o \
	utils/clm_io_config.o \
	clm1d/clm1d_varpar.o\
	clm1d/clm1d_varcon.o\
	clm1d/clm1d_type.o\
	clm1d/clm1d_typini.o\
	clm1d/clm1d_balchk.o\
	clm1d/clm1d_lake.o\
	clm1d/clm1d_subdiv.o\
	clm1d/clm1d_combin.o\
	clm1d/clm1d_leaftem.o\
	clm1d/clm1d_surfalb.o\
	clm1d/clm1d_combo.o\
	clm1d/clm1d_main.o\
	clm1d/clm1d_surfrad.o\
	clm1d/clm1d_compact.o\
	clm1d/clm1d_meltfreeze.o\
	clm1d/clm1d_thermal.o\
	clm1d/clm1d_condch.o\
	clm1d/clm1d_obuini.o\
	clm1d/clm1d_thermalk.o\
	clm1d/clm1d_condcq.o\
	clm1d/clm1d_obult.o\
	clm1d/clm1d_coszen.o\
	clm1d/clm1d_psi.o\
	clm1d/clm1d_tridia.o\
	clm1d/clm1d_dynvegpar.o\
	clm1d/clm1d_qsadv.o\
	clm1d/clm1d_twostream.o\
	clm1d/clm1d_hydro_canopy.o\
	clm1d/clm1d_snowage.o\
	clm1d/clm1d_hydro_irrig.o\
	clm1d/clm1d_snowalb.o\
	clm1d/clm1d_hydro_snow.o\
	clm1d/clm1d_soilalb.o\
	clm1d/clm1d_hydro_soil.o\
	clm1d/clm1d_hydro_wetice.o\
	clm1d/clm1d_stomata.o \
	drv/grid_type.o \
	drv/tile_type.o \
	drv/io_type.o \
	drv/drv_type.o \
	drv/drv_almaout.o \
	drv/drv_g2clm.o \
	drv/drv_readclmin.o \
	drv/drv_t2g.o \
	drv/drv_astp.o \
	drv/drv_getforce.o \
	drv/drv_readvegpf.o \
	drv/drv_tick.o \
	drv/drv_1dout.o \
	drv/drv_clm2g.o \
	drv/drv_gridave.o \
	drv/drv_g2tile.o \
	drv/drv_2dout.o \
	drv/drv_clmini.o \
	drv/drv_pout.o \
	drv/drv_restart.o \
	$(HOST_SOURCE).o \
	clm_type.o \
	clm_host_transfer.o

all: $(OBJS)
	$(AR) cr clm.a $(OBJS)

clean:
	@$(CLM_CLEAN)

allclean: clean
	@cd utils; $(CLM_CLEAN)
	@cd clm1d; $(CLM_CLEAN)
	@cd drv; $(CLM_CLEAN)

