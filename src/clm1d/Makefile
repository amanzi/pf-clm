SRCDIR	      =	.

OBJECT        = ..


default: all

include $(OBJECT)/config/Makefile.config
FCFLAGS += -I ../$(CLM_UTILS_DIR)

# Order is important due to .mod file creation
OBJS = \
	clm1d_varpar.o\
	clm1d_varcon.o\
	clm1d_type.o\
	clm1d_typini.o\
	clm1d_balchk.o\
	clm1d_lake.o\
	clm1d_subdiv.o\
	clm1d_combin.o\
	clm1d_leaftem.o\
	clm1d_surfalb.o\
	clm1d_combo.o\
	clm1d_main.o\
	clm1d_surfrad.o\
	clm1d_compact.o\
	clm1d_meltfreeze.o\
	clm1d_thermal.o\
	clm1d_condch.o\
	clm1d_obuini.o\
	clm1d_thermalk.o\
	clm1d_condcq.o\
	clm1d_obult.o\
	clm1d_coszen.o\
	clm1d_psi.o\
	clm1d_tridia.o\
	clm1d_dynvegpar.o\
	clm1d_qsadv.o\
	clm1d_twostream.o\
	clm1d_hydro_canopy.o\
	clm1d_snowage.o\
	clm1d_hydro_irrig.o\
	clm1d_snowalb.o\
	clm1d_hydro_snow.o\
	clm1d_soilalb.o\
	clm1d_hydro_soil.o\
	clm1d_hydro_wetice.o\
	clm1d_stomata.o

all: $(OBJS)
#	$(AR) rus clm1d.a $(OBJS)

clean:
	$(CLM_CLEAN)

allclean: clean

#include $(SRCDIR)/Makefile.depend



