SRCDIR	      =	.
OBJECT        = ..

default: all

include $(OBJECT)/config/Makefile.config
FCFLAGS += -I ../$(CLM_UTILS_DIR) -I ../$(CLM_CLM1D_DIR)

# Order is important due to .mod file creation
OBJS = \
	grid_type.o \
	tile_type.o \
	io_type.o \
	drv_type.o \
	drv_almaout.o \
	drv_g2clm.o \
	drv_readclmin.o \
	drv_t2g.o \
	drv_astp.o \
	drv_getforce.o \
	drv_readvegpf.o \
	drv_tick.o \
	drv_1dout.o \
	drv_clm2g.o \
	drv_gridave.o \
	drv_readvegtf.o \
	drv_2dout.o \
	drv_clmini.o \
	drv_pout.o \
	drv_restart.o \
	drv_open_files.o \
	drv_close_files.o


all: $(OBJS)
	$(AR) rus drv.a $(OBJS)

clean:
	$(CLM_CLEAN)

allclean: clean

#include $(SRCDIR)/Makefile.depend
