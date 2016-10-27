SRCDIR	      =	.
OBJECT        = .


HOST_SOURCE   = clm_host_structured
#HOST_SOURCE   = clm_host_unstructured

default: all

include $(OBJECT)/config/Makefile.config
FCFLAGS += -I $(CLM_UTILS_DIR) -I $(CLM_CLM1D_DIR) -I $(CLM_DRV_DIR)

# Order is important due to .mod file creation
OBJS = \
	$(HOST_SOURCE).o \
	clm_type.o \
	clm_host_transfer.o
#	clm.o


all: subdirs $(OBJS)

subdirs:
	@cd utils; make
	@cd clm1d; make
	@cd drv; make

clean:
	$(CLM_CLEAN)

allclean: clean
	@cd utils; make allclean
	@cd clm1d; make allclean
	@cd drv; make allclean


#include $(SRCDIR)/Makefile.depend
