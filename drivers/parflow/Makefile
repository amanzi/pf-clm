SRCDIR	      =	../../src/
OBJECT        = ../../src/

default: all

include $(OBJECT)/config/Makefile.config
FCFLAGS += -I$(SRCDIR) -I $(SRCDIR)$(CLM_UTILS_DIR) -I $(SRCDIR)$(CLM_CLM1D_DIR) -I $(SRCDIR)$(CLM_DRV_DIR)

# Order is important due to .mod file creation
OBJS = \
	parflow_open_files.o \
	parflow_close_files.o \
	parflow_read_ground.o \
	parflow_check_mass_balance.o \
	parflow_driver.o


all: subdirs $(OBJS)
	$(AR) cr parflow_clm.a $(OBJS) ../../src/clm.a


subdirs:
	@cd ../../src; make

clean:
	@$(CLM_CLEAN)

allclean: clean
	@cd ../../src; make allclean


#include $(SRCDIR)/Makefile.depend
