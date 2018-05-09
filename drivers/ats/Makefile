SRCDIR	      =	../../src/
OBJECT        = ../../src/

default: all

include $(OBJECT)/config/Makefile.config
FCFLAGS += -I$(SRCDIR) -I $(SRCDIR)$(CLM_UTILS_DIR) -I $(SRCDIR)$(CLM_CLM1D_DIR) -I $(SRCDIR)$(CLM_DRV_DIR) -fbounds-check -g3

# Order is important due to .mod file creation
OBJS = \
	ats_clm.o

all: subdirs $(OBJS)
	$(AR) cr libats_clm.a $(OBJS)

subdirs:
	@cd ../../src; make

clean:
	@$(CLM_CLEAN)

allclean: clean
	@cd ../../src; make allclean


#include $(SRCDIR)/Makefile.depend
