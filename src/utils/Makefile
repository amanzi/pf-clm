SRCDIR	      =	.
OBJECT        = ..

default: all

include $(OBJECT)/config/Makefile.config

# Order is important due to .mod file creation
OBJS = \
	clm_precision.o \
	clm_infnan.o \
	clm_abort.o \
	clm_io_config.o


all: $(OBJS)

clean:
	$(CLM_CLEAN)

allclean: clean

#include $(SRCDIR)/Makefile.depend



