TRILINOS_DIR = /Users/uec/codes/ats/amanzi-tpls/install-dev/trilinos-12-12-1/

SRCDIR	      =	../../src/
OBJECT        = ../../src/

#CXX=/usr/local/bin/g++-7
#CC=/usr/local/bin/g++-7
CXX=`which mpicxx`
CC=`which mpicxx`
FC=`which mpif90`
CXXFLAGS += -std=c++11 -g -I$(TRILINOS_DIR)/include

default: all


include $(OBJECT)/config/Makefile.config
FCFLAGS += -I$(SRCDIR) -I $(SRCDIR)$(CLM_UTILS_DIR) -I $(SRCDIR)$(CLM_CLM1D_DIR) -I $(SRCDIR)$(CLM_DRV_DIR) -fbounds-check -g3 -funderscoring


# Order is important due to .mod file creation
OBJS = \
	ats_clm.o \
	ats_clm_interface.o

all: subdirs $(OBJS)
	$(AR) cr libats_clm.a $(OBJS)

subdirs:
	@cd ../../src; make

clean:
	@$(CLM_CLEAN)

allclean: clean
	@cd ../../src; make allclean

test: all main.o 
	$(FC) $(FFLAGS) -o run_test main.o -L. -lats_clm -L../../src -lclm -L$(TRILINOS_DIR)/lib -lepetra -lteuchoscore -lteuchoscomm -lteuchoskokkoscompat -lteuchosparameterlist -lkokkoscore -lblas -lc++

#include $(SRCDIR)/Makefile.depend
