GSLROOT = /usr/local/other/SLES11.1/gsl/1.16/intel-13.0.1.117
#GSLROOT = /opt/local
BAMBIROOT = /discover/nobackup/jgbaker/sw/bambi/
#BAMBIROOT = $(HOME)/build/bambi
CC = icc
CPP = icpc
LD = $(CPP)
LDFLAGS =
#Uncomment this for MPI and specify your needed MPI libraries
CC += -DPARALLEL
CPP += -DPARALLEL
MPILIBS = -lmpi -lmpi_cxx -lmpi_mpifh

#Need to add/change the following to build on discover
GSLROOT = /usr/local/other/SLES11.1/gsl/1.15/intel-12.1.0.233
FC = mpif90 -DPARALLEL
CC = mpicc -DPARALLEL
CXX = mpiicpc -DPARALLEL
LAPACKLIB = -lmkl_intel_lp64 -lmkl_intel_thread -liomp5 -lpthread -lm -lmkl_core -lmkl_lapack95_lp64
#LAPACKLIB = -L/usr/local/other/SLES11/lapack/3.3.1/intel-12.1.0.233/lib -llapack -lblas -lm
MPILIBS += $(LAPACKLIB)
#MPILIBS += -lmpi -lintlc -limf -lsvml -lmkl_sequential -lmkl_intel_ilp64 -lmkl_core
LD = mpif90
LDFLAGS += -cxxlib -nofor_main -g -traceback -C
#end discover

GSLINC = $(GSLROOT)/include
CFLAGS += -g -C -O0 -traceback -std=c99 -I$(GSLINC) -I./tools -I./EOBNRv2HMROM -I./integration -I./LISAsim -I./LLVsim -I./LLVinference

SUBDIRS = tools EOBNRv2HMROM LISAsim LLVsim integration LISAinference LLVinference
SUBCLEAN = $(addsuffix .clean,$(SUBDIRS))

export CC CPP GSLROOT BAMBIROOT MPILIBS LD LDFLAGS

.PHONY: all clean subdirs $(SUBDIRS)

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

EOBNRv2HMROM: tools

LISAsim: tools integration EOBNRv2HMROM

LLVsim: tools integration EOBNRv2HMROM

LISAinference: tools integration EOBNRv2HMROM LISAsim

LLVinference: tools integration EOBNRv2HMROM LLVsim

phaseSNR: tools integration EOBNRv2HMROM LLVsim
	$(MAKE) -C LLVinference phaseSNR

all: subdirs

clean: $(SUBCLEAN)

$(SUBCLEAN): %.clean:
	$(MAKE) -C $* clean
