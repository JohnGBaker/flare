MACHINE=discover
ifeq ($MACHINE,sylvainsmac)
  GSLROOT = /opt/local
  BAMBIROOT = $(HOME)/build/bambi
  CC = gcc
  CPP = g++
  LD = $(CPP)	
  LDFLAGS=
  #Uncomment this for MPI and specify your needed MPI libraries
  CC += -DPARALLEL
  CPP += -DPARALLEL
  MPILIBS = -lmpi -lmpi_cxx -lmpi_mpifh
else if ($MACHINE,discover) 
  #based on modules:
  #module load comp/intel-15.0.3.187 lib/mkl-15.0.3.187 mpi/impi-5.0.3.048
  GSLROOT = /usr/local/other/SLES11.1/gsl/1.16/intel-13.0.1.117
  BAMBIROOT = /discover/nobackup/jgbaker/sw/bambi/
  FC = mpif90 -DPARALLEL
  CC = mpicc -DPARALLEL
  CPP = mpiicpc -DPARALLEL
  LAPACKLIB = -lmkl_intel_lp64 -lmkl_intel_thread -liomp5 -lpthread -lm -lmkl_core -lmkl_lapack95_lp64
  MPILIBS += $(LAPACKLIB)
  LD = mpif90
  LDFLAGS = -cxxlib -nofor_main -g -traceback -C
endif

GSLINC = $(GSLROOT)/include
CFLAGS += -O2 -std=c99 -I$(GSLINC) -I./tools -I./EOBNRv2HMROM -I./integration -I./LISAsim -I./LLVsim -I./LLVinference

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
