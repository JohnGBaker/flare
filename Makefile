MESSAGE="Specify which machine to compile for in the Makefile."
MACHINE="johnsmac"

ifeq ($(MACHINE),"sylvainsmac")
  MESSAGE="Compiling for Sylvain's Mac"
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
else ifeq ($(MACHINE),"johnsmac")
  MESSAGE="Compiling for John's Mac"
  GSLROOT = /opt/local
  BAMBIROOT = ../BAMBI
  CC = gcc-mp-4.7 -g
  CXX = g++-mp-4.7
  CPP = g++-mp-4.7
  LD = $(CPP)
  #LD = gfortran-mp-4.7	
  LDFLAGS= -fopenmp -L/opt/local/lib -L/opt/local/lib/mpich-mp -lgfortran -llapack -latlas -lblas
  #Uncomment this for MPI and specify your needed MPI libraries
  #CC += -DPARALLEL
  #CPP += -DPARALLEL
  MPILIBS = -lmpi -lmpicxx -lmpifort
  CFLAGS += -g -fopenmp -I/opt/local/include/mpich-mp
  CPPFLAGS += -g -I/opt/local/include/mpich-mp
  CXXFLAGS = -g -fopenmp
  PTMCMC=$(PWD)/ptmcmc
else ifeq ($(MACHINE),"discover") 
  #based on modules:
  #module load comp/intel-15.0.3.187 lib/mkl-15.0.3.187 mpi/impi-5.0.3.048
  MESSAGE="Compiling for Discover at NCCS"
  GSLROOT = /usr/local/other/SLES11.1/gsl/1.16/intel-13.0.1.117
  BAMBIROOT = /discover/nobackup/jgbaker/sw/bambi/
  FC = mpif90 -DPARALLEL
  CC = mpicc -DPARALLEL
  CPP = mpiicpc -DPARALLEL
  LAPACKLIB = -lmkl_intel_lp64 -lmkl_intel_thread -liomp5 -lpthread -lm -lmkl_core -lmkl_lapack95_lp64
  MPILIBS += $(LAPACKLIB)
  LD = mpif90
  LDFLAGS = -cxxlib -nofor_main -g -traceback -C
else ifeq ($(MACHINE),"datura") 
  #based on modules:
  #module add Compiler/intel/ips_xe_2015/ips_xe_2015_intel15 mpi/openmpi/1.10.0-intel15 hdf5/1.8.13-intel15 gsl/1.15
  #environment:
  #AEI_GSL_HOME=/cluster/gsl/SL6/1.15 AEI_MKLROOT=/cluster/Compiler/Intel/ips_xe_2015/composer_xe_2015.1.133/mkl
  #export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$AEI_MKLROOT/lib/intel64
  MESSAGE="Compiling for Datura at AEI"
  GSLROOT = $(AEI_GSL_HOME)
  MKLROOT = $(AEI_MKLROOT)
  MKLINC = $(MKLROOT)/include/
  #GSLINC = 3D -I$(GSLROOT)/include -DHAVE_INLINE -DGSL_C99_INLINE -DGSL_RANGE_CHECK_OFF
  BAMBIROOT = $(HOME)/build/bambi
  FC = mpif90 -DPARALLEL
  CC = mpicc -DPARALLEL
  CPP = mpic++ -DPARALLEL
  MPILIBS = -lmpi -lmpi_cxx -lmpi_mpifh
  IFORTLIB = -lifcore
  LAPACKLIB = -L$(AEI_MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -liomp5 -lpthread -lm -lmkl_core -lmkl_lapack95_lp64
  MPILIBS += $(IFORTLIB)
  MPILIBS += $(LAPACKLIB)
  LD = mpif90
  LDFLAGS += -cxxlib -nofor_main -g -traceback -C
  CFLAGS += -I$(MKLINC)
  CPPFLAGS += -I$(MKLINC)
endif

GSLINC = $(GSLROOT)/include
BAMBIINC = $(BAMBIROOT)/include
BAMBILIB = $(BAMBIROOT)/lib
CFLAGS += -std=c99
CPPFLAGS += -O2 -I$(GSLINC)
CXXFLAGS += -std=c++11

SUBDIRS = tools EOBNRv2HMROM LISAsim LLVsim integration
ifdef PTMCMC
  SUBDIRS += ptmcmc
  CXXFLAGS+= -I$(PTMCMC)/include
endif
SUBDIRS +=  LISAinference LLVinference

SUBCLEAN = $(addsuffix .clean,$(SUBDIRS))

export CC CPP CXX GSLROOT GSLINC BAMBIROOT BAMBIINC BAMBILIB MPILIBS CFLAGS CPPFLAGS LD LDFLAGS PTMCMC CXXFLAGS


.PHONY: all clean message subdirs $(SUBDIRS)


all: message subdirs

message:
	@echo $(MESSAGE)

subdirs: $(SUBDIRS)
	@echo Making in subdirs:

$(SUBDIRS):
	$(MAKE) -C $@

EOBNRv2HMROM: tools

LISAsim: tools integration EOBNRv2HMROM

LLVsim: tools integration EOBNRv2HMROM

LISAinference: tools integration EOBNRv2HMROM LISAsim ptmcmc

LLVinference: tools integration EOBNRv2HMROM LLVsim

phaseSNR: tools integration EOBNRv2HMROM LLVsim
	$(MAKE) -C LLVinference phaseSNR

ifdef PTMCMC
.ptmcmc-version: $(PTMCMC)/lib/libptmcmc.a $(PTMCMC)/lib/libprobdist.a
	cd ptmcmc;git rev-parse HEAD > ../.ptmcmc-version;git status >> ../.ptmcmc-version;git diff >> ../.ptmcmc-version
ptmcmc:
	@echo "Do we need to check out ptmcmc from github?:";\
	if [ \! -d ptmcmc ]; then git clone https://github.com/JohnGBaker/ptmcmc.git; fi;
	$(MAKE) CFLAGS="$(CPPFLAGS) $(CXXFLAGS)" -C ptmcmc
endif

clean: $(SUBCLEAN)

$(SUBCLEAN): %.clean:
	$(MAKE) -C $* clean
