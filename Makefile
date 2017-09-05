MESSAGE="Specify which machine to compile for in the Makefile."
MACHINE="sylvainsmac"
#MACHINE="discover"
#MACHINE="johnsmac"

ifeq ($(MACHINE),"sylvainsmac")
  MESSAGE="Compiling for Sylvain's Mac"
  GSLROOT = /opt/local
  BAMBIROOT = $(HOME)/build/bambi
  CC = gcc-mp-5
  CPP = g++-mp-5
	CXX = g++-mp-5
  LD = $(CPP)
  LDFLAGS=  -L$(GSLROOT)/lib
  #Uncomment this for MPI and specify your needed MPI libraries
	CFLAGS += -I/usr/local/include -I/opt/local/include
	CPPFLAGS += -I/usr/local/include -I/opt/local/include
	#Required for openmp with gcc
	CFLAGS += -fopenmp
	CPPFLAGS += -fopenmp
	LDFLAGS += -fopenmp
  CC += -DPARALLEL
  CPP += -DPARALLEL
  MPILIBS = -lmpi -lmpi_cxx -lmpi_mpifh
	PTMCMC=$(PWD)/ptmcmc
else ifeq ($(MACHINE),"johnsmac")
  MESSAGE="Compiling for John's Mac"
  GSLROOT = /opt/local
  BAMBIROOT = ../BAMBI
  CC = gcc-mp-4.7 -g
  CXX = g++-mp-4.7
  CPP = g++-mp-4.7
  LD = $(CPP)
  #LD = gfortran-mp-4.7
  LDFLAGS= -fopenmp -L/opt/local/lib -L/opt/local/lib/mpich-gcc47 -lgfortran -llapack -latlas -lblas
  #LDFLAGS= -fopenmp -L/opt/local/lib -L/opt/local/lib/mpich-mp -lgfortran -llapack -latlas -lblas
  #Uncomment this for MPI and specify your needed MPI libraries
  MPILIBS = -lmpi -lmpicxx -lmpifort
  #CFLAGS += -g -fopenmp -I/opt/local/include/mpich-mp
  #CPPFLAGS += -g -I/opt/local/include/mpich-mp
  CFLAGS += -g -O3 -fopenmp -I/opt/local/include/mpich-gcc47
  CPPFLAGS += -g -O3 -I/opt/local/include/mpich-gcc47
  CXXFLAGS = -g -O3 -fopenmp
  #Uncomment this for MPI and specify your needed MPI libraries
  CFLAGS += -DPARALLEL
  CXXFLAGS += -DPARALLEL
  PTMCMC=$(PWD)/ptmcmc
else ifeq ($(MACHINE),"discover")
  #based on modules:
  #module load comp/intel-15.0.3.187 lib/mkl-15.0.3.187 mpi/impi-5.0.3.048
  MESSAGE="Compiling for Discover at NCCS"
  GSLROOT = /usr/local/other/SLES11.1/gsl/1.16/intel-13.0.1.117
  BAMBIROOT = /discover/nobackup/jgbaker/sw/bambi/
  FFTWROOT = /usr/local/other/SLES11.1/fftw/3.3.3/intel-12.1.0.233/
  FC = mpif90 -DPARALLEL -traceback
  CC = mpicc -DPARALLEL -traceback
  CXX = mpiicpc -DPARALLEL -traceback
  CPP = mpiicpc -DPARALLEL -traceback
  LAPACKLIB = -lmkl_intel_lp64 -lmkl_intel_thread -liomp5 -lpthread -lm -lmkl_core -lmkl_lapack95_lp64
  CFLAGS += -g -O3 -fopenmp -I $(GSLROOT)/include -I $(FFTWROOT)/include -L$(FFTWROOT)/lib -L$(GSLROOT)/lib
  CXXFLAGS = -g -O3 -fopenmp
  MPILIBS += $(LAPACKLIB)
  LD = mpiicpc
  LDFLAGS = -cxxlib -g -traceback -C -fopenmp -L$(FFTWROOT)/lib -L$(GSLROOT)/lib
  #LDFLAGS += -nofor_main
  LDFLAGS += -lifcore
  PTMCMC=$(PWD)/ptmcmc
  #LD = mpif90
  #LDFLAGS = -cxxlib -nofor_main -g -traceback -C
else ifeq ($(MACHINE),"datura")
  #based on modules:
  #module add Compiler/intel/ips_xe_2015/ips_xe_2015_intel15 mpi/openmpi/1.10.0-intel15 hdf5/1.8.13-intel15 gsl/1.15
  #environment:
  #AEI_GSL_HOME=/cluster/gsl/SL6/1.15 AEI_MKLROOT=/cluster/Compiler/Intel/ips_xe_2015/composer_xe_2015.1.133/mkl
  #AEI_FFTW_HOME=/cluster/fftw/3.3
  #export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$AEI_MKLROOT/lib/intel64
  MESSAGE="Compiling for Datura at AEI"
  GSLROOT = $(AEI_GSL_HOME)
  MKLROOT = $(AEI_MKLROOT)
  MKLINC = $(MKLROOT)/include/
	FFTWROOT = $(AEI_FFTW_HOME)
	FFTWINC = $(FFTWROOT)/include
	FFTWLIB = $(FFTWROOT)/lib
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
	LDFLAGS += -L$(FFTWLIB)
  LDFLAGS += -cxxlib -nofor_main -g -traceback -C -fopenmp
  CFLAGS += -I$(MKLINC) -I$(FFTWINC) -fopenmp
  CPPFLAGS += -I$(MKLINC) -I$(FFTWINC) -fopenmp
endif

GSLINC = $(GSLROOT)/include
#FFTWINC = $(FFTWROOT)/include
BAMBIINC = $(BAMBIROOT)/include
BAMBILIB = $(BAMBIROOT)/lib
CFLAGS += -std=c99
CPPFLAGS += -O3 -I$(GSLINC)
CXXFLAGS += -g -std=c++11 -O3 -I$(GSLINC)

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

ifdef PTMCMC
LISAinference: tools integration EOBNRv2HMROM LISAsim ptmcmc
else
LISAinference: tools integration EOBNRv2HMROM LISAsim
endif

LLVinference: tools integration EOBNRv2HMROM LLVsim

phaseSNR: tools integration EOBNRv2HMROM LLVsim
	$(MAKE) -C LLVinference phaseSNR

ifdef PTMCMC
.ptmcmc-version: $(PTMCMC)/lib/libptmcmc.a $(PTMCMC)/lib/libprobdist.a
	cd ptmcmc;git rev-parse HEAD > ../.ptmcmc-version;git status >> ../.ptmcmc-version;git diff >> ../.ptmcmc-version
ptmcmc:
	@echo "Do we need to check out ptmcmc from github?:";\
	if [ \! -d ptmcmc ]; then git clone https://github.com/JohnGBaker/ptmcmc.git; fi;
	@echo "INCLUDE="$(INCLUDE)
	@$(MAKE) CFLAGS="$(CPPFLAGS) $(CXXFLAGS)" INCLUDE="$(PTMCMC)/include" -C ptmcmc
endif

clean: $(SUBCLEAN)

ifdef PTMCMC
$(SUBCLEAN): %.clean:
	$(MAKE) -C $* clean
	$(MAKE) -C ptmcmc INCLUDE="$(PTMCMC)/include" clean
else
$(SUBCLEAN): %.clean:
	$(MAKE) -C $* clean
endif
