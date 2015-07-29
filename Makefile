GSLROOT = /usr
BAMBIROOT = ${HOME}/build/bambi
CC = mpicc
CPP = mpicxx
#Uncomment this for MPI and specify your needed MPI libraries
CC += -DPARALLEL
CPP += -DPARALLEL
MPILIBS = -lmpi -lmpi_f77

GSLINC = $(GSLROOT)/include
CFLAGS += -O2 -std=c99 -I$(GSLINC) -I./tools -I./EOBNRv2HMROM -I./integration -I./LISAsim -I./LLVsim -I./LLVinference

SUBDIRS = tools EOBNRv2HMROM LISAsim LLVsim integration LISAinference LLVinference
SUBCLEAN = $(addsuffix .clean,$(SUBDIRS))

export CC CPP GSLROOT BAMBIROOT MPILIBS

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

findDist: tools integration EOBNRv2HMROM LLVsim
	$(MAKE) -C LLVinference findDist

all: subdirs

clean: $(SUBCLEAN)

$(SUBCLEAN): %.clean:
	$(MAKE) -C $* clean
