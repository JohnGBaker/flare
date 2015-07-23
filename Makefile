GSLROOT = /opt/local
BAMBIROOT = $(HOME)/build/bambi
CC = gcc
CPP = g++
### Uncomment this for MPI and specify your needed MPI libraries
#CC += -DPARALLEL
#CPP += -DPARALLEL
#MPILIBS = -lmpi -lmpi_cxx -lmpi_mpifh

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

all: subdirs

clean: $(SUBCLEAN)

$(SUBCLEAN): %.clean:
	$(MAKE) -C $* clean
