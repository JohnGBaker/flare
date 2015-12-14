GSLROOT = /opt/local
GSLINC = $(GSLROOT)/include
BAMBIROOT = $(HOME)/build/bambi
BAMBIINC = $(BAMBIROOT)/include
BAMBILIB = $(BAMBIROOT)/lib
CC = gcc
CPP = g++
#Uncomment this for MPI and specify your needed MPI libraries
CC += -DPARALLEL
CPP += -DPARALLEL
MPILIBS = -lmpi -lmpi_cxx -lmpi_mpifh

CFLAGS += -O2 -std=c99 -I$(GSLINC)
CPPFLAGS += -O2 -I$(GSLINC)

SUBDIRS = tools EOBNRv2HMROM LISAsim LLVsim integration LISAinference LLVinference
SUBCLEAN = $(addsuffix .clean,$(SUBDIRS))

export CC CPP GSLROOT GSLINC BAMBIROOT BAMBIINC BAMBILIB MPILIBS CFLAGS CPPFLAGS

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
