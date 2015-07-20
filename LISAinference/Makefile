GSLROOT = /opt/local
GSLINC = $(GSLROOT)/include
BAMBIROOT = $(HOME)/build/bambi
BAMBILIB = $(BAMBIROOT)/lib
BAMBIINCLUDE = $(BAMBIROOT)/include
CC = gcc
CPP = g++
CFLAGS := -O2 -std=c99 -I$(GSLINC) -I../tools -I../integration -I../EOBNRv2HMROM -I../LISAsim -I../LLVsim -I../LISAinference
CPPFLAGS := -O2 -I$(GSLINC) -I../tools -I../integration -I../EOBNRv2HMROM -I../LISAsim -I../LLVsim -I../LISAinference

OBJ = LISAinference.o LISAutils.o bambi.o


all: $(OBJ) LISAinference

LISAutils.o: LISAutils.c LISAutils.h
	$(CC) -c $(CFLAGS) LISAutils.c

bambi.o: bambi.cc bambi.h
	mpic++ -c $(CPPFLAGS) -I$(BAMBIINCLUDE) -D PARALLEL bambi.cc

LISAinference.o: LISAinference.c LISAinference.h ../LISAsim/LISAFDresponse.h ../LISAsim/LISAnoise.h ../LISAsim/LISAgeometry.h ../tools/constants.h ../tools/struct.h ../tools/likelihood.h ../EOBNRv2HMROM/EOBNRv2HMROM.h ../EOBNRv2HMROM/EOBNRv2HMROMstruct.h ../integration/wip.h
	mpicc -c $(CFLAGS) -I$(BAMBIINCLUDE) -D PARALLEL LISAinference.c

LISAinference: LISAinference.o LISAutils.o bambi.o ../LISAsim/LISAFDresponse.o ../LISAsim/LISAnoise.o ../LISAsim/LISAgeometry.o ../tools/struct.o ../tools/timeconversion.o ../tools/likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o
	mpic++ $(CPPFLAGS) -o LISAinference LISAinference.o LISAutils.o bambi.o ../LISAsim/LISAFDresponse.o ../LISAsim/LISAnoise.o ../LISAsim/LISAgeometry.o ../tools/struct.o ../tools/timeconversion.o ../tools/likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o -L$(GSLROOT)/lib -L$(BAMBILIB) -lgsl -lgslcblas -lm -lbambi-1.2

clean:
	-rm *.o
