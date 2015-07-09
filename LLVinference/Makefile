GSLROOT = /opt/local
GSLINC = $(GSLROOT)/include
BAMBILIB = $(HOME)/build/bambi/lib
BAMBIINCLUDE = $(HOME)/build/bambi/include
CC = gcc
CPP = g++
CFLAGS := -g -O0 -std=c99 -I$(GSLINC) -I../tools -I../integration -I../EOBNRv2HMROM -I../LISAsim -I../LLVsim -I../LLVinference
CPPFLAGS := -g -O0 -I$(GSLINC) -I../tools -I../integration -I../EOBNRv2HMROM -I../LISAsim -I../LLVsim -I../LLVinference

OBJ = LLVinference.o LLVInit.o bambi.o


all: $(OBJ) LLVinference

LLVInit.o: LLVInit.c LLVInit.h
	$(CC) -c $(CFLAGS) LLVInit.c

bambi.o: bambi.cc bambi.h
	$(CPP) -c $(CPPFLAGS) -I$(BAMBIINCLUDE) bambi.cc

LLVinference.o: LLVinference.c LLVinference.h ../LLVsim/LLVFDresponse.h ../LLVsim/LLVnoise.h ../LLVsim/LLVgeometry.h ../tools/constants.h ../tools/struct.h ../tools/likelihood.h ../EOBNRv2HMROM/EOBNRv2HMROM.h ../EOBNRv2HMROM/EOBNRv2HMROMstruct.h ../integration/wip.h
	$(CC) -c $(CFLAGS) -I$(BAMBIINCLUDE) LLVinference.c

LLVinference: LLVinference.o LLVInit.o bambi.o ../LLVsim/LLVFDresponse.o ../LLVsim/LLVnoise.o ../tools/struct.o ../tools/timeconversion.o ../tools/likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o
	$(CPP) $(CPPFLAGS) -o LLVinference LLVinference.o LLVInit.o bambi.o ../LLVsim/LLVFDresponse.o ../LLVsim/LLVnoise.o ../tools/struct.o ../tools/timeconversion.o ../tools/likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o -L$(GSLROOT)/lib -L$(BAMBILIB) -lgsl -lgslcblas -lm -lbambi-1.2

clean:
	-rm *.o
