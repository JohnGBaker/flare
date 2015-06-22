GSLROOT = /opt/local
GSLINC = $(GSLROOT)/include
CC = gcc
CFLAGS := -g -O0 -std=c99 -I$(GSLINC) -I../tools -I../integration -I../EOBNRv2HMROM -I../LISAsim -I../LLVsim -I../LLVinference

OBJ = LLVinference.o LLVInit.o


all: $(OBJ) LLVinference

LLVInit.o: LLVInit.c LLVInit.h
	$(CC) -c $(CFLAGS) LLVInit.c

LLVinference.o: LLVinference.c LLVinference.h ../LLVsim/LLVFDresponse.h ../LLVsim/LLVnoise.h ../LLVsim/LLVgeometry.h ../tools/constants.h ../tools/struct.h ../tools/likelihood.h ../EOBNRv2HMROM/EOBNRv2HMROM.h ../EOBNRv2HMROM/EOBNRv2HMROMstruct.h ../integration/wip.h
	$(CC) -c $(CFLAGS) LLVinference.c

LLVinference: LLVinference.o ../LLVsim/LLVFDresponse.o ../LLVsim/LLVnoise.o ../tools/struct.o ../tools/likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o
	$(CC) $(CFLAGS) -lgsl -lgslcblas -lm -o LLVinference LLVinference.o ../LLVsim/LLVFDresponse.o ../LLVsim/LLVnoise.o ../tools/struct.o ../tools/likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o -lgsl -L$(GSLROOT)/lib

clean:
	-rm *.o
