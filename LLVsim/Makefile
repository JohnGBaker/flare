GSLROOT = /opt/local
GSLINC = $(GSLROOT)/include
CC = gcc
CFLAGS := -g -std=c99 -I$(GSLINC) -I../tools -I../integration -I../EOBNRv2HMROM -I../LISAsim -I../LLVsim

OBJ = LLVFDresponse.o LLVnoise.o LLVexampleSNR.o


all: $(OBJ) LLVexampleSNR

LLVFDresponse.o: LLVFDresponse.c  LLVFDresponse.h LLVgeometry.h ../tools/constants.h ../tools/struct.h 
	$(CC) -c $(CFLAGS) LLVFDresponse.c

LLVnoise.o: LLVnoise.c LLVnoise.h ../tools/constants.h 
	$(CC) -c $(CFLAGS) LLVnoise.c

LLVexampleSNR.o: LLVexampleSNR.c LLVFDresponse.h LLVnoise.h LLVgeometry.h ../tools/constants.h ../tools/struct.h ../tools/likelihood.h ../EOBNRv2HMROM/EOBNRv2HMROM.h ../EOBNRv2HMROM/EOBNRv2HMROMstruct.h ../integration/wip.h
	$(CC) -c $(CFLAGS) LLVexampleSNR.c

LLVexampleSNR: LLVexampleSNR.o LLVFDresponse.o LLVnoise.o ../tools/struct.o ../tools/likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o
	$(CC) $(CFLAGS) -lgsl -lgslcblas -lm -o LLVexampleSNR LLVexampleSNR.o LLVFDresponse.o LLVnoise.o ../tools/struct.o ../tools/likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o -lgsl -L$(GSLROOT)/lib

clean:
	-rm *.o
