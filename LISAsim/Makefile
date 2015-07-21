GSLINC = $(GSLROOT)/include
CFLAGS := -O2 -std=c99 -I$(GSLINC) -I../tools -I../integration -I../EOBNRv2HMROM -I../LISAsim -I../LLVsim -I../LLVinference

OBJ = LISAFDresponse.o LISAgeometry.o #LISAexampleSNR.o


all: $(OBJ) #LISAexampleSNR

LISAgeometry.o: LISAgeometry.c LISAgeometry.h ../tools/constants.h 
	$(CC) -c $(CFLAGS) LISAgeometry.c

LISAFDresponse.o: LISAFDresponse.c  LISAFDresponse.h LISAgeometry.h ../tools/constants.h ../tools/struct.h 
	$(CC) -c $(CFLAGS) LISAFDresponse.c

#LISAexampleSNR.o: LISAexampleSNR.c LISAFDresponse.h LISAgeometry.h ../tools/constants.h ../tools/struct.h ../tools/likelihood.h ../EOBNRv2HMROM/EOBNRv2HMROM.h ../EOBNRv2HMROM/EOBNRv2HMROMstruct.h ../integration/wip.h
#	$(CC) -c $(CFLAGS) LISAexampleSNR.c

#LISAexampleSNR: LISAexampleSNR.o LISAFDresponse.o LISAgeometry.o ../tools/struct.o ../tools/likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o
#	$(CC) $(CFLAGS) -lgsl -lgslcblas -lm -o LISAexampleSNR LISAexampleSNR.o LISAFDresponse.o LISAgeometry.o ../tools/struct.o ../tools/likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o -lgsl -L$(GSLROOT)/lib

clean:
	-rm *.o
