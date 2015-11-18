GSLINC = $(GSLROOT)/include

BAMBILIB = $(BAMBIROOT)/lib
BAMBIINCLUDE = $(BAMBIROOT)/include

CFLAGS := -O2 -std=c99 -I$(GSLINC) -I../tools -I../integration -I../EOBNRv2HMROM -I../LISAsim -I../LLVsim -I../LLVinference
CPPFLAGS := -O2 -I$(GSLINC) -I../tools -I../integration -I../EOBNRv2HMROM -I../LISAsim -I../LLVsim -I../LLVinference

OBJ = LLVinference.o LLVutils.o bambi.o


all: $(OBJ) LLVinference

LLVutils.o: LLVutils.c LLVutils.h
	$(CC) -c $(CFLAGS) LLVutils.c

bambi.o: bambi.cc bambi.h
	$(CPP) -c $(CPPFLAGS) -I$(BAMBIINCLUDE) bambi.cc

LLVinference.o: LLVinference.c LLVinference.h ../LLVsim/LLVFDresponse.h ../LLVsim/LLVnoise.h ../LLVsim/LLVgeometry.h ../tools/constants.h ../tools/struct.h ../tools/likelihood.h ../EOBNRv2HMROM/EOBNRv2HMROM.h ../EOBNRv2HMROM/EOBNRv2HMROMstruct.h ../integration/wip.h
	$(CC) -c $(CFLAGS) -I$(BAMBIINCLUDE) LLVinference.c

LLVinference: LLVinference.o LLVutils.o bambi.o ../LLVsim/LLVFDresponse.o ../LLVsim/LLVnoise.o ../tools/struct.o ../tools/timeconversion.o ../tools/likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o
	$(LD) $(LDFLAGS) -o LLVinference LLVinference.o LLVutils.o bambi.o ../LLVsim/LLVFDresponse.o ../LLVsim/LLVnoise.o ../tools/struct.o ../tools/timeconversion.o ../tools/likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o -L$(GSLROOT)/lib -L$(BAMBILIB) -lgsl -lgslcblas -lm -lbambi-1.2 $(MPILIBS)

phaseSNR.o: phaseSNR.c LLVinference.h ../LLVsim/LLVFDresponse.h ../LLVsim/LLVnoise.h ../LLVsim/LLVgeometry.h ../tools/constants.h ../tools/struct.h ../tools/likelihood.h ../EOBNRv2HMROM/EOBNRv2HMROM.h ../EOBNRv2HMROM/EOBNRv2HMROMstruct.h ../integration/wip.h
	$(CC) -c $(CFLAGS) phaseSNR.c

phaseSNR: phaseSNR.o LLVutils.o ../LLVsim/LLVFDresponse.o ../LLVsim/LLVnoise.o ../tools/struct.o ../tools/timeconversion.o ../tools/likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o
	$(CC) $(CPPFLAGS) -o phaseSNR phaseSNR.o LLVutils.o ../LLVsim/LLVFDresponse.o ../LLVsim/LLVnoise.o ../tools/struct.o ../tools/timeconversion.o ../tools/likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o -L$(GSLROOT)/lib -lgsl -lgslcblas -lm $(MPILIBS)	

findDist.o: findDist.c LLVinference.h ../LLVsim/LLVFDresponse.h ../LLVsim/LLVnoise.h ../LLVsim/LLVgeometry.h ../tools/constants.h ../tools/struct.h ../tools/likelihood.h ../EOBNRv2HMROM/EOBNRv2HMROM.h ../EOBNRv2HMROM/EOBNRv2HMROMstruct.h ../integration/wip.h
	$(CC) -c $(CFLAGS) findDist.c

findDist: findDist.o LLVutils.o ../LLVsim/LLVFDresponse.o ../LLVsim/LLVnoise.o ../tools/struct.o ../tools/timeconversion.o ../tools/likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o
	$(CC) $(CPPFLAGS) -o findDist findDist.o LLVutils.o ../LLVsim/LLVFDresponse.o ../LLVsim/LLVnoise.o ../tools/struct.o ../tools/timeconversion.o ../tools/likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o -L$(GSLROOT)/lib -lgsl -lgslcblas -lm $(MPILIBS)	

clean:
	-rm *.o
