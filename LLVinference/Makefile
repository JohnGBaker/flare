CFLAGS += -I../tools -I../integration -I../EOBNRv2HMROM -I../LISAsim -I../LLVsim -I../LLVinference
CPPFLAGS += -I../tools -I../integration -I../EOBNRv2HMROM -I../LISAsim -I../LLVsim -I../LLVinference

OBJ = LLVinference.o LLVutils.o bambi.o LLVlikelihood.o


all: $(OBJ) LLVinference LLVlikelihood

LLVutils.o: LLVutils.c LLVutils.h
	$(CC) -c $(CFLAGS) LLVutils.c

bambi.o: bambi.cc bambi.h
	@echo CPP=$(CPP)
	$(CPP) -c $(CPPFLAGS) -I$(BAMBIINC) bambi.cc

LLVlikelihood.o: LLVlikelihood.c LLVinference.h LLVutils.h ../LLVsim/LLVFDresponse.h ../LLVsim/LLVnoise.h ../LLVsim/LLVgeometry.h ../tools/constants.h ../tools/struct.h ../tools/likelihood.h ../EOBNRv2HMROM/EOBNRv2HMROM.h ../EOBNRv2HMROM/EOBNRv2HMROMstruct.h ../integration/wip.h
	$(CC) -c $(CFLAGS) LLVlikelihood.c

LLVinference.o: LLVinference.c LLVinference.h LLVutils.h ../LLVsim/LLVFDresponse.h ../LLVsim/LLVnoise.h ../LLVsim/LLVgeometry.h ../tools/constants.h ../tools/struct.h ../tools/likelihood.h ../EOBNRv2HMROM/EOBNRv2HMROM.h ../EOBNRv2HMROM/EOBNRv2HMROMstruct.h ../integration/wip.h
	$(CC) -c $(CFLAGS) -I$(BAMBIINC) LLVinference.c

LLVlikelihood: LLVlikelihood.o LLVutils.o ../LLVsim/LLVFDresponse.o ../LLVsim/LLVnoise.o ../tools/struct.o ../tools/waveform.o ../tools/timeconversion.o ../tools/splinecoeffs.o ../tools/fresnel.o ../tools/likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o
	$(LD) $(LDFLAGS) -o LLVlikelihood LLVlikelihood.o LLVutils.o ../LLVsim/LLVFDresponse.o ../LLVsim/LLVnoise.o ../tools/struct.o ../tools/waveform.o ../tools/timeconversion.o ../tools/splinecoeffs.o ../tools/fresnel.o ../tools/likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o -L$(GSLROOT)/lib -lgsl -lgslcblas -lm

LLVinference: LLVinference.o LLVutils.o bambi.o ../LLVsim/LLVFDresponse.o ../LLVsim/LLVnoise.o ../tools/struct.o ../tools/waveform.o ../tools/timeconversion.o ../tools/splinecoeffs.o ../tools/fresnel.o ../tools/likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o
	$(LD) $(LDFLAGS) -o LLVinference LLVinference.o LLVutils.o bambi.o ../LLVsim/LLVFDresponse.o ../LLVsim/LLVnoise.o ../tools/struct.o ../tools/waveform.o ../tools/timeconversion.o ../tools/splinecoeffs.o ../tools/fresnel.o ../tools/likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o -L$(GSLROOT)/lib -L$(BAMBILIB) -lgsl -lgslcblas -lm -lbambi-1.2 $(MPILIBS)

phaseSNR.o: phaseSNR.c LLVinference.h ../LLVsim/LLVFDresponse.h ../LLVsim/LLVnoise.h ../LLVsim/LLVgeometry.h ../tools/constants.h ../tools/struct.h ../tools/likelihood.h ../EOBNRv2HMROM/EOBNRv2HMROM.h ../EOBNRv2HMROM/EOBNRv2HMROMstruct.h ../integration/wip.h
	$(CC) -c $(CFLAGS) phaseSNR.c

phaseSNR: phaseSNR.o LLVutils.o ../LLVsim/LLVFDresponse.o ../LLVsim/LLVnoise.o ../tools/struct.o ../tools/waveform.o ../tools/timeconversion.o ../tools/likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o
	$(CC) $(CPPFLAGS) -o phaseSNR phaseSNR.o LLVutils.o ../LLVsim/LLVFDresponse.o ../LLVsim/LLVnoise.o ../tools/struct.o ../tools/waveform.o ../tools/timeconversion.o ../tools/likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o -L$(GSLROOT)/lib -lgsl -lgslcblas -lm $(MPILIBS)

findDist.o: findDist.c LLVinference.h ../LLVsim/LLVFDresponse.h ../LLVsim/LLVnoise.h ../LLVsim/LLVgeometry.h ../tools/constants.h ../tools/struct.h ../tools/likelihood.h ../EOBNRv2HMROM/EOBNRv2HMROM.h ../EOBNRv2HMROM/EOBNRv2HMROMstruct.h ../integration/wip.h
	$(CC) -c $(CFLAGS) findDist.c

findDist: findDist.o LLVutils.o ../LLVsim/LLVFDresponse.o ../LLVsim/LLVnoise.o ../tools/struct.o ../tools/waveform.o ../tools/timeconversion.o ../tools/likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o
	$(CC) $(CPPFLAGS) -o findDist findDist.o LLVutils.o ../LLVsim/LLVFDresponse.o ../LLVsim/LLVnoise.o ../tools/struct.o ../tools/waveform.o ../tools/timeconversion.o ../tools/likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o -L$(GSLROOT)/lib -lgsl -lgslcblas -lm $(MPILIBS)

clean:
	-rm *.o
