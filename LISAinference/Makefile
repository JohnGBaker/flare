CFLAGS +=-I../tools -I../integration -I../EOBNRv2HMROM -I../LISAsim -I../LLVsim -I../LISAinference
CPPFLAGS +=-I../tools -I../integration -I../EOBNRv2HMROM -I../LISAsim -I../LLVsim -I../LISAinference


ifdef PTMCMC
all: LISAinference ComputeLISASNR TestLISAData LISAlikelihood LISAinference_ptmcmc
else
all: LISAinference ComputeLISASNR TestLISAData LISAlikelihood
endif

LISAutils.o: LISAutils.c LISAutils.h ../LISAsim/LISAFDresponse.h ../LISAsim/LISAnoise.h ../LISAsim/LISAgeometry.h ../tools/constants.h ../tools/struct.h ../tools/waveform.h ../tools/fresnel.h ../tools/splinecoeffs.h ../tools/likelihood.h ../EOBNRv2HMROM/EOBNRv2HMROM.h ../EOBNRv2HMROM/EOBNRv2HMROMstruct.h ../integration/wip.h
	$(CC) -c $(CFLAGS) LISAutils.c

bambi.o: bambi.cc bambi.h
	@echo CPP=$(CPP)
	$(CPP) -c $(CPPFLAGS) -I$(BAMBIINC) bambi.cc

TestLISAData.o: TestLISAData.c ComputeLISASNR.h LISAutils.h ../LISAsim/LISAFDresponse.h ../LISAsim/LISAnoise.h ../LISAsim/LISAgeometry.h ../tools/constants.h ../tools/struct.h ../tools/waveform.h ../tools/fft.h ../tools/fresnel.h ../tools/splinecoeffs.h ../tools/likelihood.h ../EOBNRv2HMROM/EOBNRv2HMROM.h ../EOBNRv2HMROM/EOBNRv2HMROMstruct.h ../integration/wip.h
	$(CC) -c $(CFLAGS) -I$(BAMBIINC) TestLISAData.c

TestLISAData: TestLISAData.o LISAutils.o ../LISAsim/LISAFDresponse.o ../LISAsim/LISAnoise.o ../LISAsim/LISAgeometry.o ../tools/struct.o ../tools/waveform.o ../tools/fft.o ../tools/timeconversion.o ../tools/splinecoeffs.o ../tools/fresnel.o ../tools/likelihood.o ../tools/data_likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -L$(GSLROOT)/lib -lgsl -lgslcblas -lm -lfftw3

ComputeLISASNR.o: ComputeLISASNR.c ComputeLISASNR.h LISAutils.h ../LISAsim/LISAFDresponse.h ../LISAsim/LISAnoise.h ../LISAsim/LISAgeometry.h ../tools/constants.h ../tools/struct.h ../tools/waveform.h ../tools/fft.h ../tools/fresnel.h ../tools/splinecoeffs.h ../tools/likelihood.h ../EOBNRv2HMROM/EOBNRv2HMROM.h ../EOBNRv2HMROM/EOBNRv2HMROMstruct.h ../integration/wip.h
	$(CC) -c $(CFLAGS) -I$(BAMBIINC) ComputeLISASNR.c

ComputeLISASNR: ComputeLISASNR.o LISAutils.o ../LISAsim/LISAFDresponse.o ../LISAsim/LISAnoise.o ../LISAsim/LISAgeometry.o ../tools/struct.o ../tools/waveform.o ../tools/fft.o ../tools/timeconversion.o ../tools/splinecoeffs.o ../tools/fresnel.o ../tools/likelihood.o ../tools/data_likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -L$(GSLROOT)/lib -lgsl -lgslcblas -lm -lfftw3

LISAinference_common.o: LISAinference_common.c LISAutils.h ../LISAsim/LISAFDresponse.h ../LISAsim/LISAnoise.h ../LISAsim/LISAgeometry.h ../tools/constants.h ../tools/struct.h ../tools/fresnel.h ../tools/splinecoeffs.h ../tools/likelihood.h ../EOBNRv2HMROM/EOBNRv2HMROM.h ../EOBNRv2HMROM/EOBNRv2HMROMstruct.h ../integration/wip.h
		$(CC) -c $(CFLAGS) -I$(BAMBIINC) LISAinference_common.c

LISAlikelihood.o: LISAlikelihood.c LISAutils.h LISAinference_common.h ../LISAsim/LISAFDresponse.h ../LISAsim/LISAnoise.h ../LISAsim/LISAgeometry.h ../tools/constants.h ../tools/struct.h ../tools/fresnel.h ../tools/splinecoeffs.h ../tools/likelihood.h ../EOBNRv2HMROM/EOBNRv2HMROM.h ../EOBNRv2HMROM/EOBNRv2HMROMstruct.h ../integration/wip.h
		$(CC) -c $(CFLAGS) LISAlikelihood.c

LISAinference.o: LISAinference.c LISAinference.h LISAutils.h ../LISAsim/LISAFDresponse.h ../LISAsim/LISAnoise.h ../LISAsim/LISAgeometry.h ../tools/constants.h ../tools/struct.h ../tools/fresnel.h ../tools/splinecoeffs.h ../tools/likelihood.h ../EOBNRv2HMROM/EOBNRv2HMROM.h ../EOBNRv2HMROM/EOBNRv2HMROMstruct.h ../integration/wip.h
	$(CC) -c $(CFLAGS) -I$(BAMBIINC) LISAinference.c

LISAinference: LISAinference.o LISAutils.o bambi.o ../LISAsim/LISAFDresponse.o ../LISAsim/LISAnoise.o ../LISAsim/LISAgeometry.o ../tools/struct.o ../tools/waveform.o ../tools/timeconversion.o ../tools/splinecoeffs.o ../tools/fresnel.o ../tools/likelihood.o ../tools/data_likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o LISAinference_common.o
	$(LD) $(LDFLAGS) -o $@ $^ -L$(GSLROOT)/lib -L$(BAMBILIB) -lgsl -lgslcblas -lm -lbambi-1.2 $(MPILIBS)

LISAlikelihood: LISAlikelihood.o LISAutils.o ../LISAsim/LISAFDresponse.o ../LISAsim/LISAnoise.o ../LISAsim/LISAgeometry.o ../tools/struct.o ../tools/waveform.o ../tools/timeconversion.o ../tools/splinecoeffs.o ../tools/fresnel.o ../tools/likelihood.o ../tools/data_likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o LISAinference_common.o
	$(LD) $(LDFLAGS) -o $@ $^ -L$(GSLROOT)/lib -lgsl -lgslcblas -lm  $(MPILIBS)

ifdef PTMCMC
LISAinference_ptmcmc.o:  LISAinference_ptmcmc.cc  $(PTMCMC)/lib/libptmcmc.a

LISAinference_ptmcmc:  LISAinference_ptmcmc.o LISAinference_common.o LISAutils.o ../LISAsim/LISAFDresponse.o ../LISAsim/LISAnoise.o ../LISAsim/LISAgeometry.o ../tools/struct.o ../tools/waveform.o ../tools/timeconversion.o ../tools/splinecoeffs.o ../tools/fresnel.o ../tools/likelihood.o ../tools/data_likelihood.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../integration/wip.o ../integration/spline.o ../integration/Faddeeva.o $(PTMCMC)/lib/libptmcmc.a
	@echo $(LD)
	$(LD) $(LDFLAGS) -o $@ $^ -L$(GSLROOT)/lib -L$(PTMCMC)/lib -lgsl -lgslcblas -lm -lptmcmc -lprobdist -I$(GSLROOT)/include -I$(PTMCMC)/include  $(MPILIBS)
endif

clean:
	-rm *.o
