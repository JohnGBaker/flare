# Makefile for building a shared library that can be used from Python via ctypes

all: libflare.so.1.0.0

COMPILE=gcc -c -O3 -fPIC -o $@ -I../tools -I../EOBNRv2HMROM -I../integration \
		-I../LISAsim

OBJECTS=LISAutils.o LISAgeometry.o LISAFDresponse.o LISAnoise.o struct.o \
		waveform.o fresnel.o EOBNRv2HMROM.o EOBNRv2HMROMstruct.o \
		splinecoeffs.o likelihood.o wip.o Faddeeva.o spline.o

LISAutils.o: ../LISAinference/LISAutils.c
	$(COMPILE) ../LISAinference/LISAutils.c

LISAgeometry.o: ../LISAsim/LISAgeometry.c
	$(COMPILE) ../LISAsim/LISAgeometry.c

LISAFDresponse.o: ../LISAsim/LISAFDresponse.c
	$(COMPILE) ../LISAsim/LISAFDresponse.c

LISAnoise.o: ../LISAsim/LISAnoise.c
	$(COMPILE) ../LISAsim/LISAnoise.c

EOBNRv2HMROM.o: ../EOBNRv2HMROM/EOBNRv2HMROM.c
	$(COMPILE) ../EOBNRv2HMROM/EOBNRv2HMROM.c

EOBNRv2HMROMstruct.o: ../EOBNRv2HMROM/EOBNRv2HMROMstruct.c
	$(COMPILE) ../EOBNRv2HMROM/EOBNRv2HMROMstruct.c

struct.o: ../tools/struct.c
	$(COMPILE) ../tools/struct.c

waveform.o: ../tools/waveform.c
	$(COMPILE) ../tools/waveform.c

splinecoeffs.o: ../tools/splinecoeffs.c
	$(COMPILE) ../tools/splinecoeffs.c

likelihood.o: ../tools/likelihood.c
	$(COMPILE) ../tools/likelihood.c

fresnel.o: ../tools/fresnel.c
	$(COMPILE) ../tools/fresnel.c

wip.o: ../integration/wip.c
	$(COMPILE) ../integration/wip.c

Faddeeva.o: ../integration/Faddeeva.c
	$(COMPILE) ../integration/Faddeeva.c

spline.o: ../integration/spline.c
	$(COMPILE) ../integration/spline.c

libflare.so.1.0.0: $(OBJECTS)
	gcc -shared -fPIC \
		-Wl,-soname,libflare.so.1,--unresolved-symbols=report-all \
        -lgsl -lgslcblas -lm -lc -o $@ $(OBJECTS)

clean:
	rm -f *.o libflare.so.1.0.0