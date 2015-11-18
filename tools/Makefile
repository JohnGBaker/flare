GSLINC = $(GSLROOT)/include
CFLAGS := -O2 -std=c99 -I$(GSLINC) -I../tools -I../integration -I../EOBNRv2HMROM -I../LISAsim -I../LLVsim -I../LLVinference

OBJ = struct.o splinecoeffs.o fresnel.o likelihood.o timeconversion.o


all: $(OBJ)

struct.o: struct.c constants.h struct.h
	$(CC) -c $(CFLAGS) struct.c

splinecoeffs.o: splinecoeffs.c constants.h struct.h splinecoeffs.h
	$(CC) -c $(CFLAGS) splinecoeffs.c

fresnel.o: fresnel.c constants.h struct.h fresnel.h
	$(CC) -c $(CFLAGS) fresnel.c

likelihood.o: likelihood.c constants.h struct.h splinecoeffs.h fresnel.h likelihood.h ../integration/wip.h
	$(CC) -c $(CFLAGS) likelihood.c

timeconversion.o: timeconversion.c constants.h
	$(CC) -c $(CFLAGS) timeconversion.c

clean:
	-rm *.o
