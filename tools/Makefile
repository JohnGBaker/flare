CFLAGS += -I../tools -I../integration -I../EOBNRv2HMROM -I../LISAsim -I../LLVsim -I../LLVinference

OBJ = struct.o splinecoeffs.o fresnel.o likelihood.o timeconversion.o fft.o waveform.o


all: $(OBJ)

struct.o: struct.c constants.h struct.h
	$(CC) -c $(CFLAGS) struct.c

splinecoeffs.o: splinecoeffs.c constants.h struct.h splinecoeffs.h
	$(CC) -c $(CFLAGS) splinecoeffs.c

fresnel.o: fresnel.c constants.h struct.h fresnel.h
	$(CC) -c $(CFLAGS) fresnel.c

likelihood.o: likelihood.c constants.h struct.h waveform.h splinecoeffs.h fresnel.h likelihood.h ../integration/wip.h
	$(CC) -c $(CFLAGS) likelihood.c

timeconversion.o: timeconversion.c constants.h
	$(CC) -c $(CFLAGS) timeconversion.c

fft.o: fft.c fft.h struct.h constants.h
	$(CC) -c $(CFLAGS) fft.c

waveform.o: waveform.c waveform.h struct.h constants.h ../EOBNRv2HMROM/EOBNRv2HMROM.h
	$(CC) -c $(CFLAGS) waveform.c

clean:
	-rm *.o
