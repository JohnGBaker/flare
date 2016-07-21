CFLAGS += -I../tools -I../integration -I../EOBNRv2HMROM -I../LISAsim -I../LLVsim -I../LLVinference

OBJ = LLVFDresponse.o LLVnoise.o GenerateLLVFD.o


all: GenerateLLVFD $(OBJ)

LLVFDresponse.o: LLVFDresponse.c  LLVFDresponse.h LLVgeometry.h ../tools/constants.h ../tools/struct.h ../tools/waveform.h ../tools/timeconversion.h
	$(CC) -c $(CFLAGS) LLVFDresponse.c

LLVnoise.o: LLVnoise.c LLVnoise.h ../tools/constants.h
	$(CC) -c $(CFLAGS) LLVnoise.c

GenerateLLVFD.o: LLVgeometry.h LLVFDresponse.h ../tools/constants.h ../tools/struct.h ../tools/timeconversion.h ../EOBNRv2HMROM/EOBNRv2HMROM.h ../EOBNRv2HMROM/EOBNRv2HMROMstruct.h ../tools/waveform.h ../tools/fft.h
	$(CC) -c $(CFLAGS) GenerateLLVFD.c

GenerateLLVFD: GenerateLLVFD.o LLVgeometry.h LLVFDresponse.h LLVFDresponse.o ../tools/constants.h ../tools/struct.h ../tools/timeconversion.h ../EOBNRv2HMROM/EOBNRv2HMROM.h ../EOBNRv2HMROM/EOBNRv2HMROMstruct.h ../tools/waveform.h ../tools/fft.h ../tools/struct.o ../tools/timeconversion.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../tools/waveform.o ../tools/fft.o
	$(LD) $(LDFLAGS) -o GenerateLLVFD GenerateLLVFD.o LLVFDresponse.o ../tools/struct.o ../tools/timeconversion.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../tools/waveform.o ../tools/fft.o -lgsl -lgslcblas -lm -lfftw3

clean:
	-rm *.o
