CFLAGS += -I../tools -I../integration -I../EOBNRv2HMROM -I../LISAsim -I../LLVsim -I../LLVinference

OBJ = EOBNRv2HMROM.o EOBNRv2HMROMstruct.o GenerateWaveform.o GenerateWaveform


all: $(OBJ)

EOBNRv2HMROMstruct.o: EOBNRv2HMROMstruct.c ../tools/constants.h ../tools/struct.h EOBNRv2HMROMstruct.h
	$(CC) -c $(CFLAGS) EOBNRv2HMROMstruct.c

EOBNRv2HMROM.o: EOBNRv2HMROM.c ../tools/constants.h ../tools/struct.h EOBNRv2HMROM.h EOBNRv2HMROMstruct.h
	$(CC) -c $(CFLAGS) EOBNRv2HMROM.c

GenerateWaveform.o: EOBNRv2HMROM.h EOBNRv2HMROMstruct.h GenerateWaveform.h GenerateWaveform.c ../tools/constants.h ../tools/struct.h ../tools/waveform.h ../tools/fft.h
	$(CC) -c $(CFLAGS) GenerateWaveform.c

GenerateWaveform: GenerateWaveform.o EOBNRv2HMROM.h EOBNRv2HMROMstruct.h ../tools/constants.h ../tools/struct.h ../tools/waveform.h ../tools/fft.h EOBNRv2HMROM.o EOBNRv2HMROMstruct.o ../tools/struct.o ../tools/waveform.o ../tools/fft.o
	$(LD) $(LDFLAGS) -o GenerateWaveform GenerateWaveform.o EOBNRv2HMROM.o EOBNRv2HMROMstruct.o ../tools/struct.o ../tools/waveform.o ../tools/fft.o -lgsl -lgslcblas -lm -lfftw3

clean:
	-rm *.o
