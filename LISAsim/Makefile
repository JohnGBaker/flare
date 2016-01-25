CFLAGS += -I$(GSLINC) -I../tools -I../integration -I../EOBNRv2HMROM -I../LISAsim -I../LLVsim -I../LLVinference

OBJ = LISAFDresponse.o LISAgeometry.o LISAnoise.o GenerateTDITD.o GenerateTDITD GenerateTDIFD.o GenerateTDIFD #LISAexampleSNR.o


all: $(OBJ) #LISAexampleSNR

LISAgeometry.o: LISAgeometry.c LISAgeometry.h ../tools/constants.h 
	$(CC) -c $(CFLAGS) LISAgeometry.c

LISAFDresponse.o: LISAFDresponse.c  LISAFDresponse.h LISAgeometry.h ../tools/constants.h ../tools/struct.h ../tools/waveform.h 
	$(CC) -c $(CFLAGS) LISAFDresponse.c

LISANoise.o: LISANoise.c LISANoise.h ../tools/constants.h ../tools/struct.h 
	$(CC) -c $(CFLAGS) LISANoise.c

GenerateTDITD.o: LISAgeometry.h ../tools/constants.h ../tools/struct.h ../EOBNRv2HMROM/EOBNRv2HMROM.h ../EOBNRv2HMROM/EOBNRv2HMROMstruct.h ../tools/waveform.h
	$(CC) -c $(CFLAGS) GenerateTDITD.c

GenerateTDIFD.o: LISAgeometry.h ../tools/constants.h ../tools/struct.h ../EOBNRv2HMROM/EOBNRv2HMROM.h ../EOBNRv2HMROM/EOBNRv2HMROMstruct.h ../tools/waveform.h ../tools/fft.h
	$(CC) -c $(CFLAGS) GenerateTDIFD.c

GenerateTDITD: GenerateTDITD.o LISAgeometry.h LISAgeometry.o ../tools/constants.h ../tools/struct.h ../EOBNRv2HMROM/EOBNRv2HMROM.h ../EOBNRv2HMROM/EOBNRv2HMROMstruct.h ../tools/waveform.h ../tools/struct.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../tools/waveform.o
	$(LD) $(LDFLAGS) -o GenerateTDITD GenerateTDITD.o LISAgeometry.o ../tools/struct.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../tools/waveform.o -lgsl -lgslcblas -lm

GenerateTDIFD: GenerateTDIFD.o LISAgeometry.h LISAgeometry.o LISAFDresponse.h LISAFDresponse.o ../tools/constants.h ../tools/struct.h ../EOBNRv2HMROM/EOBNRv2HMROM.h ../EOBNRv2HMROM/EOBNRv2HMROMstruct.h ../tools/waveform.h ../tools/fft.h ../tools/struct.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../tools/waveform.o ../tools/fft.o
	$(LD) $(LDFLAGS) -o GenerateTDIFD GenerateTDIFD.o LISAgeometry.o LISAFDresponse.o ../tools/struct.o ../EOBNRv2HMROM/EOBNRv2HMROM.o ../EOBNRv2HMROM/EOBNRv2HMROMstruct.o ../tools/waveform.o ../tools/fft.o -lgsl -lgslcblas -lm -lfftw3

clean:
	-rm *.o
