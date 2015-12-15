CFLAGS += -I$(GSLINC) -I../tools -I../integration -I../EOBNRv2HMROM -I../LISAsim -I../LLVsim -I../LLVinference

OBJ = LISAFDresponse.o LISAgeometry.o LISAnoise.o #LISAexampleSNR.o


all: $(OBJ) #LISAexampleSNR

LISAgeometry.o: LISAgeometry.c LISAgeometry.h ../tools/constants.h 
	$(CC) -c $(CFLAGS) LISAgeometry.c

LISAFDresponse.o: LISAFDresponse.c  LISAFDresponse.h LISAgeometry.h ../tools/constants.h ../tools/struct.h 
	$(CC) -c $(CFLAGS) LISAFDresponse.c

LISANoise.o: LISANoise.c LISANoise.h ../tools/constants.h ../tools/struct.h 
	$(CC) -c $(CFLAGS) LISANoise.c

clean:
	-rm *.o
