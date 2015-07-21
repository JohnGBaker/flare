GSLINC = $(GSLROOT)/include
CFLAGS := -O2 -std=c99 -I$(GSLINC) -I../tools -I../integration -I../EOBNRv2HMROM -I../LISAsim -I../LLVsim -I../LLVinference

OBJ = EOBNRv2HMROM.o EOBNRv2HMROMstruct.o


all: $(OBJ)

EOBNRv2HMROMstruct.o: EOBNRv2HMROMstruct.c ../tools/constants.h ../tools/struct.h EOBNRv2HMROMstruct.h
	$(CC) -c $(CFLAGS) EOBNRv2HMROMstruct.c

EOBNRv2HMROM.o: EOBNRv2HMROM.c ../tools/constants.h ../tools/struct.h EOBNRv2HMROM.h EOBNRv2HMROMstruct.h
	$(CC) -c $(CFLAGS) EOBNRv2HMROM.c

clean:
	-rm *.o
