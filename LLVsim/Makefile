GSLINC = $(GSLROOT)/include
CFLAGS := -O2 -std=c99 -I$(GSLINC) -I../tools -I../integration -I../EOBNRv2HMROM -I../LISAsim -I../LLVsim -I../LLVinference

OBJ = LLVFDresponse.o LLVnoise.o


all: $(OBJ)

LLVFDresponse.o: LLVFDresponse.c  LLVFDresponse.h LLVgeometry.h ../tools/constants.h ../tools/struct.h ../tools/timeconversion.h
	$(CC) -c $(CFLAGS) LLVFDresponse.c

LLVnoise.o: LLVnoise.c LLVnoise.h ../tools/constants.h 
	$(CC) -c $(CFLAGS) LLVnoise.c

clean:
	-rm *.o
