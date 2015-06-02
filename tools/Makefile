GSLROOT = /opt/local
GSLINC = $(GSLROOT)/include
CC = gcc
CFLAGS := -g -std=c99 -I$(GSLINC) -I../tools -I../integration -I../EOBNRv2HMROM -I../LISAsim -I../LLVsim

OBJ = struct.o likelihood.o


all: $(OBJ)

struct.o: struct.c constants.h struct.h
	$(CC) -c $(CFLAGS) struct.c

likelihood.o: likelihood.c constants.h struct.h likelihood.h ../integration/wip.h
	$(CC) -c $(CFLAGS) likelihood.c

clean:
	-rm *.o
