CFLAGS += -I../tools -I../integration -I../IMRPhenomD

OBJ = IMRPhenomD.o IMRPhenomD_internals.o


all: $(OBJ)

IMRPhenomD_internals.o: IMRPhenomD_internals.c ../tools/constants.h ../tools/struct.h IMRPhenomD_internals.h
	$(CC) -c $(CFLAGS) IMRPhenomD_internals.c

IMRPhenomD.o: IMRPhenomD.c ../tools/constants.h ../tools/struct.h IMRPhenomD.h IMRPhenomD_internals.h
	$(CC) -c $(CFLAGS) IMRPhenomD.c

clean:
	-rm *.o
