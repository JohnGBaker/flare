CFLAGS += -I../tools -I../EOBNRv2HMROM 

OBJ = BBHWaveformFD.o


all: $(OBJ)

%.o : %.c 
	$(CC) -c $(CFLAGS) $<

clean:
	-rm *.o
