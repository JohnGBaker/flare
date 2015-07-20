CFLAGS= -O2 

_splinecspagh.so: spline.c spline.h spline.i
	swig -python spline.i
	$(CC) -c `python-config --cflags` spline.c spline_wrap.c
	mv spline.o splinepy.o
	$(CC) -bundle `python-config --ldflags` splinepy.o spline_wrap.o -o _splinecspagh.so


Faddeeva.o: Faddeeva.h Faddeeva.c Faddeeva.cc

wip.o: wip.c wip.h Faddeeva.o

spline.o: spline.c spline.h	

wiptest: wiptest.c wip.o spline.o Faddeeva.o

clean:
	rm -f *.o *.so
