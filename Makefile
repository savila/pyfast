CC              =       gcc
OPT             =       -O2 -fopenmp -std=c99 -fPIC -shared
MLIB            =       -lm
DEFS            =

place_halos: place_halos.c
		 $(CC) $(OPT) $(DEFS) -fPIC -shared -o place_halos.so place_halos.c $(MLIB)
clean:
	rm -f *.o *~place_halos
