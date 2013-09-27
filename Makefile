CFLAGS=-O3 -std=c99 -Wall -pedantic -DPTHREADS
LDFLAGS=-lm -lpthread

all: calc_sasa 

calc_sasa: calc_sasa.c src/*.c src/*.h
	gcc calc_sasa.c -o calc_sasa src/*.c $(CFLAGS) $(LDFLAGS)

example: example.c src/*.c src/*.h
	gcc example.c -o example src/*.c $(CFLAGS) $(LDFLAGS)

debug: calc_sasa.c src/*.c src/*.h
	gcc calc_sasa.c -o calc_sasa src/*.c $(CFLAGS) -g -p -DDEBUG $(LDFLAGS)

prof: calc_sasa.c src/*.c src/*.h
	gcc calc_sasa.c -o calc_sasa src/*.c $(CFLAGS) -g -p $(LDFLAGS)

clean:
	@if [ -e calc_sasa ] ; then rm calc_sasa; fi;	
	@if [ -e example ] ; then rm example; fi;	

