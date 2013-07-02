CCFLAGS=-O3 -std=c99 -Wall -pedantic
LDFLAGS=-lm

all: calc_sasa

calc_sasa: calc_sasa.c src/*.c src/*.h
	gcc calc_sasa.c -o calc_sasa src/*.c $(CCFLAGS) $(LDFLAGS)

debug: test.c calc_sasa.c src/*.c src/*.h
	gcc test.c -o test src/*.c $(CCFLAGS) -g -p -DDEBUG $(LDFLAGS)
	gcc calc_sasa.c -o calc_sasa src/*.c $(CCFLAGS) -g -p $(LDFLAGS)

clean:
	if [ -e calc_sasa ] ; then rm calc_sasa; fi;	
	if [ -e test ] ; then rm test; fi;	
