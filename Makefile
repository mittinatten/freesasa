CCFLAGS=-O3 -std=c99 -Wall -pedantic
LDFLAGS=-lm

all: calc_sasa

calc_sasa:
	gcc calc_sasa.c -o calc_sasa src/*.c $(CCFLAGS) -DNDEBUG $(LDFLAGS)
debug:
	gcc test.c -o test src/*.c $(CCFLAGS) -g -p -DDEBUG $(LDFLAGS)
	gcc calc_sasa.c -o calc_sasa src/*.c $(CCFLAGS) -g -p $(LDFLAGS)

clean:
	if [ -a calc_sasa ] ; then rm calc_sasa; fi;	
	if [ -a test ] ; then rm test; fi;	
