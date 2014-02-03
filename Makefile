CC=gcc
CFLAGS=-O2 -std=c99 -Wall -pedantic -DPTHREADS
LDFLAGS=-lm -lpthread
srcdir=src/

all: calc_sasa 

.PHONY: doc

doc:
	$(MAKE) -C doc/

calc_sasa: calc_sasa.c src/*.c src/*.h
	$(CC) calc_sasa.c -o calc_sasa src/*.c $(CFLAGS) $(LDFLAGS)

example: example.c src/*.c src/*.h
	$(CC) example.c -o example src/*.c $(CFLAGS) $(LDFLAGS)

debug: calc_sasa.c src/*.c src/*.h
	$(CC) calc_sasa.c -o calc_sasa src/*.c $(CFLAGS) \
	-g -p -DDEBUG $(LDFLAGS)

prof: calc_sasa.c src/*.c src/*.h
	$(CC) calc_sasa.c -o calc_sasa src/*.c $(CFLAGS) \
	-g -p -fprofile-arcs -ftest-coverage $(LDFLAGS)

test: test.c src/*.c src/*.h
	$(CC) test.c -o test tests/*.c src/*.c -Isrc/ $(CFLAGS) -g $(LDFLAGS)
	./test

clean:	
	rm -f calc_sasa example test
	rm -f *.gcda *.gcno gmon.out *.gcov *~
