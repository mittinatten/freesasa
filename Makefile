
all:
	gcc calc_sasa.c -o calc_sasa src/*.c -O3 -std=c99 -Wall -pedantic -lm
test:
	gcc test.c -o test src/*.c -g -p -DDEBUG -std=c99 -lm
	gcc calc_sasa.c -o calc_sasa src/*.c -g -p -std=c99 -lm
