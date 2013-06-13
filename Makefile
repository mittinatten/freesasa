all:
	gcc calc_sasa.c -o calc_sasa src/*.c -O3 -std=c99 -lm
test:
	gcc test.c -o test src/*.c -g -p -DDEBUG -std=c99 -lm
