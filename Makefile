all: mcsp

mcsp: mcsp.c graph.c graph.h
	gcc -g -O3 -march=native -mcmodel=medium -Wall -std=c11 -o mcsp graph.c mcsp.c

clean:
	rm -rf mcsp
