all: mcsp

mcsp: mcsp.cc graph.cc graph.hh
	g++ -Ofast -pthread -march=native -std=c++20 -o mcsp graph.cc mcsp.cc

clean:
	rm -rf mcsp
