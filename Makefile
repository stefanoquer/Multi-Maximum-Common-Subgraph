all: mcsp

mcsp: mcsp.cc graph.cc graph.hh
	g++ -O3 -pthread -march=native -std=c++17 -o mcsp graph.cc mcsp.cc

clean:
	rm -rf mcsp
