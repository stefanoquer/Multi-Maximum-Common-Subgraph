all: mcsp

mcsp: mcsp.cc graph.cc graph.hh SecureQueue.h
	g++ -O3 -pthread -march=native -std=c++17 -o mcsp graph.cc mcsp.cc

clean:
	rm -rf mcsp
