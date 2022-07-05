all: mcsp

mcsp: mcsp.cpp graph.cpp graph.h SecureQueue.h
	g++ -O3 -pthread -march=native -std=c++17 -o mcsp graph.cpp mcsp.cpp

clean:
	rm -rf mcsp
