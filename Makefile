CXX := g++
CXXFLAGS := -O3 -march=native
all: mcsp

mcsp: mcsp-mt.cpp graph.cpp graph.h
	$(CXX) $(CXXFLAGS) -Wall -std=c++11 -o mcsp graph.cpp mcsp-mt.cpp -pthread -g -ggdb3

clean:
	rm -f v0.2_trimble_par_cpp
	rm -f *.o
