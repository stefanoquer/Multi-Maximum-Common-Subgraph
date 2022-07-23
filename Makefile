CXX := g++
CXXFLAGS := -Ofast -march=native
all: mcsp

mcsp: mcsp-mt.cpp graph.cpp graph.h
	$(CXX) $(CXXFLAGS) -Wall -std=c++20 -o mcsp graph.cpp mcsp-mt.cpp -pthread -g -ggdb3

clean:
	rm -f mcsp
