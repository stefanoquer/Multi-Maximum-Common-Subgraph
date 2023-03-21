NVCC = /usr/local/cuda-11.7/bin/nvcc
NVCCFLAGS = -arch=compute_86

all:
	$(NVCC) $(NVCCFLAGS) -std=c++17 -O2 mcsp-mt.cu graph.cu -o ./build/McSplit_Multigraph

clean:
	rm -f mcsp
