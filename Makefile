NVCC = /usr/local/cuda-11.6/bin/nvcc
NVCCFLAGS = -arch=compute_86

all:
	$(NVCC) $(NVCCFLAGS) -O3 mcsp-mt.cu graph.cu -o mcsp

clean:
	rm -f mcsp
