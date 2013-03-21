CC=GCC
CXX=g++
NVCC=nvcc
NVCCFLAGS=-ccbin /usr/bin -arch=sm_30 -m64
CFLAGS=-m64 -std=c99 `root-config --cflags`
CXXFLAGS=-m64 -std=c++11 `root-config --cflags`
SRC=src
INCLUDE_CUDA=$(INCLUDE) -I$(CUDA)/include -L$(CUDA)/lib64 -L$(CUDA)/lib -lcuda -lcudart
INCLUDE=-I./ -I$(SRC)

all: main

main: Simulator.o
	$(CXX) $(CXXFLAGS) $(INCLUDE_CUDA) *.o $(SRC)/Main.cpp -o bin/simulator

Simulator.o: $(SRC)/Simulator.cpp GenerateParticles.o Propagate.o
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $(SRC)/Simulator.cpp

GenerateParticles.o: $(SRC)/GenerateParticles.cu
	$(NVCC) $(NVCCFLAGS) $(INCLUDE_CUDA) -c $(SRC)/GenerateParticles.cu

Propagate.o: $(SRC)/Propagate.cu CudaHelper.o
	$(NVCC) $(NVCCFLAGS) $(INCLUDE_CUDA) -c $(SRC)/Propagate.cu

CudaHelper.o: $(SRC)/CudaHelper.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $(SRC)/CudaHelper.cpp

clean:
	rm ./*.o
	rm bin/*

.PHONY: clean