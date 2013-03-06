CC=GCC
CXX=g++
NVCC=nvcc
NVCCFLAGS=-ccbin /usr/bin -arch=sm_30 -m64
CFLAGS=-m64 -std=c99
CXXFLAGS=-m64 -std=c++03
SRC=src
INCLUDE_CUDA=-I$(CUDA)/include -L$(CUDA)/lib64 -L$(CUDA)/lib -lcuda -lcudart
INCLUDE_ROOT=-I$(ROOTSYS)/lib
INCLUDE=-I./ -I$(SRC) $(INCLUDE_CUDA)

all: Simulator

Simulator: $(SRC)/Simulator.cpp GenerateParticles.o Propagate.o CudaHelper.o
	$(CXX) $(CXXFLAGS) $(INCLUDE) *.o $(SRC)/Simulator.cpp -o bin/simulator
	rm ./*.o

GenerateParticles.o: $(SRC)/GenerateParticles.cu
	$(NVCC) $(NVCCFLAGS) $(INCLUDE) -c $(SRC)/GenerateParticles.cu

Propagate.o: $(SRC)/Propagate.cu
	$(NVCC) $(NVCCFLAGS) $(INCLUDE) -c $(SRC)/Propagate.cu

CudaHelper.o: $(SRC)/CudaHelper.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $(SRC)/CudaHelper.cpp

clean:
	rm src/*.o
	rm bin/*