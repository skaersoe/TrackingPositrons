CC=GCC
CXX=g++
NVCC=nvcc
NVCCFLAGS=-arch=sm_30 -m64
CFLAGS=-m64 -std=c99 `root-config --cflags`
CXXFLAGS=-m64 -std=c++11 `root-config --cflags`
SRC=src
INCLUDE_CUDA=-I$(CUDA)/include -L$(CUDA)/lib64 -L$(CUDA)/lib -lcuda -lcudart
INCLUDE=-I./ -I../Geometry -I$(SRC) $(INCLUDE_CUDA)

all: main

main: Simulator.o

standalone: Simulator.o
	$(CXX) $(CXXFLAGS) $(INCLUDE) *.o $(SRC)/Main.cpp -o bin/simulator

Simulator.o: $(SRC)/Simulator.cpp GenerateTracks.o Propagate.o Geometry.o
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $(SRC)/Simulator.cpp

GenerateTracks.o: $(SRC)/GenerateTracks.cu
	$(NVCC) $(NVCCFLAGS) $(INCLUDE) -c $(SRC)/GenerateTracks.cu

Propagate.o: $(SRC)/Propagate.cu CudaHelper.o
	$(NVCC) $(NVCCFLAGS) $(INCLUDE) -c $(SRC)/Propagate.cu

CudaHelper.o: $(SRC)/CudaHelper.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $(SRC)/CudaHelper.cpp

Geometry.o:
	make -C ../Geometry

clean:
	rm ./*.o