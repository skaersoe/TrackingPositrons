CC=GCC
CXX=g++
NVCC=nvcc
NVCCFLAGS=-arch=sm_30 -m64
CFLAGS=-m64 -std=c99 `root-config --cflags`
CXXFLAGS=-m64 -std=c++11 `root-config --cflags`
SRC=src
INCLUDE_CUDA=$(INCLUDE) -I$(CUDA)/include -L$(CUDA)/lib64 -L$(CUDA)/lib -lcuda -lcudart
INCLUDE=-I./ -I$(SRC) -I../Simulation

all: main

main: Geometry.o

Geometry.o: Material.o Volumes
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $(SRC)/Geometry.cpp

Material.o:
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $(SRC)/Material.cpp

Volumes:
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $(SRC)/Sphere.cpp

clean:
	rm ./*.o

.PHONY: clean