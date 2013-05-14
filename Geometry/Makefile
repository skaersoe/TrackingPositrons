CC=g++
NVCC=nvcc
NVCCFLAGS=-arch=sm_30 -m64
CCFLAGS=-m64 -std=c++11
SRC=src
INCLUDE_CUDA=-I../Geometry -I../Simulation -I$(SRC)
INCLUDE=-I$(CUDA)/include -L$(CUDA)/lib64 -lcuda -lcudart -I../Geometry -I../Simulation -I$(SRC) `root-config --libs --cflags`

all: geometry

test: geometry simulation Main.cc
	$(CC) $(CCFLAGS) $(INCLUDE) *.o ../Simulation/*.o Main.cc -o bin/Geometry

geometry: ccobjects nvccobjects

ccobjects:
	$(CC) $(CCFLAGS) $(INCLUDE) -c $(SRC)/*.cc

nvccobjects:
	$(NVCC) $(NVCCFLAGS) $(INCLUDE_CUDA) -c $(SRC)/*.cu

simulation:
	make -C ../Simulation

clean:
	rm ./*.o