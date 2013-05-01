CC=g++
NVCC=nvcc
NVCCFLAGS=-arch=sm_30 -m64
CCFLAGS=-m64 -std=c++11
SRC=src
INCLUDE_CUDA=-I$(CUDA)/include -L$(CUDA)/lib64 -L$(CUDA)/lib -lcuda -lcudart
INCLUDE=-I../Geometry -I../Simulation -I$(SRC) $(INCLUDE_CUDA)

all: simulator

test: simulator
	$(CC) $(CCFLAGS) $(INCLUDE) ./*.o ../Geometry/*.o ./Main.cc -o bin/simulator

simulator: ccobjects nvccobjects geometry

ccobjects:
	$(CC) $(CCFLAGS) $(INCLUDE) -c $(SRC)/*.cc

nvccobjects:
	$(NVCC) $(NVCCFLAGS) $(INCLUDE) -c $(SRC)/*.cu

geometry:
	make -C ../Geometry

clean:
	rm ./*.o