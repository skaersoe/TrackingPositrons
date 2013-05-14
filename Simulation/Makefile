CC=g++
NVCC=nvcc
NVCCFLAGS=-arch=sm_30 -m64
CCFLAGS=-m64 -std=c++11
SRC=src
INCLUDE_CUDA=-I../Geometry -I../Simulation -I$(SRC)
INCLUDE=-I$(CUDA)/include -L$(CUDA)/lib64 -lcuda -lcudart -I../Geometry -I../Simulation -I$(SRC) `root-config --libs --cflags`

all: simulator

test: simulator Main.cc dlink
	$(CC) $(CCFLAGS) $(INCLUDE) ./*.o ../Geometry/*.o Main.cc -o bin/Simulator

energyloss: EnergyLoss.cc
	$(CC) $(CCFLAGS) $(INCLUDE) ./*.o ../Geometry/*.o EnergyLoss.cc -o bin/EnergyLoss

depth: Depth.cc
	$(CC) $(CCFLAGS) $(INCLUDE) ./*.o ../Geometry/*.o Depth.cc -o bin/Depth

simulator: ccobjects nvccobjects geometry

ccobjects:
	$(CC) $(CCFLAGS) $(INCLUDE) -c $(SRC)/*.cc

nvccobjects:
	$(NVCC) $(NVCCFLAGS) $(INCLUDE_CUDA) -dc $(SRC)/*.cu

dlink: nvccobjects
	$(NVCC) $(NVCCFLAGS) -dlink PropagateGPU.o TrackGPU.o

geometry:
	make -C ../Geometry

clean:
	@rm ./*.o