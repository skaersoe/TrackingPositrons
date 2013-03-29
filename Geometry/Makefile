CC=g++
NVCC=nvcc
NVCCFLAGS=-arch=sm_30 -m64
CCFLAGS=-m64 -std=c++11
SRC=src
INCLUDE_CUDA=-I$(CUDA)/include -L$(CUDA)/lib64 -L$(CUDA)/lib -lcuda -lcudart
INCLUDE=-I../Geometry -I../Simulation -I$(SRC) $(INCLUDE_CUDA)

geometry: ccobjects nvccobjects

ccobjects:
	$(CC) $(CCFLAGS) $(INCLUDE) -c $(SRC)/*.cc

nvccobjects:
	$(NVCC) $(NVCCFLAGS) $(INCLUDE) -c $(SRC)/*.cu

clean:
	rm ./*.o