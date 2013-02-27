CC=g++
NVCC=nvcc
NVCCFLAGS=-arch=sm_30 --machine=64
CFLAGS=-pedantic
SRC=src
INCLUDE=-I$(SRC) -L$(HEPMC)/lib -lHepMC

all: simulator

simulator: $(SRC)/simulator.cpp GenerateParticles
	$(CC) $(CFLAGS) $(INCLUDE) $(SRC)/simulator.cpp -o bin/simulator
	rm src/*.o

GenerateParticles: $(SRC)/GenerateParticles.cpp
	$(NVCC) $(NVCCFLAGS) --compile $(SRC)/GenerateParticles.cpp

clean:
	rm src/*.o
	rm bin/*