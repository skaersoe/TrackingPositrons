CXX=g++
NVCC=nvcc
NVCCFLAGS=-ccbin /usr/bin -arch=sm_30 -m64
CXXFLAGS=-pedantic -m64 -std=c++11
SRC=src
HEADERS=Simulation
INCLUDE=-I$(SRC) -I$(CUDA)/include -I$(HEADERS) -L$(HEPMC)/lib -L$(CUDA)/lib -lhepmc -lcuda -lcudart

all: Simulator

Simulator: $(SRC)/Simulator.cpp GenerateParticles Step
	$(CXX) $(CXXFLAGS) $(INCLUDE) *.o $(SRC)/Simulator.cpp -o bin/simulator
	rm ./*.o

GenerateParticles: $(SRC)/GenerateParticles.cu
	$(NVCC) $(NVCCFLAGS) $(INCLUDE) -c $(SRC)/GenerateParticles.cu

Step: $(SRC)/Step.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $(SRC)/Step.cpp

clean:
	rm src/*.o
	rm bin/*