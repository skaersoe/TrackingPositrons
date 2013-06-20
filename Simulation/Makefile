CC=g++
NVCC=nvcc

NVCCFLAGS=-arch=sm_30 -m64
CCFLAGS=-m64 -std=c++11

SRCDIR=src
CCOBJDIR=ccobjects
NVCCOBJDIR=nvccobjects

CCSRC=$(shell find $(SRCDIR) -name '*.cc')
NVCCSRC=$(shell find $(SRCDIR) -name '*.cu')
CCOBJ=$(patsubst src/%.cc,$(CCOBJDIR)/%.o,$(CCSRC))
NVCCOBJ=$(patsubst src/%.cu,$(NVCCOBJDIR)/%.o,$(NVCCSRC))
LINKOBJ=$(NVCCOBJDIR)/dlink.o

ROOT=`root-config --cflags --libs`
INCLUDE_NVCC=-I.
INCLUDE_CC=-I. -I$(CUDA)/include -L$(CUDA)/lib64 -lcuda -lcudart $(ROOT)
INCLUDE_LINK=$(INCLUDE_CC) $(CCOBJ) $(NVCCOBJ) $(LINKOBJ)

all: simulator

simulator: folders $(CCOBJ) $(NVCCOBJ) $(LINKOBJ)

folders:
	@-mkdir ccobjects
	@-mkdir nvccobjects

EnergyLoss: EnergyLoss.cc
	$(CC) $(CCFLAGS) $(INCLUDE_LINK) EnergyLoss.cc -o bin/EnergyLoss

Depth: Depth.cc
	$(CC) $(CCFLAGS) $(INCLUDE_LINK) Depth.cc -o bin/Depth

Benchmark: Benchmark.cc
	$(CC) $(CCFLAGS) $(INCLUDE_LINK) Benchmark.cc -o bin/Benchmark

MeanDepth: MeanDepth.cc
	$(CC) $(CCFLAGS) $(INCLUDE_LINK) MeanDepth.cc -o bin/MeanDepth

Bremsstrahlung: Bremsstrahlung.cc
	$(CC) $(CCFLAGS) $(INCLUDE_LINK) Bremsstrahlung.cc -o bin/Bremsstrahlung

BremsstrahlungBenchmark: BremsstrahlungBenchmark.cc
	$(CC) $(CCFLAGS) $(INCLUDE_LINK) BremsstrahlungBenchmark.cc -o bin/BremsstrahlungBenchmark

photonpropagate: PhotonPropagate.cc
	$(CC) $(CCFLAGS) $(INCLUDE_LINK) PhotonPropagate.cc -o bin/PhotonPropagate

stepsperlaunchbenchmark: StepsPerLaunchBenchmark.cc
	$(CC) $(CCFLAGS) $(INCLUDE_LINK) StepsPerLaunchBenchmark.cc -o bin/StepsPerLaunchBenchmark

stepsperlaunchbenchmarkdiverse: StepsPerLaunchBenchmarkDiverse.cc
	$(CC) $(CCFLAGS) $(INCLUDE_LINK) StepsPerLaunchBenchmarkDiverse.cc -o bin/StepsPerLaunchBenchmarkDiverse
	
$(CCOBJDIR)/%.o: $(SRCDIR)/%.cc
	$(CC) $(CCFLAGS) $(INCLUDE_CC) -c $< -o $@

$(NVCCOBJDIR)/%.o: $(SRCDIR)/%.cu
	$(NVCC) $(NVCCFLAGS) $(INCLUDE_NVCC) -dc $< -o $@

$(LINKOBJ): $(NVCCOBJ)
	$(NVCC) $(NVCCFLAGS) $(NVCCOBJ) -dlink -o $(LINKOBJ)

atlasstyle:
	$(CC) $(CCFLAGS) $(INCLUDE) -c ../atlasstyle/*.C

landau:
	$(NVCC) -arch=sm_30 $(INCLUDE_NVCC) $(ROOT) $(CCOBJDIR)/Landau.o Landau.cu -o bin/Landau

boxtest:
	$(CC) $(CCFLAGS) $(INCLUDE_LINK) BoxTest.cc -o bin/BoxTest

ProbabilityTest:
	$(CC) $(CCFLAGS) $(INCLUDE_LINK) ProbabilityTest.cc -o bin/ProbabilityTest

test:
	$(CC) $(CCFLAGS) $(INCLUDE_LINK) Test.cc -o bin/Test

EnergyDeposit: EnergyDeposit.cc
	$(CC) $(CCFLAGS) $(INCLUDE_LINK) EnergyDeposit.cc -o bin/EnergyDeposit

clean: cleancc cleancuda
	@-rm bin/*

cleancc:
	@-rm -rf ccobjects

cleancuda:
	@-rm -rf nvccobjects

.PHONY: clean