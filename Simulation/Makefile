CC=g++
CFLAGS=-pedantic -m64
SRC=src
INCLUDE=-I$(SRC)

all: simulator

simulator: $(SRC)/simulator.cpp particle
	$(CC) $(CFLAGS) $(SRC)/simulator.cpp $(SRC)/Particle.o -o bin/simulator

particle: $(SRC)/Particle.cpp
	$(CC) $(CFLAGS) -c $(SRC)/Particle.cpp

clean:
	rm src/*.o
	rm bin/*