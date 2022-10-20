CFLAGS = -Ofast -g -Wall -std=c++14 -fopenmp

all: orbits

orbits: orbits.cpp
	g++ $(CFLAGS) -o orbits orbits.cpp
clean:
	rm -f *.o *.txt orbits
