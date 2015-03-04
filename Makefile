CC    = g++
FLAGS = -std=c++11 -O3 -Wall -mavx

all: testit fasttestit

testit: main.cpp dualquat.hpp
	$(CC) $(FLAGS) main.cpp -o testit 
	time ./testit

fasttestit: main.cpp dualquat.hpp
	$(CC) $(FLAGS) main.cpp -o fasttestit -march=native -Ofast
	time ./fasttestit

clean:
	rm -rf testit fasttestit
