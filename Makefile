CC    = g++
FLAGS = -std=c++11 -O3 -Wall

all: testit fasttestit blending

testit: main.cpp dualquat.hpp
	$(CC) $(FLAGS) main.cpp -o testit 
	time ./testit

fasttestit: main.cpp dualquat.hpp
	$(CC) $(FLAGS) main.cpp -o fasttestit -march=native -Ofast
	time ./fasttestit

blending: blending.cpp dualquat.hpp
	$(CC) $(FLAGS) blending.cpp -o blending 
	time ./blending

clean:
	rm -rf testit fasttestit blending
