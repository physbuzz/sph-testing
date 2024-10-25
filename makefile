CPP=g++

run: sphmain
	./sphmain

vtest: FORCE #VectorND.h VectorNDTest.cpp
	$(CPP) VectorNDTest.cpp -o vtest
	./vtest

fast: main.cpp
	$(CPP) -O3 -funsafe-math-optimizations -DNDEBUG main.cpp -o sphmain

fastrun: fast
	./sphmain

sphmain: main.cpp
	$(CPP) main.cpp -o sphmain

clean: 
	-rm sphmain

FORCE: 
