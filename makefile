CPP=g++

run: sphmain
	./sphmain

vtest: FORCE #VectorND.h VectorNDTest.cpp
	$(CPP) VectorNDTest.cpp -o vtest
	./vtest


sphmain: main.cpp
	$(CPP) main.cpp -o sphmain

clean: 
	-rm sphmain

FORCE: 
