CPP=g++

run: sphmain
	./sphmain

sphmain: main.cpp
	$(CPP) main.cpp -o sphmain

clean: 
	-rm sphmain
