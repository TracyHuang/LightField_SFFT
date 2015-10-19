all:
	g++ -ggdb `pkg-config --cflags opencv` -o collectResult resultCollection.cc -lfftw3 -lm -lfftw3f -lconfig++ `pkg-config --libs opencv`
clean:
	rm collectResult
