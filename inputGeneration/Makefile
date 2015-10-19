all: 
	g++ -ggdb `pkg-config --cflags opencv` -o getInput inputGeneration.cc -lpthread -lfftw3 -lm -lfftw3f -L -lconfig++ `pkg-config  --libs opencv`
clean:
	rm  getInput

