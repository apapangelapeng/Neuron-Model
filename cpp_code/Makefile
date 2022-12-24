# Makefile for Writing Make Files Example
 
# *****************************************************
# Variables to control Makefile operation
 
CC = g++
CFLAGS = -Wall -g
 
# ****************************************************
# Targets needed to bring the executable up to date
#makes the neuron model dependencies
neuronmodel: neuronmodel.o
	$(CC) $(CFLAGS) -o neuronmodel neuronmodel.o

test: test.o
	$(CC) $(CFLAGS) -o test test.o


# The main.o target can be written more simply
 
neuronmodel.o: NeuronModel.cpp 
	$(CC) $(CFLAGS) NeuronModel.cpp -c -o neuronmodel.o

test.o: Testing.cpp 
	$(CC) $(CFLAGS) Testing.cpp -c -o test.o

.PHONY : clean
clean :
	rm -rf neuronmodel test
	rm -rf *.ll *.exe *.s *.o testall.log
