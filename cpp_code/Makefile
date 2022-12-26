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

static_ap: static_ap.o
	$(CC) $(CFLAGS) -o static_ap static_ap.o


# The main.o target can be written more simply
 
neuronmodel.o: NeuronModel.cpp 
	$(CC) $(CFLAGS) NeuronModel.cpp -c -o neuronmodel.o

static_ap.o: Static_AP.cpp 
	$(CC) $(CFLAGS) Static_AP.cpp -c -o static_ap.o

.PHONY : clean
clean :
	rm -rf neuronmodel static_ap
	rm -rf *.ll *.exe *.s *.o testall.log
	
.PHONY : clean_img
clean_img:
	rm -rf ../graphs/*.png