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


propogate_ap: propogate_ap.o
	$(CC) $(CFLAGS) -o propogate_ap propogate_ap.o


reduced_ap: reduced_ap.o
	$(CC) $(CFLAGS) -o reduced_ap reduced_ap.o


reduced_ap_stochastic: reduced_ap_stochastic.o
	$(CC) $(CFLAGS) -o reduced_ap_stochastic reduced_ap_stochastic.o


# The main.o target can be written more simply
 



propogate_ap.o: Propogate_AP.cpp 
	$(CC) $(CFLAGS) propogate_AP.cpp -c -o propogate_ap.o

reduced_ap.o: Reduced_AP.cpp
	$(CC) $(CFLAGS) Reduced_AP.cpp -c -o reduced_ap.o



reduced_ap_stochastic.o: Reduced_AP_Stochastic.cpp
	$(CC) $(CFLAGS) Reduced_AP_Stochastic.cpp -c -o reduced_ap_stochastic.o




.PHONY : clean
clean :
	rm -rf neuronmodel static_ap reduced_ap reduced_ap_stochastic Reduced_AP_Stochastic
	rm -rf *.ll *.exe *.s *.o testall.log
	
.PHONY : clean_img
clean_img:
	rm -rf ../graphs/*.png

clean_csv:
	rm -rf -rf ../data_files/*.csv