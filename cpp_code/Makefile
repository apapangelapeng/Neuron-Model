.PHONY : all
all : neuronmodel

neuronmodel:
	g++ -o neuronmodel NeuronModel.cpp

test:
	g++ -o test Testing.cpp

# "make clean" removes all generated files

.PHONY : clean
clean :
	rm -rf neuronmodel test
	rm -rf *.ll *.exe *.s *.o testall.log
