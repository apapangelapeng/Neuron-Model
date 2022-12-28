# Neuron-Model
To compile the executable run `make test` in the `cpp_code` directory \
To compile the Neuron Model Exectuable run `make neuronmodel` in the `cpp_code` directory  \
To run the python code run `python graph_test.py` in the `python_code` directory. The resulting graph will be in the `graphs` directory
To remove all the png graphs generated run `make clean_img` in the `cpp_code` directory
To regenerate all the png graphs first run `bash to_python_script` in the `python_code` directory, then run `python graph_test.py` 
# Graph Generation regarding Current
To generate all the graphs with the corresponding current run `bash generate_graph.sh $injected_current` where `$injected_current` is a number in the `cpp_code` directory

Set up VM on google cloud
