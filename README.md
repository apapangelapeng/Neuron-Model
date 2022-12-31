# Neuron-Model
To compile the executable run `make test` in the `cpp_code` directory \
To compile the Neuron Model Exectuable run `make neuronmodel` in the `cpp_code` directory  \
To run the python code run `python graph_test.py` in the `python_code` directory. The resulting graph will be in the `graphs` directory
To remove all the png graphs generated run `make clean_img` in the `cpp_code` directory
To regenerate all the png graphs first run `bash to_python_script` in the `python_code` directory, then run `python graph_test.py` 
# Graph Generation regarding Current
To generate all the graphs with the corresponding current run `bash generate_graph.sh $injected_current` where `$injected_current` is a number in the `cpp_code` directory

To compile and run the executable run `bash run <cpp_file_name_without_.cpp> in the `cpp_code` directory \

 To remove all the png graphs generated run `make clean_img` in the `cpp_code` directory

# Graphing
The `graph()` accepts the dataframe, begin_index, end_index(but does not include the end index). If no begin_index is given, then it will automatically be 0. If no end index is given the function will graph till the end of the columns.
