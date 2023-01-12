# Code Purpose 
Hello! The overall purpose of this code is to test differnet neuron models and to use them to learn about neuron behavior, particularly with regard to excitability.\\
In general, the bulk of the algorithmic load is done in c++ and then the output is graphed using python. \\
A brief rundown of the `main` branch's current files: \\
`Static_AP.cpp` uses the standard Hodgkin-Huxely model to generate the time constants (tau) and dynamical variables (m, h, n) of sodium and potassium ion channels. This is then used to generate a trace of a neuron's voltage. It can also take in some input current as desired. This is output into a `.csv` file, which is read by the Jupyter Notebook file `Static_AP.ipynb`.\\
`Reduced_AP.cpp` uses the standard Fitzhugh-Nagumo Reduction model to generate the change in voltage and blocking/inactivation of voltage (n) in response to initial parameters. This is output as a `.csv` and can be used to graph a static action potential using the Jupyter Notebook file `Reduced_AP.ipynb`. An equally important part of this file is generating a phase plane/nullcline, which is outputted and graphed in `Reduced_AP.ipynb` as well. There is another component, entitled diffusion, which has been abandoned in favor of a different file, commente on next. \\
`Propogate_AP.cpp` is meant to track a moving action potential. It uses the Fitzhugh-Nagumo Reduction model with the added spatial component. It is currently not functional, though. \\
`Static_AP_HCN.cpp` is used to test the affect of adding HCN channels into the Hodgkin-Huxely model. There are a variety of `.csv` outputs, which are graphed in a variety of `.ipynb` files. Likely the most important being `Static_Vcompare.ipynb`, which plots the voltage including and excluding HCN channels on the same graph. The other outputs are important, though, and are used to track things like current in an independent manner. \\
Okay, that's all for now. Hopefully the file `Angela.cpp` doesn't `Segmentation Fault` anymore because I finally wrote this section of the `README.md`.

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
