NOTE TO SELF: 
HAVING BUFFER CONCENTRATION CHANGE IS ACTUALLY SUPER DUMB, SINCE THE ENTIRE NEURON HOLDS A WEALTH OF BUFFER. BE LESS DUMB. 


# Model_Overview.tex

Note: To this point, this code had been an entirely "for fun" exploration of some topics in electrophysiology. The ReadMe is mostly out of date and only describes some of the original elements in the project. A more in-depth discussion can be found in the .tex file titled Model_Overview. In this file you will find some of the mathematical basis of the code, as well as general concepts in neurbiology that are implemented.  


# Compilation of Executables
To compile all the executable run `make all`.
To just compile one executable run `make <exectuable_name>` here are the available names
- static_ap 
- nullcline 
- propogate_ap 
- reduced_ap 
- static_wt_ap 
- static_ap_hcn 
- static_vcompare 
- static_hcn_ap

## Clean Up
To remove all the executables and `.o` files, run `make clean`
To remove all the `.csv` files run `make clean_csv`

# Graph Generation

The graphing code is stored in the jupyter notebook under the `python_code` directory. Each executable has its own graph generation code, simply press run all the inspect the output.

# Data Testing
The constant of each type of cpp file is stored in the corresponding header files `*.h` . For the reduce code the constant is stored in `ReduceConst.h`.
All the static action potential constants are stored in `StaticConst.h`

# Code Purpose
Hello! The overall purpose of this code is to test differnet neuron models and to use them to learn about neuron behavior, particularly with regard to excitability.\
In general, the bulk of the algorithmic load is done in c++ and then the output is graphed using python. \
A brief rundown of the `main` branch’s current files: \
`Static_AP.cpp` uses the standard Hodgkin-Huxely model to generate the time constants (tau) and dynamical variables (m, h, n) of sodium and potassium ion channels. This is then used to generate a trace of a neuron’s voltage. It can also take in some input current as desired. This is output into a `.csv` file, which is read by the Jupyter Notebook file `Static_AP.ipynb`.\
`Reduced_AP.cpp` uses the standard Fitzhugh-Nagumo Reduction model to generate the change in voltage and blocking/inactivation of voltage (n) in response to initial parameters. This is output as a `.csv` and can be used to graph a static action potential using the Jupyter Notebook file `Reduced_AP.ipynb`. An equally important part of this file is generating a phase plane/nullcline, which is outputted and graphed in `Reduced_AP.ipynb` as well. There is another component, entitled diffusion, which has been abandoned in favor of a different file, commente on next. \
`Propogate_AP.cpp` is meant to track a moving action potential. It uses the Fitzhugh-Nagumo Reduction model with the added spatial component. It is currently not functional, though. \
`Static_AP_HCN.cpp` is used to test the affect of adding HCN channels into the Hodgkin-Huxely model. There are a variety of `.csv` outputs, which are graphed in a variety of `.ipynb` files. Likely the most important being `Static_Vcompare.ipynb`, which plots the voltage including and excluding HCN channels on the same graph. The other outputs are important, though, and are used to track things like current in an independent manner. \
Okay, that’s all for now. Hopefully the file `Angela.cpp` doesn’t `Segmentation Fault` anymore because I finally wrote this section of the `README.md`.
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

# CPP and Header Files

## nullcline.cpp

This file illustrate ....

# Propogate_AP.cpp

This File illustrate

## Reduced_AP

This file illustrate

## Static_Ap .cpp

This file illustrate
Make sure to pass a $I$ value (current) value when call the executable. Like this `./static_ap 0.1`

## Static_AP_HCN.cpp

This file illustrate

## Static_WT_AP.cpp

## Static_Vcompare.cpp

## DynamicFunc.h
This header file stores the Hodgkin-Huxley model

$$C\dot{V} = I-\overbrace{\bar{g}_k n^4(V-E_k)}^{I_k}-\overbrace{\bar{g}_{Na}m^3h(V-E_{Na})}^{I_{Na}}-\overbrace{g_L(V-E_L)}^{I_L}$$


$$\begin{align*}
    \dot{n} &=  \alpha_n(V)(10n)-\beta_n(V)n\\
     \dot{m} &= \alpha_m(V)(1-m)-\beta_m(V)m\\
    \dot{h} &=  \alpha_h(V)(1-h)-\beta_h(V)h\\
\end{align*}$$




s.t

$$\begin{align*}
    \alpha_n(V) &= 0.01 \frac{10-V}{exp(\frac{10-V}{10})-1}\\
    \beta_n(V) &= 0.125 \exp(\frac{-V}{80})\\
    \alpha_m (V) &= 0.1\quad \frac{25-V}{\exp(\frac{25-V}{10})-1}\\
    \beta_m(V) &= 4\exp(\frac{-V}{18})
\end{align*}$$


This header file also stores the for loop that updates the output voltage with regard to the value of injected current $I$.

