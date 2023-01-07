# Compilation

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
To remove all the executables and `.o` files, run `make clean`
To remove all the `.csv` files run `make clean_csv`

# Graph Generation

The graphing code is stored in the jupyter notebook under the `python_code` directory. Each executable has its own graph generation code, simply press run all the inspect the output.

# Data Testing
The constant of each type of cpp file is stored in the corresponding header files `*.h` . For the reduce code the constant is stored in `ReduceConst.h`.
All the static action potential constants are stored in `StaticConst.h`

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
```math 
$$\begin{align*}
    C\dot{V} &= I-\overbrace{\bar{g}_k n^4(V-E_k)}^{I_k}-\overbrace{\bar{g}_{Na}m^3h(V-E_{Na})}^{I_{Na}}-\overbrace{g_L(V-E_L)}^{I_L}
    \dot{n} &= \alpha_n(V)(10n)-\beta_n(V)n
    \dot{m} &= \alpha_m(V)(1-m)-\beta_m(V)m
    \dot{h} &= \alpha_h(V)(1-h)-\beta_h(V)h
\end{align*}$$``` 
s.t
.
```math
$$\begin{align*}
    \alpha_n(V) &= 0.01 \frac{10-V}{exp(\frac{10-p}{10})-1}\\
    \beta_n(V &= 0.125 \exp(\frac{-V}{80})\\
    \alpha_m (V) &= 0.1\quad \frac{25-V}{\exp(\frac{25-V}{10})-1}\\
    \beta_m(V) &= 4exp(\frac{-V}{18})
\end{align*}$$ ```


This header file also stores the for loop that updates the output voltage with regard to the value of injected current $I$.

