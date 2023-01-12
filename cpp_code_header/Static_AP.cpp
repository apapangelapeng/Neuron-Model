#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <vector>
#include "StaticConst.h"
#include "StringPercision.h"
#include "DynamicFunc.h"
#include"DynamicOutput.h"

using namespace std;




const char *path1 = "../data_files/static_output.csv";
string path2;

double a_n, b_n, a_h, b_h, a_m, b_m;
double n_dynamic, m_dynamic, h_dynamic, n_inf, m_inf, h_inf;
vector<double> vec_n, vec_m, vec_h;
vector<double> vec_inf_n, vec_inf_m, vec_inf_h;
vector<double> vec_Nap, vec_Kp;

double tau_n, tau_m, tau_h;

double V_dt;
vector<double> vec_V, vec_Na_I, vec_K_I, vec_L_I;

double current;

/*
n is a kinetic equation built to track ONE kind of POTASSIUM channel's opening
it will be a proportion of channels open, and or the activation of the channel
potassium current can be generalized as: g*n^4*(V - E)
Where V is resting voltage, E is membrane potential, and g is conductance

m is the same, but for Na
h represents INACTIVATION of Na channels, so h and m compete

The formulas are taken from page 37 of the Electrophysiology textbook
*/
int reset_vecs(int x)
{
    vec_V.clear();
    vec_n.clear();
    vec_m.clear();
    vec_h.clear();
  
    vec_inf_m.clear();
  
    vec_inf_h.clear();
   
    vec_inf_n.clear();
    return (0);
}



double Static_AP(int arbitrary_variable)
{
    reset_vecs(0);
    vector<double> vec_tau_n, vec_tau_m, vec_tau_h;
    vector<double> vec_inf_n, vec_inf_m, vec_inf_h;
    double V_start = 0;

    double Na_I_temp, K_I_temp, L_I_temp;
    int x = 0;

    vec_V.push_back(V_start);
    vec_n.push_back(0);
    vec_m.push_back(0);
    vec_h.push_back(0);
    
    for (double i = 0; i <= 10; i += delta_t)
    {
        double local_cur;
        if (i <= 4 && i >= 2)
        {
            local_cur = current;
        }
        else
        {
            local_cur = 0;
        }

        // cout << "Break point 1" << endl;
        dynamical_m(vec_V[x],vec_tau_m,vec_inf_m);
        dynamical_h(vec_V[x],vec_tau_h,vec_inf_h);
        dynamical_n(vec_V[x],vec_tau_n,vec_inf_n);


        // cout << "Break point 2" << endl;

        vec_n.push_back(vec_n[x] + delta_t * ((vec_inf_n[x] - vec_n[x]) / vec_tau_n[x]));
        vec_m.push_back(vec_m[x] + delta_t * ((vec_inf_m[x] - vec_m[x]) / vec_tau_m[x]));
        vec_h.push_back(vec_h[x] + delta_t * ((vec_inf_h[x] - vec_h[x]) / vec_tau_h[x]));

        // cout << "Break point 3" << endl;

        K_I_temp = (g_k * pow(vec_n[x + 1], 4) * ((vec_V[x]) - E_k));
        Na_I_temp = (g_Na * pow(vec_m[x + 1], 3) * pow(vec_h[x + 1], 1) * ((vec_V[x]) - E_Na));
        L_I_temp = (g_l * ((vec_V[x]) - E_l));

        // cout << V << endl;

        V_dt = (current - K_I_temp - Na_I_temp - L_I_temp) / C_m;
        vec_V.push_back(vec_V[x] + delta_t * V_dt);

        vec_Na_I.push_back(Na_I_temp);
        vec_K_I.push_back(K_I_temp);
        vec_L_I.push_back(L_I_temp);
        // cout << V << endl;
        // cout << V_temp << endl;
        x += 1;
        // cout << x << endl;
    }
    return (0);
}

double Run_time(double x)
{
    ofstream create_file(path2);
    ofstream myfile;
    myfile.open(path2);

    reset_vecs(0);
    Static_AP(0);

    myfile << "V,";
    for (int i = 0; i < vec_V.size(); i++)
    {
        myfile << vec_V[i] << ",";
        // cout << vec_V[i] << endl;
    }
    myfile << "\n K_I, ";
    for (int i = 0; i < vec_K_I.size(); i++)
    {
        myfile << vec_K_I[i] << ",";
        // cout << vec_K_I[i] << endl;
    }
    myfile << "\n Na_I,";
    for (int i = 0; i < vec_Na_I.size(); i++)
    {
        myfile << vec_Na_I[i] << ",";
        // cout << vec_V[i] << endl;
    }
    myfile << "\n L_I,";
    for (int i = 0; i < vec_L_I.size(); i++)
    {
        myfile << vec_L_I[i] << ",";
        // cout << vec_V[i] << endl;
    }
    myfile.close();
    return (x);
}

int main(int argc, char *argv[])
{
    /*
    Program workflow is main -calls-> ouput_file -calls-> Porportion_open for Sodium and Potassium individually
    Proportion_open -calls-> the dynamical variables, which stores the values in vectors, and also calculates
    the actual proportion open.
    Output_file then writes then info to TestingDynamicVars.csv
    */
    if (argc < 2)
    {
        throw std::invalid_argument("Did not pass in a valid current argument");
        return 1;
    }
    try
    {
        // Block of code to try
        current = stod(argv[1]); // grab the first val and convert to double
        
        path2 = string("../data_files/testV_output_")+ to_string_with_precision(current,2)+ string(".csv");
    }
    catch (invalid_argument)
    {
        throw std::invalid_argument("The given command line input is not a double for current.");
    }

    cout << "Begin" << endl;
    output_file(path1);
    Run_time(0);
    cout << "End" << endl;
}
