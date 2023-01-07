#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <vector>
#include <algorithm>
#include "StaticConst.h"
#include "DynamicFunc.h"
#include "StaticOutput.h"

using namespace std;

/*
THIS FILE IS JUST FOR ME TO TEST THINGS
*/

const char *output_path="../data_files/static_wt_output.csv";

double a_n, b_n, a_h, b_h, a_m, b_m;
double n_dynamic, m_dynamic, h_dynamic, n_inf, m_inf, h_inf;
vector<double> vec_n, vec_m, vec_h;

double tau_n, tau_m, tau_h;

double V_dt; 
vector<double> vec_V, vec_Na_I, vec_K_I, vec_L_I, vec_HCN_I, vec_Cl_I; 
vector<double> vec_Nap, vec_Kp, vec_HCNp;
vector<double> vec_VWT, vec_VHCN1, vec_VHCN2, vec_VHCN3; 
vector<double> vec_HCN_I1, vec_HCN_I2, vec_HCN_I3; 


/*
n is a kinetic equation built to track ONE kind of POTASSIUM channel's opening
it will be a proportion of channels open, and or the activation of the channel
potassium current can be generalized as: g*n^4*(V - E)
Where V is resting voltage, E is membrane potential, and g is conductance

m is the same, but for Na 
h represents INACTIVATION of Na channels, so h and m compete

The formulas are taken from page 37 of the Electrophysiology textbook
*/


double Static_WT_AP(int arbitrary_variable){
    vector<double> vec_tau_n, vec_tau_m, vec_tau_h;
    vector<double> vec_inf_n, vec_inf_m, vec_inf_h;
    double V_start = 0;
    double current;
    double V_temp;
    double Na_I_temp, K_I_temp, L_I_temp, HCN_I_temp, Cl_I_temp = 0;
    int x = 0; 

    vec_V.push_back(V_start);
    vec_n.push_back(0);
    vec_m.push_back(0);
    vec_h.push_back(0);

    //vec_Cl_I.push_back(0);

    for(double i = 0; i <= 50; i += delta_t){
        //||  10 <= i <= 12 || 20 <= i <= 22 || 30 <= i <= 32 || 40 <= i <= 42
        if((i <= 25 && i >= 23) || (i <= 28 && i >= 26) || (i <= 31 && i >= 29) || (i <= 34 && i >= 32)){
            current = global_current;
            //cout << x << endl; 
        }
        else{
            current = 0;
        }

        //cout << "Break point 1" << endl;

        dynamical_m(vec_V[x],vec_tau_m,vec_inf_m);
        dynamical_h(vec_V[x],vec_tau_h,vec_inf_h);
        dynamical_n(vec_V[x],vec_tau_n,vec_inf_n);

        //cout << "Break point 2" << endl;

        vec_n.push_back(vec_n[x] + delta_t*((vec_inf_n[x] - vec_n[x])/vec_tau_n[x]));
        vec_m.push_back(vec_m[x] + delta_t*((vec_inf_m[x] - vec_m[x])/vec_tau_m[x]));
        vec_h.push_back(vec_h[x] + delta_t*((vec_inf_h[x] - vec_h[x])/vec_tau_h[x]));

        //cout << "Break point 3" << endl;

        K_I_temp = (g_k*pow(vec_n[x+1],4)*((vec_V[x]) - E_k));
        Na_I_temp = (g_Na*pow(vec_m[x+1],3)*pow(vec_h[x+1],1)*((vec_V[x]) - E_Na));
        L_I_temp = (g_l*((vec_V[x]) - E_l));

        //Cl_I_temp = t_Cl*((g_Cl*pow(2.71828,(vec_V[x]-E_Cl)/b_Cl))*(vec_V[x] - E_Cl) - vec_Cl_I[x]);

        //cout << "Break point 4" << endl;

        //cout << V << endl;

        V_dt = (current - K_I_temp - Na_I_temp - L_I_temp)/C_m;
        vec_V.push_back(vec_V[x] + delta_t*V_dt);

        vec_VWT.push_back(vec_V[x] + delta_t*V_dt);

        vec_Na_I.push_back(Na_I_temp);
        vec_K_I.push_back(K_I_temp);
        vec_L_I.push_back(L_I_temp);
        vec_Cl_I.push_back(Cl_I_temp);


        //cout << V << endl;
        //cout << V_temp << endl;
        x += 1; 
        //cout << x << endl;
    }
    vec_Cl_I.erase(vec_Cl_I.begin());
    output_data(output_path,vec_V,vec_K_I,vec_Na_I,vec_L_I,vec_Cl_I,vec_HCN_I,"None");
    return(0);
}
/*
double output_WT_Static_AP(double x)
{
    ofstream create_file(output_path);
    ofstream myfile;
    myfile.open(output_path);

    Static_WT_AP(0);

    vector<int> sizes;
    cout<<vec_V.size()<<endl;
    cout<<vec_N.size()<<endl;
    cout<<vec_tiny_N.size()<<endl;
    sizes.insert(sizes.begin(),vec_V.size());
    sizes.insert(sizes.begin(),vec_K_I.size());
    sizes.insert(sizes.begin(),vec_Na_I.size());
    sizes.insert(sizes.begin(),vec_L_I.size());
    sizes.insert(sizes.begin(),vec_Cl_I.size());
    sort(sizes.begin(), sizes.end());
    int max_size = sizes.back();
    
    cout << max_size << endl;

    bool V; 
    bool K_I; 
    bool Na_I; 
    bool L_I; 
    bool Cl_I;
    
    //cout << "Break point 4" << endl;

    myfile << "V,K_I,Na_I,L_I,Cl_I\n";
    for (int i = 0; i < max_size; i++)
    {
        //cout << "Break point 5" << endl;
        V = (vec_V.size() > i) ? true : false;
        K_I = (vec_K_I.size() > i) ? true : false;
        Na_I = (vec_Na_I.size() > i) ? true : false;
        L_I = (vec_L_I.size() > i) ? true : false;
        Cl_I = (vec_Cl_I.size() > i) ? true : false;

        //cout << "Break point 6" << endl;

        if(V) myfile << vec_V[i] <<"," ;
        if(!V) myfile <<"," ;
        if(K_I) myfile << vec_K_I[i] << ",";
        if(!K_I) myfile << ",";
        if(Na_I) myfile << vec_Na_I[i] << ",";
        if(!Na_I) myfile <<",";
        if(L_I) myfile << vec_L_I[i] << ",";
        if(!L_I) myfile <<",";
        if(Cl_I) myfile << vec_Cl_I[i];

        //if(!dn_V_0) myfile << vec_tiny_N[i];

        myfile << "\n";
        //cout << "Break point 6" << endl;
    }
    myfile.close();
    return (x);
}*/


int main(void) {
  cout << "Begin" << endl;
Static_WT_AP(0);

  cout << "End" << endl;
}
