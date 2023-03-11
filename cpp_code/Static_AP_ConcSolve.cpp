#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <vector>
#include <algorithm>

using namespace std;

/*
THIS FILE IS JUST FOR ME TO TEST THINGS
*/

const char *path1="../data_files/ConcSolve_output.csv";

double F = 96485.332;

double a_n, b_n, a_h, b_h, a_m, b_m;
double n_dynamic, m_dynamic, h_dynamic, n_inf, m_inf, h_inf;
vector<double> vec_n, vec_m, vec_h;
vector<double> vec_inf_n, vec_inf_m, vec_inf_h;
vector<double> vec_tau_n, vec_tau_m, vec_tau_h;
double tau_n, tau_m, tau_h;

double V_dt; 
vector<double> vec_V, vec_Na_I, vec_K_I, vec_L_I; 
vector<double> vec_Na_conc, vec_K_conc, vec_L_conc, vec_approx_conc, vec_conc_diff;
vector<double> vec_Nap, vec_Kp;
vector<double> vec_VWT;
double Na_m, K_m, L_m, approx_m; 
double Na_conc, K_conc, L_conc, approx_conc; 

double C_m = 1; 

double g_k = 36; 
double g_Na = 120;
double g_l = 0.3; 

double E_k = -12; 
double E_Na = 120; 
double E_l = 10.6; 
double delta_t = 0.001;

double global_current = -5;

/*
n is a kinetic equation built to track ONE kind of POTASSIUM channel's opening
it will be a proportion of channels open, and or the activation of the channel
potassium current can be generalized as: g*n^4*(V - E)
Where V is resting voltage, E is membrane potential, and g is conductance

m is the same, but for Na 
h represents INACTIVATION of Na channels, so h and m compete

The formulas are taken from page 37 of the Electrophysiology textbook
*/
int reset_vecs(int x){
    vec_V.clear();
    vec_n.clear();
    vec_m.clear();
    vec_h.clear();
    vec_tau_m.clear();
    vec_inf_m.clear();
    vec_tau_h.clear();
    vec_inf_h.clear();
    vec_tau_n.clear();
    vec_inf_n.clear();
    vec_Na_I.clear();
    vec_K_I.clear();
    vec_L_I.clear();
    return(0);
}

double dynamical_h(double V){
    //for some reason I decided to add in double h, so that we can store this value local to the method
    //and, h_dynamic can be shared to main. I forget why I did this lol, I probably had something in mind
    a_h = 0.07*exp(-V/20);
    b_h = ((1)/(exp((30-V)/10) + 1));
    h_inf = a_h/(a_h + b_h);
    tau_h = 1/(a_h + b_h);
    //the if statements are added to expunge NaN from the data
    //storing the values in vectors so that they can be easily written to a csv file
    vec_tau_h.push_back(tau_h);
    vec_inf_h.push_back(h_inf);
    return(0);
}

double dynamical_n(double V){
    a_n = 0.01*((10-V)/(exp((10-V)/10) - 1));
    if(V == 10){
        a_n = 0.1; // This is the Taylor approx value for when divide by 0
    }
    b_n = 0.125*exp(-V/80);
    n_inf = a_n/(a_n + b_n);
    tau_n = 1/(a_n + b_n);
    //cout << tau_n << endl;
    vec_tau_n.push_back(tau_n);
    vec_inf_n.push_back(n_inf);
    //cout << V << endl;
    return (0);
}

double dynamical_m(double V){
    a_m = 0.1*((25 - V)/(exp((25-V)/10) - 1));
    if (V == 25){
        a_m = 1; // This is the Taylor approx value for when divide by 0
    }
    b_m = 4*exp(-V/18);
    m_inf = a_m/(a_m + b_m);
    tau_m = 1/(a_m + b_m);
    //cout << tau_m << endl;
    vec_tau_m.push_back(tau_m);
    vec_inf_m.push_back(m_inf);
    return (0);
}

double Proportion_open_test(int x){
    //This scales between membrane voltages of -40 to 100, which is near the typical operating range of neurons
    for(double i = -40; i <= 100; i++){
        double V_temp = i;
        double m_temp, h_temp, n_temp;
        //cout << V_temp << endl;
        //inputs a, b, and c are dependent on the type of channel
        //Sodium channels will only utilize a and b
        //Potassium channels will only use c
            m_temp = dynamical_m(V_temp);
            h_temp = dynamical_h(V_temp);
            n_temp = dynamical_n(V_temp);
        }
  return (0);
}

int output_file(int x){
    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    //Proportion_open_test(0);

    myfile << "Tau_n, ";
    for(int i = 0; i < vec_tau_n.size(); i++){
       myfile << vec_tau_n[i] << ","; 
    }
    myfile << "\n Tau_M, ";
    for(int i = 0; i < vec_tau_m.size(); i++){
       myfile << vec_tau_m[i] << ","; 
        //cout << vec_tau_m[i] << endl;
    }
    myfile << "\n Tau_H, ";
    for(int i = 0; i < vec_tau_h.size(); i++){
        myfile << vec_tau_h[i] << ","; 
    }
    myfile << "\n Inf_n, ";
    for(int i = 0; i < vec_inf_n.size(); i++){
       myfile << vec_inf_n[i] << ","; 
    }
    myfile << "\n Inf_M, ";
    for(int i = 0; i < vec_inf_m.size(); i++){
       myfile << vec_inf_m[i] << ","; 
    }
    myfile << "\n Inf_H, ";
    for(int i = 0; i < vec_inf_h.size(); i++){
        myfile << vec_inf_h[i] << ","; 
    }
    return(0);
}

double Static_WT_AP(int arbitrary_variable){

    double V_start = 0;
    double current;
    double V_temp;
    double Na_I_temp, K_I_temp, L_I_temp = 0;
    int x = 0; 

    vec_V.push_back(V_start);
    vec_n.push_back(dynamical_n(0));
    vec_m.push_back(dynamical_m(0));
    vec_h.push_back(dynamical_h(0));

    vec_Na_conc.push_back(0);
    vec_K_conc.push_back(0);
    vec_L_conc.push_back(0);
    vec_approx_conc.push_back(0);

    for(double i = 0; i <= 20; i += delta_t){

        //cout << "Break point 1" << endl;

        dynamical_m(vec_V[x]);
        dynamical_h(vec_V[x]);
        dynamical_n(vec_V[x]);

        //cout << "Break point 2" << endl;

        vec_n.push_back(vec_n[x] + delta_t*((vec_inf_n[x] - vec_n[x])/vec_tau_n[x]));
        vec_m.push_back(vec_m[x] + delta_t*((vec_inf_m[x] - vec_m[x])/vec_tau_m[x]));
        vec_h.push_back(vec_h[x] + delta_t*((vec_inf_h[x] - vec_h[x])/vec_tau_h[x]));

        //cout << "Break point 3" << endl;

        K_I_temp = (g_k*pow(vec_n[x+1],4)*((vec_V[x]) - E_k));
        Na_I_temp = (g_Na*pow(vec_m[x+1],3)*pow(vec_h[x+1],1)*((vec_V[x]) - E_Na));
        L_I_temp = (g_l*((vec_V[x]) - E_l));

        V_dt = (current - K_I_temp - Na_I_temp - L_I_temp)/C_m;
        vec_V.push_back(vec_V[x] + delta_t*V_dt);

        vec_VWT.push_back(vec_V[x] + delta_t*V_dt);

        vec_Na_I.push_back(Na_I_temp);
        vec_K_I.push_back(K_I_temp);
        vec_L_I.push_back(L_I_temp);

        Na_m = (Na_I_temp*0.001*delta_t*pow(10,-9))/F;
        K_m = (K_I_temp*0.001*delta_t*pow(10,-9))/F;
        L_m = (L_I_temp*0.001*delta_t*pow(10,-9))/F;

        Na_conc = 1000*(Na_m/((4/3)*M_PI*pow(0.001,3)));
        K_conc = 1000*(K_m/((4/3)*M_PI*pow(0.001,3)));
        L_conc = 1000*(L_m/((4/3)*M_PI*pow(0.001,3)));

        vec_Na_conc.push_back(vec_Na_conc[x] + Na_conc);
        vec_K_conc.push_back(vec_K_conc[x] + K_conc);
        vec_L_conc.push_back(vec_L_conc[x] + L_conc);

        approx_m = (1.257*pow(10,-12)*(vec_V[x]/1000))/F;
        approx_conc = 1000*(approx_m/((4/3)*M_PI*pow(0.001,3)));
        vec_conc_diff.push_back(vec_Na_conc[x] + vec_K_conc[x] + vec_L_conc[x]);

        vec_approx_conc.push_back(approx_conc);

        //cout << V << endl;
        //cout << V_temp << endl;
        x += 1; 
        //cout << x << endl;
    }
    return(0);
}

double output_WT_Static_AP(double x)
{
    reset_vecs(0);
    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    Static_WT_AP(0);

    vector<int> sizes;
    /*cout<<vec_V.size()<<endl;
    cout<<vec_N.size()<<endl;
    cout<<vec_tiny_N.size()<<endl;*/
    sizes.insert(sizes.begin(),vec_V.size());
    sizes.insert(sizes.begin(),vec_K_I.size());
    sizes.insert(sizes.begin(),vec_Na_I.size());
    sizes.insert(sizes.begin(),vec_L_I.size());
    sort(sizes.begin(), sizes.end());
    int max_size = sizes.back();
    
    cout << max_size << endl;

    bool V; 
    bool K_I; 
    bool Na_I; 
    bool L_I; 
    bool bool_Na_conc;
    bool bool_K_conc;
    bool bool_L_conc;
    bool bool_approx_conc;
    bool bool_conc_diff;
    
    //cout << "Break point 4" << endl;

    myfile << "V,K_I,Na_I,L_I,Na_conc,K_conc,L_conc,approx_conc,conc_diff\n";
    for (int i = 0; i < max_size; i++)
    {
        //cout << "Break point 5" << endl;
        V = (vec_V.size() > i) ? true : false;
        K_I = (vec_K_I.size() > i) ? true : false;
        Na_I = (vec_Na_I.size() > i) ? true : false;
        L_I = (vec_L_I.size() > i) ? true : false;
        bool_Na_conc = (vec_Na_conc.size() > i) ? true : false;
        bool_K_conc = (vec_K_conc.size() > i) ? true : false;
        bool_L_conc = (vec_L_conc.size() > i) ? true : false;
        bool_approx_conc = (vec_approx_conc.size() > i) ? true : false;
        bool_conc_diff = (vec_conc_diff.size() > i) ? true : false;

        //cout << "Break point 6" << endl;

        if(V) myfile << vec_V[i] <<"," ;
        if(!V) myfile <<"," ;
        if(K_I) myfile << vec_K_I[i] << ",";
        if(!K_I) myfile << ",";
        if(Na_I) myfile << vec_Na_I[i] << ",";
        if(!Na_I) myfile <<",";
        if(L_I) myfile << vec_L_I[i] << ",";
        if(!L_I) myfile <<",";

        if(bool_Na_conc) myfile << vec_Na_conc[i] << ",";
        if(!bool_Na_conc) myfile <<",";
        if(bool_K_conc) myfile << vec_K_conc[i] << ",";
        if(!bool_K_conc) myfile <<",";
        if(bool_L_conc) myfile << vec_L_conc[i] << ",";
        if(!bool_L_conc) myfile <<",";
        if(bool_approx_conc) myfile << vec_approx_conc[i]<<",";;
        if(!bool_approx_conc) myfile <<",";
        if(bool_conc_diff) myfile << vec_conc_diff[i];


        //if(!dn_V_0) myfile << vec_tiny_N[i];

        myfile << "\n";
        //cout << "Break point 6" << endl;
    }
    myfile.close();
    return (x);
}

int main(void) {
  cout << "Begin" << endl;
  //output_file(0);
  cout << "break point 1" << endl;
  output_WT_Static_AP(0);
  cout << "break point 2" << endl;
  cout << "End" << endl;
}