#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <vector>

using namespace std;

/*
THIS FILE IS JUST FOR ME TO TEST THINGIES
*/

double a_n, b_n, a_h, b_h, a_m, b_m;
double n_dynamic, m_dynamic, h_dynamic, n_inf, m_inf, h_inf;
vector<double> vec_n, vec_m, vec_h;
vector<double> vec_inf_n, vec_inf_m, vec_inf_h;
vector<double> vec_Nap, vec_Kp;
vector<double> vec_tau_n, vec_tau_m, vec_tau_h;
double tau_n, tau_m, tau_h;

/*
n is a kinetic equation built to track ONE kind of POTASSIUM channel's opening
it will be a proportion of channels open, and or the activation of the channel
potassium current can be generalized as: g*n^4*(V - E)
Where V is resting voltage, E is membrane potential, and g is conductance

m is the same, but for Na 
h represents INACTIVATION of Na channels, so h and m compete

The formulas are taken from page 37 of the Electrophysiology textbook
*/

double dynamical_h(double V){
    //for some reason I decided to add in double h, so that we can store this value local to the method
    //and, h_dynamic can be shared to main. I forget why I did this lol, I probably had something in mind
    double h;
    a_h = 0.07*exp(-V/20);
    b_h = ((1)/(exp((30-V)/10) + 1));
    h_inf = a_h/(a_h + b_h);
    tau_h = 1/(a_h + b_h);
    h_dynamic = (h_inf - h)/tau_h;
    h = h_dynamic;
    //the if statements are added to expunge NaN from the data
    if (tau_h != tau_h){
        tau_h = 0; 
    }
    if (h != h){
        h = 0; 
    }
    //storing the values in vectors so that they can be easily written to a csv file
    vec_h.push_back(h);
    vec_tau_h.push_back(tau_h);
    vec_inf_h.push_back(h_inf);
    //cout << "h: " << h << endl;
    return(h_inf);
}

double dynamical_n(double V){
    double n;
    a_n = 0.01*((10 - V)/(exp((10-V)/10) - 1));
    b_n = 0.125*exp(-V/80);
    n_inf = a_n/(a_n + b_n);
    tau_n = 1/(a_n + b_n);
    n_dynamic = (n_inf - n)/tau_n;
    n = n_dynamic;
    if (tau_n != tau_n){
        tau_n = 0; 
    }
    if (n != n){
        n = 0; 
    }
    vec_n.push_back(n);
    vec_tau_n.push_back(tau_n);
    vec_inf_n.push_back(n_inf);
    //cout << "n: " << n << endl;
    return (n_inf);
}

double dynamical_m(double V){
    double m;
    a_m = 0.01*((25 - V)/(exp((25-V)/10) - 1));
    b_m = 4*exp(-V/18);
    m_inf = a_m/(a_m + b_m);
    tau_m = 1/(a_m + b_m);
    m_dynamic = (m_inf - m)/tau_m;
    m = m_dynamic;
    if (tau_m != tau_m){
        tau_m = 0; 
    }
    if (m != m){
        m = 0; 
    }
    vec_m.push_back(m);
    vec_tau_m.push_back(tau_m);
    vec_inf_m.push_back(m_inf);
    //cout << "m: " << m << endl;
    return (m_inf);
}

double Proportion_open(int a, int b, int c){
    double V_temp = 0;
    //This scales between membrane voltages of -40 to 100, which is near the typical operating range of neurons
    for(int i = -40; i <= 100; i++){
        double V_temp = i;
        //cout << V_temp << endl;
        double m_temp, h_temp, n_temp = 1;
        //inputs a, b, and c are dependent on the type of channel
        //Sodium channels will only utilize a and b
        //Potassium channels will only use c
        if(a != 0){
            m_temp = dynamical_m(V_temp);
        }
        if(b != 0){
            h_temp = dynamical_h(V_temp);
        }
        if(c != 0){
            n_temp = dynamical_n(V_temp);
        }
        double temp_P = pow(m_temp,a)*pow(h_temp,b)*pow(n_temp,c);
        //cout << temp_P << endl;
        if (temp_P != temp_P){
            temp_P = 0; 
        }
        if(a != 0){
            vec_Nap.push_back(temp_P);
        }
        if(c != 0){
            vec_Kp.push_back(temp_P);
        }
    }
  return (0);
}

int output_file(int x){
    const char *path="../data_files/test_output.csv";
     ofstream create_file(path);
    ofstream myfile;
    myfile.open(path);
    cout << path;

    Proportion_open(3,1,0); //This will call the proportion method to calculate for SODIUM!!
    Proportion_open(0,0,4); //Same but for POTASSIUM!!!
    myfile << "\n N, ";
    for(int i = 0; i < vec_n.size(); i++){
       myfile << vec_n[i] << ","; 
    }
    myfile << "\n M, ";
    for(int i = 0; i < vec_m.size(); i++){
       myfile << vec_m[i] << ","; 
    }
    myfile << "\n H, ";
    for(int i = 0; i < vec_h.size(); i++){
        myfile << vec_h[i] << ","; 
    }
    myfile << "\n Tau_n, ";
    for(int i = 0; i < vec_tau_n.size(); i++){
       myfile << vec_tau_n[i] << ","; 
    }
    myfile << "\n Tau_M, ";
    for(int i = 0; i < vec_tau_m.size(); i++){
       myfile << vec_tau_m[i] << ","; 
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
    myfile << "\n Nap, ";
    for(int i = 0; i < vec_Nap.size(); i++){
        myfile << vec_Nap[i] << ","; 
    }
    myfile << "\n Kp, ";
    for(int i = 0; i < vec_Kp.size(); i++){
        myfile << vec_Kp[i] << ","; 
        //cout << vec_Kp[i] << endl;
    }
    return(0);
}

int main(void) {
/*
Program workflow is main -calls-> ouput_file -calls-> Porportion_open for Sodium and Potassium individually
Proportion_open -calls-> the dynamical variables, which stores the values in vectors, and also calculates
the actual proportion open. 
Output_file then writes then info to TestingDynamicVars.csv
*/
  cout << "Ran" << endl;
  output_file(0);
}
