#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <vector>
#include <algorithm>
#include <random>

using namespace std;

/*
THIS FILE IS JUST FOR ME TO TEST THINGS
*/

const char *path1="../data_files/HH_propagate_output.csv";


double a_n, b_n, a_h, b_h, a_m, b_m;
double n_dynamic, m_dynamic, h_dynamic, n_inf, m_inf, h_inf;
vector<double> vec_n, vec_m, vec_h;
vector<double> vec_inf_n, vec_inf_m, vec_inf_h;
vector<double> vec_tau_n, vec_tau_m, vec_tau_h;
double tau_n, tau_m, tau_h;


vector<vector<double> > n_2d, m_2d, h_2d;
vector<vector<double> > inf_n_2d, inf_m_2d, inf_h_2d;
vector<vector<double> > tau_n_2d, tau_m_2d, tau_h_2d;


double V_dt; 

vector<double> vec_V, vec_Na_I, vec_K_I, vec_L_I;
vector<double> vec_Nap, vec_Kp;
vector<double> vec_VWT; 

vector<vector<double> > v_WT_2d;

double C_m = 1; 

double g_k = 36; 
double g_Na = 120;
double g_l = 0.3; 

double E_k = -12; 
double E_Na = 120; 
double E_l = 10.6; 

double delta_t = 0.01;
double delta_x = 0.01;

double current = 0;

double x_range = 1;

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

void output_file(vector<vector<double> > v_map){
    cout << "I HAVE BEEN SUMMONED " << endl;
    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    int col_num = x_range/delta_x;

    //myfile << "V" << V_start << "\n";
    for (int j = 0; j < v_map.size(); j++)
    {
        for (int i = 0; i < col_num; i++)
        {
            if (i == 0)
            {
                myfile << v_map[j][i];
            }
            else
            {
                myfile << "," << v_map[j][i];
            }
            
        }
        myfile << "\n";
    }
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

double Static_WT_AP(int arbitrary_variable){
    double V_start = 0;
    double V_temp;
    double Na_I_temp, K_I_temp, L_I_temp, HCN_I_temp, Cl_I_temp = 0;
    int counter_time = 0; 
    int counter_space = 0;

    for (double space = 0; space <= x_range+delta_x; space += delta_x) //THIS COUNTS SPACE!!
    { 
        vec_V.push_back(V_start);
        vec_n.push_back(0);
        vec_m.push_back(0);
        vec_h.push_back(0);
    }
    v_WT_2d.push_back(vec_V);
    n_2d.push_back(vec_n);
    m_2d.push_back(vec_m);
    h_2d.push_back(vec_h);

    for(double t = 0; t <= 10; t += delta_t){
        counter_space = 0;
        reset_vecs(0);
        for(double x = 0; x <= x_range; x += delta_x){
        //||  10 <= i <= 12 || 20 <= i <= 22 || 30 <= i <= 32 || 40 <= i <= 42

        //cout << "Break point 1" << endl;

        dynamical_m(v_WT_2d[counter_time][counter_space]);
        dynamical_h(v_WT_2d[counter_time][counter_space]);
        dynamical_n(v_WT_2d[counter_time][counter_space]);

        //cout << "Break point 2" << endl;

        vec_n.push_back(n_2d[counter_time][counter_space] + delta_t*((vec_inf_n[counter_space] - n_2d[counter_time][counter_space])/vec_tau_n[counter_space]));
        vec_m.push_back(m_2d[counter_time][counter_space] + delta_t*((vec_inf_m[counter_space] - m_2d[counter_time][counter_space])/vec_tau_m[counter_space]));
        vec_h.push_back(h_2d[counter_time][counter_space] + delta_t*((vec_inf_h[counter_space] - h_2d[counter_time][counter_space])/vec_tau_h[counter_space]));

        //cout << "Break point 3" << endl;

        K_I_temp = (g_k*pow(vec_n[counter_time+1],4)*((vec_V[counter_time]) - E_k));
        Na_I_temp = (g_Na*pow(vec_m[counter_time+1],3)*pow(vec_h[counter_time+1],1)*((vec_V[counter_time]) - E_Na));
        L_I_temp = (g_l*((vec_V[counter_time]) - E_l));


        //cout << "Break point 4" << endl;

        //cout << V << endl;

        V_dt = (current - K_I_temp - Na_I_temp - L_I_temp)/C_m;
        vec_V.push_back(v_WT_2d[counter_time][counter_space] + delta_t*V_dt);

        vec_Na_I.push_back(Na_I_temp);
        vec_K_I.push_back(K_I_temp);
        vec_L_I.push_back(L_I_temp);

        //cout << V << endl;
        //cout << V_temp << endl;
        counter_space += 1;

        if(counter_space == 1){
            //cout << "K_I_temp " << K_I_temp << endl;
            cout << "vec_n " << vec_n[counter_space] << endl;
        }

        }

        counter_time += 1; 

        n_2d.push_back(vec_n);
        m_2d.push_back(vec_m);
        h_2d.push_back(vec_h);

        inf_n_2d.push_back(vec_inf_n);
        inf_m_2d.push_back(vec_inf_m);
        inf_h_2d.push_back(vec_inf_h);

        tau_n_2d.push_back(vec_tau_n);
        tau_m_2d.push_back(vec_tau_m);
        tau_h_2d.push_back(vec_tau_h);

        v_WT_2d.push_back(vec_V);

        //cout << "time = " << counter_time << endl;

    }

    //cout << "Break point 5" << endl;

    output_file(v_WT_2d);

    return(0);
}

int main(void) {
  cout << "Begin" << endl;
  Static_WT_AP(0);
  cout << "End" << endl;
}
