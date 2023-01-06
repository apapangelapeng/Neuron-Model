#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <vector>
#include <map>
#include<algorithm>
#include <random>

using namespace std;

/*
THIS FILE IS JUST FOR ME TO TEST THINGS
*/

const char *path1 = "../data_files/reduced_stochastic_V_output.csv";


vector<double> vec_V, vec_N, vec_Space,vec_tiny_N,vec_N_Zerodt,vec_V_Zerodt;

vector<double> vec_V_Zero_dvdt, vec_N_Zero_dvdt, vec_V_Zero_dndt, vec_N_Zero_dndt;

vector<double> vec_Vdt, vec_Ndt, vec_x_Vdt, vec_y_Ndt;

vector<vector<double> > v_map_2d;

double V_dt, N_dt;

double v_threshold = 0.05;
double v_max = 1;
double V_start = 0;
double v_range = 0.5;

double N_start = 0;

double gam = 3;
double gam_range = 0;

double e = 0.005;

//double current = 0;
double current_range = 0.5;
double current_temp;

double diffusion = 0.5;
double diffusion_range = 2;

double spatial_d = 0;
double delta_x = 0.1;
double x_range = 1;

double time_d = 0;
double delta_t = 0.1;

int reset_vecs(int x)
{
    vec_N.clear();
    vec_V.clear();
    return (0);
}


void write_Reduced_AP(vector<vector<double> > & v_map){
    cout << "I HAVE BEEN SUMMONED " << endl;
    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    /*
    for (V_start = -0.5; V_start <= v_range; V_start += 0.3)
    {

        string end = (V_start >= v_range - 0.3) ? "\n" : ",";
        myfile << "V" << V_start << end;
    }
    */

    myfile << "V" << V_start << "\n";
    for (int i = 0; i < vec_V.size(); i++)
    {
        for (int j = 0; j < v_map.size(); j++)
        {
            string end = (j == v_map.size() - 1) ? "\n" : ",";
            myfile << v_map[j][i] << end;
        }

        // cout << vec_V[i] << endl;
    }

    myfile << "\n";

}

double Reduced_AP(double z)
{
    reset_vecs(0);
    int x = 0;

    V_start = 0;

    vec_V.push_back(V_start);
    vec_N.push_back(N_start);
    default_random_engine generator;
    normal_distribution<double> distribution_d(0,0.01);
    normal_distribution<double> distribution_N(0,0.001);
    for (double i = 0; i <= 50; i += delta_t)
    {
        
        if(i >= 5 && i <= 10){
            current_temp = 0.05;
        }
        else{
            current_temp = 0;
        }
        
        cout << current_temp << endl; 

        // cout << "Break Point 4" << endl;
        double noise_d = distribution_d(generator);
        double noise_N = distribution_N(generator);
        time_d = vec_V[x] * (vec_V[x] - v_threshold) * (v_max - vec_V[x]);
        V_dt = current_temp + time_d+noise_d - vec_N[x]-noise_N;

        cout << "V_DT " << V_dt << endl;
        vec_V.push_back(vec_V[x] + V_dt);
        N_dt = e * (vec_V[x] - (gam * vec_N[x]));
        vec_N.push_back(vec_N[x] + N_dt);

        // cout << "Break Point 5" << endl;

        // cout << "x = " << x << endl;

        x += 1;
        }

    vector<double> vec_V_deep = vec_V;
    v_map_2d.push_back(vec_V_deep); 

    cout << "V_start: " << V_start << endl;

    write_Reduced_AP(v_map_2d);
    return (0);
}


int main(void)
{
    cout << "Begin" << endl;

    Reduced_AP(0);
    //output__nullcline_file(0);
    //Diffusion_AP(0);

    cout << "End" << endl;
}
