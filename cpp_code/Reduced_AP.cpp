#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <vector>

using namespace std;

/*
THIS FILE IS JUST FOR ME TO TEST THINGS
*/

const char *path1="../data_files/reduced_output.csv";

vector<double> vec_V, vec_N; 

double V_dt, N_dt;

double v_threshold = 0.3;
double v_max = 1;

double gam = 1;
double gam_range = 0;

double e = 0.001;

double current_range = 0.5;
double current_temp;

double diffusion = 0; 
double diffusion_range = 2; 

double spatial_d = 0; 
double delta_x = 0.1;

double time_d = 0; 
double delta_t = 0.1;

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
    vec_N.clear();
    return(0);
}

double Reduced_AP(double z){
    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    //cout << "Break Point 1" << endl;

    double V_start = 0;
    double N_start = 0;

    //cout << "Break Point 3" << endl;

    for(diffusion = 0; diffusion <= diffusion_range; diffusion += 0.5){
        double current = 0.05; 
        reset_vecs(0);
        int x = 0;
        vec_V.push_back(V_start);
        vec_N.push_back(N_start);

        for(double i = 0; i <= 10; i += delta_t){
            if(i <= 10){
                current_temp = current;
            }
            else{
                current_temp = 0;
            }

            //cout << "Break Point 4" << endl;


            time_d = vec_V[x]*(vec_V[x] - v_threshold)*(v_max - vec_V[x]);
            spatial_d = diffusion*(vec_V[x] - (2*vec_V[x-1]))*pow(delta_x,2);
            V_dt = current + time_d - vec_N[x] + spatial_d;
            vec_V.push_back(vec_V[x] + V_dt);
            N_dt = e*(vec_V[x] - (gam*vec_N[x]));
            vec_N.push_back(vec_N[x] + N_dt);

            //cout << "Break Point 5" << endl;

            //cout << "x = " << x << endl;

            x += 1; 

        }

        myfile << "V" << diffusion << ",";
        for(int i = 0; i < vec_V.size(); i++){
            myfile << vec_V[i] << ","; 
            //cout << vec_V[i] << endl;
        }

        myfile << "\n";
        myfile << "N" << diffusion << ",";
        for(int i = 0; i < vec_N.size(); i++){
            myfile << vec_N[i] << ","; 
            //cout << vec_N[i] << endl;
            //cout << "Break Point 6" << endl;
        }
        myfile << "\n";

        cout << "Diffusion: " << diffusion << endl;
    }
    myfile.close();
    return(0);
}

int main(void) {
/*
Program workflow is main -calls-> ouput_file -calls-> Porportion_open for Sodium and Potassium individually
Proportion_open -calls-> the dynamical variables, which stores the values in vectors, and also calculates
the actual proportion open. 
Output_file then writes then info to TestingDynamicVars.csv
*/
  cout << "Begin" << endl;

  Reduced_AP(0);

  cout << "End" << endl;
}
