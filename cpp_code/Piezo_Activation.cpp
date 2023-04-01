#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <vector>
#include <algorithm>
#include <random>

using namespace std;

const char *path1="../data_files/Piezo_Activation.csv";

default_random_engine generator;
normal_distribution<double> stochastic_opening(0,4);

vector<double> vec_P_open;

double Piezo_kinetics(double i, double j, double k, int time){
    double P_inf;

    double i_temp = (70 - i)/10;
    double j_temp = (50 - i)/10;
    double k_temp = (30 - k)/10; 


    P_inf = 1/(exp(i_temp + j_temp + k_temp) + 1);

    vec_P_open.push_back(P_inf);
        
}

double Piezo_activation(double x){
    for(int stiffness = 0; stiffness <= 100; stiffness++){
        
    }
    
}

double voltage_output(double x)
{
    Piezo_activation(0);

    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    vector<int> sizes;

    sizes.insert(sizes.begin(),vec_P_open.size());

    sort(sizes.begin(), sizes.end());
    int max_size = sizes.back();
    
    cout << max_size << endl;

    bool Bool_P_open;

    //cout << "Break point 4" << endl;

    myfile << "Ca_concentration,Buffering,Piezo_current,J_ryr,J_serca,Piezo_open\n";
    for (int i = 0; i < max_size; i++)
    {
      if(i >= 1000){
        //cout << "Break point 5" << endl;
        Bool_P_open = (vec_P_open.size() > i) ? true : false;
    
        //cout << "Break point 6" << endl;

        if(Bool_P_open) myfile << vec_P_open[i] << ",";
   
        //if(!dn_V_0) myfile << vec_tiny_N[i];

        myfile << "\n";
        //cout << "Break point 6" << endl;
      }
    }
    myfile.close();
    return (x);
}

int main(void) {
  cout << "Begin" << endl;
  voltage_output(0);
  cout << "End" << endl;
}
