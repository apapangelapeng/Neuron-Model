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

vector<double> vec_P_Substrate;
vector<double> vec_P_Pressure;
vector<double> vec_P_Voltage;

double Piezo_P_Pressure(double i){
    double F_inf;
    double force; 
    
    F_inf = 1/(exp((30 - i)/10) + 1);

    vec_P_Pressure.push_back(F_inf);

    return(0);
}

double Piezo_P_Substrate(double i){
    double S_inf;

    S_inf = (1/(0.25*pow(2*M_PI,0.5))*exp(-0.5*pow((i - 0.7)/0.25,2)));

    vec_P_Substrate.push_back(S_inf);

    return(0);
}

double Piezo_P_Voltage(double i){
    double V_inf;

    V_inf = 1/(exp((100 - i)/20) + 1);

    vec_P_Voltage.push_back(V_inf);

    return(0);
}

double Piezo_activation(double x){

    for(int temp = 0; temp <= 100; temp++){
        Piezo_P_Pressure(temp);
    }

    for(double temp = 0; temp <= 3; temp += 0.01){
        Piezo_P_Substrate(temp);
    }

    for(double temp = 0; temp <= 150; temp += 1){
        Piezo_P_Voltage(temp);
    }


    // for(int temp = 1; temp <= 100; temp++){
    //     Piezo_kinetics(0, temp, 0, 0);
    // }

    // for(int temp = 1; temp <= 100; temp++){
    //     Piezo_kinetics(0, 0, temp, 0);
    // }

    return(0);   
}

double output_file(double x)
{
    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    Piezo_activation(0);

    vector<int> sizes;

    sizes.insert(sizes.begin(),vec_P_Pressure.size());
    sizes.insert(sizes.begin(),vec_P_Substrate.size());
    sizes.insert(sizes.begin(),vec_P_Voltage.size());

    sort(sizes.begin(), sizes.end());
    int max_size = sizes.back();
    
    cout << max_size << endl;

    bool Bool_P_open1;
    bool Bool_P_open2;
    bool Bool_P_open3;

    //cout << "Break point 4" << endl;

    myfile << "Pressure,Substrate,Voltage\n";

    for (int i = 0; i < max_size; i++)
    {
        //cout << "Break point 5" << endl;
        Bool_P_open1 = (vec_P_Pressure.size() > i) ? true : false;
        Bool_P_open2 = (vec_P_Substrate.size() > i) ? true : false;
        Bool_P_open3 = (vec_P_Voltage.size() > i) ? true : false;

        //cout << "Break point 6" << endl;

        if(Bool_P_open1) myfile << vec_P_Pressure[i] << ",";
        if(!Bool_P_open1) myfile <<",";
        if(Bool_P_open2) myfile << vec_P_Substrate[i] << ",";
        if(!Bool_P_open2) myfile <<",";
        if(Bool_P_open3) myfile << vec_P_Voltage[i];

        //if(!dn_V_0) myfile << vec_tiny_N[i];

        myfile << "\n";
        //cout << "Break point 6" << endl;
    }

    myfile.close();
    return (x);
}

int main(void) {
  cout << "Begin" << endl;
  output_file(0);
  cout << "End" << endl;
}
