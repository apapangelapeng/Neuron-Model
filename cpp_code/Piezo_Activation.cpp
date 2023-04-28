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
vector<double> vec_P_Total; 

double Piezo_P_Pressure(double i){
    double F_inf;
    F_inf = 1/(exp((30 - i)/6) + 1);
    vec_P_Pressure.push_back(F_inf);
    return(F_inf);
}

double Piezo_P_Substrate(double i){
    double S_inf;
    S_inf = (1/(0.25*pow(2*M_PI,0.5))*exp(-0.5*pow((i - 0.7)/0.25,2)))/1.6;
    vec_P_Substrate.push_back(S_inf);
    return(S_inf);
}

double Piezo_P_Voltage(double i){
    double V_inf;
    V_inf = 1/(exp((100 - i)/20) + 1);
    vec_P_Voltage.push_back(V_inf);
    return(V_inf);
}

double Piezo_Channel(double x, double y, double z){

    int stochastic = stochastic_opening(generator);

    int local_N_Piezo, open_local, closed_inactive_local, closed_active_local; 

    double P_opening_temp;
    double Pressure_input, Substrate_input, Voltage_input; 

    Pressure_input = x;
    Substrate_input = y;
    Voltage_input = z; 

    double P_P = Piezo_P_Pressure(Pressure_input);
    double P_S = Piezo_P_Substrate(Substrate_input);
    double P_V = Piezo_P_Voltage(Voltage_input);
    double P_total = P_P*P_S + P_V;

    P_opening_temp = 1/(exp((0.5 - P_total)/0.1) + 1);

    // TO DO: %%%%%%%
    // Add time screen, can add arbitrary force that acts at time x
    // Fix/add variable names
    // Make closing mechanism work 

    // int local_tau = 1.6/delta_T;
    
    // if(!(open_counter % local_tau)){
    //     closed_inactive_local = vec_num_open[time][x][y]*0.9048; //this is the time constant of Piezo, so it is not relevant on a micro s scale 
    // }

    vec_P_Total.push_back(P_opening_temp);

    return(0);
}

double output_file(double x)
{
    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    for(double i = 0; i <= 100; i++){
        Piezo_Channel(i, 0.7, 0);
    }
    for(double i = 0; i <= 3; i += 0.1){
        Piezo_Channel(60, i, 0);
    }
    for(double i = 0; i <= 200; i++){
        Piezo_Channel(0, 0, i);
    }

    vector<int> sizes;

    sizes.insert(sizes.begin(),vec_P_Total.size());
    sizes.insert(sizes.begin(),vec_P_Pressure.size());
    sizes.insert(sizes.begin(),vec_P_Substrate.size());
    sizes.insert(sizes.begin(),vec_P_Voltage.size());

    sort(sizes.begin(), sizes.end());
    int max_size = sizes.back();
    
    cout << max_size << endl;

    bool Bool_P_open1;
    bool Bool_P_open2;
    bool Bool_P_open3;
    bool Bool_P_total;

    //cout << "Break point 4" << endl;

    myfile << "Piezo,Pressure,Substrate,Voltage\n";

    for (int i = 0; i < max_size; i++)
    {
        //cout << "Break point 5" << endl;
        Bool_P_total = (vec_P_Total.size() > i) ? true : false;
        Bool_P_open1 = (vec_P_Pressure.size() > i) ? true : false;
        Bool_P_open2 = (vec_P_Substrate.size() > i) ? true : false;
        Bool_P_open3 = (vec_P_Voltage.size() > i) ? true : false;

        //cout << "Break point 6" << endl;

        if(Bool_P_total) myfile << vec_P_Total[i] << ",";
        if(!Bool_P_total) myfile <<",";
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
