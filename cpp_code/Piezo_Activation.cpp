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
// normal_distribution<double> stochastic_opening(0,4);


normal_distribution<double> stiffness(0.7,0.05);
normal_distribution<double> pressure(5,5);
normal_distribution<double> voltage(-70,10);

vector<double> vec_P_Substrate;
vector<double> vec_P_Pressure;
vector<double> vec_P_Voltage;
vector<double> vec_P_Total; 

vector<double> vec_closed_inactive_local;
vector<double> vec_closed_active_local;

double open_local, closed_inactive_local, closed_active_local; 

double local_N_Piezo = 100;

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

double Piezo_Channel(int time){

    double stochastic_stiffness = abs(stiffness(generator));
    double stochastic_pressure = abs(pressure(generator));
    double stochastic_voltage = voltage(generator);

    if((time > 50) && (time < 60)){
        stochastic_pressure = 30;
    }

    //cout << stochastic_stiffness << " " << stochastic_pressure << " " << stochastic_voltage << endl;

    double P_opening_temp;
    double Pressure_input, Substrate_input, Voltage_input; 

    Pressure_input = stochastic_pressure;
    Substrate_input = stochastic_stiffness;
    Voltage_input = stochastic_voltage;

    double P_P = Piezo_P_Pressure(Pressure_input);
    double P_S = Piezo_P_Substrate(Substrate_input);
    double P_V = Piezo_P_Voltage(Voltage_input);

    double P_total = P_P*P_S + P_V;

    P_opening_temp = 1/(exp((0.5 - P_total)/0.1) + 1);

    vec_closed_inactive_local.push_back(0.904837418*vec_closed_inactive_local[time] + (vec_P_Total[time] - 0.904837418*vec_P_Total[time]));
    vec_closed_active_local.push_back((1 - P_opening_temp)*vec_closed_active_local[time] + (vec_closed_inactive_local[time] - 0.904837418*vec_closed_inactive_local[time]));

    // 0.904837418 = e^-1/10

    open_local = 0.904837418*vec_P_Total[time] + P_opening_temp*vec_closed_active_local[time];
    // cout << open_local << " and " << P_opening_temp*vec_closed_active_local[time] << endl;

    // cout << P_opening_temp << endl;
    // open_local = vec_num_open[time] + P_opening_temp*closed_active_local;

    if(open_local >= local_N_Piezo){
        open_local = local_N_Piezo;
    }

    vec_P_Total.push_back(open_local);

    return(0);
}

double output_file(double x)
{
    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    vec_P_Total.push_back(0);
    vec_closed_inactive_local.push_back(0);
    vec_closed_active_local.push_back(100);

    for(double i = 0; i <= 100; i++){
        Piezo_Channel(i);
    }

    vector<int> sizes;

    sizes.insert(sizes.begin(),vec_P_Total.size());
    sizes.insert(sizes.begin(),vec_P_Pressure.size());
    sizes.insert(sizes.begin(),vec_P_Substrate.size());
    sizes.insert(sizes.begin(),vec_P_Voltage.size());
    sizes.insert(sizes.begin(),vec_closed_active_local.size());
    sizes.insert(sizes.begin(),vec_closed_inactive_local.size());

    sort(sizes.begin(), sizes.end());
    int max_size = sizes.back();
    
    cout << max_size << endl;

    bool Bool_P_open1;
    bool Bool_P_open2;
    bool Bool_P_open3;
    bool Bool_P_total;
    bool Bool_P_active;
    bool Bool_P_inactive;

    //cout << "Break point 4" << endl;

    myfile << "Piezo_Open,Piezo_Closed,Piezo_Inactive,Pressure,Substrate,Voltage\n";

    for (int i = 0; i < max_size; i++)
    {
        //cout << "Break point 5" << endl;
        Bool_P_total = (vec_P_Total.size() > i) ? true : false;
        Bool_P_open1 = (vec_P_Pressure.size() > i) ? true : false;
        Bool_P_open2 = (vec_P_Substrate.size() > i) ? true : false;
        Bool_P_open3 = (vec_P_Voltage.size() > i) ? true : false;
        Bool_P_active = (vec_closed_active_local.size() > i) ? true : false;
        Bool_P_inactive = (vec_closed_inactive_local.size() > i) ? true : false;

        //cout << "Break point 6" << endl;

        if(Bool_P_total) myfile << vec_P_Total[i] << ",";
        if(!Bool_P_total) myfile <<",";
        if(Bool_P_active) myfile << vec_closed_active_local[i] << ",";
        if(!Bool_P_active) myfile <<",";
        if(Bool_P_inactive) myfile << vec_closed_inactive_local[i] << ",";
        if(!Bool_P_inactive) myfile <<",";
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
