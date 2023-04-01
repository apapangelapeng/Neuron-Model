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

vector<double> vec_P_open1;
vector<double> vec_P_open2;
vector<double> vec_P_open3;

double Piezo_kinetics(double i, double j, double k, int time){
    double P_inf;

    double i_temp = (70 - i)/10;
    double j_temp = (1 - j)/10;
    double k_temp = (30 - k)/10; 


    P_inf = 1/(exp(i_temp + j_temp + k_temp) + 1);

    if(i > 0){vec_P_open1.push_back(P_inf);}
    // if(j > 0){vec_P_open2.push_back(P_inf);}
    // if(k > 0){vec_P_open3.push_back(P_inf);}

    return(0);
}

double Piezo_activation(double x){
    for(int temp = 1; temp <= 100; temp++){
        Piezo_kinetics(temp, 20, 10, 0);
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

    sizes.insert(sizes.begin(),vec_P_open1.size());
    sizes.insert(sizes.begin(),vec_P_open2.size());
    sizes.insert(sizes.begin(),vec_P_open3.size());

    sort(sizes.begin(), sizes.end());
    int max_size = sizes.back();
    
    cout << max_size << endl;

    bool Bool_P_open1;
    bool Bool_P_open2;
    bool Bool_P_open3;

    //cout << "Break point 4" << endl;

    myfile << "1,2,3\n";

    for (int i = 0; i < max_size; i++)
    {
        //cout << "Break point 5" << endl;
        Bool_P_open1 = (vec_P_open1.size() > i) ? true : false;
        Bool_P_open2 = (vec_P_open2.size() > i) ? true : false;
        Bool_P_open3 = (vec_P_open3.size() > i) ? true : false;

        //cout << "Break point 6" << endl;

        if(Bool_P_open1) myfile << vec_P_open1[i] << ",";
        if(!Bool_P_open1) myfile <<",";
        if(Bool_P_open2) myfile << vec_P_open2[i] << ",";
        if(!Bool_P_open2) myfile <<",";
        if(Bool_P_open3) myfile << vec_P_open3[i];

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
