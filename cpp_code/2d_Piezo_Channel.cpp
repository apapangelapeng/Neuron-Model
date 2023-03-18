#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <vector>
#include <algorithm>
#include <random>

#include "2d_Calcium_Dynamics.h"

using namespace std;

const char *path1="../data_files/2d_Piezo_Channel.csv";
const char *path2="../data_files/2d_Piezo_Channel_test.csv";

default_random_engine generator;
normal_distribution<double> stochastic_opening(0,0.6);

int reset_vecs(int x){
    vec_x.clear();
    vec_y.clear();
    vec_coords.clear();

    return(0);
}

vector<vector<double> > fill_2dvecs(int x_max, int y_max, double value){

    for(double i = 0; i <= x_max; i++){
        for(double j = 0; j <= y_max; j++){ 
            vec_y.push_back(value);
        }
    vec_coords.push_back(vec_y);
    }


    return(vec_coords);

    reset_vecs(0);
}

double PotentialE(double out, double in, int Z) { //calculated in VOLTS NOT MILLIVOLTS
  double E = (R_constant * body_temp) / (F * Z) * log(out / in); // log(x) = ln(x) in cpp
  if ((isinf(E)) || (E != E)) {
    cout << "YOUR E FUNCTION IS FAULTING! Probably, the concentration inside went to 0, or you entered z = 0." << endl;
  }
  //THIS IS USED TO TEST POTENTIAL CALCULATION
  cout << "\n THIS IS E in V: " << E << endl;
  return (E);
}

double Piezo_Channel(double potential, int time, int x, int y){

    int stochastic = stochastic_opening(generator);
    double open_temp;
    double P_open_temp;

    open_temp = vec_num_open[time][x][y] + abs(stochastic);

    if (open_temp >= N_Piezo_channels){
        open_temp = N_Piezo_channels;
    }

    int local_tau = 1.6/delta_T;
    if (!(open_counter % local_tau)){
        open_temp = open_temp*0.9048; //this is the time constant of Piezo, so it is not relevant on a micro s scale 
    }

    vec_num_open[time + 1][x][y] = open_temp;
    vec_num_closed[time + 1][x][y] = 1 - open_temp;
    
    G_Piezo_total = open_temp*G_Piezo_single;
    Piezo_current = (potential*G_Piezo_total)/2; //dividing by 2 because Ca is 2+ charge (affecting current by factor of 2)
    vec_Piezo_current[time + 1][x][y] = Piezo_current;

    return(Piezo_current);
}

double Calcium_concentration(double x){

    ofstream create_file(path2);
    ofstream myfile;
    myfile.open(path2);

    myfile << "test\n";


    myfile.close();

    cout << "high" << endl;

    double x_max = 9;
    double y_max = 9;
    double time_max = 1;
    E_Ca = PotentialE(0.0024, 0.00000012, 2);

    vec_time.push_back(fill_2dvecs(x_max, y_max, 0.00000012));
    vec_num_open.push_back(fill_2dvecs(x_max, y_max, 0));
    vec_num_closed.push_back(fill_2dvecs(x_max, y_max, N_Piezo_channels));
    vec_Piezo_current.push_back(fill_2dvecs(x_max, y_max, 0));

    cout << "break point 2" << endl;

    vec_buff_bound.push_back(0);

    double scaling_factor = 1;

    for(int time_temp = 0; time_temp <= time_max; time_temp++){

    vec_num_open.push_back(fill_2dvecs(x_max, y_max, 0));
    vec_num_closed.push_back(fill_2dvecs(x_max, y_max, 0));
    vec_time.push_back(fill_2dvecs(x_max, y_max, 0));
    vec_Piezo_current.push_back(fill_2dvecs(x_max, y_max, 0));
    
        for(int i = 0; i <= x_max; i++){
            for(int j = 0; j <= y_max; j++){

                //cout << time_temp << "," << i << "," << j << endl;

                C_cyt = vec_time[time_temp][i][j];
                

                Ca_c_dT = delta_T*(scaling_factor*Piezo_Channel(E_Ca, time_temp, i, j));
                vec_time[time_temp + 1][i][j] = vec_time[time_temp][i][j] + Ca_c_dT;

            }
        }
    }

    cout << "low" << endl;
    return(0);
}


double output_file(double x)
{
    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    if (myfile.is_open())
    {
        cout << "File successfully open" << endl;
        myfile << "File successfully open,";
    }
    else
    {
        cout << "Error opening file";
    }

    myfile << "time,x,y\n";

    //cout << "Break point 1" << endl;

    Calcium_concentration(0);

    vector<int> sizes;

    //cout << "Break point 2" << endl;

    sizes.insert(sizes.begin(),vec_time.size());

    //cout << "Break point 3" << endl;

    sort(sizes.begin(), sizes.end());
    int max_size = sizes.back();
    
    cout << max_size << endl;
    bool bool_y; 

    cout << "Break point 4" << endl;

    double x_max = 9;
    double y_max = 9;
    double time_max = 1;

    myfile << "x,y\n";

    for(int time = 0; time <= vec_time.size(); time++){
        for (int x = 0; x <= x_max; x++){
            //cout << "Break point 5" << endl;
            for (int y = 0; y <= y_max; y++){  
                if(y < y_max) {
                    double temp = vec_time[time][x][y]; 
                    myfile << vec_time[time][x][y] << ",";
                    //cout << vec_time[time][x][y] << endl;
                }  
                else {
                    myfile << vec_time[time][x][y];
                    //cout << "high to low" << endl;
                }
            myfile << "\n";
        //cout << "Break point 7" << endl;
            }
        myfile << "\n";
        }
    myfile << "\n";
    }

    myfile.close();

    return (x);
}

int main(void) {
  cout << "Begin" << endl;
  output_file(0);
  cout << "End" << endl;
}
