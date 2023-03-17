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

const char *path1="../data_files/Piezo_Channel.csv";

default_random_engine generator;
normal_distribution<double> stochastic_opening(0,4);

double PotentialE(double out, double in, int Z) { //calculated in VOLTS NOT MILLIVOLTS
  double E = (R_constant * body_temp) / (F * Z) * log(out / in); // log(x) = ln(x) in cpp
  if ((isinf(E)) || (E != E)) {
    cout << "YOUR E FUNCTION IS FAULTING! Probably, the concentration inside went to 0, or you entered z = 0." << endl;
  }
  //THIS IS USED TO TEST POTENTIAL CALCULATION
  cout << "\n THIS IS E in V: " << E << endl;
  return (E);
}

double Compute_J_on(double C_cyt){
  double scaling_factor = 100000; 
  double local_C_cyt = C_cyt*scaling_factor;
  //cout << "Compute_J_on active" << endl;
  k_buff_bind = 0.600; //for the buffer BAPTA in mM
  k_buff_unbind = 0.100;

  J_on = k_buff_bind*local_C_cyt*buff_unbound;
  J_off = k_buff_unbind*buff_bound; 
  buff_c_dT = k_buff_bind*local_C_cyt*buff_unbound - k_buff_unbind*buff_bound;
  vec_buff_bound.push_back(vec_buff_bound[buff_counter] + buff_c_dT);

  buff_diff = J_off - J_on;
  buff_diff = buff_diff*pow(scaling_factor,-1);
  vec_J_on.push_back(buff_diff);

  buff_counter++;
  return(buff_diff*0.000001); //*0.000001
}

double Compute_J_serca(double serc_local){
  //cout << "Compute_J_serca active" << endl;
  double local_C_cyt = 1000*serc_local;

  J_serca = v_serca*(pow(local_C_cyt,2)/(pow(local_C_cyt,2) + pow(K_p,2)));
  if(J_serca != J_serca){
    J_serca = 0;
  }
  vec_J_serca.push_back(J_serca*0.01); //this and ryr are scaled weirdly, I don't know why this works better - otherwise Piezo will dominate
  return(J_serca*0.01);
}

double Compute_J_ryr(double ryr_local){ // I am almost certain that there is something wrong with the kinetic equations that go beyond the paper
  //cout << "Compute_J_ryr active1" << endl;

  double local_C_cyt = 1000*ryr_local; // this is here in case we want to scale 

  w_inf = ((K_a/pow(local_C_cyt,4)) + 1 + (pow(local_C_cyt,3)/K_b))/((1/K_c) + (K_a/pow(local_C_cyt,4)) + 1 + (pow(local_C_cyt,3)/K_b)); 
  //cout << (K_a/pow(local_C_cyt,4)) << endl;

  //cout << "Compute_J_ryr active2" << endl;

  tau_w = w_inf/K_d;
  w_dt = (w_inf - vec_w[w_counter])/tau_w;

  //cout << "Compute_J_ryr active3" << endl;

  vec_w.push_back(vec_w[w_counter] + w_dt);
  P_open = (vec_w[w_counter + 1]*((1 + pow(local_C_cyt,3))/K_b))/((K_a/pow(local_C_cyt,4)) + 1 + (pow(local_C_cyt,3)/K_b));
  J_ryr = (v_rel*P_open + v_leak)*(C_er - local_C_cyt); //I do not like the C_er - local_c_cyt

  if(J_ryr != J_ryr){
    J_ryr = 0;
  }

  w_counter++;

  vec_J_ryr.push_back(J_ryr*0.01);
  return(J_ryr*0.01);
}

double Piezo_Channel(double potential){
  //cout << "Piezo_Channel active" << endl;
  int stochastic = stochastic_opening(generator);
  int open_temp;
  double P_open_temp;

  open_temp = vec_num_open[open_counter] + abs(stochastic);

  // cout << P_open_temp << endl;
  // cout << exp((delta_T*0.001)/T_Piezo) << endl;

  if (open_temp >= N_Piezo_channels){
    open_temp = N_Piezo_channels;
  }

  int local_tau = 1.6/delta_T;
  //cout << local_tau << endl;
  //cout << open_counter % local_tau << endl;
  if (!(open_counter % local_tau)){

    //cout << open_counter << endl;
    open_temp = open_temp*0.9048; //this is the time constant of Piezo, so it is not relevant on a micro s scale 
  }

   if ((open_counter >= 50000) && (open_counter <= 50010)){
     open_temp = N_Piezo_channels; //this is the time constant of Piezo, so it is not relevant on a micro s scale 
   }


  //cout << open_temp << endl;

  vec_num_open.push_back(open_temp);  
  vec_num_closed.push_back(1 - open_temp);
    
  G_Piezo_total = open_temp*G_Piezo_single;

  Piezo_current = (potential*G_Piezo_total)/2; //dividing by 2 because Ca is 2+ charge (affecting current by factor of 2)
  vec_Piezo_current.push_back(Piezo_current);

  open_counter++;
  return(Piezo_current);
}

double Calcium_concentration(int x){


    E_Ca = PotentialE(0.0024, 0.00000012, 2);
    vec_w.push_back(0);
    vec_num_open.push_back(0);
    vec_num_closed.push_back(N_Piezo_channels);
    vec_buff_bound.push_back(0);
    vec_Ca_conc.push_back(0.00000012); //supposed to be 120nM
    // cone_circumference = 2*M_PI*cone_radius;
    // cone_cross_area = M_PI*pow(cone_radius,2);

    // Ca_c_dT = D_diff_Ca*Ca_c_dT_dT + J_ipr + (cone_circumference/cone_cross_area)*(J_in - J_pm) + Compute_J_ryr(C_cyt) - Compute_J_serca(C_cyt) + Compute_J_on(C_cyt);

    double scaling_factor = 1;

    double x_max = 9;
    double y_max = 9;
    double time_max = 1;

    for(int time_temp = 0; time_temp <= time_max; time_temp++){
        for(double i = 0; i <= x_max; i++){
            for(double j = 0; j <= y_max; j++){
                vec_x.push_back(i);
                vec_y.push_back(j);
                C_cyt = vec_time[time_temp][i][j];

                //cout << C_cyt << endl;

                Ca_c_dT = delta_T*(scaling_factor*Piezo_Channel(E_Ca) + scaling_factor*Compute_J_ryr(C_cyt) - scaling_factor*Compute_J_serca(C_cyt) + Compute_J_on(C_cyt));
                vec_time[time_temp + 1][i][j] = vec_time[time_temp][i][j] + Ca_c_dT;
            }
            vec_coords.push_back(vec_x);
            vec_coords.push_back(vec_y);
            vec_x.clear();
            vec_y.clear();

        }
        vec_time.push_back(vec_coords);
        vec_coords.clear();
    }
  return(0);
}


double output_file(double x)
{

    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

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
for(int time = 0; time < vec_time.size(); time++)
    {
    for (int x = 0; x < vec_time[1].size(); x++)
    {
        //cout << "Break point 5" << endl;
        for (int y = 0; y <= y_max; y++)
        {  
            if(y < x_max) 
                myfile << vec_time[time][x][y] << ",";
            else
                myfile << vec_time[time][x][y];

        //cout << "Break point 7" << endl;
        }
        myfile << "\n";
    }
    }
    myfile.close();
    return (x);
}

int main(void) {
  cout << "Begin" << endl;
  output_file(0);
  cout << "End" << endl;
}
