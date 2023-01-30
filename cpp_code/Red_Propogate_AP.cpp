#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

/*
THIS FILE IS JUST FOR ME TO TEST THINGS
*/

const char *path1 = "../data_files/reduced2dV_output.csv";

vector<double> vec_V, vec_N, vec_Space,vec_tiny_N,vec_N_Zerodt,vec_V_Zerodt;

vector<double> vec_V_Zero_dvdt, vec_N_Zero_dvdt, vec_V_Zero_dndt, vec_N_Zero_dndt;

vector<double> vec_Vdt, vec_Ndt, vec_x_Vdt, vec_y_Ndt;

vector<vector<double> > v_map_2d;
vector<vector<double> > voltage_2d;
vector<vector<double> > space_2d;
vector<vector<double> > n_2d;

double V_dt, N_dt;

double v_threshold = 0.2;
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
double current_input;

double diffusion = 0.2;
double diffusion_range = 2;

double spatial_d = 0;
double delta_x = 0.01;
double x_range = 1; 

double time_d = 0;
double delta_t = 0.001;

int reset_vecs(int x)
{
    vec_N.clear();
    vec_V.clear();
    //vec_Space.clear();
    return (0);
}

void output_file(vector<vector<double> > v_map){
    cout << "I HAVE BEEN SUMMONED " << endl;
    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    int col_num = x_range/delta_x;

    //myfile << "V" << V_start << "\n";
    for (int j = 0.19/delta_t; j < 0.4/delta_t; j+=1)
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

double Diffusion_AP(double z)
{
    //cout << "Break Point 1" << endl;

    double V_start = 0.4;
    double N_start = 0;
    //double Space_start = 0.05;
    double current = 0; 

    int counter_space = 0;
    int counter_time = 0;

    for (double space = 0; space <= x_range+delta_x; space += delta_x) //THIS COUNTS SPACE!!
    { 
        vec_V.push_back(V_start);
        vec_N.push_back(N_start);
    }
    
    voltage_2d.push_back(vec_V);
    n_2d.push_back(vec_N);
    //space_2d.push_back(vec_Space);

    //cout << "Break Point 2" << endl;

    for (double i = 0; i <= 0.5; i += delta_t) //THIS COUNTS TIME!!
    {
        counter_space = 0;

        reset_vecs(0);

        for (double space = 0; space <= x_range+delta_x; space += delta_x) //THIS COUNTS SPACE!!
        { 

        if(space <= 0.1){
            if((i <= 0.21 && i >= 0.2)){
                current = 0.05;
            }
            else{
                current = 0;
            }
        }
        else{
            current = 0;
        }
        

        // if((i <= 0.1) || (i <= 0.6 && i >= 0.5)){
        //         cout << counter_space << " current = " << current << endl;
        //     }

        if((counter_space >= 1) && (counter_space <= 99)){
            spatial_d = diffusion*(voltage_2d[counter_time][counter_space + 1] - (2 * voltage_2d[counter_time][counter_space]) + voltage_2d[counter_time][counter_space - 1]);
            // cout << "pos1 " << voltage_2d[counter_time][counter_space] << endl;
            // cout << "pos2 " << 2 * voltage_2d[counter_time][counter_space - 1] << endl;
            // cout << "pos3 " << voltage_2d[counter_time][counter_space - 2] << endl;
            //cout << diffusion*vec_Space[counter_space - 1] << endl;
            //cout << current << endl;
            current_input = spatial_d;
        }
        else{
            current_input = 0;
        }

        time_d = voltage_2d[counter_time][counter_space] * (voltage_2d[counter_time][counter_space] - v_threshold) * (v_max - voltage_2d[counter_time][counter_space]);
        //cout << "POSITION " << counter_time << "," << counter_space << " = " << voltage_2d[counter_time][counter_space] << endl; 
        //cout << "time_d " << time_d << endl;
        V_dt = current + current_input + time_d - n_2d[counter_time][counter_space];

        // if((i <= 0.1) || (i <= 0.6 && i >= 0.5)){
        //         cout << counter_space << " " << current << endl;
        //     }

        vec_V.push_back(voltage_2d[counter_time][counter_space] + V_dt);
        //cout << "V_dt " << V_dt << endl;

        N_dt = e * (voltage_2d[counter_time][counter_space] - (gam * n_2d[counter_time][counter_space]));

        vec_N.push_back(n_2d[counter_time][counter_space] + N_dt);

        // if(counter_time <= 100){
        // cout << "N_dt " << N_dt << endl;
        // cout << "POSITION " << counter_time << "," << counter_space << " = " << n_2d[counter_time][counter_space] << endl; 
        // cout << "POSITION " << counter_space << " = " << vec_N[counter_space] << endl;
        // }

        //cout << vec_Space[2] << endl;

        // cout << counter_space << " position " << vec_Space[counter_space] << endl;

        //cout << counter_space << " spatial_d " << spatial_d << endl;

        counter_space++;

        }

        voltage_2d.push_back(vec_V);
        n_2d.push_back(vec_N);

    counter_time += 1;

    
    for (int i = 0; i < vec_V.size(); i++)
    {
        //cout << "2d " << voltage_2d[counter_time][i] << endl;
        //cout << "1d " << vec_V[i] << endl;
    }
        //cout << "Break Point 3 and: " << counter_time << endl
    }
    output_file(voltage_2d);

    return (0);
}


int main(void)
{
    cout << "Begin" << endl;

    //Reduced_AP(0);
    Diffusion_AP(0);

    cout << "End" << endl;
}
