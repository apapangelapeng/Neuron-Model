#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <vector>
#include <map>
#include<algorithm>
#include"ReduceConst.h"

using namespace std;


const char *path1 = "../data_files/reduced_V_output.csv";

void write_Reduced_AP(vector<vector<double> > & v_map){
    cout << "I HAVE BEEN SUMMONED " << endl;
    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    myfile << "V" << V_start << "\n";
    for (int i = 0; i < v_map[0].size(); i++)
    {
        for (int j = 0; j < v_map.size(); j++)
        {
            string end = (j == v_map.size() - 1) ? "\n" : ",";
            myfile << v_map[j][i] << end;
        }
    }

    myfile << "\n";

}

vector<vector<double> > Reduced_AP(double z)
{
    vector<double> vec_V, vec_N;

vector<double> vec_Vdt, vec_Ndt, vec_x_Vdt, vec_y_Ndt, vec_V_equal_n, vec_N_equal_v;



double V_dt, N_dt;
    int x = 0;

    V_start = 0;

    vec_V.push_back(V_start);
    vec_N.push_back(N_start);
    vector<vector<double> > v_map_2d;

    for (double i = 0; i <= 50; i += delta_t)
    {
        
        if(i >= 5 && i <= 10){
            current_temp = 0.05;
        }
        else{
            current_temp = 0;
        }
        
        cout << current_temp << endl; 


        time_d = vec_V[x] * (vec_V[x] - v_threshold) * (v_max - vec_V[x]);
        V_dt = current_temp + time_d - vec_N[x];

        cout << "V_DT " << V_dt << endl;
        vec_V.push_back(vec_V[x] + V_dt);
        N_dt = e * (vec_V[x] - (gam * vec_N[x]));
        vec_N.push_back(vec_N[x] + N_dt);

        x += 1;
        }

    vector<double> vec_V_deep = vec_V;
    v_map_2d.push_back(vec_V_deep); 
   
    cout << "V_start: " << V_start << endl;

   
   return v_map_2d;
}


int main(void)
{
    cout << "Reduced Action potential Running. Outputfile location: " << endl;
    cout<<path1<<endl;
    vector<vector<double> >  v_map_2d = Reduced_AP(0);
    write_Reduced_AP(v_map_2d);

    cout << "End" << endl;
}
