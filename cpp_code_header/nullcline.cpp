        #include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <vector>
#include <map>
#include<algorithm>
#include "ReduceConst.h"

using namespace std;

/*
THIS FILE IS JUST FOR ME TO TEST THINGS
*/

const char *path1 = "../data_files/nullcline_metadata.csv";
const char *path2 = "../data_files/nullcline_output.csv";
const char *path3 = "../data_files/reduced_diffusion_output.csv";

vector<double> vec_V, vec_N, vec_Space,vec_tiny_N,vec_N_Zerodt,vec_V_Zerodt;

vector<double> vec_V_Zero_dvdt, vec_N_Zero_dvdt, vec_V_Zero_dndt, vec_N_Zero_dndt;

vector<double> vec_Vdt, vec_Ndt, vec_x_Vdt, vec_y_Ndt, vec_V_equal_n, vec_N_equal_v;

vector<vector<double> > v_map_2d;

double V_dt, N_dt;



int reset_vecs(int x)
{
    vec_N.clear();
    vec_V.clear();
    return (0);
}

void output_metadata(double local_v_start,double local_v_end,double v_step,double local_n_start,double local_n_end,double n_step)
{
    ofstream create_file(path2);
    ofstream myfile;
    myfile.open(path2);
    myfile<<"local_v_start,local_v_end,v_step,local_n_start,local_n_end,n_step \n";
    myfile<<local_v_start<<",";
    myfile<<local_v_end<<",";
    myfile<<v_step<<",";
    myfile<<local_n_start<<",";
    myfile<<local_n_end<<",";
    myfile<<n_step;
    myfile.close();


}
double nullcline_generation(double x)
{
    int y;
    reset_vecs(0);
    
    double local_v_start = -0.15;
    double local_v_end	=0.45;
    double v_step	=0.0001;
   double local_n_start = -0.05;
    double local_n_end	= 0.12;
    double n_step=0.0001;
    output_metadata( local_v_start, local_v_end, v_step, local_n_start, local_n_end, n_step);
    vec_N.push_back(local_n_start);
    vec_V.push_back(local_v_start);


    double current = 0.0; 

    y = 0;
    for (double v_temp = local_v_start; v_temp <= local_v_end; v_temp += v_step)
    {
        // cout << "v temp : " << v_temp << endl;
        N_dt = e * (v_temp - (gam * vec_N[y]));
        vec_N.push_back(vec_N[y] + N_dt);
        y++;
    }

    y = 0;
    for (double n_temp = local_n_start; n_temp <= local_n_end; n_temp += n_step)
    {
        time_d = vec_V[y] * (vec_V[y] - v_threshold) * (v_max - vec_V[y]);
        V_dt = current + time_d - n_temp;
        vec_V.push_back(vec_V[y] + V_dt);
        y++;
    }
    for (double v_temp = -0.22; v_temp <= 1.1; v_temp += 0.01)
    {
        for (double n_temp = local_n_start; n_temp <= 0.20; n_temp += 0.005)
        {
            time_d = v_temp * (v_temp - v_threshold) * (v_max - v_temp);
            V_dt = current + time_d - n_temp;
            if(V_dt <= 0.01 && V_dt >= -0.01){ 
                vec_N_Zero_dvdt.push_back(n_temp);
                vec_V_Zero_dvdt.push_back(v_temp);
            }
        }
    }
    for (double v_temp = local_v_start; v_temp <= 1.1; v_temp += 0.01)
    {
        for (double n_temp = local_n_start; n_temp <= 0.15; n_temp += 0.005)
        {
            N_dt = (v_temp - (gam * n_temp));
            if(N_dt <= 0.01 && N_dt >= -0.01){ 
                vec_N_Zero_dndt.push_back(n_temp);
                vec_V_Zero_dndt.push_back(v_temp);
                //cout << "V " << v_temp << endl;
                //cout << "N " << n_temp << endl;
            }
        }
    }
    
    for (double v_temp = -0.2; v_temp <= 1.1; v_temp += 0.05)
    {
        for (double n_temp = local_n_start; n_temp <= 0.15; n_temp += 0.01)
        {
            N_dt = e*(v_temp - (gam * n_temp));
            time_d = v_temp * (v_temp - v_threshold) * (v_max - v_temp);
            V_dt = current + time_d - n_temp;
            vec_Ndt.push_back(N_dt);
            vec_Vdt.push_back(V_dt);
            vec_y_Ndt.push_back(n_temp);
            vec_x_Vdt.push_back(v_temp);
            if((N_dt <= (V_dt+0.01)) && (N_dt >= (V_dt-0.01))){
                vec_V_equal_n.push_back(n_temp);
                vec_N_equal_v.push_back(v_temp);
            }
        }
    }
    return (0);
}



double output__nullcline_file(double x)
{
    ofstream create_file(path2);
    ofstream myfile;
    myfile.open(path2);

   

    vector<int> sizes;
    /*cout<<vec_V.size()<<endl;
    cout<<vec_N.size()<<endl;
    cout<<vec_tiny_N.size()<<endl;*/
    sizes.insert(sizes.begin(),vec_V.size());
    sizes.insert(sizes.begin(),vec_N.size());
    sizes.insert(sizes.begin(),vec_N_Zero_dvdt.size());
    sizes.insert(sizes.begin(),vec_V_Zero_dvdt.size());
    sizes.insert(sizes.begin(),vec_N_Zero_dndt.size());
    sizes.insert(sizes.begin(),vec_V_Zero_dndt.size());
    sizes.insert(sizes.begin(),vec_Ndt.size());
    sizes.insert(sizes.begin(),vec_Vdt.size());
    sizes.insert(sizes.begin(),vec_y_Ndt.size());
    sizes.insert(sizes.begin(),vec_x_Vdt.size());
    sizes.insert(sizes.begin(),vec_V_equal_n.size());
    sizes.insert(sizes.begin(),vec_N_equal_v.size());
    sort(sizes.begin(), sizes.end());
    int max_size = sizes.back();
    
    cout <<max_size << endl;
    bool N ;
    bool V;
    bool dv_N_0;
    bool dv_V_0;
    bool dn_N_0;
    bool dn_V_0;
    bool Ndt;
    bool Vdt;
    bool N_y_dt;
    bool V_x_dt;
    bool Vequal;
    bool Nequal;
    
    myfile << "N,V,N_dv_0,V_dv_0,N_dn_0,V_dn_0,Ndt,Vdt,N_y_dt,V_x_dt,Nequal,Vequal\n";
    for (int i = 0; i < max_size; i++)
    {
        N = (vec_N.size()> i) ? true  : false;
        V = (vec_V.size() > i) ? true : false;
        dv_N_0 = (vec_N_Zero_dvdt.size()> i) ? true  : false;
        dv_V_0 = (vec_V_Zero_dvdt.size()> i) ? true  : false;
        dn_N_0 = (vec_N_Zero_dndt.size()> i) ? true  : false;
        dn_V_0 = (vec_V_Zero_dndt.size()> i) ? true  : false;
        Ndt = (vec_Ndt.size()> i) ? true  : false;
        Vdt = (vec_Vdt.size()> i) ? true  : false;
        N_y_dt = (vec_Ndt.size()> i) ? true  : false;
        V_x_dt = (vec_Vdt.size()> i) ? true  : false;
        Nequal = (vec_V_equal_n.size()> i) ? true  : false;
        Vequal = (vec_N_equal_v.size()> i) ? true  : false;
        if(N) myfile << vec_N[i] << "," ;
        if(!N) myfile << "," ;
        if(V) myfile << vec_V[i]<<"," ;
        if(!V) myfile <<"," ;
        if(dv_N_0) myfile << vec_N_Zero_dvdt[i] << ",";
        if(!dv_N_0) myfile <<",";
        if(dv_V_0) myfile << vec_V_Zero_dvdt[i] << ",";
        if(!dv_V_0) myfile <<",";
        if(dn_N_0) myfile << vec_N_Zero_dndt[i] << ",";
        if(!dn_N_0) myfile <<",";
        if(dn_V_0) myfile << vec_V_Zero_dndt[i] << ",";
        if(!dn_V_0) myfile <<",";
        if(Ndt) myfile << vec_Ndt[i] << ",";
        if(!Ndt) myfile <<",";
        if(Vdt) myfile << vec_Vdt[i] << ",";
        if(!Vdt) myfile <<",";
        if(N_y_dt) myfile << vec_y_Ndt[i] << ",";
        if(!N_y_dt) myfile <<",";
        if(V_x_dt) myfile << vec_x_Vdt[i] << ",";
        if(!N_y_dt) myfile <<",";
        if(Nequal) myfile << vec_N_equal_v[i] << ",";
        if(!Nequal) myfile <<",";
        if(Vequal) myfile << vec_V_equal_n[i];
        //if(!dn_V_0) myfile << vec_tiny_N[i];
        myfile<<"\n";

        //cout << "break";
    }
    myfile.close();
    return (x);
}


int main(void)
{
    cout << "Begin" << endl;

    //Reduced_AP(0);
     nullcline_generation(0);
    output__nullcline_file(0);
    //Diffusion_AP(0);

    cout << "End" << endl;
}