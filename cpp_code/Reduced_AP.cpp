#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <vector>
#include <map>
using namespace std;

/*
THIS FILE IS JUST FOR ME TO TEST THINGS
*/

const char *path1 = "../data_files/reduced_V_output.csv";
const char *path2 = "../data_files/reduced_null_output.csv";
const char *path3 = "../data_files/reduced_diffusion_output.csv";

vector<double> vec_V, vec_N, vec_Space,vec_tiny_N,vec_N_Zerodt,vec_V_Zerodt;

vector<double> vec_V_Zerodt, vec_N_Zerodt;

    vector<vector<double>> v_map_2d;

double V_dt, N_dt;

double v_threshold = 0.2;
double v_max = 1;
double V_start = 0;
double v_range = 0.5;

double N_start = 0;

double gam = 3;
double gam_range = 0;

double e = 0.005;

double current = 0;
double current_range = 0.5;
double current_temp;

double diffusion = 0.5;
double diffusion_range = 2;

double spatial_d = 0;
double delta_x = 0.1;
double x_range = 1;

double time_d = 0;
double delta_t = 0.1;

int reset_vecs(int x)
{
    vec_N.clear();
    vec_V.clear();
    return (0);
}

double nullcline_generation_a(double x)
{
    int y;
    reset_vecs(0);
    double local_v_start = -0.3;
    double local_n_start = -0.05;
    double tiny_n_start = -0.05;
    vec_N.push_back(local_n_start);
    vec_V.push_back(local_v_start);
    vec_tiny_N.push_back(local_n_start);

    y = 0;
    for (double v_temp = local_v_start; v_temp <= 1.1; v_temp += 0.0001)
    {
        // cout << "v temp : " << v_temp << endl;
        N_dt = e * (v_temp - (gam * vec_N[y]));
        vec_N.push_back(vec_N[y] + N_dt);
        y++;
    }

    y = 0;
    for (double n_temp = local_n_start; n_temp <= 0.15; n_temp += 0.0001)
    {
        time_d = vec_V[y] * (vec_V[y] - v_threshold) * (v_max - vec_V[y]);
        V_dt = current + time_d - n_temp;
        vec_V.push_back(vec_V[y] + V_dt);
        y++;
    }
    //vec_tiny_N
    y = 0;
    for (double v_temp = local_v_start; v_temp <= 0.5; v_temp += 0.000001)
    {
        // cout << "v temp : " << v_temp << endl;
        N_dt = e * (v_temp - (gam * vec_tiny_N[y]));
        vec_tiny_N.push_back(vec_tiny_N[y] + N_dt);
        y++;
    }


    return (0);
}

double nullcline_generation(double x)
{
    int y;
    reset_vecs(0);
    double local_v_start = -0.3;
    double local_n_start = -0.05;
    double tiny_n_start = -0.05;
    vec_N.push_back(local_n_start);
    vec_V.push_back(local_v_start);
    vec_tiny_N.push_back(local_n_start);

    y = 0;
    for (double v_temp = local_v_start; v_temp <= 1.1; v_temp += 0.0001)
    {
        // cout << "v temp : " << v_temp << endl;
        N_dt = e * (v_temp - (gam * vec_N[y]));
        vec_N.push_back(vec_N[y] + N_dt);
        y++;
    }

    y = 0;
    for (double n_temp = local_n_start; n_temp <= 0.15; n_temp += 0.0001)
    {
        time_d = vec_V[y] * (vec_V[y] - v_threshold) * (v_max - vec_V[y]);
        V_dt = current + time_d - n_temp;
        vec_V.push_back(vec_V[y] + V_dt);
        y++;
    }
    for (double v_temp = local_v_start; v_temp <= 1.1; v_temp += 0.01)
    {
        for (double n_temp = local_n_start; n_temp <= 0.35; n_temp += 0.005)
        {
            time_d = v_temp * (v_temp - v_threshold) * (v_max - v_temp);
            V_dt = current + time_d - n_temp;
            if(V_dt <= 0.01 && V_dt >= -0.01){ 
                vec_N_Zerodt.push_back(n_temp);
                vec_V_Zerodt.push_back(v_temp);
            }
        }
    }
    for (double v_temp = local_v_start; v_temp <= 1.1; v_temp += 0.01)
    {
        for (double n_temp = local_n_start; n_temp <= 0.15; n_temp += 0.005)
        {
            N_dt = (v_temp - (gam * n_temp));
            if(N_dt <= 0.01 && N_dt >= -0.01){ 
                vec_N_Zerodt.push_back(n_temp);
                vec_V_Zerodt.push_back(v_temp);
            }
        }
    }
    


    return (0);
}


double output__nullcline_file_a(double x)
{
    ofstream create_file(path2);
    ofstream myfile;
    myfile.open(path2);

    nullcline_generation_a(0);

   
    vector<int> sizes;
    /*cout<<vec_V.size()<<endl;
    cout<<vec_N.size()<<endl;
    cout<<vec_tiny_N.size()<<endl;*/
    sizes.insert(sizes.begin(),vec_V.size());
    sizes.insert(sizes.begin(),vec_N.size());
    sizes.insert(sizes.begin(),vec_tiny_N.size());
    sort(sizes.begin(), sizes.end());
    int max_size =sizes.back();
    
    cout <<max_size << endl;
    bool N ;
    bool V;
    bool tiny_N;
    myfile << "N,V,tiny_N\n";
    for (int i = 0; i < max_size; i++)
    {
        N = (vec_N.size()> i) ? true  : false;
        V = (vec_V.size() > i) ? true : false;
        tiny_N = (vec_tiny_N.size()> i) ? true  : false;
        if(N) myfile << vec_N[i] << "," ;
        if(!N) myfile << "," ;
        if(V) myfile << vec_V[i]<<"," ;
        if(!V) myfile <<"," ;
        if(tiny_N) myfile << vec_tiny_N[i];
        myfile<<"\n";
    }
    myfile.close();
    return (x);
}


double output__nullcline_file(double x)
{
    ofstream create_file(path2);
    ofstream myfile;
    myfile.open(path2);

    nullcline_generation(0);

   
    vector<int> sizes;
    /*cout<<vec_V.size()<<endl;
    cout<<vec_N.size()<<endl;
    cout<<vec_tiny_N.size()<<endl;*/
    sizes.insert(sizes.begin(),vec_V.size());
    sizes.insert(sizes.begin(),vec_N.size());
    sizes.insert(sizes.begin(),vec_N_Zerodt.size());
    sizes.insert(sizes.begin(),vec_V_Zerodt.size());
    sort(sizes.begin(), sizes.end());
    int max_size =sizes.back();
    
    cout <<max_size << endl;
    bool N ;
    bool V;
    bool N_0;
    bool V_0;
    myfile << "N,V,N_0,V_0\n";
    for (int i = 0; i < max_size; i++)
    {
        N = (vec_N.size()> i) ? true  : false;
        V = (vec_V.size() > i) ? true : false;
        N_0 = (vec_N_Zerodt.size()> i) ? true  : false;
        V_0 = (vec_V_Zerodt.size()> i) ? true  : false;
        if(N) myfile << vec_N[i] << "," ;
        if(!N) myfile << "," ;
        if(V) myfile << vec_V[i]<<"," ;
        if(!V) myfile <<"," ;
        if(N_0) myfile << vec_N_Zerodt[i]<<",";
        if(!N_0) myfile <<",";
        if(V_0) myfile << vec_V_Zerodt[i];
        //if(!V_0) myfile << vec_tiny_N[i];
        myfile<<"\n";
    }
    myfile.close();
    return (x);
}

void write_Reduced_AP(vector<vector<double>> & v_map){
    cout << "I HAVE BEEN SUMMONED " << endl;
    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    /*
    for (V_start = -0.5; V_start <= v_range; V_start += 0.3)
    {

        string end = (V_start >= v_range - 0.3) ? "\n" : ",";
        myfile << "V" << V_start << end;
    }
    */

    myfile << "V" << V_start << "\n";
    for (int i = 0; i < vec_V.size(); i++)
    {
        for (int j = 0; j < v_map.size(); j++)
        {
            string end = (j == v_map.size() - 1) ? "\n" : ",";
            myfile << v_map[j][i] << end;
        }

        // cout << vec_V[i] << endl;
    }

    myfile << "\n";

}

double Reduced_AP(double z)
{

    for (V_start = -0.5; V_start <= v_range; V_start += 0.3)
    {
        double current = 0;
        reset_vecs(0);
        int x = 0;
        vec_V.push_back(V_start);
        vec_N.push_back(N_start);

        for (double i = 0; i <= 50; i += delta_t)
        {
            /*
                if(i <= 10){
                    current_temp = current;
                }
                else{
                    current_temp = 0;
                }
            */

            // cout << "Break Point 4" << endl;

            time_d = vec_V[x] * (vec_V[x] - v_threshold) * (v_max - vec_V[x]);
            V_dt = current + time_d - vec_N[x];
            vec_V.push_back(vec_V[x] + V_dt);
            N_dt = e * (vec_V[x] - (gam * vec_N[x]));
            vec_N.push_back(vec_N[x] + N_dt);

            // cout << "Break Point 5" << endl;

            // cout << "x = " << x << endl;

            x += 1;
        }

        vector<double> vec_V_deep = vec_V;
        v_map_2d.push_back(vec_V_deep); 

        cout << "V_start: " << V_start << endl;
    }
    write_Reduced_AP(v_map_2d);
    return (0);
}

double Diffusion_AP(double z)
{
    ofstream create_file(path3);
    ofstream myfile;
    myfile.open(path3);

    // cout << "Break Point 1" << endl;

    double V_start = 0.4;
    double N_start = 0;
    double Space_start = 0;

    int counter = 0;

    // cout << "Break Point 3" << endl;

    vec_V.push_back(V_start);
    vec_N.push_back(N_start);
    vec_Space.push_back(Space_start);

    for (double i = 0; i <= 100; i += delta_t)
    {
        vec_Space.clear();

        /*
        if (i <= 2)
        {
            current_temp = current;
        }
        else
        {
            current_temp = 0;
        }
        */

        time_d = vec_V[counter] * (vec_V[counter] - v_threshold) * (v_max - vec_V[counter]);
        V_dt = current + time_d - vec_N[counter];
        vec_V.push_back(vec_V[counter] + V_dt);
        N_dt = e * (vec_V[counter] - (gam * vec_N[counter]));
        vec_N.push_back(vec_N[counter] + N_dt);

        /*
        cout << "volt " << vec_V[counter] << endl;
        cout << "time d = " << time_d << endl;
        cout << "v_dt = " << V_dt << endl;
        cout << "n " << vec_N[counter] << endl;
        cout << "counter " << counter << endl;
        */

        counter += 1;

        // cout << "Break Point 4" << endl;

        // cout << i << ", and ROUNDED: " << floor(i) << endl;
        if (counter % 10 == 0 || counter == 1)
        {
            //cout << "PASSED counter: " << counter << endl;
            vec_Space.push_back(vec_V[counter]);
            vec_Space.push_back(vec_V[counter]);
            vec_Space.push_back(vec_V[counter]);
            int local_counter = 2;
            for (double space = 0; space <= x_range; space += delta_x)
            {
                // cout << "Voltage at this time: " << vec_V[counter] << endl;
                // cout << "Break Point 5" << endl;

                spatial_d = (vec_Space[local_counter] - (2 * vec_Space[local_counter - 1]) + vec_Space[local_counter - 2]);
                vec_Space.push_back(diffusion * vec_Space[local_counter] - pow(delta_t, 2) * spatial_d);
                local_counter += 1;
                //cout << spatial_d << endl;
                //cout << "counter = " << counter << endl; 
            }
            myfile << "T" << i << ",";
            for (int i = 0; i < vec_Space.size(); i++)
            {
                myfile << vec_Space[i] << ",";
                // cout << "Break Point 6" << endl;
            }
            myfile << "\n";
        }

        // cout << "Break point 6" << endl;

        /*
        myfile << "V" << diffusion << ",";
        for(int i = 0; i < vec_V.size(); i++){
            myfile << vec_V[i] << ",";
            //cout << vec_V[i] << endl;
        }
        */

        // cout << "Time: " << i << endl;
    }
    /*
    myfile << "Voltage,";
    for(int i = 0; i < vec_V.size(); i++){
        myfile << vec_V[i] << ",";
        //cout << vec_N[i] << endl;
        //cout << "Break Point 6" << endl;
    }
    myfile << "\n";

    myfile << "N,";
    for(int i = 0; i < vec_N.size(); i++){
        myfile << vec_N[i] << ",";
        //cout << vec_N[i] << endl;
        //cout << "Break Point 6" << endl;
    }
    */

    myfile.close();

    vector<double> vec_V_deep = vec_V;
    v_map_2d.push_back(vec_V_deep); 
    write_Reduced_AP(v_map_2d);

    cout << "V_start: " << V_start << endl;
    return (0);
}

int main(void)
{
    cout << "Begin" << endl;

    //Reduced_AP(0);
    //output__nullcline_file(0);
    Diffusion_AP(0);

    cout << "End" << endl;
}
