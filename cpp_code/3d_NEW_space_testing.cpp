#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <vector>
#include <algorithm>
#include <random>

using namespace std;

const char *path1="../data_files/3dtest_output.csv";


vector<double> vec_x;
vector<double> vec_y; 
vector<double> vec_z; 
vector<double> vec_valid; 
vector<double> vec_temp; 
vector<vector<double> > vec_time;

// build a sphere via enterring radius, populate via finding area of sphere and dividing by 1000 or so

int reset_vecs(int x){
    vec_x.clear();
    vec_y.clear();
    vec_z.clear();
    vec_valid.clear();
    vec_time.clear();
    return(0);
}

double access_coordinates(double x, double y, double z){
    double temp;
    int counter_time;

    bool bool_x, bool_y, bool_z;
    for(double i = 0; i <= 9; i ++){
        for(double j = 0; j <= 9; j ++){
            for(double k = 0; k <= 9; k ++){
                bool_x = (vec_time[0][i] == x) ? true : false;
                bool_y = (vec_time[1][j] == y) ? true : false;
                bool_z = (vec_time[2][k] == z) ? true : false;

                if(bool_x && bool_y && bool_z){
                    cout << "found coordinates: " << x << y << z << endl;
                }
            }
        }
    }
    return(1);
}

int display_vectors(int x){
    double temp;
    // for(const vector<double>& row : vec_time){
    //     for(double col : row){
    //         cout << col << " ";
    //     }
    //     cout << endl;
    // }
    temp = access_coordinates(3,3,3);
    cout << temp << endl;
    return(0);
}

int fill_vectors(int temp){
    int x_max = 9;
    int y_max = 9;
    int z_max = 9;
    int time_max = 4;

    //for(int time_temp = 0; time_temp <= time_max; time_temp++){

    for(int i = 0; i <= x_max; i ++){
        for(int j = 0; j <= y_max; j ++){
            for(int k = 0; k <= z_max; k ++){
                vec_x.push_back(i);
                vec_y.push_back(j);
                vec_z.push_back(k);
                if((i >= 3 && i <= 6) && (j >= 3 && j <= 6) && (k >= 3 && k <= 6)){
                    vec_valid.push_back(1);
                }
                else{
                    vec_valid.push_back(0);
                }
            }
        }
    }
        vec_time.push_back(vec_x);
        vec_time.push_back(vec_y);
        vec_time.push_back(vec_z);
        vec_time.push_back(vec_valid);
        vec_x.clear();
        vec_y.clear();
        vec_z.clear();
        vec_valid.clear();
    //}

    cout << "high" << endl;

    display_vectors(0);
    return(0);
}

double output_file(double x)
{
    reset_vecs(0);
    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    cout << "Break point 1" << endl;

    fill_vectors(0);

    vector<int> sizes;

    cout << "Break point 2" << endl;

    sizes.insert(sizes.begin(),vec_time.size());

    cout << "Break point 3" << endl;

    sort(sizes.begin(), sizes.end());
    int max_size = sizes.back();
    
    //cout << max_size << endl;
    bool bool_y; 

    //cout << "Break point 4" << endl;

    myfile << "x,y\n";
    for (int i = 0; i < max_size; i++)
    {
        //cout << "Break point 5" << endl;
        for (int z = 0; z < vec_time[i].size(); z++)
        {  

        //cout << "Break point 6" << endl;

        if(z < vec_time[i].size() - 1) 
            myfile << vec_time[i][z] << ",";
        else
            myfile << vec_time[i][z];

        //cout << "Break point 7" << endl;
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