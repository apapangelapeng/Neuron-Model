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


vector<double> vec_xyz; 
vector<vector<double> > vec_time;

// build a sphere via enterring radius, populate via finding area of sphere and dividing by 1000 or so

int reset_vecs(int x){
    vec_xyz.clear();
    vec_time.clear();
    return(0);
}

int display_vectors(int x){
    for(const vector<double>& row : vec_time){
        for(double col : row){
            cout << col << " ";
        }
        cout << endl;
    }
    return(0);
}

int fill_vectors(int temp){
    int x_max = 10;
    int y_max = 10;
    int z_max = 10;
    int time_max = 4;

    for(int time_temp = 0; time_temp <= time_max; time_temp++){

    for(int i = 0; i <= x_max; i ++){
        for(int j = 0; j <= y_max; j ++){
            for(int k = 0; k <= z_max; k ++){
                vec_xyz.push_back(i);
                vec_xyz.push_back(j);
                vec_xyz.push_back(k);
            }
        }
    }
        vec_time.push_back(vec_xyz);
        vec_xyz.clear();
    }

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