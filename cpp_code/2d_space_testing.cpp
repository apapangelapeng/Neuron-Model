#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <vector>
#include <algorithm>
#include <random>

using namespace std;

const char *path1="../data_files/2dtest_output.csv";


vector<double> vec_x; 
vector<vector<double> > vec_y;

// build a sphere via enterring radius, populate via finding area of sphere and dividing by 1000 or so

int reset_vecs(int x){
    vec_x.clear();
    vec_y.clear();
    return(0);
}

int display_vectors(int x){
    for(const vector<double>& row : vec_y){
        for(double col : row){
            cout << col << " ";
        }
        cout << endl;
    }
    return(0);
}

int fill_vectors(int temp){
    int start_i = 0;
    int start_z = 0;
    int z_max = 100;
    int i_max = 100;
    for (int i = start_i; i <= i_max; i ++){
        for(int z = 0; z <= z_max; z ++){
            if(i == start_i){
                vec_x.push_back(0);
            }
            else if(start_z == z){
                vec_x.push_back(0);
            }
            else if(i_max == i){
                vec_x.push_back(0);
            }
            else if(z_max == z){
                vec_x.push_back(0);
            }
            else 
                vec_x.push_back(1);
        }
        vec_y.push_back(vec_x);
        vec_x.clear();
    }

    cout << "high" << endl;

    for (int i = start_i + 1; i <= i_max - 1; i++)
    {
        for (int z = start_z + 1; z <= z_max - 1; z++)
        { 
            if(vec_y[i - 1][z] == 0)
                vec_y[i][z] = 2;
            else if(vec_y[i + 1][z] == 0)
                vec_y[i][z] = 2;
            else if(vec_y[i][z - 1] == 0)
                vec_y[i][z] = 2;
            else if(vec_y[i][z + 1] == 0)
                vec_y[i][z] = 2;
            else 
                vec_y[i][z] = 3;
        }
    cout << "low" << endl;
    }

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

    sizes.insert(sizes.begin(),vec_y.size());

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
        for (int z = 0; z < vec_y[i].size(); z++)
        {  

        //cout << "Break point 6" << endl;

        if(z < vec_y[i].size() - 1) 
            myfile << vec_y[i][z] << ",";
        else
            myfile << vec_y[i][z];

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