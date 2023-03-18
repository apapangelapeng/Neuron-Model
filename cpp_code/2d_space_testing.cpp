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
vector<double> vec_y; 
vector<double> vec_valid; 
vector<double> vec_temp; 
vector<vector<double> > vec_coords;
vector<vector<vector<double> > > vec_time;

// build a sphere via enterring radius, populate via finding area of sphere and dividing by 1000 or so

int reset_vecs(int x){
    vec_x.clear();
    vec_y.clear();
    vec_valid.clear();
    vec_time.clear();
    return(0);
}

// double access_coordinates(double x, double y){
//     double temp;
//     int counter_time;

//     bool bool_x, bool_y;
//     for(double i = 0; i <= 9; i ++){
//         for(double j = 0; j <= 9; j ++){
//                 bool_x = (vec_time[i][y] == x) ? true : false;
//                 bool_y = (vec_time[x][j] == y) ? true : false;
//                 if(bool_x){
//                     // cout << "hi " << endl;
//                 }
//                 if(bool_x && bool_y){
//                     cout << "found coordinates: " << x << "," << y << " == " << vec_time[i][j] << endl;
//                 }
//         }
//     }
//     return(1);
// }

int display_vectors(int x){
    double temp;
    // for(const vector<double>& row : vec_time){
    //     for(double col : row){
    //         cout << col << " ";
    //     }
    //     cout << endl;
    // }
    //cout << "test test " << vec_time[1][0] << endl;

    //temp = access_coordinates(3,3);
    cout << "accessed point " << vec_time[3][3][3] << endl;
    // cout << temp << endl;
    return(0);
}

int fill_vectors(int temp){
    cout << "high" << endl;
    double x_max = 9;
    double y_max = 9;
    double time_max = 4;

    for(int time_temp = 0; time_temp <= time_max; time_temp++){
        for(double i = 0; i <= x_max; i++){
            for(double j = 0; j <= y_max; j++){
                    vec_x.push_back(i);
                    vec_y.push_back(j);
                // if((i >= 3 && i <= 6) && (j >= 3 && j <= 6)){
                //     vec_valid.push_back(1);
                // }
                // else{
                //     vec_valid.push_back(0);
                // }
            }
            vec_coords.push_back(vec_x);
            vec_coords.push_back(vec_y);
            vec_x.clear();
            vec_y.clear();
        }
        vec_time.push_back(vec_coords);
        vec_coords.clear();
    }
        //vec_time.push_back(vec_valid);
        //vec_valid.clear();
    display_vectors(0);
    cout << "low" << endl;
    return(0);
}

double output_file(double x)
{
    reset_vecs(0);
    
    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    //cout << "Break point 1" << endl;

    fill_vectors(0);

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