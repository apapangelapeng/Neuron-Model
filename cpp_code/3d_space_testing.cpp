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
vector<vector<double> > vec_y;
vector<vector<vector<double> > > vec_z;
vector<vector<vector<vector<double> > > > vec_time; 

// build a sphere via enterring radius, populate via finding area of sphere and dividing by 1000 or so

int reset_vecs(int x){
    vec_x.clear();
    vec_y.clear();
    vec_z.clear();
    vec_time.clear();
    return(0);
}

int display_vectors(int x){
    for(int i = 0; i < vec_time.size(); i++){
        cout << vec_time[i][i][i][i] << endl;
    }
    cout << "end here!" << endl;
    return(0);
}

int fill_vectors(int x_range, int y_range, int z_range, int time_range){
    for (int i = 0; i <= x_range; i++) { 
        vec_x.push_back(i);
    }
    for (int i = 0; i <= y_range; i++) { 
        vec_y.push_back(vec_x);
    }
    for (int i = 0; i <= y_range; i++) { 
        vec_z.push_back(vec_y);
    }
    for (int i = 0; i <= time_range; i++) { 
        vec_time.push_back(vec_z);
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

    fill_vectors(5, 5, 5, 5);

    vector<int> sizes;

    cout << "Break point 2" << endl;

    sizes.insert(sizes.begin(),vec_x.size());
    sizes.insert(sizes.begin(),vec_y.size());
    sizes.insert(sizes.begin(),vec_z.size());
    sizes.insert(sizes.begin(),vec_time.size());

    cout << "Break point 3" << endl;

    sort(sizes.begin(), sizes.end());
    int max_size = sizes.back();
    
    //cout << max_size << endl;

    bool bool_x; 
    bool bool_y; 
    bool bool_z; 
    bool bool_time; 

    //cout << "Break point 4" << endl;

    myfile << "x,y,z,time\n";
    for (int i = 0; i < max_size; i++)
    {
        //cout << "Break point 5" << endl;

        bool_x = (vec_x.size() > i) ? true : false;
        bool_y = (vec_y.size() > i) ? true : false;
        bool_z = (vec_z.size() > i) ? true : false;
        bool_time = (vec_time.size() > i) ? true : false;


        cout << "Break point 6" << endl;

        if(bool_x) myfile << vec_x[i] <<"," ;
        if(!bool_x) myfile << "," ;
        if(bool_y) myfile << vec_y[i][i] << ",";
        if(!bool_y) myfile << ",";
        if(bool_z) myfile << vec_z[i][i][i] << ",";
        if(!bool_z) myfile << ",";
        if(bool_time) myfile << vec_time[i][i][i][i];

        myfile << "\n";

        cout << "Break point 7" << endl;

    }
    myfile.close();
    return (x);
}

int main(void) {
    cout << "Begin" << endl;
    output_file(0);
    cout << "End" << endl;
}