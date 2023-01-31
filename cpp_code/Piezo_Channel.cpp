#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <vector>
#include <algorithm>
#include <random>

using namespace std;

/*
THIS FILE IS JUST FOR ME TO TEST THINGS
*/

const char *path1="../data_files/Piezo_Channel.csv";


double temp;

double G_Piezo = 0.000000000030; 
// Piezo1 = 29pS https://www.sciencedirect.com/science/article/pii/S0968000416301505
// The decay rate, according to this paper, is 1ms
// 25-30pS https://anatomypubs.onlinelibrary.wiley.com/doi/full/10.1002/dvdy.401
// Time constant = 6.2 ± 0.3 ms https://www.nature.com/articles/nature10812
// This paper also says conductance is 30pS and that Drosophila is closer to 3.3pS....

vector<double> temp_vec;

default_random_engine generator;
normal_distribution<double> error(1,0.025);

int reset_vecs(int x){
    return(0);
}

int Piezo_Channel(int x){
    return(0);
}


double voltage_output(double x)
{
    Piezo_Channel(0);
    reset_vecs(0);

    for (int i = 0; i < 3; i++)
    {
        cout << i << endl; 
        Piezo_Channel(i);
        reset_vecs(0);
    }

    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    vector<int> sizes;
    /*cout<<vec_V.size()<<endl;
    cout<<vec_N.size()<<endl;
    cout<<vec_tiny_N.size()<<endl;*/
    sizes.insert(sizes.begin(),temp_vec.size());

    sort(sizes.begin(), sizes.end());
    int max_size = sizes.back();
    
    cout << max_size << endl;

    bool temp_bool;

    //cout << "Break point 4" << endl;

    myfile << "Temp\n";
    for (int i = 0; i < max_size; i++)
    {
        //cout << "Break point 5" << endl;
        temp_bool = (temp_vec.size() > i) ? true : false;

        //cout << "Break point 6" << endl;

        if(temp_bool) myfile << temp_vec[i];

        //if(!dn_V_0) myfile << vec_tiny_N[i];

        myfile << "\n";
        //cout << "Break point 6" << endl;
    }
    myfile.close();
    return (x);
}

int main(void) {
  cout << "Begin" << endl;
  voltage_output(0);
  cout << "End" << endl;
}
