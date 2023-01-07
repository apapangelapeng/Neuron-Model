#if !defined(DYNAMIC_OUTPUT_H)
#define DYNAMIC_OUTPUT_H 1
#include <cmath>
#include <vector>
#include <fstream>
#include"DynamicFunc.h"


using namespace std;
int output_file(const char* path)
{
    ofstream create_file(path);
    ofstream myfile;
    myfile.open(path);

    
    vector< vector <double> > v_map = Proportion_open_test(0);

    myfile << "Tau_N,Tau_M, Tau_H,  Inf_n,Inf_M,  Inf_H \n";

        for (int i = 0; i < v_map[0].size(); i++)
    {
        for (int j = 0; j < v_map.size(); j++)
        {
            string end = (j == v_map.size() - 1) ? "\n" : ",";
            myfile << v_map[j][i] << end;
        }
    }
   
    return (0);
}

#endif