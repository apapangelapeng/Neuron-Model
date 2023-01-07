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

const char *path1="../data_files/static_vcompare_output.csv";


double a_n, b_n, a_h, b_h, a_m, b_m;
double n_dynamic, m_dynamic, h_dynamic, n_inf, m_inf, h_inf;
vector<double> vec_n, vec_m, vec_h;
vector<double> vec_inf_n, vec_inf_m, vec_inf_h;
vector<double> vec_tau_n, vec_tau_m, vec_tau_h;
double tau_n, tau_m, tau_h;

double V_dt; 
vector<double> vec_V, vec_Na_I, vec_K_I, vec_L_I, vec_HCN_I, vec_Cl_I; 
vector<double> vec_Nap, vec_Kp, vec_HCNp;
vector<double> vec_VWT, vec_VHCN1, vec_VHCN2, vec_VHCN3; 
vector<double> vec_HCN_I1, vec_HCN_I2, vec_HCN_I3; 

double C_m = 1; 

double g_k = 36; 
double g_Na = 120;
double g_l = 0.3; 
double g_HCN = 0.1; 
double g_Cl; //https://link.springer.com/article/10.1007/s11538-017-0289-y

double g_cl_max = 0; 

double E_k = -12; 
double E_Na = 120; 
double E_l = 10.6; 
double E_HCN3 = -43; 
double E_HCN1 = -1; 
double E_HCN2 = -21; 
double E_Cl = -10; 

double b_Cl = 50; 
double t_Cl = 0.02; 

double delta_t = 0.001;

double global_current = 5;

default_random_engine generator;

normal_distribution<double> error(1,0.025);

int reset_vecs(int x){
    vec_V.clear();
    vec_n.clear();
    vec_m.clear();
    vec_h.clear();
    vec_tau_m.clear();
    vec_inf_m.clear();
    vec_tau_h.clear();
    vec_inf_h.clear();
    vec_tau_n.clear();
    vec_inf_n.clear();
    vec_Na_I.clear();
    vec_K_I.clear();
    vec_L_I.clear();
    vec_HCN_I.clear();
    vec_Cl_I.clear();
    return(0);
}

double dynamical_h(double V){
    //for some reason I decided to add in double h, so that we can store this value local to the method
    //and, h_dynamic can be shared to main. I forget why I did this lol, I probably had something in mind
    a_h = 0.07*exp(-V/20);
    b_h = ((1)/(exp((30-V)/10) + 1));
    h_inf = a_h/(a_h + b_h);
    tau_h = 1/(a_h + b_h);
    //the if statements are added to expunge NaN from the data
    //storing the values in vectors so that they can be easily written to a csv file
    vec_tau_h.push_back(tau_h);
    vec_inf_h.push_back(h_inf);
    return(0);
}

double dynamical_n(double V){
    a_n = 0.01*((10-V)/(exp((10-V)/10) - 1));
    if(V == 10){
        a_n = 0.1; // This is the Taylor approx value for when divide by 0
    }
    b_n = 0.125*exp(-V/80);
    n_inf = a_n/(a_n + b_n);
    tau_n = 1/(a_n + b_n);
    //cout << tau_n << endl;
    vec_tau_n.push_back(tau_n);
    vec_inf_n.push_back(n_inf);
    //cout << V << endl;
    return (0);
}

double dynamical_m(double V){
    a_m = 0.1*((25 - V)/(exp((25-V)/10) - 1));
    if (V == 25){
        a_m = 1; // This is the Taylor approx value for when divide by 0
    }
    b_m = 4*exp(-V/18);
    m_inf = a_m/(a_m + b_m);
    tau_m = 1/(a_m + b_m);
    //cout << tau_m << endl;
    vec_tau_m.push_back(tau_m);
    vec_inf_m.push_back(m_inf);
    return (0);
}

double Proportion_open_test(int x){
    //This scales between membrane voltages of -40 to 100, which is near the typical operating range of neurons
    for(double i = -40; i <= 100; i++){
        double V_temp = i;
        double m_temp, h_temp, n_temp;
        //cout << V_temp << endl;
        //inputs a, b, and c are dependent on the type of channel
        //Sodium channels will only utilize a and b
        //Potassium channels will only use c
            m_temp = dynamical_m(V_temp);
            h_temp = dynamical_h(V_temp);
            n_temp = dynamical_n(V_temp);
        }
  return (0);
}

int output_file(int x){
    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    //Proportion_open_test(0);

    myfile << "Tau_n, ";
    for(int i = 0; i < vec_tau_n.size(); i++){
       myfile << vec_tau_n[i] << ","; 
    }
    myfile << "\n Tau_M, ";
    for(int i = 0; i < vec_tau_m.size(); i++){
       myfile << vec_tau_m[i] << ","; 
        //cout << vec_tau_m[i] << endl;
    }
    myfile << "\n Tau_H, ";
    for(int i = 0; i < vec_tau_h.size(); i++){
        myfile << vec_tau_h[i] << ","; 
    }
    myfile << "\n Inf_n, ";
    for(int i = 0; i < vec_inf_n.size(); i++){
       myfile << vec_inf_n[i] << ","; 
    }
    myfile << "\n Inf_M, ";
    for(int i = 0; i < vec_inf_m.size(); i++){
       myfile << vec_inf_m[i] << ","; 
    }
    myfile << "\n Inf_H, ";
    for(int i = 0; i < vec_inf_h.size(); i++){
        myfile << vec_inf_h[i] << ","; 
    }
    return(0);
}

double Static_HCN_AP(int HCN_num){
    

    double V_start = 0;
    double current;
    double V_temp;
    double Na_I_temp, K_I_temp, L_I_temp, HCN_I_temp, Cl_I_temp = 0;
    int x = 0; 

    vec_V.push_back(V_start);
    vec_n.push_back(0);
    vec_m.push_back(0);
    vec_h.push_back(0);

    vec_Cl_I.push_back(0);

    for(double i = 0; i <= 100; i += delta_t){
        //||  10 <= i <= 12 || 20 <= i <= 22 || 30 <= i <= 32 || 40 <= i <= 42
        if(i <= 50 && i >= 45){
            current = global_current;
            //cout << x << endl; 
        }
        else{
            current = 0;
        }

        if(i <= 55 && i >= 45){
            g_Cl = g_cl_max;
        }
        else{
            g_Cl = g_cl_max;
        }


        double error_applied = error(generator);

        //cout << "Break point 1" << endl;

        dynamical_m(vec_V[x]);
        dynamical_h(vec_V[x]);
        dynamical_n(vec_V[x]);

        //cout << "Break point 2" << endl;

        vec_n.push_back(vec_n[x] + delta_t*((vec_inf_n[x] - vec_n[x])/vec_tau_n[x]));
        vec_m.push_back(vec_m[x] + delta_t*((vec_inf_m[x] - vec_m[x])/vec_tau_m[x]));
        vec_h.push_back(vec_h[x] + delta_t*((vec_inf_h[x] - vec_h[x])/vec_tau_h[x]));

        //cout << "Break point 3" << endl;

        K_I_temp = (error_applied)*(g_k*pow(vec_n[x+1],4)*((vec_V[x]) - E_k));
        Na_I_temp = (error_applied)*(g_Na*pow(vec_m[x+1],3)*pow(vec_h[x+1],1)*((vec_V[x]) - E_Na));
        L_I_temp = (error_applied)*(g_l*((vec_V[x]) - E_l));

        Cl_I_temp = (error_applied)*t_Cl*((g_Cl*pow(2.71828,(vec_V[x]-E_Cl)/b_Cl))*(vec_V[x] - E_Cl) - vec_Cl_I[x]);

        //cout << "Break point 4" << endl;
        if(HCN_num == 0){
            HCN_I_temp = (error_applied)*(g_HCN*pow(vec_h[x+1],1)*(vec_V[x] - E_HCN1));
        }
        if(HCN_num == 1){
            HCN_I_temp = (error_applied)*(g_HCN*pow(vec_h[x+1],1)*(vec_V[x] - E_HCN2));
        }
        if(HCN_num == 2){
            HCN_I_temp = (error_applied)*(g_HCN*pow(vec_h[x+1],1)*(vec_V[x] - E_HCN3));
        }


        //cout << V << endl;

        V_dt = (current - K_I_temp - Na_I_temp - L_I_temp - HCN_I_temp - Cl_I_temp)/C_m;
        vec_V.push_back(vec_V[x] + delta_t*V_dt);

        if(HCN_num == 0){
            vec_VHCN1.push_back(vec_V[x] + delta_t*V_dt);
            vec_HCN_I1.push_back(HCN_I_temp);
        }
        if(HCN_num == 1){
            vec_VHCN2.push_back(vec_V[x] + delta_t*V_dt);
            vec_HCN_I2.push_back(HCN_I_temp);
        }
        if(HCN_num == 2){
            vec_VHCN3.push_back(vec_V[x] + delta_t*V_dt);
            vec_HCN_I3.push_back(HCN_I_temp);
        }


        vec_Na_I.push_back(Na_I_temp);
        vec_K_I.push_back(K_I_temp);
        vec_L_I.push_back(L_I_temp);
        vec_Cl_I.push_back(Cl_I_temp);

        vec_HCN_I.push_back(HCN_I_temp);


        //cout << V << endl;
        //cout << V_temp << endl;
        x += 1; 
        //cout << x << endl;
    }
    vec_Cl_I.erase(vec_Cl_I.begin());
    return(0);
}


double Static_WT_AP(int arbitrary_variable){
    double V_start = 0;
    double current;
    double V_temp;
    double Na_I_temp, K_I_temp, L_I_temp, HCN_I_temp, Cl_I_temp = 0;
    int x = 0; 

    vec_V.push_back(V_start);
    vec_n.push_back(0);
    vec_m.push_back(0);
    vec_h.push_back(0);

    vec_Cl_I.push_back(0);

    for(double i = 0; i <= 100; i += delta_t){
        //||  10 <= i <= 12 || 20 <= i <= 22 || 30 <= i <= 32 || 40 <= i <= 42
        if(i <= 55 && i >= 45){
            current = global_current;
            //cout << x << endl; 
        }
        else{
            current = 0;
        }
    
        if(i <= 55 && i >= 45){
            g_Cl = g_cl_max;
        }
        else{
            g_Cl = g_cl_max;
        }

        double error_applied = error(generator);

        //cout << "Break point 1" << endl;

        dynamical_m(vec_V[x]);
        dynamical_h(vec_V[x]);
        dynamical_n(vec_V[x]);

        //cout << "Break point 2" << endl;

        vec_n.push_back(vec_n[x] + delta_t*((vec_inf_n[x] - vec_n[x])/vec_tau_n[x]));
        vec_m.push_back(vec_m[x] + delta_t*((vec_inf_m[x] - vec_m[x])/vec_tau_m[x]));
        vec_h.push_back(vec_h[x] + delta_t*((vec_inf_h[x] - vec_h[x])/vec_tau_h[x]));

        //cout << "Break point 3" << endl;

        K_I_temp = (error_applied)*(g_k*pow(vec_n[x+1],4)*((vec_V[x]) - E_k));
        Na_I_temp = (error_applied)*(g_Na*pow(vec_m[x+1],3)*pow(vec_h[x+1],1)*((vec_V[x]) - E_Na));
        L_I_temp = (error_applied)*(g_l*((vec_V[x]) - E_l));

        Cl_I_temp = (error_applied)*t_Cl*((g_Cl*pow(2.71828,(vec_V[x]-E_Cl)/b_Cl))*(vec_V[x] - E_Cl) - vec_Cl_I[x]);

        //cout << "Break point 4" << endl;

        //cout << V << endl;

        V_dt = (current - K_I_temp - Na_I_temp - L_I_temp - Cl_I_temp)/C_m;
        vec_V.push_back(vec_V[x] + delta_t*V_dt);

        vec_VWT.push_back(vec_V[x] + delta_t*V_dt);

        vec_Na_I.push_back(Na_I_temp);
        vec_K_I.push_back(K_I_temp);
        vec_L_I.push_back(L_I_temp);
        vec_Cl_I.push_back(Cl_I_temp);


        //cout << V << endl;
        //cout << V_temp << endl;
        x += 1; 
        //cout << x << endl;
    }
    vec_Cl_I.erase(vec_Cl_I.begin());
    return(0);
}


double voltage_output(double x)
{
    Static_WT_AP(0);
    reset_vecs(0);

    for (int i = 0; i < 3; i++)
    {
        cout << i << endl; 
        Static_HCN_AP(i);
        reset_vecs(0);
    }

    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    vector<int> sizes;
    /*cout<<vec_V.size()<<endl;
    cout<<vec_N.size()<<endl;
    cout<<vec_tiny_N.size()<<endl;*/
    sizes.insert(sizes.begin(),vec_VWT.size());
    sizes.insert(sizes.begin(),vec_VHCN1.size());
    sizes.insert(sizes.begin(),vec_VHCN2.size());
    sizes.insert(sizes.begin(),vec_VHCN3.size());
    sizes.insert(sizes.begin(),vec_HCN_I1.size());
    sizes.insert(sizes.begin(),vec_HCN_I2.size());
    sizes.insert(sizes.begin(),vec_HCN_I3.size());

    sort(sizes.begin(), sizes.end());
    int max_size = sizes.back();
    
    cout << max_size << endl;

    bool VWT; 
    bool VHCN1; 
    bool VHCN2; 
    bool VHCN3; 
    bool HCN_I1; 
    bool HCN_I2; 
    bool HCN_I3; 
    
    //cout << "Break point 4" << endl;

    myfile << "V_WT,V_HCN1,V_HCN2,V_HCN3,HCN_I1,HCN_I2,HCN_I3\n";
    for (int i = 0; i < max_size; i++)
    {
        //cout << "Break point 5" << endl;
        VWT = (vec_VWT.size() > i) ? true : false;
        VHCN1 = (vec_VHCN1.size() > i) ? true : false;
        VHCN2 = (vec_VHCN2.size() > i) ? true : false;
        VHCN3 = (vec_VHCN3.size() > i) ? true : false;
        HCN_I1 = (vec_VHCN1.size() > i) ? true : false;
        HCN_I2 = (vec_VHCN2.size() > i) ? true : false;
        HCN_I3 = (vec_VHCN3.size() > i) ? true : false;

        //cout << "Break point 6" << endl;

        if(VWT) myfile << vec_VWT[i] << "," ;
        if(!VWT) myfile << "," ;
        if(VHCN1) myfile << vec_VHCN1[i] << ",";
        if(!VHCN1) myfile <<"," ;
        if(VHCN2) myfile << vec_VHCN2[i] << ",";
        if(!VHCN2) myfile << "," ;
        if(VHCN3) myfile << vec_VHCN3[i] << ",";
        if(!VHCN3) myfile << "," ;
        if(HCN_I1) myfile << vec_HCN_I1[i] << ",";
        if(!HCN_I1) myfile <<"," ;
        if(HCN_I2) myfile << vec_HCN_I2[i] << ",";
        if(!HCN_I2) myfile << "," ;
        if(HCN_I3) myfile << vec_HCN_I3[i];


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
