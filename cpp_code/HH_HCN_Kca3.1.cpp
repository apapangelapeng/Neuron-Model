#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <vector>
#include <algorithm>

using namespace std;

/*
THIS FILE IS JUST FOR ME TO TEST THINGS
*/

const char *path1="../data_files/pacemaker_output.csv";

double a_n, b_n, a_h, b_h, a_m, b_m;
double n_dynamic, m_dynamic, h_dynamic, n_inf, m_inf, h_inf;
vector<double> vec_n, vec_m, vec_h;
vector<double> vec_inf_n, vec_inf_m, vec_inf_h;
vector<double> vec_tau_n, vec_tau_m, vec_tau_h;
double tau_n, tau_m, tau_h;

double V_dt; 
vector<double> vec_V, vec_K_I, vec_L_I, vec_HCN_I; 
vector<double> vec_Kp, vec_HCNp;
vector<double> vec_VWT, vec_VHCN1, vec_VHCN2, vec_VHCN3; 
vector<double> vec_HCN_I1, vec_HCN_I2, vec_HCN_I3; 

double C_m = 1; 

double g_k = 36; 
double g_Na = 120;
double g_l = 0.3; 
double g_HCN = 1; 
double g_Cl = 10; //https://link.springer.com/article/10.1007/s11538-017-0289-y

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

double global_current = 0;

double Cone, Ctwo, Cthree, Cfour, Cfive, Vzero;

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
    vec_K_I.clear();
    vec_L_I.clear();
    vec_HCN_I.clear();
    return(0);
}

double dynamical_m(double V){

    a_m = (Cone*exp((V-Vzero)/Ctwo)+Cthree*(V-Vzero))/(1+Cfour*exp((V-Vzero)/Cfive));
    b_m = 
    m_inf = a_m/(a_m + b_m);
    tau_m = 1/(a_m + b_m);
    vec_tau_m.push_back(tau_m);
    vec_inf_m.push_back(m_inf);
    return(0);
}

double dynamical_h(double V){
    a_m =
    b_m = 
    m_inf = a_m/(a_m + b_m);
    tau_m = 1/(a_m + b_m);
    vec_tau_m.push_back(tau_m);
    vec_inf_m.push_back(m_inf);
    return(0);
}

double dynamical_p(double V){
    a_m =
    b_m = 
    m_inf = a_m/(a_m + b_m);
    tau_m = 1/(a_m + b_m);
    vec_tau_m.push_back(tau_m);
    vec_inf_m.push_back(m_inf);
    return(0);
}

double dynamical_d(double V){
    a_m =
    b_m = 
    m_inf = a_m/(a_m + b_m);
    tau_m = 1/(a_m + b_m);
    vec_tau_m.push_back(tau_m);
    vec_inf_m.push_back(m_inf);
    return(0);
}

double dynamical_f(double V){
    a_m =
    b_m = 
    m_inf = a_m/(a_m + b_m);
    tau_m = 1/(a_m + b_m);
    vec_tau_m.push_back(tau_m);
    vec_inf_m.push_back(m_inf);
    return(0);
}

double Static_HCN_AP(int HCN_num){

    double V_start = 0;
    double current = 0;
    double V_temp;
    double K_I_temp, L_I_temp, HCN_I_temp = 0;
    int x = 0; 

    vec_V.push_back(V_start);
    vec_n.push_back(0);
    vec_m.push_back(0);
    vec_h.push_back(0);

    for(double i = 0; i <= 100; i += delta_t){

        //cout << "Break point 1" << endl;

        dynamical_m(vec_V[x]);
        dynamical_h(vec_V[x]);
        dynamical_n(vec_V[x]);

        //cout << "Break point 2" << endl;

        vec_n.push_back(vec_n[x] + delta_t*((vec_inf_n[x] - vec_n[x])/vec_tau_n[x]));
        vec_m.push_back(vec_m[x] + delta_t*((vec_inf_m[x] - vec_m[x])/vec_tau_m[x]));
        vec_h.push_back(vec_h[x] + delta_t*((vec_inf_h[x] - vec_h[x])/vec_tau_h[x]));

        //cout << "Break point 3" << endl;

        K_I_temp = (g_k*pow(vec_n[x+1],4)*((vec_V[x]) - E_k));
        L_I_temp = (g_l*((vec_V[x]) - E_l));

        //cout << "Break point 4" << endl;
        if(HCN_num == 0){
            HCN_I_temp = (g_HCN*pow(vec_h[x+1],1)*(vec_V[x] - E_HCN1));
        }
        if(HCN_num == 1){
            HCN_I_temp = (g_HCN*pow(vec_h[x+1],1)*(vec_V[x] - E_HCN2));
        }
        if(HCN_num == 2){
            HCN_I_temp = (g_HCN*pow(vec_h[x+1],1)*(vec_V[x] - E_HCN3));
        }


        //cout << V << endl;

        V_dt = (current - K_I_temp - L_I_temp - HCN_I_temp)/C_m;
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

        vec_K_I.push_back(K_I_temp);
        vec_L_I.push_back(L_I_temp);
        vec_HCN_I.push_back(HCN_I_temp);


        //cout << V << endl;
        //cout << V_temp << endl;
        x += 1; 
        //cout << x << endl;
    }
    return(0);
}

double output_HCN_Static_AP(double x)
{
    reset_vecs(0);
    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    for(int i = 0; i <= 2; i++){
    Static_HCN_AP(i);
    }

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
    bool anode_break;

    //cout << "Break point 4" << endl;

    myfile << "V_HCN1,V_HCN2,V_HCN3\n";
    for (int i = 0; i < max_size; i++)
    {
        //cout << "Break point 5" << endl;
        HCN_I1 = (vec_VHCN1.size() > i) ? true : false;
        HCN_I2 = (vec_VHCN2.size() > i) ? true : false;
        HCN_I3 = (vec_VHCN3.size() > i) ? true : false;
        //cout << "Break point 6" << endl;

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

  output_HCN_Static_AP(0);

  cout << "End" << endl;
}