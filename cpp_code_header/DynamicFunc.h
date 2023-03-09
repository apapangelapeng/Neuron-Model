#if !defined(DYNAMIC_FUNC_H)
#define DYNAMIC_FUNC_H 1
#include <cmath>
#include <vector>


using namespace std;
void dynamical_h(double V,vector<double>& vec_tau_h,vector<double>& vec_inf_h)
{
    // for some reason I decided to add in double h, so that we can store this value local to the method
    // and, h_dynamic can be shared to main. I forget why I did this lol, I probably had something in mind
    double a_h = 0.07 * exp(-V / 20);
    double b_h = ((1) / (exp((30 - V) / 10) + 1));
    double h_inf = a_h / (a_h + b_h);
    double tau_h = 1 / (a_h + b_h);
    // the if statements are added to expunge NaN from the data
    // storing the values in vectors so that they can be easily written to a csv file
    vec_tau_h.push_back(tau_h);
    vec_inf_h.push_back(h_inf);

}

void dynamical_n(double V,vector<double>& vec_tau_n,vector<double>& vec_inf_n)
{
    double a_n = 0.01 * ((10 - V) / (exp((10 - V) / 10) - 1));
    if (V == 10)
    {
        a_n = 0.1; // This is the Taylor approx value for when divide by 0
    }
    double b_n = 0.125 * exp(-V / 80);
    double n_inf = a_n / (a_n + b_n);
    double tau_n = 1 / (a_n + b_n);
    // cout << tau_n << endl;
    vec_tau_n.push_back(tau_n);
    vec_inf_n.push_back(n_inf);
    // cout << V << endl;
  
}

void dynamical_m(double V,vector<double>& vec_tau_m,vector<double>& vec_inf_m)
{
    double a_m = 0.1 * ((25 - V) / (exp((25 - V) / 10) - 1));
    if (V == 25)
    {
        a_m = 1; // This is the Taylor approx value for when divide by 0
    }
    double b_m = 4 * exp(-V / 18);
    double m_inf = a_m / (a_m + b_m);
    double tau_m = 1 / (a_m + b_m);
    // cout << tau_m << endl;
    vec_tau_m.push_back(tau_m);
    vec_inf_m.push_back(m_inf);
}

vector< vector <double> > Proportion_open_test(int x)
{
    vector<double> vec_tau_n, vec_tau_m, vec_tau_h;
    vector<double> vec_inf_n, vec_inf_m, vec_inf_h;
    // This scales between membrane voltages of -40 to 100, which is near the typical operating range of neurons
    for (double i = -40; i <= 100; i++)
    {
        double V_temp = i;

        // cout << V_temp << endl;
        // inputs a, b, and c are dependent on the type of channel
        // Sodium channels will only utilize a and b
        // Potassium channels will only use c
        dynamical_m(V_temp,vec_tau_m,vec_inf_m);
        dynamical_h(V_temp,vec_tau_h,vec_inf_h);
        dynamical_n(V_temp,vec_tau_n,vec_inf_n);
    }
   vector< vector <double> > dynamic_vecs ;
   dynamic_vecs.push_back(vec_tau_n);
   dynamic_vecs.push_back(vec_tau_m);
   dynamic_vecs.push_back(vec_tau_h);
   dynamic_vecs.push_back(vec_inf_n);
   dynamic_vecs.push_back(vec_inf_m);
   dynamic_vecs.push_back(vec_inf_h);
   return dynamic_vecs;
    
}

#endif