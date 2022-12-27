#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <vector>

using namespace std;

/*
THIS FILE IS JUST FOR ME TO TEST THINGS
*/

#include <sstream>

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

const char *path1 = "../data_files/test_output.csv";
string path2;

double a_n, b_n, a_h, b_h, a_m, b_m;
double n_dynamic, m_dynamic, h_dynamic, n_inf, m_inf, h_inf;
vector<double> vec_n, vec_m, vec_h;
vector<double> vec_inf_n, vec_inf_m, vec_inf_h;
vector<double> vec_Nap, vec_Kp;
vector<double> vec_tau_n, vec_tau_m, vec_tau_h;
double tau_n, tau_m, tau_h;

double V_dt;
vector<double> vec_V, vec_Na_I, vec_K_I, vec_L_I;

double C_m = 1;

double g_k = 36;
double g_Na = 120;
double g_l = 0.3;

double E_k = -12;
double E_Na = 120;
double E_l = 10.6;

double delta_t = 0.1;
double current;

/*
n is a kinetic equation built to track ONE kind of POTASSIUM channel's opening
it will be a proportion of channels open, and or the activation of the channel
potassium current can be generalized as: g*n^4*(V - E)
Where V is resting voltage, E is membrane potential, and g is conductance

m is the same, but for Na
h represents INACTIVATION of Na channels, so h and m compete

The formulas are taken from page 37 of the Electrophysiology textbook
*/
int reset_vecs(int x)
{
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
    return (0);
}

double dynamical_h(double V)
{
    // for some reason I decided to add in double h, so that we can store this value local to the method
    // and, h_dynamic can be shared to main. I forget why I did this lol, I probably had something in mind
    a_h = 0.07 * exp(-V / 20);
    b_h = ((1) / (exp((30 - V) / 10) + 1));
    h_inf = a_h / (a_h + b_h);
    tau_h = 1 / (a_h + b_h);
    // the if statements are added to expunge NaN from the data
    // storing the values in vectors so that they can be easily written to a csv file
    vec_tau_h.push_back(tau_h);
    vec_inf_h.push_back(h_inf);
    return (0);
}

double dynamical_n(double V)
{
    a_n = 0.01 * ((10 - V) / (exp((10 - V) / 10) - 1));
    if (V == 10)
    {
        a_n = 0.1; // This is the Taylor approx value for when divide by 0
    }
    b_n = 0.125 * exp(-V / 80);
    n_inf = a_n / (a_n + b_n);
    tau_n = 1 / (a_n + b_n);
    // cout << tau_n << endl;
    vec_tau_n.push_back(tau_n);
    vec_inf_n.push_back(n_inf);
    // cout << V << endl;
    return (0);
}

double dynamical_m(double V)
{
    a_m = 0.1 * ((25 - V) / (exp((25 - V) / 10) - 1));
    if (V == 25)
    {
        a_m = 1; // This is the Taylor approx value for when divide by 0
    }
    b_m = 4 * exp(-V / 18);
    m_inf = a_m / (a_m + b_m);
    tau_m = 1 / (a_m + b_m);
    // cout << tau_m << endl;
    vec_tau_m.push_back(tau_m);
    vec_inf_m.push_back(m_inf);
    return (0);
}

double Proportion_open_test(int x)
{
    // This scales between membrane voltages of -40 to 100, which is near the typical operating range of neurons
    for (double i = -40; i <= 100; i++)
    {
        double V_temp = i;
        double m_temp, h_temp, n_temp;
        // cout << V_temp << endl;
        // inputs a, b, and c are dependent on the type of channel
        // Sodium channels will only utilize a and b
        // Potassium channels will only use c
        m_temp = dynamical_m(V_temp);
        h_temp = dynamical_h(V_temp);
        n_temp = dynamical_n(V_temp);
    }
    return (0);
}

int output_file(int x)
{
    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    Proportion_open_test(0);

    myfile << "Tau_n, ";
    for (int i = 0; i < vec_tau_n.size(); i++)
    {
        myfile << vec_tau_n[i] << ",";
    }
    myfile << "\n Tau_M, ";
    for (int i = 0; i < vec_tau_m.size(); i++)
    {
        myfile << vec_tau_m[i] << ",";
        // cout << vec_tau_m[i] << endl;
    }
    myfile << "\n Tau_H, ";
    for (int i = 0; i < vec_tau_h.size(); i++)
    {
        myfile << vec_tau_h[i] << ",";
    }
    myfile << "\n Inf_n, ";
    for (int i = 0; i < vec_inf_n.size(); i++)
    {
        myfile << vec_inf_n[i] << ",";
    }
    myfile << "\n Inf_M, ";
    for (int i = 0; i < vec_inf_m.size(); i++)
    {
        myfile << vec_inf_m[i] << ",";
    }
    myfile << "\n Inf_H, ";
    for (int i = 0; i < vec_inf_h.size(); i++)
    {
        myfile << vec_inf_h[i] << ",";
    }
    return (0);
}

double Static_AP(int arbitrary_variable)
{
    reset_vecs(0);

    double V_start = 0;
    double V_temp;
    double local_cur;
    double Na_I_temp, K_I_temp, L_I_temp;
    int x = 0;

    vec_V.push_back(V_start);
    vec_n.push_back(0);
    vec_m.push_back(0);
    vec_h.push_back(0);

    for (double i = 0; i <= 10; i += delta_t)
    {
        if (i <= 4 && i >= 2)
        {
            local_cur = current;
        }
        else
        {
            local_cur = 0;
        }

        // cout << "Break point 1" << endl;

        dynamical_m(vec_V[x]);
        dynamical_h(vec_V[x]);
        dynamical_n(vec_V[x]);

        // cout << "Break point 2" << endl;

        vec_n.push_back(vec_n[x] + delta_t * ((vec_inf_n[x] - vec_n[x]) / vec_tau_n[x]));
        vec_m.push_back(vec_m[x] + delta_t * ((vec_inf_m[x] - vec_m[x]) / vec_tau_m[x]));
        vec_h.push_back(vec_h[x] + delta_t * ((vec_inf_h[x] - vec_h[x]) / vec_tau_h[x]));

        // cout << "Break point 3" << endl;

        K_I_temp = (g_k * pow(vec_n[x + 1], 4) * ((vec_V[x]) - E_k));
        Na_I_temp = (g_Na * pow(vec_m[x + 1], 3) * pow(vec_h[x + 1], 1) * ((vec_V[x]) - E_Na));
        L_I_temp = (g_l * ((vec_V[x]) - E_l));

        // cout << V << endl;

        V_dt = (local_cur - K_I_temp - Na_I_temp - L_I_temp) / C_m;
        vec_V.push_back(vec_V[x] + delta_t * V_dt);

        vec_Na_I.push_back(Na_I_temp);
        vec_K_I.push_back(K_I_temp);
        vec_L_I.push_back(L_I_temp);
        // cout << V << endl;
        // cout << V_temp << endl;
        x += 1;
        // cout << x << endl;
    }
    return (0);
}

double Run_time(double x)
{
    ofstream create_file(path2);
    ofstream myfile;
    myfile.open(path2);

    reset_vecs(0);
    Static_AP(0);

    myfile << "V,";
    for (int i = 0; i < vec_V.size(); i++)
    {
        myfile << vec_V[i] << ",";
        // cout << vec_V[i] << endl;
    }
    myfile << "\n K_I, ";
    for (int i = 0; i < vec_K_I.size(); i++)
    {
        myfile << vec_K_I[i] << ",";
        // cout << vec_K_I[i] << endl;
    }
    myfile << "\n Na_I,";
    for (int i = 0; i < vec_Na_I.size(); i++)
    {
        myfile << vec_Na_I[i] << ",";
        // cout << vec_V[i] << endl;
    }
    myfile << "\n L_I,";
    for (int i = 0; i < vec_L_I.size(); i++)
    {
        myfile << vec_L_I[i] << ",";
        // cout << vec_V[i] << endl;
    }
    myfile.close();
    return (x);
}

int main(int argc, char *argv[])
{
    /*
    Program workflow is main -calls-> ouput_file -calls-> Porportion_open for Sodium and Potassium individually
    Proportion_open -calls-> the dynamical variables, which stores the values in vectors, and also calculates
    the actual proportion open.
    Output_file then writes then info to TestingDynamicVars.csv
    */
    if (argc < 2)
    {
        throw std::invalid_argument("Did not pass in a valid current argument");
        return 1;
    }
    try
    {
        // Block of code to try
        current = stod(argv[1]); // grab the first val and convert to double
        
        path2 = string("../data_files/testV_output_")+ to_string_with_precision(current,2)+ string(".csv");
    }
    catch (invalid_argument)
    {
        throw std::invalid_argument("The given command line input is not a double for current.");
    }

    cout << to_string_with_precision(current,2)<< endl;
    output_file(0);
    Run_time(0);
}
