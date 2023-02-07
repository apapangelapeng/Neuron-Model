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


//CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double log_convert = 2.303; // to convert from ln to log10
double R_constant = 8.1345;
double F = 96845;
double body_temp = 310.15;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double temp;

// calcium concentration: 2.4mM insid, 100nM outside https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3553253/#:~:text=Extracellular%20calcium%2C%20and%20particularly%20the,calcium%20of%201.1%E2%80%931.4%20mM
double Ca_in, Ca_out;
double E_Ca = 131.373; // this is for humans, i.e., body temp of 310K etc. Unsure what it is for Drosophila


// Piezo Kinetics %%%%%%%%%%%%%%%%%%%%%%%%%
double G_Piezo = 0.000000000030; 
// Piezo1 = 29pS https://www.sciencedirect.com/science/article/pii/S0968000416301505
// The decay rate, according to this paper, is 1ms
// 25-30pS https://anatomypubs.onlinelibrary.wiley.com/doi/full/10.1002/dvdy.401
// Time constant = 6.2 ± 0.3 ms https://www.nature.com/articles/nature10812
// This paper also says conductance is 30pS and that Drosophila is closer to 3.3pS....
// Piezo1 (16ish ms) vs Piezo2 (7ish ms) inactiation kinetics https://www.science.org/doi/full/10.1126/science.1193270?casa_token=WnOGIrz8PbAAAAAA%3AdA422PZ9tigaVHLD6IxvFlWpQer0KOFm_jvCrohrhgnSzzdEPOXhZo_koo-cZeWbLeNSTFPtSJ3bdtA
// Force gating for Piezo1 is about –30 mmHg or 1.4 mN/m https://www.sciencedirect.com/science/article/pii/S0968000421000220#bb0270
// Piezo1 half-maximal activation T50 = 2.7 ± 0.1 mN/m https://elifesciences.org/articles/12088
// Piezo1 P50 of –28.1 ± 2.8 and –31.2 ± 3.5 mmHg https://www.science.org/doi/full/10.1126/science.1193270?casa_token=gKZRiU1R_vAAAAAA%3AAuhBFZP-QjXckn7G9aBL6A_ZJfCjsRrUIqoXCbBA1887i29wPWtzlMBfdwShr45kBM7Pj-N4NYVlfEQ


// Piezo dimensions %%%%%%%%%%%%%%%%%%%%%%%
// Blade = 200A https://www.nature.com/articles/nature25453
// SA, then, is 17320.51 A^2, or 173.20 nm^2, or 72658.95
// If we assume the growth cone of a sphere is 1um, then the SA is 12.57 um^2
// I think that means that 72658.95 Piezo proteins can be on the surface


// Overall Definitions %%%%%%%%%%%%%%%%%%%%%
double Ca_c_dT; //derivative of concentration of calcium per time
double D_diff; //coefficient of diffusion
double Ca_c_dT_dT; //2nd derivative of Ca concentration per unit time
double J_ipr; //IPR is the IP3 receptor, and is supposedly the biggest calcium leaker
// Note: it does not seem that Piezo upregulating PLC/IP3 has been shown
double cone_circumference; //circumference of the growth cone
double cone_cross_area; //cross sectional area of the growth cone
double J_in; //flux of Ca2+ into the cell 
double J_pm; //flux pumped out of the cell
double J_ryr; //flux from RyR receptor
double J_serca; //flux from serca pump, maybe includes good SERCA model: https://link.springer.com/article/10.1140/epjp/s13360-023-03691-1
double J_on; //function of binding of Ca2+ to buffers, Page 309 of mathematical physiology seems good
double J_off; //function of unbinding of Ca2+ from buffers

// Ca_diff_dT = D_c*Ca_diff_dT_dT + J_ipr + (cone_circumference/cone_cross_area)*(J_in - J_pm) + J_ryr - J_serca - J_on + J_off; 

// RyR Definitions %%%%%%%%%%%%%%%%%%%%%%%%
// Very good reference for calcium dynamics: https://www.frontiersin.org/articles/10.3389/fphys.2012.00114/full#F1
double vol_D; //dyadic space volume
double vol_ER; //ER volume
double g_ryr; //RyR channel conductance
double c_j; //concentration inside the ER 
double c_d; //concentration in the cytoplasm
double N_ryr; //stochastic number of RyR channels

// J_ryr = N_ryr*(g_ryr/vol_D)*(c_j - c_d); //This kind of decribes the local movement due to gradient

// Buffering Definitions %%%%%%%%%%%%%%%%%
// Reference includes a list of models published by year: https://www.frontiersin.org/articles/10.3389/fncom.2018.00014/full
double buff_unbound; //concentration of unbound buffer
double buff_bound; //concentration of bound buffer
double k_buff_bind; //binding affinity/Kon of buffer
double k_buff_unbind; //unbinding affinity/Koff of buffer
double buff_c_dT; //derivative of concentration of buffer with respect to time
double buff_c_dT_dT; //second derivative of concentration of buffer with respect to time


vector<double> temp_vec;

default_random_engine generator;
normal_distribution<double> error(1,0.025);

double PotentialE(double out, double in, int Z) {
  double E = 1000 * (R_constant * body_temp) / (F * Z) * log(out / in); // log(x) = ln(x) in cpp
  if ((isinf(E)) || (E != E)) {
    cout << "YOUR E FUNCTION IS FAULTING! Probably, the concentration inside went to 0, or you entered z = 0." << endl;
  }
  //THIS IS USED TO TEST POTENTIAL CALCULATION
  cout << "\n THIS IS E in mV: " << E << endl;
  return (E);
}



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
  PotentialE(0.0024, 0.0000001, 2);
  cout << "End" << endl;
}
