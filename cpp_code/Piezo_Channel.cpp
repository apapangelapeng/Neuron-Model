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
double F = 96485.3321;
double body_temp = 310.15;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//General use %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double delta_T = 0.01;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


// calcium concentration: 2.4mM outside, 100nM inside https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3553253/#:~:text=Extracellular%20calcium%2C%20and%20particularly%20the,calcium%20of%201.1%E2%80%931.4%20mM
double Ca_in, Ca_out;
double E_Ca; // 131.373 --> this is for humans, i.e., body temp of 310K etc. Unsure what it is for Drosophila


// Piezo Kinetics %%%%%%%%%%%%%%%%%%%%%%%%%
double G_Piezo_single = 0.000000000030; 
double G_Piezo_total;
int N_Piezo_channels = 2000;
double p_open = 0; 
vector<int> vec_num_open;
double p_closed = 1;
vector<int> vec_num_closed;
double Piezo_current;
vector<double> vec_Piezo_current;
double placeholder_opening_function;
int open_counter = 0;
// Piezo1 = 29pS https://www.sciencedirect.com/science/article/pii/S0968000416301505
// The decay rate, according to this paper, is 1ms
// 25-30pS https://anatomypubs.onlinelibrary.wiley.com/doi/full/10.1002/dvdy.401
// Time constant = 6.2 ± 0.3 ms https://www.nature.com/articles/nature10812
// This paper also says conductance is 30pS and that Drosophila is closer to 3.3pS....
// Piezo1 (16ish ms) vs Piezo2 (7ish ms) inactiation kinetics https://www.science.org/doi/full/10.1126/science.1193270?casa_token=WnOGIrz8PbAAAAAA%3AdA422PZ9tigaVHLD6IxvFlWpQer0KOFm_jvCrohrhgnSzzdEPOXhZo_koo-cZeWbLeNSTFPtSJ3bdtA
// Force gating for Piezo1 is about –30 mmHg or 1.4 mN/m https://www.sciencedirect.com/science/article/pii/S0968000421000220#bb0270
// Piezo1 half-maximal activation T50 = 2.7 ± 0.1 mN/m https://elifesciences.org/articles/12088
// Piezo1 P50 of –28.1 ± 2.8 and –31.2 ± 3.5 mmHg https://www.science.org/doi/full/10.1126/science.1193270?casa_token=gKZRiU1R_vAAAAAA%3AAuhBFZP-QjXckn7G9aBL6A_ZJfCjsRrUIqoXCbBA1887i29wPWtzlMBfdwShr45kBM7Pj-N4NYVlfEQ
double T_Piezo = 0.016; //16 ms tau

// Piezo dimensions %%%%%%%%%%%%%%%%%%%%%%%
// Blade = 200A https://www.nature.com/articles/nature25453
// SA, then, is 17320.51 A^2, or 173.20 nm^2, or 72658.95
// If we assume the growth cone of a sphere is 1um, then the SA is 12.57 um^2
// I think that means that 72658.95 Piezo proteins can be on the surface


// Overall Definitions %%%%%%%%%%%%%%%%%%%%%
vector<double> vec_Ca_conc; //calcium concentration
int Ca_counter = 0;
double Ca_c_dT; //derivative of concentration of calcium per time in the cytoplasm
double D_diff_Ca; //coefficient of diffusion
double Ca_c_dT_dT; //2nd derivative of Ca concentration per unit time
double J_ipr; //IPR is the IP3 receptor, and is supposedly the biggest calcium leaker
// Note: it does not seem that Piezo upregulating PLC/IP3 has been shown
double cone_radius = 0.000175; //this is in ... centimeters ... for some reason
double cone_circumference; //circumference of the growth cone
double cone_cross_area; //cross sectional area of the growth cone
double J_in; //flux of Ca2+ into the cell 
double J_pm; //flux pumped out of the cell
double J_ryr; //flux from RyR receptor
double J_serca; //flux from serca pump, maybe includes good SERCA model: https://link.springer.com/article/10.1140/epjp/s13360-023-03691-1
double J_on; //function of binding of Ca2+ to buffers, Page 309 of mathematical physiology seems good
double J_off; //function of unbinding of Ca2+ from buffers
// Ca_c_dT = D_c*Ca_c_dT_dT + J_ipr + (cone_circumference/cone_cross_area)*(J_in - J_pm) + J_ryr - J_serca - J_on + J_off; 

// RyR Definitions %%%%%%%%%%%%%%%%%%%%%%%%
// Model 1
// Very good reference for calcium dynamics: https://www.frontiersin.org/articles/10.3389/fphys.2012.00114/full#F1
double vol_D; //dyadic space volume
double vol_ER; //ER volume
double g_ryr; //RyR channel conductance
double C_er = 0.0005; //concentration inside the ER; units of M; in this case, we are assuming that C_er is relatively constant
double C_cyt; //concentration in the cytoplasm
double N_ryr; //stochastic number of RyR channels
// J_ryr = N_ryr*(g_ryr/vol_D)*(C_er - C_cyt); //This kind of decribes the local movement due to gradient

// RyR Definitions %%%%%%%%%%%%%%%%%%%%%%%%
// Model 2
// From this paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4648987/pdf/10827_2015_Article_579.pdf
double v_rel = 5; // no clue what this is, just some constant, units of s^-1
double v_leak = 0.12; // still, no clue, just some constant, units of s^-1
// c_cyt and c_er are still used 
double P_open; // probability/proportion that channels are open
double w; // dynamical variable -- I don't want to talk about it
double w_inf; // infinite version of w; w_inf is a state-function and does not need a vector
double w_dt; // w derivative with respect to time
double tau_w; // time constant of w
double K_a = 0.0000000000000000000000000192; // kinetics a, units of M^4
double K_b = 0.0000000000000000002573; // kinetics b, units of M^3
double K_c = 0.0571; // kinetics c, unitless
double K_d = 0.1; // kinetics d, units of s^-1
int w_counter = 0;
vector<double> vec_J_ryr;
// J_ryr = (v_rel*P_open + v_leak)(C_er - C_cyt);
// P_open = (w*((1 + C_cyt^3)/K_b))/((K_a/C_Cyt^4) + 1 + (C_cyt^3/K_b));
// w_inf = ((K_a/C_Cyt^4) + 1 + (C_cyt^3/K_b))/((1/K_c) + (K_a/C_Cyt^4) + 1 + (C_cyt^3/K_b));
// w_dt = (w_inf - w)/tau_w; 
// tau_w = w_inf/K_d;
// w = w + w_dt;

// SERCA Definitions %%%%%%%%%%%%%%%%%%%%%%%%
// Most SERCA models seem to simply be the Hill function
double v_serca = 0.00012; // some constant, units of M / s
double K_p = 0.0000003; // kinetics p, units of uM
vector<double> vec_J_serca;
// TECHNICALLY THiS MODEL IS OUTDATED, CHECK PAGE 284 OF MATHEMATICAL PHYS
// J_serca = v_serca*((C_cyt^2)/(C_cyt^2 + K_p^2));

// IPR Definitions %%%%%%%%%%%%%%%%%%%%%%%%
// Taking this mostly from Mathematical Physiology page 293. This is a simplified IPR model
// It is difficult to know if this will be sufficient 


// Buffering Definitions %%%%%%%%%%%%%%%%%
// Reference includes a list of models published by year: https://www.frontiersin.org/articles/10.3389/fncom.2018.00014/full
double buff_unbound = 0.0000000001; //concentration of unbound buffer, which we are taking to be b_total
double buff_bound = 0.00000001; //concentration of bound buffer
vector<double> vec_buff_bound;
int buff_counter = 0;
double k_buff_bind = 600000000; //binding affinity/Kon of buffer
double k_buff_unbind = 100; //unbinding affinity/Koff of buffer
double buff_c_dT; //derivative of concentration of buffer with respect to time
double buff_c_dT_dT; //second derivative of concentration of buffer with respect to time, we do not really ned this unless we are measuring the diffusion of the buffer
double buff_diff;
vector<double> vec_J_on;
// J_on = k_buff_bind*C_cyt*buff_unbound;
// J_off = k_buff_unbind*buff_bound; 
// buff_c_dT = k_buff_bind*C_cyt*buff_unbound - k_buff_unbind*buff_bound;
// buff_bound = buff_bound + buff_c_dT;


vector<double> temp_vec;
vector<double> vec_channel_number;
vector<double> vec_time;
vector<double> vec_w;

default_random_engine generator;
normal_distribution<double> stochastic_opening(0,1);

double PotentialE(double out, double in, int Z) {
  double E = 1000 * (R_constant * body_temp) / (F * Z) * log(out / in); // log(x) = ln(x) in cpp
  if ((isinf(E)) || (E != E)) {
    cout << "YOUR E FUNCTION IS FAULTING! Probably, the concentration inside went to 0, or you entered z = 0." << endl;
  }
  //THIS IS USED TO TEST POTENTIAL CALCULATION
  cout << "\n THIS IS E in mV: " << E << endl;
  return (E);
}

double Compute_J_on(double C_cyt){
  double scaling_factor = 1; 
  double local_C_cyt = C_cyt*scaling_factor;
  //cout << "Compute_J_on active" << endl

  J_on = k_buff_bind*local_C_cyt*buff_unbound;
  J_off = k_buff_unbind*buff_bound; 
  buff_c_dT = k_buff_bind*local_C_cyt*buff_unbound - k_buff_unbind*buff_bound;
  vec_buff_bound.push_back(vec_buff_bound[buff_counter] + buff_c_dT);

  buff_diff = J_off - J_on;
  buff_diff = buff_diff*pow(scaling_factor,-1);
  vec_J_on.push_back(buff_diff);

  buff_counter++;
  return(buff_diff);
}

double Compute_J_serca(double C_cyt){
  //cout << "Compute_J_serca active" << endl;
  J_serca = v_serca*(pow(C_cyt,2)/(pow(C_cyt,2) + pow(K_p,2)));
  if(J_serca != J_serca){
    J_serca = 0;
  }
  vec_J_serca.push_back(J_serca);
  return(J_serca);
}

double Compute_J_ryr(double C_Cyt){
  //cout << "Compute_J_ryr active1" << endl;

  double local_C_cyt = C_Cyt; // this is here in case we want to scale 

  w_inf = ((K_a/pow(local_C_cyt,4)) + 1 + (pow(local_C_cyt,3)/K_b))/((1/K_c) + (K_a/pow(local_C_cyt,4)) + 1 + (pow(local_C_cyt,3)/K_b)); 
  cout << ((K_a/pow(local_C_cyt,4)) + 1 + (pow(local_C_cyt,3)/K_b)) << "/" << ((1/K_c) + (K_a/pow(local_C_cyt,4)) + 1 + (pow(local_C_cyt,3)/K_b)) << endl;

  //cout << "Compute_J_ryr active2" << endl;

  tau_w = w_inf/K_d;
  w_dt = (w_inf - vec_w[w_counter])/tau_w;

  //cout << "Compute_J_ryr active3" << endl;

  vec_w.push_back(vec_w[w_counter] + w_dt);
  P_open = (vec_w[w_counter + 1]*((1 + pow(local_C_cyt,3))/K_b))/((K_a/pow(local_C_cyt,4)) + 1 + (pow(local_C_cyt,3)/K_b));
  J_ryr = (v_rel*P_open + v_leak)*(C_er - local_C_cyt);

  if(J_ryr != J_ryr){
    J_ryr = 0;
  }

  // cout << w_inf << endl;
  // cout << w << endl;
  // cout << vec_w[w_counter] << endl;
  // cout << P_open << endl;

  w_counter++;

  vec_J_ryr.push_back(1000*J_ryr);
  return(J_ryr);
}

// int reset_vecs(int x){
//     return(0);
// }

// double Piezo_screen(double channel_number, double time_max){
//   E_Ca = PotentialE(0.0024, 0.0000001, 2); //taking it to be 100nM inside, 2.4mM outside, also this outputs in millivolts
//   double Charge_temp = 0;
//   double Charge_temp_dt = 0;
//   double Conc_temp = 0;
//   for(double x = 0; x <= channel_number; x += 10){
//     for(double y = 0; y <= time_max; y += 0.001){ // time is in milliseconds, so step with 1/1000 of ms seems good
//       Charge_temp_dt = x*(E_Ca/1000)*G_Piezo_single*0.000001; // the 0.0000001 was added to convert from seconds to 0.001 ms 
//       Charge_temp += Charge_temp_dt;
      
//       Conc_temp = (Charge_temp/F)/((4/3)*M_PI*pow(0.00175,3)); //the radius, in this case, is in centimeters. I do not know why, but this seems to work better
//       //cout << x*(E_Ca/1000)*G_Piezo_single*0.001 << endl;
//       //cout << ((4/3)*M_PI*pow(0.00000175,3)) << endl;
//       //cout << Conc_temp << endl;
//       //cout << Charge_temp/F << endl;

//       // C = C + C_dt
//       if((Conc_temp >= 0.0000000022) && (Conc_temp <= 0.0000000023)){
//         vec_time.push_back(y);
//         vec_channel_number.push_back(x);
//         //cout << "occured" << endl;
//       }
//     }
//     Charge_temp = 0;
//     Conc_temp = 0; 
//     Charge_temp_dt = 0;
//   }
//   return(0);
// }

double Piezo_Channel(double potential){
  //cout << "Piezo_Channel active" << endl;
  int stochastic = stochastic_opening(generator);
  int open_temp;
  double P_open_temp;

  //P_open_temp = (vec_num_open[open_counter]*(1/exp((delta_T*0.001)/T_Piezo)));
  
  //cout << 1/(exp((delta_T*0.001)/T_Piezo)) << endl;

  open_temp = vec_num_open[open_counter] + abs(stochastic);

  // cout << P_open_temp << endl;
  // cout << exp((delta_T*0.001)/T_Piezo) << endl;

  if (open_temp >= N_Piezo_channels){
    open_temp = N_Piezo_channels;
  }

  int local_tau = 1.6/delta_T;
  //cout << local_tau << endl;
  //cout << open_counter % local_tau << endl;
  if (!(open_counter % local_tau)){

    //cout << open_counter << endl;
    open_temp = open_temp*0.9048; //this is the time constant of Piezo, so it is not relevant on a micro s scale 
  }

  // if ((open_counter >= 5000) && (open_counter <= 5005)){
  //   open_temp = N_Piezo_channels; //this is the time constant of Piezo, so it is not relevant on a micro s scale 
  // }


  //cout << open_temp << endl;

  vec_num_open.push_back(open_temp);  
  vec_num_closed.push_back(1 - open_temp);
    
  G_Piezo_total = open_temp*G_Piezo_single;

  Piezo_current = (potential*G_Piezo_total)/2; //dividing by 2 because Ca is 2+ charge (affecting current by factor of 2)
  vec_Piezo_current.push_back(Piezo_current);

  open_counter++;
  return(Piezo_current);
}

double Calcium_concentration(double time_range, double delta_T){
  E_Ca = PotentialE(0.0024, 0.0000001, 2);
  vec_w.push_back(0);
  vec_num_open.push_back(0);
  vec_num_closed.push_back(N_Piezo_channels);
  vec_buff_bound.push_back(0);
  vec_Ca_conc.push_back(0.0012); //supposed to be 120nM
  // cone_circumference = 2*M_PI*cone_radius;
  // cone_cross_area = M_PI*pow(cone_radius,2);

  // Ca_c_dT = D_diff_Ca*Ca_c_dT_dT + J_ipr + (cone_circumference/cone_cross_area)*(J_in - J_pm) + Compute_J_ryr(C_cyt) - Compute_J_serca(C_cyt) + Compute_J_on(C_cyt);

  double scaling_factor = 1;

  for(double i = 0; i <= time_range; i += delta_T){
    C_cyt = vec_Ca_conc[Ca_counter]; 

    //cout << C_cyt << endl;

    Ca_c_dT = delta_T*(scaling_factor*Piezo_Channel(E_Ca) + scaling_factor*Compute_J_ryr(C_cyt) - scaling_factor*Compute_J_serca(C_cyt) + Compute_J_on(C_cyt));
    vec_Ca_conc.push_back(vec_Ca_conc[Ca_counter] + Ca_c_dT);
    Ca_counter++;
  }

  return(0);
}

double voltage_output(double x)
{
    //reset_vecs(0);
    //Piezo_screen(1000, 10);
    Calcium_concentration(10, delta_T);

    // for (int i = 0; i < 3; i++)
    // {
    //     cout << i << endl; 
    //     Piezo_Channel(i);
    //     reset_vecs(0);
    // }

    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    vector<int> sizes;

    sizes.insert(sizes.begin(),vec_Ca_conc.size());
    sizes.insert(sizes.begin(),vec_Piezo_current.size());
    sizes.insert(sizes.begin(),vec_J_on.size());
    sizes.insert(sizes.begin(),vec_J_ryr.size());
    sizes.insert(sizes.begin(),vec_J_serca.size());
    sizes.insert(sizes.begin(),vec_num_open.size());

    sort(sizes.begin(), sizes.end());
    int max_size = sizes.back();
    
    cout << max_size << endl;

    bool bool_Ca_conc;
    bool bool_Piezo_current;
    bool bool_J_on;
    bool bool_J_ryr;
    bool bool_J_serca;
    bool bool_num_open;

    //cout << "Break point 4" << endl;

    myfile << "Ca_concentration,Buffering,Piezo_current,J_ryr,J_serca,Piezo_open\n";
    for (int i = 0; i < max_size; i++)
    {
        //cout << "Break point 5" << endl;
        bool_Ca_conc = (vec_Ca_conc.size() > i) ? true : false;
        bool_Piezo_current = (vec_Piezo_current.size() > i) ? true : false;
        bool_J_on = (vec_J_on.size() > i) ? true : false;
        bool_J_ryr = (vec_J_ryr.size() > i) ? true : false;
        bool_J_serca = (vec_J_serca.size() > i) ? true : false;
        bool_num_open = (vec_num_open.size() > i) ? true : false;

        //cout << "Break point 6" << endl;

        if(bool_Ca_conc) myfile << 100000*vec_Ca_conc[i] << ",";
        if(!bool_Ca_conc) myfile << ",";
        if(bool_J_on) myfile << 100000*vec_J_on[i] << ",";
        if(!bool_J_on) myfile << ",";
        if(bool_Piezo_current) myfile << 100000*vec_Piezo_current[i] << ",";
        if(!bool_Piezo_current) myfile << ",";
        if(bool_J_ryr) myfile << 100000*vec_J_ryr[i] << ",";
        if(!bool_J_ryr) myfile << ",";
        if(bool_J_serca) myfile << -1*100000*vec_J_serca[i]<< ",";;
        if(!bool_J_serca) myfile << ",";
        if(bool_num_open) myfile << vec_num_open[i];


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
