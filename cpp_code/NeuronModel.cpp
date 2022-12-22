#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <vector>

using namespace std;

/*
IDEAS:
1) Adding dependence on lipids would be very cool. I.e., sphingomyelin is neutral and on 
extracellcular side, some lipids have negative charge and reside on the cytoplasmic side. 
*/

vector<double> vec1, vec2, vec3, vec4, vecE, vecG;
vector<double> vec_K_I, vec_Na_I, vec_V_dt;
double E, I, V, R = 0;
vector<int> vec_channel_density;

//CONSTANTS %%%%%%%%%%%%%%%%
double log_convert = 2.303; // to convert from ln to log10
double R_constant = 8.135;
double F = 96840;
double body_temp = 310.15;
//%%%%%%%%%%%%%%%%%%%%%%%%%%

double radius = 0.0001; //100um
double resistivity = 100; //this is in Ohms apparently

int Z = 1;
double G1 = 0;
double channel_resistance, channel_current, channel_voltage;
double V_rest, V_disp, R_inp, G_inp;
double a_n, b_n, a_h, b_h, a_m, b_m;
double n_dynamic, m_dynamic, h_dynamic, n_inf, m_inf, h_inf;
vector<double> vec_n, vec_m, vec_h;
double tau_n, tau_m, tau_h;
double K_I, Na_I;
double G_sum, C_sum;

double injected_current;

double Na_G = 120;
double K_G = 36; 
double Leak_G = 0.3;
double Leak_E = 10.6;

vector<string> grid_vec1;
vector<double> grid_vec2;
int grid_rows;

int Read_file_col(int x) {
  ifstream inFile;
  inFile.open("Gridtest.csv");
  string test1, testa, testb;
  /*
  if (inFile.is_open()) {
  //THIS IS USED TO ENSURE FILE IS OPENED PROPERLY
    cout << "\n Grid file has been opened" << endl;
  } else {
    cout << "NO GRID FILE HAS BEEN OPENED" << endl;
  }
  */
  while (!inFile.eof()) {
      while (getline(inFile, test1, '\n')) {
        grid_vec1.push_back(test1);
        stringstream ss(grid_vec1[0]);
        while (ss.good()) {
          getline(ss, testb, ',');
          //grid_vec2.push_back(stod(testb));
        }
      }
    }
/*
cout << "grid_vec1 size:" << grid_vec1.size() << endl;
for (size_t i = 0; i < grid_vec1.size(); i++) {
    cout << "grid_vec1 " << grid_vec1[i] << "," << endl;
  }
*/
//grid_rows = grid_vec2.size()/grid_vec1.size();
//cout << "THIS IS GRID ROW : " << grid_rows; 
  return (0);
}

int Questions(int x){
  //cout << "COMMANDLINE TEST!!" << endl; //Update me if you wish to test CommandLine
  cout << "\n Questions to ask Professor Mori:" << endl;
  cout << "1. Do you consider 3D structure? See comment for more. For instance:" << endl;
  cout << "---00--------------     ---00---------00---" << endl;
  cout << "      case 1                   case 2" << endl;
  cout << "--------------00---     -------------------" << endl;
  cout << "In case 1 the first ion channel likely can not activate the second. \n" << endl;
  //Presumably if the charge difference is only very slight and just at the membrane, if another VGCC is on the opposite end of the cross section, maybe it will not be activated by the action potential?
  //Secondly, having a larger diameter will increase the "wire resistance." This seems to be linked to the comment above.
  cout << "2. Dynamically changing with time?" << endl;
  //I am wondering if people just average all of the sodium channels, potassium channels, etc. 
  cout << "3. Wire model does not consider spacing or density of ion channels? Voltage divider does not care about distance." << endl;
  cout << "4. Conductance and R_input needs to be re-written by the proportion of channels open, correct? So why does the model not include this?" << endl;
  return(0);
}

int Open_file(int x) {
  ifstream inFile;
  inFile.open("Input.csv");

  if (inFile.is_open()) {
  //THIS IS USED TO ENSURE FILE IS OPENED PROPERLY
    cout << "\n File has been opened" << endl;
  } else {
    cout << "NO FILE HAS BEEN OPENED" << endl;
  }
  return (0);
}

int Wire_resistance(double x){
  double wire_resistance, wire_length, cross_sectionA;
  wire_length = 0.001; //1mm
  cross_sectionA = M_PI*radius*radius;
  wire_resistance = (resistivity*wire_length)/cross_sectionA;
  //cout << "\n THE WIRE RESISTANCE IS: " << wire_resistance << endl; //show wire resistance
  return(wire_resistance);
}

int Read_file(int x) {
  //THIS CHECKS THAT THE FILE HAS BEEN READ IN PROPERLY
  //cout << "\n" << "READING INTO VECTORS" << endl;
  ifstream inFile;
  inFile.open("Input.csv");
  string test1, test2, test3;
  string testa, testb, testc;

  while (!inFile.eof()) {

    // READING IN LINE ONE
    // EVIDENTLY THE ROWS NEED TO BE SAME LENGTH
    while (getline(inFile, test1, '\n')) {
      //cout << "THIS IS TEST 1: " << test1 << endl;
      stringstream ss(test1);
      while (ss.good()) {
        getline(ss, testa, ',');
        //cout << testa << endl;
        vec1.push_back(stod(testa));
      }
      /*To describe what the above code is doing, it takes the 
      input from inFile as the variable test1 (which is a string),
      with each input separate by a comma then and converts it 
      to a "stringstream" which can be separated via the commas 
      into the variable testa. Testa is then fed into the vector*/

      /*Also, test1 is fed in per each row of the csv file, which 
      is denoted by the \n. Not so sure how to make this less 
      ugly though, as each line needs to read in in a nested manner
      in the current file. Which is quite annoying. Certainly, 
      there is a better way.*/

      // READING IN LINE TWO
      // successive lines seem to need to be imbeded within the preceeding While loop
      
      while (getline(inFile, test2, '\n')) {
        stringstream ss(test2);
        while (ss.good()) {
          getline(ss, testb, ',');
          vec2.push_back(stod(testb));
        }

        // READING IN LINE THREE
        while (getline(inFile, test3, '\n')) {
          stringstream ss(test3);
          while (ss.good()) {
            getline(ss, testc, ',');
            vec3.push_back(stod(testc));
          }
        }
        // add more vectors here
      }
    }
  }
//DO NO DELETE!!!
//THIS SECTION CAN BE USED TO CHECK VECTORS
/*
  for (size_t i = 0; i < vec1.size(); i++) {
    cout << vec1[i] << ",";
  }
  cout << endl;
  for (size_t i = 0; i < vec2.size(); i++) {
    cout << vec2[i] << ",";
  }
  cout << endl;
  for (size_t i = 0; i < vec3.size(); i++) {
    cout << vec3[i] << ",";
  }
  cout << endl;
*/
  return (0);
}

//THIS METHOD IS USED TO TEST CALCULATIONS WITH VECTORS 
//BEFORE OFFICIALLY INTEGRATING INTO CODE
/*
int Test_calculations(int x) {
  cout << endl;
  double temp = 0;
  temp = vec1[0] + vec2[0];
  cout << "TEST FUNCTION " << temp << endl;
  for (int i = 0; i < vec1.size(); i++) {
    vec4.push_back(vec1[i] + vec2[i]);
    cout << "NEW VEC pos: " << i << "=" << vec4[i] << endl;
  }

  return (0);
}
*/

double PotentialE(double out, double in) {
  double E = 1000 * log_convert * (R_constant * body_temp) / (F * Z) * log10(out / in);
  if (isinf(E)) {
    E = 0;
  }
  //THIS IS USED TO TEST POTENTIAL CALCULATION
  //cout << "THIS IS E in mV: " << E << endl;
  return (E);
}

double C_sum_G_sum(int x){
for(int i = 0; i < vecG.size(); i++){
      C_sum += (vecG[i]*vecE[i]);
      G_sum += (vecG[i]);
      //cout << vecG[i] << " times " << vecE[i] << endl;
    }
    return(0);
}

double Resting_voltage(int x){
    V_rest = C_sum/G_sum;
    //cout << "C sum : " << C_sum << endl; 
    //cout << "RESTING V: " << V_rest << endl;
  return(V_rest);
}

double Input_G_R(int x){
  G_inp = 0;
    for(int i = 0; i < vecG.size(); i++){
      G_inp += vecG[i];
    }
    R_inp = 1/G_inp;
    return (R_inp);
}

double Voltage_displacement(double x, double injected){
    double magnitude;
    double V_temp_rest = Resting_voltage(0);
    double R_temp_inp = Input_G_R(0);
    V_disp = V_temp_rest + (injected*R_temp_inp);
    //cout << "INJECTED : " << injected << endl;
    //cout << "V_TEMP_REST: " << V_temp_rest << endl;
    magnitude = injected*R_temp_inp; 
    //cout << "R = " << R_temp_inp << endl;
    //cout << "Displaced voltage: " << V_disp << endl;
    //cout << "Magnitude of displacement: " << magnitude << endl;
    return (V_disp);
}

double dynamical_h(double V){
  double h;
  a_h = 0.07*exp(-V/20);
  b_h = ((1)/(exp((30-V)/10) - 1));
  h_inf = a_h/(a_h + b_h);
  tau_h = 1/(a_h + b_h);
  h_dynamic = (h_inf - h)/tau_h;
  h = h_dynamic;
  vec_h.push_back(h);
  //cout << "h: " << h << endl;
  return(h_dynamic);
}

double dynamical_n(double V){
  double n;
  a_n = 0.01*((10 - V)/(exp((10-V)/10) - 1));
  b_n = 0.125*exp(-V/80);
  n_inf = a_n/(a_n + b_n);
  tau_n = 1/(a_n + b_n);
  n_dynamic = (n_inf - n)/tau_n;
  n = n_dynamic;
  vec_n.push_back(n);
  //cout << "n: " << n << endl;
  return (n_dynamic);
}

double dynamical_m(double V){
  double m;
  a_m = 0.01*((25 - V)/(exp((25-V)/10) - 1));
  b_m = 4*exp(-V/18);
  m_inf = a_m/(a_m + b_m);
  tau_m = 1/(a_m + b_m);
  m_dynamic = (m_inf - m)/tau_m;
  m = m_dynamic;
  vec_m.push_back(m);
  //cout << "m: " << m << endl;
  return (m_dynamic);
}

double Proportion_open(int a, int b, int c, double injected){
  double V_temp = Voltage_displacement(0, injected);
  double m_temp, h_temp, n_temp = 1;
  if(a != 0){
    m_temp = dynamical_m(V_temp);
  }
  if(b != 0){
    h_temp = dynamical_h(V_temp);
  }
  if(c != 0){
    n_temp = dynamical_n(V_temp);
  }
  double temp_P = pow(m_temp,a)*pow(h_temp,b)*pow(n_temp,c);
  return (temp_P);
}

double Generic_CurrentI(double G, double E, int a, int b, int c, double injected){
    //V is the voltage of the membrane
    //G is max conductance of channel
    //E is potential of channel
    //P is proportion of channels open
    double P_temp = Proportion_open(a,b,c, injected);
    double V_temp = Voltage_displacement(0, injected);
    channel_current = P_temp*G*(V_temp - E);
    //cout << "P: " << P_temp << endl;
    //cout << "V: " << V_temp << endl;
    //cout << "G: " << G << endl;
    //cout << "E: " << E << endl;
    //cout << "\n";
    //cout << "Calculated current: " << channel_current << endl;
    return(channel_current);
}

double Potassium_I(double x, double injected){
  double I_temp, E_temp;
  double Out_temp = 6;
  double In_temp = 120;
  E_temp = PotentialE(Out_temp, In_temp);
  I_temp = Generic_CurrentI(K_G, E_temp, 0, 0, 4, injected);
  //cout << "I_temp: " << I_temp << endl;
  return(I_temp);
}

double Sodium_I(double x, double injected){
  double I_temp, E_temp;
  double Out_temp = 120;
  double In_temp = 6;
  E_temp = PotentialE(Out_temp, In_temp);
  I_temp = Generic_CurrentI(Na_G, E_temp, 3, 1, 0, injected);
  //cout << "I_temp: " << I_temp << endl;
  return(I_temp);
}

double Writing_vectors(int x){
  //I am just putting this here to organize, but in the future this may be useful
  vecG.push_back(Na_G);
  vecG.push_back(K_G);
  vecG.push_back(Leak_G);
  vecE.push_back(PotentialE(120,10));
  vecE.push_back(PotentialE(10,120));
  //cout << "Csum unit " << PotentialE(120,10) << endl;;
  vecE.push_back(Leak_E);
  C_sum_G_sum(0);
  return(0);
}

double current_injector(double x){
  ofstream myfile;
  myfile.open("Currents.csv");
  myfile << "\n Potassium,";
  //cout << vecE[0] << " AND " << vecE[1] << endl;
  double I_temp;
  for(int i = 0; i < vec1.size(); i++){
    injected_current = vec1[i];
    //cout << "HERE IS THE INJECTED: " << injected_current << endl;
    I_temp = Potassium_I(0, injected_current);
    myfile << I_temp << ",";
    vec_K_I.push_back(I_temp);
  }
  myfile << "\n Sodium,";
  for(int i = 0; i < vec1.size(); i++){
    injected_current = vec1[i];
    I_temp = Sodium_I(0, injected_current);
    myfile << I_temp << ",";
    vec_Na_I.push_back(I_temp);
  }
  myfile.close();
  myfile.open("DynamicVariables.csv");
  myfile << "\n N, ";
  for(int i = 0; i < vec_n.size(); i++){
    myfile << vec_n[i] << ","; 
  }
  myfile << "\n M, ";
  for(int i = 0; i < vec_m.size(); i++){
    myfile << vec_m[i] << ","; 
  }
  myfile << "\n H, ";
  for(int i = 0; i < vec_h.size(); i++){
    myfile << vec_h[i] << ","; 
  }
return(0);
}

double Writing_wire_dimensions(double density, int spacing, int wire_length){
  //density in units of ion channels per um, lets say
  //nodes of ranvier are 1um wide, for reference
  ofstream myfile;
  myfile.open("Output.csv");
  for(int i = 1; i <= wire_length; i++){
        if(i%spacing != 0){
          vec_channel_density.push_back(density);
        }
        else{
          vec_channel_density.push_back(0);
        }
  }
  myfile.close();
  return(0);
}

int Output_file(int x) {
  //cout << endl << "OUTPUT FILE TESTING" << endl;
  ofstream myfile;
  myfile.open("Output.csv");
  myfile << "Test3.\n";
  myfile << "Potential,";
  for (int i = 0; i < vec1.size(); i++) {
    vecE.push_back(PotentialE(vec1[i], vec2[i]));
    myfile << vecE[i] << ",";
    //cout << "\n VECE I =" << vecE[i] << endl;
  }
  myfile << "\n";
  //testing contents of vec
  myfile << "\n Testing wire dimensions,";
  for(int i = 0; i < vec_channel_density.size(); i++){
    myfile << vec_channel_density[i] << ",";
  }
  myfile << "\n";
  myfile.close();
  return 0;
}

double Propagating_AP(double x){
  double V_dt = 0; 
  current_injector(0);
  for(int i = 0; i < vec1.size(); i++){
  V_dt = ((radius/(2*resistivity)) - vec1[i] - vec_K_I[0] - vec_Na_I[0])/C_sum; 
  vec_V_dt.push_back(V_dt);
  cout << "C SUM " << C_sum << endl;
  cout << "V dt : " << vec_V_dt[i] << endl;
  //WHAT IS THE TIME DEPENDENCE IN THIS CASE?
  }
  return(0); 
}

int main(void) {
  cout << "Ran" << endl;
  //Questions(0);
  Open_file(0);
  Read_file(0);
  Writing_vectors(0);
  //Wire_resistance(0);
  //Test_calculations(0);
  //Writing_wire_dimensions(10,3,100);
  //Read_file_col(0);
  Propagating_AP(0);
}
