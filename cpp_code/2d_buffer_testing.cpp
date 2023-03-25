#include "2d_Calcium_Dynamics.h"

const char *path1="../data_files/2d_Piezo_Channel.csv";
const char *path2="../data_files/2d_Piezo_Channel_avg.csv";

default_random_engine generator;
normal_distribution<double> stochastic_opening(0,0.025);

int reset_vecs(int x){
    vec_time.clear();
    vec_buff_bound.clear();
    vec_w.clear();

    return(0);
}

double Compute_J_diffusion(int time, int x, int y) { 
    double R = 1;
    double total_diffusion, x_diffusion, y_diffusion;
    if((x >= 1) && (x <= x_max - 1)){
        x_diffusion = (1/R)*(vec_time[time][y][x + 1] - (2 * vec_time[time][y][x]) + vec_time[time][y][x - 1]);
    }
    if((y >= 1) && (y <= y_max - 1)){
        y_diffusion = (1/R)*(vec_time[time][y + 1][x] - (2 * vec_time[time][y][x]) + vec_time[time][y - 1][x]);
    }
    
    if((y == 0)){
        y_diffusion = (1/R)*(vec_time[time][y + 1][x] - (vec_time[time][y][x]));
    }
    if((y == y_max)){
        y_diffusion = (1/R)*(vec_time[time][y - 1][x] - (vec_time[time][y][x]));
    }

    if((x == 0)){
        x_diffusion = (1/R)*(vec_time[time][y][x + 1] - (vec_time[time][y][x]));
    }
    if((x == x_max)){
        x_diffusion = (1/R)*(vec_time[time][y][x - 1] - (vec_time[time][y][x]));
    }

    total_diffusion = x_diffusion + y_diffusion;
  return(total_diffusion);
}


double Compute_J_on(double C_cyt, int time, int x, int y){
    buff_bound = vec_buff_bound[time][x][y]; 
    buff_unbound = vec_buff_unbound[time][x][y]; 

    double local_C_cyt = C_cyt*10000000000;

    k_buff_bind = 0.600; //for the buffer BAPTA in mM
    k_buff_unbind = 0.100;

    J_on = k_buff_bind*local_C_cyt*buff_unbound;
    J_off = k_buff_unbind*buff_bound; 
    buff_c_dT = k_buff_bind*local_C_cyt*buff_unbound - k_buff_unbind*buff_bound;
    
    //buff_c_dT = -(k_buff_unbind*buff_bound);

    vec_buff_bound[time + 1][x][y] = (vec_buff_bound[time][x][y] + buff_c_dT);
    vec_buff_unbound[time + 1][x][y] = (vec_buff_unbound[time][x][y] - buff_c_dT);

    buff_diff = J_off - J_on;
    //buff_diff = J_off;
    //buff_diff = buff_diff*pow(scaling_factor,-1);
    //vec_J_on.push_back(buff_diff);

    return(buff_diff); //*0.000001
}

double Calcium_concentration(double x){

    cout << "  High" << endl;
    
    double divs = (x_max + 1)*(y_max + 1);

    double mols_divs = 0.0000000012/divs;

    for(int i = 0; i <= x_max; i++){
        for(int j = 0; j <= y_max; j++){
            vec_time[0][i][j] = mols_divs;
            vec_buff_bound[0][i][j] = buff_bound;
            vec_buff_unbound[0][i][j] = buff_unbound;
        }
    }

    //cout << "break point 2" << endl;

    double scaling_factor;
    double divide = (y_max + 1) * (x_max + 1);
    scaling_factor = 1/divide; 
    double avg_temp, avg_buff_temp;

    for(int time_temp = 0; time_temp <= time_max; time_temp++){
        for(int i = 0; i <= x_max; i++){
            for(int j = 0; j <= y_max; j++){

                C_cyt = vec_time[time_temp][i][j];

                Ca_c_dT = delta_T*(scaling_factor*Compute_J_diffusion(time_temp, j, i) + scaling_factor*Compute_J_on(C_cyt, time_temp, i, j));
                
                vec_time[time_temp + 1][i][j] = vec_time[time_temp][i][j] + 40*Ca_c_dT;

                // what will be better is to solve for the number of moles in each cube, then use that to calculate overall concentration
                avg_temp += vec_time[time_temp][i][j];
                avg_buff_temp += vec_buff_bound[time_temp][i][j];
                //cout << avg_temp << endl;
            }
        }
        vec_average.push_back(avg_temp/(divide));
        vec_buff_average.push_back(avg_buff_temp/(divide));
        //cout << avg_temp/(x_max*y_max) << endl;
        avg_temp = 0;
        avg_buff_temp = 0;
    }

    cout << "    to Low" << endl;

    return(0);
}


double output_file(double x)
{
    ofstream create_file(path1);
    ofstream myfile;
    myfile.open(path1);

    Calcium_concentration(0);

    vector<int> sizes;

    sizes.insert(sizes.begin(),vec_time.size());

    sort(sizes.begin(), sizes.end());
    int max_size = sizes.back();
    
    cout << "Vector length: " << max_size << endl;

    for(int time = 0; time <= time_max; time++){
        for (int x = 0; x <= x_max; x++){
            for (int y = 0; y <= y_max; y++){  
                if(y < y_max) {
                    myfile << vec_time[time][x][y] << ",";
                }  
                else {
                    myfile << vec_time[time][x][y];
                }
            }
        myfile << "\n";
        }
    myfile << "\n";
    }

    reset_vecs(0);

    myfile.close();

    return (x);
}

double output_avg_file(double x)
{
    ofstream create_file(path2);
    ofstream myfile;
    myfile.open(path2);

    vector<int> sizes;

    sizes.insert(sizes.begin(),vec_average.size());
    sizes.insert(sizes.begin(),vec_buff_average.size());

    sort(sizes.begin(), sizes.end());
    int max_size = sizes.back();
    
    cout << "Vector length: " << max_size << endl;

    bool bool_average;
    bool bool_buff_average;

    //cout << "Break point 4" << endl;

    myfile << "Average,buff_average\n";

    for (int i = 0; i < max_size; i++)
    {
        //cout << "Break point 5" << endl;
        bool_average = (vec_average.size() > i) ? true : false;
        bool_average = (vec_buff_average.size() > i) ? true : false;

        //cout << "Break point 6" << endl;

        if(bool_average) myfile << vec_average[i] << ",";
        if(!bool_average) myfile << ",";
        if(bool_average) myfile << vec_buff_average[i];



        myfile << "\n";
        //cout << "Break point 6" << endl;
    }
    reset_vecs(0);

    myfile.close();

    return (x);
}

int main(void) {
  cout << "Begin" << endl;
  output_file(0);
  output_avg_file(0);
  cout << "End" << endl;
}
