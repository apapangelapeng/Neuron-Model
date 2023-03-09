#include <fstream>
#include <vector>

using namespace std;
void output_data(const char* path, vector<double> vec_V, vector<double> vec_K_I, vector<double> vec_Na_I, vector<double> vec_L_I, vector<double> vec_Cl_I, vector<double> vec_HCN_I,char* vec_name)
{

    ofstream create_file(path);
    ofstream myfile;
    myfile.open(path);



    vector<int> sizes;
    /*cout<<vec_V.size()<<endl;
    cout<<vec_N.size()<<endl;
    cout<<vec_tiny_N.size()<<endl;*/
    sizes.insert(sizes.begin(),vec_V.size());
    sizes.insert(sizes.begin(),vec_K_I.size());
    sizes.insert(sizes.begin(),vec_Na_I.size());
    sizes.insert(sizes.begin(),vec_L_I.size());
    sizes.insert(sizes.begin(),vec_Cl_I.size());
    sizes.insert(sizes.begin(),vec_HCN_I.size());
    sort(sizes.begin(), sizes.end());
    int max_size = sizes.back();
    

    bool V; 
    bool K_I; 
    bool Na_I; 
    bool L_I; 
    bool Cl_I;
    bool HCN_I;
    
    //cout << "Break point 4" << endl;

    myfile << "V,K_I,Na_I,L_I,Cl_I,"<<vec_name<<"\n";
    for (int i = 0; i < max_size; i++)
    {
    
        //cout << "Break point 5" << endl;
        V = (vec_V.size() > i) ? true : false;
        K_I = (vec_V.size() > i) ? true : false;
        Na_I = (vec_V.size() > i) ? true : false;
        L_I = (vec_V.size() > i) ? true : false;
        HCN_I = (vec_V.size() > i) ? true : false;

        //cout << "Break point 6" << endl;

        if(V) myfile << vec_V[i] <<"," ;
        if(!V) myfile <<"," ;
        if(K_I) myfile << vec_K_I[i] << ",";
        if(!K_I) myfile << ",";
        if(Na_I) myfile << vec_Na_I[i] << ",";
        if(!Na_I) myfile <<",";
        if(L_I) myfile << vec_L_I[i] << ",";
        if(!L_I) myfile <<",";
        if(Cl_I) myfile << vec_Cl_I[i] << ",";
        if(!Cl_I) myfile <<",";
        if (vec_name =="HCN_I"){
        if(HCN_I) myfile << vec_HCN_I[i];}

        //if(!dn_V_0) myfile << vec_tiny_N[i];

        myfile << "\n";
        //cout << "Break point 6" << endl;

    }
    myfile.close();

}