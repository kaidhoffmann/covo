// Kai Hoffmann
// Institute for Computational Science (ICS)
// University of Zurich
// 2019

#include <string>
#include <iostream>
#include <fstream>
#include <sstream> 
#include <algorithm>
#include <vector>
#include "parameters.h"



//===========================================================
// read parameters from file
//===========================================================
void parameters::read(std::string fname){
    
    verbose = std::stoi(get_param(fname, "verbose"));
    
    fname_cat_1 = get_param(fname, "fname_cat_1");
    fname_cat_2 = get_param(fname, "fname_cat_2");

    delim_in = get_param(fname, "delim_in");   
    delim_out = get_param(fname, "delim_out");
    
    if(delim_in == "tab") delim_in = "\t";
    if(delim_out == "tab") delim_out = "\t";

    if(delim_in == "space") delim_in = " ";
    if(delim_out == "space") delim_out = " ";

    
    dir_out = get_param(fname, "dir_out");
    fname_out = get_param(fname, "fname_out");
    fname_out_prefix = get_param(fname, "fname_out_prefix");
    fname_out_suffix = get_param(fname, "fname_out_suffix");
    fname_out_extention = get_param(fname, "fname_out_extention");

    header_out = (get_param(fname, "header_out") == "true");
    
    //column numbers with vector components in input catalogue
    cols_pos = extract_numbers_int( get_param(fname, "cols_pos") );
    cols_vec_a = extract_numbers_int( get_param(fname, "cols_vec_a") );
    cols_vec_b = extract_numbers_int( get_param(fname, "cols_vec_b") );
    
    type_subsample = get_param(fname, "type_subsample");
    
    //number of jk samples per axis
    numb_jk = extract_numbers_int( get_param(fname, "numb_jk") );
  
    
    //binning variables
    numb_bin = std::stoi( get_param(fname, "numb_bin") );
    r_min = std::stod( get_param(fname, "r_min") );
    r_max = std::stod( get_param(fname, "r_max") );
    lg_bins = (get_param(fname, "lg_bins") == "true");

  
    //vector cobinations to compute
    r12_v1a = (get_param(fname, "r12_v1a") == "true");
    r12_v1b = (get_param(fname, "r12_v1b") == "true");
    
    r12_v2a = (get_param(fname, "r12_v2a") == "true");
    r12_v2b = (get_param(fname, "r12_v2b") == "true");
    
    v1a_v2a = (get_param(fname, "v1a_v2a") == "true");
    v1b_v2b = (get_param(fname, "v1b_v2b") == "true");
    
    v1a_v2b = (get_param(fname, "v1a_v2b") == "true");
    v1b_v2a = (get_param(fname, "v1b_v2a") == "true");
    
    
    //variable for generating random catalogue
    numb_rand = std::stoi( get_param(fname, "numb_rand") );
    xlim_rand = extract_numbers_double( get_param(fname, "xlim_rand") );
    ylim_rand = extract_numbers_double( get_param(fname, "ylim_rand") );
    zlim_rand = extract_numbers_double( get_param(fname, "zlim_rand") );
  
    make_rand = (get_param(fname, "make_rand") == "true");
    rand_seed = std::stod( get_param(fname, "rand_seed") );
    fname_rand = get_param(fname, "fname_rand");
};



//===========================================================
// health checks
//===========================================================
bool parameters::check(){
    
    bool eishockey = true;
    
    int dim = cols_pos.size();

    if(cols_vec_a.size() !=dim){
        std::cerr<<"# ##### ERROR: " << cols_vec_a.size() << " columns for vector a, but " << dim << " columns for position #####"<<std::endl;
        eishockey = false;
    }

    if(cols_vec_b.size() !=dim){
        std::cerr<<"# ##### ERROR: " << cols_vec_b.size() << " columns for vector b, but " << dim << " columns for position #####"<<std::endl;
        eishockey = false;
    }
    
    if(numb_jk.size() !=dim){
        std::cerr<<"# ##### ERROR: " << numb_jk.size() << " dimensions for numb_jk sampling, but " << dim << " columns for position #####"<<std::endl;
        eishockey = false;
    }

    if(lg_bins && r_min <=0){
        std::cerr<<"# ##### ERROR: set r_min > 0 when lg_bins = true" << std::endl;
        eishockey = false;        
    }
    
    if(delim_in.size() > 1){
        std::cerr<<"# ##### ERROR: delim_in must be one character" << std::endl;
        eishockey = false;        
    }
    
    if(type_subsample !="cartesian" && type_subsample !="healpix"){
        std::cerr<<"# ##### ERROR: type_subsample must be cartesian or healpix" << std::endl;
        eishockey = false;        
    }
    
        std::cout<<std::endl;

    
    return eishockey;
};



//===========================================================
// print parameters
//===========================================================
void parameters::print(){
    
    std::cout << std::endl << "# ================ PARAMETERS ================" << std::endl;

    
    std::cout << "# fname_cat_1: " << fname_cat_1 << std::endl;
    std::cout << "# fname_cat_2: " << fname_cat_2 << std::endl;
    
    std::cout<< "# fname_out_prefix: " << fname_out_prefix << std::endl;
    std::cout<< "# fname_out_suffix: " << fname_out_suffix << std::endl;
    std::cout<< "# fname_out_extention: " << fname_out_extention << std::endl;
    
    std::cout << "# dir_out: " << dir_out << std::endl;
    
    std::cout << "# delim_in: " << delim_in << std::endl;
    std::cout << "# delim_out: " << delim_out << std::endl;

    std::cout << "# header_out: " << header_out << std::endl;

    std::cout<<std::endl;

    
    int Njk=1;
    for(int i = 0; i < numb_jk.size(); i++){ Njk *= numb_jk[i]; }
    
    
    std::cout << "# type_subsample: ";
    std::cout << type_subsample << std::endl;
    std::cout<<std::endl;
    
    
    std::cout << "# JK samples per axis: ";
    for(int i = 0; i < numb_jk.size(); i++){ std::cout<<numb_jk[i]<<" "; }
    std::cout << " => " << Njk << " samples in total" << std::endl;
    std::cout<<std::endl;


    std::cout<<"# columns in input catalogue: "<< std::endl;

    std::cout<<"# - position vector: ";
    for(int i = 0; i < cols_pos.size(); i++){ std::cout<<cols_pos[i]<<" "; }
    std::cout<<std::endl;
    
    std::cout<<"# - vector a: ";
    for(int i = 0; i < cols_vec_a.size(); i++){ std::cout<<cols_vec_a[i]<<" "; }
    std::cout<<std::endl;
    
    std::cout<<"# - vector b: ";
    for(int i = 0; i < cols_vec_b.size(); i++){ std::cout<<cols_vec_b[i]<<" "; }
    std::cout<<std::endl;
    std::cout<<std::endl;
    
    
    std::cout<<"# r_min = "<<r_min<<std::endl;
    std::cout<<"# r_max = "<<r_max<<std::endl;
    std::cout<<"# numb_bin = "<<numb_bin<<std::endl;
    std::cout<<"# lg_bins = "<<lg_bins<<std::endl;
    std::cout<<std::endl;
    
    
    std::cout<<"# r12_v1a = "<<r12_v1a<<std::endl;
    std::cout<<"# r12_v1b = "<<r12_v1b<<std::endl;
    std::cout<<"# r12_v2a = "<<r12_v2a<<std::endl;
    std::cout<<"# r12_v2b = "<<r12_v2b<<std::endl;
    
    std::cout<<"# v1a_v2a = "<<v1a_v2a<<std::endl;
    std::cout<<"# v1b_v2b = "<<v1b_v2b<<std::endl;
    std::cout<<"# v1a_v2b = "<<v1a_v2b<<std::endl;
    std::cout<<"# v1b_v2a = "<<v1b_v2a<<std::endl;
    std::cout<<std::endl;
    
    
    std::cout<<"# make_rand = "<<make_rand<<std::endl;
    std::cout<<"# rand_seed = "<<rand_seed<<std::endl;
    std::cout<<"# fname_rand = "<<fname_rand<<std::endl;
    
    
    std::cout<<"# numb_rand = "<<numb_rand<<std::endl;
    
    std::cout<<"# xlim_rand: ";
    for(int i = 0; i < xlim_rand.size(); i++){ std::cout<<xlim_rand[i]<<" "; }
    std::cout<<std::endl;

    std::cout<<"# ylim_rand: ";
    for(int i = 0; i < ylim_rand.size(); i++){ std::cout<<ylim_rand[i]<<" "; }
    std::cout<<std::endl;

    std::cout<<"# zlim_rand: ";
    for(int i = 0; i < zlim_rand.size(); i++){ std::cout<<zlim_rand[i]<<" "; }
    std::cout<<std::endl;

    std::cout << "# ============================================" << std::endl;
};



//===========================================================
// read string variables from parameter file
//===========================================================
std::string parameters::get_param(std::string fname_param_file, std::string param_name) {
    
    std::ifstream cFile (fname_param_file);

    if(cFile.is_open()){
        
        bool param_found = false;
        std::string param = "nan";
        std::string line;
        
        while(getline(cFile, line)){
            
            line.erase(std::remove_if(line.begin(), line.end(), isspace),line.end());

            if(line[0] == '#' || line.empty())
                continue;
            auto delimiterPos = line.find("=");
            auto var = line.substr(0, delimiterPos);
            auto val = line.substr(delimiterPos + 1);
            //std::cout << var << " = " << val << std::endl;
            
            if(param_name == var){
                param = val;
                param_found = true;
            }
            
        }
        
        if(param_found){
            return param;
        }else{
            std::cerr << "# #### ERROR [parameters::get_param]: parameter " <<param_name << " not found in "<< fname_param_file <<" ####" << std::endl;
            exit (EXIT_FAILURE);
        }
        
    }else {
        std::cerr << "# #### ERROR [parameters::get_param]: Couldn't open parameter file: " <<fname_param_file << " ####" << std::endl;
        exit (EXIT_FAILURE);
    }
    
};



//===========================================================
// get vector with int numbers from string
//===========================================================
std::vector<int> parameters::extract_numbers_int(std::string str){
//from https://stackoverflow.com/questions/41820639/c-extracting-the-integer-in-a-string-with-multiple-delimiters

    std::vector<int> numbers;
    std::stringstream ss(str);
    
    if( str.length() != 0 ){
        while( !ss.eof() ){
            int number;
            ss>>number;
            numbers.push_back(number);
            ss.get();
        }
    }
 
    return numbers;
}



//===========================================================
// get vector with double numbers from string
//===========================================================
std::vector<double> parameters::extract_numbers_double(std::string str){
//from https://stackoverflow.com/questions/41820639/c-extracting-the-integer-in-a-string-with-multiple-delimiters

    std::vector<double> numbers;
    std::stringstream ss(str);
    
    if( str.length() != 0 ){
        while( !ss.eof() ){
            double number;
            ss>>number;
            numbers.push_back(number);
            ss.get();
        }
    }
 
    return numbers;
}
