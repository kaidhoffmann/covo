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
#include <cmath>
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
    
    auto_limits = (get_param(fname, "auto_limits") == "true");
    
    type_subsample = get_param(fname, "type_subsample");
  
    r_lim = extract_numbers_double( get_param(fname, "r_lim") );
    theta_lim = extract_numbers_double( get_param(fname, "theta_lim") );
    phi_lim = extract_numbers_double( get_param(fname, "phi_lim") );

    x_lim = extract_numbers_double( get_param(fname, "x_lim") );
    y_lim = extract_numbers_double( get_param(fname, "y_lim") );
    z_lim = extract_numbers_double( get_param(fname, "z_lim") );
    
    //convert angular limits from degree to radians
    for(int i=0; i < 2 ;i++){
        theta_lim[i] *= M_PI/180.;
        phi_lim[i] *= M_PI/180.;
    }

  
    // parameters for healpix sampling
    nside = std::stoi( get_param(fname, "nside") );
    nrad = std::stoi( get_param(fname, "nrad") );
    
    //number of jk samples per axis for cartesian subsampling
    numb_jk_cart = extract_numbers_int( get_param(fname, "numb_jk_cart") );
  
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
    make_rand = (get_param(fname, "make_rand") == "true");
    rand_seed = std::stod( get_param(fname, "rand_seed") );
    fname_rand = get_param(fname, "fname_rand");
    
    x_lim_rand = x_lim;
    y_lim_rand = y_lim;
    z_lim_rand = z_lim;
    
    r_lim_rand = r_lim;
    theta_lim_rand = theta_lim;
    phi_lim_rand = phi_lim;

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
    
    if(numb_jk_cart.size() !=dim){
        std::cerr<<"# ##### ERROR: " << numb_jk_cart.size() << " dimensions for numb_jk_cart sampling, but " << dim << " columns for position #####"<<std::endl;
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
    
    
    if(nside < 1){
        //TODO: check order
        std::cerr<<"# ##### ERROR: nside must be power of 2 and > 0" << std::endl;
        eishockey = false;
    }
    
    
    if(nrad < 1){
        std::cerr<<"# ##### ERROR: nrad must be > 0" << std::endl;
        eishockey = false;        
    }
    
    
    if(!auto_limits){
        if(theta_lim[0] < 0 || theta_lim[1] < 0 || theta_lim[0] > M_PI || theta_lim[1] > M_PI){
            std::cerr<<"# ##### ERROR: theta_lim must be in range [0, 180] degree" << std::endl;
            eishockey = false;        
        }

        if(r_lim[0] < 0 || r_lim[1] < 0){
            std::cerr<<"# ##### ERROR: r_lim must be > 0" << std::endl;
            eishockey = false;
        }

    }    
    
    //TODO:
    //check if limits are set correctly (min-max, not max, min) to avoid negative dr, dphi, dtheta
    
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
    for(int i = 0; i < numb_jk_cart.size(); i++){ Njk *= numb_jk_cart[i]; }
    
    std::cout << "# type_subsample: ";
    std::cout << type_subsample << std::endl;
       
    std::cout << "# auto_limits: " << auto_limits << std::endl;

    if(type_subsample=="cartesian"){

        if(!auto_limits){
            
            std::cout<<"# x_lim: ";
            for(int i = 0; i < x_lim.size(); i++){ std::cout<<x_lim[i]<<" "; }
            std::cout<<std::endl;

            std::cout<<"# y_lim: ";
            for(int i = 0; i < y_lim.size(); i++){ std::cout<<y_lim[i]<<" "; }
            std::cout<<std::endl;

            std::cout<<"# z_lim: ";
            for(int i = 0; i < z_lim.size(); i++){ std::cout<<z_lim[i]<<" "; }
            std::cout<<std::endl;
        }
        
        std::cout << "# JK samples per axis: ";
        for(int i = 0; i < numb_jk_cart.size(); i++){ std::cout<<numb_jk_cart[i]<<" "; }
        std::cout << " => " << Njk << " samples in total" << std::endl;        
        std::cout<<std::endl;
    }
    
    
    if(type_subsample=="healpix"){
        
        if(!auto_limits){
            
            std::cout<<"# r_lim: ";
            for(int i = 0; i < r_lim.size(); i++){ std::cout<<r_lim[i]<<" "; }
            std::cout<<std::endl;

            std::cout<<"# theta_lim: ";
            for(int i = 0; i < theta_lim.size(); i++){ std::cout<<theta_lim[i] * 180./M_PI <<" "; }
            std::cout<<std::endl;

            std::cout<<"# phi_lim: ";
            for(int i = 0; i < phi_lim.size(); i++){ std::cout<<phi_lim[i] * 180./M_PI<<" "; }
            std::cout<<std::endl;
        }
        
        std::cout << "# nside: " << nside << std::endl;
        std::cout << "# nrad: " << nrad << std::endl;
        std::cout<<std::endl;
    }
    

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
    
    if(make_rand){
        
        std::cout<<"# rand_seed = "<<rand_seed<<std::endl;
        std::cout<<"# fname_rand = "<<fname_rand<<std::endl;
        
        std::cout<<"# numb_rand = "<<numb_rand<<std::endl;
        
        std::cout<<"# x_lim_rand: ";
        for(int i = 0; i < x_lim_rand.size(); i++){ std::cout<<x_lim_rand[i]<<" "; }
        std::cout<<std::endl;

        std::cout<<"# y_lim_rand: ";
        for(int i = 0; i < y_lim_rand.size(); i++){ std::cout<<y_lim_rand[i]<<" "; }
        std::cout<<std::endl;

        std::cout<<"# z_lim_rand: ";
        for(int i = 0; i < z_lim_rand.size(); i++){ std::cout<<z_lim_rand[i]<<" "; }
        std::cout<<std::endl;

    }
    
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
