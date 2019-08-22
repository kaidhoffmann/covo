// Kai Hoffmann
// Institute for Computational Science (ICS)
// University of Zurich
// 2019

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <vector>
#include <sstream> 

class parameters{
    
    std::string get_param(std::string fname_param_file, std::string param_name);

    std::vector<int> extract_numbers_int(std::string str);
    std::vector<double> extract_numbers_double(std::string str);
    
    public:
    
    short verbose;
    
    std::string fname_cat_1, fname_cat_2, dir_out, delim_in, delim_out;
    
    std::string fname_out, fname_out_prefix, fname_out_suffix, fname_out_extention;
    
    std::vector <int> cols_pos, cols_vec_a, cols_vec_b;
    
    std::vector <double> x_lim, y_lim, z_lim;
    std::vector <double> r_lim, theta_lim, phi_lim;
    
    bool auto_limits;
    
    std::string type_subsample;

    double nside, nrad;
    
    std::vector <int> numb_jk_cart;
    
    void read(std::string fname);
    
    void print();
    
    bool check();
        
    double r_min, r_max;
    int numb_bin;
    
    bool lg_bins;
    
    bool
    r12_v1a, r12_v1b,
    r12_v2a, r12_v2b,
    v1a_v2a, v1b_v2b,
    v1a_v2b, v1b_v2a;
    
    bool header_out;
    
    bool make_rand;
    int numb_rand;
    double rand_seed;    
    std::string fname_rand;
    std::vector <double> x_lim_rand, y_lim_rand, z_lim_rand;
    std::vector <double> r_lim_rand, theta_lim_rand, phi_lim_rand;
    
};

#endif // PARAMETERS_H
