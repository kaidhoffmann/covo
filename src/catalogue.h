// Kai Hoffmann
// Institute for Computational Science (ICS)
// University of Zurich
// 2019

#ifndef CATALOGUE_H
#define CATALOGUE_H

#include <string>
#include <vector>
#include "parameters.h"

class catalogue{
    
                      
    int pos_to_ID_cart(
        const std::vector < int > & numb_jk_cart,
        const std::vector < double > & pos,
        const std::vector < std::vector < double > > & pos_limits,
        const std::vector < double > & Lcell);
    
    std::vector<double> rand_vec_sphere(double radius);
    std::vector<double> rand_vec_shell(const std::vector<double> & intrsq_lim, const std::vector<double> & cos_theta_lim, const std::vector<double> & phi_lim);
    std::vector<double> rand_vec_box(std::vector<double> x_lim, std::vector<double> y_lim, std::vector<double> z_lim);

    public:

    struct object {
        std::vector <double> pos;
        std::vector <double> vec_a;
        std::vector <double> vec_b;
    };
        
    struct sample {
        long hp_ID;
        int bin_rad;
        std::vector < object > obj;
        std::vector < double > cent;
        std::vector < std::vector < double > > edge;
    };
    
   
    std::vector < std::vector < double > > find_limits_cart(sample smp);
    std::vector < std::vector < double > > find_limits_sphere(sample smp);
    
    void get_pos_limits(const parameters p);

    void show_pos_limits();

    void cut_input(const parameters p);
    
    std::vector < std::vector < double > > pos_limits_cart, pos_limits_sphere;
    
    
    std::vector < sample > samp;
    
    sample input;
    
    sample random;
    
    void read( const parameters p, const std::string filename);
    
    void normalize_vectors();
   
    void make_samples_cart(const parameters p);

    void make_samples_healpix(const parameters p);
    
    void delete_input();

    void make_random_box(const parameters p);

    void make_random_shell(const parameters p);

    void write_input(const parameters p, const std::string filename);

};


void cut_overlap(const parameters p, catalogue & cat_1, catalogue & cat_2);


#endif // CATALOGUE_H
