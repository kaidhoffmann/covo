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
        const std::vector < int > & numb_jk,
        const std::vector < double > & pos,
        const std::vector < std::vector < double > > & pos_limits,
        const std::vector < double > & Lcell);
    
   
    
    std::vector<double>rand_vec_sphere(double radius);
    std::vector<double> rand_vec_box(std::vector<double> x_lim, std::vector<double> y_lim, std::vector<double> z_lim);

    public:

    struct object {
        std::vector <double> pos;
        std::vector <double> vec_a;
        std::vector <double> vec_b;
    };
        
    struct sample {
        std::vector < object > obj;
        std::vector < double > cent;
        std::vector < std::vector < double > > edge;
    };
    
   
    std::vector < std::vector < double > > limits(sample smp);
    
    std::vector < sample > samp;
    
    sample input;
    
    sample random;
    
    void read( const parameters p, const std::string filename);
    
    void normalize_vectors();
   
    void make_samples_cart(std::vector < int > & numb_jk);
    
    void delete_input();

    void make_random(const parameters p);

    void write_input(const parameters p, const std::string filename);

};

#endif // CATALOGUE_H
