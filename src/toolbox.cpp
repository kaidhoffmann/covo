// Kai Hoffmann
// Institute for Computational Science (ICS)
// University of Zurich
// 2019

#include <vector>
#include <cmath>
#include "parameters.h"
#include "toolbox.h"
//#include <functional>   // std::divides
#include <algorithm>

// ==========================================================
// absolute value of vector
// ==========================================================
double vabs(const std::vector<double> vec) {

    double abs_val = 0;

    for(int i=0; i<vec.size(); i++){
        abs_val += pow(vec[i],2);
    }
    
    return sqrt(abs_val);
}

// ==========================================================
// distance of two Nd positions
// ==========================================================
double distance(const std::vector <double> & pos_1,  const std::vector <double> & pos_2){
    
    double dist=0;
    
    for(int i=0; i<pos_1.size();i++){
        dist += pow(pos_2[i] - pos_1[i],2);
    }
    
    return sqrt(dist);
};


// ==========================================================
// distance of two Nd positions
// ==========================================================
double distance_new(const std::vector <double> & pos_1,  const std::vector <double> & pos_2){
    
    const std::vector <double> dist_vec;// = subtract_vectors_double(pos_2, pos_1);

    return vabs(dist_vec);
    
};


// ==========================================================
// convert cartesian to spehrical coordinates
// ==========================================================
std::vector<double> cart_to_sphere(const std::vector<double> pos_cart) {

    std::vector<double> pos_sphere;
    
    double r = vabs(pos_cart);
    
    pos_sphere.push_back( r );
    pos_sphere.push_back( atan2(pos_cart[1],  pos_cart[0]) );
    pos_sphere.push_back( acos(pos_cart[2] / r) );
    
    return pos_sphere;
}


// ==========================================================
// convert spehrical coordinates to cartesian
// ==========================================================
std::vector<double> sphere_to_cart(const std::vector<double> pos_sphere) {

    std::vector<double> pos_cart;
    
    pos_cart.push_back( pos_sphere[0] * sin(pos_sphere[2]) * cos(pos_sphere[1]) );
    pos_cart.push_back( pos_sphere[0] * sin(pos_sphere[2]) * sin(pos_sphere[1]) );
    pos_cart.push_back( pos_sphere[0] * cos(pos_sphere[2]) );
    
    
    return pos_cart;
}




// ==========================================================
// add integer vector "to_add" to integer vector "total_sum"
// ==========================================================
void sum_vectors_int(std::vector < int > & total_sum, std::vector < int > & to_add){
    
    std::transform (
            total_sum.begin(), total_sum.end(),
            to_add.begin(),total_sum.begin(),
            std::plus<int>());
        
}



// ==========================================================
// add double vector "to_add" to double vector "total_sum"
// ==========================================================
void sum_vectors_double(std::vector < double > & total_sum, std::vector < double > & to_add){
    
    std::transform (
            total_sum.begin(), total_sum.end(),
            to_add.begin(),total_sum.begin(),
            std::plus<double>());
        
}



// ==========================================================
// subtract integer vector "to_subtract" from integer vector "total_sum" 
// ==========================================================
void subtract_vectors_int(std::vector < int > & total_sum, std::vector < int > & to_subtract){
    
    std::transform (
            total_sum.begin(), total_sum.end(),
            to_subtract.begin(),total_sum.begin(),
            std::minus<int>());
        
}



// ==========================================================
// subtract double vector "to_subtract" from double vector "total_sum" 
// ==========================================================
void subtract_vectors_double(std::vector < double > & total_sum, std::vector < double > & to_subtract){
    
    std::transform (
            total_sum.begin(), total_sum.end(),
            to_subtract.begin(),total_sum.begin(),
            std::minus<double>());
}



// ==========================================================
// devide elements in double vectors "numerator" and "denominator"
// ==========================================================
void divide_vectors_double(std::vector < double > & numerator, std::vector < double > & denominator){
    std::transform (
            numerator.begin(), numerator.end(),
            denominator.begin(),numerator.begin(),
            std::divides<double>());
}




