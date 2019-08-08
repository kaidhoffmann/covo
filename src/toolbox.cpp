// Kai Hoffmann
// Institute for Computational Science (ICS)
// University of Zurich
// 2019

#include <vector>
#include <cmath>
#include "toolbox.h"


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

