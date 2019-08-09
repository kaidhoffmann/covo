// Kai Hoffmann
// Institute for Computational Science (ICS)
// University of Zurich
// 2019

#ifndef TOOLBOX_H
#define TOOLBOX_H


#include <vector>
#include <cmath>
#include <algorithm>


// ==========================================================
// subtract vector "to_subtract" from vector "total_sum" 
// ==========================================================
template <typename T> inline void subtract_vectors(std::vector < T > & total_sum, std::vector < T > & to_subtract){
    std::transform (
        total_sum.begin(), total_sum.end(),
        to_subtract.begin(),total_sum.begin(),
        std::minus<T>());
}


// ==========================================================
//add vector "to_add" to vector "total_sum"
// ==========================================================
template <typename T> inline void sum_vectors(std::vector < T > & total_sum, std::vector < T > & to_add){
    
    std::transform (
            total_sum.begin(), total_sum.end(),
            to_add.begin(),total_sum.begin(),
            std::plus<T>());
        
}


// ==========================================================
// devide elements in vectors "numerator" and "denominator"
// ==========================================================
template <typename T> inline void divide_vectors(std::vector < T > & numerator, std::vector < T > & denominator){
    std::transform (
            numerator.begin(), numerator.end(),
            denominator.begin(),numerator.begin(),
            std::divides<T>());
}


// ==========================================================
// absolute value of vector
// ==========================================================
template <typename T> inline double vabs(const std::vector<T> vec) {

    T abs_val = 0;

    for(int i=0; i<vec.size(); i++){
        abs_val += pow(vec[i],2);
    }
    
    return sqrt(abs_val);
}


// ==========================================================
// distance between two Nd positions
// ==========================================================
template <typename T> inline double distance(const std::vector <T> & pos_1,  const std::vector <T> & pos_2){
    
    T dist=0;
    
    for(int i=0; i<pos_1.size();i++){
        dist += pow(pos_2[i] - pos_1[i],2);
    }
    
    return sqrt(dist);
};



// ==========================================================
// convert spehrical coordinates to cartesian
// ==========================================================
template <typename T> inline std::vector<T> sphere_to_cart(const std::vector<T> pos_sphere) {

    std::vector<T> pos_cart;
    
    pos_cart.push_back( pos_sphere[0] * sin(pos_sphere[2]) * cos(pos_sphere[1]) );
    pos_cart.push_back( pos_sphere[0] * sin(pos_sphere[2]) * sin(pos_sphere[1]) );
    pos_cart.push_back( pos_sphere[0] * cos(pos_sphere[2]) );
    
    return pos_cart;
}



// ==========================================================
// convert cartesian to spehrical coordinates
// ==========================================================
template <typename T> inline std::vector<T> cart_to_sphere(const std::vector<T> pos_cart) {

    std::vector<T> pos_sphere;
    
    T r = vabs(pos_cart);
    
    pos_sphere.push_back( r );
    pos_sphere.push_back( atan2(pos_cart[1],  pos_cart[0]) );
    pos_sphere.push_back( acos(pos_cart[2] / r) );
    
    return pos_sphere;
}




#endif // TOOLBOX_H
