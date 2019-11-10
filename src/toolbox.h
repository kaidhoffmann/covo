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
// 1D distance in periodic box
// ==========================================================
template <typename T> inline T periodic_distance(T dr, T Lbox, T r_max){
    if(dr >= r_max){
        dr -= Lbox;
    }else{
        if(dr <= -r_max){
            dr += Lbox;
        }
    }
    return dr;
}



// ==========================================================
// convert spehrical coordinates to cartesian
// ==========================================================
template <typename T> inline std::vector<T> sphere_to_cart(const std::vector<T> pos_sphere) {

    std::vector<T> pos_cart;
    
    pos_cart.push_back( pos_sphere[0] * sin(pos_sphere[1]) * cos(pos_sphere[2]) );
    pos_cart.push_back( pos_sphere[0] * sin(pos_sphere[1]) * sin(pos_sphere[2]) );
    pos_cart.push_back( pos_sphere[0] * cos(pos_sphere[1]) );
    
    return pos_cart;
}



// ==========================================================
// convert cartesian to spehrical coordinates
// ==========================================================
template <typename T> inline std::vector<T> cart_to_sphere(const std::vector<T> pos_cart) {

    
    T r = vabs(pos_cart);
    T theta = acos(pos_cart[2] / r);
    T phi = atan2(pos_cart[1],  pos_cart[0]);
    if(phi < 0){ phi += 2*M_PI; }
    
    std::vector<T> pos_sphere {r, theta, phi};
    
    return pos_sphere;
}


// ==========================================================
// check if x is power of y
// source: https://www.geeksforgeeks.org/check-if-a-number-is-power-of-another-number/
// ==========================================================
template <typename Ta, typename Tb > inline bool isPower(Ta x, Tb y){
    
    // logarithm function to calculate value 
    int res1 = int(log(y) / log(x)); 
    double res2 = double(log(y) / log(x)); // Note : this is double 
  
    // compare to the result1 or result2 both are equal 
    return (res1 == res2); 
} 



#endif // TOOLBOX_H
