// Kai Hoffmann
// Institute for Computational Science (ICS)
// University of Zurich
// 2019

#ifndef TOOLBOX_H
#define TOOLBOX_H


#include <vector>
#include <cmath>


//commit massage: "move sum_vectors_int, sum_vectors_double, subtract_vectors_int, subtract_vectors_double, divide_vectors_double" from correlation.h to toolbox.h


// absolute value of vector
double vabs(const std::vector<double> vec);

//add integer vector "to_add" to integer vector "total_sum"
void sum_vectors_int(const parameters p, std::vector < int > & total_sum, std::vector < int > & to_add);

// add double vector "to_add" to double vector "total_sum"
void sum_vectors_double(const parameters p, std::vector < double > & total_sum, std::vector < double > & to_add);

// subtract integer vector "to_subtract" from integer vector "total_sum"
void subtract_vectors_int(const parameters p, std::vector < int > & total_sum, std::vector < int > & to_subtract);

// subtract double vector "to_subtract" from double vector "total_sum" 
void subtract_vectors_double(const parameters p, std::vector < double > & total_sum, std::vector < double > & to_subtract);

// devide elements in double vectors "numerator" and "denominator"
void divide_vectors_double(const parameters p, std::vector < double > & numerator, std::vector < double > & denominator);


// distance of two Nd positions
double distance(const std::vector <double> & pos_1,  const std::vector <double> & pos_2);
double distance_new(const std::vector <double> & pos_1,  const std::vector <double> & pos_2);

// convert cartesian to spehrical coordinates
std::vector<double> cart_to_sphere(const std::vector<double> pos_cart);


// convert spehrical coordinates to cartesian
std::vector<double> sphere_to_cart(const std::vector<double> pos_sphere);



#endif // TOOLBOX_H
